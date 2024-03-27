from doctest import set_unittest_reportflags
import os
from typing import List, Iterable, Dict, Tuple, Collection, Sequence
from copy import deepcopy
import numpy as np
import pandas as pd
from Bio import SeqIO
import bean as be
from bean.framework.Edit import Edit, Allele
from bean.framework.AminoAcidEdit import (
    AminoAcidEdit,
    CodingNoncodingAllele,
)
from bean.annotate.utils import get_cds_seq_pos_from_gene_name, find_overlap, revcomp
import logging
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n",
    datefmt="%a, %d %b %Y %H:%M:%S",
    stream=sys.stderr,
    filemode="w",
)
error = logging.critical
warn = logging.warning
debug = logging.debug
info = logging.info

BASE_SET = {"A", "C", "T", "G"}
reverse_map = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "-": "-"}

codon_map = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


class RefBaseMismatchException(Exception):
    pass


# This function is adopted from https://github.com/gpp-rnd/be-validation-pipeline/blob/main/notebooks/01_BEV_allele_frequencies.ipynb
def _translate(seq: str, codon_map: Dict[str, str]):
    """Translate seq as is."""
    #     if not seq: # if remove_introns returned False -> possible splice site mutation
    if seq == "UTR":
        return "Possible UTR mutation"
    if seq == "intron":
        return "Possible intronic mutation"
    aa = ""
    i = 0
    frame = 1
    while i < len(seq):
        substring = ""
        while frame <= 3:
            if i >= len(seq):
                return aa + ">"
            if seq[i] != "-":
                substring += seq[i]
                frame += 1
            i += 1
        if len(substring) == 3:
            aa = aa + "-" if "N" in substring else aa + codon_map[substring]
        frame = 1  # reset frame
    return aa


def _parse_range(chr_range: str) -> Tuple[str, int, int]:
    """Parse 'chr:start-end' into chr,start,end tuple"""
    chrom, gen_range = chr_range.split(":")
    start, end = gen_range.split("-")
    start = int(start)
    end = int(end)
    return (chrom, start, end)


def _parse_description(desc_str: str):
    """Parse description line of fasta file to get range."""
    sstr = desc_str.split(" ")
    strand = None
    genome_range = None
    for s in sstr:
        if s.startswith("range="):
            genome_range = _parse_range(s[6:])
        if s.startswith("strand="):
            strand = s[-1]
            assert strand in ["+", "-"]
    if strand is None:
        strand = "+"
    return genome_range, strand


def get_cds_seq_pos_from_fasta(fasta_file_name: str) -> Tuple[List[str], List[int]]:
    """Obtain tuple of lists, first with nucleotide and second with genomic position of the nucleotide."""
    exons = list(SeqIO.parse(fasta_file_name, "fasta"))
    translated_seq = []
    genomic_pos = []
    for exon in exons:
        (exon_chrom, exon_start, exon_end), strand = _parse_description(
            exon.description
        )
        for i, nt in enumerate(str(exon.seq)):
            if nt.islower():
                continue
            translated_seq.append(nt)
            genomic_pos.append(exon_start + i)
    return (exon_chrom, translated_seq, genomic_pos, strand)


def _translate_single_codon(
    codon: List[str],
) -> str:  # nt_seq_string: str, aa_pos: int) -> str:
    """Translate `aa_pos`-th codon of `nt_seq_string`."""
    if len(codon) != 3:
        print("reached the end of CDS, frameshift.")
        return "/"
    try:
        codon = "".join(codon)
        return codon_map[codon]
    except KeyError:
        if codon[-1] == "N" and codon[0] in BASE_SET and codon[1] in BASE_SET:
            aa_set = {codon_map[codon[:2] + N] for N in BASE_SET}
            if len(aa_set) == 1:
                return next(iter(aa_set))
            raise ValueError(f"warning: no matching aa with codon {codon}")
        else:
            raise ValueError(f"Cannot translate codon due to ambiguity: {codon}")


class CDS:
    gene_info_dict = {}
    mismatch_tuples = []  # (gene, pos, ref, edit)

    def __init__(self):
        self.edited_aa_pos = set()
        self.edits_noncoding = set()
        self.edited_nt: List[str] = []
        self.nt: List[str] = []
        self.strand: int = 1
        self.gene_name: str = ""
        self.chrom: str = ""
        self.translated_seq: List[str] = []
        self.pos: np.ndarray = None
        self.genomic_pos: Sequence = []

    @classmethod
    def from_fasta(cls, fasta_file_name, gene_name, suppressMessage=True):
        cds = cls()
        # if fasta_file_name is not None:
        #     if not suppressMessage:
        #         raise ValueError("No fasta file provided as reference: using LDLR")
        #     fasta_file_name = os.path.dirname(be.__file__) + "/annotate/ldlr_exons.fa"
        if gene_name not in cls.gene_info_dict:
            chrom, translated_seq, genomic_pos, strand = get_cds_seq_pos_from_fasta(
                fasta_file_name
            )
            cls.gene_info_dict[gene_name] = {
                "chrom": chrom,
                "translated_seq": translated_seq,
                "genomic_pos": genomic_pos,
                "strand": strand,
            }
        cds.gene_name = gene_name
        cds.strand = cls.gene_info_dict[gene_name]["strand"]
        cds.chrom = cls.gene_info_dict[gene_name]["chrom"]
        cds.translated_seq = deepcopy(cls.gene_info_dict[gene_name]["translated_seq"])
        cds.genomic_pos = cls.gene_info_dict[gene_name]["genomic_pos"]
        cds.nt = cds.gene_info_dict[gene_name]["translated_seq"]
        cds.pos = np.array(cds.genomic_pos)
        cds.edited_nt = deepcopy(cds.nt)
        return cds

    @classmethod
    def from_gene_name(cls, gene_name, ref_version: str = "GRCh38"):
        cds = cls()
        if gene_name not in cls.gene_info_dict:
            chrom, translated_seq, genomic_pos, strand = get_cds_seq_pos_from_gene_name(
                gene_name, ref_version
            )
            cls.gene_info_dict[gene_name] = {
                "chrom": chrom,
                "translated_seq": translated_seq,
                "genomic_pos": genomic_pos,
                "strand": strand,
            }
        cds.gene_name = gene_name
        cds.strand = cls.gene_info_dict[gene_name]["strand"]
        cds.chrom = cls.gene_info_dict[gene_name]["chrom"]
        cds.translated_seq = deepcopy(cls.gene_info_dict[gene_name]["translated_seq"])
        cds.genomic_pos = cls.gene_info_dict[gene_name]["genomic_pos"]
        cds.nt = cds.gene_info_dict[gene_name]["translated_seq"]  # in sense strand
        # print(cds.gene_name + ":" + "".join(cds.nt))
        # if cds.strand == -1:
        #     cds.nt = revcomp(cds.nt)
        cds.pos = np.array(cds.genomic_pos)
        cds.edited_nt = deepcopy(cds.nt)
        return cds

    # def translate(self):
    #     if self.strand == -1:
    #         self.edited_nt = revcomp(self.edited_nt)
    #     self.aa = _translate(self.edited_nt, codon_map)

    def _get_relative_nt_pos(self, absolute_pos):
        """0-based relative position"""
        nt_relative_pos = np.where(self.pos == absolute_pos)[0]
        assert len(nt_relative_pos) <= 1, nt_relative_pos
        return nt_relative_pos.astype(int).item() if nt_relative_pos else -1

    def _edit_pos_to_aa_pos(self, edit_pos):
        """0-based nt position. Adds in sense direction, needs to be reversed for antisense gene"""
        nt_relative_pos = self._get_relative_nt_pos(edit_pos)
        if nt_relative_pos != -1:
            self.edited_aa_pos.add(nt_relative_pos // 3)

        return nt_relative_pos

    def edit_single(self, edit_str):
        """Add a mutation induced by a single `edit_str` to CDS.
        For the negative CDS, nt and edited_nt are antisense."""
        edit = Edit.from_str(edit_str)
        rel_pos = self._edit_pos_to_aa_pos(edit.pos)
        if edit.strand == "-":
            ref_base = reverse_map[edit.ref_base]
            alt_base = reverse_map[edit.alt_base]
        else:
            ref_base = edit.ref_base
            alt_base = edit.alt_base
        if rel_pos == -1:  # position outside CDS
            self.edits_noncoding.add(edit)
        elif self.nt[rel_pos] != ref_base:
            if ref_base != "-":
                CDS.mismatch_tuples.append(
                    (
                        (self.gene_name if hasattr(self, "gene_name") else ""),
                        rel_pos,
                        self.nt[rel_pos],
                        edit,
                    )
                )
                raise RefBaseMismatchException(
                    f"{self.gene_name + ';' if hasattr(self, 'gene_name') else ''}ref:{self.nt[rel_pos]} at pos {rel_pos}, got edit {edit}. Gene sequence: {''.join(self.nt)}"
                )
        else:
            self.edited_nt[rel_pos] = alt_base
        if alt_base == "-":  # frameshift
            # self.edited_nt.pop(rel_pos)
            if self.strand == 1:
                self.edited_aa_pos.update(list(range(rel_pos, len(self.nt) // 3)))
            else:
                self.edited_aa_pos.update(list(range(0, rel_pos)))

    def edit_allele(self, allele_str):
        if isinstance(allele_str, Allele):
            edit_strs = allele_str.edits
        else:
            edit_strs = allele_str.split(",")
        for edit_str in edit_strs:
            self.edit_single(edit_str)
        if "-" in self.edited_nt:
            self.edited_nt = [nt for nt in self.edited_nt if nt != "-"]
        if self.strand == -1:
            # Reverse logged edited positions as it was in the sense direction.
            # rev_pos = self.edited_aa_pos[::-1]
            self.edited_aa_pos = {len(self.nt) // 3 - 1 - r for r in self.edited_aa_pos}
            self.nt = revcomp(self.nt)
            self.edited_nt = revcomp(self.edited_nt)

    def get_aa_change(
        self, allele_str, include_synonymous=True
    ) -> CodingNoncodingAllele:
        """1-based amino acid editing result"""
        self.edit_allele(allele_str)
        mutations = CodingNoncodingAllele()
        mutations.nt_allele.update(self.edits_noncoding)
        for edited_aa_pos in self.edited_aa_pos:
            try:
                # print("".join(self.nt))
                # print((3 * edited_aa_pos), (3 * edited_aa_pos + 3))
                ref_aa = _translate_single_codon(
                    self.nt[(3 * edited_aa_pos) : (3 * edited_aa_pos + 3)]
                )
            except ValueError as e:
                print(
                    f"Translation mismatch in translating ref for {allele_str}: {e}. Check the input .fasta or genome version used for the reporter. Ignoring this allele."
                )
                continue
            try:
                # print("".join(self.nt))
                print(self.edited_nt[(3 * edited_aa_pos) : (3 * edited_aa_pos + 3)])
                mt_aa = _translate_single_codon(
                    self.edited_nt[(3 * edited_aa_pos) : (3 * edited_aa_pos + 3)]
                )
            except ValueError as e:
                print(
                    f"Translation mismatch in translating mutated gene for {allele_str}: {e}. Check the input .fasta or genome version used for the reporter. Ignoring this allele."
                )
                continue
            except IndexError as e:
                print(f"End of gene reached by frameshift {allele_str}: {e}")
                mt_aa = ">"
            if not include_synonymous and ref_aa == mt_aa:
                continue
            mutations.aa_allele.add(
                AminoAcidEdit(edited_aa_pos + 1, ref_aa, mt_aa, gene=self.gene_name)
            )
        return mutations


def get_mismatch_df():
    return pd.DataFrame.from_records(
        CDS.mismatch_tuples, columns=["gene", "rel_pos", "ref", "edit"]
    ).drop_duplicates()


def get_cds_dict(gene_names, ref_version: str = "GRCh38"):
    return {gname: CDS.from_gene_name(gname, ref_version) for gname in gene_names}


def export_gene_info_to_json(gene_dict, write_path=".tmp_gene_info.csv"):
    gene_dfs = []
    for k, v in gene_dict.items():
        gene_df = pd.DataFrame(
            {"translated_seq": v["translated_seq"], "genomic_pos": v["genomic_pos"]}
        )
        gene_df["gene"] = k
        gene_df["chrom"] = v["chrom"]
        gene_dfs.append(gene_df)
    pd.concat(gene_dfs).to_csv(write_path)


class CDSCollection:
    """
    Represents a collection of coding sequences (CDS) for multiple genes.
    """

    unedited_cds_dict = {}

    def __init__(
        self,
        gene_names: List[str] = None,
        fasta_file_names: List[str] = None,
        suppressMessage=True,
    ):
        if fasta_file_names is None:
            if not CDSCollection.unedited_cds_dict:
                CDSCollection.unedited_cds_dict = get_cds_dict(gene_names)
        elif len(gene_names) != len(fasta_file_names):
            raise ValueError("gene_names and fasta_file_names have different lengths")
        else:
            if not CDSCollection.unedited_cds_dict:
                for gid, fasta_file in zip(fasta_file_names, gene_names):
                    CDSCollection.unedited_cds_dict[gid] = CDS.from_fasta(
                        fasta_file, fasta_file_id=gid
                    )
        CDSCollection.cds_ranges = CDSCollection.get_cds_ranges()
        export_gene_info_to_json(CDS.gene_info_dict)

    @classmethod
    def get_cds_ranges(cls):
        gids = []
        seqnames = []
        starts = []
        ends = []
        for gene_id, cds in cls.unedited_cds_dict.items():
            gids.append(gene_id)
            seqnames.append(cds.chrom)
            starts.append(cds.genomic_pos[0])
            ends.append(cds.genomic_pos[-1])
        starts = np.array(starts)
        ends = np.array(ends)
        genomic_start = np.minimum(starts, ends)
        genomic_end = np.maximum(starts, ends)
        return pd.DataFrame(
            {
                "chrom": seqnames,
                "start": genomic_start,
                "end": genomic_end,
            },
            index=gids,
        )

    def get_aa_change(
        self, allele: Allele, include_synonymous: bool = True
    ) -> CodingNoncodingAllele:  # sourcery skip: use-named-expression
        """Finds overlapping CDS and call the same function for the CDS, else return CodingNonCodingAllele with no translated allele."""
        if len(allele.edits) == 0:
            return CodingNoncodingAllele.from_alleles(nt_allele=allele)
        chrom, start, end = allele.get_range()
        overlapping_cds = find_overlap(chrom, start, end, self.cds_ranges)
        if overlapping_cds:
            return deepcopy(self.unedited_cds_dict[overlapping_cds]).get_aa_change(
                allele, include_synonymous
            )
        else:
            return CodingNoncodingAllele.from_alleles(nt_allele=allele)


def get_allele_aa_change_single_gene(
    allele: Allele,
    gene_name: str = None,
    fasta_file: str = None,
    fasta_file_id: str = None,
    include_synonymous: bool = True,
    ref_version: str = "GRCh38",
):
    """
    Obtain amino acid changes
    """
    if gene_name is None and fasta_file is not None:
        cds = CDS.from_fasta(fasta_file, fasta_file_id)
    elif gene_name is not None and fasta_file is None:
        cds = CDS.from_gene_name(gene_name, ref_version)
    else:
        raise ValueError("Only one of `gene_name` or `fasta_file` should be provided.")
    return cds.get_aa_change(allele, include_synonymous)


def get_allele_aa_change_multi_genes(
    allele: Allele, gene_names=None, fasta_file_dict=None, include_synonymous=True
):
    """
    Obtain amino acid changes for multiple provided gene exon sequences
    """
    if gene_names is None and fasta_file_dict is not None:
        cdss = CDSCollection(
            gene_names=fasta_file_dict.keys(), fasta_file_list=fasta_file_dict.values()
        )
    elif gene_names is not None and fasta_file_dict is None:
        cdss = CDSCollection(gene_names=gene_names)
    else:
        raise ValueError(
            "Only one of `gene_names` or `fasta_file_dict` should be provided."
        )
    return cdss.get_aa_change(allele, include_synonymous)


def translate_allele(
    allele: Allele,
    include_synonymous=True,
    allow_ref_mismatch=True,
    gene_name: str = None,
    gene_names: List[str] = None,
    fasta_file: str = None,
    fasta_file_dict: Dict[str, str] = None,
):
    try:
        if gene_name:
            if (
                gene_names is not None
                or (fasta_file is not None)
                or (fasta_file_dict is not None)
            ):
                raise ValueError(
                    "Only one of `gene_name`, `gene_names`, `fasta_file`, `fasta_file_dict` can be provided as argument."
                )
            return get_allele_aa_change_single_gene(
                allele,
                gene_name=gene_name,
                include_synonymous=include_synonymous,
            )
        elif gene_names:
            if (fasta_file is not None) or (fasta_file_dict is not None):
                raise ValueError(
                    "Only one of `gene_name`, `gene_names`, `fasta_file`, `fasta_file_dict` can be provided as argument."
                )
            return get_allele_aa_change_multi_genes(
                allele,
                gene_names=gene_names,
                include_synonymous=include_synonymous,
            )
        elif fasta_file:
            if fasta_file_dict is not None:
                raise ValueError(
                    "Only one of `gene_name`, `gene_names`, `fasta_file`, `fasta_file_dict` can be provided as argument."
                )
            return get_allele_aa_change_single_gene(
                allele,
                fasta_file=fasta_file,
                include_synonymous=include_synonymous,
            )
        elif fasta_file_dict:
            return get_allele_aa_change_multi_genes(
                allele,
                fasta_file_dict=fasta_file_dict,
                include_synonymous=include_synonymous,
            )
        else:
            fasta_file_name = os.path.dirname(be.__file__) + "/annotate/ldlr_exons.fa"
            return get_allele_aa_change_single_gene(
                allele,
                fasta_file=fasta_file_name,
                include_synonymous=include_synonymous,
            )
    except RefBaseMismatchException as e:
        if not allow_ref_mismatch:
            raise e
        print(e)
        return "ref mismatch"


def translate_allele_df(
    allele_df,
    include_synonymous=True,
    allow_ref_mismatch=True,
    gene_name=None,
    gene_names=None,
    fasta_file=None,
    fasta_file_dict: Dict[str, str] = None,
):
    if fasta_file and fasta_file_dict:
        raise ValueError("Both fasta file and fasta file dict provided.")
    allele_df = allele_df.copy()
    translated_alleles = allele_df.allele.map(
        lambda a: be.translate_allele(
            a,
            include_synonymous=include_synonymous,
            allow_ref_mismatch=allow_ref_mismatch,
            gene_name=gene_name,
            gene_names=gene_names,
            fasta_file=fasta_file,
            fasta_file_dict=fasta_file_dict,
        )
    )
    allele_df.insert(2, "cn_allele", translated_alleles)
    allele_df["cn_allele"] = allele_df.cn_allele.map(str)
    allele_df = allele_df.loc[
        ~allele_df.cn_allele.map(str).isin(
            ["ref mismatch", "|", "", "translation error"]
        ),
        :,
    ]
    allele_df = allele_df.drop("allele", axis=1).groupby(["guide", "cn_allele"]).sum()
    allele_df = allele_df.reset_index().rename(columns={"cn_allele": "aa_allele"})
    allele_df.aa_allele = allele_df.aa_allele.map(
        lambda s: be.CodingNoncodingAllele.from_str(s)
    )
    return allele_df


def filter_nt_allele(cn_allele: CodingNoncodingAllele, pos_include: Iterable[int]):
    """
    For CodingNoncodingAllele object, retain all amino acid mutation while filtering the nt_allele based on edit position.
    """
    cn_allele = deepcopy(cn_allele)
    edit_list = [e for e in cn_allele.nt_allele.edits if e.pos in pos_include]
    cn_allele.nt_allele = be.Allele(edit_list)
    return cn_allele


def filter_nt_alleles(cn_allele_df: pd.DataFrame, pos_include: Iterable[int]):
    """
    For CodingNoncodingAllele object, retain all amino acid mutation while filtering the nt_allele based on edit position.
    Arguments
    -- cn_allele_df (pd.DataFrame): Allele dataframe that has 'guide', 'aa_allele' as columns.
    """
    splice_only = cn_allele_df.aa_allele.map(
        lambda a: be.an.translate_allele.filter_nt_allele(a, pos_include)
    )
    alleles = cn_allele_df.copy()
    alleles = alleles.drop("aa_allele", axis=1, inplace=False)
    alleles.insert(1, "aa_allele", splice_only)
    alleles = alleles.groupby(["guide", "aa_allele"]).sum().reset_index()
    alleles = alleles.loc[alleles.aa_allele.map(bool), :]
    return alleles


def strsplit_edit(edit_str):
    if len(edit_str.split(":")) == 3:
        chrom, pos, transition = edit_str.split(":")
    elif len(edit_str.split(":")) == 2:
        pos, transition = edit_str.split(":")
        chrom = None
    else:
        raise ValueError(f"{edit_str} is not in the correct format.")
    ref, alt = transition.split(">")
    return chrom, pos, ref, alt


def annotate_edit(
    edit_info: pd.DataFrame,
    edit_col: str = "edit",
    control_tag: str = "CONTROL",
    splice_sites: Collection[
        int
    ] = None,  # TODO: may be needed to extended into multi-chromosome case
):
    """Classify edit strings into coding/noncoding and synonymous/misseinse[splicing/trunc]

    Args
    edit_info: pd.DataFrame with at least 1 column of 'edit_col', which has 'Edit' format.
    control_tag: String tag identifying non-targeting control guides so their variant signal are not aggregated.
    splice_sites: Collection of integer splice site positions. If the edit position matches the positions, it will be annotated as 'splicing'.

    """
    edit_info = edit_info.copy()
    edit_info["group"] = ""
    edit_info["int_pos"] = -1
    if "pos" not in edit_info.columns:
        edit_info["chrom"], edit_info["pos"], edit_info["ref"], edit_info["alt"] = zip(
            *(edit_info[edit_col].map(strsplit_edit))
        )
    edit_info["coding"] = ""
    edit_info.loc[edit_info.pos.map(lambda s: s.startswith("A")), "coding"] = "coding"
    edit_info.loc[
        edit_info.pos.map(lambda s: not s.startswith("A")), "coding"
    ] = "noncoding"
    if control_tag is not None:
        edit_info.loc[
            edit_info.pos.map(lambda s: control_tag in s), "group"
        ] = "negctrl"
        edit_info.loc[
            edit_info.pos.map(lambda s: control_tag in s), "coding"
        ] = "negctrl"
    edit_info.loc[
        (edit_info.coding == "noncoding") & (edit_info.group != "negctrl"), "int_pos"
    ] = edit_info.loc[
        (edit_info.coding == "noncoding") & (edit_info.group != "negctrl"), "pos"
    ].map(
        int
    )

    edit_info.loc[
        (edit_info.alt != edit_info.ref) & (edit_info.coding == "coding"), "group"
    ] = "missense"
    edit_info.loc[edit_info.alt == "*", "group"] = "trunc"
    edit_info.loc[
        (edit_info.alt == edit_info.ref) & (edit_info.coding == "coding"), "group"
    ] = "syn"
    if splice_sites is not None:
        edit_info.loc[
            edit_info.pos.isin(splice_sites.astype(str)), "group"
        ] = "splicing"
    edit_info.loc[edit_info.int_pos < -100, "group"] = "negctrl"
    edit_info.loc[edit_info.int_pos < -100, "coding"] = "negctrl"
    return edit_info
