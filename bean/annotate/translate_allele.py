import os
from typing import List, Iterable, Dict, Tuple, Collection, Optional
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
from bean.annotate.utils import get_cds_seq_pos_from_gene_name
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
    """Translate seq"""
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
    for s in sstr:
        if s.startswith("range="):
            return _parse_range(s[6:])


def _get_seq_pos_from_fasta(fasta_file_name: str) -> Tuple[List[str], List[int]]:
    """Obtain tuple of lists, first with nucleotide and second with genomic position of the nucleotide."""
    exons = list(SeqIO.parse(fasta_file_name, "fasta"))
    translated_seq = []
    genomic_pos = []
    for exon in exons:
        exon_chrom, exon_start, exon_end = _parse_description(exon.description)
        for i, nt in enumerate(str(exon.seq)):
            if nt.islower():
                continue
            translated_seq.append(nt)
            genomic_pos.append(exon_start + i)
    return (exon_chrom, translated_seq, genomic_pos)


def _translate_single_codon(nt_seq_string: str, aa_pos: int) -> str:
    """Translate `aa_pos`-th codon of `nt_seq_string`."""
    codon = "".join(nt_seq_string[aa_pos * 3 : (aa_pos * 3 + 3)])
    if len(codon) != 3:
        print("reached the end of CDS, frameshift.")
        return "/"
    try:
        return codon_map[codon]
    except KeyError:
        if codon[-1] == "N" and codon[0] in BASE_SET and codon[1] in BASE_SET:
            aa_set = {codon_map[codon[:2] + N] for N in BASE_SET}
            if len(aa_set) == 1:
                return next(iter(aa_set))
            print(f"warning: no matching aa with codon {codon}")
        else:
            print(f"Cannot translate codon due to ambiguity: {codon}")

        return "_"


class CDS:
    def __init__(self):
        self.edited_aa_pos = set()
        self.edits_noncoding = set()

    @classmethod
    def from_fasta(cls, fasta_file_name, suppressMessage=True):
        cds = cls()
        if fasta_file_name is None:
            if not suppressMessage:
                print("No fasta file provided as reference: using LDLR")
            fasta_file_name = os.path.dirname(be.__file__) + "/annotate/ldlr_exons.fa"
        try:
            cds.set_exon_fasta_name(fasta_file_name)
        except FileNotFoundError as e:
            print(os.getcwd())
            raise e
        cds.edited_nt = cds.nt.copy()
        return cds

    @classmethod
    def from_gene_name(cls, gene_name):
        cds = cls()
        cds.chrom, cds.translated_seq, cds.genomic_pos = get_cds_seq_pos_from_gene_name(
            gene_name
        )
        cds.nt = cds.translated_seq
        cds.pos = np.array(cds.genomic_pos)

    def set_exon_fasta_name(self, fasta_file_name: str):
        self.fasta_file_name = fasta_file_name
        self.chrom, self.translated_seq, self.genomic_pos = _get_seq_pos_from_fasta(
            fasta_file_name
        )
        self.nt = self.translated_seq
        self.pos = np.array(self.genomic_pos)

    def translate(self):
        self.aa = _translate(self.edited_nt, codon_map)

    def _get_relative_nt_pos(self, absolute_pos):
        nt_relative_pos = np.where(self.pos == absolute_pos)[0]
        assert len(nt_relative_pos) <= 1, nt_relative_pos
        return nt_relative_pos.astype(int).item() if nt_relative_pos else -1

    def _edit_pos_to_aa_pos(self, edit_pos):
        nt_relative_pos = self._get_relative_nt_pos(edit_pos)
        if nt_relative_pos != -1:
            self.edited_aa_pos.add(nt_relative_pos // 3)
        return nt_relative_pos

    def edit_single(self, edit_str):
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
                raise RefBaseMismatchException(
                    f"ref:{self.nt[rel_pos]} at pos {rel_pos}, got edit {edit}"
                )
        else:
            self.edited_nt[rel_pos] = alt_base
        if alt_base == "-":  # frameshift
            # self.edited_nt.pop(rel_pos)
            self.edited_aa_pos.update(list(range(rel_pos, len(self.nt) // 3)))

    def edit_allele(self, allele_str):
        if isinstance(allele_str, Allele):
            edit_strs = allele_str.edits
        else:
            edit_strs = allele_str.split(",")
        for edit_str in edit_strs:
            self.edit_single(edit_str)
        if "-" in self.edited_nt:
            self.edited_nt.remove("-")

    def get_aa_change(
        self, allele_str, include_synonymous=True
    ) -> CodingNoncodingAllele:
        """1-based amino acid editing result"""
        self.edit_allele(allele_str)
        mutations = CodingNoncodingAllele()
        mutations.nt_allele.update(self.edits_noncoding)
        for edited_aa_pos in self.edited_aa_pos:
            ref_aa = _translate_single_codon(self.nt, edited_aa_pos)
            mt_aa = _translate_single_codon(self.edited_nt, edited_aa_pos)
            if mt_aa == "_":
                return "translation error"
            if not include_synonymous and ref_aa == mt_aa:
                continue
            mutations.aa_allele.add(AminoAcidEdit(edited_aa_pos + 1, ref_aa, mt_aa))
        return mutations


def get_cds_dict(gene_names):
    return {gname: CDS.from_gene_name(gname) for gname in gene_names}


class CDSCollection:
    """
    Represents a collection of coding sequences (CDS) for multiple genes.
    """

    def __init__(
        self,
        gene_ids: List[str] = None,
        fasta_file_names: List[str] = None,
        suppressMessage=True,
    ):
        if fasta_file_names is None:
            self.cds_dict = get_cds_dict(gene_ids)
        elif len(gene_ids) != len(fasta_file_names):
            raise ValueError("gene_ids and fasta_file_names have different lengths")
        self.cds_dict = {}
        for gid, fasta_file in zip(fasta_file_names, gene_ids):
            self.cds_dict[gid] = CDS(fasta_file)
        self.cds_ranges = self.get_cds_ranges()

    def get_cds_ranges(self):
        gids = []
        seqnames = []
        starts = []
        ends = []
        for gene_id, cds in self.cds_dict.items():
            gids.append(gene_id)
            seqnames.append(cds.chrom)
            starts.append(cds.start)
            ends.append(cds.end)

        return pd.DataFrame(
            {
                "chrom": seqnames,
                "start": starts,
                "end": ends,
            },
            index=gids,
        )

    def find_overlap(
        self, chrom: str, start: int, end: int, range_df: pd.DataFrame
    ) -> Optional[str]:
        """Find overlap between query range and range_df and return ID of overlapping region in range_df."""
        if not chrom in range_df.chrom:
            return None
        overlap = range_df.loc[
            (range_df.chrom == chrom)
            & ((start <= range_df.end) | (end >= range_df.start))
        ]
        if len(overlap) == 0:
            return None
        if len(overlap) > 2:
            raise ValueError("We cannot handle overlapping genes.")
        return overlap.index.tolist()[0]

    def get_aa_change(
        self, allele: Allele, include_synonymous: bool = True
    ) -> CodingNoncodingAllele:
        """Finds overlapping CDS and call the same function for the CDS, else return CodingNonCodingAllele with no translated allele."""
        chrom, start, end = allele.get_range()
        overlapping_cds = find_overlap(chrom, start, end, self.cds_ranges)
        if overlapping_cds:
            return self.cds_dict[overlapping_cds].get_aa_change(
                allele, include_synonymous
            )
        else:
            return CodingNonCodingAllele.from_alleles(nt_allele=allele)


def get_allele_aa_change_single_gene(
    allele: Allele, gene_name=None, fasta_file=None, include_synonymous=True
):
    """
    Obtain amino acid changes
    """
    if gene_name is None and fasta_file is not None:
        cds = CDS.from_fasta(fasta_file)
    elif gene_name is not None and fasta_file is None:
        cds = CDS.from_gene_name(gene_name)
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


def annotate_edit(
    edit_info: pd.DataFrame,
    edit_col="edit",
    splice_sites: Collection[
        int
    ] = None,  # TODO: may be needed to extended into multi-chromosome case
):
    """Classify edit strings into coding/noncoding and synonymous/misseinse[splicing/trunc]

    Args
    edit_info: pd.DataFrame with at least 1 column of 'edit_col', which has 'Edit' format.
    splice_sites: Collection of integer splice site positions. If the edit position matches the positions, it will be annotated as 'splicing'.

    """
    edit_info = edit_info.copy()
    edit_info["group"] = ""
    edit_info["int_pos"] = -1
    if "pos" not in edit_info.columns:
        edit_info["pos"], transition = zip(*(edit_info[edit_col].str.split(":")))
        edit_info["ref"], edit_info["alt"] = zip(
            *(pd.Series(transition).str.split(">"))
        )
    edit_info.loc[edit_info.pos.map(lambda s: s.startswith("A")), "coding"] = "coding"
    edit_info.loc[
        edit_info.pos.map(lambda s: not s.startswith("A")), "coding"
    ] = "noncoding"
    edit_info.loc[edit_info.pos.map(lambda s: "CONTROL" in s), "group"] = "negctrl"
    edit_info.loc[edit_info.pos.map(lambda s: "CONTROL" in s), "coding"] = "negctrl"
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
