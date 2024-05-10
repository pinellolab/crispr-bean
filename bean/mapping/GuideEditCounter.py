from typing import Tuple, Optional, Sequence, Union
import gzip
import logging
from copy import deepcopy
import os
import sys
import random

import numpy as np
import pandas as pd
from bean import Allele, ReporterScreen
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqPhredIterator
if sys.stderr.isatty():
    # Output into terminal
    from tqdm import tqdm
else:
    # Writing into file
    def tqdm(iterable, **kwargs):
        return iterable

from ._supporting_fn import (
    _base_edit_to_from,
    _get_edited_allele_crispresso,
    _get_fastq_handle,
    _read_count_match,
    _read_is_good_quality,
    _write_alignment_matrix,
    _write_paired_end_reads,
    revcomp,
)

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


class InputFileError(Exception):
    """Raised when the input file is not valid."""


class NoReadsAfterQualityFiltering(Exception):
    """Raised when no reads are remaining after quality filtering."""


strand_str_to_int = {"neg": -1, "pos": 1, "-": -1, "+": 1}
base_revcomp = {"A": "T", "T": "A", "G": "C", "C": "G"}


def _get_stranded_guide_offset(strand: int, start_pos: int, guide_len: int) -> int:
    if strand == -1:
        offset = start_pos + 32 - 6 - 1
    elif strand == 1:
        offset = start_pos - (32 - 6 - guide_len)
    return offset


class GuideEditCounter:
    def __init__(self, **kwargs):
        self.R1_filename = kwargs["R1"]
        self.R2_filename = kwargs["R2"]
        self.target_base_edits = {
            k: _base_edit_to_from(k) for k in kwargs["edited_base"]
        }
        self.min_average_read_quality = kwargs["min_average_read_quality"]
        self.min_single_bp_quality = kwargs["min_single_bp_quality"]

        self.guides_info_df = pd.read_csv(kwargs["sgRNA_filename"])
        self.guides_has_strands = (
            "strand" in self.guides_info_df.columns.tolist()
        )  # @TODO: fix this?
        if self.guides_has_strands:
            info("Considering `strand` column in sgRNA info file")
            if "start_pos" not in self.guides_info_df.columns.tolist():
                raise InputFileError("Guide `start_pos` not in sgRNA info file.")
            if "guide_len" not in self.guides_info_df.columns.tolist():
                self.guides_info_df["guide_len"] = self.guides_info_df.sequence.map(len)
        else:
            info("Ignoring guide strands, all guides are considered positive")

        self.qstart_R1 = kwargs["qstart_R1"]
        self.qend_R1 = kwargs["qend_R1"]
        self.qstart_R2 = kwargs["qstart_R2"]
        self.qend_R2 = kwargs["qend_R2"]
        self.name = kwargs["name"]

        self.count_only_bcmatched = False
        self.sgRNA_filename = kwargs["sgRNA_filename"]
        self._set_sgRNA_df()

        self.database_id = self._get_database_name()
        self.output_dir = os.path.join(
            os.path.abspath(kwargs["output_folder"]),
            f"bean_count_{self.database_id}",
        )
        self._write_start_log()

        self.screen = ReporterScreen(
            X=np.zeros((len(self.guides_info_df), 1)),
            X_edit=np.zeros((len(self.guides_info_df), 1)),
            X_bcmatch=np.zeros((len(self.guides_info_df), 1)),
            guides=self.guides_info_df,
            samples=pd.DataFrame(index=[self.database_id]),
            target_base_changes=",".join(
                [
                    f"{base_edited_from}>{base_edited_to}"
                    for base_edited_from, base_edited_to in self.target_base_edits.items()
                ]
            ),
            tiling=kwargs["tiling"],
        )
        self.screen.guides["guide_len"] = self.screen.guides.sequence.map(len)
        self.screen.uns["reporter_length"] = kwargs["reporter_length"]
        self.screen.uns["reporter_right_flank_length"] = kwargs["reporter_length"] - kwargs["gstart_reporter"] - self.screen.guides["guide_len"].max()
        self.count_guide_edits = kwargs["count_guide_edits"]
        if self.count_guide_edits:
            self.screen.uns["guide_edit_counts"] = {}
        self.count_reporter_edits = (
            kwargs["count_reporter"] or kwargs["count_guide_reporter_alleles"]
        )
        if self.count_reporter_edits:
            self.screen.uns["edit_counts"] = {}
            self.gstart_reporter = kwargs["gstart_reporter"]
            self.screen.uns["allele_counts"] = {}
        self.count_guide_reporter_alleles = kwargs["count_guide_reporter_alleles"]
        if self.count_guide_reporter_alleles:
            self.guide_to_guide_reporter_allele = {}
        self.align_score_threshold = 80
        self.target_pos_col = kwargs["target_pos_col"]

        self.guide_start_seq = kwargs["guide_start_seq"].upper()
        self.guide_end_seq = kwargs["guide_end_seq"].upper()
        self.barcode_start_seq = kwargs["barcode_start_seq"].upper()
        if not self.guide_start_seq == "":
            info(
                f"{self.name}: Using guide_start_seq={self.guide_start_seq} for {self.output_dir}"
            )
        assert (
            self.guide_start_seq == "" or self.guide_end_seq == ""
        ), "Doesn't support both start & end seq matching"
        self.guide_bc = kwargs["guide_bc"]
        if self.guide_bc:
            self.guide_bc_len = kwargs["guide_bc_len"]

        self.offset = kwargs["offset"]
        if self.count_reporter_edits:
            self.reporter_length = kwargs["reporter_length"]
        self.guide_to_allele = {}
        self.n_total_reads = _read_count_match(self.R1_filename, self.R2_filename)

        self.objectify_allele = not kwargs["string_allele"]
        if not self.objectify_allele:
            info(f"{self.name}: Storing allele as strings.")
        self.keep_intermediate = kwargs["keep_intermediate"]
        self.skip_filtering = kwargs["skip_filtering"]
        self.id_number = random.randint(0, int(1e6))
        self.semimatch = 0
        self.bcmatch = 0
        self.nomatch = 0
        self.duplicate_match = 0
        self.duplicate_match_wo_barcode = 0

    def mask_sequence(self, seq, revcomp=False):
        seq_masked = seq
        for base_edited_from, base_edited_to in self.target_base_edits.items():
            if revcomp:
                base_edited_from = base_revcomp[base_edited_from]
                base_edited_to = base_revcomp[base_edited_to]
            seq_masked = seq_masked.replace(base_edited_from, base_edited_to)
        return seq_masked

    def masked_equal(self, seq1, seq2):
        """Tests if two sequences are equal, ignoring the allowed base transition."""
        return self.mask_sequence(seq1) == self.mask_sequence(seq2)

    def _set_sgRNA_df(self):
        """set gRNA info dataframe"""
        with open(self.sgRNA_filename) as infile:
            sgRNA_df = pd.read_csv(infile)
            if "name" not in sgRNA_df.columns or "sequence" not in sgRNA_df.columns:
                raise InputFileError(
                    "Input gRNA info file doesn't have the column 'name' or 'sequence'."
                )
            if self.count_only_bcmatched and "barcode" not in sgRNA_df.columns:
                raise InputFileError(
                    "Input gRNA info file doesn't have the column 'barcode'."
                )
            sgRNA_df = sgRNA_df.set_index("name")
        self.guides_info_df = sgRNA_df
        self.guides_info_df["masked_sequence"] = self.guides_info_df.sequence.map(
            self.mask_sequence
        )
        self.guides_info_df["masked_barcode"] = self.guides_info_df.barcode.map(
            self.mask_sequence
        )
        self.guide_lengths = sgRNA_df.sequence.map(lambda s: len(s)).unique()

    def check_filter_fastq(self):
        """Checks if the quality filtered fastq files already exists,
        and use them if the do."""
        if self.skip_filtering:
            self.filtered_R1_filename = self.R1_filename
            self.filtered_R2_filename = self.R2_filename
            self.n_reads_after_filtering = self.n_total_reads
        else:
            self.filtered_R1_filename = self._jp(
                os.path.basename(self.R1_filename)
                .replace(".fastq", "")
                .replace(".gz", "")
                + "_filtered.fastq.gz"
            )
            self.filtered_R2_filename = self._jp(
                os.path.basename(self.R2_filename)
                .replace(".fastq", "")
                .replace(".gz", "")
                + "_filtered.fastq.gz"
            )
            self._check_names_filter_fastq()
        # if (
        #     path.exists(self.filtered_R1_filename)
        #     and path.exists(self.filtered_R2_filename)
        #     and self.reuse_filtered_reads
        # ):
        #     warn("Using preexisting filtered file")
        # else:
        #     self._check_names_filter_fastq()

    def get_counts(self):
        self.nomatch_R1_filename = self.R1_filename.replace(".fastq", "_nomatch.fastq")
        self.nomatch_R2_filename = self.R2_filename.replace(".fastq", "_nomatch.fastq")
        self.semimatch_R1_filename = self.R1_filename.replace(
            ".fastq", "_semimatch.fastq"
        )
        self.semimatch_R2_filename = self.R2_filename.replace(
            ".fastq", "_semimatch.fastq"
        )
        if self.count_reporter_edits or self.count_guide_edits:
            _write_alignment_matrix(
                self.target_base_edits,
                self.output_dir + "/.aln_mat.txt",
            )

        if self.count_only_bcmatched:
            # count X
            raise NotImplementedError
        else:  # count both bc matched & unmatched guide counts
            self._get_guide_counts_bcmatch_semimatch()

        if self.count_reporter_edits:
            self.screen.uns["edit_counts"] = pd.DataFrame.from_dict(
                self.screen.uns["edit_counts"]
            )

        if self.count_reporter_edits:
            self._write_reporter_alleles()
        if self.count_guide_edits:
            self._write_guide_self_alleles()
        if self.count_guide_reporter_alleles:
            self._write_guide_reporter_alleles()
        count_stat_path = self._jp("mapping_stats.txt")
        count_stat_file = open(count_stat_path, "w")
        count_stat_file.write("Read count with \n")
        count_stat_file.write(
            f"Unique guide match without barcode:\t{self.semimatch}\n"
        )
        count_stat_file.write(f"Unique guide match with barcode:\t{self.bcmatch}\n")
        count_stat_file.write(f"No match:\t{self.nomatch}\n")
        count_stat_file.write(
            f"Duplicate match with barcode:\t{self.duplicate_match}\n"
        )
        count_stat_file.write(
            f"No match with barcode & Duplicate match w/o barcode:\t{self.duplicate_match_wo_barcode}\n"
        )

    def _write_guide_reporter_alleles(self):
        guides = []
        guide_alleles = []
        reporter_alleles = []
        counts = []
        for guide, allele_to_count in self.guide_to_guide_reporter_allele.items():
            if len(allele_to_count.keys()) == 0:
                continue
            guides.extend([guide] * len(allele_to_count.keys()))
            for key, count in allele_to_count.items():
                reporter_allele, guide_allele = key
                guide_alleles.append(guide_allele)
                reporter_alleles.append(reporter_allele)
                counts.append(count)
        if not (
            len(guides) == len(reporter_alleles) == len(guide_alleles) == len(counts)
        ):
            raise ValueError(
                f"Guides:{len(guides)}, guide_alleles:{len(guide_alleles)}, reporter_alleles: {len(reporter_alleles)}, counts:{len(counts)}"
            )
        self.screen.uns["guide_reporter_allele_counts"] = pd.DataFrame(
            {
                "guide": guides,
                "guide_allele": guide_alleles,
                "reporter_allele": reporter_alleles,
                self.database_id: counts,
            }
        )

        if "guide" in self.screen.uns["guide_reporter_allele_counts"].columns:
            self.screen.uns["guide_reporter_allele_counts"].guide = (
                self.screen.guides.index[
                    self.screen.uns["guide_reporter_allele_counts"].guide
                ]
            )

    def _write_reporter_alleles(self):
        guides = []
        alleles = []
        counts = []
        for guide, allele_to_count in self.guide_to_allele.items():
            if len(allele_to_count.keys()) == 0:
                continue
            guides.extend([guide] * len(allele_to_count.keys()))
            for allele, count in allele_to_count.items():
                alleles.append(allele)
                counts.append(count)
        if not (len(guides) == len(alleles) == len(counts)):
            raise ValueError(
                f"Guides:{len(guides)}, alleles:{len(alleles)}, counts:{len(counts)}"
            )
        self.screen.uns["allele_counts"] = pd.DataFrame(
            {"guide": guides, "allele": alleles, self.database_id: counts}
        )

        if "guide" in self.screen.uns["allele_counts"].columns:
            self.screen.uns["allele_counts"].guide = self.screen.guides.index[
                self.screen.uns["allele_counts"].guide
            ]

    def _write_guide_self_alleles(self):
        guides = []
        alleles = []
        counts = []
        guide_edits = deepcopy(self.screen.uns["guide_edit_counts"])
        for guide, allele_to_count in guide_edits.items():
            if len(allele_to_count.keys()) == 0:
                continue
            guides.extend([guide] * len(allele_to_count.keys()))
            for allele, count in allele_to_count.items():
                alleles.append(allele)
                counts.append(count)
        if not (len(guides) == len(alleles) == len(counts)):
            raise ValueError(
                f"Guides:{len(guides)}, alleles:{len(alleles)}, counts:{len(counts)}"
            )
        self.screen.uns["guide_edit_counts"] = pd.DataFrame(
            {"guide": guides, "allele": alleles, self.database_id: counts}
        )

        if "guide" in self.screen.uns["guide_edit_counts"].columns:
            self.screen.uns["guide_edit_counts"].guide = self.screen.guides.index[
                self.screen.uns["guide_edit_counts"].guide
            ]

    def _count_guide_edits(
        self, matched_guide_idx, R1_record: SeqIO.SeqRecord, single_base_qual_cutoff=30
    ):
        if self.guides_has_strands:
            strand = self.screen.guides.strand.iloc[matched_guide_idx]
            guide_strand = strand_str_to_int.get(strand, 1)
        else:
            guide_strand = 1
        ref_guide_seq = self.screen.guides.sequence.iloc[matched_guide_idx]
        if "chrom" in self.screen.guides.columns.tolist():
            chrom = self.screen.guides.chrom.iloc[matched_guide_idx]
        else:
            if self.screen.tiling:
                raise ValueError(
                    "'chrom' column denoting gRNA chromosome location column not present in the input. Please check the input file."
                )
            else:
                chrom = None
        read_guide_seq, read_guide_qual = self.get_guide_seq_qual(
            R1_record, len(ref_guide_seq)
        )
        guide_edit_allele, score = _get_edited_allele_crispresso(
            ref_seq=ref_guide_seq,
            query_seq=read_guide_seq,
            target_base_edits=self.target_base_edits,
            aln_mat_path=self.output_dir + "/.aln_mat.txt",
            offset=0,
            strand=guide_strand,
            chrom=chrom,
            start_pos=0,
            end_pos=len(ref_guide_seq),
            positionwise_quality=np.array(read_guide_qual),
            quality_thres=single_base_qual_cutoff,
            objectify_allele=self.objectify_allele,
        )

        if not self.count_reporter_edits:
            self._update_counted_guide_allele(matched_guide_idx, guide_edit_allele)
        return (guide_edit_allele, score)

    def _get_strand_offset_from_guide_index(self, guide_idx: int) -> Tuple[int, int]:
        """Returns guide starnd and offset for a given guide index."""
        if self.guides_has_strands:
            strand = self.screen.guides.strand.iloc[guide_idx]
            if strand in strand_str_to_int:
                guide_strand = strand_str_to_int[strand]
                offset = _get_stranded_guide_offset(
                    strand=guide_strand,
                    start_pos=self.screen.guides.iloc[guide_idx].start_pos,
                    guide_len=self.screen.guides.iloc[guide_idx].guide_len,
                )
            else:
                guide_strand = 1
                offset = 0

        else:
            guide_strand = 1
            if self.target_pos_col in self.screen.guides.columns:
                offset = -(self.screen.guides.iloc[guide_idx][self.target_pos_col] - 1)
            else:
                offset = 0
        return (guide_strand, offset)

    def _update_counted_allele(
        self, guide_idx: int, allele: Union[Allele, str]
    ) -> None:
        """Add allele count to self.guide_to_allele dictionary."""
        if guide_idx in self.guide_to_allele.keys():
            if allele in self.guide_to_allele[guide_idx].keys():
                self.guide_to_allele[guide_idx][allele] += 1
            else:
                self.guide_to_allele[guide_idx][allele] = 1
        else:
            self.guide_to_allele[guide_idx] = {allele: 1}

    def _update_counted_guide_allele(
        self, guide_idx: int, allele: Union[Allele, str]
    ) -> None:
        """Add allele count to self.guide_to_allele dictionary."""
        if guide_idx in self.screen.uns["guide_edit_counts"].keys():
            if allele in self.screen.uns["guide_edit_counts"][guide_idx].keys():
                self.screen.uns["guide_edit_counts"][guide_idx][allele] += 1
            else:
                self.screen.uns["guide_edit_counts"][guide_idx][allele] = 1
        else:
            self.screen.uns["guide_edit_counts"][guide_idx] = {allele: 1}

    def _update_counted_allele_and_guideAllele(
        self,
        guide_idx: int,
        allele: Union[Allele, str],
        guide_allele: Union[str, Allele],
    ) -> None:
        """Add count of (guide allele, reporter allele) combination to self.guide_reporter_allele dictionary."""
        if guide_idx in self.guide_to_guide_reporter_allele.keys():
            if (allele, guide_allele) in self.guide_to_guide_reporter_allele[
                guide_idx
            ].keys():
                self.guide_to_guide_reporter_allele[guide_idx][
                    (allele, guide_allele)
                ] += 1
            else:
                self.guide_to_guide_reporter_allele[guide_idx][
                    (allele, guide_allele)
                ] = 1
        else:
            self.guide_to_guide_reporter_allele[guide_idx] = {(allele, guide_allele): 1}

    def _count_reporter_edits(
        self,
        matched_guide_idx: int,
        R1_seq: str,
        R2_record: SeqIO.SeqRecord,
        R2_start: int = 0,
        single_base_qual_cutoff: int = 30,
        guide_allele: Optional[Union[str, Allele]] = None,
    ):
        """
        Count edits in a single read to save as allele.

        Args
        --
        matched_guide_idx: index of guides in self.screen.guides to get information from
        R1_seq: Read1 sequence
        R2_record: Read2 sequence record with quality
        single_base_qual_cutoff: Ignore this base if the Phread quality score is less than this threshold
        guide_allele: Allele from baseedit in gRNA spacer sequence when paired with guide allele.
        """
        ref_reporter_seq = self.screen.guides.reporter.iloc[matched_guide_idx]
        read_reporter_seq, read_reporter_qual = self.get_reporter_seq_qual(
            R2_record, R2_start
        )

        guide_strand, offset = self._get_strand_offset_from_guide_index(
            matched_guide_idx
        )
        if "chrom" in self.screen.guides.columns.tolist():
            chrom = self.screen.guides.chrom.iloc[matched_guide_idx]
        else:
            chrom = None
        allele, score = _get_edited_allele_crispresso(
            ref_seq=ref_reporter_seq,
            query_seq=read_reporter_seq,
            target_base_edits=self.target_base_edits,
            aln_mat_path=self.output_dir + "/.aln_mat.txt",
            offset=offset,
            strand=guide_strand,
            chrom=chrom,
            positionwise_quality=np.array(read_reporter_qual),
            quality_thres=single_base_qual_cutoff,
            objectify_allele=self.objectify_allele,
        )

        if score < self.align_score_threshold:
            self.semimatch += 1
            self.bcmatch -= 1
            return

        if self.count_reporter_edits:
            self._update_counted_allele(matched_guide_idx, allele)

        if self.count_guide_reporter_alleles and (guide_allele is not None):
            self._update_counted_allele_and_guideAllele(
                matched_guide_idx, allele, guide_allele
            )

    def _get_guide_counts_bcmatch_semimatch(
        self, bcmatch_layer="X_bcmatch", semimatch_layer="X_semimatch"
    ):
        self.screen.layers[semimatch_layer] = np.zeros_like(self.screen.X)
        R1_iter, R2_iter = self._get_fastq_iterators(
            self.filtered_R1_filename, self.filtered_R2_filename
        )
        if self.keep_intermediate:
            outfile_R1_nomatch, outfile_R2_nomatch = self._get_fastq_handle("nomatch")
            outfile_R1_semimatch, outfile_R2_semimatch = self._get_fastq_handle(
                "semimatch"
            )
            outfile_R1_dup_wo_bc, outfile_R2_dup_wo_bc = self._get_fastq_handle(
                "duplicate_wo_barcode"
            )
            outfile_R1_dup, outfile_R2_dup = self._get_fastq_handle("duplicate")
        tqdm_reads= tqdm(
            enumerate(zip(R1_iter, R2_iter)),
            total=self.n_reads_after_filtering,
            postfix=f"n_read={self.bcmatch}",
        )
        for i, (r1, r2) in tqdm_reads:
            R1_seq = str(r1.seq)
            R2_seq = str(r2.seq)

            (
                bc_match,
                semimatch,
                R2_start,
            ) = self._match_read_to_sgRNA_bcmatch_semimatch(R1_seq, R2_seq)
            if len(bc_match) == 0:
                if len(semimatch) == 0:  # no guide match
                    if self.keep_intermediate:
                        _write_paired_end_reads(
                            r1, r2, outfile_R1_nomatch, outfile_R2_nomatch
                        )
                    self.nomatch += 1
                elif len(semimatch) >= 2:  # Duplicate match if w/o barcode
                    if self.keep_intermediate:
                        _write_paired_end_reads(
                            r1, r2, outfile_R1_dup_wo_bc, outfile_R2_dup_wo_bc
                        )
                    self.duplicate_match_wo_barcode += 1
                else:  # guide match with no barcode match
                    matched_guide_idx = semimatch[0]
                    self.screen.layers[semimatch_layer][matched_guide_idx, 0] += 1
                    if self.count_guide_edits:
                        guide_allele, _ = self._count_guide_edits(
                            matched_guide_idx, r1
                        )
                    self.semimatch += 1

            elif len(bc_match) >= 2:  # duplicate mapping
                if self.keep_intermediate:
                    _write_paired_end_reads(r1, r2, outfile_R1_dup, outfile_R2_dup)
                self.duplicate_match += 1

            else:  # unique barcode match
                matched_guide_idx = bc_match[0]
                self.screen.layers[bcmatch_layer][matched_guide_idx, 0] += 1
                self.bcmatch += 1
                if self.count_guide_edits or self.count_guide_reporter_alleles:
                    guide_allele, _ = self._count_guide_edits(matched_guide_idx, r1)
                if self.count_reporter_edits:
                    # TODO: what if reporter seq doesn't match barcode & guide?
                    if self.count_guide_reporter_alleles:
                        self._count_reporter_edits(
                            matched_guide_idx,
                            R1_seq,
                            r2,
                            R2_start=R2_start,
                            guide_allele=guide_allele,
                        )
                    else:
                        self._count_reporter_edits(
                            matched_guide_idx,
                            R1_seq,
                            r2,
                            R2_start=R2_start,
                        )
            if sys.stderr.isatty():
                tqdm_reads.postfix = f"n_read={self.bcmatch}"
                tqdm_reads.update()

        self.screen.X = (
            self.screen.layers[semimatch_layer] + self.screen.layers[bcmatch_layer]
        )

    def _write_allele(self, guide_idx: int, allele: Allele, uns_key="allele_counts"):
        if (guide_idx, str(allele)) in self.screen.uns[uns_key].keys():
            self.screen.uns[uns_key][(guide_idx, str(allele))] += 1
        else:
            self.screen.uns[uns_key][(guide_idx, str(allele))] = 1

    def _write_edits(self, guide_idx: int, allele: Allele, uns_key="edit_counts"):
        for edit in allele.edits:
            if (guide_idx, str(edit)) in self.screen.uns[uns_key].keys():
                self.screen.uns[uns_key][(guide_idx, str(edit))] += 1
            else:
                self.screen.uns[uns_key][(guide_idx, str(edit))] = 1

    def _write_guide_reporter_allele(
        self,
        guide_idx: int,
        allele: Allele,
        guide_allele: Allele,
        uns_key="guide_reporter_allele_counts",
    ):
        if len(allele.edits) == 0 and len(guide_allele.edits) == 0:
            return
        if (guide_idx, str(allele), str(guide_allele)) in self.screen.uns[
            uns_key
        ].keys():
            self.screen.uns[uns_key][(guide_idx, str(allele), str(guide_allele))] += 1
        else:
            self.screen.uns[uns_key][(guide_idx, str(allele), str(guide_allele))] = 1

    def get_guide_seq(self, R1_seq, R2_seq, guide_length):
        """This can be edited by user based on the read construct."""
        if self.guide_end_seq == "":
            guide_start_idx = self.mask_sequence(R1_seq).find(
                self.mask_sequence(self.guide_start_seq)
            )
            if guide_start_idx == -1:
                return None
            if guide_start_idx + guide_length >= len(R1_seq):
                return None
            guide_start_idx = guide_start_idx + len(self.guide_start_seq)
            gRNA_seq = R1_seq[guide_start_idx : (guide_start_idx + guide_length)]
        else:
            guide_end_idx = self.mask_sequence(R1_seq).find(
                self.mask_sequence(self.guide_end_seq)
            )
            if guide_end_idx == -1:
                return None
            if guide_end_idx - guide_length < 0:
                return None
            gRNA_seq = R1_seq[(guide_end_idx - guide_length) : guide_end_idx]

        return None if len(gRNA_seq) != guide_length else gRNA_seq

    def get_guide_seq_qual(
        self, R1_record: SeqIO.SeqRecord, guide_length: int
    ) -> Tuple[str, Sequence]:
        R1_seq = R1_record.seq
        guide_start_idx = self.mask_sequence(R1_seq).find(
            self.mask_sequence(self.guide_start_seq)
        )
        guide_start_idx = guide_start_idx + len(self.guide_start_seq)
        seq = R1_record[guide_start_idx : (guide_start_idx + guide_length)]
        return (str(seq.seq), seq.letter_annotations["phred_quality"])

    def get_reporter_seq(self, R1_seq, R2_seq):
        """This can be edited by user based on the read construct."""
        return revcomp(
            R2_seq[self.guide_bc_len : (self.guide_bc_len + self.reporter_length)]
        )

    def get_reporter_seq_qual(self, R2_record: SeqIO.SeqRecord, R2_start=0):
        seq = R2_record[
            (R2_start + self.guide_bc_len) : (
                R2_start + self.guide_bc_len + self.reporter_length
            )
        ].reverse_complement()
        return (str(seq.seq), seq.letter_annotations["phred_quality"])

    def get_barcode(self, R1_seq, R2_seq):
        if self.barcode_start_seq != "":
            barcode_start_idx = self.mask_sequence(R2_seq, revcomp=True).find(
                self.mask_sequence(revcomp(self.barcode_start_seq), revcomp=True)
            )
            if barcode_start_idx == -1:
                return -1, ""
            barcode_start_idx += len(self.barcode_start_seq)
        else:
            barcode_start_idx = 0
        return barcode_start_idx, revcomp(
            R2_seq[barcode_start_idx : (barcode_start_idx + self.guide_bc_len)]
        )

    def _match_read_to_sgRNA_bcmatch_semimatch(self, R1_seq: str, R2_seq: str):
        # This should be adjusted for each experimental recipes.
        bc_start_idx, guide_barcode = self.get_barcode(R1_seq, R2_seq)
        bc_match_idx = np.array([])
        semimatch_idx = np.array([])
        for guide_length in self.guide_lengths:
            seq = self.get_guide_seq(R1_seq, R2_seq, guide_length)
            if seq is None:
                continue

            _seq_match = np.where(
                self.mask_sequence(seq) == self.screen.guides.masked_sequence
            )[0]
            _bc_match = np.where(
                self.mask_sequence(guide_barcode) == self.screen.guides.masked_barcode
            )[0]

            bc_match_idx = np.append(
                bc_match_idx, np.intersect1d(_seq_match, _bc_match)
            )
            semimatch_idx = np.append(semimatch_idx, _seq_match)

        return bc_match_idx.astype(int), semimatch_idx.astype(int), bc_start_idx

    def _get_guide_position_seq_of_read(self, seq):
        guide_start_idx = self._get_guide_start_idx(seq)
        if guide_start_idx == -1:
            return None

        return [
            seq[guide_start_idx : (guide_start_idx + guide_length)]
            for guide_length in self.guide_lengths
        ]

    def _get_guide_start_idx(self, seq):
        start_seq_idx = seq.find(self.guide_start_seq)
        return -1 if start_seq_idx == -1 else start_seq_idx + len(self.guide_start_seq)

    def get_gRNA_barcode(self, R1_seq, R2_seq):
        # This can be adjusted for different construct design.
        return revcomp(R2_seq[: self.guide_bc_len])

    def _get_fastq_handle(
        self,
        out_type: Optional[str] = None,
    ):
        assert out_type in {
            "semimatch",
            "nomatch",
            "duplicate_wo_barcode",
            "duplicate",
        }
        R1_filename = self.R1_filename.replace(".fastq", f"_{out_type}.fastq")
        R2_filename = self.R2_filename.replace(".fastq", f"_{out_type}.fastq")
        R1_handle = _get_fastq_handle(R1_filename, "w")
        R2_handle = _get_fastq_handle(R2_filename, "w")

        return (R1_handle, R2_handle)

    def _get_fastq_iterators(self, R1_filename=None, R2_filename=None):
        if R1_filename is None:
            R1_filename = self.R1_filename
        if R2_filename is None:
            R2_filename = self.R2_filename
        R1_handle = _get_fastq_handle(R1_filename)
        R2_handle = _get_fastq_handle(R2_filename)

        R1_iterator = FastqPhredIterator(R1_handle)
        R2_iterator = FastqPhredIterator(R2_handle)

        return (R1_iterator, R2_iterator)

    def _get_seq_handles(self):
        R1_handle = _get_fastq_handle(self.R1_filename)
        R2_handle = _get_fastq_handle(self.R2_filename)
        return (R1_handle, R2_handle)

    def _check_names_filter_fastq(self, filter_by_qual=False):
        if self.min_average_read_quality > 0 or self.min_single_bp_quality > 0:
            info(
                f"Filtering reads with average bp quality < {self.min_average_read_quality} and single bp quality < {self.min_single_bp_quality} ..."
            )

        if self.qend_R1 > 0 or self.qend_R2 > 0:
            info(
                f"In the filtering, bases up to position {self.qend_R1} of R1 and {self.qend_R2} of R2 are considered."
            )

        R1_iter, R2_iter = self._get_fastq_iterators()
        (self.n_reads_after_filtering) = self._check_readname_match_and_filter_quality(
            R1_iter, R2_iter, filter_by_qual
        )

        if self.n_reads_after_filtering == 0:
            raise NoReadsAfterQualityFiltering(
                "No reads in input or no reads survived the average or single bp quality filtering."
            )
        else:
            info(
                "Number of reads in input:%d\tNumber of reads after filtering:%d\n"
                % (self.n_total_reads, self.n_reads_after_filtering)
            )

    def _check_readname_match_and_filter_quality(
        self, R1_iter, R2_iter, filter_by_qual=False
    ) -> int:
        R1_filtered = gzip.open(self.filtered_R1_filename, "wt+")
        R2_filtered = gzip.open(self.filtered_R2_filename, "wt+")

        n_reads_after_filtering = 0
        for R1_record, R2_record in tqdm(
            zip(R1_iter, R2_iter), total=self.n_total_reads
        ):
            if R1_record.name != R2_record.name:
                raise InputFileError(
                    "R1 and R2 read discordance in read {} and {}".format(
                        R1_record.name, R2_record.name
                    )
                )
            if filter_by_qual:
                R1_quality_pass = _read_is_good_quality(
                    R1_record,
                    self.min_average_read_quality,
                    self.min_single_bp_quality,
                    self.qend_R1,
                )
                R2_quality_pass = _read_is_good_quality(
                    R2_record,
                    self.min_average_read_quality,
                    self.min_single_bp_quality,
                    self.qend_R2,
                )
            else:
                R1_quality_pass = True
                R2_quality_pass = True
            if R1_quality_pass and R2_quality_pass:
                n_reads_after_filtering += 1
                R1_filtered.write(R1_record.format("fastq"))
                R2_filtered.write(R2_record.format("fastq"))
        R1_filtered.close()
        R2_filtered.close()
        return n_reads_after_filtering

    def _write_start_log(self):
        try:
            os.makedirs(self.output_dir)
            info(f"Creating Folder {self.output_dir}")
        except OSError:
            info(f"Folder {self.output_dir} already exists.")
        self.log_filename = self._jp("beanCount_RUNNING_LOG.txt")
        logging.getLogger().addHandler(logging.FileHandler(self.log_filename))
        with open(self.log_filename, "w+") as outfile:
            outfile.write(
                "[Command used]:\nCRISPRessoCount %s\n\n[Execution log]:\n"
                % " ".join(sys.argv)
            )

    def _jp(self, filename):
        return os.path.join(self.output_dir, filename)

    def _get_database_name(self):
        get_name_from_fasta = (
            lambda x: os.path.basename(x)
            .replace(".fastq", "")
            .replace(".gz", "")
            .replace("_R1", "")
        )

        return self.name or f"{get_name_from_fasta(self.R1_filename)}"
