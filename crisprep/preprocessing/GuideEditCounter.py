from typing import Set, Literal, Iterable
import os
from os import path
import sys
import logging
from tqdm import tqdm
import pdb
import gzip
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from crisprep import ReporterScreen, Edit, Allele
#from crisprep.framework.Edit import Allele

from ._supporting_fn import (
    _base_edit_to_from,
    _read_count_match,
    _get_fastq_handle,
    _write_paired_end_reads,
    _get_edited_allele,
    _check_readname_match,
    _read_is_good_quality,
    revcomp
)

logging.basicConfig(level=logging.INFO,
                     format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                     datefmt='%a, %d %b %Y %H:%M:%S',
                     stream=sys.stderr,
                     filemode="w"
                     )
error   = logging.critical
warn    = logging.warning
debug   = logging.debug
info    = logging.info



class NTException(Exception):
    pass


class InputFileError(Exception):
    pass


class NoReadsAfterQualityFiltering(Exception):
    pass


class GuideEditCounter:
    def __init__(self, **kwargs):
        self.R1_filename = kwargs["R1"]
        self.R2_filename = kwargs["R2"]
        self.base_edited_from = kwargs["edited_base"]
        self.base_edited_to = _base_edit_to_from(self.base_edited_from)
        self.min_average_read_quality = kwargs["min_average_read_quality"]
        self.min_single_bp_quality = kwargs["min_single_bp_quality"]

        self.guides_info_df = pd.read_csv(kwargs["sgRNA_filename"])
        self.guides_has_strands = ("Strand" in self.guides_info_df.columns) or \
        ("strand" in self.guides_info_df.columns)
        if self.guides_has_strands:
            info("Considering strand information of guides")
            assert "gRNA_start" in self.guides_info_df.columns
            assert "gRNA_end" in self.guides_info_df.columns
        else:
            info("Ignoring guide strands, all guides are considered positive")

        self.qstart_R1 = kwargs["qstart_R1"]
        self.qend_R1 = kwargs["qend_R1"]
        self.qstart_R2 = kwargs["qstart_R2"]
        self.qend_R2 = kwargs["qend_R2"]
        self.name = kwargs["name"]

        self.sgRNA_filename = kwargs["sgRNA_filename"]
        self.count_only_bcmatched = False
        self._set_sgRNA_df()

        self.database_id = self._get_database_name()
        self.output_dir = os.path.join(
            os.path.abspath(kwargs["output_folder"]),
            "CRISPRessoCount_on_%s" % self.database_id,
        )
        self._write_start_log()

        self.screen = ReporterScreen(
            X=np.zeros((len(self.guides_info_df), 1)),
            X_edit=np.zeros((len(self.guides_info_df), 1)),
            X_bcmatch=np.zeros((len(self.guides_info_df), 1)),
            guides=self.guides_info_df,
            condit=pd.DataFrame(index=[self.database_id]),
        )

        self.count_reporter_edits = kwargs["count_reporter"]
        if self.count_reporter_edits:
            self.screen.uns["edit_counts"] = dict()
            # pd.DataFrame(
            #     columns = ["guide_name", "edit_id", "count"]).set_index(
            #     ["guide_name", "edit_id"], drop = True)
            # #self.screen.layers["edits"] = np.zeros_like(self.screen.X)
            self.gstart_reporter = kwargs["gstart_reporter"]

        self.count_edited_alleles = kwargs["count_allele"]
        if self.count_edited_alleles:
            self.screen.uns["allele_counts"] = dict()
            # pd.DataFrame(
            #     columns = ["guide_name", "allele_id", "count"]).set_index(
            #     ["guide_name", "allele_id"], drop = True)
            self.screen.uns["alleles"] = pd.DataFrame(columns = ["allele_id", "edit_id"]) 
            # Dict[guide_name -> List[edit_ID]]
            
        self.guide_start_seq = kwargs["guide_start_seq"]
        self.guide_bc = kwargs["guide_bc"]
        if self.guide_bc:
            self.guide_bc_len = kwargs["guide_bc_len"]

        self.offset = kwargs["offset"]
        if self.count_reporter_edits:
            self.reporter_length = kwargs["reporter_length"]

        self.n_total_reads = _read_count_match(self.R1_filename, self.R2_filename)
        

        self.keep_intermediate = kwargs["keep_intermediate"]
        self.nomatch = 0
        self.duplicate_match = 0
        self.duplicate_match_wo_barcode = 0
        

    def _set_sgRNA_df(self):
        with open(self.sgRNA_filename) as infile:
            sgRNA_df = pd.read_csv(infile)
            if not ("name" in sgRNA_df.columns and "sequence" in sgRNA_df.columns):
                raise InputFileError(
                    "Input gRNA info file doesn't have the column 'name' or 'sequence'."
                )
            if self.count_only_bcmatched and "barcode" not in sgRNA_df.columns:
                raise InputFileError(
                    "Input gRNA info file doesn't have the column 'barcode'."
                )
            sgRNA_df = sgRNA_df.set_index("name")
        self.guides_info_df = sgRNA_df
        self.guide_lengths = sgRNA_df.sequence.map(lambda s: len(s)).unique()

    def check_filter_fastq(self):
        self.filtered_R1_filename = self._jp(
            os.path.basename(self.R1_filename).replace(".fastq", "").replace(".gz", "")
            + "_filtered.fastq.gz"
        )
        self.filtered_R2_filename = self._jp(
            os.path.basename(self.R2_filename).replace(".fastq", "").replace(".gz", "")
            + "_filtered.fastq.gz"
        )
        if path.exists(self.filtered_R1_filename) and path.exists(
            self.filtered_R2_filename
        ):
            warn("Using preexisting filtered file")
        else:
            self._check_names_filter_fastq()

    def get_counts(self):
        infile_R1 = _get_fastq_handle(self.R1_filename)
        infile_R2 = _get_fastq_handle(self.R2_filename)

        self.nomatch_R1_filename = self.R1_filename.replace(".fastq", "_nomatch.fastq")
        self.nomatch_R2_filename = self.R2_filename.replace(".fastq", "_nomatch.fastq")
        self.semimatch_R1_filename = self.R1_filename.replace(
            ".fastq", "_semimatch.fastq"
        )
        self.semimatch_R2_filename = self.R2_filename.replace(
            ".fastq", "_semimatch.fastq"
        )

        if self.count_only_bcmatched:
            # count X
            self._get_guide_counts_bcmatch()
        else:  # count both bc matched & unmatched guide counts
            self._get_guide_counts_bcmatch_semimatch()

        if self.count_reporter_edits:
            mi = pd.MultiIndex.from_tuples(self.screen.uns["edit_counts"].keys(), names = ["guide", "edit"])
            self.screen.uns["edit_counts"] = pd.DataFrame.from_dict(
                self.screen.uns["edit_counts"], orient = "index", columns = [self.database_id])
            self.screen.uns["edit_counts"].index = mi
            self.screen.uns["edit_counts"].reset_index(inplace = True)
            self.screen.uns["edit_counts"].rename(columns = {"level_0":"guide", "level_1":"edit"}, inplace= True)
            self.screen.uns["edit_counts"].guide = self.screen.guides.index[self.screen.uns["edit_counts"].guide.to_numpy(dtype = int)]
            #self.screen.uns["edit_counts"].set_index(["guide", "edit"], drop = True, inplace = True)

        if self.count_edited_alleles:
            mi = pd.MultiIndex.from_tuples(self.screen.uns["allele_counts"].keys(), names = ["guide", "allele"])
            self.screen.uns["allele_counts"] = pd.DataFrame.from_dict(self.screen.uns["allele_counts"], orient = "index", columns = [self.database_id])
            self.screen.uns["allele_counts"].index = mi
            self.screen.uns["allele_counts"].reset_index(inplace = True)
            self.screen.uns["allele_counts"].rename(columns = {"level_0":"guide", "level_1":"allele"}, inplace = True)
            self.screen.uns["allele_counts"].guide = self.screen.guides.index[self.screen.uns["allele_counts"].guide.to_numpy(dtype = int)]
            #self.screen.uns["allele_counts"].set_index(["guide", "allele"], drop = True, inplace = True)

        info("Read count with \nno match:\t{}\nduplicate match:\t{}\nduplicate match wo barcode:\t{}\n".format(
            self.nomatch, self.duplicate_match, self.duplicate_match_wo_barcode))

    def _gRNA_eq(self, guide:str, observed:str):
        if len(observed) != len(guide):
            return False
        for observed_nt, guide_nt in zip(observed, guide):
            if (observed_nt == guide_nt) or (
                (guide_nt, observed_nt) == (self.base_edited_from, self.base_edited_to)
            ):
                continue
            else:
                return False
        return True

    def _get_guide_counts_bcmatch(self):
        NotImplemented

    def _get_guide_counts_bcmatch_semimatch(self, bcmatch_layer = "X_bcmatch", semimatch_layer = "X"):
        strand_str_to_int = {"neg":-1, "pos":1}

        self.screen.layers[semimatch_layer] = np.zeros_like((self.screen.X))
        R1_iter, R2_iter = self._get_fastq_iterators()

        outfile_R1_nomatch, outfile_R2_nomatch = self._get_fastq_handle("nomatch")
        outfile_R1_semimatch, outfile_R2_semimatch = self._get_fastq_handle("semimatch")
        outfile_R1_dup_wo_bc, outfile_R2_dup_wo_bc = self._get_fastq_handle(
            "duplicate_wo_barcode"
        )
        outfile_R1_dup, outfile_R2_dup = self._get_fastq_handle("duplicate")

        for i, (r1, r2) in tqdm(
            enumerate(zip(R1_iter, R2_iter)), total=self.n_reads_after_filtering
        ):
            _, R1_seq, _ = r1
            _, R2_seq, _ = r2

            bc_match, semimatch = self._match_read_to_sgRNA_bcmatch_semimatch(
                R1_seq, R2_seq
            )

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
                    try:
                        self.screen.layers[semimatch_layer][semimatch[0], 0] += 1
                    except:
                        print(semimatch)

            elif len(bc_match) >= 2:
                if self.keep_intermediate: _write_paired_end_reads(r1, r2, outfile_R1_dup, outfile_R2_dup)
                self.duplicate_match += 1
            else:  # unique barcode match
                matched_guide_idx = bc_match[0]
                self.screen.layers[bcmatch_layer][matched_guide_idx, 0] += 1

                if self.count_reporter_edits:
                    # TBD: what if reporter seq doesn't match barcode & guide?
                    ref_reporter_seq = self.screen.guides.Reporter[matched_guide_idx]
                    read_reporter_seq = self.get_reporter_seq(R1_seq, R2_seq)

                    if self.guides_has_strands:
                        try:
                            guide_strand = strand_str_to_int[self.screen.guides.Strand[matched_guide_idx]]
                            if guide_strand == -1: 
                                offset = self.screen.guides.start_pos[matched_guide_idx] + \
                                    self.screen.guides.guide_len[matched_guide_idx] - 1
                            if guide_strand == 1: 
                                offset = self.screen.guides.start_pos[matched_guide_idx]
                        except KeyError: #control guides
                            guide_strand = 1
                            offset = 0
                    else:
                        guide_strand = 1
                        offset = self.screen.guides.offset[matched_guide_idx]
                    
                    # Set relative start/end position in reporter to count edits
                    start_pos = len(read_reporter_seq) - self.gstart_reporter - len(self.screen.guides.sequence[matched_guide_idx])
                    end_pos = len(read_reporter_seq) - self.gstart_reporter
                    assert (start_pos == 7 and len(self.screen.guides.sequence[matched_guide_idx])==19) or \
                        (start_pos == 6 and len(self.screen.guides.sequence[matched_guide_idx])==20)
                    allele = _get_edited_allele(
                        ref_seq = ref_reporter_seq,
                        query_seq = read_reporter_seq,
                        offset = offset,
                        strand = guide_strand,
                        start_pos = start_pos,
                        end_pos = end_pos
                    )

                    if self.count_edited_alleles:
                        self._write_allele(matched_guide_idx, allele)
                    self._write_edits(matched_guide_idx, allele)

        self.screen.X = self.screen.layers[semimatch_layer] + self.screen.layers[bcmatch_layer]


    def _write_allele(self, guide_idx: int, allele: Allele):
        if len(allele.edits) == 0: return
        if (guide_idx, str(allele)) in self.screen.uns["allele_counts"].keys():
            self.screen.uns["allele_counts"][(guide_idx, str(allele))] += 1
        else:
            self.screen.uns["allele_counts"][(guide_idx, str(allele))] = 1


    def _write_edits(self, guide_idx: int, allele: Allele):
        for edit in allele.edits:
            if (guide_idx, str(edit)) in self.screen.uns["edit_counts"].keys():
                self.screen.uns["edit_counts"][(guide_idx, str(edit))] += 1
            else: 
                self.screen.uns["edit_counts"][(guide_idx, str(edit))] = 1


    def get_guide_seq(self, R1_seq, R2_seq, guide_length):
        """This can be edited by user based on the read construct."""
        guide_start_idx = R1_seq.find(self.guide_start_seq)
        if guide_start_idx == -1:
            return None
        guide_start_idx = guide_start_idx + len(self.guide_start_seq)
        gRNA_seq = R1_seq[guide_start_idx : guide_start_idx + guide_length]
        return gRNA_seq

    def get_reporter_seq(self, R1_seq, R2_seq):
        """This can be edited by user based on the read construct."""
        reporter_seq = revcomp(
            R2_seq[self.guide_bc_len:(self.guide_bc_len + self.reporter_length)]
        )
        return reporter_seq

    def get_barcode(self, R1_seq, R2_seq):
        """This can be edited by user based on the read construct."""
        return(revcomp(R2_seq[:self.guide_bc_len]))


    def _match_read_to_sgRNA_bcmatch_semimatch(self, R1_seq: str, R2_seq: str):
        # This should be adjusted for each experimental recipes.'
        guide_barcode = self.get_barcode(R1_seq, R2_seq)
        bc_match_idx = np.array([])
        semimatch_idx = np.array([])
        for guide_length in self.guide_lengths:
            seq = self.get_guide_seq(R1_seq, R2_seq, guide_length)
            if seq is None: continue

            _seq_match = np.array(list(map(lambda x: self._gRNA_eq(x, seq), self.guides_info_df.sequence)))
            assert len(_seq_match) == len(self.guides_info_df.sequence)
            _bc_match = np.array(list(map(lambda x: x == guide_barcode, self.guides_info_df.barcode)))
            assert len(_bc_match) == len(self.guides_info_df.barcode)
            bc_match_idx = np.append(bc_match_idx, np.where(_seq_match & _bc_match)[0])
            semimatch_idx = np.append(semimatch_idx, np.where(_seq_match & np.invert(_bc_match))[0])

        return (bc_match_idx.astype(int), semimatch_idx.astype(int))

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
        if start_seq_idx == -1:
            return -1
        return start_seq_idx + len(self.guide_start_seq)

    def get_gRNA_barcode(self, R1_seq, R2_seq):
        # This can be adjusted for different construct design.
        return revcomp(R2_seq[:self.guide_bc_len])

    def _get_fastq_handle(
        self,
        out_type: Literal[
            "semimatch", "nomatch", "duplicate_wo_barcode", "duplicate"
        ] = None,
    ):
        R1_filename = self.R1_filename.replace(".fastq", "_{}.fastq".format(out_type))
        R2_filename = self.R2_filename.replace(".fastq", "_{}.fastq".format(out_type))
        R1_handle = _get_fastq_handle(R1_filename, "w")
        R2_handle = _get_fastq_handle(R2_filename, "w")

        return (R1_handle, R2_handle)

    def _get_fastq_iterators(self):
        R1_handle = _get_fastq_handle(self.R1_filename)
        R2_handle = _get_fastq_handle(self.R2_filename)

        R1_iterator = FastqGeneralIterator(R1_handle)
        R2_iterator = FastqGeneralIterator(R2_handle)

        return (R1_iterator, R2_iterator)

    def _get_seq_records(self):
        R1_handle = _get_fastq_handle(self.R1_filename)
        R2_handle = _get_fastq_handle(self.R2_filename)
        R1 = list(SeqIO.parse(R1_handle, "fastq"))
        R2 = list(SeqIO.parse(R2_handle, "fastq"))
        R1_handle.close()
        R2_handle.close()
        return (R1, R2)

    def _check_names_filter_fastq(self, filter_by_qual=False):
        if self.min_average_read_quality > 0 or self.min_single_bp_quality > 0:
            info(
                "Filtering reads with average bp quality < {} \
                and single bp quality < {} ... in up to position {} of R1 and {} of R2".format(
                    self.min_average_read_quality, self.min_single_bp_quality
                )
            )

        if self.qend_R1 > 0 or self.qend_R2 > 0:
            info(
                "In the filtering, bases up to position {} of R1 and {} of R2 are considered.".format(
                    self.qend_R1, self.qend_R2
                )
            )

        R1, R2 = self._get_seq_records()

        _check_readname_match(R1, R2)
        if filter_by_qual:
            self.n_reads_after_filtering = self._filter_read_quality(R1, R2)
            if self.n_reads_after_filtering == 0:
                raise NoReadsAfterQualityFiltering(
                    "No reads in input or no reads survived the average or single bp quality filtering."
                )
            else:
                info(
                    "Number of reads in input:%d\tNumber of reads after filtering:%d\n"
                    % (self.n_total_reads, self.n_reads_after_filtering)
                )
        else:
            self.n_reads_after_filtering = self.n_total_reads

    def _filter_read_quality(self, R1=None, R2=None) -> int:
        R1_filtered = gzip.open(self.filtered_R1_filename, "w+")
        R2_filtered = gzip.open(self.filtered_R2_filename, "w+")

        if R1 is None or R2 is None:
            R1, R2 = self._get_seq_records()

        n_reads_after_filtering = 0
        for i in range(len(R1)):
            R1_record = R1[i]
            R2_record = R2[i]

            if self.filter_by_qual:
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

                if R1_quality_pass and R2_quality_pass:
                    n_reads_after_filtering += 1
                    R1_filtered.write(R1.format("fastq"))
                    R2_filtered.write(R2.format("fastq"))
        return n_reads_after_filtering

    def _write_start_log(self):
        try:
            os.makedirs(self.output_dir)
            print("Creating Folder %s" % self.output_dir)
        except:
            print("Folder %s already exists." % self.output_dir)
        self.log_filename = self._jp("CRISPRessoCount_RUNNING_LOG.txt")
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

        if not self.name:
            database_id = "%s" % get_name_from_fasta(self.R1_filename)
        else:
            database_id = self.name
        return database_id
