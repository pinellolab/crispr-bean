from typing import List, Union, Dict, Optional, Tuple, Iterable
import subprocess as sb
import numpy as np
import pandas as pd
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from bean.framework.Edit import Allele, Edit
from bean.mapping.CRISPResso2Align import read_matrix, global_align_base_editor


class InputFileError(Exception):
    pass


def _base_edit_to_from(start_base: str = "A"):
    try:
        base_map = {"A": "G", "C": "T"}
    except KeyError:
        raise ValueError("Only A/C are supported for base to be edited.")
    return base_map[start_base]


def _read_count_match(R1_filename: str, R2_filename: str) -> int:
    R1_count = _get_n_reads_fastq(R1_filename)
    R2_count = _get_n_reads_fastq(R2_filename)
    if R1_count != R2_count:
        raise InputFileError(
            "Paired end read numbers are different in R1({}) and R2({})".format(
                R1_count, R2_count
            )
        )
    return R1_count


def _get_n_reads_fastq(fastq_filename: str):
    p = sb.Popen(
        ("z" if fastq_filename.endswith(".gz") else "")
        + "cat < %s | wc -l" % fastq_filename,
        shell=True,
        stdout=sb.PIPE,
    )
    return int(float(p.communicate()[0]) / 4.0)


def _get_fastq_handle(fastq_filename: str, mode="r"):
    if fastq_filename.endswith(".gz"):
        if mode == "r":
            mode = "rt"
        elif mode == "w":
            mode = "wb"
        fastq_handle = gzip.open(fastq_filename, mode)
    else:
        fastq_handle = open(fastq_filename, mode)
    return fastq_handle


def _read_is_good_quality(
    record: SeqIO.SeqRecord, min_bp_quality=0, min_single_bp_quality=0, qend=-1
):
    mean_quality_pass = (
        np.array(record.letter_annotations["phred_quality"])[:qend].mean()
        >= min_bp_quality
    )
    min_quality_pass = (
        np.array(record.letter_annotations["phred_quality"])[:qend].min()
        >= min_single_bp_quality
    )
    return mean_quality_pass and min_quality_pass


def _check_readname_match(R1: List[SeqIO.SeqRecord], R2: List[SeqIO.SeqRecord]):
    if len(R1) != len(R2):
        raise ValueError("The number of reads in R1 and R2 file does not match.")

        for i in range(len(R1)):
            R1_record = R1[i]
            R2_record = R2[i]

            if R1_record.name != R2_record.name:
                raise InputFileError(
                    "R1 and R2 read discordance in read {} and {}".format(
                        R1_record.name, R2_record.name
                    )
                )


def _get_guide_to_reporter_df(sgRNA_filename: str) -> pd.DataFrame:
    """Returns a gRNA name to reporter sequence mapping."""

    with open(sgRNA_filename) as infile:
        sgRNA_df = pd.read_csv(infile)
        if not ("name" in sgRNA_df.columns and "Reporter" in sgRNA_df.columns):
            raise InputFileError(
                "Input gRNA file doesn't have the column 'gRNA' or 'gRNA_barcode'."
            )
        sgRNA_df.set_index("name", inplace=True)
        return sgRNA_df


def revcomp(seq: Union[Seq, str]) -> str:
    if isinstance(seq, str):
        seq = Seq(seq)
    return str(seq.reverse_complement())


def _fastq_iter_to_text(record: SeqIO.SeqRecord):
    t = (record.id,)
    seq = record.seq
    q = record.letter_annotations["phred_quality"]
    return "{}\n{}\n+\n{}\n".format(t, seq, q)


def _write_paired_end_reads(R1_record, R2_record, R1_out_handle, R2_out_handle):
    try:
        R1_out_handle.write(_fastq_iter_to_text(R1_record))
        R2_out_handle.write(_fastq_iter_to_text(R2_record))
    except TypeError:
        R1_out_handle.write(_fastq_iter_to_text(R1_record).encode())
        R2_out_handle.write(_fastq_iter_to_text(R2_record).encode())


def _get_edited_allele(
    ref_seq: str,
    query_seq: str,
    offset: int,
    strand: int = 1,
    start_pos: int = 0,
    end_pos: int = 100,
):
    allele = Allele()

    assert len(ref_seq) == len(query_seq), "reference and query seq length mismatch"

    for i, (ref_nt, sample_nt) in enumerate(zip(ref_seq, query_seq)):
        if i < start_pos or i >= end_pos:
            continue
        if ref_nt == sample_nt:
            continue
        else:
            edit = Edit(i - start_pos, ref_nt, sample_nt, offset=offset, strand=strand)
            allele.add(edit)

    return allele


def _write_alignment_matrix(
    target_base_edits: Dict[str, str], path, allow_complementary=False
):
    """
    Writes base substitution matrix
    """
    bases = ["A", "C", "T", "G"]
    for ref_base, alt_base in target_base_edits.items():
        if ref_base not in bases or alt_base not in bases:
            raise ValueError(
                f"Specified ref base '{ref_base}' or alt base '{alt_base}' isn't valid"
            )
    mat = np.ones((4, 4), dtype=int) * -4
    np.fill_diagonal(mat, 5)
    aln_df = pd.DataFrame(mat, index=bases, columns=bases)
    for ref_base, alt_base in target_base_edits.items():
        aln_df.loc[ref_base, alt_base] = 0
        if allow_complementary:
            comp_map = {"A": "T", "C": "G", "T": "A", "G": "C"}
            aln_df.loc[comp_map[ref_base], comp_map[alt_base]] = 0
    aln_df.to_csv(path, sep=" ")


def _get_allele_from_alignment(
    ref_aligned: str,
    query_aligned: str,
    offset: int,
    strand: int,
    start_pos: int,
    end_pos: int,
    chrom: Optional[str] = None,
    positionwise_quality: Optional[np.ndarray] = None,
    quality_thres: float = -1,
):
    assert len(ref_aligned) == len(query_aligned)
    allele = Allele()
    ref_gaps = 0
    alt_gaps = 0
    alt_seq_len = len(query_aligned) - query_aligned.count("-")
    if positionwise_quality is None:
        # alt_position_is_good_quality = np.ones(alt_seq_len, dtype=bool)
        alt_position_is_good_quality = np.array(
            [c != "N" for c in query_aligned.replace("-", "")]
        )
    else:
        assert len(positionwise_quality) == alt_seq_len
        alt_position_is_good_quality = positionwise_quality > quality_thres
    for i in range(len(ref_aligned)):
        if ref_aligned[i] == query_aligned[i]:
            continue
        ref_base = ref_aligned[i]
        alt_base = query_aligned[i]
        if alt_base != "-":
            alt_base_is_good_quality = alt_position_is_good_quality[i - alt_gaps]
        else:
            alt_base_is_good_quality = True
        if ref_base == "-":
            ref_gaps += 1
        elif alt_base == "-":
            alt_gaps += 1
        ref_pos = i - ref_gaps
        if alt_base_is_good_quality and ref_pos >= start_pos and ref_pos < end_pos:
            allele.add(
                Edit(
                    chrom=chrom,
                    rel_pos=ref_pos,
                    ref_base=ref_base,
                    alt_base=alt_base,
                    offset=offset,
                    strand=strand,
                )
            )
    return allele


def _get_edited_allele_crispresso(
    ref_seq: str,
    query_seq: str,
    target_base_edits: Dict[str, str],
    aln_mat_path: str,
    offset: int,
    strand: int = 1,
    chrom: Optional[str] = None,
    start_pos: int = 0,
    end_pos: int = 100,
    positionwise_quality: Optional[np.ndarray] = None,
    quality_thres: float = 30,
    objectify_allele=True,
) -> Tuple[Union[Allele, str], float]:
    aln_matrix = read_matrix(aln_mat_path)
    assert strand in [-1, +1]
    gap_incentive = np.zeros(len(ref_seq) + 1, dtype=int)
    query_aligned, ref_aligned, score = global_align_base_editor(
        query_seq,
        ref_seq,
        target_base_edits,
        aln_matrix,
        gap_incentive,
        gap_open=-20,
        gap_extend=-10,
    )
    if objectify_allele:
        allele = _get_allele_from_alignment(
            ref_aligned,
            query_aligned,
            offset,
            strand,
            start_pos,
            end_pos,
            chrom,
            positionwise_quality,
            quality_thres,
        )
        for e in allele.edits:
            if e.ref_base == "-":
                continue
            assert ref_seq[e.rel_pos] == e.ref_base, (
                "relative position mismatch: ref pos {}: {} vs {}".format(
                    e.rel_pos, ref_seq[e.rel_pos], e.ref_base
                )
                + "\nallele {}, \nrefseq {}, \naltseq {}, \n".format(
                    allele, ref_seq, query_seq
                )
                + "ref_align {}, \nalt_align {}".format(ref_aligned, query_aligned)
            )
    else:
        allele = _string_filter_basewise_quality(
            ref_aligned, query_aligned, positionwise_quality, quality_thres
        )
    return (allele, score)


def _string_filter_basewise_quality(
    ref_seq, query_seq, positionwise_quality, quality_thres
):
    query_seq_chars = list(query_seq)
    for i in range(len(positionwise_quality)):
        if positionwise_quality[i] < quality_thres:
            query_seq_chars[i] = ref_seq[i]
    return "".join(query_seq_chars)


def _multiindex_dict_to_df(input_dict, key_column_names, value_column_name):
    if not isinstance(key_column_names, list):
        key_column_names = [key_column_names]
    mi = pd.MultiIndex.from_tuples(
        input_dict.keys(), names=["guide"] + key_column_names
    )
    df = pd.DataFrame.from_dict(
        input_dict,
        orient="index",
        columns=[value_column_name],
    )
    df.index = mi
    df.reset_index(inplace=True)
    if len(key_column_names) == 1:
        df.rename(
            columns={"level_0": "guide", "level_1": key_column_names[0]}, inplace=True
        )
    elif len(key_column_names) == 2:
        df.rename(
            columns={
                "level_0": "guide",
                "level_1": key_column_names[0],
                "level_2": key_column_names[1],
            },
            inplace=True,
        )
    return df


def hamming_distance(
    seq1: str,
    seq2: str,
    match_score: float = 0,
    mismatch_penalty: float = 1,
    allowed_substitutions_penalties: Dict[Tuple[str, str], float] = {
        ("A", "G"): 0.2,
        ("T", "C"): 0.2,
    },
):
    """
    Calculates the Hamming distance between two DNA sequences with different penalties.

    Args:
        seq1 (str): The first DNA sequence.
        seq2 (str): The second DNA sequence.
        match_score (int): Score for a match (default: 0).
        mismatch_penalty (int): Penalty for a mismatch (default: 1).

    Returns:
        int: The Hamming distance.
    """

    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length.")

    distance = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            distance += match_score
        elif (seq1[i], seq2[i]) in allowed_substitutions_penalties:
            distance += allowed_substitutions_penalties[(seq1[i], seq2[i])]
        else:
            distance += mismatch_penalty
    return distance


def find_closest_sequence_index(
    query_seqs: Iterable[str],
    ref_seqs: Iterable[str],
    hamming_distance_threshold: int = 0.1,
    match_score: float = 0,
    mismatch_penalty: float = 1,
    allowed_substitutions_penalties: Dict[Tuple[str, str], float] = {
        ("A", "G"): 0.2,
        ("T", "C"): 0.2,
    },
) -> int:
    """
    Find the closest sequence to a query sequence.

    Args:
        query_seq (str): The query sequence.
        ref_seqs (Iterable[str]): The reference sequences.
        hamming_distance_threshold (int, optional): The maximum allowed hamming distance. Defaults to 3.
        match_score (int, optional): Score for a match (default: 0). Defaults to 0.
        mismatch_penalty (int, optional): Penalty for a mismatch (default: 1). Defaults to 1.
        allowed_substitutions_penalties (Dict[Tuple[str, str], float], optional): Allowed substitutions and penalties. Defaults to [("A", "G"):0.2, ("T", "C"):0.2].

    Returns:
        int: The closest sequence pair's index in the input ref_seqs list.
    """
    min_distance = float("inf")
    closest_seq_index = None
    for i, (query_seq, ref_seq) in enumerate(zip(query_seqs, ref_seqs)):
        distance = hamming_distance(
            query_seq,
            ref_seq,
            match_score,
            mismatch_penalty,
            allowed_substitutions_penalties,
        ) / len(query_seq)
        if distance < min_distance and distance <= hamming_distance_threshold:
            min_distance = distance
            closest_seq_index = i
    return closest_seq_index
