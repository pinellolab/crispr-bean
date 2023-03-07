import argparse
import gzip

import pandas as pd
from Bio import SeqIO
from slugify import slugify


class InputFileError(Exception):
    """Raised when the input file is not valid."""

    pass


def _check_file(filename):
    """Check if file exists and is readable."""
    try:
        with open(filename, "r"):  # pylint: disable=unspecified-encoding
            pass
    except IOError as exc:
        raise InputFileError(f"I cannot open the file: {filename}") from exc


def _check_library(library_name):
    """Check if library exists."""
    try:
        return __import__(library_name)
    except ImportError as exc:
        raise ImportError(
            f"You need to install {library_name} module to use CRISPRessoCount!"
        ) from exc


def _get_input_parser():
    """Get the input data"""
    print(
        r"""
    _ _       
  /  \ '\                       _   
  |   \  \      __ ___ _  _ _ _| |_ 
   \   \  |    / _/ _ \ || | ' \  _|
    `.__|/     \__\___/\_,_|_||_\__|
    """
    )

    parser = argparse.ArgumentParser(
        description="CRISPRessoCount parameters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-b",
        "--edited-base",
        type=str,
        required=True,
        help="For base editors, the base that should be ignored when matching the gRNA sequence",
    )
    parser.add_argument(
        "-f",
        "--sgRNA-filename",
        type=str,
        required=True,
        help="""sgRNA description file. The format requires three columns: gRNA, Reporter, gRNA_barcode.""",
    )

    # optional
    parser.add_argument(
        "--guide-start-seq",
        type=str,
        help="Guide starts after this sequence in R1",
        default="",
    )
    parser.add_argument(
        "--guide-end-seq",
        type=str,
        help="Guide starts after this sequence in R1",
        default="",
    )
    parser.add_argument(
        "-r", "--count-reporter", help="Count reporter edits.", action="store_true"
    )
    parser.add_argument(
        "-q",
        "--min_average_read_quality",
        type=int,
        help="Minimum average quality score (phred33) to keep a read",
        default=30,
    )
    parser.add_argument(
        "-s",
        "--min-single-bp-quality",
        type=int,
        help="Minimum single bp score (phred33) to keep a read",
        default=0,
    )
    parser.add_argument("-n", "--name", help="Output name", default="")
    parser.add_argument("-o", "--output-folder", help="", default="")
    parser.add_argument(
        "-l", "--reporter-length", type=int, help="length of the reporter", default=32
    )
    parser.add_argument(
        "--keep-intermediate",
        help="Keep all the  intermediate files",
        action="store_true",
    )
    parser.add_argument(
        "--qstart-R1",
        help="Start position of the read when filtering for quality score of the read 1",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--qend-R1",
        help="End position of the read when filtering for quality score of the read 1",
        type=int,
        default=47,
    )
    parser.add_argument(
        "--qstart-R2", help="Same as qstart_R1, for read 2 fastq file", default=0
    )
    parser.add_argument(
        "--qend-R2", help="Same as qstart_R2, for read 2 fastq file", default=36
    )
    parser.add_argument(
        "--gstart-reporter",
        help="Start position of the guide sequence in the reporter",
        type=int,
        default=6,
    )
    parser.add_argument(
        "--match-target-pos",
        help="Count the edit in the exact target position.",
        action="store_true",
    )
    parser.add_argument(
        "--target-pos-col",
        help="Column name specifying the relative target position within reporter sequence.",
        default="target_pos",
    )

    parser.add_argument("--guide-bc", help="Construct has guide barcode", default=True)
    parser.add_argument(
        "--guide-bc-len",
        help="Guide barcode sequence length at the beginning of the R2",
        type=str,
        default=4,
    )
    parser.add_argument(
        "--offset",
        help="Guide file has offest column that will be added to the relative position of reporters.",
        action="store_true",
    )
    parser.add_argument(
        "--align-fasta",
        help="gRNA is aligned to this sequence to infer the offset. Can be used when the exact offset is not provided.",
        type=str,
        default="",
    )

    parser.add_argument(
        "-as",
        "--string-allele",
        help="Store allele as quality filtered string instead of Allele object",
        action="store_true",
    )
    parser.add_argument(
        "-g",
        "--count-guide-edits",
        help="count the self editing of guides",
        action="store_true",
    )
    parser.add_argument(
        "-m",
        "--count-guide-reporter-alleles",
        help="count the matched allele of guide and reporter edit",
        action="store_true",
    )

    return parser


def _get_first_read_length(fastq_filename: str) -> int:
    """Obtain length of first read entry of fastq file"""
    if fastq_filename.split(".")[-1] == "gz":
        handle = gzip.open(fastq_filename, "rt")
    else:
        handle = fastq_filename

    for record in SeqIO.parse(handle, "fastq"):
        return len(record)
    raise InputFileError("Provided R1 file doesn't have any read to parse")


def _check_read_length(args: argparse.Namespace, read_length: int, warn_logger):
    """Checks if quality filter positions are compatible with provided read length."""
    if args.qstart_R1 >= read_length or args.qstart_R2 >= read_length:
        raise ValueError(
            f"The start position of base quality filter is not nonnegative ({args.qstart_R1} for R1, {args.qstart_R2} for R2 provided)"
        )
    if args.qend_R1 >= read_length or args.qend_R2 >= read_length:
        raise ValueError(
            f"The start position of base quality filter is not nonnegative ({args.qstart_R1} for R1, {args.qstart_R2} for R2 provided)"
        )
    if args.qend_R2 != args.guide_bc_len + args.reporter_length:
        warn_logger(
            f"Quality of R2 checked up until {args.qend_R2}bp, while the length of guide barcode and reporter combined is {args.guide_bc_len + args.reporter_length}bp."
        )


def _check_arguments(args, info_logger, warn_logger, error_logger):
    """Check the argument validity of the ArgumentParser"""

    if args.sgRNA_filename:
        _check_file(args.sgRNA_filename)

    # Edited base should be one of A/C/T/G
    if args.edited_base.upper() not in ["A", "C", "T", "G"]:
        raise ValueError(
            f"The edited base should be one of A/C/T/G, {args.edited_base} provided."
        )

    edited_base = args.edited_base.upper()
    info_logger(f"Using specified edited base: {edited_base}")
    info_logger(
        f"Using guide barcode length {args.guide_bc_len}, guide start '{args.guide_start_seq}'"
    )
    # normalize name and remove not allowed characters
    # if args.name:
    #     clean_name = slugify(args.name)
    #     if args.name != clean_name:
    #         warn_logger(
    #             f"The specified name {args.name} contained characters not allowed and was changed to: {clean_name}"
    #         )
    #         args.name = clean_name
    sgRNA_info_tbl = pd.read_csv(args.sgRNA_filename)

    def _check_sgrna_info_table(args, sgRNA_info_tbl):
        if args.offset:
            if "offset" not in sgRNA_info_tbl.columns:
                raise InputFileError(
                    "Offset option is set but the input file doesn't contain the offset column."
                )
            if len(args.align_fasta) > 0:
                error_logger("Can't have --offset and --align_fasta option together.")
        if (
            args.match_target_pos
            and args.target_pos_col is not None
            and args.target_pos_col not in sgRNA_info_tbl.columns
        ):
            raise InputFileError(
                f"Specified target position column '{args.target_pos_col}' not in the input file {args.sgRNA_filename}."
            )

    _check_sgrna_info_table(args, sgRNA_info_tbl)

    if args.match_target_pos and (args.target_pos_col not in sgRNA_info_tbl.columns):
        raise InputFileError(
            "Target position option is set but the input file doesn't contain the target position column."
        )

    info_logger("Done checking input arguments.")

    return args
