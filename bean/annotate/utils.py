from typing import Tuple
import os
import sys
from pathlib import Path
import requests
from typing import Optional, List
import argparse
import numpy as np
import pandas as pd
from itertools import chain
import logging

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

complement_base = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"}


def revcomp(nt_list: List[str]):
    rev_list = nt_list[::-1]
    return list(map(lambda c: complement_base[c], rev_list))


def fast_flatten(input_list):
    return list(chain.from_iterable(input_list))


def fast_concat(df_list: List[pd.DataFrame]):
    """Faster concatenation of many dataframes from
    https://gist.github.com/TariqAHassan/fc77c00efef4897241f49e61ddbede9e
    """
    colnames = df_list[0].columns
    df_dict = dict.fromkeys(colnames, [])
    for col in colnames:
        extracted = (df[col] for df in df_list)
        # Flatten and save to df_dict
        df_dict[col] = fast_flatten(extracted)
    return pd.DataFrame.from_dict(df_dict)[colnames]


def find_overlap(
    chrom: str, start: int, end: int, range_df: pd.DataFrame
) -> Optional[str]:
    """Find overlap between query range and range_df and return ID of overlapping region in range_df."""
    if chrom not in range_df.chrom.tolist():
        return None
    overlap = range_df.loc[
        (range_df.chrom == chrom)
        & (
            ((start <= range_df.end) & (end >= range_df.start))
            | ((end >= range_df.start) & (start <= range_df.end))
        )
    ]
    if len(overlap) == 0:
        return None
    if len(overlap) > 2:
        raise ValueError(
            f"We cannot handle overlapping genes for {chrom}, {start}, {end} with : {overlap}"
        )
    return overlap.index.tolist()[0]


def get_mane_transcript_id(gene_name: str):
    """
    Retrieves the MANE transcript ID and version for a given gene name.

    Args:
        gene_name (str): The gene name for which to retrieve the MANE transcript ID.

    Returns:
        tuple: A tuple containing the MANE transcript ID and version.
    """
    api_url = "http://tark.ensembl.org/api/transcript/manelist/"
    response = requests.get(api_url, headers={"Content-Type": "application/json"})
    mane_json = response.json()
    mane_df = pd.DataFrame.from_records(mane_json)
    try:
        mane_transcript_id = mane_df.loc[
            mane_df.ens_gene_name == gene_name, "ens_stable_id"
        ].values[0]
        id_version = mane_df.loc[
            mane_df.ens_gene_name == gene_name, "ens_stable_id_version"
        ].values[0]
    except IndexError as e:
        print(
            f"Cannot find {gene_name} from MANE database: check http://tark.ensembl.org/api/transcript/manelist/ or use custom fasta."
        )
        print(gene_name)
        print(e)
        exit(1)
    return mane_transcript_id, id_version


def get_exons_from_transcript_id(
    transcript_id: str, id_version: int, ref_version: str = "GRCh38"
):
    """
    Retrieves the exons and the start position of the coding sequence (CDS) for a given transcript ID and version.

    Args:
        transcript_id (str): The transcript ID for which to retrieve the exons.
        id_version (int): The version of the transcript ID.

    Returns:
        tuple: A tuple containing the exons and the start and end (inclusive) position of the CDS.
    """
    api_url = f"http://tark.ensembl.org/api/transcript/?stable_id={transcript_id}&stable_id_version={id_version}&expand=exons"
    response = requests.get(api_url, headers={"Content-Type": "application/json"})
    transcript_json = response.json()
    if transcript_json["count"] != 1:
        if transcript_json["count"] > 1:
            api_url = f"http://tark.ensembl.org/api/transcript/?stable_id={transcript_id}&stable_id_version={id_version}&assembly_name={ref_version}&expand=exons"
            response = requests.get(
                api_url, headers={"Content-Type": "application/json"}
            )
            transcript_json = response.json()
            if transcript_json["count"] != 1:
                raise ValueError(
                    f"Non-unique entry for transcript ID , version {id_version} and assembly {ref_version}:\n{transcript_json}"
                )
        else:
            raise ValueError(
                f"No entry found for transcript ID {transcript_id}, version {id_version} and assembly {ref_version}"
            )
    transcript_record = transcript_json["results"][0]
    exons_list = transcript_record["exons"]
    strand = transcript_record["loc_strand"]  # +1/-1
    if strand == 1:
        cds_start = transcript_record["five_prime_utr_end"] + 1
        cds_end = transcript_record["three_prime_utr_start"] - 1
    else:
        assert strand == -1, f"Invalid loc_strand returned from {api_url}: {strand}"
        cds_start = transcript_record["three_prime_utr_start"] + 1
        cds_end = transcript_record["five_prime_utr_end"] - 1

    return exons_list, cds_start, cds_end, strand


def get_exon_pos_seq(exon_id, id_version, cds_start, cds_end, ref_version="GRCh38"):
    seq = []
    genomic_pos = []
    api_url = f"http://tark.ensembl.org/api/exon/?stable_id={exon_id}&stable_id_version={id_version}&expand=sequence"
    response = requests.get(api_url, headers={"Content-Type": "application/json"})
    exon_json = response.json()
    if exon_json["count"] != 1:
        if exon_json["count"] > 1:
            exon_record = None
            for record in exon_json["results"]:
                if record["assembly"] == ref_version:
                    if exon_record is not None:
                        raise ValueError(
                            f"Non-unique entry for exon ID and version:\n{exon_json}"
                        )
                    exon_record = record
            if exon_record is None or exon_json["count"] == 0:
                raise ValueError(
                    f"Non-unique entry for exon ID and version:\n{exon_json}"
                )
    else:
        exon_record = exon_json["results"][0]
    sequence = exon_record["sequence"]["sequence"]
    strand = exon_record["loc_strand"]
    start_pos = exon_record["loc_start"]
    end_pos = exon_record["loc_end"]
    chrom = exon_record["loc_region"]
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    if strand == -1:
        sequence = revcomp([*sequence])

    if cds_start > end_pos or cds_end < start_pos:
        return chrom, [], []
    if cds_start > start_pos and cds_start < end_pos:
        if cds_start > end_pos:
            warn(f"Exon {exon_id} doesn't have coding sequence.")
            return chrom, seq, genomic_pos
        else:
            sequence = sequence[cds_start - start_pos :]
            start_pos = cds_start
    if cds_end > start_pos and cds_end < end_pos:
        sequence = sequence[: (cds_end - end_pos)]
        end_pos = cds_end

    assert len(sequence) == end_pos - start_pos + 1, (
        len(sequence),
        cds_start,
        cds_end,
        end_pos,
        start_pos,
        (end_pos - start_pos + 1),
        sequence,
    )
    return chrom, list(sequence), list(range(start_pos, end_pos + 1))


def get_cds_seq_pos_from_gene_name(gene_name: str, ref_version: str = "GRCh38"):
    transcript_id, id_version = get_mane_transcript_id(gene_name)
    print(f"MANE transcript ID {transcript_id} for {gene_name} will be used.")
    exons_list, cds_start, cds_end, strand = get_exons_from_transcript_id(
        transcript_id, id_version, ref_version
    )
    if strand == -1:
        exons_list = exons_list[::-1]
    cds_seq = []
    cds_pos = []
    for exon_dict in exons_list:
        cds_chrom, _cds_seq, _cds_pos = get_exon_pos_seq(
            exon_dict["stable_id"], exon_dict["stable_id_version"], cds_start, cds_end
        )
        cds_seq.extend(_cds_seq)
        cds_pos.extend(_cds_pos)
    return cds_chrom, cds_seq, cds_pos, strand


def get_splice_positions_from_gene_name(
    gene_name: str, ref_version: str = "GRCh38"
) -> Tuple[str, np.ndarray, np.ndarray]:
    transcript_id, id_version = get_mane_transcript_id(gene_name)
    print(f"MANE transcript ID {transcript_id} for {gene_name} will be used.")
    exons_list, cds_start, cds_end, strand = get_exons_from_transcript_id(
        transcript_id, id_version, ref_version
    )
    if strand == -1:
        exons_list = exons_list[::-1]
    splice_donor_pos = []
    splice_acceptor_pos = []
    cds_positions = []
    for exon_dict in exons_list:
        cds_chrom, _cds_seq, _cds_pos = get_exon_pos_seq(
            exon_dict["stable_id"], exon_dict["stable_id_version"], cds_start, cds_end
        )
        if len(_cds_pos) > 0:
            cds_positions.append(_cds_pos)
    for i, _cds_pos in enumerate(cds_positions):
        if i > 0:
            splice_acceptor_pos.append(_cds_pos[0])
        if i < len(exons_list):
            splice_donor_pos.append(_cds_pos[-1])
    return cds_chrom, np.array(splice_donor_pos), np.array(splice_acceptor_pos)


def get_splice_parser(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            "Get splice site position",
            usage="Get splice site position from exon fasta and target editing base.",
        )

    parser.add_argument(
        "exon_fa_path", help="File path to fasta file with exon position information."
    )
    parser.add_argument(
        "--gene-name",
        action="store_true",
        help="File path to fasta file with exon position information.",
    )
    parser.add_argument("edited_base", help="Edited base, either A or C.")
    parser.add_argument("output_path", help="output path of the splice site csv file.")

    return parser


def parse_args(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            prog="allele_filter",
            description="Filter alleles based on edit position in spacer and frequency across samples.",
        )
    parser.add_argument(
        "bdata_path",
        type=str,
        help="Input ReporterScreen file of which allele will be filtered out.",
    )
    parser.add_argument(
        "--output-prefix",
        "-o",
        type=str,
        default=None,
        help="Output prefix for log and ReporterScreen file with allele assignment",
    )
    parser.add_argument(
        "--plasmid-path",
        "-p",
        type=str,
        default=None,
        help="Plasmid ReporterScreen object path. If provided, alleles are filtered based on if a nucleotide edit is more significantly enriched in sample compared to the plasmid data. Negative control data where no edit is expected can be fed in instead of plasmid library.",
    )
    parser.add_argument(
        "--reporter-length",
        type=int,
        default=32,
        help="Length of reporter sequence in the construct.",
    )
    parser.add_argument(
        "--reporter-right-flank-length",
        type=int,
        default=6,
        help="Length of the right-flanking nucleotides of protospacer in the reporter.",
    )
    parser.add_argument(
        "--edit-start-pos",
        "-s",
        type=int,
        default=2,
        help="0-based start posiiton (inclusive) of edit relative to the start of guide spacer.",
    )
    parser.add_argument(
        "--edit-end-pos",
        "-e",
        type=int,
        default=7,
        help="0-based end position (exclusive) of edit relative to the start of guide spacer.",
    )
    parser.add_argument(
        "--jaccard-threshold",
        "-j",
        type=float,
        help="Jaccard Index threshold when the alleles are mapped to the most similar alleles. In each filtering step, allele counts of filtered out alleles will be mapped to the most similar allele only if they have Jaccard Index of shared edit higher than this threshold.",
        default=0.3,
    )
    parser.add_argument(
        "--filter-spacer",
        help="Only consider edit within protospacer positions of reporter.",
        action="store_true",
    )
    parser.add_argument(
        "--filter-window",
        "-w",
        help="Only consider edit within window provided by (edit-start-pos, edit-end-pos). If this flag is not provided, `--edit-start-pos` and `--edit-end-pos` flags are ignored.",
        action="store_true",
    )
    parser.add_argument(
        "--keep-indels",
        "-i",
        help="Include indels.",
        action="store_true",
    )
    parser.add_argument(
        "--filter-target-basechange",
        "-b",
        help="Only consider target edit (stored in bdata.uns['target_base_changes'])",
        action="store_true",
    )
    parser.add_argument(
        "--translate", "-t", help="Translate alleles", action="store_true"
    )
    parser.add_argument(
        "--translate-fasta",
        "-f",
        type=str,
        help="fasta file path with exon positions. If not provided, LDLR hg19 coordinates will be used.",
        default=None,
    )
    parser.add_argument(
        "--translate-fastas-csv",
        "-fs",
        type=str,
        help=".csv with two columns with gene IDs and FASTA file path corresponding to each gene.",
        default=None,
    )
    parser.add_argument(
        "--translate-gene",
        "-g",
        type=str,
        help="Gene symbol if a gene is tiled. If not provided, LDLR hg19 coordinates will be used.",
        default=None,
    )
    parser.add_argument(
        "--translate-genes-list",
        "-gs",
        type=str,
        help="File with gene symbols, one per line, if multiple genes are tiled.",
        default=None,
    )
    parser.add_argument(
        "--filter-allele-proportion",
        "-ap",
        type=float,
        default=0.05,
        help="If provided, alleles that exceed `filter_allele_proportion` in `filter-sample-proportion` will be retained.",
    )
    parser.add_argument(
        "--filter-allele-count",
        "-ac",
        type=int,
        default=5,
        help="If provided, alleles that exceed `filter_allele_proportion` AND `filter_allele_count` in `filter-sample-proportion` will be retained.",
    )
    parser.add_argument(
        "--filter-sample-proportion",
        "-sp",
        type=float,
        default=0.2,
        help="If `filter_allele_proportion` is provided, alleles that exceed `filter_allele_proportion` in `filter-sample-proportion` will be retained.",
    )
    parser.add_argument(
        "--load-tmp",
        action="store_true",
        help="Load temporary file and work from there.",
    )
    return parser


def check_args(args):
    if args.output_prefix is None:
        args.output_prefix = args.bdata_path.rsplit(".h5ad", 1)[0] + "_alleleFiltered"
    info(f"Saving results to {args.output_prefix}")
    Path(os.path.dirname(args.output_prefix)).mkdir(parents=True, exist_ok=True)
    if args.filter_window:
        if args.edit_start_pos is None and args.edit_end_pos is None:
            raise ValueError(
                "Invalid arguments: --filter-window option set but none of --edit-start-pos and --edit-end-pos specified."
            )
        if args.edit_start_pos is None:
            warn(
                "--filter-window option set but none of --edit-start-pos not provided. Using 0 as its value."
            )
            args.edit_start_pos = 0
        if args.edit_end_pos is None:
            warn(
                "--filter-window option set but none of --edit-end-pos not provided. Using 20 as its value."
            )
            args.edit_end_pos = 20
    if args.filter_allele_proportion is not None and (
        args.filter_allele_proportion < 0 or args.filter_allele_proportion > 1
    ):
        raise ValueError(
            "Invalid arguments: filter-allele-proportion should be in range [0, 1]."
        )
    if args.filter_sample_proportion < 0 or args.filter_sample_proportion > 1:
        raise ValueError(
            "Invalid arguments: filter-sample-proportion should be in range [0, 1]."
        )
    if args.translate and (
        int(args.translate_fasta is not None)
        + int(args.translate_fastas_csv is not None)
        + int(args.translate_gene is not None)
        + int(args.translate_genes_list is not None)
        > 1
    ):
        raise ValueError(
            "Invalid arguments: You should specify exactly one of --translate-fasta, --translate-fastas-csv, --translate-gene, translate-genes-list to translate alleles."
        )
    if args.translate_genes_list is not None:
        args.translate_genes_list = (
            pd.read_csv(args.translate_genes_list, header=None).values[:, 0].tolist()
        )
        info(f"Using {args.translate_genes_list} as genes for translation.")
    if args.translate_fastas_csv:
        tbl = pd.read_csv(
            args.translate_fastas_csv,
            header=None,
        )
        if len(tbl) == 0 or len(tbl.columns != 2):
            raise ValueError(
                "Invalid arguments: Table should have two columns and more than 0 entry"
            )
        for path in tbl.iloc[:, 2].tolist():
            if not os.path.isfile(path):
                raise FileNotFoundError(
                    f"Invalid input file: {path} does not exist. Check your input in {args.translate_fastas_csv}"
                )

    return args
