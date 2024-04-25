from typing import Tuple
import argparse
import re
import numpy as np
import pandas as pd
from bean.annotate.utils import get_splice_positions_from_gene_name


def get_parser(parser=None):
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


def get_splice_positions(exon_fa_path) -> Tuple[str, np.ndarray, np.ndarray]:
    splice_donor_pos = []
    splice_acceptor_pos = []
    p = re.compile("range=chr(\d+):(\d+)-(\d+)")
    with open(exon_fa_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                result = p.search(line)
                chrom = result[1]
                exon_start = int(result[2])
                exon_end = int(result[3])
                splice_donor_pos.append(exon_end)
                splice_acceptor_pos.append(exon_start)
    splice_donor_pos = np.array(splice_donor_pos)
    splice_acceptor_pos = np.array(splice_acceptor_pos)
    splice_donor_pos = splice_donor_pos[:-1]
    splice_acceptor_pos = splice_acceptor_pos[1:]
    return chrom, splice_donor_pos, splice_acceptor_pos


def get_targetable_splice_positions(
    chrom: str,
    splice_donor_pos: np.ndarray,
    splice_acceptor_pos: np.ndarray,
    edited_base: str = "A",
) -> pd.DataFrame:
    """
    Splice donor: GT
    Splice acceptor: AG
    """
    splice_site_dfs = []
    revcomp_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    splice_donor_consensus = {"G": 1, "T": 2}
    splice_acceptor_consensus = {"A": -2, "G": -1}
    for splice_site_pos, rel_basepos_map, label in zip(
        [splice_donor_pos, splice_acceptor_pos],
        [splice_donor_consensus, splice_acceptor_consensus],
        ["SD", "SA"],
    ):
        splice_site_dfs.extend(
            pd.DataFrame(
                {
                    "chrom": chrom,
                    "pos": splice_site_pos + rel_basepos_map[base],
                    "type": label,
                    "target_base": base,
                }
            )
            for base in [edited_base, revcomp_map[edited_base]]
            if base in rel_basepos_map
        )
    return pd.concat(splice_site_dfs)


def main(args):
    if args.gene_name is None:
        chrom, sd_pos, sa_pos = get_splice_positions(args.exon_fa_path)
    else:
        chrom, sd_pos, sa_pos = get_splice_positions_from_gene_name(
            args.exon_fa_path.strip()
        )
    splice_target_df = get_targetable_splice_positions(
        chrom, sd_pos, sa_pos, args.edited_base
    )
    splice_target_df.to_csv(args.output_path)
