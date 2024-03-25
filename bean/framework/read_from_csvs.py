import argparse

from perturb_tools._readwrite._funcs._read_screen_from_csvs import read_csvs
from perturb_tools import Screen


def get_input_parser() -> argparse.Namespace:
    """Add multi-sample specific arguments to the base parser."""
    print(
        r"""
    _ _       
  /  \ '\                      _        
  |   \  \     __ _ _ ___ __ _| |_ ___ 
   \   \  |   / _| '_/ -_) _` |  _/ -_)
    `.__|/    \__|_| \___\__,_|\__\___|
    """
    )
    parser = argparse.ArgumentParser(
        description="bean-create-screen parameters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "gRNA_info_table_csv",
        type=str,
        help="Path to gRNA info .csv table, with index at first column and column names at the first row.",
    )
    parser.add_argument(
        "samples_info_table_csv",
        type=str,
        help="Path to sample info .csv table, with index at first column and column names at the first row.",
    )
    parser.add_argument(
        "gRNA_counts_table_csv",
        type=str,
        help="Path to gRNA counts .csv table, with index at first column and column names at the first row.",
    )
    parser.add_argument(
        "-e",
        "--edits",
        type=str,
        help="Path to edit counts .csv table, with index at first column and column names at the first row.",
        default=None,
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        type=str,
        help="Output file prefix (output will be saved as `output_prefix.h5ad`). If not provided, `gRNA_counts_table_csv` file prefix is used.",
        default=None,
    )
    return parser


def create_screen(args: argparse.Namespace) -> Screen:
    return read_csvs(
        args.gRNA_counts_table_csv,
        args.gRNA_info_table_csv,
        args.samples_info_table_csv,
        layers_filenames_dict={"edits": args.edits} if args.edits else None,
    )
