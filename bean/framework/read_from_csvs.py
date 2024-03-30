import argparse

from perturb_tools._readwrite._funcs._read_screen_from_csvs import read_csvs
from perturb_tools import Screen




def create_screen(args: argparse.Namespace) -> Screen:
    return read_csvs(
        args.gRNA_counts_table_csv,
        args.gRNA_info_table_csv,
        args.samples_info_table_csv,
        layers_filenames_dict={"edits": args.edits} if args.edits else None,
    )
