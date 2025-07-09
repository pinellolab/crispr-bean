import argparse

from perturb_tools._readwrite._funcs._read_screen_from_csvs import read_csvs
from perturb_tools import Screen




def create_screen(args: argparse.Namespace) -> Screen:
    layers_filenames_dict = None
    if args.edits or args.X_bcmatch:
        layers_filenames_dict = {}
        if args.edits:
            layers_filenames_dict["edits"] = args.edits
        if args.X_bcmatch:
            layers_filenames_dict["X_bcmatch"] = args.X_bcmatch
    return read_csvs(
        args.gRNA_counts_table_csv,
        args.gRNA_info_table_csv,
        args.samples_info_table_csv,
        layers_filenames_dict=layers_filenames_dict,
        matched_indices = False,
    )
