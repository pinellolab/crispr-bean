#!/usr/bin/env python
"""Count guides, optionally with reporter and alleles for multiple samples that share the same guide library."""


import argparse
import logging
from copy import deepcopy
import os
import sys
from multiprocessing import Pool

import bean
import pandas as pd
from bean.mapping.utils import (
    InputFileError,
    _check_arguments,
    _get_first_read_length,
    _check_read_length,
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


def count_sample(R1: str, R2: str, sample_id: str, args: argparse.Namespace):
    """Count single sample given R1 and R2 paths.
    Arguments are modified accordingly to the provided sample_id before being passed to GuideEditCounter.
    """
    args_dict = deepcopy(vars(args))
    args_dict["R1"] = R1
    args_dict["R2"] = R2
    args_dict["name"] = sample_id
    args_dict["output_folder"] = os.path.join(args.output_folder, sample_id)

    base_editing_map = {"A": "G", "C": "T"}
    try:
        target_base_edits = {k: base_editing_map[k] for k in args_dict["edited_base"]}
    except KeyError as e:
        raise KeyError(args_dict["edited_base"]) from e

    match_target_pos = args_dict["match_target_pos"]
    if (
        "guide_start_seqs_tbl" in args_dict
        and args_dict["guide_start_seqs_tbl"] is not None
    ):
        args_dict["guide_start_seq"] = str(args_dict["guide_start_seqs_tbl"][sample_id])
    if (
        "guide_end_seqs_tbl" in args_dict
        and args_dict["guide_end_seqs_tbl"] is not None
    ):
        args_dict["guide_end_seq"] = args_dict["guide_end_seqs_tbl"][sample_id]
    if (
        "barcode_start_seqs_tbl" in args_dict
        and args_dict["barcode_start_seqs_tbl"] is not None
    ):
        args_dict["barcode_start_seq"] = str(
            args_dict["barcode_start_seqs_tbl"][sample_id]
        )
    counter = bean.mp.GuideEditCounter(**args_dict)
    if os.path.exists(f"{counter.output_dir}.h5ad") and not args_dict["rerun"]:
        screen = bean.read_h5ad(f"{counter.output_dir}.h5ad")
        if counter.count_reporter_edits and match_target_pos:
            screen.uns["allele_counts"] = screen.uns["allele_counts"].loc[
                screen.uns["allele_counts"].allele.map(str) != "", :
            ]
            try:
                screen.get_edit_from_allele("allele_counts", "allele")
            except ValueError as exc:
                raise ValueError(
                    f"File {counter.output_dir}.h5ad doesn't have alllele information stored."
                ) from exc
            screen.get_edit_mat_from_uns(target_base_edits, match_target_pos)
        info(
            f"Reading already existing data for {sample_id} from \n\
            {counter.output_dir}.h5ad"
        )

    else:
        info(f"Counting {sample_id}")
        counter.check_filter_fastq()
        counter.get_counts()
        counter.screen.write(f"{counter.output_dir}.h5ad")
        screen = counter.screen
        if screen.X.max() == 0:
            warn(f"Nothing counted for {sample_id}. Check your input.")
        if counter.count_reporter_edits:
            screen.uns["allele_counts"] = screen.uns["allele_counts"].loc[
                screen.uns["allele_counts"].allele.map(str) != "", :
            ]
            screen.get_edit_from_allele("allele_counts", "allele")
            if match_target_pos:
                screen.get_edit_mat_from_uns(target_base_edits, match_target_pos)
            else:
                screen.get_edit_mat_from_uns(target_base_edits)
        info(
            f"Done for {sample_id}. \n\
            Output written at {counter.output_dir}.h5ad"
        )
    del counter
    return screen


def check_arguments(args: argparse.Namespace) -> argparse.Namespace:
    """Checks the validity of the argument."""
    args = _check_arguments(args, info, warn, error)
    sample_tbl = pd.read_csv(args.sample_list)
    if len(sample_tbl.iloc[:, 2].unique()) != len(sample_tbl.iloc[:, 2]):
        raise InputFileError(
            f"Sample ID not unique. Please check your input file {args.sample_list}."
        )
    first_read_lengths = [
        _get_first_read_length(fastq_R1) for fastq_R1 in sample_tbl.iloc[:, 0]
    ]
    for read_length in first_read_lengths:
        _check_read_length(args, read_length, warn)

    def _check_return_guide_seqs_tbl(guide_seqs_file, sample_tbl, label):
        """Check if the provided `guide_[start,end]_seqs_file` contains information about all samples in `sample_tbl`."""
        guide_seqs_tbl = pd.read_csv(guide_seqs_file, header=None, dtype=str).fillna("")
        if len(guide_seqs_tbl.columns) == 1:
            info("No guide start seq provided. Ignoring the file.")
            return None
        sample_has_seq = sample_tbl.iloc[:, 2].isin(guide_seqs_tbl[0])
        if not sample_has_seq.all():
            raise InputFileError(
                f"Sample ID(s) {sample_tbl.iloc[:,2][~sample_has_seq]} not in {label}_seqs_file {guide_seqs_tbl}"
            )
        guide_seqs_tbl.columns = ["sample", "seq"]
        return guide_seqs_tbl.set_index("sample")["seq"]

    if args.guide_start_seqs_file is not None:
        args.guide_start_seqs_tbl = _check_return_guide_seqs_tbl(
            args.guide_start_seqs_file, sample_tbl, "guide_start"
        )
    if args.guide_end_seqs_file is not None:
        args.guide_end_seqs_tbl = _check_return_guide_seqs_tbl(
            args.guide_end_seqs_file, sample_tbl, "guide_end"
        )
    if args.barcode_start_seqs_file is not None:
        args.barcode_start_seqs_tbl = _check_return_guide_seqs_tbl(
            args.barcode_start_seqs_file, sample_tbl, "barcode_start"
        )
    return args


def main(args):
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
    args = check_arguments(args)
    sample_tbl = pd.read_csv(args.sample_list)  # R1_filepath, R2_filepath, sample_name
    sample_tbl_input = sample_tbl.iloc[:, :3]
    sample_info_tbl = sample_tbl.iloc[:, 2:].set_index(sample_tbl.columns[2])
    with Pool(processes=args.threads, maxtasksperchild=1) as p:
        result = p.starmap(
            count_sample,
            [
                list(tup) + [args]
                for tup in list(sample_tbl_input.to_records(index=False))
            ],
        )
        # result = p.starmap(count_sample, sample_tbl[0], sample_tbl[1], sample_tbl[2])

    screen = bean.concat(result, axis=1)
    database_id = args.name or os.path.basename(args.sample_list).split(".")[0]
    output_path = os.path.join(
        os.path.abspath(args.output_folder), f"bean_count_{database_id}"
    )
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    try:
        screen.samples = screen.samples.join(sample_info_tbl, how="left")
    except TypeError:
        print(screen.samples)
        print(sample_info_tbl)
    # screen.guides = result[0].guides.loc[screen.guides.index, :]
    screen.write(f"{output_path}.h5ad")
    screen.to_Excel(f"{output_path}.xlsx")

    info("All Done!")
    print(
        r"""
    _ _       
  /  \ '\                       _   
  |   \  \      __ ___ _  _ _ _| |_ 
   \   \  |    / _/ _ \ || | ' \  _|
    `.__|/     \__\___/\_,_|_||_\__|
    """
    )
