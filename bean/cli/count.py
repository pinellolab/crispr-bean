#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Count guides, optionally with reporter and alleles of a single sequencing sample."""

import logging
import os
import sys

import bean
from bean.mapping.utils import (
    _check_arguments,
    _check_file,
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


_ROOT = os.path.abspath(os.path.dirname(__file__))


def check_arguments(args, info_logger, warn_logger, error_logger):
    args = _check_arguments(
        args, info_logger=info, warn_logger=warn, error_logger=error
    )
    _check_file(args.R1)
    _check_file(args.R2)
    read_length = _get_first_read_length(args.R1)
    _check_read_length(args, read_length, warn)
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
    args = check_arguments(args, info_logger=info, warn_logger=warn, error_logger=error)
    args_dict = vars(args)

    base_editing_map = {"A": "G", "C": "T"}
    try:
        target_base_edits = {k: base_editing_map[k] for k in args_dict["edited_base"]}
    except KeyError as e:
        raise KeyError(args_dict["edited_base"]) from e
    match_target_pos = args_dict["match_target_pos"]

    counter = bean.mp.GuideEditCounter(**args_dict)
    counter.check_filter_fastq()

    counter.get_counts()
    if counter.count_reporter_edits:
        counter.screen.uns["allele_counts"] = counter.screen.uns["allele_counts"].loc[
            counter.screen.uns["allele_counts"].allele.map(str) != "", :
        ]
        if match_target_pos:
            counter.screen.get_edit_mat_from_uns(target_base_edits, match_target_pos)
        else:
            counter.screen.get_edit_mat_from_uns(target_base_edits)
    counter.screen.write(f"{counter.output_dir}.h5ad")
    counter.screen.to_Excel(f"{counter.output_dir}.xlsx")
    info(f"Output written at:\n {counter.output_dir}.h5ad,\n {counter.output_dir}.xlsx")
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
