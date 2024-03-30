#!/usr/bin/env python
"""Create Screen object (AnnData) from gRNA metadata, sample metadata, and count tables."""

import os
import sys
import logging
from bean.framework.read_from_csvs import create_screen

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


def main(args):
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
    screen = create_screen(args)
    info(f"Done obtaining screen:\n{screen}\nWriting result...")
    output_path = f"{args.output_prefix if args.output_prefix else os.path.splitext(args.gRNA_counts_table_csv)[0]}.h5ad"
    screen.write(output_path)
    info(f"Done writing screen object to {output_path}.")
