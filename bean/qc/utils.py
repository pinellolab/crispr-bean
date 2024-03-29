import distutils
from typing import Union, List
import numpy as np
import pandas as pd
from copy import deepcopy
import argparse
from ..framework.ReporterScreen import ReporterScreen, concat


def parse_args():
    print("  \n~~~BEANQC~~~")
    print("-Check guide/sample level quality and mask / discard-")
    print(
        r"""
    _ _       
  /  \ '\        ___   ___ 
  |   \  \      / _ \ / __|
   \   \  |    | (_) | (__ 
    `.__|/      \__\_\\___|
    """
    )
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "bdata_path", help="Path to the ReporterScreen object to run QC on", type=str
    )
    parser.add_argument(
        "-o",
        "--out-screen-path",
        help="Path where quality-filtered ReporterScreen object to be written to",
        type=str,
    )
    parser.add_argument(
        "-i",
        "--ignore-missing-samples",
        help="If the flag is not provided, if the ReporterScreen object does not contain all condiitons for each replicate, make fake empty samples. If the flag is provided, don't add dummy samples.",
        action="store_true",
    )
    parser.add_argument(
        "-r",
        "--out-report-prefix",
        help="Output prefix of qc report (prefix.html, prefix.ipynb)",
        type=str,
    )
    parser.add_argument(
        "--tiling",
        dest="tiling",
        type=lambda x: bool(distutils.util.strtobool(x)),
        help="Specify that the guide library is tiling library without 'n guides per target' design",
    )
    parser.add_argument(
        "--replicate-label",
        help="Label of column in `bdata.samples` that describes replicate ID.",
        type=str,
        default="rep",
    )
    parser.add_argument(
        "--sample-covariates",
        help="Comma-separated list of column names in `bdata.samples` that describes non-selective experimental condition. (drug treatment, etc.)",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--condition-label",
        help="Label of column in `bdata.samples` that describes experimental condition. (sorting bin, time, etc.)",
        type=str,
        default="condition",
    )
    parser.add_argument(
        "--no-editing",
        help="Ignore QC about editing. Can be used for QC of other editing modalities.",
        action="store_true",
    )
    parser.add_argument(
        "--target-pos-col",
        help="Target position column in `bdata.guides` specifying target edit position in reporter",
        type=str,
        default="target_pos",
    )
    parser.add_argument(
        "--rel-pos-is-reporter",
        help="Specifies whether `edit_start_pos` and `edit_end_pos` are relative to reporter position. If `False`, those are relative to spacer position.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--edit-start-pos",
        help="Edit start position to quantify editing rate on, 0-based inclusive.",
        default=2,
    )
    parser.add_argument(
        "--edit-end-pos",
        help="Edit end position to quantify editing rate on, 0-based exclusive.",
        default=7,
    )
    parser.add_argument(
        "--count-correlation-thres",
        help="Correlation threshold to mask out.",
        type=float,
        default=0.7,
    )
    parser.add_argument(
        "--edit-rate-thres",
        help="Mean editing rate threshold per sample to mask out.",
        type=float,
        default=0.1,
    )
    parser.add_argument(
        "--posctrl-col",
        help="Column name in ReporterScreen.guides DataFrame that specifies guide category. To use all gRNAs, feed empty string ''.",
        type=str,
        default="target_group",
    )
    parser.add_argument(
        "--posctrl-val",
        help="Value in ReporterScreen.guides[`posctrl_col`] that specifies guide will be used as the positive control in calculating log fold change.",
        type=str,
        default="PosCtrl",
    )
    parser.add_argument(
        "--lfc-thres",
        help="Positive guides' correlation threshold to filter out.",
        type=float,
        default=-0.1,
    )
    parser.add_argument(
        "--lfc-conds",
        help="Values in of column in `ReporterScreen.samples[condition_label]` for LFC will be calculated between, delimited by comma",
        type=str,
        default="top,bot",
    )
    parser.add_argument(
        "--ctrl-cond",
        help="Values in of column in `ReporterScreen.samples[condition_label]` for guide-level editing rate to be calculated",
        type=str,
        default="bulk",
    )
    parser.add_argument(
        "--recalculate-edits",
        help="Even when ReporterScreen.layers['edit_count'] exists, recalculate the edit counts from ReporterScreen.uns['allele_count'].",
        action="store_true",
    )
    args = parser.parse_args()
    if args.out_screen_path is None:
        args.out_screen_path = f"{args.bdata_path.rsplit('.h5ad', 1)[0]}.filtered.h5ad"
    if args.out_report_prefix is None:
        args.out_report_prefix = f"{args.bdata_path.rsplit('.h5ad', 1)[0]}.qc_report"
    return args


def check_args(args):
    lfc_conds = args.lfc_conds.split(",")
    if not len(lfc_conds) == 2:
        raise ValueError(
            f"lfc_conds must be two condition labels delimited by comma. Provided {args.lfc_conds}"
        )
    args.lfc_cond1 = lfc_conds[0]
    args.lfc_cond2 = lfc_conds[1]
    if args.sample_covariates is not None:
        if "," in args.sample_covariates:
            args.sample_covariates = args.sample_covariates.split(",")
            args.replicate_label = [args.replicate_label] + args.sample_covariates
        else:
            args.replicate_label = [args.replicate_label, args.sample_covariates]
    if args.no_editing:
        args.base_edit_data = False
    else:
        args.base_edit_data = True
    return args


def _add_dummy_sample(
    bdata, rep, cond, condition_label: str, replicate_label: Union[str, List[str]]
):
    sample_id = f"{rep}_{cond}"
    cond_df = deepcopy(bdata.samples)
    # cond_df = cond_df.drop_duplicates()
    cond_row = cond_df.loc[cond_df[condition_label] == cond, :]
    if len(cond_row) != 1:
        cond_row = cond_row.iloc[[0], :]
    cond_row.index = [sample_id]
    cond_row[replicate_label] = rep
    dummy_sample_bdata = ReporterScreen(
        X=np.zeros((bdata.n_obs, 1)),
        X_bcmatch=np.zeros((bdata.n_obs, 1)),
        guides=bdata.guides,
        samples=cond_row,
    )
    for k in bdata.uns.keys():
        if isinstance(bdata.uns[k], pd.DataFrame):
            dummy_sample_bdata.uns[k] = pd.DataFrame(
                columns=bdata.uns[k].columns.tolist()[:2] + [sample_id]
            )
        else:
            dummy_sample_bdata.uns[k] = bdata.uns[k]
    bdata = concat([bdata, dummy_sample_bdata])
    return bdata


def fill_in_missing_samples(
    bdata, condition_label: str, replicate_label: Union[str, List[str]]
):
    """If not all condition exists for every replicate in bdata, fill in fake sample"""
    added_dummy = False
    if isinstance(replicate_label, str):
        rep_list = bdata.samples[replicate_label].unique()
    else:
        rep_list = (
            bdata.samples[replicate_label].drop_duplicates().to_records(index=False)
        )
        # print(rep_list)
    for rep in rep_list:
        for cond in bdata.samples[condition_label].unique():
            if isinstance(replicate_label, str):
                rep_samples = bdata.samples[replicate_label] == rep
            else:
                rep = list(rep)
                rep_samples = (bdata.samples[replicate_label] == rep).all(axis=1)
            if (
                len(np.where(rep_samples & (bdata.samples[condition_label] == cond))[0])
                != 1
            ):
                print(f"Adding dummy samples for {rep}, {cond}")
                bdata = _add_dummy_sample(
                    bdata, rep, cond, condition_label, replicate_label
                )
                if not added_dummy:
                    added_dummy = True
    if added_dummy:
        if isinstance(replicate_label, str):
            sort_labels = [replicate_label, condition_label]
        else:
            sort_labels = replicate_label + [condition_label]
            bdata = bdata[:, bdata.samples.sort_values(sort_labels).index]

    return bdata
