from typing import Union, List
import numpy as np
import pandas as pd
from copy import deepcopy
from bean.framework.ReporterScreen import ReporterScreen, concat
import bean as be


def check_args(args):
    bdata = be.read_h5ad(args.bdata_path)
    if args.replicate_col not in bdata.samples.columns:
        raise ValueError(
            f"Specified --replicate-col `{args.replicate_col}` does not exist in ReporterScreen.samples.columns ({bdata.samples.columns}). Please check your input."
        )
    if args.condition_col not in bdata.samples.columns:
        raise ValueError(
            f"Specified --condition-col `{args.condition_col}` does not exist in ReporterScreen.samples.columns ({bdata.samples.columns}). Please check your input."
        )
    if (
        (
            (args.tiling is not None and args.tiling == False)
            or (args.tiling is None and not bdata.tiling)
        )
        and (
            ("target_base_changes" in bdata.uns and bdata.uns["target_base_changes"])
            or ("target_base_change" in bdata.uns and bdata.uns["target_base_change"])
        )
        and args.target_pos_col not in bdata.guides.columns
    ):
        raise ValueError(
            f"Specified --target-pos-col `{args.target_pos_col}` does not exist in ReporterScreen.guides.columns ({bdata.guides.columns}). Please check your input. (--tiling {args.tiling}, ReporterScreen.tiling: {bdata.tiling})"
        )
    if args.posctrl_col != "":
        if args.posctrl_col and args.posctrl_col not in bdata.guides.columns:
            raise ValueError(
                f"Specified --posctrl-col `{args.posctrl_col}` does not exist in ReporterScreen.guides.columns ({bdata.guides.columns}). Please check your input."
            )
        else:
            bdata.guides[args.posctrl_col] = bdata.guides[args.posctrl_col].astype(str)
    if (
        args.posctrl_col
        and args.posctrl_val not in bdata.guides[args.posctrl_col].tolist()
    ):
        raise ValueError(
            f"Specified --posctrl-val `{args.posctrl_val}` does not exist in ReporterScreen.guides[{args.posctrl_col}] ({bdata.guides[args.posctrl_col]}). Please check your input. To proceed without positive control, please provide --posctrl-col='' argument."
        )
    if args.control_condition not in bdata.samples[args.condition_col].tolist():
        raise ValueError(
            f"Specified --control-condition `{args.control_condition}` does not exist in ReporterScreen.samples[{args.condition_col}] :\n{bdata.samples[args.condition_col]}.\n Please check your input. \nFeed the condition where the editing rate would be quantified as the --control-condition argument, usually the condition with the least selection. (Closest to T0 for survival, pre-sort or bulk for sorting screens)."
        )

    lfc_conds = args.lfc_conds.split(",")
    if not len(lfc_conds) == 2:
        raise ValueError(
            f"lfc_conds must be two condition labels delimited by comma. Provided {args.lfc_conds}"
        )
    args.lfc_cond1 = lfc_conds[0]
    args.lfc_cond2 = lfc_conds[1]
    if args.lfc_cond1 not in bdata.samples[args.condition_col].tolist():
        raise ValueError(
            f"Specified --lfc-conds `{args.lfc_cond1}` does not exist in ReporterScreen.samples[{args.condition_col}]:\n{bdata.samples[args.condition_col]}. Please check your input."
        )
    if args.lfc_cond2 not in bdata.samples[args.condition_col].tolist():
        raise ValueError(
            f"Specified --lfc-conds `{args.lfc_cond2}` does not exist in ReporterScreen.samples[{args.condition_col}]:\n{bdata.samples[args.condition_col]}. Please check your input."
        )
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
    if args.reporter_length is None:
        if "reporter_length" in bdata.uns:
            args.reporter_length = bdata.uns["reporter_length"]
        else:
            args.reporter_length = 32
    if args.reporter_right_flank_length is None:
        if "reporter_right_flank_length" in bdata.uns:
            args.reporter_right_flank_length = bdata.uns["reporter_right_flank_length"]
        else:
            args.reporter_right_flank_length = 6
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
