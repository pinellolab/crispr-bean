from typing import Union, List
import numpy as np
import pandas as pd
from copy import deepcopy
from bean.framework.ReporterScreen import ReporterScreen, concat


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
