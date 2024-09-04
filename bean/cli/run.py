#!/usr/bin/env python
import os
import sys
import logging
import warnings
from functools import partial
from copy import deepcopy
import numpy as np
import pandas as pd

import torch
import pyro
import pyro.infer
import pyro.optim
import pickle as pkl

import bean.model.model as m
from bean.model.readwrite import write_result_table
from bean.preprocessing.data_class import DATACLASS_DICT
from bean.preprocessing.utils import (
    prepare_bdata,
    _obtain_effective_edit_rate,
    _obtain_n_guides_alleles_per_variant,
    _obtain_n_cooccurring_variants,
)

import bean as be
from bean.model.run import (
    run_inference,
    _get_guide_info,
    _get_guide_to_variant_df,
    _get_guide_target_info,
    check_args,
    identify_model_guide,
    identify_negctrl_model_guide,
    _check_prior_params,
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

pyro.set_rng_seed(101)

warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    message=r".*is_categorical_dtype is deprecated and will be removed in a future version.*",
)

warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    message=r".*FutureWarning: The default of observed=False is deprecated and will be changed to True in a future version of pandas.*",
)


def main(args, return_data=False):
    print(
        r"""
    _ _       
  /  \ '\                 
  |   \  \      _ _ _  _ _ _  
   \   \  |    | '_| || | ' \ 
    `.__|/     |_|  \_,_|_||_|
    """
    )
    print("bean-run: Run model to identify targeted variants and their impact.")
    bdata = be.read_h5ad(args.bdata_path)
    args, bdata = check_args(args, bdata)
    prefix = (
        args.outdir
        + "/bean_run_result."
        + os.path.basename(args.bdata_path).rsplit(".", 1)[0]
    )
    os.makedirs(prefix, exist_ok=True)
    file_logger = logging.FileHandler(f"{prefix}/bean_run.log")
    file_logger.setLevel(logging.INFO)
    logging.getLogger().addHandler(file_logger)
    info(f"Running: {' '.join(sys.argv[:])}")
    if args.cuda:
        torch.set_default_tensor_type(torch.cuda.FloatTensor)
    else:
        torch.set_default_tensor_type(torch.FloatTensor)

    model_label, model, guide = identify_model_guide(args)
    info("Done loading data. Preprocessing...")

    # Format bdata into data structure compatible with Pyro model
    bdata = prepare_bdata(bdata, args, warn, prefix)
    negctrl_idx = np.where(
        bdata.guides[args.negctrl_col].map(lambda s: s.lower())
        == args.negctrl_col_value.lower()
    )[0]
    ndata = DATACLASS_DICT[args.selection][model_label](
        screen=bdata,
        device=args.device,
        repguide_mask=args.repguide_mask,
        sample_mask_column=args.sample_mask_col,
        accessibility_col=args.acc_col,
        accessibility_bw_path=args.acc_bw_path,
        use_const_pi=args.const_pi,
        condition_column=args.condition_col,
        time_column=args.time_col,
        control_condition=args.control_condition,
        control_can_be_selected=~args.exclude_control_condition_for_inference,
        allele_df_key=args.allele_df_key,
        control_guide_tag=args.control_guide_tag,
        target_col=args.target_col,
        shrink_alpha=args.shrink_alpha,
        popt=args.popt,
        replicate_col=args.replicate_col,
        use_bcmatch=(not args.ignore_bcmatch),
        negctrl_guide_idx=negctrl_idx,
    )
    guide_index = ndata.screen.guides.index.copy()
    assert len(guide_index) == bdata.n_obs, (len(guide_index), bdata.n_obs)
    if return_data:
        return ndata
    # Build variant dataframe
    adj_negctrl_idx = None

    if args.library_design == "variant":
        if not args.uniform_edit:
            if "edit_rate" not in ndata.screen.guides.columns:
                ndata.screen.get_edit_from_allele()
                ndata.screen.get_edit_mat_from_uns(rel_pos_is_reporter=True)
                ndata.screen.get_guide_edit_rate(
                    unsorted_condition_label=args.control_condition
                )
        target_info_df = _get_guide_target_info(
            ndata.screen, args, cols_include=[args.negctrl_col]
        )
        if args.adjust_confidence_by_negative_control:
            adj_negctrl_idx = np.where(
                target_info_df[args.negctrl_col].map(lambda s: s.lower())
                == args.negctrl_col_value.lower()
            )[0]
    else:
        if "edit_rate_norm" not in ndata.screen.guides.columns:
            ndata.screen.get_edit_from_allele()
            ndata.screen.get_edit_mat_from_uns(rel_pos_is_reporter=True)
            ndata.screen.get_guide_edit_rate(
                unsorted_condition_label=args.control_condition
            )
        if args.splice_site_path is not None:
            splice_site = pd.read_csv(args.splice_site_path).pos
        target_info_df = be.an.translate_allele.annotate_edit(
            pd.DataFrame(pd.Series(ndata.edit_index))
            .reset_index()
            .rename(columns={"index": "edit"}),
            control_tag=args.control_guide_tag,
            splice_sites=None if args.splice_site_path is None else splice_site,
        )
        _guide_idxs, _per_guide_edit_rates, eff_editing_rates = (
            _obtain_effective_edit_rate(ndata, count_thres=0)
        )
        # Summary stats
        # 1. guide_target_group: Aggregated `target_group` columns of guides that generated the variant.
        # 2. effective_edit_rate: Sum of normalized editing rate per variant.
        # 3. editing_guides: guide IDs that generated the variant.
        # 4. per_guide_editing_rates: comma-separated editing rates, ordered as the same as `editing_guides``.
        # 5. n_guides: # of guides generated the variant.
        # 6. n_coocc: # of unique co-occurring variants that appears together in allele that generated the variant.
        target_info_df["guide_target_group"] = [
            ",".join(
                np.unique(
                    bdata.guides["target_group"]
                    .values[_gidx.cpu().numpy().flatten()]
                    .tolist()
                )
            )
            for _gidx in _guide_idxs
        ]
        target_info_df["effective_edit_rate"] = eff_editing_rates.cpu()
        target_info_df["editing_guides"] = [
            ",".join(guide_index.values[_gidx.cpu().numpy().flatten()].tolist())
            for _gidx in _guide_idxs
        ]
        target_info_df["per_guide_editing_rates"] = [
            ",".join([f"{p:.3g}" for p in pger]) for pger in _per_guide_edit_rates
        ]
        target_info_df["n_guides"] = _obtain_n_guides_alleles_per_variant(ndata).cpu()
        target_info_df["n_coocc"] = _obtain_n_cooccurring_variants(ndata)
        if args.adjust_confidence_by_negative_control:
            adj_negctrl_idx = np.where(
                (target_info_df.ref == target_info_df.alt)
                & (target_info_df.coding == "coding")
            )[0]
            info(
                f"Using {len(adj_negctrl_idx)} synonymous variants to adjust confidence."
            )
    guide_info_df = _get_guide_info(
        ndata.screen, args, guide_lfc_pseudocount=args.guide_lfc_pseudocount
    )
    if args.library_design == "tiling":
        # add variant into to guide_info_tbl
        guide_info_df = pd.concat(
            [
                guide_info_df,
                _get_guide_to_variant_df(target_info_df).reindex(guide_info_df.index),
            ],
            axis=1,
        )
    # Add user-defined prior.
    if args.prior_params is not None:
        prior_params = _check_prior_params(args.prior_params, ndata)
        model = partial(model, prior_params=prior_params)

    # Run the inference steps
    info(f"Running inference for {model_label}...")
    if args.load_existing:
        with open(f"{prefix}/{model_label}.result.pkl", "rb") as handle:
            param_history_dict = pkl.load(handle)
    else:
        save_dict = dict()
        if args.fit_negctrl:
            negctrl_model, negctrl_guide = identify_negctrl_model_guide(
                args, "X_bcmatch" in bdata.layers
            )
            print(f"Using {negctrl_model} & {negctrl_guide}")
            negctrl_idx = np.where(
                ndata.screen.guides[args.negctrl_col].map(lambda s: s.lower())
                == args.negctrl_col_value.lower()
            )[0]
            info(
                f"Using {len(negctrl_idx)} negative control elements to adjust phenotypic effect sizes..."
            )
            ndata_negctrl = ndata[negctrl_idx]
            print(
                f"ndata size factor: {ndata.size_factor}, {ndata_negctrl.size_factor}\esf"
            )
            param_history_dict_negctrl, save_dict["negctrl"] = deepcopy(
                run_inference(
                    negctrl_model, negctrl_guide, ndata_negctrl, num_steps=args.n_iter
                )
            )
            if args.selection == "survival":
                print(
                    f"Feeding mu_negctrl={param_history_dict_negctrl['mu_loc']} into model..."
                )
                model = partial(model, mu_negctrl=param_history_dict_negctrl["mu_loc"])
        else:
            param_history_dict_negctrl = None
        save_dict["data"] = ndata
        param_history_dict, save_dict_model = deepcopy(
            run_inference(model, guide, ndata, num_steps=args.n_iter)
        )
        for k, v in save_dict_model.items():
            save_dict[k] = v
    # Save results
    outfile_path = (
        f"{prefix}/bean_element[sgRNA]_result.{model_label}{args.result_suffix}.csv"
    )
    info(f"Done running inference. Writing result at {outfile_path}...")
    if not os.path.exists(prefix):
        os.makedirs(prefix)
    if args.save_raw:
        with open(
            f"{prefix}/{model_label}.result{args.result_suffix}.pkl", "wb"
        ) as handle:
            pkl.dump(save_dict, handle)
    write_result_table(
        target_info_df,
        guide_info_df,
        param_history_dict,
        negctrl_params=param_history_dict_negctrl,
        model_label=model_label,
        prefix=f"{prefix}/",
        suffix=args.result_suffix,
        guide_acc=(
            ndata.guide_accessibility.cpu().numpy()
            if hasattr(ndata, "guide_accessibility")
            and ndata.guide_accessibility is not None
            else None
        ),
        adjust_confidence_by_negative_control=args.adjust_confidence_by_negative_control,
        adjust_confidence_negatives=adj_negctrl_idx,
        sd_is_fitted=(args.selection == "sorting"),
        sample_covariates=(
            ndata.sample_covariates if hasattr(ndata, "sample_covariates") else None
        ),
    )
    info("Done!")
