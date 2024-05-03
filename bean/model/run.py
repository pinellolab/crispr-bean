import os
import sys

from tqdm import tqdm
import pickle as pkl
import pandas as pd
import logging
from functools import partial
import pyro
import bean.model.model as sorting_model
import bean.model.survival_model as survival_model

# import bean.model.survival_model as survival_model

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


def check_args(args, bdata):
    args.adjust_confidence_by_negative_control = (
        not args.no_negative_control
    )
    if args.scale_by_acc:
        if args.acc_col is None and args.acc_bw_path is None:
            raise ValueError(
                "--scale-by-acc not accompanied by --acc-col nor --acc-bw-path to use. Pass either one."
            )
        elif args.acc_col is not None and args.acc_bw_path is not None:
            warn(
                "Both --acc-col and --acc-bw-path is specified. --acc-bw-path is ignored."
            )
            args.acc_bw_path = None
    if args.outdir is None:
        args.outdir = os.path.dirname(args.bdata_path)
    if args.library_design == "variant":
        if args.adjust_confidence_by_negative_control and (args.negctrl_col not in bdata.guides.columns):
            raise ValueError(f"--negctrl-col argument '{args.negctrl_col}' not in ReporterScreen.guides.columns {bdata.guides.columns}. Please check the input or provide --no-negative-control flag if you don't have the negative controls.")
    elif args.library_design == "tiling":
        if args.allele_df_key is None:
            key_to_use = "allele_counts"
            n_alleles = len(bdata.uns[key_to_use])
            for key, df in bdata.uns.items():
                if "allele_counts" not in key or not isinstance(df, pd.DataFrame):
                    continue
                if len(df) < n_alleles:
                    key_to_use = key
                    n_alleles = len(df)
            warn(
                f"--allele-df-key not provided for tiling screen. Using the most filtered allele counts with {n_alleles} alleles stored in {key_to_use}."
            )
            args.allele_df_key = key_to_use
        elif args.allele_df_key not in bdata.uns:
            raise ValueError(
                f"--allele-df-key {args.allele_df_key} not in ReporterScreen.uns. Check your input."
            )
    else:
        raise ValueError(
            "Invalid library_design provided. Select either 'variant' or 'tiling'."
        )  # TODO: change this into discrete modes via argparse
    if args.fit_negctrl:
        n_negctrl = (
            bdata.guides[args.negctrl_col].map(lambda s: s.lower())
            == args.negctrl_col_value.lower()
        ).sum()
        if not n_negctrl >= 10:
            raise ValueError(
                f"Not enough negative control guide in the input data: {n_negctrl}. Check your input arguments."
            )
    if args.repguide_mask is not None and args.repguide_mask not in bdata.uns.keys():
        bdata.uns[args.repguide_mask] = pd.DataFrame(
            index=bdata.guides.index, columns=bdata.samples[args.replicate_col].unique()
        ).fillna(1)
        warn(
            f"{args.bdata_path} does not have replicate x guide outlier mask. All guides are included in analysis."
        )
    if args.sample_mask_col is not None:
        if args.sample_mask_col not in bdata.samples.columns.tolist():
            raise ValueError(
                f"{args.bdata_path} does not have specified sample mask column {args.sample_mask_col} in .samples"
            )
    if args.condition_col not in bdata.samples.columns:
        raise ValueError(
            f"Condition column set by `--condition-col` {args.condition_col} not in ReporterScreen.samples.columns:{bdata.samples.columns}. Check your input."
        )
    if args.control_condition not in bdata.samples[args.condition_col].tolist():
        raise ValueError(
            f"No sample has control label (set by `--control-condition`) {args.control_condition} in ReporterScreen.samples[{args.condition_col}]: {bdata.samples[args.condition_col]}. Check your input."
        )
    if args.replicate_col not in bdata.samples.columns:
        raise ValueError(
            f"Condition column set by `--replicate-col` {args.replicate_col} not in ReporterScreen.samples.columns:{bdata.samples.columns}. Check your input."
        )
    if args.control_guide_tag is not None:
        if args.library_design == "variant":
            raise ValueError(
                "`--control-guide-tag` is not used for the variant mode. Make sure you provide the separate `target` column for negative control guide that targets different negative control variant."
            )
        elif not bdata.guides.index.map(lambda s: args.control_guide_tag in s).any():
            raise ValueError(
                f"Negative control guide label {args.control_guide_tag} provided by `--control-guide-tag` doesn't appear in any of the guide names. Check your input."
            )
    if args.alpha_if_overdispersion_fitting_fails is not None:
        try:
            b0, b1 = args.alpha_if_overdispersion_fitting_fails.split(",")
            args.popt = (float(b0), float(b1))
        except TypeError as e:
            raise ValueError(
                f"Input --alpha-if-overdispersion-fitting-fails {args.alpha_if_overdispersion_fitting_fails} is malformatted! Provide [float].[float] format."
            )
    else:
        args.popt = None
    return args, bdata


def _get_guide_target_info(bdata, args, cols_include=[]):
    guide_info = bdata.guides.copy()
    target_info = (
        guide_info[
            [args.target_col]
            + [
                col
                for col in guide_info.columns
                if (
                    (
                        (col.startswith("target_"))
                        and len(guide_info[[args.target_col, col]].drop_duplicates())
                        == len(guide_info[args.target_col].unique())
                    )
                    or col in cols_include
                )
                and col != args.target_col
            ]
        ]
        .drop_duplicates()
        .set_index(args.target_col, drop=True)
    )
    target_info["n_guides"] = guide_info.groupby("target").size()

    if "edit_rate" in guide_info.columns.tolist():
        edit_rate_info = (
            guide_info[[args.target_col, "edit_rate"]]
            .groupby(args.target_col, sort=False)
            .agg({"edit_rate": ["mean", "std"]})
        )
        edit_rate_info.columns = edit_rate_info.columns.get_level_values(1)
        edit_rate_info = edit_rate_info.rename(
            columns={"mean": "edit_rate_mean", "std": "edit_rate_std"}
        )
        target_info = target_info.join(edit_rate_info)
    return target_info


def run_inference(
    model, guide, data, initial_lr=0.01, gamma=0.1, num_steps=2000, autoguide=False
):
    pyro.clear_param_store()
    lrd = gamma ** (1 / num_steps)
    svi = pyro.infer.SVI(
        model=model,
        guide=guide,
        optim=pyro.optim.ClippedAdam({"lr": initial_lr, "lrd": lrd}),
        loss=pyro.infer.Trace_ELBO(),
    )
    losses = []
    try:
        for t in tqdm(range(num_steps)):
            loss = svi.step(data)
            if t % 100 == 0:
                print(f"loss {loss} @ iter {t}")
            losses.append(loss)
    except ValueError as exc:
        error(
            "Error occurred during fitting. Saving temporary output at tmp_result.pkl."
        )
        with open("tmp_result.pkl", "wb") as handle:
            pkl.dump({"param": pyro.get_param_store()}, handle)

        raise ValueError(
            f"Fitting halted for command: {' '.join(sys.argv)} with following error: \n {exc}"
        )
    return pyro.get_param_store(), {
        "loss": losses,
        "params": pyro.get_param_store().get_state(),
    }


def identify_model_guide(args):
    if args.selection == "sorting":
        m = sorting_model
    else:
        m = survival_model
    if args.library_design == "tiling":
        info("Using Mixture Normal model...")
        return (
            f"MultiMixtureNormal{'+Acc' if args.scale_by_acc else ''}",
            partial(
                m.MultiMixtureNormalModel,
                scale_by_accessibility=args.scale_by_acc,
                use_bcmatch=(not args.ignore_bcmatch,),
            ),
            partial(
                m.MultiMixtureNormalGuide,
                scale_by_accessibility=args.scale_by_acc,
                fit_noise=~args.dont_fit_noise,
            ),
        )
    if args.uniform_edit:
        if args.guide_activity_col is not None:
            raise ValueError(
                "Can't use the guide activity column while constraining uniform edit."
            )
        info("Using Normal model...")
        return (
            "Normal",
            partial(m.NormalModel, use_bcmatch=(not args.ignore_bcmatch)),
            m.NormalGuide,
        )
    elif args.const_pi:
        if args.guide_activity_col is not None:
            raise ValueError(
                "--guide-activity-col to be used as constant pi is not provided."
            )
        info("Using Mixture Normal model with constant weight ...")
        return (
            "MixtureNormalConstPi",
            partial(m.MixtureNormalConstPiModel, use_bcmatch=(not args.ignore_bcmatch)),
            m.MixtureNormalGuide,
        )
    else:
        info(
            f"Using Mixture Normal model {'with accessibility normalization' if args.scale_by_acc else ''}..."
        )
        return (
            f"{'_' if args.dont_fit_noise else ''}MixtureNormal{'+Acc' if args.scale_by_acc else ''}",
            partial(
                m.MixtureNormalModel,
                scale_by_accessibility=args.scale_by_acc,
                use_bcmatch=(not args.ignore_bcmatch,),
            ),
            partial(
                m.MixtureNormalGuide,
                scale_by_accessibility=args.scale_by_acc,
                fit_noise=(not args.dont_fit_noise),
            ),
        )


def identify_negctrl_model_guide(args, data_has_bcmatch):
    if args.selection == "sorting":
        m = sorting_model
    else:
        m = survival_model
    negctrl_model = partial(
        m.ControlNormalModel,
        use_bcmatch=(not args.ignore_bcmatch and data_has_bcmatch),
    )

    negctrl_guide = partial(
        m.ControlNormalGuide,
        use_bcmatch=(not args.ignore_bcmatch and data_has_bcmatch),
    )
    return negctrl_model, negctrl_guide
