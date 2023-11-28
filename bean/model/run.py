import os
import sys
import argparse
from tqdm import tqdm
import pickle as pkl
import pandas as pd
import logging
from functools import partial
import pyro
import bean.model.model as sorting_model

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


def none_or_str(value):
    if value == "None":
        return None
    return value


def parse_args():
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
    parser = argparse.ArgumentParser(description="Run model on data.")
    parser.add_argument(
        "selection",
        type=str,
        choices=["sorting", "survival"],
        help="Screen selection type whether cells are sorted based on continuous phenotype ('sorting') or proliferated based on their viability ('survival').",
    )
    parser.add_argument(
        "library_design",
        type=str,
        choices=["variant", "tiling"],
        help="Library design type whether to run variant or tiling screen model.\nVariant library design assumes gRNA has specific target variant and bystander edits are ignored. Tiling library design considers all alleles generated by gRNA in reporter.",
    )
    parser.add_argument("bdata_path", type=str, help="Path of an ReporterScreen object")
    parser.add_argument(
        "--rep-pi",
        "-r",
        action="store_true",
        default=False,
        help="Fit replicate specific scaling factor. Recommended to set as True if you expect variable editing activity across biological replicates.",
    )
    parser.add_argument(
        "--uniform-edit",
        "-p",
        action="store_true",
        default=False,
        help="Assume uniform editing rate for all guides.",
    )
    parser.add_argument(
        "--scale-by-acc",
        action="store_true",
        default=False,
        help="Scale guide editing efficiency by the target loci accessibility",
    )
    parser.add_argument(
        "--acc-bw-path",
        type=str,
        default=None,
        help="Accessibility .bigWig file to be used to assign accessibility of guides.",
    )
    parser.add_argument(
        "--acc-col",
        type=str,
        default=None,
        help="Column name in bdata.guides that specify raw ATAC-seq signal.",
    )
    parser.add_argument(
        "--const-pi",
        default=False,
        action="store_true",
        help="Use constant pi provided in --guide-activity-col (instead of fitting from reporter data)",
    )
    parser.add_argument(
        "--shrink-alpha",
        default=False,
        action="store_true",
        help="Instead of using the trend-fitted alpha values, use estimated alpha values for each gRNA that are shrunk towards the fitted trend.",
    )
    parser.add_argument(
        "--condition-col",
        default="bin",
        type=str,
        help="Column key in `bdata.samples` that describes experimental condition.",
    )
    parser.add_argument(
        "--time-col",
        default="time",
        type=str,
        help="Column key in `bdata.samples` that describes time elapsed.",
    )
    parser.add_argument(
        "--control-condition-label",
        default="bulk",
        type=str,
        help="Value in `bdata.samples[condition_col]` that indicates control experimental condition.",
    )
    parser.add_argument(
        "--include-control-condition-for-inference",
        "-ic",
        default=False,
        action="store_true",
        help="Include control conditions for inference. Currently only supported for survival screens.",
    )
    parser.add_argument(
        "--replicate-col",
        default="rep",
        type=str,
        help="Column key in `bdata.samples` that describes experimental replicates.",
    )
    parser.add_argument(
        "--target-col",
        default="target",
        type=str,
        help="Column key in `bdata.guides` that describes the target element of each guide.",
    )
    parser.add_argument(
        "--guide-activity-col",
        "-a",
        type=str,
        default=None,
        help="Column in ReporterScreen.guides DataFrame showing the editing rate estimated via external tools",
    )
    parser.add_argument(
        "--outdir",
        "-o",
        default=".",
        type=str,
        help="Directory to save the run result.",
    )
    parser.add_argument(
        "--result-suffix",
        default="",
        type=str,
        help="Suffix of the output files",
    )
    parser.add_argument(
        "--sorting-bin-upper-quantile-col",
        "-uq",
        help="Column name with upper quantile values of each sorting bin in [Reporter]Screen.samples (or AnnData.var)",
        default="upper_quantile",
    )
    parser.add_argument(
        "--sorting-bin-lower-quantile-col",
        "-lq",
        help="Column name with lower quantile values of each sorting bin in [Reporter]Screen.samples (or AnnData var)",
        default="lower_quantile",
    )
    parser.add_argument(
        "--alpha-if-overdispersion-fitting-fails",
        "-af",
        default=None,
        type=str,
        help="Comma-separated regression coefficient (b0, b1) of log(a0) ~ log(q) that will be used if fitting dispersion on the data fails.",
    )
    parser.add_argument("--cuda", action="store_true", default=False, help="run on GPU")
    parser.add_argument(
        "--sample-mask-col",
        type=str,
        default=None,
        help="Name of the column indicating the sample mask in [Reporter]Screen.samples (or AnnData.var). Sample is ignored if the value in this column is 0. This can be used to mask out low-quality samples.",
    )
    parser.add_argument(
        "--fit-negctrl",
        action="store_true",
        default=False,
        help="Fit the shared negative control distribution to normalize the fitted parameters",
    )
    parser.add_argument(
        "--negctrl-col",
        type=str,
        default="target_group",
        help="Column in bdata.obs specifying if a guide is negative control. If the `bdata.guides[negctrl_col].lower() == negctrl_col_value`, it is treated as negative control guide.",
    )
    parser.add_argument(
        "--negctrl-col-value",
        type=str,
        default="negctrl",
        help="Column value in bdata.guides specifying if a guide is negative control. If the `bdata.guides[negctrl_col].lower() == negctrl_col_value`, it is treated as negative control guide.",
    )
    parser.add_argument(
        "--repguide-mask",
        type=none_or_str,
        default="repguide_mask",
        help="n_replicate x n_guide mask to mask the outlier guides. screen.uns[repguide_mask] will be used.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="Optionally use GPU if provided valid GPU device name (ex. cuda:0)",
    )
    parser.add_argument(
        "--ignore-bcmatch",
        action="store_true",
        default=False,
        help="If provided, even if the screen object has .X_bcmatch, ignore the count when fitting.",
    )
    parser.add_argument(
        "--allele-df-key",
        type=str,
        default=None,
        help="screen.uns[allele_df_key] will be used as the allele count.",
    )
    parser.add_argument(
        "--splice-site-path",
        type=str,
        default=None,
        help="Path to splicing site",
    )
    parser.add_argument(
        "--control-guide-tag",
        type=none_or_str,
        default="CONTROL",
        help="If this string is in guide name, treat each guide separately not to mix the position. Used for negative controls.",
    )
    parser.add_argument(
        "--dont-fit-noise",  # TODO: add check args
        action="store_true",
    )
    parser.add_argument(
        "--dont-adjust-confidence-by-negative-control",
        action="store_true",
        help="Adjust confidence by negative controls. For variant library_design, this uses negative control variants. For tiling library_design, adjusts confidence by synonymous edits.",
    )
    parser.add_argument(
        "--n-iter",  # TODO: add check args
        type=int,
        default=2000,
        help="# of SVI steps taken for inference.",
    )
    parser.add_argument(
        "--load-existing",  # TODO: add check args
        action="store_true",
        help="Load existing .pkl file if present.",
    )

    return parser.parse_args()


def check_args(args, bdata):
    args.adjust_confidence_by_negative_control = (
        not args.dont_adjust_confidence_by_negative_control
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
        pass
    elif args.library_design == "tiling":
        if args.allele_df_key is None:
            raise ValueError(
                "--allele-df-key not provided for tiling screen. Feed in the key then allele counts in screen.uns[allele_df_key] will be used."
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
        if not n_negctrl >= 20:
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
    if args.alpha_if_overdispersion_fitting_fails is not None:
        try:
            b0, b1 = args.alpha_if_overdispersion_fitting_fails.split(",")
            args.popt = (float(b0), float(b1))
        except TypeError as e:
            raise e(
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
    return {
        "loss": losses,
        "params": pyro.get_param_store(),
    }


def identify_model_guide(args):
    if args.selection == "sorting":
        m = sorting_model
    # else:
    #     m = survival_model
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