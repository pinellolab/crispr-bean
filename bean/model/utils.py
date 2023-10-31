import os
import sys
import argparse
from tqdm import tqdm
import logging
import pickle as pkl
import pandas as pd
import torch
import torch.distributions as tdist
import pyro
import pyro.distributions as dist
import pyro.distributions.constraints as constraints

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
        "mode",
        type=str,
        help="[variant, tiling]- Screen type whether to run variant or tiling screen model.",
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
        "--control-condition-label",
        default="bulk",
        type=str,
        help="Value in `bdata.samples[condition_col]` that indicates control experimental condition.",
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
        help="Adjust confidence by negative controls. For variant mode, this uses negative control variants. For tiling mode, adjusts confidence by synonymous edits.",
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
    if args.mode == "variant":
        pass
    elif args.mode == "tiling":
        if args.allele_df_key is None:
            raise ValueError(
                "--allele-df-key not provided for tiling screen. Feed in the key then allele counts in screen.uns[allele_df_key] will be used."
            )
    else:
        raise ValueError(
            "Invalid mode provided. Select either 'variant' or 'tiling'."
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

    return args, bdata


def get_alpha(
    expected_guide_p, size_factor, sample_mask, a0, epsilon=1e-5, normalize_by_a0=True
):
    p = (
        expected_guide_p.permute(0, 2, 1) * size_factor[:, None, :]
    )  # (n_reps, n_guides, n_bins)
    if normalize_by_a0:
        a = (
            (p + epsilon / p.shape[-1])
            / (p.sum(axis=-1)[:, :, None] + epsilon)
            * a0[None, :, None]
        )
        a = (a * sample_mask[:, None, :]).clamp(min=epsilon)
        return a
    a = (p * sample_mask[:, None, :]).clamp(min=epsilon)
    return a


def get_std_normal_prob(
    upper_quantile: torch.Tensor,
    lower_quantile: torch.Tensor,
    mu: torch.Tensor,
    sd: torch.Tensor,
    mask=None,
) -> torch.Tensor:
    """
    Returns the probability that the normal distribution with mu and sd will
    lie between upper_quantile and lower_quantile of normal distribution
    centered at 0 and has scale sd.
    Arguments
    - mask: ignore index if 0 (=False).
    """
    inf_mask = upper_quantile == 1.0
    ninf_mask = lower_quantile == 0.0

    upper_thres = torch.ones_like(upper_quantile)
    lower_thres = torch.zeros_like(lower_quantile)
    upper_thres[~inf_mask] = tdist.Normal(0, 1).icdf(upper_quantile[~inf_mask])
    lower_thres[~ninf_mask] = tdist.Normal(0, 1).icdf(lower_quantile[~ninf_mask])

    if mask is not None:
        # Temporarily add mask value
        sd = sd + (~mask).long() * 100
    assert (sd > 0).all(), f"sd > 0 not met at {sd[~(sd > 0)]}"
    cdf_upper = torch.ones_like(upper_quantile)
    cdf_lower = torch.zeros_like(lower_quantile)
    assert (
        cdf_upper.shape == inf_mask.shape
    ), f"mask:{inf_mask.shape}, cdf:{cdf_upper.shape}"

    cdf_upper[~inf_mask] = tdist.Normal(mu[~inf_mask], sd[~inf_mask]).cdf(
        upper_thres[~inf_mask]
    )
    cdf_lower[~ninf_mask] = tdist.Normal(mu[~ninf_mask], sd[~ninf_mask]).cdf(
        lower_thres[~ninf_mask]
    )
    res = cdf_upper - cdf_lower
    if mask is not None:
        res[~mask] = 0

    return res


def _scale_edited_pi(
    pi: torch.Tensor,
    guide_accessibility: torch.Tensor,
    a: float = 0.2513,
    b: float = -1.9458,
):
    """Scale editied pi by its accessibility.

    Data fitted through relationship observed from data (updated 1/17/2023).
    Transformation derived from linear regression of log(endo/reporter) ~ a*log(atac_signal + 1)+b.
    If pi of multiple alleles are provided, pi is scaled so that total scaled pi wouldn't exceed 1.0.

    Args
    pi: Editing rate
    guide_accessibility: raw accessibility score
    """
    assert pi.shape[-2] == guide_accessibility.shape[-1], (
        pi.shape,
        guide_accessibility.shape,
    )
    return (
        pi
        * torch.exp(torch.tensor(b))
        * torch.pow(guide_accessibility, a).unsqueeze(-1)
    )


def scale_pi_by_accessibility(
    pi: torch.Tensor,
    guide_accessibility: torch.Tensor,
    a: float = 0.2513,
    b: float = -1.9458,
    fit_noise: bool = False,
):
    """Scale pi by its accessibility.

    Data fitted through relationship observed from data (updated 1/17/2023).
    Transformation derived from linear regression of log(endo/reporter) ~ a*log(atac_signal + 1)+b.
    If pi of multiple alleles are provided, pi is scaled so that total scaled pi wouldn't exceed 1.0.

    Args
    pi: Editing rate (Has control editing rate at 0th index of last dimension)
    guide_accessibility: raw accessibility score for each guide. Shape (n_guides,)
    """
    pi_shape = pi.shape
    scaled_pi = _scale_edited_pi(pi[..., 1:], guide_accessibility)
    ctrl_pi = torch.ones(pi[..., 0].shape) - scaled_pi.sum(axis=-1)
    pi = torch.concat([ctrl_pi.unsqueeze(-1), scaled_pi], axis=-1)
    pi = pi / pi.sum(axis=-1).clamp(min=1.0)[..., None]
    assert pi.shape == pi_shape, (pi.shape, scaled_pi.shape, pi_shape)
    pi = add_noise_to_pi(pi, fit_noise=fit_noise)
    return pi


def add_noise_to_pi(pi: torch.Tensor, pi_noise_sd: float = 0.655, fit_noise=False):
    """Add noise to pi.

    Gaussian noise is added on the logit scale. Scale of the noise is fitted from the data (updated 1/17/2023).

    Args
    pi: Editing rate of shape (n_reps, 1, n_guides, n_alleles). Has control editing rate at 0th index of last dimension.
    """
    pi_shape = pi.shape
    n_reps, _, n_guides, n_alleles = pi_shape
    logit_pi = torch.logit(pi[..., 1:].clamp(min=1e-3, max=1 - 1e-3))
    if fit_noise:
        with pyro.plate("guide_plate_noise", n_guides):
            noise_loc = pyro.param("noise_loc", torch.zeros((n_guides,)))
            noise_scale = pyro.param(
                "noise_scale",
                torch.ones((n_guides,)) * pi_noise_sd,
                constraint=constraints.positive,
            )
            logit_pi_noise = pyro.sample(
                "logit_pi_noise",
                dist.Normal(noise_loc, noise_scale),  # infer={"is_auxiliary": True}
            )
    else:
        with pyro.plate("guide_plate_noise", n_guides):
            logit_pi_noise = pyro.sample(
                "logit_pi_noise",
                dist.Normal(0, pi_noise_sd),  # infer={"is_auxiliary": True}
            )
    logit_pi += (
        logit_pi_noise.unsqueeze(0)
        .unsqueeze(0)
        .unsqueeze(-1)
        .expand((n_reps, 1, -1, n_alleles - 1))
    )
    exp_pi_noised = torch.exp(logit_pi)
    pi_noised = (exp_pi_noised / (1 + exp_pi_noised)).clamp(min=1e-3, max=1 - 1e-3)
    pi = torch.concat(
        [
            (torch.ones(pi[:, :, :, 0].shape) - pi_noised.sum(axis=-1)).unsqueeze(-1),
            pi_noised,
        ],
        axis=-1,
    )
    assert pi.shape == pi_shape, pi.shape
    return pi
