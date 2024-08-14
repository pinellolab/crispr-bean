from typing import Union, Sequence, Optional, List
import numpy as np
import pandas as pd
from statistics import NormalDist
from scipy.special import logit, expit
from scipy.stats import norm


def get_novl(df, mu_col, mu_sd_col):
    return df.apply(
        lambda row: 1
        - NormalDist(mu=row[mu_col], sigma=row[mu_sd_col]).overlap(
            NormalDist(mu=0, sigma=1)
        ),
        axis=1,
    )


def get_quantile(mu, sd, q):
    dist = norm(mu, sd)
    return dist.ppf(q)


def add_credible_interval(df, mu_col, mu_sd_col, alpha=0.05):
    """Add credible interval to df"""
    df = df.copy()
    df[f"CI[{alpha/2}"] = get_quantile(df[mu_col], df[mu_sd_col], alpha / 2)
    df[f"{1-alpha/2}]"] = get_quantile(df[mu_col], df[mu_sd_col], 1 - alpha / 2)
    return df


def adjust_normal_params_by_control(
    param_df: pd.DataFrame,
    sd0: float,
    suffix: str = "_adj",
    mu_adjusted_col="mu",
    mu_sd_adjusted_col="mu_sd",
    mu0: float = 0.0,
) -> pd.DataFrame:
    """Adjust the z-score by scaling by the standard deviation of negative variants sd0."""
    param_df[f"mu{suffix}"] = param_df[mu_adjusted_col] - mu0
    param_df[f"mu_sd{suffix}"] = param_df[mu_sd_adjusted_col] * sd0
    param_df[f"mu_z{suffix}"] = param_df[f"mu{suffix}"] / param_df[f"mu_sd{suffix}"]
    param_df[f"novl{suffix}"] = get_novl(param_df, f"mu{suffix}", f"mu_sd{suffix}")

    return param_df


def write_result_table(
    target_info_df: pd.DataFrame,
    guide_info_df: pd.DataFrame,
    param_hist_dict,
    model_label: str,
    prefix: str = "",
    suffix: str = "",
    negctrl_params=None,
    adjust_confidence_by_negative_control: bool = True,
    adjust_confidence_negatives: Optional[np.ndarray] = None,
    guide_acc: Optional[Sequence] = None,
    sd_is_fitted: bool = True,
    sample_covariates: Optional[List[str]] = None,
    return_result: bool = False,
) -> Union[pd.DataFrame, None]:
    """Combine target information and scores to write result table to a csv file or return it."""
    if param_hist_dict["mu_loc"].dim() == 2:
        mu = param_hist_dict["mu_loc"].detach()[:, 0].cpu().numpy()
        mu_sd = param_hist_dict["mu_scale"].detach()[:, 0].cpu().numpy()
        if sd_is_fitted:
            sd = param_hist_dict["sd_loc"].detach().exp()[:, 0].cpu().numpy()
    elif param_hist_dict["mu_loc"].dim() == 1:
        mu = param_hist_dict["mu_loc"].detach().cpu().numpy()
        mu_sd = param_hist_dict["mu_scale"].detach().cpu().numpy()
        if sd_is_fitted:
            sd = param_hist_dict["sd_loc"].detach().exp().cpu().numpy()
    else:
        raise ValueError(
            f'`mu_loc` has invalid shape of {param_hist_dict["mu_loc"].shape}'
        )
    param_dict = {
        "mu": mu,
        "mu_sd": mu_sd,
        "mu_z": mu / mu_sd,
    }
    if sd_is_fitted:
        param_dict["sd"] = sd
    if sample_covariates is not None:
        assert (
            "mu_cov_loc" in param_hist_dict and "mu_cov_scale" in param_hist_dict
        ), param_hist_dict.keys()
        for i, sample_cov in enumerate(sample_covariates):
            param_dict[f"mu_{sample_cov}"] = (
                mu + param_hist_dict["mu_cov_loc"].detach().cpu().numpy()[i]
            )
            param_dict[f"mu_sd_{sample_cov}"] = np.sqrt(
                mu_sd**2
                + param_hist_dict["mu_cov_scale"].detach().cpu().numpy()[i] ** 2
            )
            param_dict[f"mu_z_{sample_cov}"] = (
                param_dict[f"mu_{sample_cov}"] / param_dict[f"mu_sd_{sample_cov}"]
            )

    fit_df = pd.DataFrame(param_dict)
    if negctrl_params is not None:
        print("Normalizing with common negative control distribution")
        mu0 = negctrl_params["mu_loc"].detach().cpu().numpy()
        if sd_is_fitted:
            sd0 = negctrl_params["sd_loc"].detach().exp().cpu().numpy()
        else:
            sd0 = 1.0
        print(f"Fitted mu0={mu0}" + (f", sd0={sd0}." if sd_is_fitted else ""))
        fit_df["mu_scaled"] = (mu - mu0) / sd0
        fit_df["mu_sd_scaled"] = mu_sd / sd0
        fit_df["mu_z_scaled"] = fit_df.mu_scaled / fit_df.mu_sd_scaled
        if sd_is_fitted:
            fit_df["sd_scaled"] = sd / sd0
        fit_df["novl_scaled"] = get_novl(fit_df, "mu_scaled", "mu_sd_scaled")
        if sample_covariates is not None:
            for i, sample_cov in enumerate(sample_covariates):
                fit_df[f"mu_{sample_cov}_scaled"] = (
                    fit_df[f"mu_{sample_cov}"] - mu0
                ) / sd0
                fit_df[f"mu_sd_{sample_cov}_scaled"] = (
                    fit_df[f"mu_sd_{sample_cov}"] / sd0
                )
                fit_df[f"mu_z_{sample_cov}_scaled"] = (
                    fit_df[f"mu_{sample_cov}_scaled"] / fit_df["mu_sd_scaled"]
                )

    fit_df = pd.concat(
        [target_info_df.reset_index(), fit_df.reset_index(drop=True)], axis=1
    )

    if adjust_confidence_by_negative_control:
        assert adjust_confidence_negatives is not None
        if len(adjust_confidence_negatives) < 10:
            print(
                f"Cannot adjust confidence by negative control due to too small number ({len(adjust_confidence_negatives)}) of negatives."
            )
            fit_df = add_credible_interval(fit_df, "mu", "mu_sd")
            fit_df = fit_df.iloc[(-fit_df.mu_z.abs()).argsort()]
        else:
            ncvar = fit_df.iloc[adjust_confidence_negatives]
            if "mu_z_scaled" in ncvar.columns:
                print("Using mu_z_scaled for normalization input..")
                _, std = norm.fit(ncvar.mu_z_scaled, floc=0)
            else:
                _, std = norm.fit(ncvar.mu_z, floc=0)
            fit_df = adjust_normal_params_by_control(
                fit_df,
                std,
                suffix="_adj",
                mu_adjusted_col=(
                    "mu_scaled" if "negctrl" in param_hist_dict.keys() else "mu"
                ),
                mu_sd_adjusted_col=(
                    "mu_sd_scaled" if "negctrl" in param_hist_dict.keys() else "mu_sd"
                ),
            )
            fit_df = add_credible_interval(fit_df, "mu_adj", "mu_sd_adj")
            fit_df = fit_df.iloc[(-fit_df.mu_z_adj.abs()).argsort()]
            if sample_covariates is not None:
                for i, sample_cov in enumerate(sample_covariates):
                    fit_df = adjust_normal_params_by_control(
                        fit_df,
                        std,
                        suffix=f"_{sample_cov}_adj",
                        mu_adjusted_col=(
                            f"mu_{sample_cov}_scaled"
                            if "negctrl" in param_hist_dict.keys()
                            else f"mu_{sample_cov}"
                        ),
                        mu_sd_adjusted_col=(
                            f"mu_sd_{sample_cov}_scaled"
                            if "negctrl" in param_hist_dict.keys()
                            else f"mu_sd_{sample_cov}"
                        ),
                    )
                    fit_df = add_credible_interval(
                        fit_df, f"mu_{sample_cov}_adj", f"mu_sd_{sample_cov}_adj"
                    )
    else:
        fit_df = add_credible_interval(fit_df, "mu", "mu_sd")
        fit_df = fit_df.iloc[(-fit_df.mu_z.abs()).argsort()]
    # write sgRNA info
    if "alpha_pi" not in param_hist_dict.keys():
        pi = 1.0
    else:
        a_fitted = param_hist_dict["alpha_pi"].detach().cpu().numpy()
        pi = a_fitted[..., 1:].sum(axis=1) / a_fitted.sum(axis=1)
    # guide_info_df["edit_rate_fitted"] = pi
    guide_info_df = guide_info_df[
        ["edit_rate"]
        + [col for col in guide_info_df.columns if col != "Mid" and col != "edit_rate"]
    ]
    if guide_acc is not None:
        guide_info_df.insert(1, "accessibility", guide_acc)
        scaled_pi = _scale_pi(
            pi,
            guide_acc,
            fitted_noise_logit=(
                param_hist_dict["noise_scale"].detach().cpu().numpy()
                if "noise_scale" in param_hist_dict.keys()
                else None
            ),
        )
        guide_info_df.insert(2, "scaled_edit_eff", scaled_pi)

    guide_info_df.to_csv(f"{prefix}bean_sgRNA_result.{model_label}{suffix}.csv")
    if return_result:
        return fit_df
    fit_df.to_csv(f"{prefix}bean_element_result.{model_label}{suffix}.csv")


def _scale_edited_pi(
    pi: np.ndarray,
    guide_accessibility: np.ndarray,
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
    return pi * np.exp(b) * guide_accessibility**a


def _add_noise_to_pi(pi: np.ndarray, fitted_noise_logit: np.ndarray):
    logit_pi = logit(pi.clip(min=1e-3, max=1 - 1e-3)) + fitted_noise_logit
    return expit(logit_pi).clip(min=1e-3, max=1 - 1e-3)


def _scale_pi(pi: np.ndarray, guide_acc: np.ndarray, fitted_noise_logit=None):
    scaled_pi = _scale_edited_pi(pi, guide_acc)
    if fitted_noise_logit is None:
        return scaled_pi
    return _add_noise_to_pi(scaled_pi, fitted_noise_logit)
