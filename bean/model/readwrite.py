from typing import Union, Sequence, Optional
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.special import logit, expit
from statsmodels.stats.multitest import fdrcorrection


def get_fdr(mu_z, plot=False):
    p_dec = norm.cdf(mu_z)
    p_inc = norm.cdf(-mu_z)
    _, fdr_dec = fdrcorrection(p_dec)
    _, fdr_inc = fdrcorrection(p_inc)
    fdr = np.minimum(fdr_dec, fdr_inc)
    return (fdr_dec, fdr_inc, fdr)


def write_result_table(
    target_info_df: pd.DataFrame,
    param_hist_dict,
    model_label: str,
    prefix: str = "",
    suffix: str = "",
    write_fitted_eff: bool = True,
    guide_index: Optional[Sequence[str]] = None,
    guide_acc: Optional[Sequence] = None,
    return_result: bool = False,
) -> Union[pd.DataFrame, None]:
    """Combine target information and scores to write result table to a csv file or return it."""
    if param_hist_dict["params"]["mu_loc"].dim() == 2:
        mu = param_hist_dict["params"]["mu_loc"].detach()[:, 0].cpu().numpy()
        mu_sd = param_hist_dict["params"]["mu_scale"].detach()[:, 0].cpu().numpy()
        sd = param_hist_dict["params"]["sd_loc"].detach().exp()[:, 0].cpu().numpy()
    elif param_hist_dict["params"]["mu_loc"].dim() == 1:
        mu = param_hist_dict["params"]["mu_loc"].detach().cpu().numpy()
        mu_sd = param_hist_dict["params"]["mu_scale"].detach().cpu().numpy()
        sd = param_hist_dict["params"]["sd_loc"].detach().exp().cpu().numpy()
    else:
        raise ValueError(
            f'`mu_loc` has invalid shape of {param_hist_dict["params"]["mu_loc"].shape}'
        )
    mu_z = mu / mu_sd
    fdr_dec, fdr_inc, fdr = get_fdr(mu_z)
    fit_df = pd.DataFrame(
        {
            "mu": mu,
            "mu_sd": mu_sd,
            "mu_z": mu / mu_sd,
            "sd": sd,
            "fdr_dec": fdr_dec,
            "fdr_inc": fdr_inc,
            "fdr": fdr,
        }
    )
    if "negctrl" in param_hist_dict.keys():
        print("Normalizing with common negative control distribution")
        mu0 = param_hist_dict["negctrl"]["params"]["mu_loc"].detach().cpu().numpy()
        sd0 = (
            param_hist_dict["negctrl"]["params"]["sd_loc"].detach().exp().cpu().numpy()
        )
        print(f"Fitted mu0={mu0}, sd0={sd0}.")
        fit_df["mu_adj"] = (mu - mu0) / sd0
        fit_df["mu_sd_adj"] = mu_sd / sd0
        fit_df["mu_z_adj"] = fit_df.mu_adj / fit_df.mu_sd_adj
        fit_df["sd_adj"] = sd / sd0
        fdr_dec, fdr_inc, fdr = get_fdr(fit_df.mu_z_adj)
        fit_df["fdr_dec_adj"], fit_df["fdr_inc_adj"], fit_df["fdr_adj"] = (
            fdr_dec,
            fdr_inc,
            fdr,
        )
    fit_df = pd.concat(
        [target_info_df.reset_index(), fit_df.reset_index(drop=True)], axis=1
    )

    if write_fitted_eff or guide_acc is not None:
        if "alpha_pi" not in param_hist_dict["params"].keys():
            pi = 1.0
        else:
            a_fitted = param_hist_dict["params"]["alpha_pi"].detach().cpu().numpy()
            pi = a_fitted[..., 1:].sum(axis=1) / a_fitted.sum(axis=1)
        sgRNA_df = pd.DataFrame({"edit_eff": pi}, index=guide_index)
        if guide_acc is not None:
            sgRNA_df["accessibility"] = guide_acc
            sgRNA_df["scaled_edit_eff"] = _scale_pi(
                pi,
                guide_acc,
                fitted_noise_logit=param_hist_dict["params"]["noise_scale"]
                .detach()
                .cpu()
                .numpy(),
            )
        sgRNA_df.to_csv(f"{prefix}bean_sgRNA_result.{model_label}{suffix}.csv")

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
        return scaled_pi, None
    return _add_noise_to_pi(scaled_pi, fitted_noise_logit)
