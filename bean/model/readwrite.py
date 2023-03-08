from typing import Union, Sequence, Optional
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection


def get_fdr(mu_z, plot=False):
    p_dec = norm.cdf(mu_z)
    p_inc = 1 - norm.cdf(mu_z)
    _, fdr_dec = fdrcorrection(p_dec)
    _, fdr_inc = fdrcorrection(p_inc)
    fdr = np.minimum(fdr_dec, fdr_inc)
    return (fdr_dec, fdr_inc, fdr)


def write_result_table(
    target_info_df: pd.DataFrame,
    param_hist_dict,
    prefix: str = "",
    write_fitted_eff: bool = True,
    guide_index: Optional[Sequence[str]] = None,
    guide_acc: Optional[Sequence] = None,
    return_result: bool = False,
) -> Union[pd.DataFrame, None]:
    """Combine target information and scores to write result table to a csv file or return it."""
    if param_hist_dict["params"]["mu_loc"].dim() == 2:
        mu = param_hist_dict["params"]["mu_loc"].detach()[:, 0].numpy()
        mu_sd = param_hist_dict["params"]["mu_scale"].detach()[:, 0].numpy()
        sd = param_hist_dict["params"]["sd_loc"].detach().exp()[:, 0].numpy()
    elif param_hist_dict["params"]["mu_loc"].dim() == 1:
        mu = param_hist_dict["params"]["mu_loc"].detach().numpy()
        mu_sd = param_hist_dict["params"]["mu_scale"].detach().numpy()
        sd = param_hist_dict["params"]["sd_loc"].detach().exp().numpy()
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
        mu0 = param_hist_dict["negctrl"]["params"]["mu_loc"].detach().numpy()
        sd0 = param_hist_dict["negctrl"]["params"]["sd_loc"].detach().exp().numpy()
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
            a_fitted = param_hist_dict["params"]["alpha_pi"].detach().numpy()
            pi = a_fitted[..., 0] / a_fitted.sum(axis=1)
        sgRNA_df = pd.DataFrame({"edit_eff": pi}, index=guide_index)
        if guide_acc is not None:
            sgRNA_df["accessibility"] = guide_acc
        sgRNA_df.to_csv(f"{prefix}CRISPRbean_sgRNA_result.csv")

    if return_result:
        return fit_df
    fit_df.to_csv(f"{prefix}CRISPRbean_element_result.csv")
