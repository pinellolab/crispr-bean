import numpy as np
import torch
from scipy.optimize import curve_fit

"""
Returns fitted alpha 0 (sum of concentration) of Dirichlet-Multinomial
For the notation, see DESeq paper
"""


def get_valid_idx(x: torch.Tensor, y: torch.Tensor):
    return ~(torch.isnan(x) | torch.isnan(y) | torch.isinf(x) | torch.isinf(y))


def get_valid_vals(x, y):
    valid = get_valid_idx(x, y)
    return (x[valid].detach().cpu().numpy(), y[valid].detach().cpu().numpy())


def get_q(X, sample_size_factors, sample_mask=None):
    """
    Obtain depth-normalized sample mean
    """
    n_reps, n_condits, n_guides, n_alleles = X.shape
    if sample_mask is None:
        sample_mask = torch.ones((n_reps, n_condits))
    q_condit_guide = (
        (X / sample_size_factors[:, :, None, None]) * sample_mask[:, :, None, None]
    ).sum(axis=0) / sample_mask.sum(axis=0)[:, None, None]
    assert q_condit_guide.shape == (
        n_condits,
        n_guides,
        n_alleles,
    ), q_condit_guide.shape
    return q_condit_guide[0, ...]


def get_w(X, sample_size_factors, sample_mask=None):
    """
    Obtain depth-normalized sample variance
    """
    n_reps, n_condits, n_guides, n_alleles = X.shape
    if sample_mask is None:
        sample_mask = torch.ones((n_reps, n_condits), device="cpu")
    q = get_q(X, sample_size_factors, sample_mask)
    se = (X / sample_size_factors[:, :, None, None] - q) ** 2
    w = (se * sample_mask[:, :, None, None]).sum(axis=0) / (sample_mask.sum(axis=0))[
        :, None, None
    ]
    assert w.shape == (n_condits, n_guides, n_alleles), w.shape
    return w[0, ...]


def get_a0(q, w):
    n = q.sum(-1)
    p = q / q.sum(axis=-1)[:, None]
    r = (w - q) / (n[:, None] * p * (1 - p))
    a0 = torch.nanmean(((n[:, None] - 1) / (r - 1 + 1 / (1 - p)) - 1), axis=-1)
    return a0


def linear(x, b0, b1):
    return b0 + b1 * x


def estimate_variance(y, y_est):
    return (y - y_est) ** 2 / len(y)


def shrink_normal_normal(y, y_est, shrink_prior_var: float = 1.0):
    var = estimate_variance(y, y_est)
    prior_weight = var / (var + shrink_prior_var)
    data_weight = 1 - prior_weight
    post_est = prior_weight * y_est + data_weight * y
    return post_est


def get_fitted_alpha0(
    X: torch.Tensor,
    sample_size_factors: torch.Tensor,
    sample_mask: torch.Tensor = None,
    fit: bool = True,
    fit_quantile: float = None,
    shrink: bool = False,
    shrink_prior_var: float = 1.0,
) -> (
    torch.Tensor
):  # TODO: set up for variable bin size of X? or scenario when alpha is strained to 0
    """Fits sum of concentration of DirichletMultinomial distribution.

    Args:
        X: allele counts tensor with replicate in the first dimension and allele in the last dimension.
        fit_quantile: if not None, alpha is fitted conservatively with lowest `fit_quantile`
            guides.
        shrink: Return shrinked alpha based on fitted trend, instead of fitted value itself.
    """
    n_reps, n_condits, n_guides, n_alleles = X.shape
    if sample_mask is None:
        sample_mask = torch.ones((n_reps, n_condits), device="cpu")
    w = get_w(X + 1, sample_size_factors, sample_mask=sample_mask)
    q = get_q(X + 1, sample_size_factors, sample_mask=sample_mask)
    n = q.sum(-1)
    a0 = get_a0(q, w)
    assert a0.shape == (n_guides,), a0.shape
    if not fit:
        return (a0, np.nan)

    x, y = get_valid_vals(n.log(), a0.log())
    popt, pcov = curve_fit(linear, x, y)
    print("Linear fit of log(pi_a0) ~ log(q): [b0, b1]={}, cov={}".format(popt, pcov))
    if not fit_quantile is None:
        print(f"Using lowest {fit_quantile} residual to fit again")
        sel_idx = np.where(
            (y - linear(x, *popt)) < np.quantile((y - linear(x, *popt)), fit_quantile)
        )[0]
        if len(sel_idx) < 5:
            print(
                "Too few points available for conservative re-fitting. Consider increasing `fit_quantile`."
            )
        else:
            popt, pcov = curve_fit(linear, x[sel_idx], y[sel_idx])
            print(
                "Adjusted conservative linear fit of log(pi_a0) ~ log(q): [b0, b1]={}, cov={}".format(
                    popt, pcov
                )
            )
            # TODO: Check if b1 differs a lot: if so, retrain with the same b1
    log_a0_est = torch.as_tensor(linear(n.log().cpu().numpy(), *popt))
    if shrink:
        y = torch.as_tensor(a0.log())
        y[torch.isnan(y)] = log_a0_est[torch.isnan(y)]
        log_a0_est = shrink_normal_normal(
            y, log_a0_est, shrink_prior_var=shrink_prior_var
        )
    return (torch.exp(log_a0_est), popt)


def get_pred_alpha0(X, sample_size_factors, popt, sample_mask=None):
    n_reps, n_condits, n_guides, n_alleles = X.shape
    if sample_mask is None:
        sample_mask = torch.ones((n_reps, n_condits), device="cpu")
    q = get_q(X + 1, sample_size_factors, sample_mask=sample_mask)
    a0_est = np.exp(linear(q.sum(axis=-1).log().cpu().numpy(), *popt))
    return a0_est
