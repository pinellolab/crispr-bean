from typing import Optional, Tuple
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
    n_reps, n_condits, n_guides = X.shape
    if sample_mask is None:
        sample_mask = torch.ones((n_reps, n_condits), device="cpu")
    else:
        sample_size_factors[
            torch.logical_and(sample_mask == 0, sample_size_factors == 0)
        ] = 1.0
    q_condit_guide = (
        (X / sample_size_factors[:, :, None]) * sample_mask[:, :, None]
    ).sum(axis=0) / sample_mask.sum(axis=0)[:, None]
    assert q_condit_guide.shape == (
        n_condits,
        n_guides,
    ), f"shape mismatch: {(n_condits, n_guides)}"
    return q_condit_guide


def get_w(X, sample_size_factors, sample_mask=None):
    """
    Obtain depth-normalized sample variance
    """
    n_reps, n_condits, n_guides = X.shape
    if sample_mask is None:
        sample_mask = torch.ones((n_reps, n_condits), device="cpu")
    q = get_q(X, sample_size_factors, sample_mask)
    se = (X / sample_size_factors[:, :, None] - q) ** 2
    w = (se * sample_mask[:, :, None]).sum(axis=0) / (sample_mask.sum(axis=0))[:, None]
    return w


def linear(x, b0, b1):
    return b0 + b1 * x


def estimate_variance(y, y_est):
    return ((y - y_est) ** 2).sum() / (len(y) - 1)


def shrink_normal_normal(y, y_est, shrink_prior_var: float = 1.0):
    var = estimate_variance(y, y_est)
    prior_weight = var / (var + shrink_prior_var)
    post_est = prior_weight * y_est + (1 - prior_weight) * y
    return post_est


def get_fitted_alpha0(
    X,
    sample_size_factors,
    sample_mask=None,
    fit_quantile: float = None,
    shrink: bool = False,
    shrink_prior_var: float = 1.0,
    popt: Optional[Tuple[float, float]] = None,
):
    """Fits sum of concentration of DirichletMultinomial distribution.

    Args:
        fit: if False, return the raw value
        fit_quantile: if not None, alpha is fitted conservatively with lowest `fit_quantile`
            guides.
        popt: Regression coefficient (b0, b1) of log(a0) ~ log(q) that will be used if fitting dispersion on the data fails
    """
    n_reps, n_condits, n_guides = X.shape
    if sample_mask is None:
        sample_mask = torch.ones((n_reps, n_condits), device="cpu")
    elif (sample_mask.sum(axis=0) == 0).any():
        raise ValueError("Some bins have no data.")
    print(sample_size_factors)
    w = get_w(X + 1, sample_size_factors, sample_mask=sample_mask)
    q = get_q(X + 1, sample_size_factors, sample_mask=sample_mask)
    n = (
        torch.nanmean(q, axis=0) * q.shape[0]
    )  # depth-normalized total counts across bins
    p = q / n[None, :]  # depth-normalized p for Multinomial
    multinom_var = n[None, :] * p * (1 - p)  # theoretical Multinomial variance
    r = (w - q) / multinom_var
    a0 = ((n - 1) / (r - 1 + 1 / (1 - p)) - 1).mean(axis=0)

    x, y = get_valid_vals(n.log(), a0.log())
    if len(y) < 5:
        if popt is None:
            popt = (-1.510, 0.7861)
        print(
            f"Cannot fit log(a0) ~ log(q): data too sparse ({len(y)} valid values)! Using pre-fitted values [b0, b1]={popt}"
        )
        print(n)
        print(a0)
    else:
        popt, pcov = curve_fit(linear, x, y)
        print("Linear fit of log(a0) ~ log(q): [b0, b1]={}, cov={}".format(popt, pcov))
    log_a0_est = torch.as_tensor(linear(n.log().cpu().numpy(), *popt))
    if shrink:
        y = torch.as_tensor(a0.log())
        y[torch.isnan(y)] = log_a0_est[torch.isnan(y)]
        log_a0_est = shrink_normal_normal(
            y, log_a0_est, shrink_prior_var=shrink_prior_var
        )
    return (torch.exp(log_a0_est), popt)


def get_pred_alpha0(X, sample_size_factors, popt, sample_mask=None):
    n_reps, n_condits, n_guides = X.shape
    if sample_mask is None:
        sample_mask = torch.ones((n_reps, n_condits), device="cpu")
    q = get_q(X + 1, sample_size_factors, sample_mask=sample_mask)
    a0_est = np.exp(linear(q.sum(axis=0).log().cpu().numpy(), *popt))
    return a0_est
