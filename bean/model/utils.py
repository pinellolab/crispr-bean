import torch
import torch.distributions as tdist
import pyro
import pyro.distributions as dist
import pyro.distributions.constraints as constraints


def get_alpha(
    expected_guide_p, size_factor, sample_mask, a0, epsilon=1e-5, normalize_by_a0=True
):
    try:
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
    except:
        print(size_factor.shape)
        print(expected_guide_p.shape)
        print(a0.shape)
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
