from typing import Optional
import torch
import pyro
from pyro import poutine
import pyro.distributions as dist
import torch.distributions.constraints as constraints
from .utils import get_alpha, scale_pi_by_accessibility, MAX_LOGPI
from ..preprocessing.data_class import (
    VariantSurvivalScreenData,
    VariantSurvivalReporterScreenData,
    TilingSurvivalReporterScreenData,
)


def NormalModel(
    data: VariantSurvivalScreenData,
    mask_thres: int = 10,
    use_bcmatch: bool = True,
    prior_params: Optional[dict] = None,
):
    """
    Fit only on guide counts

    Args:
        data: input data
        mask_thres: threshold for masking guide counts for stability. Defaults to 10.
        use_bcmatch: whether to use barcode-matched counts. Defaults to True.
        sd_scale: scale for prior standard deviation. Defaults to 0.01 for improved identifiability.
        prior_params: prior parameters. If provided, specified prior parameters will be used.
    """
    replicate_plate = pyro.plate("rep_plate", data.n_reps, dim=-3)
    replicate_plate2 = pyro.plate("rep_plate2", data.n_reps, dim=-2)
    time_plate = pyro.plate("time_plate", data.n_condits, dim=-2)
    guide_plate = pyro.plate("guide_plate", data.n_guides, dim=-1)

    mu_dist = dist.Laplace(0, 1)
    initial_abundance = torch.ones(data.n_guides) / data.n_guides
    if prior_params is not None:
        if "mu_loc" in prior_params or "mu_scale" in prior_params:
            mu_loc = 0.0
            mu_scale = 1.0
            if "mu_loc" in prior_params:
                mu_loc = prior_params["mu_loc"]
            if "mu_scale" in prior_params:
                mu_scale = prior_params["mu_scale"]
            mu_dist = dist.Normal(mu_loc, mu_scale)
        if "initial_abundance" in prior_params:
            initial_abundance = prior_params["initial_abundance"]

    # Set the prior for phenotype means
    with pyro.plate("guide_plate0", 1):
        with pyro.plate("target_plate", data.n_targets):
            # In survival analysis, fitted effect size is not
            mu_targets = pyro.sample("mu_targets", mu_dist)

    mu_center = mu_targets
    mu = torch.repeat_interleave(mu_center, data.target_lengths, dim=0)
    if hasattr(data, "negctrl_guide_idx"):
        mu[data.negctrl_guide_idx, :] = 0.0
    r = torch.exp(mu)
    assert r.shape == (data.n_guides, 1)

    with pyro.plate("replicate_plate0", data.n_reps, dim=-1):
        q_0 = pyro.sample(
            "initial_guide_abundance",
            dist.Dirichlet(initial_abundance.unsqueeze(0).expand(data.n_reps, -1)),
        )
    with replicate_plate:
        with time_plate as t:
            time = data.timepoints[t]
            assert time.shape == (data.n_condits,)
            with guide_plate:
                alleles_p_time = torch.pow(
                    r.unsqueeze(0).expand((data.n_condits, -1, -1)),
                    time.unsqueeze(-1).unsqueeze(-1).expand((-1, data.n_guides, 1)),
                )
                assert alleles_p_time.shape == (data.n_condits, data.n_guides, 1)

            expected_allele_p = alleles_p_time.unsqueeze(0).expand(
                data.n_reps, -1, -1, -1
            ) * q_0.unsqueeze(1).unsqueeze(-1).expand((-1, data.n_condits, -1, -1))
            expected_guide_p = expected_allele_p.sum(axis=-1)
            assert expected_guide_p.shape == (
                data.n_reps,
                data.n_condits,
                data.n_guides,
            ), expected_guide_p.shape

    with replicate_plate2:
        with pyro.plate("guide_plate3", data.n_guides, dim=-1):
            a = get_alpha(expected_guide_p, data.size_factor, data.sample_mask, data.a0)
            assert data.X.shape == (
                data.n_reps,
                data.n_condits,
                data.n_guides,
            )
            with poutine.mask(
                mask=torch.logical_and(
                    data.X_masked.permute(0, 2, 1).sum(axis=-1) > mask_thres,
                    data.repguide_mask,
                )
            ):
                pyro.sample(
                    "guide_counts",
                    dist.DirichletMultinomial(a, validate_args=False),
                    obs=data.X_masked.permute(0, 2, 1),
                )
            if use_bcmatch:
                a_bcmatch = get_alpha(
                    expected_guide_p,
                    data.size_factor_bcmatch,
                    data.sample_mask,
                    data.a0_bcmatch,
                )
                with poutine.mask(
                    mask=torch.logical_and(
                        data.X_bcmatch_masked.permute(0, 2, 1).sum(axis=-1)
                        > mask_thres,
                        data.repguide_mask,
                    )
                ):
                    assert (
                        a_bcmatch > 0
                    ).all(), f"{torch.where(a_bcmatch < 1e-6)}\n{torch.where(torch.isnan(a_bcmatch))}"
                    pyro.sample(
                        "guide_bcmatch_counts",
                        dist.DirichletMultinomial(a_bcmatch, validate_args=False),
                        obs=data.X_bcmatch_masked.permute(0, 2, 1),
                    )
    return alleles_p_time


def ControlNormalModel(data, mask_thres=10, use_bcmatch=True):
    """
    Fit shared control distribution
    """
    replicate_plate = pyro.plate("rep_plate", data.n_reps, dim=-3)
    replicate_plate2 = pyro.plate("rep_plate2", data.n_reps, dim=-2)
    time_plate = pyro.plate("time_plate", data.n_condits, dim=-2)
    guide_plate = pyro.plate("guide_plate", data.n_guides, dim=-1)

    # Set the prior for phenotype means
    mu_targets = pyro.sample("mu_targets", dist.Laplace(0, 1))
    mu = mu_targets.repeat(data.n_guides).unsqueeze(-1)
    r = torch.exp(mu)
    with pyro.plate("rep_plate1", data.n_reps, dim=-1):
        q_0 = pyro.sample(
            "initial_guide_abundance",
            dist.Dirichlet(torch.ones((data.n_reps, data.n_guides))),
        )
    with replicate_plate:
        with time_plate as t:
            time = data.timepoints[t]
            assert time.shape == (data.n_condits,)
            with guide_plate:
                alleles_p_time = r.unsqueeze(0).expand(
                    (data.n_condits, -1, -1)
                ) ** time.unsqueeze(-1).unsqueeze(-1).expand((-1, data.n_guides, 1))
                assert alleles_p_time.shape == (data.n_condits, data.n_guides, 1)

            expected_allele_p = alleles_p_time.unsqueeze(0).expand(
                data.n_reps, -1, -1, -1
            ) * q_0.unsqueeze(1).unsqueeze(-1).expand((-1, data.n_condits, -1, -1))
            expected_guide_p = expected_allele_p.sum(axis=-1)
            assert expected_guide_p.shape == (
                data.n_reps,
                data.n_condits,
                data.n_guides,
            ), expected_guide_p.shape

    with replicate_plate2:
        with pyro.plate("guide_plate3", data.n_guides, dim=-1):
            a = get_alpha(expected_guide_p, data.size_factor, data.sample_mask, data.a0)
            assert data.X.shape == (
                data.n_reps,
                data.n_condits,
                data.n_guides,
            ), (
                data.X.shape,
                (
                    data.n_reps,
                    data.n_condits,
                    data.n_guides,
                ),
            )
            with poutine.mask(
                mask=torch.logical_and(
                    data.X_masked.permute(0, 2, 1).sum(axis=-1) > mask_thres,
                    data.repguide_mask,
                )
            ):
                pyro.sample(
                    "guide_counts",
                    dist.DirichletMultinomial(a, validate_args=False),
                    obs=data.X_masked.permute(0, 2, 1),
                )
            if use_bcmatch:
                a_bcmatch = get_alpha(
                    expected_guide_p,
                    data.size_factor_bcmatch,
                    data.sample_mask,
                    data.a0_bcmatch,
                )
                with poutine.mask(
                    mask=torch.logical_and(
                        data.X_bcmatch_masked.permute(0, 2, 1).sum(axis=-1)
                        > mask_thres,
                        data.repguide_mask,
                    )
                ):
                    pyro.sample(
                        "guide_bcmatch_counts",
                        dist.DirichletMultinomial(a_bcmatch, validate_args=False),
                        obs=data.X_bcmatch_masked.permute(0, 2, 1),
                    )

    return alleles_p_time


def MixtureNormalModel(
    data: VariantSurvivalReporterScreenData,
    alpha_prior: float = 1,
    use_bcmatch: bool = True,
    use_all_timepoints_for_pi: bool = True,
    sd_scale: float = 0.01,
    scale_by_accessibility: bool = False,
    fit_noise: bool = False,
    mask_thres: int = 10,
    prior_params: Optional[dict] = None,
):
    """
    Using the reporter outcome, phenotype of cells with a guide will be modeled as mixture of two normal distributions of edited and unedited cells.

    Args:
        data: Input data of type VariantSortingReporterScreenData.
        alpha_prior: Prior parameter for controlling the concentration of the Dirichlet process. Defaults to 1.
        use_bcmatch: Flag indicating whether to use barcode-matched counts. Defaults to True.
        use_all_timepoints_for_pi: Use all available timepoints instead of the `--control-condition` timepoint.
        sd_scale: Scale for the prior standard deviation. Defaults to 0.01.
        scale_by_accessibility: If True, pi fitted from reporter data is scaled by accessibility.
        fit_noise: Valid only when scale_by_accessibility is True. If True, parametrically fit noise of endo ~ reporter + noise.
        prior_params: Optional dictionary of prior parameters. If provided, specified prior parameters will be used.
    """
    torch.autograd.set_detect_anomaly(True)
    replicate_plate = pyro.plate("rep_plate", data.n_reps, dim=-3)
    replicate_plate2 = pyro.plate("rep_plate2", data.n_reps, dim=-2)
    time_plate = pyro.plate("time_plate", data.n_condits, dim=-2)
    guide_plate = pyro.plate("guide_plate", data.n_guides, dim=-1)

    mu_dist = dist.Laplace(0, 1)
    initial_abundance = torch.ones(data.n_guides) / data.n_guides
    if prior_params is not None:
        if "mu_loc" in prior_params or "mu_scale" in prior_params:
            mu_loc = 0.0
            mu_scale = 1.0
            if "mu_loc" in prior_params:
                mu_loc = prior_params["mu_loc"]
            if "mu_scale" in prior_params:
                mu_scale = prior_params["mu_scale"]
            mu_dist = dist.Normal(mu_loc, mu_scale)
        if "initial_abundance" in prior_params:
            initial_abundance = prior_params["initial_abundance"]

    # Set the prior for phenotype means
    with pyro.plate("guide_plate0", 1):
        with pyro.plate("guide_plate1", data.n_targets):
            mu_targets = pyro.sample("mu_targets", mu_dist)
    # with pyro.plate("negctrl_plate", len(data.negctrl_guide_idx)):
    #     mu_negctrl = pyro.sample("mu_negctrl", dist.Normal(0, 1))
    mu_center = torch.cat([torch.zeros((data.n_targets, 1)), mu_targets], axis=-1)
    mu = torch.repeat_interleave(mu_center, data.target_lengths, dim=0)
    # Fix negative control's mu to be 0
    if hasattr(data, "negctrl_guide_idx"):
        mu[data.negctrl_guide_idx, :] = 0.0  # mu_negctrl[:, None]
    assert mu.shape == (data.n_guides, 2)
    r = torch.exp(mu)

    # with pyro.plate("replicate_plate0", data.n_reps, dim=-1):
    #     q_0 = pyro.sample(
    #         "initial_guide_abundance",
    #         dist.Dirichlet(initial_abundance.unsqueeze(0).expand(data.n_reps, -1)),
    #     )
    alpha_pi = pyro.param(
        "alpha_pi",
        torch.ones(
            (
                data.n_guides,
                2,
            )
        )
        * alpha_prior,
        constraint=constraints.positive,
    )
    pi_a_scaled = alpha_pi / alpha_pi.sum(axis=-1)[:, None] * data.pi_a0[:, None]
    assert alpha_pi.shape == (
        data.n_guides,
        2,
    ), alpha_pi.shape

    with replicate_plate:
        with guide_plate, poutine.mask(mask=data.repguide_mask.unsqueeze(1)):
            pi = pyro.sample(
                "pi",
                dist.Dirichlet(
                    pi_a_scaled.unsqueeze(0).unsqueeze(0).expand(data.n_reps, 1, -1, -1)
                ),
            )
            assert pi.shape == (
                data.n_reps,
                1,
                data.n_guides,
                2,
            ), pi.shape
        with pyro.plate("time_plate0", len(data.control_timepoint), dim=-2) as t:
            with guide_plate, poutine.mask(mask=data.repguide_mask.unsqueeze(1)):
                # if use_all_timepoints_for_pi:
                #     time_pi = data.timepoints
                #     expanded_allele_p = pi * r.expand(
                #         data.n_reps, len(data.timepoints), -1, -1
                #     ) ** time_pi.unsqueeze(0).unsqueeze(-1).unsqueeze(-1).expand(
                #         data.n_reps, len(data.timepoints), -1, -1
                #     )
                #     pyro.sample(
                #         "allele_count",
                #         dist.Multinomial(probs=expanded_allele_p, validate_args=False),
                #         obs=data.allele_counts,
                #     )
                # else:
                time_pi = data.control_timepoint[t]
                # If pi is sampled in later timepoint, account for the selection.
                expanded_allele_p = pi * torch.pow(
                    r.expand(data.n_reps, 1, -1, -1), time_pi
                )
                pyro.sample(
                    "control_allele_count",
                    dist.Multinomial(probs=expanded_allele_p, validate_args=False),
                    obs=data.allele_counts_control,
                )
    if scale_by_accessibility:
        # Endogenous target site editing rate may be different
        pi = scale_pi_by_accessibility(
            pi, data.guide_accessibility, fit_noise=fit_noise
        )
    with replicate_plate:
        with time_plate as t:
            time = data.timepoints[t]
            assert time.shape == (data.n_condits,)

            with guide_plate:
                alleles_p_time = torch.pow(
                    r.unsqueeze(0).expand((data.n_condits, -1, -1)),
                    time.unsqueeze(-1).unsqueeze(-1).expand((-1, data.n_guides, 1)),
                )
                # alleles_p_time = torch.clamp(
                #     time.unsqueeze(-1).unsqueeze(-1).expand((-1, data.n_guides, 1))
                #     * torch.log(r).unsqueeze(0).expand((data.n_condits, -1, -1)),
                #     max=MAX_LOGPI,
                # ).exp()
                # negctrl_abundance = pyro.param(
                #     "negctrl_abundance",
                #     torch.ones((data.n_condits,)),
                #     constraint=constraints.positive,
                # )
                # alleles_p_time = (
                #     alleles_p_time / negctrl_abundance.clamp(min=1e-5)[:, None, None]
                # )
                assert alleles_p_time.shape == (data.n_condits, data.n_guides, 2)

            expected_allele_p = (
                pi.expand(-1, data.n_condits, -1, -1) * alleles_p_time[None, :, :, :]
            )  # * q_0.unsqueeze(1).unsqueeze(-1).expand((-1, data.n_condits, -1, -1))
            expected_guide_p = expected_allele_p.sum(axis=-1)
            assert expected_guide_p.shape == (
                data.n_reps,
                data.n_condits,
                data.n_guides,
            ), expected_guide_p.shape

    with replicate_plate2:
        with pyro.plate("guide_plate3", data.n_guides, dim=-1):
            a = get_alpha(expected_guide_p, data.size_factor, data.sample_mask, data.a0)

            assert (
                data.X.shape
                == data.X_bcmatch_masked.shape
                == (
                    data.n_reps,
                    data.n_condits,
                    data.n_guides,
                )
            )
            with poutine.mask(
                mask=torch.logical_and(
                    data.X_masked.permute(0, 2, 1).sum(axis=-1) > mask_thres,
                    data.repguide_mask,
                )
            ):
                pyro.sample(
                    "guide_counts",
                    dist.DirichletMultinomial(a, validate_args=False),
                    obs=data.X_masked.permute(0, 2, 1),
                )
            if use_bcmatch:
                a_bcmatch = get_alpha(
                    expected_guide_p,
                    data.size_factor_bcmatch,
                    data.sample_mask,
                    data.a0_bcmatch,
                )
                with poutine.mask(
                    mask=torch.logical_and(
                        data.X_bcmatch_masked.permute(0, 2, 1).sum(axis=-1)
                        > mask_thres,
                        data.repguide_mask,
                    )
                ):
                    pyro.sample(
                        "guide_bcmatch_counts",
                        dist.DirichletMultinomial(a_bcmatch, validate_args=False),
                        obs=data.X_bcmatch_masked.permute(0, 2, 1),
                    )


def MultiMixtureNormalModel(
    data: TilingSurvivalReporterScreenData,
    alpha_prior=1,
    use_bcmatch=True,
    use_all_timepoints_for_pi: bool = True,
    sd_scale=0.01,
    norm_pi=False,
    scale_by_accessibility=False,
    fit_noise: bool = False,
    prior_params: Optional[dict] = None,
    epsilon=1e-5,
):
    """
    Using the reporter outcome, phenotype of cells with a guide will be modeled as mixture of normal distributions of all major alleles (including WT) produced by the guide.

    Args:
        data: Input data of type VariantSortingReporterScreenData.
        alpha_prior: Prior parameter for controlling the concentration of the Dirichlet process. Defaults to 1.
        use_bcmatch: Flag indicating whether to use barcode-matched counts. Defaults to True.
        use_all_timepoints_for_pi: Use all available timepoints instead of the `--control-condition` timepoint.
        sd_scale: Scale for the prior standard deviation. Defaults to 0.01.
        scale_by_accessibility: If True, pi fitted from reporter data is scaled by accessibility.
        fit_noise: Valid only when scale_by_accessibility is True. If True, parametrically fit noise of endo ~ reporter + noise.
        prior_params: Optional dictionary of prior parameters. If provided, specified prior parameters will be used.
        epsilon: Small value to avoid division by zero, assigned as Dirichlet parameters for non-existing alleles.
    """

    replicate_plate = pyro.plate("rep_plate", data.n_reps, dim=-3)
    replicate_plate2 = pyro.plate("rep_plate2", data.n_reps, dim=-2)
    time_plate = pyro.plate("time_plate", data.n_condits, dim=-2)
    guide_plate = pyro.plate("guide_plate", data.n_guides, dim=-1)

    mu_dist = dist.Laplace(0, 1)
    initial_abundance = torch.ones(data.n_guides) / data.n_guides
    if prior_params is not None:
        if "mu_loc" in prior_params or "mu_scale" in prior_params:
            mu_loc = 0.0
            mu_scale = 1.0
            if "mu_loc" in prior_params:
                mu_loc = prior_params["mu_loc"]
            if "mu_scale" in prior_params:
                mu_scale = prior_params["mu_scale"]
            mu_dist = dist.Normal(mu_loc, mu_scale)
        if "initial_abundance" in prior_params:
            initial_abundance = prior_params["initial_abundance"]

    # Set the prior for phenotype means
    with pyro.plate("guide_plate1", data.n_edits):
        mu_edits = pyro.sample("mu_targets", mu_dist)
    assert mu_edits.shape == (data.n_edits,)
    assert data.allele_to_edit.shape == (
        data.n_guides,
        data.n_max_alleles - 1,
        data.n_edits,
    )
    mu_targets = torch.matmul(data.allele_to_edit, mu_edits)
    assert mu_targets.shape == (data.n_guides, data.n_max_alleles - 1)

    mu = torch.cat([torch.zeros((data.n_guides, 1)), mu_targets], axis=-1)
    r = torch.exp(mu)

    with pyro.plate("replicate_plate0", data.n_reps, dim=-1):
        q_0 = pyro.sample(
            "initial_guide_abundance",
            dist.Dirichlet(initial_abundance.unsqueeze(0).expand(data.n_reps, -1)),
        )
    # The pi should be Dirichlet distributed instead of independent betas
    alpha_pi0 = (
        torch.ones(
            (
                data.n_guides,
                data.n_max_alleles,
            )
        )
        * alpha_prior
    )
    # Effectively remove alphas for non-existing alleles
    assert data.allele_mask.shape == (data.n_guides, data.n_max_alleles)
    alpha_pi0[~data.allele_mask] = epsilon
    alpha_pi = pyro.param("alpha_pi", alpha_pi0, constraint=constraints.positive)
    alpha_pi[~data.allele_mask] = epsilon
    pi_a_scaled = (
        (alpha_pi + epsilon / alpha_pi.shape[-1])
        / (alpha_pi.sum(axis=-1)[:, None] + epsilon)
        * data.pi_a0[:, None]
    )
    pi_a_scaled[pi_a_scaled < epsilon] = epsilon
    if torch.isnan(pi_a_scaled).any():
        print(torch.where(alpha_pi.isnan()))
        print(torch.where(alpha_pi < 0))
        exit(1)
    if (pi_a_scaled <= 0).any():
        print(torch.where(alpha_pi < 0))

    with replicate_plate:
        with guide_plate, poutine.mask(mask=data.repguide_mask.unsqueeze(1)):
            pi = pyro.sample(
                "pi",
                dist.Dirichlet(
                    pi_a_scaled.unsqueeze(0).unsqueeze(0).expand(data.n_reps, 1, -1, -1)
                ),
            )
        with pyro.plate("time_plate0", len(data.control_timepoint), dim=-2):
            with guide_plate, poutine.mask(mask=data.repguide_mask.unsqueeze(1)):
                # if use_all_timepoints_for_pi:
                #     time_pi = data.timepoints
                #     # If pi is sampled in later timepoint, account for the selection.
                #     expanded_allele_p = pi * r.expand(
                #         data.n_reps, 1, -1, -1
                #     ) ** time_pi.unsqueeze(0).unsqueeze(-1).unsqueeze(-1).expand(
                #         data.n_reps, len(data.timepoints), -1, -1
                #     )
                #     pyro.sample(
                #         "allele_count",
                #         dist.Multinomial(probs=expanded_allele_p, validate_args=False),
                #         obs=data.allele_counts,
                #     )
                # else:
                time_pi = data.control_timepoint
                # If pi is sampled in later timepoint, account for the selection.
                expanded_allele_p = pi * r.expand(data.n_reps, 1, -1, -1) ** time_pi
                pyro.sample(
                    "control_allele_count",
                    dist.Multinomial(probs=expanded_allele_p, validate_args=False),
                    obs=data.allele_counts_control,
                )
    if scale_by_accessibility:
        # Endogenous target site editing rate may be different
        pi = scale_pi_by_accessibility(
            pi, data.guide_accessibility, fit_noise=fit_noise
        )

    with replicate_plate:
        with time_plate as t:
            time = data.timepoints[t]
            assert time.shape == (data.n_condits,)

            with guide_plate, poutine.mask(mask=data.repguide_mask.unsqueeze(1)):
                alleles_p_time = torch.clamp(
                    time.unsqueeze(-1).unsqueeze(-1).expand((-1, data.n_guides, 1))
                    * torch.log(r).unsqueeze(0).expand((data.n_condits, -1, -1)),
                    max=MAX_LOGPI,
                ).exp()

                mask = data.allele_mask.unsqueeze(0).expand((data.n_condits, -1, -1))
                alleles_p_time = alleles_p_time * mask

                assert alleles_p_time.shape == (
                    data.n_condits,
                    data.n_guides,
                    data.n_max_alleles,
                )
            expected_allele_p = (
                pi.expand(data.n_reps, data.n_condits, -1, -1)
                * alleles_p_time[None, :, :, :]
                * q_0.unsqueeze(1).unsqueeze(-1).expand((-1, data.n_condits, -1, -1))
            )
            expected_guide_p = expected_allele_p.sum(axis=-1)
            assert expected_guide_p.shape == (
                data.n_reps,
                data.n_condits,
                data.n_guides,
            ), expected_guide_p.shape
    try:
        with replicate_plate2:
            with pyro.plate("guide_plate3", data.n_guides, dim=-1):
                a = get_alpha(
                    expected_guide_p, data.size_factor, data.sample_mask, data.a0
                )
                a_bcmatch = get_alpha(
                    expected_guide_p,
                    data.size_factor_bcmatch,
                    data.sample_mask,
                    data.a0_bcmatch,
                )
                # assert a.shape == a_bcmatch.shape == (data.n_reps, data.n_guides, data.n_condits)
                assert (
                    data.X.shape
                    == data.X_bcmatch_masked.shape
                    == (
                        data.n_reps,
                        data.n_condits,
                        data.n_guides,
                    )
                )
                with poutine.mask(
                    mask=torch.logical_and(
                        data.X_masked.permute(0, 2, 1).sum(axis=-1) > 10,
                        data.repguide_mask,
                    )
                ):
                    pyro.sample(
                        "guide_counts",
                        dist.DirichletMultinomial(a, validate_args=False),
                        obs=data.X_masked.permute(0, 2, 1),
                    )
                if use_bcmatch:
                    with poutine.mask(
                        mask=torch.logical_and(
                            data.X_bcmatch_masked.permute(0, 2, 1).sum(axis=-1) > 10,
                            data.repguide_mask,
                        )
                    ):
                        pyro.sample(
                            "guide_bcmatch_counts",
                            dist.DirichletMultinomial(a_bcmatch, validate_args=False),
                            obs=data.X_bcmatch_masked.permute(0, 2, 1),
                        )
    except ValueError as e:
        print(f"ERROR a is 0 at {torch.sum(a.sum(axis=-1) ==0)}")
        print(
            f"ERROR expected_guide_p is 0 at {torch.sum(expected_guide_p.sum(axis=1) ==0)}"
        )
        print(f"ERROR a is NaN at {torch.where(a.isnan().any(axis=-1))}")
        print(
            f"ERROR data.size_factor is NaN at {torch.where(data.size_factor.isnan())}"
        )
        print(
            f"ERROR expected_guide_p is NaN at {torch.where(expected_guide_p.isnan().any(axis=1))}"
        )
        print(f"ERROR a0 is NaN at {torch.where(data.a0.isnan())}")
        raise e


def NormalGuide(data):
    initial_abundance = pyro.param(
        "initial_abundance",
        torch.ones(data.n_guides) / data.n_guides,
        constraint=constraints.positive,
    )
    with pyro.plate("replicate_plate0", data.n_reps, dim=-1):
        q_0 = pyro.sample(
            "initial_guide_abundance",
            dist.Dirichlet(initial_abundance),
        )
    with pyro.plate("guide_plate0", 1):
        with pyro.plate("guide_plate1", data.n_targets):
            mu_loc = pyro.param("mu_loc", torch.zeros((data.n_targets, 1)))
            mu_scale = pyro.param(
                "mu_scale",
                torch.ones((data.n_targets, 1)),
                constraint=constraints.positive,
            )
            pyro.sample("mu_targets", dist.Normal(mu_loc, mu_scale))


def MixtureNormalGuide(
    data,
    alpha_prior: float = 1,
    use_bcmatch: bool = True,
    scale_by_accessibility: bool = False,
    fit_noise: bool = False,
):
    replicate_plate = pyro.plate("rep_plate", data.n_reps, dim=-3)
    guide_plate = pyro.plate("guide_plate", data.n_guides, dim=-1)
    initial_abundance = pyro.param(
        "initial_abundance",
        torch.ones(data.n_guides) / data.n_guides,
        constraint=constraints.positive,
    )
    with pyro.plate("replicate_plate0", data.n_reps, dim=-1):
        q_0 = pyro.sample(
            "initial_guide_abundance",
            dist.Dirichlet(initial_abundance),
        )
    # Set the prior for phenotype means
    mu_loc = pyro.param("mu_loc", torch.zeros((data.n_targets, 1)))
    mu_scale = pyro.param(
        "mu_scale", torch.ones((data.n_targets, 1)), constraint=constraints.positive
    )
    with pyro.plate("guide_plate0", 1):
        with pyro.plate("guide_plate1", data.n_targets):
            mu_targets = pyro.sample("mu_targets", dist.Normal(mu_loc, mu_scale))
    mu_center = torch.cat([torch.zeros((data.n_targets, 1)), mu_targets], axis=-1)
    mu = torch.repeat_interleave(mu_center, data.target_lengths, dim=0)
    assert mu.shape == (data.n_guides, 2)

    # The pi should be Dirichlet distributed instead of independent betas
    alpha_pi = pyro.param(
        "alpha_pi",
        torch.ones(
            (
                data.n_guides,
                2,
            )
        )
        * alpha_prior,
        constraint=constraints.positive,
    )
    assert alpha_pi.shape == (
        data.n_guides,
        2,
    ), alpha_pi.shape
    pi_a_scaled = alpha_pi / alpha_pi.sum(axis=-1)[:, None] * data.pi_a0[:, None]

    with replicate_plate:
        with guide_plate:
            pi = pyro.sample(
                "pi",
                dist.Dirichlet(
                    pi_a_scaled.unsqueeze(0)
                    .unsqueeze(0)
                    .expand(data.n_reps, 1, -1, -1)
                    .clamp(1e-5)
                ),
            )
            assert pi.shape == (
                data.n_reps,
                1,
                data.n_guides,
                2,
            ), pi.shape
    if scale_by_accessibility:
        # Endogenous target site editing rate may be different
        pi = scale_pi_by_accessibility(
            pi, data.guide_accessibility, fit_noise=fit_noise
        )


def ControlNormalGuide(data, mask_thres=10, use_bcmatch=True):
    """
    Fit shared mean
    """
    # Set the prior for phenotype means
    mu_loc = pyro.param("mu_loc", torch.tensor(0.0))
    mu_scale = pyro.param(
        "mu_scale", torch.tensor(1.0), constraint=constraints.positive
    )
    pyro.sample("mu_targets", dist.Normal(mu_loc, mu_scale))


def MultiMixtureNormalGuide(
    data,
    alpha_prior=1,
    use_bcmatch=True,
    epsilon=1e-5,
    scale_by_accessibility: bool = False,
    fit_noise: bool = False,
):
    """
    Guide for model C14
    """
    replicate_plate = pyro.plate("rep_plate", data.n_reps, dim=-3)
    guide_plate = pyro.plate("guide_plate", data.n_guides, dim=-1)

    initial_abundance = pyro.param(
        "initial_abundance",
        torch.ones(data.n_guides) / data.n_guides,
        constraint=constraints.positive,
    )
    with pyro.plate("replicate_plate0", data.n_reps, dim=-1):
        q_0 = pyro.sample(
            "initial_guide_abundance",
            dist.Dirichlet(initial_abundance),
        )
    # Set the prior for phenotype means
    mu_loc = pyro.param("mu_loc", torch.zeros((data.n_edits,)))
    mu_scale = pyro.param(
        "mu_scale", torch.ones((data.n_edits,)), constraint=constraints.positive
    )
    with pyro.plate("guide_plate1", data.n_edits):
        mu_edits = pyro.sample("mu_targets", dist.Normal(mu_loc, mu_scale))
    mu_targets = torch.matmul(data.allele_to_edit, mu_edits)
    assert mu_targets.shape == (data.n_guides, data.n_max_alleles - 1), (
        mu_targets.shape,
        data.n_max_alleles,
        data.n_edits,
    )

    mu = torch.cat([torch.zeros((data.n_guides, 1)), mu_targets], axis=-1)
    assert mu.shape == (data.n_guides, data.n_max_alleles), (mu.shape,)
    # The pi should be Dirichlet distributed instead of independent betas
    alpha_pi0 = (
        torch.ones(
            (
                data.n_guides,
                data.n_max_alleles,
            )
        )
        * alpha_prior
    )
    # Effectively remove alphas for non-existing alleles
    assert data.allele_mask.shape == (data.n_guides, data.n_max_alleles)
    alpha_pi0[~data.allele_mask] = epsilon
    alpha_pi = pyro.param("alpha_pi", alpha_pi0, constraint=constraints.positive)
    alpha_pi[~data.allele_mask] = epsilon
    pi_a_scaled = alpha_pi / alpha_pi.sum(axis=-1)[:, None] * data.pi_a0[:, None]

    with replicate_plate:
        with guide_plate, poutine.mask(mask=data.repguide_mask.unsqueeze(1)):
            pi = pyro.sample(
                "pi",
                dist.Dirichlet(
                    pi_a_scaled.unsqueeze(0)
                    .unsqueeze(0)
                    .expand(data.n_reps, 1, -1, -1)
                    .clamp(1e-5)
                ),
            )

            assert pi.shape == (
                data.n_reps,
                1,
                data.n_guides,
                data.n_max_alleles,
            ), pi.shape
    if scale_by_accessibility:
        # Endogenous target site editing rate may be different
        pi = scale_pi_by_accessibility(
            pi, data.guide_accessibility, fit_noise=fit_noise
        )
