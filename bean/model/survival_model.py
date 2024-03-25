import torch
import pyro
from pyro import poutine
import pyro.distributions as dist
import torch.distributions.constraints as constraints
from .utils import (
    get_alpha,
    scale_pi_by_accessibility,
)
from ..preprocessing.data_class import (
    VariantSurvivalScreenData,
    VariantSurvivalReporterScreenData,
    TilingSurvivalReporterScreenData,
)


def NormalModel(
    data: VariantSurvivalScreenData, mask_thres: int = 10, use_bcmatch: bool = True
):
    """
    Fit only on guide counts
    Args
    --
    data: input data
    """
    replicate_plate = pyro.plate("rep_plate", data.n_reps, dim=-3)
    replicate_plate2 = pyro.plate("rep_plate2", data.n_reps, dim=-2)
    time_plate = pyro.plate("time_plate", data.n_condits, dim=-2)
    guide_plate = pyro.plate("guide_plate", data.n_guides, dim=-1)

    # Set the prior for phenotype means
    with pyro.plate("guide_plate0", 1):
        with pyro.plate("target_plate", data.n_targets):
            # In survival analysis, fitted effect size is not
            mu_alleles = pyro.sample("mu_alleles", dist.Laplace(0, 1))

    mu_center = mu_alleles
    mu = torch.repeat_interleave(mu_center, data.target_lengths, dim=0)
    r = torch.exp(mu)
    assert r.shape == (data.n_guides, 1)

    with pyro.plate("replicate_plate0", data.n_reps, dim=-1):
        q_0 = pyro.sample(
            "initial_guide_abundance",
            dist.Dirichlet(torch.ones((data.n_reps, data.n_guides))),
        )
    with replicate_plate:
        with time_plate as t:
            time = data.timepoints[t]
            assert time.shape == (data.n_condits,)
            # with guide_plate, poutine.mask(mask=(data.allele_counts.sum(axis=-1) == 0)):
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
    mu_alleles = pyro.sample("mu_alleles", dist.Laplace(0, 1))
    mu = mu_alleles.repeat(data.n_guides).unsqueeze(-1)
    r = torch.exp(mu)
    with replicate_plate:
        with pyro.plate("guide_plate2", data.n_guides):
            q_0 = pyro.sample(
                "initial_guide_abundance",
                dist.Dirichlet(torch.ones((data.n_reps, data.n_guides))),
            )
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

            assert (
                data.X.shape
                == data.X_bcmatch.shape
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

    return alleles_p_time


def MixtureNormalModel(
    data: VariantSurvivalReporterScreenData,
    alpha_prior: float = 1,
    use_bcmatch: bool = True,
    sd_scale: float = 0.01,
    scale_by_accessibility: bool = False,
    fit_noise: bool = False,
):
    """
    Args:
        scale_by_accessibility: If True, pi fitted from reporter data is scaled by accessibility
    """
    torch.autograd.set_detect_anomaly(True)
    replicate_plate = pyro.plate("rep_plate", data.n_reps, dim=-3)
    replicate_plate2 = pyro.plate("rep_plate2", data.n_reps, dim=-2)
    time_plate = pyro.plate("time_plate", data.n_condits, dim=-2)
    guide_plate = pyro.plate("guide_plate", data.n_guides, dim=-1)

    # Set the prior for phenotype means
    with pyro.plate("guide_plate0", 1):
        with pyro.plate("guide_plate1", data.n_targets):
            mu_alleles = pyro.sample("mu_alleles", dist.Laplace(0, 1))
    mu_center = torch.cat([torch.zeros((data.n_targets, 1)), mu_alleles], axis=-1)
    mu = torch.repeat_interleave(mu_center, data.target_lengths, dim=0)
    assert mu.shape == (data.n_guides, 2)
    r = torch.exp(mu)

    with pyro.plate("replicate_plate0", data.n_reps, dim=-1):
        q_0 = pyro.sample(
            "initial_guide_abundance",
            dist.Dirichlet(torch.ones((data.n_reps, data.n_guides))),
        )
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
    pi_a_scaled = alpha_pi / alpha_pi.sum(axis=-1)[:, None] * data.pi_a0[:, None]
    assert alpha_pi.shape == (
        data.n_guides,
        2,
    ), alpha_pi.shape
    with replicate_plate:
        with guide_plate, poutine.mask(mask=data.repguide_mask.unsqueeze(1)):
            # Accounting for sample specific overall edit rate across all guides.
            # P(allele | guide, bin=bulk)
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
            pyro.sample(
                "bulk_allele_count",
                dist.Multinomial(probs=pi, validate_args=False),
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

            # with guide_plate, poutine.mask(mask=(data.allele_counts.sum(axis=-1) == 0)):
            with guide_plate:
                alleles_p_time = r.unsqueeze(0).expand(
                    (data.n_condits, -1, -1)
                ) ** time.unsqueeze(-1).unsqueeze(-1).expand((-1, data.n_guides, 1))
                assert alleles_p_time.shape == (data.n_condits, data.n_guides, 2)

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

    with replicate_plate2:
        with pyro.plate("guide_plate3", data.n_guides, dim=-1):
            a = get_alpha(expected_guide_p, data.size_factor, data.sample_mask, data.a0)
            a_bcmatch = get_alpha(
                expected_guide_p,
                data.size_factor_bcmatch,
                data.sample_mask,
                data.a0_bcmatch,
            )
            # a_bcmatch = get_alpha(expected_guide_p, data.size_factor_bcmatch, data.sample_mask, data.a0_bcmatch)
            # assert a.shape == a_bcmatch.shape == (data.n_reps, data.n_guides, data.n_condits)
            assert (
                data.X.shape
                == data.X_bcmatch.shape
                == (
                    data.n_reps,
                    data.n_condits,
                    data.n_guides,
                )
            )
            with poutine.mask(
                mask=torch.logical_and(
                    data.X_masked.permute(0, 2, 1).sum(axis=-1) > 10, data.repguide_mask
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


def NormalGuide(data):
    with pyro.plate("guide_plate0", 1):
        with pyro.plate("guide_plate1", data.n_targets):
            mu_loc = pyro.param("mu_loc", torch.zeros((data.n_targets, 1)))
            mu_scale = pyro.param(
                "mu_scale",
                torch.ones((data.n_targets, 1)),
                constraint=constraints.positive,
            )
            pyro.sample("mu_alleles", dist.Normal(mu_loc, mu_scale))


def MixtureNormalGuide(
    data,
    alpha_prior: float = 1,
    use_bcmatch: bool = True,
    scale_by_accessibility: bool = False,
    fit_noise: bool = False,
):
    """
    Guide for MixtureNormal model
    """

    replicate_plate = pyro.plate("rep_plate", data.n_reps, dim=-3)
    guide_plate = pyro.plate("guide_plate", data.n_guides, dim=-1)

    # Set the prior for phenotype means
    mu_loc = pyro.param("mu_loc", torch.zeros((data.n_targets, 1)))
    mu_scale = pyro.param(
        "mu_scale", torch.ones((data.n_targets, 1)), constraint=constraints.positive
    )
    with pyro.plate("guide_plate0", 1):
        with pyro.plate("guide_plate1", data.n_targets):
            mu_alleles = pyro.sample("mu_alleles", dist.Normal(mu_loc, mu_scale))
    mu_center = torch.cat([torch.zeros((data.n_targets, 1)), mu_alleles], axis=-1)
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
    pyro.sample("mu_alleles", dist.Normal(mu_loc, mu_scale))


def MultiMixtureNormalModel(
    data: TilingSurvivalReporterScreenData,
    alpha_prior=1,
    use_bcmatch=True,
    sd_scale=0.01,
    norm_pi=False,
    scale_by_accessibility=False,
    epsilon=1e-5,
    fit_noise: bool = False,
):
    """Tiling version of MixtureNormalModel"""

    replicate_plate = pyro.plate("rep_plate", data.n_reps, dim=-3)
    replicate_plate2 = pyro.plate("rep_plate2", data.n_reps, dim=-2)
    time_plate = pyro.plate("time_plate", data.n_condits, dim=-2)
    guide_plate = pyro.plate("guide_plate", data.n_guides, dim=-1)

    # Set the prior for phenotype means
    with pyro.plate("guide_plate1", data.n_edits):
        mu_edits = pyro.sample("mu_alleles", dist.Laplace(0, 1))
        sd_edits = pyro.sample(
            "sd_alleles",
            dist.LogNormal(
                torch.zeros((data.n_edits,)),
                torch.ones(
                    data.n_edits,
                )
                * sd_scale,
            ),
        )
    assert mu_edits.shape == sd_edits.shape == (data.n_edits,)
    assert data.allele_to_edit.shape == (
        data.n_guides,
        data.n_max_alleles - 1,
        data.n_edits,
    )
    mu_alleles = torch.matmul(data.allele_to_edit, mu_edits)
    assert mu_alleles.shape == (data.n_guides, data.n_max_alleles - 1)

    mu = torch.cat([torch.zeros((data.n_guides, 1)), mu_alleles], axis=-1)
    r = torch.exp(mu)
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
            pyro.sample(
                "bulk_allele_count",
                dist.Multinomial(probs=pi, validate_args=False),
                obs=data.allele_counts_control,
            )
    if scale_by_accessibility:
        # Endogenous target site editing rate may be different
        pi = scale_pi_by_accessibility(
            pi, data.guide_accessibility, fit_noise=fit_noise
        )

    with replicate_plate:
        with pyro.plate("guide_plate2", data.n_guides):
            q_0 = pyro.sample(
                "initial_guide_abundance",
                dist.Dirichlet(torch.ones((data.n_reps, data.n_guides))),
            )
        with time_plate as t:
            time = data.timepoints[t]
            assert time.shape == (data.n_condits,)
            with guide_plate, poutine.mask(mask=data.repguide_mask.unsqueeze(1)):
                alleles_p_time = r.unsqueeze(0).expand(
                    (data.n_condits, -1, -1)
                ) ** time.unsqueeze(-1).unsqueeze(-1).expand((-1, data.n_guides, 1))
                mask = data.allele_mask.unsqueeze(0).expand((data.n_condits, -1, -1))
                alleles_p_time[~mask] = 0.0

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

    with replicate_plate2:
        with pyro.plate("guide_plate3", data.n_guides, dim=-1):
            a = get_alpha(expected_guide_p, data.size_factor, data.sample_mask, data.a0)
            a_bcmatch = get_alpha(
                expected_guide_p,
                data.size_factor_bcmatch,
                data.sample_mask,
                data.a0_bcmatch,
            )
            # a_bcmatch = get_alpha(expected_guide_p, data.size_factor_bcmatch, data.sample_mask, data.a0_bcmatch)
            # assert a.shape == a_bcmatch.shape == (data.n_reps, data.n_guides, data.n_condits)
            assert (
                data.X.shape
                == data.X_bcmatch.shape
                == (
                    data.n_reps,
                    data.n_condits,
                    data.n_guides,
                )
            )
            with poutine.mask(
                mask=torch.logical_and(
                    data.X_masked.permute(0, 2, 1).sum(axis=-1) > 10, data.repguide_mask
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

    # Set the prior for phenotype means
    mu_loc = pyro.param("mu_loc", torch.zeros((data.n_edits,)))
    mu_scale = pyro.param(
        "mu_scale", torch.ones((data.n_edits,)), constraint=constraints.positive
    )
    with pyro.plate("guide_plate1", data.n_edits):
        mu_edits = pyro.sample("mu_alleles", dist.Normal(mu_loc, mu_scale))
    mu_alleles = torch.matmul(data.allele_to_edit, mu_edits)
    assert mu_alleles.shape == (data.n_guides, data.n_max_alleles - 1), (
        mu_alleles.shape,
        data.n_max_alleles,
        data.n_edits,
    )

    mu = torch.cat([torch.zeros((data.n_guides, 1)), mu_alleles], axis=-1)
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
                    # .clamp(1e-5)
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
