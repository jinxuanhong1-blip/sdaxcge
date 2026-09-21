#!/usr/bin/env python3
"""Algebra checks for partial Spearman and rank-OLS mediation. No TCGA download."""

from __future__ import annotations

import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from bootstrap import OUTCOMES, RHO_MODELS, bootstrap_primary
from io_tcga import patient_id
from stats import (
    attenuation_fraction,
    bh_fdr,
    fisher_ci,
    partial_spearman,
    product_mediation,
    random_effects_meta,
)


def test_spearman_without_covariates_matches_scipy():
    rng = np.random.default_rng(1)
    x = rng.normal(size=80)
    y = 0.4 * x + rng.normal(size=80)
    rho, p, n = partial_spearman(x, y, ())
    ref = stats.spearmanr(x, y)
    assert n == 80
    assert abs(rho - ref.statistic) < 1e-10
    assert abs(p - ref.pvalue) < 1e-8


def test_partial_spearman_matches_residual_pearson():
    rng = np.random.default_rng(2)
    n = 90
    z1 = rng.normal(size=n)
    z2 = rng.normal(size=n)
    x = 0.7 * z1 + rng.normal(size=n)
    y = -0.5 * z1 + 0.4 * z2 + 0.25 * x + rng.normal(size=n)
    rho, p, nobs = partial_spearman(x, y, [z1, z2])
    assert nobs == n and np.isfinite(p)
    xr, yr = stats.rankdata(x), stats.rankdata(y)
    design = np.column_stack([np.ones(n), stats.rankdata(z1), stats.rankdata(z2)])
    ex = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    ey = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    ref, _ = stats.pearsonr(ex, ey)
    assert abs(rho - ref) < 1e-10


def test_mediation_identity_and_linear_proportion():
    rng = np.random.default_rng(3)
    n = 4000
    z = rng.normal(size=(n, 2))
    x = 0.3 * z[:, 0] + rng.normal(size=n)
    a_true, b_true, direct = 1.15, -0.72, -0.18
    mediator = a_true * x + z @ np.array([0.2, -0.1]) + rng.normal(size=n) * 0.35
    y = direct * x + b_true * mediator + z @ np.array([0.05, 0.08]) + rng.normal(size=n) * 0.35
    raw = product_mediation(x, mediator, y, z, rank=False)
    assert abs(raw["identity_gap"]) < 1e-8
    truth = (a_true * b_true) / (direct + a_true * b_true)
    assert abs(raw["proportion"] - truth) < 0.04
    ranked = product_mediation(x, mediator, y, z, rank=True)
    assert abs(ranked["identity_gap"]) < 1e-8
    assert ranked["proportion"] > 0.5


def test_noise_mediator_is_near_zero_and_true_mediator_is_not():
    rng = np.random.default_rng(4)
    n = 2500
    x = rng.normal(size=n)
    mediator = 1.3 * x + rng.normal(size=n) * 0.4
    y = -0.9 * mediator + rng.normal(size=n) * 0.45
    noise = rng.normal(size=n)
    pm_true = product_mediation(x, mediator, y, rank=True)["proportion"]
    pm_noise = product_mediation(x, noise, y, rank=True)["proportion"]
    assert pm_true > 0.8
    assert abs(pm_noise) < 0.12


def test_attenuation_fraction_sign():
    # -0.10 moving to -0.04 is 60% of the association accounted for.
    assert abs(attenuation_fraction(-0.10, -0.04) - 0.60) < 1e-12
    # Moving farther from zero is a negative share.
    assert attenuation_fraction(-0.10, -0.16) < 0
    assert np.isnan(attenuation_fraction(0.0, -0.02))


def test_random_effects_identical_cohorts():
    out = random_effects_meta([0.2, 0.2, 0.2, 0.2], [220, 220, 220, 220], k=0)
    assert out["n_cohorts"] == 4
    assert abs(out["rho"] - 0.2) < 1e-10
    assert out["I2"] < 1e-8
    lo, hi = fisher_ci(0.2, 220, 0)
    assert lo < 0.2 < hi
    negative = random_effects_meta([-0.2] * 6, [300] * 6, k=3)
    assert negative["rho"] < 0
    assert negative["p"] < 1e-6


def test_bh_fdr_orders_and_caps():
    q = bh_fdr([0.01, 0.04, np.nan, 0.03])
    assert np.isnan(q[2])
    assert q[0] <= q[2] or np.isnan(q[2])
    finite = q[np.isfinite(q)]
    assert np.all(finite <= 1.0 + 1e-12)
    assert np.all(finite >= 0.0)


def test_patient_id_keeps_primary_tumor_only():
    assert patient_id("TCGA-05-4244-01A-01R-1107-07") == "TCGA-05-4244"
    assert patient_id("TCGA.05.4244.01A.01R") == "TCGA-05-4244"
    assert patient_id("TCGA-05-4244-11A-01R") is None
    assert patient_id("TCGA-05-4244-06A-01R") is None


def test_bootstrap_is_deterministic_and_separates_mediator_from_control():
    rng = np.random.default_rng(5)
    n = 350
    keratin = rng.normal(size=n)
    x = 0.35 * keratin + rng.normal(size=n)
    cldn4 = 1.0 * x + 0.2 * keratin + rng.normal(size=n) * 0.45
    cldn7 = rng.normal(size=n)
    epcam = 0.15 * keratin + rng.normal(size=n)
    cd8 = -0.85 * cldn4 + 0.1 * keratin + rng.normal(size=n) * 0.5
    cd3 = -0.4 * cldn4 + rng.normal(size=n)
    data = {
        "TACSTD2": x,
        "CLDN4": cldn4,
        "CLDN7": cldn7,
        "EPCAM": epcam,
        "KRT8": keratin,
        "KRT18": keratin + rng.normal(size=n) * 0.1,
        "KRT19": keratin + rng.normal(size=n) * 0.1,
        "CD8_score": cd8,
        "CD3_score": cd3,
    }
    first = bootstrap_primary(data, n_boot=40, seed=11)
    second = bootstrap_primary(data, n_boot=40, seed=11)
    assert np.allclose(first["rho"]["CD8_score"]["keratin"], second["rho"]["CD8_score"]["keratin"])
    assert first["full_pm"]["CD8_score"]["CLDN4"] > 0.6
    assert abs(first["full_pm"]["CD8_score"]["CLDN7"]) < 0.2
    # The fast full-sample rho matches the public partial Spearman.
    rho, _, nobs = partial_spearman(
        data["TACSTD2"], data["CD8_score"], [data[name] for name in RHO_MODELS["keratin"]]
    )
    assert nobs == n
    assert abs(rho - first["full_rho"]["CD8_score"]["keratin"]) < 1e-8
    assert set(first["full_rho"]) == set(OUTCOMES)


if __name__ == "__main__":
    test_spearman_without_covariates_matches_scipy()
    test_partial_spearman_matches_residual_pearson()
    test_mediation_identity_and_linear_proportion()
    test_noise_mediator_is_near_zero_and_true_mediator_is_not()
    test_attenuation_fraction_sign()
    test_random_effects_identical_cohorts()
    test_bh_fdr_orders_and_caps()
    test_patient_id_keeps_primary_tumor_only()
    test_bootstrap_is_deterministic_and_separates_mediator_from_control()
    print("all stats tests passed")
