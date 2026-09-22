"""Unit checks for partial correlation, Cliff's delta, and the selection rule."""

from __future__ import annotations

import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats import (  # noqa: E402
    bh_fdr,
    cliffs_batch,
    cliffs_delta,
    dl_random_effects,
    fisher_random_effects,
    matrix_partial_corr,
    partial_corr,
    select_abs,
    select_max,
    select_positive,
)


def test_partial_removes_shared_covariate():
    rng = np.random.default_rng(0)
    n = 400
    z = rng.normal(size=n)
    x = z + rng.normal(scale=0.1, size=n)
    y = z + rng.normal(scale=0.1, size=n)
    raw, raw_p, n_raw = partial_corr(x, y, [], rank=False)
    adj, adj_p, n_adj = partial_corr(x, y, [z], rank=False)
    assert n_raw == n and n_adj == n
    assert raw > 0.8
    assert abs(adj) < 0.15
    assert adj_p > 0.01 or abs(adj) < 0.2


def test_spearman_partial_matches_known_monotone():
    rng = np.random.default_rng(1)
    n = 300
    x = rng.normal(size=n)
    y = x + rng.normal(scale=0.5, size=n)
    rho, p, n_obs = partial_corr(x, y, [], rank=True)
    expect = stats.spearmanr(x, y).statistic
    assert n_obs == n
    assert abs(rho - expect) < 1e-10
    assert p < 1e-6


def test_cliff_delta_sign_and_mannwhitney_identity():
    rng = np.random.default_rng(2)
    low = rng.normal(loc=1.0, size=80)
    high = rng.normal(loc=0.0, size=90)
    delta, p, var, n1, n2 = cliffs_delta(high, low)
    assert n1 == 90 and n2 == 80
    assert delta < 0
    assert var > 0
    assert p < 0.01
    u = stats.mannwhitneyu(high, low, alternative="two-sided")
    mw_delta = 2 * u.statistic / (n1 * n2) - 1
    assert abs(delta - mw_delta) < 1e-10


def test_dl_recovers_common_effect():
    effects = np.array([-0.2, -0.2, -0.2])
    variances = np.array([0.01, 0.02, 0.04])
    pooled = dl_random_effects(effects, variances)
    assert pooled["n_cohorts"] == 3
    assert abs(pooled["effect"] + 0.2) < 1e-10
    assert pooled["I2"] == 0.0
    assert abs(pooled["fe"] + 0.2) < 1e-10


def test_fisher_pool_sign():
    rhos = np.array([-0.2, -0.25, -0.15, -0.22])
    ns = np.array([200, 180, 220, 160])
    pooled = fisher_random_effects(rhos, ns, k=4)
    assert pooled["n_cohorts"] == 4
    assert pooled["effect"] < 0
    assert pooled["ci_high"] < 0
    assert 0 <= pooled["I2"] <= 1


def test_bh_orders_and_caps():
    q = bh_fdr(np.array([0.01, 0.04, np.nan, 0.5]))
    assert np.isnan(q[2])
    assert q[0] <= q[1] <= q[3]
    assert q[3] <= 1


def test_cliffs_batch_matches_scalar():
    rng = np.random.default_rng(4)
    high = rng.normal(size=(40, 3))
    low = rng.normal(loc=0.5, size=(35, 3))
    high[0, 2] = np.nan
    delta, p, var = cliffs_batch(high, low)
    assert np.isnan(delta[2])
    for j in range(2):
        d, pj, v, _, _ = cliffs_delta(high[:, j], low[:, j])
        assert abs(d - delta[j]) < 1e-10
        assert abs(v - var[j]) < 1e-10
        assert abs(pj - p[j]) < 1e-10


def test_matrix_partial_matches_scalar():
    rng = np.random.default_rng(3)
    n = 120
    cov = rng.normal(size=(n, 2))
    pred = rng.normal(size=(n, 2)) + cov @ np.array([[0.4, 0.1], [0.2, 0.5]])
    out = rng.normal(size=(n, 2)) + cov @ np.array([[0.3, -0.2], [0.1, 0.4]])
    rho, p = matrix_partial_corr(pred, out, [cov[:, 0], cov[:, 1]], df=n - 2 - 2)
    for i in range(2):
        for j in range(2):
            r, pj, n_obs = partial_corr(pred[:, i], out[:, j], [cov[:, 0], cov[:, 1]], rank=False)
            assert n_obs == n
            assert abs(r - rho[i, j]) < 1e-8
            assert abs(pj - p[i, j]) < 1e-8


def test_select_max_respects_seven_of_eight_and_sign():
    effects = np.array([
        [-0.30, -0.30, -0.30],  # large but only 6/8 for the first predictor
        [-0.12, -0.14, -0.16],  # eligible, smaller
        [-0.20, -0.05, 0.02],   # not all negative
        [-0.19, -0.18, -0.15],  # eligible, larger minimum
    ])
    neg = np.array([
        [6, 8, 8],
        [7, 8, 8],
        [8, 7, 4],
        [7, 8, 8],
    ])
    worst_p = np.array([1e-8, 0.01, 0.2, 0.002])
    mean_abs = np.mean(np.abs(effects), axis=1)
    idx = select_max(effects, neg, worst_p, mean_abs)
    assert idx == 3
    none = select_max(effects[:1], neg[:1], worst_p[:1], mean_abs[:1])
    assert none == -1


def test_select_abs_requires_six_negative_and_picks_largest():
    effects = np.array([-0.40, -0.25, 0.50, -0.30, -0.25])
    neg = np.array([5, 6, 8, 7, 8])
    cohorts = np.array([8, 6, 8, 8, 8])
    pvals = np.array([1e-6, 0.01, 1e-9, 0.02, 0.001])
    # -0.40 has only 5 negative cohorts. +0.50 is the wrong direction.
    # -0.30 (7/8) beats -0.25 (8/8) on absolute effect.
    assert select_abs(effects, neg, cohorts, pvals) == 3
    # Tie on |effect|: more negative cohorts wins, then smaller p.
    effects2 = np.array([-0.20, -0.20, -0.20])
    neg2 = np.array([6, 8, 7])
    cohorts2 = np.array([8, 8, 8])
    pvals2 = np.array([0.001, 0.2, 0.0001])
    assert select_abs(effects2, neg2, cohorts2, pvals2) == 1
    assert select_positive(effects, cohorts, pvals) == 2
    assert select_positive(np.array([-0.2, -0.3]), np.array([8, 8]), np.array([0.1, 0.1])) == -1


if __name__ == "__main__":
    test_partial_removes_shared_covariate()
    test_spearman_partial_matches_known_monotone()
    test_cliff_delta_sign_and_mannwhitney_identity()
    test_dl_recovers_common_effect()
    test_fisher_pool_sign()
    test_bh_orders_and_caps()
    test_cliffs_batch_matches_scalar()
    test_matrix_partial_matches_scalar()
    test_select_max_respects_seven_of_eight_and_sign()
    test_select_abs_requires_six_negative_and_picks_largest()
    print("ok")
