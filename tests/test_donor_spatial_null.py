"""Tests for the donor-level within-FOV spatial null.

These checks do not need the CosMx matrices. They lock three properties:
labels stay inside their FOV, donors are weighted equally, and a within-FOV
signal can produce a p-value below the 8-section sign-test floor.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from donor_spatial_null import (  # noqa: E402
    DONOR_SIGN_FLOOR_ONE,
    SAMPLE_SIGN_FLOOR_ONE,
    FovBlock,
    analyze,
    beats_donor_sign_floor,
    beats_sample_sign_floor,
    binned_means,
    perm_p,
    sign_floor,
    sign_test,
    within_fov_delta,
    within_fov_spearman,
)


NAMES = ("count",)


def _block(donor, sample, fov, cldn4, y) -> FovBlock:
    return FovBlock(
        donor=donor,
        sample=sample,
        fov=str(fov),
        cldn4=np.asarray(cldn4, dtype=float),
        outcomes=np.asarray(y, dtype=float).reshape(-1, 1),
        outcome_names=NAMES,
    )


def test_sign_floors_match_binomial():
    assert sign_floor(8, "greater") == SAMPLE_SIGN_FLOOR_ONE == 1 / 256
    assert sign_floor(5, "greater") == DONOR_SIGN_FLOOR_ONE == 1 / 32
    assert sign_floor(8, "two-sided") == 1 / 128
    assert sign_floor(5, "two-sided") == 1 / 16
    assert binomtest(8, 8, 0.5, alternative="greater").pvalue == SAMPLE_SIGN_FLOOR_ONE
    assert binomtest(5, 5, 0.5, alternative="greater").pvalue == DONOR_SIGN_FLOOR_ONE
    got = sign_test(np.array([-1.0, -0.2, -0.3, -4.0, -0.5, -0.7, -0.8, -0.1]))
    assert got["n_neg"] == 8 and got["p"] == SAMPLE_SIGN_FLOOR_ONE


def test_perm_p_includes_observed_draw():
    null = np.array([0.0, 0.1, -0.2])
    # Nothing in the null is <= -1, so p = 1/4.
    assert perm_p(-1.0, null, "less") == 0.25
    # Every null draw is <= 1, so p = 1.
    assert perm_p(1.0, null, "less") == 1.0


def test_within_fov_delta_and_spearman_direction():
    x = np.arange(10, dtype=float)
    y = -x
    assert within_fov_delta(x, y, min_arm=2) < 0
    assert within_fov_spearman(x, y) == -1.0
    # Flat mark is not a contrast.
    assert np.isnan(within_fov_delta(np.ones(10), y, min_arm=2))
    assert np.isnan(within_fov_spearman(np.ones(10), y))


def test_bins_are_low_to_high_cldn4():
    x = np.arange(20, dtype=float)
    y = -x
    means = binned_means(x, y, n_bins=5)
    assert np.all(np.diff(means) < 0)


def test_between_fov_geography_is_not_credited():
    """High CLDN4 sits in an immune-cold FOV. Shuffling inside FOVs cannot move it.

    The section-pooled contrast is large and negative, but it is invariant under
    the within-FOV null, so the permutation p is 1. The within-FOV contrast is 0.
    """
    cold = _block("D1", "S1", 1, cldn4=np.arange(10, 16), y=np.zeros(6))
    hot = _block("D1", "S1", 2, cldn4=np.arange(0, 6), y=np.full(6, 10.0))
    # Four extra donors so the donor mean is defined on more than one tissue.
    extras = []
    for i, donor in enumerate(("D2", "D3", "D4", "D5")):
        extras.append(_block(donor, f"S{i+2}", 1, cldn4=np.arange(10, 16), y=np.zeros(6)))
        extras.append(_block(donor, f"S{i+2}", 2, cldn4=np.arange(0, 6), y=np.full(6, 10.0)))
    out = analyze([cold, hot, *extras], n_perm=30, n_boot=20, seed=1, min_arm=2)
    pooled = out["section_pooled"]["delta"]
    within = out["within_fov"]["delta"]
    assert pooled["T"][0] == -10.0
    assert np.allclose(pooled["null_T"][:, 0], -10.0)
    assert pooled["p_one"][0] == 1.0
    assert within["T"][0] == 0.0
    assert within["p_one"][0] == 1.0
    # CLDN4 multisets never needed to leave the FOV for this conclusion.
    assert not beats_sample_sign_floor(pooled["p_one"][0])


def test_within_fov_signal_beats_sample_sign_floor():
    rng = np.random.default_rng(7)
    blocks = []
    for donor in ("D1", "D2", "D3", "D4", "D5"):
        for fov in (1, 2, 3):
            x = np.arange(24, dtype=float)
            y = -x + rng.normal(0, 0.01, size=x.size)
            blocks.append(_block(donor, donor, fov, x, y))
    out = analyze(blocks, n_perm=400, n_boot=30, seed=11, min_arm=2)
    spear = out["within_fov"]["spearman"]
    delta = out["within_fov"]["delta"]
    assert spear["T"][0] < -0.9
    assert delta["T"][0] < 0
    assert spear["p_one"][0] == 1 / 401
    assert beats_sample_sign_floor(spear["p_one"][0])
    assert beats_donor_sign_floor(spear["p_one"][0])
    assert spear["p_one"][0] < SAMPLE_SIGN_FLOOR_ONE
    # Five donors all negative: the sign test is stuck at its floor.
    assert spear["donor_sign"][0]["n_neg"] == 5
    assert spear["donor_sign"][0]["p"] == DONOR_SIGN_FLOOR_ONE
    # FOV-block interval for a near-deterministic effect excludes 0.
    assert spear["ci_block_hi"][0] < 0


def test_donors_are_weighted_equally_not_by_fov_count():
    blocks = []
    # One donor with four FOVs, each contrast -50.
    for fov in range(4):
        blocks.append(_block("D1", "D1", fov, [0, 1, 2, 3], [50, 50, 0, 0]))
    # Four donors with one FOV each, contrast +10.
    # cldn4 [0,1,2,3], y low=0, y high=10 -> delta = 10 - 0 = +10.
    for donor in ("D2", "D3", "D4", "D5"):
        blocks.append(_block(donor, donor, 1, [0, 1, 2, 3], [0, 0, 10, 10]))
    out = analyze(blocks, n_perm=10, n_boot=15, seed=3, min_arm=2)
    delta = out["within_fov"]["delta"]
    # Equal donor weight: (-50 + 10 + 10 + 10 + 10) / 5 = -2.
    assert delta["T"][0] == -2.0
    donor_map = dict(zip(out["donors"], delta["donor"][:, 0]))
    assert donor_map["D1"] == -50.0
    assert donor_map["D2"] == 10.0
    # A FOV-count weight would have been (-50*4 + 10*4) / 8 = -20.
    assert delta["T"][0] != -20.0


def test_null_p_is_not_automatically_small():
    rng = np.random.default_rng(99)
    ps = []
    for rep in range(12):
        blocks = []
        local = np.random.default_rng(1000 + rep)
        for donor in ("D1", "D2", "D3", "D4", "D5"):
            for fov in (1, 2):
                x = local.normal(size=30)
                y = local.normal(size=30)
                blocks.append(_block(donor, donor, fov, x, y))
        out = analyze(blocks, n_perm=40, n_boot=10, seed=50 + rep, min_arm=2)
        ps.append(out["within_fov"]["spearman"]["p_one"][0])
    ps = np.asarray(ps)
    assert np.median(ps) > 0.1
    assert np.mean(ps < 0.05) <= 0.25


def test_cell_weighted_sensitivity_differs_when_fovs_differ_in_size():
    # Large FOV with a strong negative contrast, tiny FOV with a positive one,
    # same donor and section. Cell weight and equal-FOV weight disagree.
    big = _block("D1", "S", 1, np.arange(20), -np.arange(20))
    tiny_x = np.array([0.0, 1.0, 2.0, 3.0])
    tiny_y = np.array([0.0, 0.0, 5.0, 5.0])  # delta +5
    tiny = _block("D1", "S", 2, tiny_x, tiny_y)
    out = analyze([big, tiny], n_perm=15, n_boot=10, seed=4, min_arm=2)
    equal = out["within_fov"]["delta"]["T"][0]
    weighted = out["within_fov"]["delta"]["T_cell_weighted"][0]
    assert equal > weighted  # the large negative FOV pulls the weighted mean down
    assert weighted < 0 < equal or weighted < equal
