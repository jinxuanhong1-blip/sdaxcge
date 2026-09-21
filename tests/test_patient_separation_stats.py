"""Unit checks for patient-level separation. No CosMx matrix required."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from patient_separation_stats import (  # noqa: E402
    bootstrap_median_ci,
    cohort_spearman_dl,
    corrected_log2_ratio,
    exclusion_magnitude,
    labels_within_groups,
    pairwise_median_diff,
    perm_p_ge,
    permute_within_groups,
    pick_winner,
    pooled_median_diff,
    rollup_section_means,
)


ROOT = Path(__file__).resolve().parents[1]
PATIENTS = ROOT / "inputs" / "concordant4_patient_units.tsv"


def test_five_patient_bootstrap_interval_is_the_sample_range():
    rng = np.random.default_rng(0)
    values = np.array([-1.6, -1.1, -0.9, -0.7, -0.35])
    med, lo, hi = bootstrap_median_ci(values, rng, n_boot=50)
    assert med == np.median(values)
    assert lo == -1.6
    assert hi == -0.35


def test_exclusion_magnitude_requires_the_interval_below_zero():
    assert exclusion_magnitude(-0.8, -0.1) == 0.8
    assert exclusion_magnitude(-0.8, 0.05) == 0.0
    assert exclusion_magnitude(0.4, 0.9) == 0.0
    assert exclusion_magnitude(float("nan"), -0.2) == 0.0


def test_winner_is_the_larger_median_not_the_tighter_small_effect():
    rows = [
        {"median": -0.05, "lo": -0.06, "hi": -0.04, "eligible": True},
        {"median": -1.20, "lo": -1.80, "hi": -0.40, "eligible": True},
        {"median": -2.00, "lo": -2.50, "hi": 0.10, "eligible": True},
    ]
    winner = pick_winner(rows)
    assert winner["median"] == -1.20


def test_winner_tie_prefers_the_narrower_interval():
    rows = [
        {"median": -1.0, "lo": -2.0, "hi": -0.2, "eligible": True},
        {"median": -1.0, "lo": -1.3, "hi": -0.5, "eligible": True},
    ]
    winner = pick_winner(rows)
    assert winner["hi"] == -0.5


def test_detected_and_quartile_labels():
    values = np.array([0, 0, 0, 2, 5, 9])
    groups = np.zeros(6, dtype=int)
    ties = np.arange(6)
    det = labels_within_groups(values, groups, "detected", ties)
    assert list(det) == [0, 0, 0, 1, 1, 1]
    q = labels_within_groups(values, groups, "quartile", ties)
    # k = 1. Lowest value is low, highest value is high.
    assert q[0] == 0
    assert q[-1] == 1
    assert int(np.sum(q == 1)) == 1
    assert int(np.sum(q == 0)) == 1


def test_section_means_are_averaged_inside_the_patient():
    # Two sections of patient 0 must not count as two patients.
    section_high = np.array([[2.0], [4.0], [10.0]])
    section_low = np.array([[4.0], [8.0], [10.0]])
    sample_patient = np.array([0, 0, 1])
    high, low = rollup_section_means(section_high, section_low, sample_patient, n_patients=2)
    assert high[0, 0] == 3.0
    assert low[0, 0] == 6.0
    assert high[1, 0] == 10.0
    ratio = corrected_log2_ratio(high[:, 0], low[:, 0], eps=0.0)
    med, lo, hi = bootstrap_median_ci(ratio, np.random.default_rng(1), n_boot=20)
    assert med == np.median(ratio)
    assert hi == max(ratio)
    assert lo == min(ratio)


def test_within_group_permutation_preserves_group_counts():
    rng = np.random.default_rng(2)
    lab = np.array([1, 1, 0, 0, 1, 0])
    groups = np.array([0, 0, 0, 1, 1, 1])
    shuffled = permute_within_groups(lab, groups, rng)
    assert int(shuffled[:3].sum()) == int(lab[:3].sum())
    assert int(shuffled[3:].sum()) == int(lab[3:].sum())


def test_pairwise_median_does_not_cross_cohorts():
    # Cohort 0: high T/NK is lower. Cohort 1: same direction, different level.
    y = np.array([0.10, 0.20, 0.50, 0.40])
    lab = np.array([1, 0, 1, 0])
    cohort = np.array([0, 0, 1, 0])
    # Cohort 1 has only one high and no low: the contrast is undefined.
    assert np.isnan(pairwise_median_diff(y, lab, cohort))
    y2 = np.array([0.10, 0.30, 0.40, 0.50, 0.20, 0.25, 0.60, 0.70])
    lab2 = np.array([1, 1, 0, 0, 1, 1, 0, 0])
    cohort2 = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    pair0 = (np.array([0.10, 0.30])[:, None] - np.array([0.40, 0.50])[None, :]).ravel()
    pair1 = (np.array([0.20, 0.25])[:, None] - np.array([0.60, 0.70])[None, :]).ravel()
    # Every within-cohort pair is negative. A cross-cohort contrast would mix
    # the low arm of cohort 0 (0.40) with the high arm of cohort 1 (0.20).
    delta = pairwise_median_diff(y2, lab2, cohort2)
    assert delta < 0
    assert delta == np.median(np.concatenate([pair0, pair1]))
    assert pooled_median_diff(y, lab) == np.median([0.10, 0.50]) - np.median([0.20, 0.40])


def test_permutation_p_includes_the_observed_draw():
    assert perm_p_ge(np.array([0.1, 0.2, 0.3]), 0.3) == 2 / 4
    assert perm_p_ge(np.array([0.1, 0.2]), 5.0) == 1 / 3


def test_locked_concordant4_spearman_is_reproduced():
    import csv

    rows = list(csv.DictReader(PATIENTS.open(), delimiter="\t"))
    assert len(rows) == 65
    cohort_name = np.array([r["dataset"] for r in rows])
    codes, names = np.unique(cohort_name, return_inverse=True)
    score = np.array([float(r["mal_CLDN4_pct"]) for r in rows])
    outcome = np.array([float(r["frac_tnk"]) for r in rows])
    dl = cohort_spearman_dl(score, outcome, names)
    by_name = {str(codes[i]): dl["cohort_rho"][i] for i in range(len(codes))}
    assert round(by_name["GSE123902"], 3) == -0.659
    assert round(by_name["GSE131907"], 3) == -0.522
    assert round(by_name["GSE205335"], 3) == -0.435
    assert round(by_name["GSE189357"], 3) == -0.600
    assert round(dl["rho"], 3) == -0.531
    assert round(dl["ci_lo"], 3) == -0.697
    assert round(dl["ci_hi"], 3) == -0.312
    assert abs(dl["p"] - 1.65e-5) < 1e-6
    assert dl["i2"] == 0.0


def pytest_approx(value, rel=1e-4):
    class _A:
        def __eq__(self, other):
            if value == 0:
                return abs(other) <= rel
            return abs(other - value) <= rel * abs(value)

        def __repr__(self):
            return f"approx({value})"

    return _A()
