"""Unit checks for the Kras-lung Cldn4 sweep. No download."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "methods" / "kras_lung_icb_cldn4_sweep"))

import analyze
from analyze import (
    bh,
    extreme_split,
    gse155972_samples,
    gse246922_kp,
    gse246922_llc,
    safe_spearman,
)


def test_bh_orders_and_caps():
    q = bh([0.01, 0.04, 0.03, float("nan")])
    assert abs(q[0] - 0.03) < 1e-12
    assert q[1] == 0.04
    assert q[2] == 0.04
    assert np.isnan(q[3])
    assert max(v for v in q if v == v) <= 1


def test_spearman_constant_is_nan():
    rho, p = safe_spearman(np.ones(6), np.arange(6))
    assert np.isnan(rho) and np.isnan(p)


def test_perfect_spearman_p_is_not_zero():
    x = np.arange(6, dtype=float)
    rho, p = safe_spearman(x, x)
    assert rho == 1
    # 2 / 6! for a two-sided perfect rank correlation without ties.
    assert abs(p - (2 / 720)) < 1e-12


def test_quartile_tails_do_not_overlap():
    values = np.arange(12, dtype=float)
    high, low = extreme_split(values, "quartile")
    assert high.sum() == 3 and low.sum() == 3
    assert not np.any(high & low)
    assert values[high].min() > values[low].max()


def test_kp_relapse_mapping_is_not_kl():
    samples = gse246922_kp(
        ["KP_1", "KPy_1", "ResKP2_1", "ResResKP_1", "ResKPlate_1"]
    )
    by = {s.name: s for s in samples}
    assert by["KP_1"].role == "reference" and by["KP_1"].ordinal == 0
    assert by["KPy_1"].role == "other"
    assert by["ResKP2_1"].group == "relapse_1" and by["ResKP2_1"].role == "resistant"
    assert by["ResResKP_1"].group == "relapse_2"
    assert by["ResKPlate_1"].group == "relapse_late" and by["ResKPlate_1"].ordinal == 3
    llc = gse246922_llc(["LLC1_1", "LLC1y_1", "ResResLLC1_1"])
    assert llc[0].role == "reference" and llc[2].role == "resistant" and llc[1].role == "other"


def test_tismo_llc_response_is_setdb1_genotype():
    cols = [
        "llc_control_sg5_UnTx_rep1_quant",
        "llc_setdb1_sg4_UnTx_rep1_quant",
        "llc_control_sg5_Tx_rep1_quant",
        "llc_setdb1_sg4_Tx_rep1_quant",
    ]
    samples = gse155972_samples(cols)
    by = {s.name: s for s in samples}
    assert by["llc_control_sg5_Tx_rep1_quant"].role == "resistant"
    assert by["llc_setdb1_sg4_Tx_rep1_quant"].role == "reference"
    assert by["llc_control_sg5_UnTx_rep1_quant"].role == "other"
    assert by["llc_setdb1_sg4_UnTx_rep1_quant"].role == "other"
    # Untreated must not be read as treated via the _tx_ substring.
    assert all(s.role == "other" for s in samples if "UnTx" in s.name)


def test_inventory_does_not_call_llc_kl():
    for row in analyze.inventory_rows():
        if row["model"] in {"LLC", "LLC1", "LL/2"}:
            assert "not" in row["note"].lower() or "KrasG12C" in row["kras"]
            assert row["model"] != "KL"
    assert any(r["model"] == "CMT-167" and r["scored"] == "no" for r in analyze.inventory_rows())
    assert any("not rerun" in r["note"].lower() or "Not rerun" in r["note"] for r in analyze.inventory_rows())


def test_synthetic_resistant_high_cldn4_is_aligned():
    resistant = np.array([5.0, 6.0, 7.0])
    reference = np.array([1.0, 2.0, 1.5])
    effect = analyze.median_diff(resistant, reference)
    assert effect > 0
    # Exact Mann-Whitney on 3 vs 3 cannot go below 0.1 even with full separation.
    assert analyze.safe_mwu(resistant, reference) == 0.1
