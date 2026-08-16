#!/usr/bin/env python3
"""Unit + smoke tests for hier_meta.py. No real cohort findings."""

from __future__ import annotations

import math
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import hier_meta as hm  # noqa: E402


def _assert(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def test_empty_results_dir_is_quiet() -> None:
    tmp = tempfile.mkdtemp()
    try:
        rc = hm.main(["--results-dir", os.path.join(tmp, "nope"), "--outdir", os.path.join(tmp, "out")])
        _assert(rc == 0, f"missing dir should exit 0, got {rc}")
        empty = os.path.join(tmp, "empty")
        os.makedirs(empty)
        rc = hm.main(["--results-dir", empty, "--outdir", os.path.join(tmp, "out2")])
        _assert(rc == 0, f"empty dir should exit 0, got {rc}")
    finally:
        shutil.rmtree(tmp)


def test_open_only_drops_ega() -> None:
    rows = [
        {"cohort_id": "GSE000001", "gene": "TACSTD2", "endpoint": "DCB", "unit": "hedges_g",
         "yi": "0.2", "sei": "0.4", "access": "open", "__file__": "a.csv", "__line__": 2},
        {"cohort_id": "GSE000002", "gene": "TACSTD2", "endpoint": "DCB", "unit": "hedges_g",
         "yi": "-0.1", "sei": "0.3", "access": "geo", "__file__": "a.csv", "__line__": 3},
        {"cohort_id": "EGAS00001005013", "gene": "TACSTD2", "endpoint": "DCB", "unit": "hedges_g",
         "yi": "0.9", "sei": "0.15", "access": "ega", "__file__": "a.csv", "__line__": 4},
        {"cohort_id": "phs002822", "gene": "TACSTD2", "endpoint": "DCB", "unit": "hedges_g",
         "yi": "0.7", "sei": "0.12", "access": "dbgap", "__file__": "a.csv", "__line__": 5},
    ]
    cfg = hm.deep_update(hm.DEFAULT_CONFIG, {"gene": "TACSTD2", "primary_endpoint": "DCB",
                                             "effect_family": "smd", "access_policy": "open_only"})
    kept, audit = hm.harmonise(rows, cfg)
    ids = {c.cohort_id for c in kept}
    _assert(ids == {"GSE000001", "GSE000002"}, f"open_only kept {ids}")
    dropped = [a for a in audit if a["kept"] == "no"]
    _assert(any("open_only" in a["reason"] for a in dropped), "restricted rows must cite open_only")


def test_gse_without_access_column_is_open() -> None:
    rows = [{"cohort_id": "GSE126044", "gene": "CLDN4", "endpoint": "DCB", "unit": "hedges_g",
             "yi": "0.1", "sei": "0.5", "__file__": "a.csv", "__line__": 2}]
    cfg = hm.deep_update(hm.DEFAULT_CONFIG, {"gene": "CLDN4", "primary_endpoint": "DCB",
                                             "effect_family": "smd"})
    kept, _ = hm.harmonise(rows, cfg)
    _assert(len(kept) == 1 and kept[0].access == "open", "GSE* defaults to open")


def test_conjugate_known_mean() -> None:
    """Equal-precision observations of a common value: posterior mean of mu ~ that value."""
    y = np.array([0.30, 0.30, 0.30, 0.30, 0.30])
    V = np.diag(np.full(5, 0.04))
    prior = {"mu_mean": 0.0, "mu_sd": 100.0, "tau_family": "half_normal", "tau_scale": 0.05}
    fit = hm.fit_bayes(y, V, prior, {"max": 1.0, "n": 401})
    m, sd = fit.mu.moments()
    _assert(abs(m - 0.30) < 0.02, f"posterior mean {m} drifted from 0.30")
    _assert(sd < 0.12, f"posterior sd {sd} unexpectedly wide")


def test_noisier_cohorts_shrink_more() -> None:
    y = np.array([1.0, 1.0, 1.0])
    V = np.diag([0.01, 0.25, 1.0])
    prior = {"mu_mean": 0.0, "mu_sd": 1.0, "tau_family": "half_normal", "tau_scale": 0.4}
    fit = hm.fit_bayes(y, V, prior, {"max": 2.0, "n": 401})
    _assert(fit.shrinkage[0] < fit.shrinkage[1] < fit.shrinkage[2],
            f"shrinkage not monotone in SE: {fit.shrinkage}")


def test_dl_zero_when_homogeneous() -> None:
    y = np.array([0.1, 0.1, 0.1, 0.1])
    v = np.array([0.04, 0.04, 0.04, 0.04])
    _assert(hm.tau2_dl(y, v) == 0.0, "DL tau2 should be 0 for identical effects")


def test_demo_end_to_end() -> None:
    tmp = tempfile.mkdtemp()
    try:
        rc = hm.main(["--demo", "--genes", "TACSTD2,CLDN4", "--outdir", tmp])
        _assert(rc == 0, f"demo rc={rc}")
        for gene in ("TACSTD2", "CLDN4"):
            pooled = os.path.join(tmp, gene, "pooled_estimates.csv")
            _assert(os.path.isfile(pooled), f"missing {pooled}")
            shrink = os.path.join(tmp, gene, "cohort_shrinkage.csv")
            text = open(shrink, encoding="utf-8").read()
            _assert("SYNTH_EGA_01" not in text, f"{gene} pooled a restricted row")
            _assert("SYNTH_OPEN_01" in text, f"{gene} dropped an open row")
            audit = open(os.path.join(tmp, gene, "input_audit.csv"), encoding="utf-8").read()
            _assert("open_only" in audit, "audit must record restricted drops")
            miss = os.path.join(tmp, gene, "missingness_summary.csv")
            _assert(os.path.isfile(miss), "missingness summary missing")
    finally:
        shutil.rmtree(tmp)


def test_smd_does_not_exp() -> None:
    y = np.array([0.2, -0.1, 0.0])
    v = np.array([0.16, 0.16, 0.16])
    r = hm.random_effects(y, v, 0.0)
    _assert(math.isfinite(r["estimate"]), "RE failed")


if __name__ == "__main__":
    tests = [
        test_empty_results_dir_is_quiet,
        test_open_only_drops_ega,
        test_gse_without_access_column_is_open,
        test_conjugate_known_mean,
        test_noisier_cohorts_shrink_more,
        test_dl_zero_when_homogeneous,
        test_smd_does_not_exp,
        test_demo_end_to_end,
    ]
    failed = 0
    for fn in tests:
        try:
            fn()
            print(f"  ok  {fn.__name__}")
        except Exception as exc:
            failed += 1
            print(f" FAIL {fn.__name__}: {exc}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    raise SystemExit(failed)
