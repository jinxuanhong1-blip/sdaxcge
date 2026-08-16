"""Synthetic-data checks for the B1_BRCA ranking helpers."""
from __future__ import annotations

import numpy as np
import pandas as pd

from run_analysis import (
    _bh_fdr,
    _corr_rows_vs_vector,
    one_per_patient,
    partial_spearman,
    rank_partners,
    sample_type,
)


def test_one_per_patient_keeps_first_vial():
    cols = [
        "TCGA-A1-AAAA-01A-11R",
        "TCGA-A1-AAAA-01B-11R",
        "TCGA-A1-BBBB-01A-11R",
    ]
    got = one_per_patient(cols)
    assert got == ["TCGA-A1-AAAA-01A-11R", "TCGA-A1-BBBB-01A-11R"]


def test_sample_type():
    assert sample_type("TCGA-A1-AAAA-01A-11R") == "01"
    assert sample_type("TCGA-A1-AAAA-11A-11R") == "11"


def test_rank_partners_recovers_known_order():
    rng = np.random.default_rng(0)
    n = 80
    samples = [f"S{i:03d}" for i in range(n)]
    anchor = rng.normal(size=n)
    # gene A = noisy copy of anchor → should rank #1
    a = anchor + rng.normal(scale=0.05, size=n)
    # gene B = weaker copy
    b = anchor + rng.normal(scale=1.5, size=n)
    # gene C = independent
    c = rng.normal(size=n)
    expr = pd.DataFrame(
        {"TACSTD2": anchor, "GENE_A": a, "GENE_B": b, "GENE_C": c, "CLDN4": b},
        index=samples,
    ).T
    df = rank_partners(expr, ["GENE_A", "GENE_B", "GENE_C", "CLDN4"], samples)
    assert df.iloc[0]["gene"] == "GENE_A"
    assert set(df["gene"]) == {"GENE_A", "GENE_B", "GENE_C", "CLDN4"}
    assert df.loc[df["gene"] == "GENE_C", "spearman_r"].iloc[0] < df.loc[df["gene"] == "GENE_A", "spearman_r"].iloc[0]


def test_corr_and_fdr_shapes():
    mat = np.array([[1.0, 2, 3, 4], [4, 3, 2, 1], [1, 1, 1, 1]])
    vec = np.array([1.0, 2, 3, 4])
    r = _corr_rows_vs_vector(mat, vec)
    assert r.shape == (3,)
    assert r[0] > 0.99
    assert r[1] < -0.99
    p = np.array([0.001, 0.02, 0.8])
    q = _bh_fdr(p)
    assert q[0] <= q[1] <= q[2]


def test_partial_spearman_removes_shared_confounder():
    rng = np.random.default_rng(1)
    z = rng.normal(size=200)
    x = z + rng.normal(scale=0.2, size=200)
    y = z + rng.normal(scale=0.2, size=200)
    r_part, _ = partial_spearman(x, y, z)
    assert abs(r_part) < 0.25


if __name__ == "__main__":
    test_one_per_patient_keeps_first_vial()
    test_sample_type()
    test_rank_partners_recovers_known_order()
    test_corr_and_fdr_shapes()
    test_partial_spearman_removes_shared_confounder()
    print("ok")
