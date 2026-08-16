"""Unit tests that do not need GEO/TCGA downloads."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import pytest
except ImportError:  # allow running the file as a script
    class pytest:  # type: ignore
        @staticmethod
        def raises(exc, match=None):
            class _Ctx:
                def __enter__(self):
                    return self
                def __exit__(self, et, ev, tb):
                    if et is None:
                        raise AssertionError(f"did not raise {exc}")
                    if not issubclass(et, exc):
                        return False
                    if match and match not in str(ev):
                        raise AssertionError(f"{ev} did not match {match!r}")
                    return True
            return _Ctx()

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bulkimmune.genes import collapse_duplicates, strip_ensembl_version
from bulkimmune.preprocess import cpm, detect_scale, looks_log_scale
from bulkimmune.ssgsea import rank_matrix, ssgsea
from bulkimmune.stats import adjust_pvalues, spearman_table


def test_rank_average_int_matches_r_truncation():
    expr = pd.DataFrame({"s1": [1.0, 1.0, 3.0, 4.0]}, index=list("abcd"))
    ranks = rank_matrix(expr, tie_method="average_int")
    # average ranks of the two 1.0s are 1.5, 1.5 -> trunc -> 1, 1
    assert list(ranks[:, 0]) == [1.0, 1.0, 3.0, 4.0]


def test_ssgsea_monotone_within_sample():
    rng = np.random.default_rng(0)
    expr = pd.DataFrame(rng.normal(size=(80, 5)), index=[f"g{i}" for i in range(80)],
                        columns=[f"s{i}" for i in range(5)])
    sets = {"A": [f"g{i}" for i in range(10)]}
    a = ssgsea(expr, sets, normalize="none")
    b = ssgsea(np.log2(expr - expr.min().min() + 1), sets, normalize="none")
    # strictly monotone per-column transform must leave ranks (hence scores) unchanged
    np.testing.assert_allclose(a.to_numpy(), b.to_numpy(), rtol=0, atol=1e-8)


def test_ssgsea_empty_overlap_dropped():
    expr = pd.DataFrame(np.arange(20).reshape(10, 2), index=[f"g{i}" for i in range(10)],
                        columns=["s1", "s2"])
    with pytest.raises(ValueError, match="no gene set"):
        ssgsea(expr, {"ghost": ["NOT_A_GENE"]}, min_size=1)


def test_cpm_and_scale_detection():
    counts = pd.DataFrame({"s1": [10, 0, 90], "s2": [20, 20, 60]}, index=list("abc"))
    out = cpm(counts)
    assert np.allclose(out.sum(), [1e6, 1e6])
    assert "1e6" in detect_scale(out)["guess"] or detect_scale(out)["guess"].startswith("TPM")
    assert looks_log_scale(np.log2(out + 1))


def test_collapse_max_mean():
    expr = pd.DataFrame({"s1": [1.0, 5.0, 2.0]}, index=["A", "A", "B"])
    out = collapse_duplicates(expr, method="max_mean")
    assert list(out.index) == ["A", "B"]
    assert out.loc["A", "s1"] == 5.0


def test_strip_ensembl():
    assert strip_ensembl_version(["ENSG0001.10", "TACSTD2"]) == ["ENSG0001", "TACSTD2"]


def test_spearman_and_fdr_families():
    rng = np.random.default_rng(1)
    n = 40
    t = pd.DataFrame({"TACSTD2": rng.normal(size=n)}, index=[f"p{i}" for i in range(n)])
    s = pd.DataFrame(
        {
            "Exclusion": t["TACSTD2"] + rng.normal(scale=0.2, size=n),
            "IFNG_Ayers6": rng.normal(size=n),
            "something_else": rng.normal(size=n),
        },
        index=t.index,
    )
    table = spearman_table(t, s)
    table["family"] = table["score"].map(
        lambda x: "exclusion" if x == "Exclusion" else "inflamed" if x == "IFNG_Ayers6" else "other"
    )
    adj = adjust_pvalues(table, group_cols=("target",), family_col="family")
    excl = adj[adj["score"] == "Exclusion"].iloc[0]
    assert excl["spearman_r"] > 0.8
    assert excl["p_adj"] < 0.05
