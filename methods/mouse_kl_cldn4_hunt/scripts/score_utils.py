"""Shared mouse-level scoring helpers. CLDN4-only. Honest n."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats


def present_genes(available: Iterable[str], wanted: list[str]) -> list[str]:
    have = set(available)
    return [g for g in wanted if g in have]


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Z-score genes (rows) across samples. Constant genes become 0."""
    mu = df.mean(axis=1)
    sd = df.std(axis=1, ddof=0)
    sd = sd.replace(0, np.nan)
    out = df.sub(mu, axis=0).div(sd, axis=0)
    return out.fillna(0.0)


def mean_score(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Mean of present genes; NaN if none present."""
    use = [g for g in genes if g in expr.index]
    if not use:
        return pd.Series(np.nan, index=expr.columns, name="score")
    return expr.loc[use].mean(axis=0)


def spearman_safe(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    d = pd.concat([x, y], axis=1).dropna()
    n = int(len(d))
    if n < 3:
        return (float("nan"), float("nan"), n)
    r, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    return (float(r), float(p), n)


def median_split_test(score: pd.Series, group_on: pd.Series, label_high="Cldn4_high") -> dict:
    """Compare score in Cldn4-high vs low mice (median of group_on). Honest n."""
    d = pd.concat({"x": group_on, "y": score}, axis=1).dropna()
    n = int(len(d))
    out = {
        "n_mice": n,
        "n_high": 0,
        "n_low": 0,
        "median_split_value": float("nan"),
        "mean_high": float("nan"),
        "mean_low": float("nan"),
        "delta_high_minus_low": float("nan"),
        "mannwhitney_u": float("nan"),
        "mannwhitney_p": float("nan"),
        "welch_t": float("nan"),
        "welch_p": float("nan"),
        "note": "",
    }
    if n < 4:
        out["note"] = f"n={n} too small for median-split test"
        return out
    # Zero-inflated Cldn4 (common in bulk LLC): detected vs undetected.
    n_zero = int((d["x"] == 0).sum())
    n_pos = int((d["x"] > 0).sum())
    if n_zero >= 2 and n_pos >= 2 and n_zero >= max(2, n // 4):
        hi = d.loc[d["x"] > 0, "y"]
        lo = d.loc[d["x"] == 0, "y"]
        med = 0.0
        out["note"] = "zero-inflated: Cldn4 detected vs undetected"
    else:
        med = float(d["x"].median())
        hi = d.loc[d["x"] >= med, "y"]
        lo = d.loc[d["x"] < med, "y"]
        # If ties at median empty one arm, use strictly greater / less and drop equals.
        if len(hi) == 0 or len(lo) == 0:
            hi = d.loc[d["x"] > med, "y"]
            lo = d.loc[d["x"] < med, "y"]
            out["note"] = "ties at median dropped"
    out["n_high"] = int(len(hi))
    out["n_low"] = int(len(lo))
    out["median_split_value"] = med
    if len(hi) == 0 or len(lo) == 0:
        out["note"] = (out["note"] + "; one arm empty after split").strip("; ")
        return out
    out["mean_high"] = float(hi.mean())
    out["mean_low"] = float(lo.mean())
    out["delta_high_minus_low"] = out["mean_high"] - out["mean_low"]
    try:
        u, p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
        out["mannwhitney_u"] = float(u)
        out["mannwhitney_p"] = float(p)
    except ValueError:
        out["note"] = (out["note"] + "; MWU failed").strip("; ")
    if len(hi) >= 2 and len(lo) >= 2 and hi.var(ddof=1) + lo.var(ddof=1) > 0:
        t, p = stats.ttest_ind(hi, lo, equal_var=False)
        out["welch_t"] = float(t)
        out["welch_p"] = float(p)
    return out


def association_table(mouse_df: pd.DataFrame, cldn4_col: str, endpoints: list[str]) -> pd.DataFrame:
    rows = []
    for ep in endpoints:
        r, p, n = spearman_safe(mouse_df[cldn4_col], mouse_df[ep])
        split = median_split_test(mouse_df[ep], mouse_df[cldn4_col])
        rows.append(
            {
                "endpoint": ep,
                "n_mice_spearman": n,
                "spearman_r": r,
                "spearman_p": p,
                **{f"split_{k}": v for k, v in split.items()},
            }
        )
    return pd.DataFrame(rows)


def fmt_p(p: float) -> str:
    if p is None or (isinstance(p, float) and (math.isnan(p) or math.isinf(p))):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"
