"""Honest multiple-testing helpers.

Families are declared before looking at p-values. BH is applied inside a
family; Bonferroni across the combinatorial primary grid is reported as a
sensitivity column, not as the sole decision rule.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(n, dtype=float)
    tmp[order] = q
    out[ok] = tmp
    return out


def add_fdr(df: pd.DataFrame, p_col: str = "p_value", family_col: str = "fdr_family") -> pd.DataFrame:
    out = df.copy()
    out["q_bh"] = np.nan
    out["n_tests_in_family"] = 0
    for fam, sub in out.groupby(family_col, dropna=False):
        idx = sub.index
        out.loc[idx, "q_bh"] = bh_fdr(sub[p_col].to_numpy())
        out.loc[idx, "n_tests_in_family"] = int(sub[p_col].notna().sum())
    # Bonferroni on the primary combinatorial grid only (filled later if present).
    if "bonferroni_primary" not in out.columns:
        out["bonferroni_primary"] = np.nan
    return out


def primary_recovery_flag(row: pd.Series, q_cut: float = 0.10) -> str:
    """T/NK or TLS down in TACSTD2-high. Requires negative effect and BH q < cut."""
    if row.get("compartment") not in {"TNK", "TLS"}:
        return "not_primary"
    if not np.isfinite(row.get("effect", np.nan)) or not np.isfinite(row.get("q_bh", np.nan)):
        return "not_estimable"
    if row["q_bh"] < q_cut and row["effect"] < 0:
        return "recovered_down"
    if row["q_bh"] < q_cut and row["effect"] > 0:
        return "opposite_up"
    return "not_recovered"
