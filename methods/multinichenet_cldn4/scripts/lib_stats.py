"""Patient-level tests. Cells are never the inferential unit."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def mannwhitney(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    na, nb = len(a), len(b)
    if na == 0 or nb == 0:
        return {
            "n_a": na,
            "n_b": nb,
            "mean_a": np.nan,
            "mean_b": np.nan,
            "median_a": np.nan,
            "median_b": np.nan,
            "U": np.nan,
            "p": np.nan,
            "rank_biserial": np.nan,
        }
    res = stats.mannwhitneyu(a, b, alternative="two-sided", method="auto")
    # r = 1 - 2U_a,b/(n_a n_b); r < 0 when group b (Q4 / high) is lower than a (Q1 / low)
    r = 1.0 - (2.0 * float(res.statistic) / (na * nb))
    return {
        "n_a": na,
        "n_b": nb,
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "U": float(res.statistic),
        "p": float(res.pvalue),
        "rank_biserial": float(r),
    }


def wilcoxon_paired(high: np.ndarray, low: np.ndarray) -> dict:
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    m = np.isfinite(high) & np.isfinite(low)
    high, low = high[m], low[m]
    n = int(len(high))
    delta = high - low
    if n < 3 or np.allclose(high, low):
        return {
            "n": n,
            "median_delta": float(np.median(delta)) if n else np.nan,
            "mean_delta": float(np.mean(delta)) if n else np.nan,
            "p": np.nan,
        }
    try:
        p = float(stats.wilcoxon(high, low, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        p = np.nan
    return {
        "n": n,
        "median_delta": float(np.median(delta)),
        "mean_delta": float(np.mean(delta)),
        "p": p,
    }


def spearman_safe(x: pd.Series | np.ndarray, y: pd.Series | np.ndarray) -> dict:
    x = pd.Series(x, dtype=float)
    y = pd.Series(y, dtype=float)
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


def fisher_z_pool(rows: list[dict]) -> dict:
    """Inverse-variance Fisher-z pool of Spearman ρ. Honest n = sum of patient n."""
    use = [r for r in rows if np.isfinite(r.get("rho", np.nan)) and r.get("n", 0) >= 4]
    if not use:
        return {"k": 0, "N": 0, "rho": np.nan, "p": np.nan, "I2": np.nan}
    zs, ws, ns = [], [], []
    for r in use:
        n = int(r["n"])
        rho = float(np.clip(r["rho"], -0.999999, 0.999999))
        z = np.arctanh(rho)
        w = n - 3
        zs.append(z)
        ws.append(w)
        ns.append(n)
    zs = np.asarray(zs, dtype=float)
    ws = np.asarray(ws, dtype=float)
    zbar = float(np.sum(ws * zs) / np.sum(ws))
    se = float(1.0 / np.sqrt(np.sum(ws)))
    p = float(2 * stats.norm.sf(abs(zbar / se)))
    # DerSimonian–Laird I² on the z scale
    q = float(np.sum(ws * (zs - zbar) ** 2))
    df = len(use) - 1
    i2 = float(max(0.0, (q - df) / q) * 100) if q > 0 and df > 0 else 0.0
    return {
        "k": len(use),
        "N": int(sum(ns)),
        "rho": float(np.tanh(zbar)),
        "p": p,
        "I2": i2,
    }


def quartile_tails(values: pd.Series) -> tuple[pd.Index, pd.Index]:
    v = values.dropna()
    q1 = v.quantile(0.25)
    q3 = v.quantile(0.75)
    low = v[v <= q1].index
    high = v[v >= q3].index
    return low, high


def rank_biserial_q4q1(score: pd.Series, group_key: pd.Series) -> dict:
    low_idx, high_idx = quartile_tails(group_key)
    a = score.reindex(low_idx).dropna().to_numpy()
    b = score.reindex(high_idx).dropna().to_numpy()
    out = mannwhitney(a, b)
    out["n_Q1"] = out.pop("n_a")
    out["n_Q4"] = out.pop("n_b")
    out["median_Q1"] = out.pop("median_a")
    out["median_Q4"] = out.pop("median_b")
    out["mean_Q1"] = out.pop("mean_a")
    out["mean_Q4"] = out.pop("mean_b")
    out["n_compared"] = int(out["n_Q1"] + out["n_Q4"])
    return out


def minmax_scale(s: pd.Series) -> pd.Series:
    x = s.astype(float)
    lo, hi = x.min(), x.max()
    if not np.isfinite(lo) or not np.isfinite(hi) or np.isclose(hi, lo):
        return pd.Series(np.nan, index=s.index)
    return (x - lo) / (hi - lo)


def fdr_bh(p: pd.Series) -> pd.Series:
    from statsmodels.stats.multitest import multipletests

    p = p.astype(float)
    out = pd.Series(np.nan, index=p.index)
    m = p.notna()
    if m.sum() == 0:
        return out
    out.loc[m] = multipletests(p[m], method="fdr_bh")[1]
    return out
