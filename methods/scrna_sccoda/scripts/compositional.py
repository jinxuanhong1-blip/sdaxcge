"""Compositional statistics: CLR, ILR, Dirichlet-multinomial.

scCODA was not installable in this environment (pip pulls rpy2, which
requires a system R). Tests below are the pre-specified fallbacks.
Patient/sample is the unit. Cells are never treated as replicates.
"""
from __future__ import annotations

from itertools import combinations
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from scipy.special import gammaln


LINEAGE_COLS = ["Epithelial", "T_NK", "B", "Myeloid", "Stromal", "Other"]
FOCUS = ("T_NK", "Epithelial")


def clr_matrix(counts: np.ndarray, pc: float = 0.5) -> np.ndarray:
    """Centered log-ratio of compositions. counts: n x D non-negative."""
    x = np.asarray(counts, dtype=float) + pc
    x = x / x.sum(axis=1, keepdims=True)
    g = np.exp(np.mean(np.log(x), axis=1, keepdims=True))
    return np.log(x / g)


def ilr_one_vs_rest(counts: np.ndarray, idx: int, pc: float = 0.5) -> np.ndarray:
    """Pivot ILR coordinate of part `idx` versus the geometric mean of the rest."""
    x = np.asarray(counts, dtype=float) + pc
    x = x / x.sum(axis=1, keepdims=True)
    rest = np.delete(x, idx, axis=1)
    g_rest = np.exp(np.mean(np.log(rest), axis=1))
    d = x.shape[1]
    return np.sqrt((d - 1) / d) * np.log(x[:, idx] / g_rest)


def fractions(counts: np.ndarray) -> np.ndarray:
    x = np.asarray(counts, dtype=float)
    n = x.sum(axis=1, keepdims=True)
    n[n == 0] = np.nan
    return x / n


def _exact_mwu_p(a: np.ndarray, b: np.ndarray, alternative: str) -> float:
    """Exact two-sample Wilcoxon by enumerating assignments when n is small."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return float("nan")
    if n1 + n2 > 16 or (n1 + n2 > 12 and n1 not in (n2, n2 - 1, n2 + 1) and n1 * n2 > 84):
        return float(stats.mannwhitneyu(a, b, alternative=alternative).pvalue)
    pooled = np.concatenate([a, b])
    obs = stats.mannwhitneyu(a, b, alternative=alternative).statistic
    labels = np.zeros(n1 + n2, dtype=int)
    labels[:n1] = 1
    extreme = 0
    total = 0
    for idx in combinations(range(n1 + n2), n1):
        mask = np.zeros(n1 + n2, dtype=bool)
        mask[list(idx)] = True
        u = stats.mannwhitneyu(pooled[mask], pooled[~mask], alternative=alternative).statistic
        total += 1
        if alternative == "two-sided":
            # compare two-sided p via |U - null|
            null = n1 * n2 / 2.0
            if abs(u - null) + 1e-12 >= abs(obs - null):
                extreme += 1
        elif alternative == "less":
            if u <= obs + 1e-12:
                extreme += 1
        else:  # greater
            if u >= obs - 1e-12:
                extreme += 1
    return extreme / total


def mwu(a: Iterable[float], b: Iterable[float], alternative: str = "two-sided") -> dict:
    a = np.asarray(list(a), dtype=float)
    b = np.asarray(list(b), dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(np.median(a)) if len(a) else float("nan"),
        "median_b": float(np.median(b)) if len(b) else float("nan"),
        "mean_a": float(np.mean(a)) if len(a) else float("nan"),
        "mean_b": float(np.mean(b)) if len(b) else float("nan"),
        "delta_median_a_minus_b": float("nan"),
        "U": float("nan"),
        "p": float("nan"),
        "p_method": "none",
        "alternative": alternative,
    }
    if len(a) < 2 or len(b) < 2:
        return out
    u = float(stats.mannwhitneyu(a, b, alternative=alternative).statistic)
    p_asymp = float(stats.mannwhitneyu(a, b, alternative=alternative).pvalue)
    if len(a) + len(b) <= 16:
        p = _exact_mwu_p(a, b, alternative)
        method = "exact_enumeration"
    else:
        p = p_asymp
        method = "asymptotic_mwu"
    out.update(
        {
            "U": u,
            "p": float(p),
            "p_asymptotic": p_asymp,
            "p_method": method,
            "delta_median_a_minus_b": float(np.median(a) - np.median(b)),
        }
    )
    return out


def spearman(x: Iterable[float], y: Iterable[float]) -> dict:
    x = np.asarray(list(x), dtype=float)
    y = np.asarray(list(y), dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    out = {"n": int(len(x)), "rho": float("nan"), "p": float("nan")}
    if len(x) < 4:
        return out
    rho, p = stats.spearmanr(x, y)
    out["rho"] = float(rho)
    out["p"] = float(p)
    return out


def dm_loglik(alpha: np.ndarray, counts: np.ndarray) -> float:
    alpha = np.asarray(alpha, dtype=float)
    counts = np.asarray(counts, dtype=float)
    A = float(alpha.sum())
    N = counts.sum(axis=1)
    ll = np.sum(gammaln(A) - gammaln(A + N))
    ll += np.sum(gammaln(counts + alpha) - gammaln(alpha))
    return float(ll)


def fit_dm(counts: np.ndarray) -> tuple[np.ndarray, float]:
    counts = np.asarray(counts, dtype=float)
    p = counts / np.maximum(counts.sum(axis=1, keepdims=True), 1.0)
    mu = np.clip(p.mean(axis=0), 1e-6, 1.0)
    start = np.log(mu * max(float(counts.sum() / max(len(counts), 1) / 10.0), 1.0))

    def nll(log_a):
        a = np.exp(log_a)
        if np.any(~np.isfinite(a)) or np.any(a <= 0):
            return 1e12
        return -dm_loglik(a, counts)

    res = minimize(nll, start, method="L-BFGS-B")
    alpha = np.exp(res.x)
    return alpha, dm_loglik(alpha, counts)


def dm_two_group(counts: np.ndarray, group: np.ndarray, n_perm: int = 2000, seed: int = 1) -> dict:
    """Dirichlet-multinomial two-group LRT with permutation p.

    H0: one shared alpha. H1: group-specific alpha. Statistic = 2 (ll1-ll0).
    """
    counts = np.asarray(counts, dtype=float)
    group = np.asarray(group)
    g0, g1 = np.unique(group)
    if len(g0) == 0:
        return {"note": "empty"}
    if len(np.unique(group)) != 2:
        return {"note": "need exactly two groups", "n": int(len(group))}
    a0 = group == g0
    a1 = group == g1
    if a0.sum() < 2 or a1.sum() < 2:
        return {"note": "need >=2 per group", "n0": int(a0.sum()), "n1": int(a1.sum())}

    _, ll0 = fit_dm(counts)
    _, ll_a = fit_dm(counts[a0])
    _, ll_b = fit_dm(counts[a1])
    ll1 = ll_a + ll_b
    stat = 2.0 * (ll1 - ll0)
    df = counts.shape[1]
    p_chi = float(stats.chi2.sf(max(stat, 0.0), df))

    rng = np.random.default_rng(seed)
    extreme = 0
    n_ok = 0
    labels = group.copy()
    for _ in range(n_perm):
        rng.shuffle(labels)
        m0 = labels == g0
        m1 = labels == g1
        try:
            _, p_ll0 = fit_dm(counts)
            _, p_lla = fit_dm(counts[m0])
            _, p_llb = fit_dm(counts[m1])
            pstat = 2.0 * (p_lla + p_llb - p_ll0)
        except Exception:
            continue
        n_ok += 1
        if pstat + 1e-9 >= stat:
            extreme += 1
    p_perm = (extreme + 1) / (n_ok + 1) if n_ok else float("nan")
    return {
        "n0": int(a0.sum()),
        "n1": int(a1.sum()),
        "group0": str(g0),
        "group1": str(g1),
        "lrt_stat": float(stat),
        "df": int(df),
        "p_chi2": p_chi,
        "p_perm": float(p_perm),
        "n_perm": int(n_ok),
        "ll_null": float(ll0),
        "ll_alt": float(ll1),
    }


def attach_transforms(df: pd.DataFrame, count_cols: list[str] | None = None, pc: float = 0.5) -> pd.DataFrame:
    count_cols = count_cols or LINEAGE_COLS
    out = df.copy()
    C = out[count_cols].to_numpy(dtype=float)
    F = fractions(C)
    L = clr_matrix(C, pc=pc)
    for i, name in enumerate(count_cols):
        out[f"frac_{name}"] = F[:, i]
        out[f"clr_{name}"] = L[:, i]
        out[f"ilr_{name}"] = ilr_one_vs_rest(C, i, pc=pc)
    out["n_cells_comp"] = C.sum(axis=1)
    return out


def cohort_center(df: pd.DataFrame, cols: list[str], by: str = "cohort") -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        out[f"{c}__cc"] = out[c] - out.groupby(by)[c].transform("median")
    return out
