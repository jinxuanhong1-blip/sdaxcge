"""Partial Spearman and DerSimonian–Laird helpers.

Partial Spearman is the Pearson correlation of rank residuals after OLS on
an intercept plus the rank-transformed covariates. Both sides are residualized.
The t reference uses df = n - 2 - k.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import stats


def rank_average(x: np.ndarray) -> np.ndarray:
    return stats.rankdata(np.asarray(x, dtype=float), method="average").astype(float)


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(rank_average(x), rank_average(y))


def partial_spearman(x: np.ndarray, y: np.ndarray, covariates: list[np.ndarray]) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    covs = [np.asarray(z, dtype=float) for z in covariates]
    n = int(x.size)
    k = len(covs)
    if n < k + 5:
        return float("nan")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        return float("nan")
    for z in covs:
        if z.shape != x.shape or not np.isfinite(z).all():
            return float("nan")
    xr = rank_average(x)
    yr = rank_average(y)
    if k == 0:
        return pearson(xr, yr)
    design = np.column_stack([np.ones(n)] + [rank_average(z) for z in covs])
    bx, *_ = np.linalg.lstsq(design, xr, rcond=None)
    by, *_ = np.linalg.lstsq(design, yr, rcond=None)
    rx = xr - design @ bx
    ry = yr - design @ by
    return pearson(rx, ry)


def spearman_p(rho: float, n: int, k_cov: int = 0) -> float:
    df = n - 2 - k_cov
    if not math.isfinite(rho) or df <= 0:
        return float("nan")
    if abs(rho) >= 1.0 - 1e-15:
        return 0.0
    tstat = rho * math.sqrt(df / (1.0 - rho * rho))
    return float(2 * stats.t.sf(abs(tstat), df))


def fisher_ci(rho: float, n: int, k: int) -> tuple[float, float]:
    if not math.isfinite(rho) or n - 3 - k <= 0:
        return float("nan"), float("nan")
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / math.sqrt(n - 3 - k)
    zcrit = 1.959963984540054
    return float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se))


def dl_meta(rhos: list[float], ns: list[int], k_cov: int = 0) -> dict:
    """DerSimonian–Laird random-effects pool on Fisher z. Var(z) = 1/(n-3-k)."""
    pairs = [
        (float(r), int(n))
        for r, n in zip(rhos, ns)
        if math.isfinite(r) and int(n) - 3 - k_cov > 1
    ]
    empty = {
        "rho": float("nan"),
        "p": float("nan"),
        "I2": float("nan"),
        "ci_lo": float("nan"),
        "ci_hi": float("nan"),
        "k": 0,
        "N": 0,
        "tau2": float("nan"),
    }
    if len(pairs) < 2:
        if len(pairs) == 1:
            r, n = pairs[0]
            lo, hi = fisher_ci(r, n, k_cov)
            empty.update(
                {
                    "rho": r,
                    "p": spearman_p(r, n, k_cov),
                    "I2": 0.0,
                    "ci_lo": lo,
                    "ci_hi": hi,
                    "k": 1,
                    "N": n,
                    "tau2": 0.0,
                }
            )
        return empty
    rhos_a = np.array([r for r, _ in pairs], dtype=float)
    ns_a = np.array([n for _, n in pairs], dtype=float)
    z = np.arctanh(np.clip(rhos_a, -0.999999, 0.999999))
    var_z = 1.0 / (ns_a - 3.0 - k_cov)
    w = 1.0 / var_z
    zbar = float(np.sum(w * z) / np.sum(w))
    q = float(np.sum(w * (z - zbar) ** 2))
    k = len(pairs)
    dfree = k - 1
    cdenom = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - dfree) / cdenom) if cdenom > 0 else 0.0
    wstar = 1.0 / (var_z + tau2)
    zre = float(np.sum(wstar * z) / np.sum(wstar))
    se = math.sqrt(1.0 / float(np.sum(wstar)))
    p = float(2 * stats.norm.sf(abs(zre / se)))
    i2 = max(0.0, (q - dfree) / q) if q > 0 else 0.0
    return {
        "rho": float(np.tanh(zre)),
        "p": p,
        "I2": float(i2),
        "ci_lo": float(np.tanh(zre - 1.959963984540054 * se)),
        "ci_hi": float(np.tanh(zre + 1.959963984540054 * se)),
        "k": k,
        "N": int(np.sum(ns_a)),
        "tau2": float(tau2),
    }


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q.tolist()
    pv = p[ok]
    m = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    qv = ranked * m / np.arange(1, m + 1)
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    qv = np.clip(qv, 0.0, 1.0)
    out = np.empty(m, dtype=float)
    out[order] = qv
    q[ok] = out
    return q.tolist()
