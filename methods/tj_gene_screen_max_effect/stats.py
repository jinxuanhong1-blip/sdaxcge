"""Rank-based partial correlation helpers for the TCGA Trop2 analysis.

Partial Spearman is Pearson correlation of rank residuals after ordinary
least squares on an intercept plus the rank-transformed covariates.
The t reference distribution uses df = n - 2 - k, matching the
purity-partial Spearman used in the earlier Xena LUAD script.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def bh_fdr(p):
    """Benjamini-Hochberg q-values. NaN inputs stay NaN and are not counted."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
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
    return q


def _partial_corr(x, y, covariates, rank: bool):
    """Partial correlation. rank=True is partial Spearman; rank=False is partial Pearson."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    covs = [np.asarray(z, dtype=float) for z in covariates]
    mask = np.isfinite(x) & np.isfinite(y)
    for z in covs:
        mask &= np.isfinite(z)
    x = x[mask]
    y = y[mask]
    covs = [z[mask] for z in covs]
    n = int(x.size)
    k = len(covs)
    if n < k + 5:
        return np.nan, np.nan, n
    if rank:
        x = stats.rankdata(x).astype(float)
        y = stats.rankdata(y).astype(float)
        covs = [stats.rankdata(z).astype(float) for z in covs]
    if k:
        design = np.column_stack([np.ones(n)] + covs)
        bx, *_ = np.linalg.lstsq(design, x, rcond=None)
        by, *_ = np.linalg.lstsq(design, y, rcond=None)
        rx = x - design @ bx
        ry = y - design @ by
    else:
        rx, ry = x, y
    if np.std(rx) < 1e-12 or np.std(ry) < 1e-12:
        return np.nan, np.nan, n
    r, _ = stats.pearsonr(rx, ry)
    df = n - 2 - k
    r = float(r)
    if df <= 0 or not np.isfinite(r):
        return r, np.nan, n
    if abs(r) >= 1.0 - 1e-15:
        return r, 0.0, n
    tstat = r * np.sqrt(df / (1.0 - r * r))
    p = float(2 * stats.t.sf(abs(tstat), df))
    return r, p, n


def partial_spearman(x, y, covariates=()):
    """Spearman rho of x vs y controlling for covariates.

    Returns rho, two-sided p, and the number of complete observations.
    """
    return _partial_corr(x, y, covariates, rank=True)


def partial_pearson(x, y, covariates=()):
    """Pearson partial correlation on the supplied values (no rank transform)."""
    return _partial_corr(x, y, covariates, rank=False)


def fisher_ci(rho, n, k, alpha=0.05):
    """Fisher-z confidence interval. k is the number of covariates."""
    if not np.isfinite(rho) or n - 3 - k <= 0:
        return np.nan, np.nan
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3 - k)
    zcrit = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se))


def spearman(x, y):
    """Marginal Spearman. Returns rho, two-sided p."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if x.size < 5 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan"), float("nan")
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p)


def dl_meta_spearman(rhos, ns, k_cov=0):
    """DL meta on Spearman rhos. Returns dict with ci_lo/ci_hi aliases."""
    out = random_effects_meta(rhos, ns, k_cov)
    return {
        "n_cohorts": out["n_cohorts"],
        "rho": out["rho"],
        "p": out["p"],
        "I2": out["I2"],
        "ci_lo": out["ci_low"],
        "ci_hi": out["ci_high"],
        "tau2": out["tau2"],
        "Q": out["Q"],
        "Q_p": out["Q_p"],
    }


def random_effects_meta(rhos, ns, k):
    """DerSimonian-Laird random-effects meta-analysis on Fisher z.

    Variance of z is 1/(n - 3 - k). Cohorts with non-finite rho or
    n - 3 - k <= 1 are dropped.
    """
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    keep = np.isfinite(rhos) & np.isfinite(ns) & ((ns - 3 - k) > 1)
    rhos = rhos[keep]
    ns = ns[keep]
    m = int(rhos.size)
    empty = {
        "n_cohorts": 0,
        "rho": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "p": np.nan,
        "I2": np.nan,
        "tau2": np.nan,
        "Q": np.nan,
        "Q_p": np.nan,
    }
    if m < 2:
        return empty
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    se2 = 1.0 / (ns - 3 - k)
    w = 1.0 / se2
    z_fe = np.sum(w * z) / np.sum(w)
    Q = float(np.sum(w * (z - z_fe) ** 2))
    df = m - 1
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (se2 + tau2)
    z_re = float(np.sum(w_re * z) / np.sum(w_re))
    se_re = float(np.sqrt(1.0 / np.sum(w_re)))
    I2 = max(0.0, (Q - df) / Q) if Q > 0 else 0.0
    p = float(2 * stats.norm.sf(abs(z_re / se_re)))
    zcrit = stats.norm.ppf(0.975)
    return {
        "n_cohorts": m,
        "rho": float(np.tanh(z_re)),
        "ci_low": float(np.tanh(z_re - zcrit * se_re)),
        "ci_high": float(np.tanh(z_re + zcrit * se_re)),
        "p": p,
        "I2": float(I2),
        "tau2": float(tau2),
        "Q": Q,
        "Q_p": float(stats.chi2.sf(Q, df)),
    }
