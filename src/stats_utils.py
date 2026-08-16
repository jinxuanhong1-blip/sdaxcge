"""Small statistics helpers shared by the A11 bulk and single-cell analyses."""

import numpy as np
import pandas as pd
from scipy import stats


def fisher_ci(r, n, alpha=0.05):
    """Fisher z confidence interval for a correlation coefficient."""
    if n < 5 or not np.isfinite(r) or abs(r) >= 1:
        return (np.nan, np.nan)
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    crit = stats.norm.ppf(1 - alpha / 2)
    return tuple(np.tanh([z - crit * se, z + crit * se]))


def spearman(x, y):
    """Spearman rho with n, p and 95% CI, on pairwise-complete observations."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(ok.sum())
    if n < 5:
        return dict(rho=np.nan, p=np.nan, n=n, lo=np.nan, hi=np.nan)
    r, p = stats.spearmanr(x[ok], y[ok])
    lo, hi = fisher_ci(r, n)
    return dict(rho=float(r), p=float(p), n=n, lo=lo, hi=hi)


def partial_spearman(x, y, covars):
    """Spearman correlation of x and y after linearly removing covariates.

    Ranks are taken first, then covariates (also ranked) are regressed out of
    both variables; the correlation of the residuals is reported. This is the
    standard rank-based partial correlation and is robust to the strongly
    non-normal expression distributions seen here.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    covars = np.column_stack([np.asarray(c, float) for c in covars])

    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(covars).all(axis=1)
    n = int(ok.sum())
    k = covars.shape[1]
    if n < 5 + k:
        return dict(rho=np.nan, p=np.nan, n=n, lo=np.nan, hi=np.nan)

    rx = stats.rankdata(x[ok])
    ry = stats.rankdata(y[ok])
    rc = np.column_stack([stats.rankdata(covars[ok, j]) for j in range(k)])
    design = np.column_stack([np.ones(n), rc])

    def resid(v):
        beta, *_ = np.linalg.lstsq(design, v, rcond=None)
        return v - design @ beta

    ex, ey = resid(rx), resid(ry)
    if ex.std() == 0 or ey.std() == 0:
        return dict(rho=np.nan, p=np.nan, n=n, lo=np.nan, hi=np.nan)

    r = float(np.corrcoef(ex, ey)[0, 1])
    df = n - 2 - k
    if abs(r) >= 1:
        p = 0.0
    else:
        t = r * np.sqrt(df / (1 - r ** 2))
        p = float(2 * stats.t.sf(abs(t), df))
    lo, hi = fisher_ci(r, n - k)
    return dict(rho=r, p=p, n=n, lo=lo, hi=hi)


def cliffs_delta(a, b):
    """Cliff's delta via the Mann-Whitney U statistic.

    delta = P(a > b) - P(a < b); ranges -1..1, 0 = no separation.
    Reported instead of Cohen's d because expression is not Gaussian.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return dict(delta=np.nan, p=np.nan, n1=len(a), n2=len(b))
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    delta = 2.0 * u / (len(a) * len(b)) - 1.0
    return dict(delta=float(delta), p=float(p), n1=len(a), n2=len(b))


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values; NaN p-values pass through as NaN."""
    p = np.asarray(pvals, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(pv)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    q[ok] = out
    return q


def zscore_module(df, genes):
    """Mean of per-gene z-scores across the samples present in df."""
    present = [g for g in genes if g in df.columns]
    if not present:
        return pd.Series(np.nan, index=df.index), []
    sub = df[present]
    z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, np.nan)
    return z.mean(axis=1), present


def stars(q):
    if not np.isfinite(q):
        return ""
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""
