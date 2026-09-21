"""Effect-size helpers for the TCGA KRT5/6 specification sweep.

Partial Spearman is Pearson correlation of rank residuals after OLS on an
intercept plus the rank-transformed covariates. Partial Pearson skips the
rank transform. Degrees of freedom are n - 2 - k, with k the number of
covariate columns.

Cliff's delta is P(high > low) - P(high < low). Negative means the upper
predictor group has the lower immune score. Its sampling variance is the
Cliff 1993 estimator. Random-effects pooling of deltas is DerSimonian-Laird
on the delta scale. Random-effects pooling of correlations is
DerSimonian-Laird on the Fisher z scale, with z variance 1/(n - 3 - k).
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def bh_fdr(p):
    """Benjamini-Hochberg q-values. Non-finite inputs stay non-finite and are not counted."""
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


def residualize(values, covariates):
    """OLS residuals of each column on an intercept plus covariates.

    values: (n,) or (n, p). covariates: list of (n,) arrays, or an empty list.
    """
    values = np.asarray(values, dtype=float)
    single = values.ndim == 1
    if single:
        values = values[:, None]
    n = values.shape[0]
    cols = [np.asarray(z, dtype=float) for z in covariates]
    if cols:
        design = np.column_stack([np.ones(n)] + cols)
        beta, *_ = np.linalg.lstsq(design, values, rcond=None)
        out = values - design @ beta
    else:
        out = values.copy()
    return out[:, 0] if single else out


def _corr_p(rx, ry, df):
    if np.std(rx) < 1e-12 or np.std(ry) < 1e-12:
        return np.nan, np.nan
    r = float(np.corrcoef(rx, ry)[0, 1])
    if df <= 0 or not np.isfinite(r):
        return r, np.nan
    if abs(r) >= 1.0 - 1e-15:
        return r, 0.0
    tstat = r * np.sqrt(df / (1.0 - r * r))
    p = float(2 * stats.t.sf(abs(tstat), df))
    return r, p


def partial_corr(x, y, covariates, rank: bool):
    """Partial correlation. Returns rho, two-sided p, and complete-case n."""
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
    rx = residualize(x, covs)
    ry = residualize(y, covs)
    r, p = _corr_p(rx, ry, n - 2 - k)
    return r, p, n


def matrix_partial_corr(predictors, outcomes, covariates, df):
    """Pearson correlations of residual columns.

    predictors, outcomes, and covariates are already on the analysis scale
    (ranks, winsorized values, or raw). Returns rho (n_pred, n_out) and p.
    """
    pred = residualize(predictors, covariates)
    out = residualize(outcomes, covariates)
    pred = pred - pred.mean(axis=0)
    out = out - out.mean(axis=0)
    sp = np.sqrt(np.sum(pred * pred, axis=0))
    so = np.sqrt(np.sum(out * out, axis=0))
    denom = sp[:, None] * so[None, :]
    numer = pred.T @ out
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = numer / denom
    rho[(sp < 1e-12)[:, None] | (so < 1e-12)[None, :]] = np.nan
    p = np.full(rho.shape, np.nan)
    if df <= 0:
        return rho, p
    finite = np.isfinite(rho)
    r = np.clip(rho, -1.0 + 1e-15, 1.0 - 1e-15)
    tstat = r * np.sqrt(df / (1.0 - r * r))
    p[finite] = 2 * stats.t.sf(np.abs(tstat[finite]), df)
    p[finite & (np.abs(rho) >= 1.0 - 1e-15)] = 0.0
    return rho, p


def cliffs_delta(high, low):
    """Cliff's delta of high versus low, plus Cliff 1993 variance and a z p-value.

    delta = P(high > low) - P(high < low). Ties contribute 0.
    Returns delta, two-sided p, variance, n_high, n_low.
    """
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    high = high[np.isfinite(high)]
    low = low[np.isfinite(low)]
    n1 = int(high.size)
    n2 = int(low.size)
    if n1 < 2 or n2 < 2:
        return np.nan, np.nan, np.nan, n1, n2
    dom = np.sign(high[:, None] - low[None, :])
    delta = float(dom.mean())
    di = dom.mean(axis=1)
    dj = dom.mean(axis=0)
    s1 = float(np.sum((di - delta) ** 2) / (n1 - 1))
    s2 = float(np.sum((dj - delta) ** 2) / (n2 - 1))
    inter = dom - di[:, None] - dj[None, :] + delta
    s12 = float(np.sum(inter ** 2) / ((n1 - 1) * (n2 - 1)))
    var = ((n1 - 1) * s1 + (n2 - 1) * s2 + s12) / (n1 * n2)
    null_var = (n1 + n2 + 1) / (3.0 * n1 * n2)
    if not np.isfinite(var) or var <= 0:
        var = float(null_var)
    se = np.sqrt(var)
    if se <= 0 or not np.isfinite(delta):
        p = np.nan
    elif abs(delta) < 1e-15:
        p = 1.0
    else:
        p = float(2 * stats.norm.sf(abs(delta) / se))
    return delta, p, float(var), n1, n2


def cliffs_batch(high, low):
    """Cliff's delta for every column of high and low.

    high: (n1, p), low: (n2, p). Columns with a non-finite value come back non-finite.
    Returns delta, p, var, each shape (p,).
    """
    high = np.asarray(high, dtype=float)
    low = np.asarray(low, dtype=float)
    if high.ndim == 1:
        high = high[:, None]
        low = low[:, None]
    n1, p = high.shape
    n2 = low.shape[0]
    delta = np.full(p, np.nan)
    pval = np.full(p, np.nan)
    var = np.full(p, np.nan)
    if n1 < 2 or n2 < 2:
        return delta, pval, var
    ok = np.isfinite(high).all(axis=0) & np.isfinite(low).all(axis=0)
    if not np.any(ok):
        return delta, pval, var
    h = high[:, ok]
    l = low[:, ok]
    dom = np.sign(h[:, None, :] - l[None, :, :])
    d = dom.mean(axis=(0, 1))
    di = dom.mean(axis=1)
    dj = dom.mean(axis=0)
    s1 = np.sum((di - d) ** 2, axis=0) / (n1 - 1)
    s2 = np.sum((dj - d) ** 2, axis=0) / (n2 - 1)
    inter = dom - di[:, None, :] - dj[None, :, :] + d
    s12 = np.sum(inter ** 2, axis=(0, 1)) / ((n1 - 1) * (n2 - 1))
    v = ((n1 - 1) * s1 + (n2 - 1) * s2 + s12) / (n1 * n2)
    null_var = (n1 + n2 + 1) / (3.0 * n1 * n2)
    v = np.where(np.isfinite(v) & (v > 0), v, null_var)
    se = np.sqrt(v)
    pz = np.ones(d.shape, dtype=float)
    nonzero = se > 0
    pz[nonzero] = 2 * stats.norm.sf(np.abs(d[nonzero]) / se[nonzero])
    delta[ok] = d
    var[ok] = v
    pval[ok] = pz
    return delta, pval, var


def dl_random_effects(effects, variances):
    """DerSimonian-Laird random-effects pool on the effect scale."""
    effects = np.asarray(effects, dtype=float)
    variances = np.asarray(variances, dtype=float)
    keep = np.isfinite(effects) & np.isfinite(variances) & (variances > 0)
    effects = effects[keep]
    variances = variances[keep]
    m = int(effects.size)
    empty = {
        "n_cohorts": 0,
        "effect": np.nan,
        "fe": np.nan,
        "se": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "p": np.nan,
        "I2": np.nan,
        "tau2": np.nan,
        "Q": np.nan,
    }
    if m < 2:
        return empty
    w = 1.0 / variances
    fe = float(np.sum(w * effects) / np.sum(w))
    q = float(np.sum(w * (effects - fe) ** 2))
    df = m - 1
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (variances + tau2)
    re = float(np.sum(w_re * effects) / np.sum(w_re))
    se = float(np.sqrt(1.0 / np.sum(w_re)))
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    p = float(2 * stats.norm.sf(abs(re / se))) if se > 0 else np.nan
    zcrit = stats.norm.ppf(0.975)
    return {
        "n_cohorts": m,
        "effect": re,
        "fe": fe,
        "se": se,
        "ci_low": float(re - zcrit * se),
        "ci_high": float(re + zcrit * se),
        "p": p,
        "I2": float(i2),
        "tau2": float(tau2),
        "Q": q,
    }


def fisher_random_effects(rhos, ns, k):
    """DerSimonian-Laird random-effects pool of correlations on the Fisher z scale.

    k is the covariate count, shared by every cohort in the specification.
    Cohorts with non-finite rho or n - 3 - k <= 1 are dropped.
    """
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    keep = np.isfinite(rhos) & np.isfinite(ns) & ((ns - 3 - k) > 1)
    rhos = rhos[keep]
    ns = ns[keep]
    m = int(rhos.size)
    empty = {
        "n_cohorts": 0,
        "effect": np.nan,
        "fe": np.nan,
        "se": np.nan,
        "ci_low": np.nan,
        "ci_high": np.nan,
        "p": np.nan,
        "I2": np.nan,
        "tau2": np.nan,
        "Q": np.nan,
    }
    if m < 2:
        return empty
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    variances = 1.0 / (ns - 3 - k)
    pooled = dl_random_effects(z, variances)
    if pooled["n_cohorts"] == 0 or not np.isfinite(pooled["effect"]):
        return empty

    def back(z_value):
        return float(np.tanh(z_value))

    return {
        "n_cohorts": pooled["n_cohorts"],
        "effect": back(pooled["effect"]),
        "fe": back(pooled["fe"]),
        "se": pooled["se"],
        "ci_low": back(pooled["ci_low"]),
        "ci_high": back(pooled["ci_high"]),
        "p": pooled["p"],
        "I2": pooled["I2"],
        "tau2": pooled["tau2"],
        "Q": pooled["Q"],
    }


def fisher_ci(rho, n, k, alpha=0.05):
    """Fisher-z confidence interval. k is the number of covariates."""
    if not np.isfinite(rho) or n - 3 - k <= 0:
        return np.nan, np.nan
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3 - k)
    zcrit = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se))


def select_max(effects, neg_counts, worst_p, mean_abs):
    """Index of the eligible specification with the largest minimum |effect|.

    effects: (m, 3) pooled effects for TACSTD2, CLDN4, CLDN7.
    A row is eligible when every pooled effect is finite and negative and
    every predictor is negative in at least 7 cohorts.
    Tie-break is larger mean |effect|, then smaller worst p.
    Returns -1 when nothing is eligible.
    """
    effects = np.asarray(effects, dtype=float)
    neg_counts = np.asarray(neg_counts, dtype=float)
    worst_p = np.asarray(worst_p, dtype=float)
    mean_abs = np.asarray(mean_abs, dtype=float)
    min_abs = np.min(np.abs(effects), axis=1)
    eligible = (
        np.all(np.isfinite(effects), axis=1)
        & np.all(effects < 0, axis=1)
        & np.all(neg_counts >= 7, axis=1)
        & np.isfinite(min_abs)
    )
    if not np.any(eligible):
        return -1
    score_p = np.where(np.isfinite(worst_p), worst_p, np.inf)
    # lexsort uses the last key as primary. Want max min_abs, then max mean_abs, then min p.
    order = np.lexsort((score_p, -mean_abs, -min_abs))
    for idx in order:
        if eligible[idx]:
            return int(idx)
    return -1
