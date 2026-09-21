"""Partial Spearman and rank-OLS mediation.

Partial Spearman is the Pearson correlation of average-rank residuals after
ordinary least squares on an intercept plus the rank-transformed covariates.
The t reference uses df = n - 2 - k. Fisher-z intervals use variance
1/(n - 3 - k). That is the keratin-adjusted correlation used in the earlier
TCGA Trop2 scripts.

Rank-OLS mediation uses those same ranks. With an intercept and the same
covariates in the total, mediator, and direct models, the OLS algebra gives

    c - c' = a * b

exactly. The proportion mediated is a*b / c. It is an accounting identity
on the observational ranks, not a causal effect.
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


def covariate_matrix(covariates, n: int) -> np.ndarray:
    """Stack covariates into an (n, k) array. None or () yields k = 0."""
    if covariates is None:
        return np.zeros((n, 0), dtype=float)
    if isinstance(covariates, np.ndarray):
        cov = np.asarray(covariates, dtype=float)
        if cov.ndim == 1:
            if cov.size != n:
                raise ValueError(f"covariate length {cov.size} != n {n}")
            return cov.reshape(n, 1)
        if cov.ndim != 2 or cov.shape[0] != n:
            raise ValueError(f"covariate shape {cov.shape} incompatible with n {n}")
        return cov
    cols = []
    for z in covariates:
        col = np.asarray(z, dtype=float).reshape(-1)
        if col.size != n:
            raise ValueError(f"covariate length {col.size} != n {n}")
        cols.append(col)
    if not cols:
        return np.zeros((n, 0), dtype=float)
    return np.column_stack(cols)


def _complete_mask(arrays) -> np.ndarray:
    mask = np.ones(arrays[0].shape[0], dtype=bool)
    for arr in arrays:
        if arr.ndim == 1:
            mask &= np.isfinite(arr)
        else:
            mask &= np.isfinite(arr).all(axis=1)
    return mask


def _partial_pearson(x: np.ndarray, y: np.ndarray, Z: np.ndarray) -> float:
    """Pearson correlation of residuals. x, y, Z are already transformed."""
    n = int(x.size)
    if Z.size:
        design = np.column_stack([np.ones(n), Z])
        bx, *_ = np.linalg.lstsq(design, x, rcond=None)
        by, *_ = np.linalg.lstsq(design, y, rcond=None)
        rx = x - design @ bx
        ry = y - design @ by
    else:
        rx = x - np.mean(x)
        ry = y - np.mean(y)
    rx = rx - np.mean(rx)
    ry = ry - np.mean(ry)
    denom = np.sqrt(np.dot(rx, rx) * np.dot(ry, ry))
    if denom < 1e-12:
        return np.nan
    return float(np.dot(rx, ry) / denom)


def partial_spearman(x, y, covariates=()):
    """Spearman rho of x vs y controlling for covariates.

    Returns rho, two-sided p, and the number of complete observations.
    """
    x = np.asarray(x, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    n_in = int(x.size)
    if y.size != n_in:
        raise ValueError("x and y length mismatch")
    Z = covariate_matrix(covariates, n_in)
    mask = _complete_mask([x, y, Z])
    x = x[mask]
    y = y[mask]
    Z = Z[mask]
    n = int(x.size)
    k = int(Z.shape[1])
    if n < k + 5:
        return np.nan, np.nan, n
    xr = stats.rankdata(x).astype(float)
    yr = stats.rankdata(y).astype(float)
    if k:
        Z = np.column_stack([stats.rankdata(Z[:, j]).astype(float) for j in range(k)])
    r = _partial_pearson(xr, yr, Z)
    df = n - 2 - k
    if not np.isfinite(r) or df <= 0:
        return r, np.nan, n
    if abs(r) >= 1.0 - 1e-15:
        return float(r), 0.0, n
    tstat = r * np.sqrt(df / (1.0 - r * r))
    p = float(2 * stats.t.sf(abs(tstat), df))
    return float(r), p, n


def fisher_ci(rho, n, k, alpha=0.05):
    """Fisher-z confidence interval. k is the number of covariates."""
    if not np.isfinite(rho) or n - 3 - k <= 0:
        return np.nan, np.nan
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3 - k)
    zcrit = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se))


def random_effects_meta(rhos, ns, k):
    """DerSimonian-Laird random-effects meta-analysis on Fisher z.

    Variance of z is 1/(n - 3 - k). Cohorts with non-finite rho or
    n - 3 - k <= 1 are dropped. I2 is returned as a fraction in [0, 1].
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


def _mediation_ranked(x, m, y, Z):
    """Product-of-coefficients mediation on values that are already ranked.

    Returns the unstandardized OLS coefficients. c - c_prime equals a * b
    up to numerical error.
    """
    n = int(x.size)
    k = int(Z.shape[1]) if Z.size else 0
    ones = np.ones(n)
    if k:
        design = np.column_stack([ones, x, Z])
        direct = np.column_stack([ones, x, m, Z])
    else:
        design = np.column_stack([ones, x])
        direct = np.column_stack([ones, x, m])
    c = float(np.linalg.lstsq(design, y, rcond=None)[0][1])
    a = float(np.linalg.lstsq(design, m, rcond=None)[0][1])
    beta = np.linalg.lstsq(direct, y, rcond=None)[0]
    c_prime = float(beta[1])
    b = float(beta[2])
    indirect = float(a * b)
    proportion = float(indirect / c) if abs(c) > 1e-12 else np.nan
    return {
        "a": a,
        "b": b,
        "c": c,
        "c_prime": c_prime,
        "indirect": indirect,
        "proportion": proportion,
        "n": n,
        "k": k,
        "identity_gap": float((c - c_prime) - indirect),
    }


def product_mediation(x, mediator, y, covariates=(), rank: bool = True):
    """Proportion of the X→Y coefficient accounted for by a mediator.

    rank=True transforms every variable to average ranks before OLS, matching
    the partial-Spearman scale. rank=False runs OLS on the supplied values.
    """
    x = np.asarray(x, dtype=float).reshape(-1)
    m = np.asarray(mediator, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    n_in = int(x.size)
    if m.size != n_in or y.size != n_in:
        raise ValueError("x, mediator, and y length mismatch")
    Z = covariate_matrix(covariates, n_in)
    mask = _complete_mask([x, m, y, Z])
    x, m, y = x[mask], m[mask], y[mask]
    Z = Z[mask]
    n = int(x.size)
    k = int(Z.shape[1])
    empty = {
        "a": np.nan,
        "b": np.nan,
        "c": np.nan,
        "c_prime": np.nan,
        "indirect": np.nan,
        "proportion": np.nan,
        "n": n,
        "k": k,
        "identity_gap": np.nan,
    }
    if n < k + 5:
        return empty
    if rank:
        x = stats.rankdata(x).astype(float)
        m = stats.rankdata(m).astype(float)
        y = stats.rankdata(y).astype(float)
        if k:
            Z = np.column_stack([stats.rankdata(Z[:, j]).astype(float) for j in range(k)])
    elif k == 0:
        Z = np.zeros((n, 0), dtype=float)
    return _mediation_ranked(x, m, y, Z)


def attenuation_fraction(rho_before, rho_after):
    """Share of rho_before removed by the later model: (before - after) / before.

    Positive when the later correlation moves toward zero from a nonzero
    starting value. Greater than 1 when the later correlation crosses zero.
    Negative when the later correlation moves farther from zero.
    """
    if not np.isfinite(rho_before) or not np.isfinite(rho_after) or abs(rho_before) < 1e-12:
        return np.nan
    return float((rho_before - rho_after) / rho_before)


def coef_given(y, x, Z) -> float:
    """OLS coefficient of x in y ~ 1 + x + Z. Arrays are used as supplied."""
    n = int(x.size)
    if Z.size:
        design = np.column_stack([np.ones(n), x, Z])
    else:
        design = np.column_stack([np.ones(n), x])
    return float(np.linalg.lstsq(design, y, rcond=None)[0][1])
