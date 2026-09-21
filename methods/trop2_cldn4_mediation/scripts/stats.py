"""Association and product-method mediation helpers.

Partial Spearman is the classical formula on the three Spearman correlations.
The rank-residual Pearson check regresses ranks on the control rank and
correlates the residuals. With ties the two can differ; the formula is the
reported partial, and a sign disagreement is recorded.

ACME is the Imai product estimator on variables standardized within the
analysis unit (mean 0, sample SD 1): a from M ~ X, b from Y ~ X + M,
ACME = a*b. The bootstrap resamples rows. It is an observational product
method, not a proof of a causal mechanism.
"""

from __future__ import annotations

import numpy as np
from scipy import stats


def spearman(x, y) -> tuple[float, float]:
    r, p = stats.spearmanr(x, y)
    return float(r), float(p)


def _partial_from_corrs(rxy: float, rxz: float, ryz: float) -> float:
    den = np.sqrt(max(0.0, (1.0 - rxz * rxz) * (1.0 - ryz * ryz)))
    if den <= 0:
        return float("nan")
    return float((rxy - rxz * ryz) / den)


def partial_spearman(x, y, z) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    n = int(len(x))
    rxy, _ = spearman(x, y)
    rxz, _ = spearman(x, z)
    ryz, _ = spearman(y, z)
    r = _partial_from_corrs(rxy, rxz, ryz)
    if not np.isfinite(r) or abs(r) >= 1.0 or n <= 4:
        p = float("nan")
    else:
        tstat = r * np.sqrt((n - 3) / (1.0 - r * r))
        p = float(2 * stats.t.sf(abs(tstat), n - 3))
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    design = np.column_stack([np.ones(n), rz])
    bx, *_ = np.linalg.lstsq(design, rx, rcond=None)
    by, *_ = np.linalg.lstsq(design, ry, rcond=None)
    rr, rp = stats.pearsonr(rx - design @ bx, ry - design @ by)
    return {
        "n": n,
        "rho": float(r),
        "p": p,
        "rho_xy": float(rxy),
        "rho_xz": float(rxz),
        "rho_yz": float(ryz),
        "rank_residual_rho": float(rr),
        "rank_residual_p": float(rp),
    }


def attenuation_percent(raw: float, partial: float, floor: float = 0.05) -> float:
    """Percent of the raw Spearman removed by the partial.

    Positive means the partial is closer to zero than the raw coefficient
    when both are on the same side, or the partial has crossed zero.
    Undefined when |raw| is below the floor.
    """
    if not np.isfinite(raw) or not np.isfinite(partial) or abs(raw) < floor:
        return float("nan")
    return float(100.0 * (raw - partial) / raw)


def _standardize(a: np.ndarray) -> np.ndarray | None:
    a = np.asarray(a, dtype=float)
    sd = float(a.std(ddof=1))
    if not np.isfinite(sd) or sd == 0.0:
        return None
    return (a - float(a.mean())) / sd


def acme_point(x, m, y) -> dict | None:
    xz = _standardize(np.asarray(x, dtype=float))
    mz = _standardize(np.asarray(m, dtype=float))
    yz = _standardize(np.asarray(y, dtype=float))
    if xz is None or mz is None or yz is None:
        return None
    a = float(np.polyfit(xz, mz, 1)[0])
    design = np.column_stack([np.ones(len(xz)), xz, mz])
    coef, *_ = np.linalg.lstsq(design, yz, rcond=None)
    direct = float(coef[1])
    b = float(coef[2])
    indirect = a * b
    return {
        "a": a,
        "b": b,
        "acme": float(indirect),
        "ade": direct,
        "total": float(direct + indirect),
    }


def acme_bootstrap(x, m, y, n_boot: int = 4000, seed: int = 4639) -> dict:
    x = np.asarray(x, dtype=float)
    m = np.asarray(m, dtype=float)
    y = np.asarray(y, dtype=float)
    point = acme_point(x, m, y)
    if point is None:
        return {"acme": np.nan, "ade": np.nan, "total": np.nan, "a": np.nan, "b": np.nan,
                "acme_lo": np.nan, "acme_hi": np.nan, "acme_se": np.nan, "acme_p": np.nan}
    n = len(x)
    rng = np.random.default_rng(seed)
    acc = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        ix = rng.integers(0, n, n)
        est = acme_point(x[ix], m[ix], y[ix])
        acc[i] = np.nan if est is None else est["acme"]
    ok = acc[np.isfinite(acc)]
    lo, hi = np.quantile(ok, [0.025, 0.975])
    # two-sided bootstrap p, percentile method, bounded at 1
    p = 2.0 * min(float(np.mean(ok <= 0.0)), float(np.mean(ok >= 0.0)))
    point.update({
        "acme_lo": float(lo),
        "acme_hi": float(hi),
        "acme_se": float(ok.std(ddof=1)),
        "acme_p": float(min(p, 1.0)),
        "n_boot": int(len(ok)),
    })
    return point


def fisher_dl(rhos, ns, partial: bool = False) -> dict:
    """DerSimonian–Laird pool of Spearman (or partial Spearman) coefficients.

    Fisher z variance is 1/(n-3) for a raw Spearman and 1/(n-4) for a
    partial Spearman with one control.
    """
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    df = ns - (4.0 if partial else 3.0)
    if np.any(df <= 0):
        raise ValueError(f"non-positive Fisher df: {df.tolist()}")
    var = 1.0 / df
    w = 1.0 / var
    zbar = float(np.sum(w * z) / np.sum(w))
    q = float(np.sum(w * (z - zbar) ** 2))
    k = int(len(z))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    wr = 1.0 / (var + tau2)
    zp = float(np.sum(wr * z) / np.sum(wr))
    se = float(np.sqrt(1.0 / np.sum(wr)))
    rho = float(np.tanh(zp))
    p = float(2 * stats.norm.sf(abs(zp / se))) if se > 0 else float("nan")
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
    return {
        "rho": rho,
        "p": p,
        "I2": float(i2),
        "tau2": float(tau2),
        "ci_lo": float(np.tanh(zp - 1.96 * se)),
        "ci_hi": float(np.tanh(zp + 1.96 * se)),
        "k": k,
        "N": int(ns.sum()),
    }


def dl_means(estimates, ses) -> dict:
    estimates = np.asarray(estimates, dtype=float)
    ses = np.asarray(ses, dtype=float)
    var = ses ** 2
    if np.any(~np.isfinite(var) | (var <= 0)):
        raise ValueError("ACME SEs must be positive")
    w = 1.0 / var
    mu = float(np.sum(w * estimates) / np.sum(w))
    q = float(np.sum(w * (estimates - mu) ** 2))
    k = int(len(estimates))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    wr = 1.0 / (var + tau2)
    m = float(np.sum(wr * estimates) / np.sum(wr))
    se = float(np.sqrt(1.0 / np.sum(wr)))
    p = float(2 * stats.norm.sf(abs(m / se))) if se > 0 else float("nan")
    i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
    return {
        "estimate": m,
        "se": se,
        "p": p,
        "I2": float(i2),
        "ci_lo": float(m - 1.96 * se),
        "ci_hi": float(m + 1.96 * se),
        "k": k,
    }


def spearman_ci(rho: float, n: int) -> tuple[float, float]:
    if n <= 3 or not np.isfinite(rho):
        return float("nan"), float("nan")
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - 3)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def verdict_row(rho_x, rho_x_lo, rho_x_hi, rho_partial, rho_partial_p,
                acme, acme_lo, acme_hi, rho_m_partial) -> str:
    """Call for 'TROP2-immune depends partly on CLDN4'.

    Immune-cold is a negative association. Intervals are 95%.

    opposite: raw TACSTD2–immune interval is entirely positive, or the
    TACSTD2|CLDN4 partial is positive with p < 0.05 while the raw
    association is not significantly immune-cold.
    support: raw interval is entirely negative, Spearman attenuation is
    between 20% and 100%, ACME is negative and its interval excludes 0,
    and CLDN4|TACSTD2 stays negative.
    partial: raw ρ ≤ −0.20, the partial is less negative than the raw
    coefficient, and ACME is negative, but the support intervals are not
    all met.
    null: the remaining rows. Attenuation is not interpreted when |raw| < 0.05.
    """
    raw_hot_sig = rho_x_lo > 0 and rho_x > 0
    raw_cold_sig = rho_x_hi < 0 and rho_x < 0
    residual_hot_sig = rho_partial > 0 and rho_partial_p < 0.05 and not raw_cold_sig
    att = attenuation_percent(rho_x, rho_partial)
    shrunk = rho_x < 0 and rho_partial > rho_x
    acme_neg_sig = acme_hi < 0 and acme < 0
    acme_neg = acme < 0
    m_still_neg = rho_m_partial < 0
    att_in_band = np.isfinite(att) and 20.0 <= att <= 100.0
    if raw_hot_sig or residual_hot_sig:
        return "opposite"
    if raw_cold_sig and shrunk and att_in_band and acme_neg_sig and m_still_neg:
        return "support"
    if rho_x <= -0.20 and shrunk and acme_neg:
        return "partial"
    return "null"
