"""Patient-level Spearman, Cliff δ, and DerSimonian–Laird pooling."""
from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from scipy import stats


def spearman(x: Iterable[float], y: Iterable[float]) -> tuple[float, float, int]:
    xa = np.asarray(list(x), dtype=float)
    ya = np.asarray(list(y), dtype=float)
    mask = np.isfinite(xa) & np.isfinite(ya)
    xa, ya = xa[mask], ya[mask]
    n = int(xa.size)
    if n < 4 or np.nanstd(xa) == 0 or np.nanstd(ya) == 0:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(xa, ya)
    return float(rho), float(p), n


def cliffs_delta(high: Iterable[float], low: Iterable[float]) -> float:
    """Cliff δ = P(high>low) − P(high<low). Negative ⇒ CLDN4-high is immune-lower."""
    a = np.asarray(list(high), dtype=float)
    b = np.asarray(list(low), dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n1, n2 = int(a.size), int(b.size)
    if n1 < 1 or n2 < 1:
        return float("nan")
    # vectorized pairwise
    diff = a[:, None] - b[None, :]
    gt = float(np.sum(diff > 0))
    lt = float(np.sum(diff < 0))
    return (gt - lt) / (n1 * n2)


def mwu(high: Iterable[float], low: Iterable[float]) -> tuple[float, float, int, int]:
    a = np.asarray(list(high), dtype=float)
    b = np.asarray(list(low), dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n1, n2 = int(a.size), int(b.size)
    if n1 < 2 or n2 < 2:
        return float("nan"), float("nan"), n1, n2
    res = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(res.statistic), float(res.pvalue), n1, n2


def fisher_z(rho: float) -> float:
    r = float(np.clip(rho, -0.999999, 0.999999))
    return float(np.arctanh(r))


def random_effects_dl(rhos: list[float], ns: list[int]) -> dict:
    """DerSimonian–Laird RE on Fisher-z of Spearman ρ; back-transform to ρ."""
    z = np.array([fisher_z(r) for r in rhos], dtype=float)
    v = np.array([1.0 / (n - 3) if n > 3 else float("nan") for n in ns], dtype=float)
    ok = np.isfinite(z) & np.isfinite(v) & (v > 0)
    z, v, n_ok = z[ok], v[ok], np.array(ns, dtype=float)[ok]
    k = int(z.size)
    if k == 0:
        return {"k": 0}
    w = 1.0 / v
    z_fe = float(np.sum(w * z) / np.sum(w))
    q = float(np.sum(w * (z - z_fe) ** 2))
    df = k - 1
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if k > 1 else float("nan")
    tau2 = max(0.0, (q - df) / c) if k > 1 and c > 0 else 0.0
    w_re = 1.0 / (v + tau2)
    z_re = float(np.sum(w_re * z) / np.sum(w_re))
    se = float(1.0 / math.sqrt(float(np.sum(w_re))))
    z_stat = z_re / se if se > 0 else float("nan")
    p = float(2 * stats.norm.sf(abs(z_stat))) if math.isfinite(z_stat) else float("nan")
    lo, hi = z_re - 1.96 * se, z_re + 1.96 * se
    i2 = max(0.0, (q - df) / q) * 100.0 if k > 1 and q > 0 else 0.0
    return {
        "k": k,
        "n_patients_total": int(n_ok.sum()),
        "method": "DerSimonian-Laird random-effects on Fisher-z(Spearman ρ)",
        "pooled_rho": float(np.tanh(z_re)),
        "se_z": se,
        "ci95_rho": [float(np.tanh(lo)), float(np.tanh(hi))],
        "p": p,
        "I2": i2,
        "Q": q,
        "tau2": tau2,
    }
