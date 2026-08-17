#!/usr/bin/env python3
"""Patient-level Spearman helpers and random-effects / Stouffer / Fisher pooling.

Copied from methods/scrna_meta_tnk/lib_stats.py (origin/cursor/scrna-meta-tnk-4809).
"""
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
    if n < 3:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(xa, ya)
    return float(rho), float(p), n


def fisher_z(rho: float) -> float:
    r = float(np.clip(rho, -0.999999, 0.999999))
    return float(np.arctanh(r))


def fisher_z_var(n: int) -> float:
    if n <= 3:
        return float("nan")
    return 1.0 / (n - 3)


def implied_spearman_p(rho: float, n: int) -> float:
    """Two-sided p from the given Spearman ρ and n (t approximation).

    Used only to attach a p to a *given* ρ for Stouffer / forest labels.
    Does not re-estimate ρ from the matrix.
    """
    if n <= 2 or not math.isfinite(rho) or abs(rho) >= 1:
        return float("nan")
    t = rho * math.sqrt((n - 2) / (1.0 - rho * rho))
    return float(2 * stats.t.sf(abs(t), n - 2))


def random_effects_dl(rhos: list[float], ns: list[int]) -> dict:
    """DerSimonian–Laird RE on Fisher-z of Spearman ρ; back-transform to ρ."""
    z = np.array([fisher_z(r) for r in rhos], dtype=float)
    v = np.array([fisher_z_var(n) for n in ns], dtype=float)
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
    h2 = q / df if k > 1 and df > 0 else float("nan")
    return {
        "k": k,
        "n_patients_total": int(n_ok.sum()),
        "method": "DerSimonian-Laird random-effects on Fisher-z(Spearman ρ)",
        "pooled_z": z_re,
        "pooled_rho": float(np.tanh(z_re)),
        "se_z": se,
        "ci95_rho": [float(np.tanh(lo)), float(np.tanh(hi))],
        "p": p,
        "Q": q,
        "df": df,
        "tau2": tau2,
        "I2": i2,
        "H2": h2,
        "fixed_z": z_fe,
        "fixed_rho": float(np.tanh(z_fe)),
    }


def stouffer(rhos: list[float], ps: list[float], ns: list[int]) -> dict:
    """Signed Stouffer combination. Weight = sqrt(n-3) (Fisher-z information)."""
    zs, ws = [], []
    for rho, p, n in zip(rhos, ps, ns):
        if n <= 3 or not math.isfinite(p) or p <= 0:
            continue
        z_abs = abs(float(stats.norm.ppf(min(max(p / 2.0, 1e-300), 1 - 1e-16))))
        sign = 1.0 if rho >= 0 else -1.0
        zs.append(sign * z_abs)
        ws.append(math.sqrt(n - 3))
    if not zs:
        return {"k": 0}
    z = np.array(zs)
    w = np.array(ws)
    z_comb = float(np.sum(w * z) / math.sqrt(float(np.sum(w**2))))
    p = float(2 * stats.norm.sf(abs(z_comb)))
    return {
        "k": int(z.size),
        "method": "Stouffer weighted by sqrt(n-3), signed by ρ",
        "z": z_comb,
        "p": p,
    }


def fisher_combine(ps: list[float]) -> dict:
    p_ok = [p for p in ps if math.isfinite(p) and 0 < p <= 1]
    k = len(p_ok)
    if k == 0:
        return {"k": 0}
    x2 = float(-2.0 * sum(math.log(p) for p in p_ok))
    p = float(stats.chi2.sf(x2, 2 * k))
    return {
        "k": k,
        "method": "Fisher combined p (direction-agnostic)",
        "X2": x2,
        "df": 2 * k,
        "p": p,
    }


def bh_fdr(ps: list[float]) -> list[float]:
    """Benjamini–Hochberg q-values (independent). NaN stays NaN."""
    arr = np.asarray(ps, dtype=float)
    q = np.full(arr.shape, np.nan)
    ok = np.isfinite(arr) & (arr >= 0) & (arr <= 1)
    if not ok.any():
        return q.tolist()
    p = arr[ok]
    order = np.argsort(p)
    ranked = p[order]
    m = ranked.size
    raw = ranked * m / np.arange(1, m + 1)
    # monotone from the tail
    adj = np.minimum.accumulate(raw[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(m)
    out[order] = adj
    q[ok] = out
    return q.tolist()
