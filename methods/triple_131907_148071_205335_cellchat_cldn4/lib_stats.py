#!/usr/bin/env python3
"""Patient-level Spearman helpers and DerSimonian–Laird / Stouffer pooling."""
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


def spearman_p_from_rho(rho: float, n: int) -> float:
    """Two-sided Spearman p from a locked ρ and n (t approximation)."""
    if n <= 2 or not math.isfinite(rho) or abs(rho) >= 1:
        return float("nan")
    t = rho * math.sqrt((n - 2) / (1.0 - rho * rho))
    return float(2 * stats.t.sf(abs(t), n - 2))


def fisher_z(rho: float) -> float:
    r = float(np.clip(rho, -0.999999, 0.999999))
    return float(np.arctanh(r))


def fisher_z_var(n: int) -> float:
    if n <= 3:
        return float("nan")
    return 1.0 / (n - 3)


def random_effects_dl(rhos: list[float], ns: list[int]) -> dict:
    """DerSimonian–Laird RE on Fisher-z of Spearman ρ (or rank-biserial r)."""
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
    return {
        "k": k,
        "n_patients_total": int(n_ok.sum()),
        "method": "DerSimonian-Laird random-effects on Fisher-z",
        "pooled_z": z_re,
        "pooled_rho": float(np.tanh(z_re)),
        "se_z": se,
        "ci95_rho": [float(np.tanh(lo)), float(np.tanh(hi))],
        "p": p,
        "Q": q,
        "df": df,
        "tau2": tau2,
        "I2": i2,
        "fixed_z": z_fe,
        "fixed_rho": float(np.tanh(z_fe)),
    }


def stouffer(rhos: list[float], ps: list[float], ns: list[int]) -> dict:
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
    return {"k": int(z.size), "z": z_comb, "p": p}


def q4_vs_q1(cldn4, immune) -> dict | None:
    """Within-cohort CLDN4 quartiles; MWU on immune. r_rb < 0 = Q4 immune lower."""
    import pandas as pd

    frame = pd.DataFrame(
        {"c": np.asarray(cldn4, dtype=float), "i": np.asarray(immune, dtype=float)}
    )
    frame = frame[np.isfinite(frame["c"]) & np.isfinite(frame["i"])].copy()
    n = int(len(frame))
    if n < 6:
        return None
    ranks = frame["c"].rank(method="average")
    try:
        qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    except ValueError:
        return None
    if qs.nunique() < 4:
        return None
    q1 = frame.loc[qs == "Q1", "i"]
    q4 = frame.loc[qs == "Q4", "i"]
    n1, n4 = int(len(q1)), int(len(q4))
    if n1 < 2 or n4 < 2:
        return None
    u, p = stats.mannwhitneyu(q4.values, q1.values, alternative="two-sided")
    r_rb = (2.0 * float(u)) / (n4 * n1) - 1.0
    return {
        "n": n,
        "n_q1": n1,
        "n_q4": n4,
        "n_compared": n1 + n4,
        "median_q1": float(q1.median()),
        "median_q4": float(q4.median()),
        "delta_median": float(q4.median() - q1.median()),
        "mwu_u": float(u),
        "p": float(p),
        "r_rb": float(r_rb),
        "thin": n < 8 or n1 < 3 or n4 < 3,
        "quartiles": qs.astype(str).to_numpy(),
    }
