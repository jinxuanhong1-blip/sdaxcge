"""Shared statistics for histology-stratified scRNA Spearman meta-analysis."""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from scipy import stats


HISTOLOGY_MAP = {
    "luad": "LUAD",
    "adeno": "LUAD",
    "adenocarcinoma": "LUAD",
    "adc": "LUAD",
    "lung adenocarcinoma": "LUAD",
    "lusc": "LUSC",
    "squamous": "LUSC",
    "sq": "LUSC",
    "sqcc": "LUSC",
    "lung squamous": "LUSC",
    "lung squamous cell carcinoma": "LUSC",
    "asc": "ASC",
    "adeno-squamous": "ASC",
    "adenosquamous": "ASC",
    "sclc": "OTHER",
    "nut": "OTHER",
    "nsclc": "UNKNOWN",
}


def map_histology(raw: str | None) -> str:
    if raw is None or (isinstance(raw, float) and math.isnan(raw)):
        return "UNKNOWN"
    key = str(raw).strip().lower()
    return HISTOLOGY_MAP.get(key, "UNKNOWN")


def spearman_rho_p(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(x.size)
    if n < 4:
        return (float("nan"), float("nan"), n)
    rho, p = stats.spearmanr(x, y)
    return (float(rho), float(p), n)


def fisher_z(rho: float) -> float:
    rho = float(np.clip(rho, -0.999999, 0.999999))
    return float(np.arctanh(rho))


def fisher_z_se(n: int) -> float:
    if n <= 3:
        return float("nan")
    return 1.0 / math.sqrt(n - 3)


def tanh_rho(z: float) -> float:
    return float(np.tanh(z))


def meta_fisher_z(rows: Iterable[dict], min_n: int = 6) -> dict:
    """Fixed- and DerSimonian–Laird random-effects meta of Spearman ρ via Fisher z."""
    usable = []
    for r in rows:
        n = int(r["n"])
        rho = r.get("rho")
        if n < min_n or rho is None or not np.isfinite(rho):
            continue
        z = fisher_z(float(rho))
        se = fisher_z_se(n)
        if not np.isfinite(se):
            continue
        usable.append({**r, "z": z, "se": se, "w": 1.0 / (se * se)})
    k = len(usable)
    if k == 0:
        return {
            "k": 0,
            "n_total": 0,
            "rho_fe": float("nan"),
            "rho_re": float("nan"),
            "ci_fe_low": float("nan"),
            "ci_fe_high": float("nan"),
            "ci_re_low": float("nan"),
            "ci_re_high": float("nan"),
            "p_fe": float("nan"),
            "p_re": float("nan"),
            "tau2": float("nan"),
            "I2": float("nan"),
            "Q": float("nan"),
            "Q_p": float("nan"),
        }
    w = np.array([u["w"] for u in usable], dtype=float)
    z = np.array([u["z"] for u in usable], dtype=float)
    n_total = int(sum(int(u["n"]) for u in usable))
    z_fe = float(np.sum(w * z) / np.sum(w))
    se_fe = float(1.0 / math.sqrt(np.sum(w)))
    p_fe = float(2 * stats.norm.sf(abs(z_fe / se_fe)))
    Q = float(np.sum(w * (z - z_fe) ** 2))
    df = k - 1
    Q_p = float(stats.chi2.sf(Q, df)) if df > 0 else float("nan")
    C = float(np.sum(w) - np.sum(w**2) / np.sum(w)) if df > 0 else float("nan")
    tau2 = max(0.0, (Q - df) / C) if df > 0 and C > 0 else 0.0
    I2 = max(0.0, (Q - df) / Q) * 100.0 if Q > 0 and df > 0 else 0.0
    w_re = 1.0 / (np.array([u["se"] ** 2 for u in usable]) + tau2)
    z_re = float(np.sum(w_re * z) / np.sum(w_re))
    se_re = float(1.0 / math.sqrt(np.sum(w_re)))
    p_re = float(2 * stats.norm.sf(abs(z_re / se_re)))
    zcrit = 1.959963984540054
    return {
        "k": k,
        "n_total": n_total,
        "rho_fe": tanh_rho(z_fe),
        "rho_re": tanh_rho(z_re),
        "ci_fe_low": tanh_rho(z_fe - zcrit * se_fe),
        "ci_fe_high": tanh_rho(z_fe + zcrit * se_fe),
        "ci_re_low": tanh_rho(z_re - zcrit * se_re),
        "ci_re_high": tanh_rho(z_re + zcrit * se_re),
        "p_fe": p_fe,
        "p_re": p_re,
        "tau2": tau2,
        "I2": I2,
        "Q": Q,
        "Q_p": Q_p,
        "z_fe": z_fe,
        "z_re": z_re,
        "se_fe": se_fe,
        "se_re": se_re,
    }


def rho_ci(rho: float, n: int) -> tuple[float, float]:
    if n <= 3 or not np.isfinite(rho):
        return (float("nan"), float("nan"))
    z = fisher_z(rho)
    se = fisher_z_se(n)
    zcrit = 1.959963984540054
    return (tanh_rho(z - zcrit * se), tanh_rho(z + zcrit * se))
