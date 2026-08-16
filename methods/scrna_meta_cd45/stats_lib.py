"""Patient-level two-group tests and a tiny inverse-variance meta.

Effect direction: MPR minus non-MPR. Negative means lower leak in MPR.
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from scipy import stats


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's δ: P(a>b) - P(a<b). Positive => a higher than b."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0 or b.size == 0:
        return float("nan")
    gt = 0
    lt = 0
    for x in a:
        gt += int(np.sum(x > b))
        lt += int(np.sum(x < b))
    return (gt - lt) / (a.size * b.size)


def hedges_g(a: np.ndarray, b: np.ndarray) -> dict:
    """Hedges' g for a minus b, with Hedges-Olkin variance for IV meta."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    n1, n2 = a.size, b.size
    if n1 < 2 or n2 < 2:
        return {
            "hedges_g": float("nan"),
            "se": float("nan"),
            "ci_lo": float("nan"),
            "ci_hi": float("nan"),
            "n1": n1,
            "n2": n2,
        }
    m1, m2 = float(a.mean()), float(b.mean())
    s1, s2 = float(a.std(ddof=1)), float(b.std(ddof=1))
    df = n1 + n2 - 2
    sp2 = ((n1 - 1) * s1 * s1 + (n2 - 1) * s2 * s2) / df
    if sp2 <= 0:
        d = 0.0
    else:
        d = (m1 - m2) / math.sqrt(sp2)
    j = 1.0 - 3.0 / (4.0 * df - 1.0) if df > 0 else 1.0
    g = j * d
    # Hedges & Olkin sampling variance of g
    var = (n1 + n2) / (n1 * n2) + (g * g) / (2.0 * (n1 + n2))
    se = math.sqrt(var)
    z = 1.959963984540054
    return {
        "hedges_g": g,
        "se": se,
        "ci_lo": g - z * se,
        "ci_hi": g + z * se,
        "n1": n1,
        "n2": n2,
        "mean_a": m1,
        "mean_b": m2,
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
    }


def mannwhitney(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0 or b.size == 0:
        return {"U": float("nan"), "p": float("nan"), "cliffs_delta": float("nan")}
    res = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "U": float(res.statistic),
        "p": float(res.pvalue),
        "cliffs_delta": cliffs_delta(a, b),
    }


def two_group(a: Iterable[float], b: Iterable[float], label_a="MPR", label_b="non-MPR") -> dict:
    a = np.asarray(list(a), dtype=float)
    b = np.asarray(list(b), dtype=float)
    out = {
        "label_a": label_a,
        "label_b": label_b,
        "n_a": int(a.size),
        "n_b": int(b.size),
    }
    out.update({f"mw_{k}": v for k, v in mannwhitney(a, b).items()})
    out.update({f"g_{k}": v for k, v in hedges_g(a, b).items()})
    return out


def iv_meta(effects: list[dict]) -> dict:
    """Fixed-effect IV plus DerSimonian-Laird random-effects on hedges_g/se."""
    usable = [e for e in effects if np.isfinite(e.get("hedges_g", np.nan)) and np.isfinite(e.get("se", np.nan)) and e["se"] > 0]
    if not usable:
        return {"k": 0, "fixed_g": float("nan"), "random_g": float("nan")}
    ys = np.array([e["hedges_g"] for e in usable], dtype=float)
    ses = np.array([e["se"] for e in usable], dtype=float)
    w = 1.0 / (ses ** 2)
    yfix = float(np.sum(w * ys) / np.sum(w))
    se_fix = float(math.sqrt(1.0 / np.sum(w)))
    q = float(np.sum(w * (ys - yfix) ** 2))
    k = len(usable)
    df = k - 1
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w)) if k > 1 else float("nan")
    tau2 = max(0.0, (q - df) / c) if k > 1 and c > 0 else 0.0
    wr = 1.0 / (ses ** 2 + tau2)
    yran = float(np.sum(wr * ys) / np.sum(wr))
    se_ran = float(math.sqrt(1.0 / np.sum(wr)))
    z = 1.959963984540054
    i2 = max(0.0, (q - df) / q) if q > 0 and k > 1 else 0.0
    return {
        "k": k,
        "fixed_g": yfix,
        "fixed_se": se_fix,
        "fixed_ci_lo": yfix - z * se_fix,
        "fixed_ci_hi": yfix + z * se_fix,
        "fixed_p": float(2 * stats.norm.sf(abs(yfix / se_fix))) if se_fix > 0 else float("nan"),
        "random_g": yran,
        "random_se": se_ran,
        "random_ci_lo": yran - z * se_ran,
        "random_ci_hi": yran + z * se_ran,
        "random_p": float(2 * stats.norm.sf(abs(yran / se_ran))) if se_ran > 0 else float("nan"),
        "Q": q,
        "I2": i2,
        "tau2": tau2,
        "n_total_mpr": int(sum(e.get("n1", 0) for e in usable)),
        "n_total_nonmpr": int(sum(e.get("n2", 0) for e in usable)),
    }
