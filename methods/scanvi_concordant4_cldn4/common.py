"""Shared definitions for the concordant-4 scVI/scANVI analysis."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import stats

QC_MIN_GENES = 200
QC_MIN_UMI = 500
QC_MAX_MT = 20.0
MAL_CAP = 220
TNK_CAP = 140
OTHER_CAP = 40

GSE131907_UNITS = [
    "BRONCHO_11",
    "BRONCHO_58",
    "EBUS_06",
    "EBUS_10",
    "EBUS_12",
    "EBUS_13",
    "EBUS_15",
    "EBUS_19",
    "EBUS_28",
    "EBUS_49",
    "EBUS_51",
    "NS_02",
    "NS_03",
    "NS_04",
    "NS_06",
    "NS_07",
    "NS_12",
    "NS_13",
    "NS_16",
    "NS_17",
    "NS_19",
]

EPI_GENES = ["EPCAM", "KRT8", "KRT18", "KRT19"]
TNK_GENES = ["CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1"]
AUTHOR_DATASETS = {"GSE131907", "GSE205335"}


def dl_spearman(rhos, ns) -> dict:
    """DerSimonian–Laird random-effects pool of Spearman rhos on the Fisher-z scale."""
    rhos = np.asarray(rhos, dtype=float)
    ns = np.asarray(ns, dtype=float)
    ok = np.isfinite(rhos) & np.isfinite(ns) & (ns > 3)
    rhos, ns = rhos[ok], ns[ok]
    z = np.arctanh(np.clip(rhos, -0.999999, 0.999999))
    var_z = 1.0 / (ns - 3.0)
    w = 1.0 / var_z
    zbar = np.sum(w * z) / np.sum(w)
    q = float(np.sum(w * (z - zbar) ** 2))
    k = int(len(rhos))
    dfree = k - 1
    cdenom = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - dfree) / cdenom) if dfree > 0 and cdenom > 0 else 0.0
    wstar = 1.0 / (var_z + tau2)
    zre = float(np.sum(wstar * z) / np.sum(wstar))
    se = math.sqrt(1.0 / float(np.sum(wstar)))
    p = float(2 * stats.norm.sf(abs(zre / se)))
    i2 = max(0.0, (q - dfree) / q) if q > 0 else 0.0
    return {
        "rho": float(math.tanh(zre)),
        "p": p,
        "I2": float(i2),
        "ci_lo": float(math.tanh(zre - 1.96 * se)),
        "ci_hi": float(math.tanh(zre + 1.96 * se)),
        "k": k,
        "N": int(ns.sum()),
        "tau2": float(tau2),
    }


def spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4:
        return float("nan"), float("nan"), n
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), n


def within_quartile(x: pd.Series) -> pd.Series:
    r = x.rank(method="average")
    breaks = np.unique(np.nanquantile(r.to_numpy(dtype=float), [0, 0.25, 0.5, 0.75, 1]))
    if len(breaks) < 3:
        return pd.Series(np.nan, index=x.index)
    labels = ["Q1", "Q2", "Q3", "Q4"][: len(breaks) - 1]
    return pd.cut(r, bins=breaks, include_lowest=True, labels=labels)


def rank_biserial(q4, q1) -> dict:
    q4 = np.asarray(q4, dtype=float)
    q1 = np.asarray(q1, dtype=float)
    q4, q1 = q4[np.isfinite(q4)], q1[np.isfinite(q1)]
    n4, n1 = int(len(q4)), int(len(q1))
    if n4 < 1 or n1 < 1:
        return {"r": float("nan"), "p": float("nan"), "n_q1": n1, "n_q4": n4}
    res = stats.mannwhitneyu(q4, q1, alternative="two-sided", method="asymptotic")
    u = float(res.statistic)
    return {"r": 2 * u / (n4 * n1) - 1, "p": float(res.pvalue), "n_q1": n1, "n_q4": n4}


def assert_locked_dl() -> None:
    """The published four-cohort Spearmans must still pool to ρ≈−0.531."""
    rhos = [-0.659340659340659, -0.522077922077922, -0.435347261434218, -0.6]
    ns = [13, 21, 22, 9]
    out = dl_spearman(rhos, ns)
    if abs(out["rho"] - (-0.5311678045689988)) > 1e-9:
        raise AssertionError(f"DL self-check failed: {out}")
    if abs(out["p"] - 1.6462232944602917e-05) > 1e-12:
        raise AssertionError(f"DL p self-check failed: {out}")
