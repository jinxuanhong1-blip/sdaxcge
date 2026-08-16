"""Shared signatures and stats for leftover lung spatial analyses."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

EPI = ["CLDN4", "TACSTD2"]
BROAD_EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
T_GENES = ["CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "TRAC", "CD247", "IL7R"]
B_GENES = ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD22"]
MYELOID = ["CD68", "CD14", "CSF1R", "LYZ", "FCGR3A"]
IMMUNE = T_GENES + B_GENES + ["PTPRC"]


def present(names, universe):
    u = {str(x).upper() for x in universe}
    return [g for g in names if g.upper() in u]


def spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 8:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(rho), "p": float(p)}


def partial_spearman(x, y, z):
    """Rank-based partial Spearman of x vs y controlling for z."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 12:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rx = stats.rankdata(x[m])
    ry = stats.rankdata(y[m])
    rz = stats.rankdata(z[m])
    rz = (rz - rz.mean()) / (rz.std() + 1e-12)
    bx = np.polyfit(rz, rx, 1)
    by = np.polyfit(rz, ry, 1)
    ex = rx - (bx[0] * rz + bx[1])
    ey = ry - (by[0] * rz + by[1])
    rho, p = stats.spearmanr(ex, ey)
    return {"n": n, "rho": float(rho), "p": float(p)}


def wilcoxon_signed(values):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    n = int(v.size)
    if n < 3:
        return {"n": n, "median": float(np.median(v)) if n else np.nan, "p": np.nan, "n_neg": int((v < 0).sum())}
    try:
        p = float(stats.wilcoxon(v, alternative="two-sided").pvalue)
    except ValueError:
        p = np.nan
    return {
        "n": n,
        "median": float(np.median(v)),
        "p": p,
        "n_neg": int((v < 0).sum()),
        "n_pos": int((v > 0).sum()),
    }


def score_mean(mat, genes, axis=0):
    """Mean of available genes. mat is samples x genes or genes x samples."""
    cols = [g for g in genes if g in mat.columns]
    if not cols:
        return pd.Series(np.nan, index=mat.index)
    return mat[cols].mean(axis=1)


def log_cp10k(counts):
    lib = np.asarray(counts.sum(axis=1)).ravel()
    lib = np.clip(lib, 1, None)
    x = counts.multiply(1e4 / lib[:, None]) if hasattr(counts, "multiply") else counts.div(lib, axis=0) * 1e4
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.log1p(np.asarray(x))


def hex_neighbors(row, col):
    return [
        (row, col - 2),
        (row, col + 2),
        (row - 1, col - 1),
        (row - 1, col + 1),
        (row + 1, col - 1),
        (row + 1, col + 1),
    ]


def ring_neighbors(row, col, k):
    """Hex graph distance k (self excluded)."""
    if k <= 0:
        return []
    seen = {(row, col)}
    frontier = {(row, col)}
    for _ in range(k):
        nxt = set()
        for r, c in frontier:
            for nb in hex_neighbors(r, c):
                if nb not in seen:
                    nxt.add(nb)
                    seen.add(nb)
        frontier = nxt
    return list(frontier)


def dump_json(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, default=_json_default) + "\n")


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if pd.isna(o):
        return None
    raise TypeError(type(o))
