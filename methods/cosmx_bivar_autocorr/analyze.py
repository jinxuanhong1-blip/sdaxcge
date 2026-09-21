#!/usr/bin/env python3
"""Bivariate spatial autocorrelation of CLDN4 vs CD8A on CosMx NSCLC.

Dataset: figshare 25976224, He et al. 2022 CosMx SMI NSCLC, CellCharter object
(765,771 cells, 8 sections, 5 patients). Expression is log1p(CP10k) from
layers/counts and obs/n_counts. The object X matrix is not used.

Pre-specified primary estimand
------------------------------
Partial bivariate Moran's I of CLDN4 with the spatial lag of CD8A, after
within-unit OLS residualization of both genes on KRT8 and EPCAM.
Primary neighbor graph: radius 50 µm (0.18 µm/px). Co-primary: kNN k=15.
Both are row-standardized, no self-loops, isolates dropped.

Also reported, same graphs
--------------------------
- Lee's L (Lee 2001), raw and partial
- CLDN4-only residual (CD8A left raw)
- KRT8–CD8A and EPCAM–CD8A Moran (epithelial baselines; same CD8A permutations)
- Reverse Moran: CD8A with the lag of CLDN4
- Within-tumor cross-correlation: tumor-cell CLDN4 (raw and residualized
  inside the tumor set) vs the lag of CD8A. Null shuffles CLDN4 among tumor
  cells and holds the CD8A field fixed.
- Getis-Ord Gi* (binary star weights) for raw expression and for residuals,
  plus CLDN4-hot / CD8A-cold overlap

Inference
---------
199 permutations inside each FOV or section (seeded, order-independent).
FOV inverse-variance meta treats FOVs as independent and is descriptive.
Confirmatory summaries are the 8 sections and the 5 patients.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import json
import math
import time
import zlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.sparse import csr_matrix, diags
from scipy.spatial import Delaunay, cKDTree

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PX_TO_UM = 0.18  # CosMx SMI NSCLC prototype, He et al. 2022
GI_Z = 1.96
MIN_USED = 80
MIN_TUMOR = 30
CACHE_VERSION = 1
GENES = ["CLDN4", "CD8A", "KRT8", "EPCAM"]
GRAPHS = (
    "knn_k6",
    "knn_k15",
    "knn_k30",
    "radius_20um",
    "radius_50um",
    "radius_100um",
    "delaunay",
)
PRIMARY = ("radius_50um", "knn_k15")
PATIENT_ORDER = ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]
PATIENT_COLOR = {
    "Lung5": "#4C78A8",
    "Lung6": "#F58518",
    "Lung9": "#54A24B",
    "Lung12": "#E45756",
    "Lung13": "#B279A2",
}

XY = None
LOGX = None
TUMOR = None
N_PERM = 199
SEED = 20260921


def _safe_spearman(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.std() < 1e-10 or b.std() < 1e-10:
        return float("nan")
    res = stats.spearmanr(a, b)
    rho = getattr(res, "statistic", None)
    if rho is None:
        rho = res[0]
    return float(rho)


def zscore(a: np.ndarray, ddof: int) -> np.ndarray | None:
    a = np.asarray(a, dtype=np.float64)
    sd = float(a.std(ddof=ddof))
    if not np.isfinite(sd) or sd < 1e-10:
        return None
    return (a - a.mean()) / sd


def residualize(y: np.ndarray, Z: np.ndarray) -> tuple[np.ndarray, float]:
    y = np.asarray(y, dtype=np.float64)
    Z = np.asarray(Z, dtype=np.float64)
    if Z.ndim == 1:
        Z = Z[:, None]
    X = np.column_stack([np.ones(len(y)), Z])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    ss = float(np.sum((y - y.mean()) ** 2))
    r2 = float(1.0 - np.sum(resid**2) / ss) if ss > 0 else float("nan")
    return resid, r2


def pack_perm(obs: float, sims: np.ndarray) -> dict:
    sims = np.asarray(sims, dtype=np.float64)
    nperm = int(sims.size)
    if nperm < 2:
        return {
            "est": float(obs),
            "se": float("nan"),
            "z": float("nan"),
            "p_less": float("nan"),
            "p_two": float("nan"),
        }
    se = float(sims.std(ddof=1))
    mu = float(sims.mean())
    z = float((obs - mu) / se) if se > 0 else float("nan")
    p_less = float((1.0 + np.sum(sims <= obs)) / (nperm + 1.0))
    larger = int(np.sum(sims >= obs))
    if (nperm - larger) < larger:
        larger = nperm - larger
    p_two = float((larger + 1.0) / (nperm + 1.0))
    return {"est": float(obs), "se": se, "z": z, "p_less": p_less, "p_two": p_two}


def shuffled_apply(zx: np.ndarray, zy: np.ndarray, W, nperm: int, rng: np.random.Generator):
    """Return observed Moran pieces and WY for columns that are shuffled zy."""
    n = zy.shape[0]
    obs_lag = np.asarray(W @ zy).ravel()
    obs = float(zx @ obs_lag) / (n - 1.0)
    if nperm < 2:
        return obs, np.empty(0), None
    Y = np.empty((n, nperm), dtype=np.float64)
    v = zy.copy()
    for i in range(nperm):
        rng.shuffle(v)
        Y[:, i] = v
    WY = W @ Y
    sims = (zx @ WY) / (n - 1.0)
    return obs, sims, WY


def moran_from_WY(zx: np.ndarray, zy: np.ndarray, W, WY) -> tuple[float, np.ndarray]:
    n = zy.shape[0]
    obs = float(zx @ np.asarray(W @ zy).ravel()) / (n - 1.0)
    sims = (zx @ WY) / (n - 1.0)
    return obs, sims


def lee_apply(zx: np.ndarray, zy: np.ndarray, W, nperm: int, rng: np.random.Generator):
    """Lee L for row-standardized W with no zero rows: (Wzx · Wzy) / n."""
    n = zy.shape[0]
    wzx = np.asarray(W @ zx).ravel()
    obs = float(wzx @ np.asarray(W @ zy).ravel()) / n
    if nperm < 2:
        return obs, np.empty(0), None, wzx
    Y = np.empty((n, nperm), dtype=np.float64)
    v = zy.copy()
    for i in range(nperm):
        rng.shuffle(v)
        Y[:, i] = v
    WY = W @ Y
    sims = (wzx @ WY) / n
    return obs, sims, WY, wzx


def lee_from_WY(zx: np.ndarray, zy: np.ndarray, W, WY, wzx_cd8_side=None) -> tuple[float, np.ndarray]:
    n = zy.shape[0]
    wzx = np.asarray(W @ zx).ravel()
    obs = float(wzx @ np.asarray(W @ zy).ravel()) / n
    sims = (wzx @ WY) / n
    return obs, sims


def knn_binary(xy: np.ndarray, k: int):
    n = xy.shape[0]
    if n <= k + 1:
        return None
    _, ix = cKDTree(xy).query(xy, k=k + 1)
    ix = np.atleast_2d(ix)
    if ix.shape[1] < k + 1 or not np.all(ix[:, 0] == np.arange(n)):
        rows = []
        cols = []
        for i in range(n):
            nbrs = [int(j) for j in np.atleast_1d(ix[i]) if int(j) != i][:k]
            rows.extend([i] * len(nbrs))
            cols.extend(nbrs)
        if not rows:
            return None
        return csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    cols = ix[:, 1:].reshape(-1)
    rows = np.repeat(np.arange(n), k)
    return csr_matrix((np.ones(n * k, dtype=np.float64), (rows, cols)), shape=(n, n))


def radius_binary(xy: np.ndarray, radius_px: float):
    n = xy.shape[0]
    pairs = cKDTree(xy).query_pairs(radius_px, output_type="ndarray")
    if len(pairs) == 0:
        return csr_matrix((n, n), dtype=np.float64)
    rows = np.concatenate([pairs[:, 0], pairs[:, 1]])
    cols = np.concatenate([pairs[:, 1], pairs[:, 0]])
    return csr_matrix((np.ones(len(rows), dtype=np.float64), (rows, cols)), shape=(n, n))


def delaunay_binary(xy: np.ndarray):
    n = xy.shape[0]
    if n < 4:
        return None
    pts = np.asarray(xy, dtype=np.float64)
    try:
        tri = Delaunay(pts)
    except Exception:
        pts = pts + np.random.default_rng(0).normal(scale=1e-6, size=pts.shape)
        try:
            tri = Delaunay(pts)
        except Exception:
            return None
    simplices = tri.simplices
    a = np.concatenate([simplices[:, 0], simplices[:, 1], simplices[:, 0]])
    b = np.concatenate([simplices[:, 1], simplices[:, 2], simplices[:, 2]])
    lo = np.minimum(a, b)
    hi = np.maximum(a, b)
    edges = np.unique(np.stack([lo, hi], axis=1), axis=0)
    dist = np.linalg.norm(pts[edges[:, 0]] - pts[edges[:, 1]], axis=1)
    if dist.size == 0:
        return None
    edges = edges[dist <= np.quantile(dist, 0.99)]
    if edges.size == 0:
        return None
    rows = np.concatenate([edges[:, 0], edges[:, 1]])
    cols = np.concatenate([edges[:, 1], edges[:, 0]])
    return csr_matrix((np.ones(len(rows), dtype=np.float64), (rows, cols)), shape=(n, n))


def build_graphs(xy: np.ndarray) -> dict[str, csr_matrix | None]:
    return {
        "knn_k6": knn_binary(xy, 6),
        "knn_k15": knn_binary(xy, 15),
        "knn_k30": knn_binary(xy, 30),
        "radius_20um": radius_binary(xy, 20.0 / PX_TO_UM),
        "radius_50um": radius_binary(xy, 50.0 / PX_TO_UM),
        "radius_100um": radius_binary(xy, 100.0 / PX_TO_UM),
        "delaunay": delaunay_binary(xy),
    }


def row_standardize_kept(Wbin: csr_matrix):
    deg = np.asarray(Wbin.sum(axis=1)).ravel()
    keep = deg > 0
    if int(keep.sum()) < MIN_USED:
        return None
    idx = np.flatnonzero(keep)
    sub = Wbin[idx][:, idx].tocsr().astype(np.float64)
    rs = np.asarray(sub.sum(axis=1)).ravel()
    if np.any(rs <= 0):
        return None
    W = (diags(1.0 / rs) @ sub).tocsr()
    return W, keep, deg


def binary_star(Wbin: csr_matrix) -> csr_matrix:
    n = Wbin.shape[0]
    coo = Wbin.tocoo()
    rows = np.concatenate([coo.row, np.arange(n)])
    cols = np.concatenate([coo.col, np.arange(n)])
    data = np.ones(len(rows), dtype=np.float64)
    W = csr_matrix((data, (rows, cols)), shape=(n, n))
    W.sum_duplicates()
    W.data[:] = 1.0
    return W.tocsr()


def gi_star_z(y: np.ndarray, Wstar: csr_matrix) -> np.ndarray:
    """Ord/Getis Gi* z with binary star weights.

    Matches esda.G_Local(..., transform='B', star=True) when mean(y) > 0,
    and stays defined for mean-centered residuals.
    """
    y = np.asarray(y, dtype=np.float64)
    n = y.shape[0]
    ybar = float(y.mean())
    S = math.sqrt(float(np.sum((y - ybar) ** 2) / n))
    z = np.full(n, np.nan)
    if S < 1e-12:
        return z
    wy = np.asarray(Wstar @ y).ravel()
    card = np.asarray(Wstar.sum(axis=1)).ravel()
    inside = (n * card - card**2) / (n - 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        den = S * np.sqrt(inside)
        z = (wy - card * ybar) / den
    z[~np.isfinite(z)] = np.nan
    return z


def fisher_greater(hot: np.ndarray, cold: np.ndarray) -> tuple[float, float]:
    a = int(np.sum(hot & cold))
    b = int(np.sum(hot & ~cold))
    c = int(np.sum(~hot & cold))
    d = int(np.sum(~hot & ~cold))
    if min(a + b, c + d, a + c, b + d) == 0:
        return float("nan"), float("nan")
    oddsr, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
    return float(oddsr), float(p)


def overlap_block(z_cl: np.ndarray, z_cd: np.ndarray, prefix: str) -> dict:
    ok = np.isfinite(z_cl) & np.isfinite(z_cd)
    out = {}
    keys = [
        "n_hot",
        "n_cold",
        "n_overlap",
        "frac_overlap",
        "expected",
        "OR",
        "fisher_p",
        "frac_overlap10",
        "expected10",
        "OR10",
        "fisher_p10",
        "cd8_min",
        "cd8_q10",
    ]
    if ok.sum() < MIN_USED:
        for k in keys:
            out[f"{prefix}_{k}"] = float("nan")
        return out
    zc = z_cl[ok]
    zd = z_cd[ok]
    hot = zc > GI_Z
    cold = zd < -GI_Z
    q10 = float(np.quantile(zd, 0.10))
    cold10 = zd <= q10
    or_z, p_z = fisher_greater(hot, cold)
    or10, p10 = fisher_greater(hot, cold10)
    out[f"{prefix}_n_hot"] = int(hot.sum())
    out[f"{prefix}_n_cold"] = int(cold.sum())
    out[f"{prefix}_n_overlap"] = int(np.sum(hot & cold))
    out[f"{prefix}_frac_overlap"] = float(np.mean(hot & cold))
    out[f"{prefix}_expected"] = float(np.mean(hot) * np.mean(cold))
    out[f"{prefix}_OR"] = or_z
    out[f"{prefix}_fisher_p"] = p_z
    out[f"{prefix}_frac_overlap10"] = float(np.mean(hot & cold10))
    out[f"{prefix}_expected10"] = float(np.mean(hot) * np.mean(cold10))
    out[f"{prefix}_OR10"] = or10
    out[f"{prefix}_fisher_p10"] = p10
    out[f"{prefix}_cd8_min"] = float(np.min(zd))
    out[f"{prefix}_cd8_q10"] = q10
    return out


def _pearson(a: np.ndarray, b: np.ndarray) -> float | None:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a = a - a.mean()
    b = b - b.mean()
    da = float(np.sqrt(np.dot(a, a)))
    db = float(np.sqrt(np.dot(b, b)))
    if da < 1e-12 or db < 1e-12:
        return None
    return float(np.dot(a, b) / (da * db))


def tumor_cross(x: np.ndarray, lag: np.ndarray, nperm: int, rng: np.random.Generator) -> dict:
    """Pearson r of tumor-cell CLDN4 with the CD8A spatial lag. Null shuffles CLDN4."""
    empty = {
        "est": float("nan"),
        "se": float("nan"),
        "z": float("nan"),
        "p_less": float("nan"),
        "p_two": float("nan"),
        "n": int(len(x)),
    }
    if len(x) < MIN_TUMOR:
        return empty
    obs = _pearson(x, lag)
    if obs is None:
        return empty
    if nperm < 2:
        empty["est"] = obs
        return empty
    sims = np.empty(nperm)
    v = np.asarray(x, dtype=np.float64).copy()
    for i in range(nperm):
        rng.shuffle(v)
        val = _pearson(v, lag)
        sims[i] = np.nan if val is None else val
    sims = sims[np.isfinite(sims)]
    packed = pack_perm(obs, sims)
    packed["n"] = int(len(x))
    return packed


def nan_row(meta: dict, graph: str, status: str) -> dict:
    row = dict(meta)
    row["graph"] = graph
    row["status"] = status
    return row


def analyze_graph(xy_expr: np.ndarray, tumor: np.ndarray, Wbin: csr_matrix, graph: str, rng: np.random.Generator) -> dict:
    """xy_expr columns: CLDN4, CD8A, KRT8, EPCAM (log1p CP10k)."""
    prepared = row_standardize_kept(Wbin)
    giW = binary_star(Wbin)
    base = {
        "graph": graph,
        "n_cells": int(xy_expr.shape[0]),
        "median_degree": float(np.median(np.asarray(Wbin.sum(axis=1)).ravel())),
    }
    base.update(tumor_quartile_block(xy_expr, tumor, Wbin, N_PERM, rng))
    z_cl_raw = gi_star_z(xy_expr[:, 0], giW)
    z_cd_raw = gi_star_z(xy_expr[:, 1], giW)
    resid_cl, r2_cl = residualize(xy_expr[:, 0], xy_expr[:, 2:4])
    resid_cd, r2_cd = residualize(xy_expr[:, 1], xy_expr[:, 2:4])
    base["r2_cldn4"] = r2_cl
    base["r2_cd8a"] = r2_cd
    z_cl_r = gi_star_z(resid_cl, giW)
    z_cd_r = gi_star_z(resid_cd, giW)
    base.update(overlap_block(z_cl_raw, z_cd_raw, "gi_raw"))
    base.update(overlap_block(z_cl_r, z_cd_r, "gi_resid"))
    if prepared is None:
        base["status"] = "too_few_neighbors"
        return base
    W, keep, deg = prepared
    base["n_used"] = int(keep.sum())
    base["n_isolates"] = int((~keep).sum())
    base["median_degree"] = float(np.median(deg))
    E = xy_expr[keep]
    tmask = tumor[keep]
    cldn, cd8, krt, epc = E[:, 0], E[:, 1], E[:, 2], E[:, 3]
    resid_cl, r2_cl = residualize(cldn, np.column_stack([krt, epc]))
    resid_cd, r2_cd = residualize(cd8, np.column_stack([krt, epc]))
    base["r2_cldn4"] = r2_cl
    base["r2_cd8a"] = r2_cd
    base["spearman_raw"] = _safe_spearman(cldn, cd8)
    base["spearman_partial"] = _safe_spearman(resid_cl, resid_cd)

    zx = zscore(cldn, 1)
    zy = zscore(cd8, 1)
    zx0 = zscore(cldn, 0)
    zy0 = zscore(cd8, 0)
    zx_r = zscore(resid_cl, 1)
    zy_r = zscore(resid_cd, 1)
    zx_r0 = zscore(resid_cl, 0)
    zy_r0 = zscore(resid_cd, 0)
    zx_k = zscore(krt, 1)
    zx_e = zscore(epc, 1)
    if zx is None or zy is None or zx0 is None or zy0 is None:
        base["status"] = "zero_variance"
        return base

    nperm = N_PERM
    obs, sims, WY = shuffled_apply(zx, zy, W, nperm, rng)
    mor = pack_perm(obs, sims)
    put(base, "moran", mor)
    if WY is not None and zx_r is not None:
        o, s = moran_from_WY(zx_r, zy, W, WY)
        put(base, "cldn_only_moran", pack_perm(o, s))
    if WY is not None and zx_k is not None:
        o, s = moran_from_WY(zx_k, zy, W, WY)
        put(base, "krt8_moran", pack_perm(o, s))
    if WY is not None and zx_e is not None:
        o, s = moran_from_WY(zx_e, zy, W, WY)
        put(base, "epcam_moran", pack_perm(o, s))

    obs_l, sims_l, WY_l, _ = lee_apply(zx0, zy0, W, nperm, rng)
    put(base, "lee", pack_perm(obs_l, sims_l))
    if WY_l is not None and zx_r0 is not None:
        o, s = lee_from_WY(zx_r0, zy0, W, WY_l)
        put(base, "cldn_only_lee", pack_perm(o, s))

    if zx_r is not None and zy_r is not None:
        obs_p, sims_p, _ = shuffled_apply(zx_r, zy_r, W, nperm, rng)
        put(base, "partial_moran", pack_perm(obs_p, sims_p))
    if zx_r0 is not None and zy_r0 is not None:
        obs_pl, sims_pl, _, _ = lee_apply(zx_r0, zy_r0, W, nperm, rng)
        put(base, "partial_lee", pack_perm(obs_pl, sims_pl))

    # Reverse: CD8A with lag of CLDN4. Null shuffles CLDN4.
    obs_rev, sims_rev, _ = shuffled_apply(zy, zx, W, nperm, rng)
    put(base, "moran_rev", pack_perm(obs_rev, sims_rev))

    # Within-tumor cross-correlation. Lag of globally standardized CD8A.
    if int(tmask.sum()) >= MIN_TUMOR:
        lag_all = np.asarray(W @ zy).ravel()
        lag_t = lag_all[tmask]
        raw_t = tumor_cross(cldn[tmask], lag_t, nperm, rng)
        put(base, "tumor_raw", raw_t)
        base["tumor_n"] = raw_t["n"]
        rt, _ = residualize(cldn[tmask], np.column_stack([krt[tmask], epc[tmask]]))
        part_t = tumor_cross(rt, lag_t, nperm, rng)
        put(base, "tumor_partial", part_t)
    base["status"] = "ok"
    return base


def put(row: dict, prefix: str, packed: dict) -> None:
    name = {
        "moran": "moran_I",
        "lee": "lee_L",
        "cldn_only_moran": "cldn_only_moran_I",
        "krt8_moran": "krt8_moran_I",
        "epcam_moran": "epcam_moran_I",
        "cldn_only_lee": "cldn_only_lee_L",
        "partial_moran": "partial_moran_I",
        "partial_lee": "partial_lee_L",
        "moran_rev": "moran_rev_I",
        "tumor_raw": "tumor_raw_r",
        "tumor_partial": "tumor_partial_r",
        "tumor_q_diff": "tumor_q_diff",
    }[prefix]
    se_key, z_key, pl_key, pt_key = _keys(name)
    row[name] = packed["est"]
    row[se_key] = packed["se"]
    row[z_key] = packed["z"]
    row[pl_key] = packed["p_less"]
    row[pt_key] = packed.get("p_two", float("nan"))


def _se_name(name: str) -> str:
    return _keys(name)[0]


def _keys(name: str) -> tuple[str, str, str, str]:
    for suf in ("_I", "_L", "_r"):
        if name.endswith(suf):
            stem = name[: -len(suf)]
            return f"{stem}_se", f"{stem}_z", f"{stem}_p_less", f"{stem}_p_two"
    return f"{name}_se", f"{name}_z", f"{name}_p_less", f"{name}_p_two"


def _hi_lo_diff(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 40:
        return float("nan")
    q1, q3 = np.quantile(x, [0.25, 0.75])
    hi = x >= q3
    lo = x <= q1
    if hi.sum() == 0 or lo.sum() == 0:
        return float("nan")
    return float(y[hi].mean() - y[lo].mean())


def tumor_quartile_block(expr: np.ndarray, tumor: np.ndarray, Wbin: csr_matrix, nperm: int, rng: np.random.Generator) -> dict:
    """Tumor-cell CLDN4 top-minus-bottom quartile difference in mean neighbor CD8A."""
    deg = np.asarray(Wbin.sum(axis=1)).ravel()
    cd8 = np.asarray(expr[:, 1], dtype=np.float64)
    cldn = np.asarray(expr[:, 0], dtype=np.float64)
    lag = np.full(len(cd8), np.nan)
    nz = deg > 0
    if int(nz.sum()) == 0:
        return {}
    lag[nz] = np.asarray(Wbin @ cd8).ravel()[nz] / deg[nz]
    m = tumor & nz & np.isfinite(lag)
    out = {"tumor_q_n": int(m.sum())}
    if int(m.sum()) < 40:
        return out
    x = cldn[m]
    y = lag[m]
    obs = _hi_lo_diff(x, y)
    rt, _ = residualize(x, np.asarray(expr[m][:, 2:4], dtype=np.float64))
    obs_r = _hi_lo_diff(rt, y)
    out["tumor_q_hi"] = float(y[x >= np.quantile(x, 0.75)].mean())
    out["tumor_q_lo"] = float(y[x <= np.quantile(x, 0.25)].mean())
    if nperm < 2 or not np.isfinite(obs):
        out["tumor_q_diff"] = obs
        out["tumor_q_diff_partial"] = obs_r
        return out
    sims = np.empty(nperm)
    sims_r = np.empty(nperm)
    v = x.copy()
    vr = rt.copy()
    for i in range(nperm):
        rng.shuffle(v)
        sims[i] = _hi_lo_diff(v, y)
        rng.shuffle(vr)
        sims_r[i] = _hi_lo_diff(vr, y)
    packed = pack_perm(obs, sims[np.isfinite(sims)])
    packed_r = pack_perm(obs_r, sims_r[np.isfinite(sims_r)])
    put(out, "tumor_q_diff", packed)
    # partial difference uses its own column name, not the put() prefix map
    out["tumor_q_diff_partial"] = packed_r["est"]
    out["tumor_q_diff_partial_se"] = packed_r["se"]
    out["tumor_q_diff_partial_z"] = packed_r["z"]
    out["tumor_q_diff_partial_p_less"] = packed_r["p_less"]
    out["tumor_q_diff_partial_p_two"] = packed_r["p_two"]
    return out


def job_rng(unit_id: str, graph: str) -> np.random.Generator:
    h = zlib.adler32(f"{unit_id}|{graph}".encode())
    return np.random.default_rng(np.random.SeedSequence([SEED, h]))


def run_unit(expr: np.ndarray, tumor: np.ndarray, xy: np.ndarray, meta: dict) -> list[dict]:
    graphs = build_graphs(xy)
    rows = []
    for graph in GRAPHS:
        Wbin = graphs[graph]
        if Wbin is None:
            row = dict(meta)
            row.update(nan_row({}, graph, "graph_failed"))
            rows.append(row)
            continue
        rng = job_rng(meta["unit_id"], graph)
        try:
            stats_row = analyze_graph(expr, tumor, Wbin, graph, rng)
        except Exception as exc:
            stats_row = {"graph": graph, "status": "error", "error": repr(exc)}
        row = dict(meta)
        row.update(stats_row)
        rows.append(row)
    return rows


def work(job: dict) -> list[dict]:
    idx = job["idx"]
    meta = {k: v for k, v in job.items() if k != "idx"}
    try:
        return run_unit(LOGX[idx], TUMOR[idx], XY[idx], meta)
    except Exception as exc:
        return [{**meta, "graph": "all", "status": "error", "error": repr(exc)}]


def decode_arr(a) -> np.ndarray:
    out = []
    for x in a:
        if isinstance(x, bytes):
            out.append(x.decode())
        else:
            out.append(str(x))
    return np.asarray(out)


def extract(h5ad: Path, cache: Path) -> dict:
    import h5py

    print(f"[extract] {h5ad}", flush=True)
    f = h5py.File(h5ad, "r")
    genes = decode_arr(f["var/_index"][:])
    missing = [g for g in GENES if g not in set(genes)]
    if missing:
        raise SystemExit(f"panel missing {missing}")
    gix = [int(np.where(genes == g)[0][0]) for g in GENES]
    indptr = f["layers/counts/indptr"][:]
    indices = f["layers/counts/indices"][:]
    data = f["layers/counts/data"][:]
    n = len(indptr) - 1
    counts = np.zeros((n, 4), dtype=np.float32)
    for j, gi in enumerate(gix):
        m = indices == gi
        rows = np.searchsorted(indptr, np.flatnonzero(m), side="right") - 1
        counts[rows, j] = data[m]
    nc = f["obs/n_counts"][:].astype(np.float32)
    scale = np.zeros(n, dtype=np.float32)
    np.divide(1e4, nc, out=scale, where=nc > 0)
    logx = np.log1p(counts * scale[:, None]).astype(np.float32)
    xy = f["obsm/spatial"][:].astype(np.float64)

    def cats(path):
        names = decode_arr(f[path + "/categories"][:])
        codes = f[path + "/codes"][:]
        return names[codes]

    sample = cats("obs/sample")
    patient = cats("obs/patient")
    fov = f["obs/fov"][:].astype(np.int32)
    ct = cats("obs/cell_type")
    tumor = np.isin(ct, ["tumor 5", "tumor 6", "tumor 9", "tumor 12", "tumor 13"])
    f.close()
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache,
        version=np.array([CACHE_VERSION]),
        xy=xy,
        logx=logx,
        tumor=tumor.astype(np.uint8),
        sample=sample.astype("U32"),
        patient=patient.astype("U16"),
        fov=fov,
        n_counts=nc,
        genes=np.array(GENES),
    )
    print(f"[cache] {cache} cells={n}", flush=True)
    return load_cache(cache)


def load_cache(cache: Path) -> dict:
    z = np.load(cache, allow_pickle=False)
    if int(z["version"][0]) != CACHE_VERSION:
        raise SystemExit(f"stale cache {cache}")
    return {
        "xy": z["xy"],
        "logx": z["logx"],
        "tumor": z["tumor"].astype(bool),
        "sample": z["sample"].astype(str),
        "patient": z["patient"].astype(str),
        "fov": z["fov"],
        "n_counts": z["n_counts"],
    }


def histology_of(sample: str) -> str:
    return "LUSC" if "LUSC" in sample else "LUAD"


def build_jobs(arr: dict, max_fovs: int, skip_sample: bool) -> tuple[list[dict], list[dict]]:
    sample, patient, fov = arr["sample"], arr["patient"], arr["fov"]
    frame = pd.DataFrame({"sample": sample, "patient": patient, "fov": fov})
    fov_jobs = []
    for (s, fv), idx in frame.groupby(["sample", "fov"], sort=True).groups.items():
        ix = np.asarray(idx, dtype=np.int64)
        pat = str(patient[ix[0]])
        fov_jobs.append(
            {
                "level": "fov",
                "unit_id": f"{s}__fov{int(fv)}",
                "sample": str(s),
                "patient": pat,
                "histology": histology_of(str(s)),
                "fov": int(fv),
                "n_cells_unit": int(ix.size),
                "n_cldn4_pos": int(np.sum(arr["logx"][ix, 0] > 0)),
                "n_cd8a_pos": int(np.sum(arr["logx"][ix, 1] > 0)),
                "n_tumor": int(np.sum(arr["tumor"][ix])),
                **_detect_fracs(arr, ix),
                "idx": ix,
            }
        )
    if max_fovs and max_fovs < len(fov_jobs):
        fov_jobs = fov_jobs[:max_fovs]
    sample_jobs = []
    if not skip_sample:
        for s, idx in frame.groupby("sample", sort=True).groups.items():
            ix = np.asarray(idx, dtype=np.int64)
            pat = str(patient[ix[0]])
            sample_jobs.append(
                {
                    "level": "sample",
                    "unit_id": str(s),
                    "sample": str(s),
                    "patient": pat,
                    "histology": histology_of(str(s)),
                    "fov": -1,
                    "n_cells_unit": int(ix.size),
                "n_cldn4_pos": int(np.sum(arr["logx"][ix, 0] > 0)),
                "n_cd8a_pos": int(np.sum(arr["logx"][ix, 1] > 0)),
                "n_tumor": int(np.sum(arr["tumor"][ix])),
                **_detect_fracs(arr, ix),
                "idx": ix,
            }
        )
    return fov_jobs, sample_jobs


def _detect_fracs(arr: dict, ix: np.ndarray) -> dict:
    logx = arr["logx"][ix]
    tumor = arr["tumor"][ix]
    cldn = logx[:, 0] > 0
    cd8 = logx[:, 1] > 0
    out = {
        "double_pos_frac": float(np.mean(cldn & cd8)),
        "cd8a_pos_frac_tumor": float(np.mean(cd8[tumor])) if tumor.any() else float("nan"),
        "cd8a_pos_frac_other": float(np.mean(cd8[~tumor])) if (~tumor).any() else float("nan"),
    }
    return out


def inverse_variance_meta(est, se) -> dict:
    i = np.asarray(est, dtype=float)
    s = np.asarray(se, dtype=float)
    ok = np.isfinite(i) & np.isfinite(s) & (s > 0)
    i, s = i[ok], s[ok]
    out = {
        "n": int(i.size),
        "I_fixed": float("nan"),
        "I_random": float("nan"),
        "se_fixed": float("nan"),
        "se_random": float("nan"),
        "z": float("nan"),
        "p": float("nan"),
        "Q": float("nan"),
        "tau2": float("nan"),
        "I2": float("nan"),
    }
    if i.size < 2:
        return out
    w = 1.0 / s**2
    Ihat = float(np.sum(w * i) / np.sum(w))
    se_f = float(math.sqrt(1.0 / np.sum(w)))
    Q = float(np.sum(w * (i - Ihat) ** 2))
    df = i.size - 1
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    w_re = 1.0 / (s**2 + tau2)
    Ire = float(np.sum(w_re * i) / np.sum(w_re))
    se_re = float(math.sqrt(1.0 / np.sum(w_re)))
    z = Ire / se_re if se_re > 0 else float("nan")
    p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else float("nan")
    out.update(
        {
            "I_fixed": Ihat,
            "I_random": Ire,
            "se_fixed": se_f,
            "se_random": se_re,
            "z": float(z),
            "p": p,
            "Q": Q,
            "tau2": float(tau2),
            "I2": float(max(0.0, (Q - df) / Q)) if Q > 0 else 0.0,
        }
    )
    return out


def bh_q(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    m = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    qv = ranked * m / (np.arange(1, m + 1))
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    qv = np.clip(qv, 0, 1)
    back = np.empty(m)
    back[order] = qv
    q[ok] = back
    return q


def wilcoxon_less(v: np.ndarray) -> dict:
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    out = {"n": int(v.size), "n_neg": int(np.sum(v < 0)), "median": float(np.median(v)) if v.size else float("nan")}
    if v.size < 5 or np.allclose(v, 0):
        out["wilcoxon_p"] = float("nan")
        out["wilcoxon_stat"] = float("nan")
        return out
    try:
        w = stats.wilcoxon(v, alternative="less", zero_method="wilcox")
        out["wilcoxon_stat"] = float(w.statistic)
        out["wilcoxon_p"] = float(w.pvalue)
    except ValueError:
        out["wilcoxon_stat"] = float("nan")
        out["wilcoxon_p"] = float("nan")
    return out


def sign_test(v: np.ndarray) -> dict:
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v) & (v != 0)]
    n = int(v.size)
    k = int(np.sum(v < 0))
    if n == 0:
        return {"n": 0, "n_neg": 0, "p": float("nan")}
    p = float(stats.binomtest(k, n, 0.5, alternative="greater").pvalue)
    return {"n": n, "n_neg": k, "p": p}


def section_patient_rollups(df: pd.DataFrame, value: str) -> dict:
    ok = df[df["status"] == "ok"].copy()
    if ok.empty or value not in ok.columns:
        return {}
    sec = ok.groupby("sample", sort=False)[value].median()
    pat_map = ok.groupby("sample")["patient"].first()
    hist_map = ok.groupby("sample")["histology"].first()
    sec_df = pd.DataFrame({"estimate": sec, "patient": pat_map, "histology": hist_map})
    pat = sec_df.groupby("patient")["estimate"].mean()
    return {
        "fov": wilcoxon_less(ok[value].to_numpy()),
        "section_values": {str(k): float(v) for k, v in sec.items()},
        "section": wilcoxon_less(sec.to_numpy()),
        "section_sign": sign_test(sec.to_numpy()),
        "patient_values": {str(k): float(v) for k, v in pat.items()},
        "patient": wilcoxon_less(pat.to_numpy()),
        "patient_sign": sign_test(pat.to_numpy()),
        "meta_fov": inverse_variance_meta(ok[value], ok[_se_name(value)]) if _se_name(value) in ok.columns else {},
        "luad_fov": wilcoxon_less(ok.loc[ok["histology"] == "LUAD", value].to_numpy()),
        "lusc_median": float(ok.loc[ok["histology"] == "LUSC", value].median()) if (ok["histology"] == "LUSC").any() else float("nan"),
    }


def add_fdr(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "partial_moran_p_less" not in df.columns:
        return df
    q = np.full(len(df), np.nan)
    for graph, idx in df.groupby("graph").groups.items():
        ix = np.asarray(list(idx))
        q[ix] = bh_q(df.loc[ix, "partial_moran_p_less"].to_numpy())
    df["partial_moran_q_less"] = q
    return df


def fmt(x, nd=4) -> str:
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(x):
        return "NA"
    return f"{x:+.{nd}f}"


def fmt_p(x) -> str:
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def md_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for r in rows:
        lines.append("| " + " | ".join(r) + " |")
    return "\n".join(lines)


def _roll_sentence(label: str, roll: dict) -> str:
    if not roll:
        return f"{label}: not available."
    fov = roll.get("fov", {})
    sec = roll.get("section", {})
    psign = roll.get("patient_sign", {})
    return (
        f"{label}: median FOV estimate {fmt(fov.get('median'))} "
        f"({fov.get('n_neg', 'NA')}/{fov.get('n', 'NA')} FOVs negative; "
        f"FOV Wilcoxon p={fmt_p(fov.get('wilcoxon_p'))}). "
        f"Section medians negative in {sec.get('n_neg', 'NA')}/{sec.get('n', 'NA')} "
        f"(Wilcoxon p={fmt_p(sec.get('wilcoxon_p'))}). "
        f"Patients negative in {psign.get('n_neg', 'NA')}/{psign.get('n', 'NA')} "
        f"(one-sided sign p={fmt_p(psign.get('p'))})."
    )


def conclusion_text(roll_partial: dict, roll_raw: dict, roll_tumor: dict) -> str:
    psign = (roll_partial or {}).get("patient_sign", {})
    sec = (roll_partial or {}).get("section", {})
    med_p = (roll_partial or {}).get("fov", {}).get("median", float("nan"))
    med_r = (roll_raw or {}).get("fov", {}).get("median", float("nan"))
    n_pat_neg = psign.get("n_neg", 0)
    n_pat = psign.get("n", 0)
    sec_p = sec.get("wilcoxon_p", 1)
    lines = []
    lines.append(_roll_sentence("Primary partial Moran I (radius 50 µm, both genes residualized on KRT8+EPCAM)", roll_partial))
    lines.append(_roll_sentence("Raw Moran I on the same graph (no residualization)", roll_raw))
    lines.append(_roll_sentence("Within-tumor partial cross-correlation on the same graph", roll_tumor))
    if n_pat and n_pat_neg == n_pat and np.isfinite(sec_p) and sec_p < 0.05 and np.isfinite(med_p) and med_p < 0:
        lines.append(
            "Partial Moran I is negative in every patient, and the section-level Wilcoxon test is in the same direction. "
            "On this public CosMx cohort, CLDN4 remains spatially anti-associated with CD8A after linear removal of KRT8 and EPCAM. "
            "This is neighbor-graph autocorrelation of the expression field. It sits beside the locked CLDN4-high versus CLDN4-low tumor-cell neighbor-count result; it does not replace that count."
        )
    elif np.isfinite(med_p) and med_p < 0 and n_pat_neg >= max(1, n_pat - 1):
        lines.append(
            "Partial Moran I leans negative, but the patient panel is not uniformly negative or the section-level test is weak. "
            "The table below is the result. A uniform CLDN4-specific anti-association, after KRT8 and EPCAM removal, is not claimed from this test."
        )
    else:
        lines.append(
            "Partial Moran I does not support a consistent CLDN4-specific spatial anti-association with CD8A after KRT8 and EPCAM are removed. "
            "A negative raw Moran I, when present, is the spatial complement of CLDN4 marking epithelium while CD8A marks T cells."
        )
    if np.isfinite(med_r) and np.isfinite(med_p):
        lines.append(
            f"Median raw I is {fmt(med_r)} and median partial I is {fmt(med_p)} on the primary radius. "
            "The gap between them is the part of the spatial association shared with KRT8 and EPCAM."
        )
    return "\n\n".join(lines)


def graph_summary_table(df: pd.DataFrame, value: str) -> str:
    headers = ["graph", "n", "n_neg", "median", "section_neg", "section_p", "patient_neg", "patient_sign_p", "RE", "I2"]
    rows = []
    for graph in GRAPHS:
        sub = df[(df["graph"] == graph) & (df["status"] == "ok")]
        if sub.empty or value not in sub.columns:
            continue
        roll = section_patient_rollups(sub, value)
        fov = roll.get("fov", {})
        sec = roll.get("section", {})
        psign = roll.get("patient_sign", {})
        meta = roll.get("meta_fov", {})
        rows.append(
            [
                graph,
                str(fov.get("n", 0)),
                str(fov.get("n_neg", 0)),
                fmt(fov.get("median")),
                f"{sec.get('n_neg', 'NA')}/{sec.get('n', 'NA')}",
                fmt_p(sec.get("wilcoxon_p")),
                f"{psign.get('n_neg', 'NA')}/{psign.get('n', 'NA')}",
                fmt_p(psign.get("p")),
                fmt(meta.get("I_random")),
                fmt(meta.get("I2"), 3) if np.isfinite(meta.get("I2", np.nan)) else "NA",
            ]
        )
    return md_table(headers, rows)


def value_table(records: dict, title_keys: list[str]) -> str:
    headers = ["unit"] + title_keys
    rows = []
    units = []
    for key in title_keys:
        units.extend((records.get(key) or {}).keys())
    # use patient or sample order from the first non-empty
    order = []
    for key in title_keys:
        for u in (records.get(key) or {}):
            if u not in order:
                order.append(u)
    for u in order:
        row = [u]
        for key in title_keys:
            row.append(fmt((records.get(key) or {}).get(u, float("nan"))))
        rows.append(row)
    return md_table(headers, rows)


def write_results(fov: pd.DataFrame, sample: pd.DataFrame, headline: dict, out: Path) -> None:
    prim = fov[(fov["graph"] == "radius_50um") & (fov["status"] == "ok")]
    lines = []
    lines.append("# RESULTS — CosMx bivariate Moran / Lee / Gi* of CLDN4 vs CD8A")
    lines.append("")
    lines.append("**Dataset.** Figshare 25976224, `cosmx_human_nsclc_clustered.h5ad` (He et al. 2022 CosMx NSCLC; CellCharter object).")
    lines.append(f"Cells in the object: **{headline['n_cells']}**. Sample–FOV units analyzed: **{headline['n_fov']}**. Sections: **{headline['n_samples']}**. Patients: **{headline['n_patients']}**.")
    lines.append("Expression: log1p(counts / n_counts × 10,000) from `layers/counts`. Coordinates: `obsm['spatial']` global pixels, 0.18 µm/px.")
    lines.append("No private 8-KL data. CLDN4 only (no TACSTD2 gate).")
    lines.append("")
    lines.append("## Estimand")
    lines.append("")
    lines.append("Primary: partial bivariate Moran's I, CLDN4 with the row-standardized lag of CD8A, both genes residualized by within-FOV OLS on KRT8 and EPCAM. Graph: cells within 50 µm. Co-primary graph: kNN k=15.")
    lines.append("Lee's L uses the same residualized fields. Gi* uses binary star weights on the same neighbor sets.")
    lines.append("Within-tumor cross-correlation residualizes CLDN4 on KRT8 and EPCAM inside the tumor-cell set (author labels tumor 5/6/9/12/13) and correlates that residual with the lag of CD8A. The null shuffles the residual among tumor cells.")
    lines.append("KRT8–CD8A and EPCAM–CD8A Moran use the same CD8A permutations as raw CLDN4–CD8A Moran and are the epithelial baselines.")
    lines.append("")
    lines.append("Permutation p-values use 199 within-unit shuffles. Section Wilcoxon and the patient sign test are the confirmatory summaries. FOV inverse-variance p-values assume independent FOVs and are reported as descriptive effect sizes (RE = DerSimonian–Laird).")
    lines.append("")
    lines.append("## Headline")
    lines.append("")
    lines.append(headline["prose"])
    lines.append("")
    lines.append("## Partial Moran I across neighbor graphs")
    lines.append("")
    lines.append(graph_summary_table(fov, "partial_moran_I"))
    lines.append("")
    lines.append("## Raw Moran I (CLDN4, lag CD8A)")
    lines.append("")
    lines.append(graph_summary_table(fov, "moran_I"))
    lines.append("")
    lines.append("## Lee's L, partial")
    lines.append("")
    lines.append(graph_summary_table(fov, "partial_lee_L"))
    lines.append("")
    lines.append("## CLDN4 residual only (CD8A not residualized)")
    lines.append("")
    lines.append(graph_summary_table(fov, "cldn_only_moran_I"))
    lines.append("")
    lines.append("## Epithelial baselines (raw Moran)")
    lines.append("")
    lines.append("KRT8 with lag CD8A:")
    lines.append("")
    lines.append(graph_summary_table(fov, "krt8_moran_I"))
    lines.append("")
    lines.append("EPCAM with lag CD8A:")
    lines.append("")
    lines.append(graph_summary_table(fov, "epcam_moran_I"))
    lines.append("")
    lines.append("## Within-tumor CLDN4 quartile vs neighbor CD8A")
    lines.append("")
    lines.append("Among tumor cells with at least one neighbor, mean neighbor CD8A in the top CLDN4 quartile minus the bottom quartile. Negative means CLDN4-high tumor cells have less CD8A around them. `tumor_q_diff_partial` uses CLDN4 residualized on KRT8 and EPCAM inside the tumor set.")
    lines.append("")
    lines.append(graph_summary_table(fov, "tumor_q_diff"))
    lines.append("")
    lines.append(graph_summary_table(fov, "tumor_q_diff_partial"))
    lines.append("")
    if len(prim) and "cd8a_pos_frac_tumor" in prim.columns:
        lines.append(
            f"CD8A detection (expression > 0), median across primary-radius FOVs: "
            f"**{float(prim['cd8a_pos_frac_tumor'].median()):.3f}** of tumor cells and "
            f"**{float(prim['cd8a_pos_frac_other'].median()):.3f}** of other cells. "
            f"Median fraction of cells positive for both CLDN4 and CD8A: **{float(prim['double_pos_frac'].median()):.3f}**."
        )
        lines.append("")
    lines.append("## Within-tumor partial cross-correlation")
    lines.append("")
    lines.append("Pearson r of tumor-cell CLDN4 with the CD8A lag. This is not a Moran I.")
    lines.append("")
    lines.append(graph_summary_table(fov, "tumor_partial_r"))
    lines.append("")
    lines.append("Within-tumor raw CLDN4 (no residualization) is in `tables/fov_graph_stats.tsv` (`tumor_raw_r`). The within-tumor number is a Pearson correlation, not a Moran I.")
    lines.append("")
    lines.append("## Primary graph by section and patient")
    lines.append("")
    lines.append("Values are medians of FOV-level estimates.")
    lines.append("")
    roll_p = headline["rolls"]["partial_moran_I"]
    roll_r = headline["rolls"]["moran_I"]
    roll_l = headline["rolls"]["partial_lee_L"]
    roll_t = headline["rolls"]["tumor_partial_r"]
    roll_k = headline["rolls"]["krt8_moran_I"]
    lines.append(
        value_table(
            {
                "partial_I": roll_p.get("section_values", {}),
                "raw_I": roll_r.get("section_values", {}),
                "partial_L": roll_l.get("section_values", {}),
                "tumor_partial_r": roll_t.get("section_values", {}),
                "krt8_I": roll_k.get("section_values", {}),
            },
            ["partial_I", "raw_I", "partial_L", "tumor_partial_r", "krt8_I"],
        )
    )
    lines.append("")
    lines.append("Patient means of those section medians:")
    lines.append("")
    lines.append(
        value_table(
            {
                "partial_I": roll_p.get("patient_values", {}),
                "raw_I": roll_r.get("patient_values", {}),
                "partial_L": roll_l.get("patient_values", {}),
                "tumor_partial_r": roll_t.get("patient_values", {}),
                "krt8_I": roll_k.get("patient_values", {}),
            },
            ["partial_I", "raw_I", "partial_L", "tumor_partial_r", "krt8_I"],
        )
    )
    lines.append("")
    r2 = float(prim["r2_cldn4"].median()) if "r2_cldn4" in prim.columns and len(prim) else float("nan")
    r2b = float(prim["r2_cd8a"].median()) if "r2_cd8a" in prim.columns and len(prim) else float("nan")
    sp = float(prim["spearman_raw"].median()) if "spearman_raw" in prim.columns and len(prim) else float("nan")
    spp = float(prim["spearman_partial"].median()) if "spearman_partial" in prim.columns and len(prim) else float("nan")
    lines.append(
        f"On the primary radius, median within-FOV R² of CLDN4 on KRT8+EPCAM is **{r2:.3f}**, "
        f"and of CD8A on KRT8+EPCAM is **{r2b:.3f}**. "
        f"Median same-cell Spearman CLDN4 vs CD8A is **{sp:+.3f}** (raw) and **{spp:+.3f}** (residuals). "
        "Same-cell Spearman is co-expression inside a cell. Moran's I is the neighbor lag."
    )
    lines.append("")
    lines.append("## Gi* overlap")
    lines.append("")
    lines.append(gi_paragraph(fov))
    lines.append("")
    lines.append("## Whole-section graphs")
    lines.append("")
    if sample is not None and len(sample):
        lines.append("The same statistics on one graph per tissue section (global coordinates, so neighbors can cross FOV borders). Each row is a section, not a FOV.")
        lines.append("")
        lines.append(graph_summary_table(sample.rename(columns={"sample": "sample"}), "partial_moran_I") if False else _sample_graph_table(sample))
        lines.append("")
        lines.append(_sample_primary_sentence(sample))
    else:
        lines.append("Section-level graphs were not run.")
    lines.append("")
    lines.append("## What this layer is")
    lines.append("")
    lines.append("An expression-field autocorrelation layer on the public He et al. 2022 CosMx cohort. The locked cell-type result (CLDN4-high tumor neighborhoods contain fewer cytotoxic cells at 50/100 µm; exclusion without loss of GZMB/PRF1/NKG7/IFNG) is a different estimand and is not recomputed here.")
    lines.append("Lung6 is the LUSC section. Lung5 and Lung9 contribute serial sections; patient summaries average those sections so technical replicates do not count as extra patients.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("bash methods/cosmx_bivar_autocorr/download.sh")
    lines.append("python3 methods/cosmx_bivar_autocorr/analyze.py --h5ad data/cosmx_human_nsclc_clustered.h5ad")
    lines.append("```")
    lines.append("")
    lines.append("Outputs: `methods/cosmx_bivar_autocorr/tables/` and `methods/cosmx_bivar_autocorr/figures/`.")
    out.write_text("\n".join(lines) + "\n")


def gi_paragraph(fov: pd.DataFrame) -> str:
    chunks = []
    for graph in PRIMARY:
        sub = fov[(fov["graph"] == graph) & (fov["status"] == "ok")]
        if sub.empty or "gi_resid_frac_overlap" not in sub.columns:
            continue
        d = (sub["gi_resid_frac_overlap"] - sub["gi_resid_expected"]).to_numpy()
        w = wilcoxon_less(-d)  # less on -d means d>0 if we flip... use greater via less on negative
        # Report Wilcoxon on (observed - expected), alternative greater, by testing -delta < 0.
        try:
            delta = sub["gi_resid_frac_overlap"] - sub["gi_resid_expected"]
            delta = delta[np.isfinite(delta)]
            ww = stats.wilcoxon(delta, alternative="greater", zero_method="wilcox") if len(delta) >= 5 else None
        except ValueError:
            ww = None
        n_sig = int(np.sum(sub["gi_resid_fisher_p"] < 0.05)) if "gi_resid_fisher_p" in sub else 0
        n_cold = int(np.sum(sub["gi_resid_n_cold"] > 0)) if "gi_resid_n_cold" in sub else 0
        med_min = float(sub["gi_resid_cd8_min"].median()) if "gi_resid_cd8_min" in sub else float("nan")
        med_frac = float(sub["gi_resid_frac_overlap"].median())
        med_exp = float(sub["gi_resid_expected"].median())
        med_or = float(sub["gi_resid_OR"].median())
        med_frac10 = float(sub["gi_raw_frac_overlap10"].median()) if "gi_raw_frac_overlap10" in sub else float("nan")
        med_exp10 = float(sub["gi_raw_expected10"].median()) if "gi_raw_expected10" in sub else float("nan")
        ptxt = fmt_p(ww.pvalue) if ww is not None else "NA"
        chunks.append(
            f"**{graph}.** Residual Gi* (CLDN4 z>1.96 and CD8A z<−1.96): median overlap fraction {med_frac:.4f} "
            f"versus independence {med_exp:.4f} (median OR {med_or:.3f}). "
            f"FOVs with any analytic CD8 cold cells: {n_cold}/{len(sub)}. "
            f"Fisher p<0.05 in {n_sig}/{len(sub)} FOVs. "
            f"Wilcoxon on (overlap − expected), alternative greater: p={ptxt}. "
            f"Median minimum residual CD8 Gi* z={med_min:.2f}. "
            f"Raw-expression rank cold (CD8 Gi* at or below the FOV 10th percentile) overlap median {med_frac10:.4f} "
            f"versus independence {med_exp10:.4f}."
        )
    if not chunks:
        return "Gi* summaries were not produced."
    chunks.append(
        "Analytic CD8 cold spots can be scarce when CD8A is sparse, because the local mean cannot fall far below a near-zero background. "
        "The rank-based cold set is the check for that limit. Overlap is a hotspot coincidence, not a contact probability."
    )
    return "\n\n".join(chunks)


def _sample_graph_table(sample: pd.DataFrame) -> str:
    # Reuse FOV table logic: each sample-graph row is one "FOV" in that helper's language.
    # section rollup then collapses to one value per sample, and patient rollup still works.
    return "\n\n".join(
        [
            "Partial Moran I:",
            graph_summary_table(sample, "partial_moran_I"),
            "Raw Moran I:",
            graph_summary_table(sample, "moran_I"),
            "Within-tumor partial Pearson r:",
            graph_summary_table(sample, "tumor_partial_r"),
        ]
    )


def _sample_primary_sentence(sample: pd.DataFrame) -> str:
    sub = sample[(sample["graph"] == "radius_50um") & (sample["status"] == "ok")]
    if sub.empty:
        return "Primary radius was not available at section level."
    roll = section_patient_rollups(sub, "partial_moran_I")
    # For section-level rows, "FOV" in the helper is actually the section.
    return (
        "On the whole-section 50 µm graph, partial Moran I "
        + _roll_sentence("section-as-unit", roll)
        + " Because each section is a single test, the FOV count in that sentence is the section count."
    )


def make_plots(fov: pd.DataFrame, sample: pd.DataFrame, arr: dict, fig_dir: Path) -> None:
    fig_dir.mkdir(parents=True, exist_ok=True)
    prim = fov[(fov["graph"].isin(PRIMARY)) & (fov["status"] == "ok")].copy()
    if prim.empty:
        return
    _strip(prim, fig_dir / "fov_partial_moran_by_patient")
    _scatter_raw_partial(fov, fig_dir / "raw_vs_partial_moran")
    _heatmap(fov, fig_dir / "graph_sensitivity_heatmap")
    _forest(fov, sample, fig_dir / "section_forest_partial_moran")
    _gi_maps(fov, arr, fig_dir)


def _strip(prim: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=True)
    rng = np.random.default_rng(1)
    for ax, graph in zip(axes, PRIMARY):
        sub = prim[prim["graph"] == graph]
        for i, pat in enumerate(PATIENT_ORDER):
            s = sub[sub["patient"] == pat]
            if s.empty:
                continue
            jit = rng.uniform(-0.18, 0.18, len(s))
            ax.scatter(
                np.full(len(s), i) + jit,
                s["partial_moran_I"],
                s=12,
                c=PATIENT_COLOR[pat],
                alpha=0.75,
                linewidths=0,
                label=pat,
            )
        ax.axhline(0, color="#333333", lw=0.8, ls="--")
        ax.set_xticks(range(len(PATIENT_ORDER)))
        ax.set_xticklabels(PATIENT_ORDER, fontsize=8)
        ax.set_title(graph)
        ax.set_xlabel("Patient")
    axes[0].set_ylabel("Partial Moran I  (CLDN4, lag CD8A | KRT8, EPCAM)")
    axes[1].legend(frameon=False, fontsize=8, loc="best")
    fig.suptitle("FOV partial bivariate Moran's I", fontsize=12)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _scatter_raw_partial(fov: pd.DataFrame, path: Path) -> None:
    sub = fov[(fov["graph"] == "radius_50um") & (fov["status"] == "ok")]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    for pat in PATIENT_ORDER:
        s = sub[sub["patient"] == pat]
        if s.empty:
            continue
        ax.scatter(s["moran_I"], s["partial_moran_I"], s=16, c=PATIENT_COLOR[pat], alpha=0.8, linewidths=0, label=pat)
    lims = np.nanpercentile(np.concatenate([sub["moran_I"], sub["partial_moran_I"]]), [1, 99])
    pad = 0.02
    lo, hi = float(lims[0] - pad), float(lims[1] + pad)
    ax.plot([lo, hi], [lo, hi], color="#888888", lw=0.8, ls=":")
    ax.axhline(0, color="#333333", lw=0.6, ls="--")
    ax.axvline(0, color="#333333", lw=0.6, ls="--")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Raw Moran I")
    ax.set_ylabel("Partial Moran I (KRT8 + EPCAM)")
    ax.set_title("Radius 50 µm, one point per FOV")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _heatmap(fov: pd.DataFrame, path: Path) -> None:
    cols = [
        ("moran_I", "CLDN4 raw"),
        ("partial_moran_I", "CLDN4 partial"),
        ("cldn_only_moran_I", "CLDN4 resid only"),
        ("krt8_moran_I", "KRT8 raw"),
        ("epcam_moran_I", "EPCAM raw"),
        ("partial_lee_L", "Lee partial"),
        ("tumor_partial_r", "Tumor partial r"),
    ]
    mat = np.full((len(GRAPHS), len(cols)), np.nan)
    for i, g in enumerate(GRAPHS):
        sub = fov[(fov["graph"] == g) & (fov["status"] == "ok")]
        for j, (col, _) in enumerate(cols):
            if col in sub.columns and sub[col].notna().any():
                mat[i, j] = float(sub[col].median())
    finite = mat[np.isfinite(mat)]
    vmax = float(np.max(np.abs(finite))) if finite.size else 0.1
    vmax = max(vmax, 0.01)
    fig, ax = plt.subplots(figsize=(9.4, 4.6))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[1] for c in cols], rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(len(GRAPHS)))
    ax.set_yticklabels(GRAPHS, fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:+.3f}", ha="center", va="center", fontsize=7, color="black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Median FOV estimate")
    ax.set_title("Neighbor-graph sensitivity")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _forest(fov: pd.DataFrame, sample: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.4), sharey=False)
    sub = fov[(fov["graph"] == "radius_50um") & (fov["status"] == "ok")]
    if not sub.empty:
        sec = sub.groupby("sample")["partial_moran_I"].median()
        q1 = sub.groupby("sample")["partial_moran_I"].quantile(0.25)
        q3 = sub.groupby("sample")["partial_moran_I"].quantile(0.75)
        order = list(sec.sort_values().index)
        y = np.arange(len(order))
        est = np.array([sec[s] for s in order])
        lo = est - np.array([q1[s] for s in order])
        hi = np.array([q3[s] for s in order]) - est
        axes[0].errorbar(est, y, xerr=[lo, hi], fmt="o", color="#222222", ecolor="#666666", ms=5)
        axes[0].axvline(0, color="#b2182b", lw=0.8, ls="--")
        axes[0].set_yticks(y)
        axes[0].set_yticklabels(order, fontsize=8)
        axes[0].set_xlabel("FOV median partial I (bar = IQR)")
        axes[0].set_title("Within-FOV 50 µm graphs")
    if sample is not None and len(sample):
        ssub = sample[(sample["graph"] == "radius_50um") & (sample["status"] == "ok")].copy()
        if not ssub.empty:
            ssub = ssub.sort_values("partial_moran_I")
            y = np.arange(len(ssub))
            axes[1].errorbar(
                ssub["partial_moran_I"],
                y,
                xerr=1.96 * ssub["partial_moran_se"],
                fmt="o",
                color="#222222",
                ecolor="#666666",
                ms=5,
            )
            axes[1].axvline(0, color="#b2182b", lw=0.8, ls="--")
            axes[1].set_yticks(y)
            axes[1].set_yticklabels(ssub["sample"], fontsize=8)
            axes[1].set_xlabel("Section partial I (bar = 1.96 permutation SE)")
            axes[1].set_title("Whole-section 50 µm graph")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=150)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def _gi_maps(fov: pd.DataFrame, arr: dict, fig_dir: Path) -> None:
    sub = fov[(fov["graph"] == "radius_50um") & (fov["status"] == "ok") & (fov["n_cd8a_pos"] >= 30)]
    if sub.empty or "partial_moran_I" not in sub.columns:
        return
    sub = sub[np.isfinite(sub["partial_moran_I"])]
    if sub.empty:
        return
    picks = []
    neg = sub.loc[sub["partial_moran_I"].idxmin()]
    pos = sub.loc[sub["partial_moran_I"].idxmax()]
    picks.append((neg, "lowest_partial_I"))
    if pos["unit_id"] != neg["unit_id"]:
        picks.append((pos, "highest_partial_I"))
    sample = arr["sample"]
    fov_id = arr["fov"]
    for rec, tag in picks:
        m = (sample == rec["sample"]) & (fov_id == int(rec["fov"]))
        xy = arr["xy"][m]
        expr = arr["logx"][m]
        tumor = arr["tumor"][m]
        Wbin = radius_binary(xy, 50.0 / PX_TO_UM)
        Wstar = binary_star(Wbin)
        resid_cl, _ = residualize(expr[:, 0], expr[:, 2:4])
        resid_cd, _ = residualize(expr[:, 1], expr[:, 2:4])
        zc = gi_star_z(resid_cl, Wstar)
        zd = gi_star_z(resid_cd, Wstar)
        hot = zc > GI_Z
        cold = zd < -GI_Z
        cat = np.zeros(len(xy), dtype=int)
        cat[hot] = 1
        cat[cold] = 2
        cat[hot & cold] = 3
        um = (xy - xy.min(axis=0)) * PX_TO_UM
        fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8))
        panels = [
            (expr[:, 0], "CLDN4 log1p CP10k", "viridis"),
            (expr[:, 1], "CD8A log1p CP10k", "magma"),
        ]
        for ax, (val, title, cmap) in zip(axes[:2], panels):
            sca = ax.scatter(um[:, 0], um[:, 1], c=val, s=2, cmap=cmap, linewidths=0)
            ax.set_title(title, fontsize=10)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
        colors = np.array(["#d9d9d9", "#d73027", "#4575b4", "#7b3294"])
        axes[2].scatter(um[:, 0], um[:, 1], c=colors[cat], s=2, linewidths=0)
        axes[2].set_title("Residual Gi*  red hot CLDN4 / blue cold CD8 / purple both", fontsize=8)
        axes[2].set_aspect("equal")
        axes[2].set_xticks([])
        axes[2].set_yticks([])
        fig.suptitle(
            f"{rec['sample']} FOV {int(rec['fov'])} ({tag}), partial I={rec['partial_moran_I']:+.3f}",
            fontsize=11,
        )
        fig.tight_layout()
        stem = fig_dir / f"gi_overlay_{tag}"
        fig.savefig(stem.with_suffix(".png"), dpi=150)
        fig.savefig(stem.with_suffix(".pdf"))
        plt.close(fig)


def collect(jobs: list[dict], workers: int) -> pd.DataFrame:
    import multiprocessing as mp

    if not jobs:
        return pd.DataFrame()
    t0 = time.time()
    ctx = mp.get_context("fork")
    rows = []
    with ctx.Pool(workers) as pool:
        done = 0
        for batch in pool.imap_unordered(work, jobs, chunksize=1):
            rows.extend(batch)
            done += 1
            if done % 10 == 0 or done == len(jobs):
                print(f"[progress] {done}/{len(jobs)} units in {time.time() - t0:.0f}s", flush=True)
    return pd.DataFrame(rows)


def build_headline(fov: pd.DataFrame, arr: dict) -> dict:
    rolls = {}
    prim = fov[fov["graph"] == "radius_50um"]
    for value in [
        "partial_moran_I",
        "moran_I",
        "partial_lee_L",
        "tumor_partial_r",
        "krt8_moran_I",
        "cldn_only_moran_I",
        "tumor_q_diff",
        "tumor_q_diff_partial",
    ]:
        rolls[value] = section_patient_rollups(prim, value)
    prose = conclusion_text(rolls.get("partial_moran_I"), rolls.get("moran_I"), rolls.get("tumor_partial_r"))
    return {
        "n_cells": int(arr["xy"].shape[0]),
        "n_fov": int(fov["unit_id"].nunique()) if len(fov) else 0,
        "n_samples": int(fov["sample"].nunique()) if len(fov) else 0,
        "n_patients": int(fov["patient"].nunique()) if len(fov) else 0,
        "px_to_um": PX_TO_UM,
        "n_perm": N_PERM,
        "seed": SEED,
        "primary_graph": "radius_50um",
        "coprimary_graph": "knn_k15",
        "rolls": rolls,
        "prose": prose,
    }


def self_test() -> None:
    rng = np.random.default_rng(10)
    xy = rng.normal(size=(90, 2))
    x = rng.normal(size=90)
    y = np.abs(rng.normal(size=90)) + 0.2
    from esda.getisord import G_Local
    from esda.moran import Moran_BV
    from libpysal.weights import KNN

    w = KNN.from_array(xy, k=6)
    w.transform = "r"
    mb = Moran_BV(x, y, w, transformation="r", permutations=0)
    W = knn_binary(xy, 6)
    prepared = row_standardize_kept(W)
    assert prepared is not None
    Wr, keep, _ = prepared
    assert keep.all()
    zx = zscore(x, 1)
    zy = zscore(y, 1)
    I = float(zx @ (Wr @ zy)) / (len(x) - 1)
    if abs(I - float(mb.I)) > 1e-8:
        raise SystemExit(f"Moran mismatch {I} vs {mb.I}")
    # Lee matches the explicit quadratic form when rows sum to 1.
    zx0 = zscore(x, 0)
    zy0 = zscore(y, 0)
    V = Wr.tocsr()
    ctc = V.T @ V
    ones = np.ones(len(x))
    den = float(ones @ (ctc @ ones))
    L_ref = float(zx0 @ (ctc @ zy0) / den)
    L_fast = float((V @ zx0) @ (V @ zy0) / len(x))
    if abs(L_ref - L_fast) > 1e-8:
        raise SystemExit(f"Lee mismatch {L_fast} vs {L_ref}")
    wb = KNN.from_array(xy, k=6)
    g = G_Local(y, wb, transform="B", permutations=0, star=True)
    gz = gi_star_z(y, binary_star(W))
    if np.nanmax(np.abs(gz - g.Zs)) > 1e-6:
        raise SystemExit(f"Gi* mismatch max {np.nanmax(np.abs(gz - g.Zs))}")
    # Residual Gi* is finite and translation-stable in the numerator sense.
    resid = y - y.mean()
    gz_r = gi_star_z(resid, binary_star(W))
    if not np.isfinite(gz_r).all():
        raise SystemExit("residual Gi* produced non-finite z")
    y2, r2 = residualize(2 * x + 0.01 * rng.normal(size=90), x)
    if r2 < 0.99:
        raise SystemExit(f"residualize R2 {r2}")
    if abs(np.dot(y2 - y2.mean(), x - x.mean())) > 1e-6:
        raise SystemExit("residual not orthogonal")
    print("[self-test] Moran, Lee, Gi*, residualize OK", flush=True)


def find_h5ad(arg: str | None) -> Path:
    cands = []
    if arg:
        cands.append(Path(arg))
    if os.environ.get("COSMX_H5AD"):
        cands.append(Path(os.environ["COSMX_H5AD"]))
    cands.append(REPO / "data" / "cosmx_human_nsclc_clustered.h5ad")
    cands.append(Path("/tmp/cosmx/cosmx_human_nsclc_clustered.h5ad"))
    for c in cands:
        if c.is_file():
            return c
    raise SystemExit("h5ad not found. Run download.sh or pass --h5ad.")


def main() -> None:
    global XY, LOGX, TUMOR, N_PERM, SEED
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=None)
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--perms", type=int, default=199)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--sample-workers", type=int, default=2)
    ap.add_argument("--max-fovs", type=int, default=0)
    ap.add_argument("--skip-sample", action="store_true")
    ap.add_argument("--seed", type=int, default=20260921)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return
    N_PERM = int(args.perms)
    SEED = int(args.seed)
    out = Path(args.out)
    tab = out / "tables"
    fig = out / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    h5ad = find_h5ad(args.h5ad)
    cache = Path("/tmp/cosmx/expr_cache_v1.npz")
    if cache.is_file():
        arr = load_cache(cache)
        print(f"[cache] hit {cache}", flush=True)
    else:
        arr = extract(h5ad, cache)
    XY = arr["xy"]
    LOGX = arr["logx"]
    TUMOR = arr["tumor"]
    fov_jobs, sample_jobs = build_jobs(arr, args.max_fovs, args.skip_sample)
    print(f"[jobs] FOVs={len(fov_jobs)} sections={len(sample_jobs)} perms={N_PERM}", flush=True)
    fov = collect(fov_jobs, args.workers)
    fov = add_fdr(fov)
    fov.to_csv(tab / "fov_graph_stats.tsv", sep="\t", index=False)
    sample = collect(sample_jobs, args.sample_workers) if sample_jobs else pd.DataFrame()
    if len(sample):
        sample = add_fdr(sample)
        sample.to_csv(tab / "sample_graph_stats.tsv", sep="\t", index=False)
    headline = build_headline(fov, arr)
    # JSON cannot hold raw numpy types nested in a messy way; default=str via a cleaner dump.
    (tab / "headline.json").write_text(json.dumps(headline, indent=2, default=_json_default))
    write_results(fov, sample, headline, out / "RESULTS.md")
    make_plots(fov, sample, arr, fig)
    print("[done]", out / "RESULTS.md", flush=True)


def _json_default(obj):
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    raise TypeError(type(obj))


if __name__ == "__main__":
    main()
