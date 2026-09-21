#!/usr/bin/env python3
"""CosMx NSCLC: Hotspot and Potts-HMRF gene modules for CLDN4 vs CD8 fields.

Official He et al. 2022 CosMx SMI 960-plex (8 sections, 5 patients).
Flat files from the NanoString public S3 release. This does not re-estimate
the locked 20 µm contact odds ratio or the 50/100 µm cytotoxic ratios.

Pre-specified before looking at results
----------------------------------------
Graph. Within each FOV, Hotspot uses a 30-NN Gaussian spatial graph
(neighborhood_factor=3, row-normalized, then Hotspot's non-redundant
weights). HMRF uses a symmetrized 8-NN Potts graph.

Hotspot (DeTomaso & Yosef, Cell Systems 2021), reimplemented:
  DANB depth model, centered local autocorrelation, BH FDR < 0.05,
  conditional pairwise local-correlation Z, average-linkage modules.
  Program cut: min_gene_threshold=20, core_only=True (package defaults).
  Fine cut: min_gene_threshold=4, core_only=False, same pairwise FDR 0.05.
  The fine cut is there because a junction-sized module cannot clear a
  20-gene minimum on a 960-plex that contains only four epithelial
  adhesion genes.

HMRF: 2-state Potts model per gene (ICM, beta=0.5, 20 synchronous
sweeps). High-state probability fields are clustered by average linkage
on correlation distance, cut at Pearson correlation 0.40.
A cut of 0.50 is above the 99th percentile of pairwise field
correlations on this panel (sparse high-states) and forms no modules;
0.40 was set from that scale check on one FOV before the cohort run.
  Fine modules: >= 4 genes. Program modules: >= 15 genes.
  Input genes: FDR-significant Hotspot genes, capped at the top 150 by
  autocorrelation Z, plus the watchlist below.

Adhesion-like module (on-panel; not a classical TJ call): contains CLDN4
and at least 2 of {CDH1, EPCAM, TACSTD2}. Keratin-mixed if it also
contains at least 2 of {KRT7, KRT8, KRT18, KRT19}. Classical TJ genes
other than CLDN4 are recorded as absent; they are not imputed.

Spatial association, among RNA-epithelial cells (marker-score argmax,
same compartment rule as the contact analysis; CLDN4 is not a marker):
  CD8 field = count of CD8 T cells (CD8A/B > 0, CD3D/E/G > 0, not
  epithelial) inside 50 µm. Pixel size 0.18 µm, as in the contact analysis.
  Primary endpoint: Spearman rho of CLDN4 (library-size log1p) vs that count.
  Module endpoint: same, using the Hotspot fine-module score where CLDN4
  is assigned. Score direction is fixed so it increases with CLDN4.
  Patient-level statistic: unweighted mean of five patient means
  (each patient mean is the unweighted mean of its FOVs).
  Null: 199 toroidal shifts of the CD8 point pattern inside the FOV box.
  One-sided p for anti-correlation (rho < 0).

Specificity, same geometry: keratin score, CLDN4 residualized on keratin,
adhesion score, adhesion without CLDN4, Hotspot program score, HMRF
scores, and CLDN4 at 100 µm.
Same-cell CLDN4–CD8 transcript correlation is a cell-type control, not
an exclusion test. Domain HMRF (5 domains, 10 PCs, beta=0.5) is a
composition description, not a permutation test.

FOV rules: modules if >= 80 QC cells (nCount>=20, nGene>=5 on the 960
panel, cell_ID!=0). Spatial tests also require >= 40 epithelial cells
and >= 5 CD8 T cells.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import zlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial import cKDTree
from scipy.spatial.distance import squareform
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "cosmx"
CACHE = ROOT / "data" / "cosmx_cache"
OUT = ROOT / "results" / "cosmx_hotspot_hmrf"

SAMPLES = [
    "Lung5_Rep1",
    "Lung5_Rep2",
    "Lung5_Rep3",
    "Lung6",
    "Lung9_Rep1",
    "Lung9_Rep2",
    "Lung12",
    "Lung13",
]
PATIENT = {
    "Lung5_Rep1": "Lung5",
    "Lung5_Rep2": "Lung5",
    "Lung5_Rep3": "Lung5",
    "Lung6": "Lung6",
    "Lung9_Rep1": "Lung9",
    "Lung9_Rep2": "Lung9",
    "Lung12": "Lung12",
    "Lung13": "Lung13",
}
PATIENT_ORDER = ["Lung5", "Lung6", "Lung9", "Lung12", "Lung13"]

PX_TO_UM = 0.18
RADIUS_UM = (50.0, 100.0)
MIN_COUNTS = 20
MIN_GENES = 5
MIN_FOV_QC = 80
MIN_EPI = 40
MIN_CD8 = 5
HOTSPOT_K = 30
HOTSPOT_FACTOR = 3
HMRF_K = 8
HMRF_BETA = 0.5
HMRF_ITERS = 20
# Pearson 0.5 is above the 99th percentile of these sparse HMRF fields
# (checked on one FOV before the cohort run) and returns no modules.
# 0.40 is that FOV's upper decile and yields modules of about 4–12 genes.
HMRF_CORR = 0.40
HMRF_TOP = 150
GENE_FDR = 0.05
PAIR_FDR = 0.05
N_PERM_DEFAULT = 199
SEED = 20260921

EPI_MARKERS = [
    "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT17", "KRT15",
    "KRT14", "KRT13", "KRT16", "KRT6A", "KRT5", "KRT7", "CEACAM6",
    "MUC1", "ELF3", "SFTPC", "SFTPB", "NAPSA", "S100A14",
]
IMM_MARKERS = [
    "PTPRC", "CD3D", "CD3E", "CD3G", "CD2", "CD8A", "CD8B", "CD4",
    "IL7R", "CCL5", "GZMA", "GZMB", "GZMK", "GNLY", "KLRB1", "KLRK1",
    "CD68", "CD163", "C1QA", "C1QB", "C1QC", "LYZ", "FCGR3A", "ITGAX",
    "ITGAM", "CD14", "CD79A", "MS4A1", "CD19", "CD37", "CD52", "CD53",
]
STR_MARKERS = [
    "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "COL6A2", "DCN",
    "LUM", "FN1", "BGN", "PDGFRB", "ACTA2", "MYH11", "PECAM1", "VWF",
    "CDH5", "CLEC14A", "FLT1", "KDR", "ESAM", "RAMP2",
]
ADHESION = ["CLDN4", "CDH1", "EPCAM", "TACSTD2"]
ADHESION_WO = ["CDH1", "EPCAM", "TACSTD2"]
KERATIN = ["KRT7", "KRT8", "KRT18", "KRT19"]
ENDO_JUNC = ["ESAM", "CDH5"]
CLASSICAL_TJ = [
    "CLDN1", "CLDN3", "CLDN4", "CLDN5", "CLDN7", "CLDN18",
    "OCLN", "TJP1", "TJP2", "TJP3", "F11R", "CGN",
    "MARVELD2", "MARVELD3", "CRB3", "PARD3", "JAM2", "JAM3",
]
WATCHLIST = ADHESION + KERATIN + ENDO_JUNC + [
    "MUC1", "CEACAM6", "ELF3", "KRT5", "KRT17",
    "CD8A", "CD8B", "CD3D", "CD3E", "PTPRC", "CD68", "GZMB", "PRF1", "NKG7",
    "COL1A1", "DCN", "PECAM1", "VWF",
]
HEATMAP = [
    "CLDN4", "CDH1", "EPCAM", "TACSTD2",
    "KRT8", "KRT18", "KRT19", "KRT7",
    "ESAM", "CDH5", "CD8A", "CD8B", "PTPRC", "COL1A1",
]
PERM_ENDPOINTS = [
    "cldn4_50",
    "hs_fine_50",
    "keratin_50",
    "resid_50",
    "adhesion_50",
    "adhesion_wo_50",
    "hs_program_50",
    "hmrf_fine_50",
    "hmrf_program_50",
    "cldn4_100",
]


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=np.float64)
    n = p.size
    if n == 0:
        return p.copy()
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out = np.empty(n, dtype=np.float64)
    out[order] = q
    return out


def find_col(columns, options):
    lower = {str(c).replace("\ufeff", "").strip().lower(): c for c in columns}
    for opt in options:
        if opt.lower() in lower:
            return lower[opt.lower()]
    return None


def make_weights_non_redundant(neighbors: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Hotspot knn.make_weights_non_redundant, including the j < i skip."""
    w = np.array(weights, dtype=np.float64, copy=True)
    n, k = neighbors.shape
    inv = []
    for i in range(n):
        inv.append({int(neighbors[i, t]): t for t in range(k)})
    for i in range(n):
        for t in range(k):
            j = int(neighbors[i, t])
            if j < i:
                continue
            slot = inv[j].get(i)
            if slot is not None:
                w[i, t] += w[j, slot]
                w[j, slot] = 0.0
    return w


def hotspot_graph(xy: np.ndarray, k: int = HOTSPOT_K, neighborhood_factor: int = HOTSPOT_FACTOR):
    n = xy.shape[0]
    kk = min(k, n - 1)
    if kk < 1:
        raise ValueError("need at least 2 cells")
    nn = NearestNeighbors(n_neighbors=kk + 1, algorithm="kd_tree")
    nn.fit(xy)
    dist, ind = nn.kneighbors(xy)
    dist = np.ascontiguousarray(dist[:, 1:], dtype=np.float64)
    ind = np.ascontiguousarray(ind[:, 1:], dtype=np.int32)
    radius_ii = int(np.ceil(kk / neighborhood_factor))
    sigma = dist[:, [radius_ii - 1]].copy()
    sigma[sigma == 0] = 1.0
    weights = np.exp(-(dist ** 2) / (sigma ** 2))
    wnorm = weights.sum(axis=1, keepdims=True)
    wnorm[wnorm == 0] = 1.0
    weights = weights / wnorm
    weights = make_weights_non_redundant(ind, weights)
    rows = np.repeat(np.arange(n, dtype=np.int32), kk)
    cols = ind.ravel()
    data = weights.ravel()
    mask = data != 0
    W = csr_matrix((data[mask], (rows[mask], cols[mask])), shape=(n, n))
    return W.tocsr(), ind, weights


def binary_knn_adj(xy: np.ndarray, k: int = HMRF_K) -> csr_matrix:
    n = xy.shape[0]
    kk = min(k, n - 1)
    nn = NearestNeighbors(n_neighbors=kk + 1, algorithm="kd_tree")
    nn.fit(xy)
    ind = nn.kneighbors(xy, return_distance=False)[:, 1:]
    rows = np.repeat(np.arange(n, dtype=np.int32), kk)
    cols = ind.ravel().astype(np.int32)
    A = csr_matrix((np.ones(rows.size, dtype=np.float64), (rows, cols)), shape=(n, n))
    A = A.maximum(A.T).tocsr()
    A.setdiag(0)
    A.eliminate_zeros()
    return A


def danb_center(counts: np.ndarray) -> np.ndarray:
    """Depth-adjusted negative binomial standardization. counts: cells x genes."""
    counts = np.asarray(counts, dtype=np.float64)
    n, g = counts.shape
    tis = counts.sum(axis=1)
    total = float(tis.sum())
    if total <= 0:
        return np.zeros_like(counts)
    tj = counts.sum(axis=0)
    mu = np.outer(tis, tj) / total
    my_rowvar = (counts - mu).var(axis=0, ddof=1)
    denom = (n - 1) * my_rowvar - tj
    numer = ((tj ** 2) / total) * ((np.dot(tis, tis)) / total)
    size = np.full(g, 1e9, dtype=np.float64)
    ok = np.isfinite(denom) & (denom > 0) & np.isfinite(numer)
    size[ok] = numer[ok] / denom[ok]
    size[~np.isfinite(size) | (size < 0)] = 1e9
    size = np.clip(size, 1e-10, 1e9)
    var = mu * (1.0 + mu / size[None, :])
    std = np.sqrt(np.maximum(var, 1e-12))
    return (counts - mu) / std


def autocorr_z(X: np.ndarray, W: csr_matrix) -> np.ndarray:
    WX = W @ X
    gstat = np.sum(X * WX, axis=0)
    wtot2 = float(W.power(2).sum())
    std = np.sqrt(wtot2) if wtot2 > 0 else 1.0
    return gstat / std, WX


def lc_z_against_one(X: np.ndarray, WX: np.ndarray, eg2: np.ndarray, g: int) -> np.ndarray:
    """Conditional local-correlation Z of every gene with column g."""
    lc = X.T @ WX[:, g] + WX.T @ X[:, g]
    std = np.sqrt(np.maximum(eg2, 1e-12))
    z_partner = lc / std
    z_self = lc / std[g]
    z = np.where(np.abs(z_partner) < np.abs(z_self), z_partner, z_self)
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    z[g] = np.nan
    return z


def pairwise_from_products(X: np.ndarray, WX: np.ndarray, eg2: np.ndarray) -> np.ndarray:
    lc = X.T @ WX
    lc = lc + lc.T
    np.fill_diagonal(lc, 0.0)
    std = np.sqrt(np.maximum(eg2, 1e-12))
    z_xy = lc / std[:, None]
    z_yx = lc / std[None, :]
    z = np.where(np.abs(z_xy) < np.abs(z_yx), z_xy, z_yx)
    np.fill_diagonal(z, 0.0)
    return np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)


def eg2_from_W(X: np.ndarray, W: csr_matrix) -> np.ndarray:
    Wsym = (W + W.T).tocsr()
    S = Wsym @ X
    return np.sum(S * S, axis=0), Wsym


def calc_mean_dists(Z, node_index, out_mean_dists):
    n = Z.shape[0] + 1
    left = int(Z[node_index, 0] - n)
    right = int(Z[node_index, 1] - n)
    if left < 0:
        left_avg, left_m = 0.0, 0.0
    else:
        left_avg, left_m = calc_mean_dists(Z, left, out_mean_dists)
    if right < 0:
        right_avg, right_m = 0.0, 0.0
    else:
        right_avg, right_m = calc_mean_dists(Z, right, out_mean_dists)
    height = float(Z[node_index, 2])
    merges = left_m + right_m + 1.0
    average = (left_avg * left_m + right_avg * right_m + height) / merges
    out_mean_dists[node_index] = average
    return average, merges


def prop_label(Z, node_index, label, labels, out_clusters):
    n = Z.shape[0] + 1
    if label == -1:
        label = labels[node_index]
    left = int(Z[node_index, 0] - n)
    right = int(Z[node_index, 1] - n)
    if left < 0:
        out_clusters[left + n] = label
    else:
        prop_label(Z, left, label, labels, out_clusters)
    if right < 0:
        out_clusters[right + n] = label
    else:
        prop_label(Z, right, label, labels, out_clusters)


def prop_label2(Z, node_index, label, labels, out_clusters):
    n = Z.shape[0] + 1
    parent_label = label
    this_label = labels[node_index]
    # Hotspot: a node label of -1 inherits the parent; otherwise it overrides.
    new_label = parent_label if this_label == -1 else this_label
    left = int(Z[node_index, 0] - n)
    right = int(Z[node_index, 1] - n)
    if left < 0:
        out_clusters[left + n] = new_label
    else:
        prop_label2(Z, left, new_label, labels, out_clusters)
    if right < 0:
        out_clusters[right + n] = new_label
    else:
        prop_label2(Z, right, new_label, labels, out_clusters)


def _remap_modules(raw: np.ndarray) -> np.ndarray:
    out = np.asarray(raw, dtype=np.int32)
    uniq = sorted({int(v) for v in out if int(v) != -1})
    mp = {u: i + 1 for i, u in enumerate(uniq)}
    mp[-1] = -1
    return np.array([mp.get(int(v), -1) for v in out], dtype=np.int32)


def assign_modules(Zlink, offset, min_threshold, z_threshold):
    n_internal = Zlink.shape[0]
    n = n_internal + 1
    labels = np.full(n_internal, -1, dtype=np.int32)
    mean_dists = np.zeros(n_internal, dtype=np.float64)
    calc_mean_dists(Zlink, n_internal - 1, mean_dists)
    clust_i = 0
    for i in range(n_internal):
        ca = int(Zlink[i, 0])
        cb = int(Zlink[i, 1])
        if ca - n < 0:
            n_a, clust_a = 1, -1
        else:
            n_a, clust_a = int(Zlink[ca - n, 3]), int(labels[ca - n])
        if cb - n < 0:
            n_b, clust_b = 1, -1
        else:
            n_b, clust_b = int(Zlink[cb - n, 3]), int(labels[cb - n])
        if Zlink[i, 2] > offset - z_threshold:
            new = -1
        elif n_a >= min_threshold and n_b >= min_threshold:
            dist_a = mean_dists[ca - n]
            dist_b = mean_dists[cb - n]
            new = clust_a if dist_a >= dist_b else clust_b
        elif n_a >= min_threshold:
            new = clust_a
        elif n_b >= min_threshold:
            new = clust_b
        elif (n_a + n_b) >= min_threshold:
            new = clust_i
            clust_i += 1
        else:
            new = -1
        labels[i] = new
    out = np.full(n, -2, dtype=np.int32)
    prop_label2(Zlink, n_internal - 1, int(labels[-1]), labels, out)
    return _remap_modules(out)


def assign_modules_core(Zlink, offset, min_threshold, z_threshold):
    n_internal = Zlink.shape[0]
    n = n_internal + 1
    labels = np.full(n_internal, -1, dtype=np.int32)
    clust_i = 0
    for i in range(n_internal):
        ca = int(Zlink[i, 0])
        cb = int(Zlink[i, 1])
        if ca - n < 0:
            n_a, clust_a = 1, -1
        else:
            n_a, clust_a = int(Zlink[ca - n, 3]), int(labels[ca - n])
        if cb - n < 0:
            n_b, clust_b = 1, -1
        else:
            n_b, clust_b = int(Zlink[cb - n, 3]), int(labels[cb - n])
        if n_a >= min_threshold and n_b >= min_threshold:
            new = -1
        elif Zlink[i, 2] > offset - z_threshold:
            new = -1
        elif n_a >= min_threshold:
            new = clust_a
        elif n_b >= min_threshold:
            new = clust_b
        elif (n_a + n_b) >= min_threshold:
            new = clust_i
            clust_i += 1
        else:
            new = -1
        labels[i] = new
    out = np.full(n, -2, dtype=np.int32)
    prop_label(Zlink, n_internal - 1, int(labels[-1]), labels, out)
    return _remap_modules(out)


def hotspot_modules(Z: np.ndarray, min_gene_threshold: int, fdr_threshold: float, core_only: bool):
    Z = np.asarray(Z, dtype=np.float64)
    Z = (Z + Z.T) / 2.0
    np.fill_diagonal(Z, 0.0)
    Z = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)
    n = Z.shape[0]
    if n < 2:
        return np.full(n, -1, dtype=np.int32)
    upper = squareform(Z, checks=False)
    all_z = np.sort(upper)
    all_q = bh_fdr(stats.norm.sf(all_z))
    ii = np.nonzero(all_q < fdr_threshold)[0]
    z_threshold = float(all_z[ii[0]]) if ii.size else float(all_z[-1] + 1.0)
    condensed = -upper
    offset = float(-condensed.min())
    condensed = condensed + offset
    floor = float(condensed.min())
    if floor < 0:
        condensed = condensed - floor
        offset = offset - floor
    link = linkage(condensed, method="average")
    if core_only:
        return assign_modules_core(link, offset, min_gene_threshold, z_threshold)
    return assign_modules(link, offset, min_gene_threshold, z_threshold)


def smooth_cols(X: np.ndarray, W: csr_matrix, lam: float = 0.9) -> np.ndarray:
    Wsym = (W + W.T).tocsr()
    deg = np.asarray(Wsym.sum(axis=1)).ravel()
    deg[deg == 0] = 1.0
    sm = (Wsym @ X) / deg[:, None]
    return lam * sm + (1.0 - lam) * X


def align_to_reference(score: np.ndarray, ref: np.ndarray) -> tuple[np.ndarray, bool]:
    s = score - np.mean(score)
    r = ref - np.mean(ref)
    agrees = True
    if np.dot(s, r) < 0:
        score = -score
        agrees = False
    return score, agrees


def hotspot_module_score(X: np.ndarray, W: csr_matrix, members: np.ndarray, ref_col: int):
    sub = smooth_cols(X[:, members], W)
    if sub.shape[1] == 1:
        score = sub[:, 0].copy()
    else:
        pca = PCA(n_components=1, random_state=0)
        score = pca.fit_transform(sub)[:, 0]
        loadings = pca.components_[0].copy()
        if float(loadings.mean()) < 0:
            score = -score
    score, agrees = align_to_reference(score, X[:, ref_col])
    return score, agrees


def hmrf_high_prob(x: np.ndarray, A: csr_matrix, beta: float = HMRF_BETA, n_iter: int = HMRF_ITERS) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    if np.unique(x).size < 2:
        return np.full(n, 0.5)
    mu = np.array([np.quantile(x, 0.25), np.quantile(x, 0.75)], dtype=np.float64)
    var = np.array([float(np.var(x)) + 1e-3, float(np.var(x)) + 1e-3])
    lab = (x >= np.median(x)).astype(np.float64)
    deg = np.asarray(A.sum(axis=1)).ravel()
    for _ in range(n_iter):
        for s in (0, 1):
            m = lab == s
            if int(m.sum()) >= 5:
                mu[s] = float(x[m].mean())
                var[s] = float(x[m].var()) + 1e-3
        if mu[1] < mu[0]:
            mu = mu[::-1].copy()
            var = var[::-1].copy()
            lab = 1.0 - lab
        v1 = A @ lab
        v0 = deg - v1
        ll1 = -0.5 * np.log(var[1]) - 0.5 * ((x - mu[1]) ** 2) / var[1]
        ll0 = -0.5 * np.log(var[0]) - 0.5 * ((x - mu[0]) ** 2) / var[0]
        lab = (ll1 + beta * v1 >= ll0 + beta * v0).astype(np.float64)
    for s in (0, 1):
        m = lab == s
        if int(m.sum()) >= 5:
            mu[s] = float(x[m].mean())
            var[s] = float(x[m].var()) + 1e-3
    if mu[1] < mu[0]:
        mu = mu[::-1].copy()
        var = var[::-1].copy()
        lab = 1.0 - lab
        v1 = A @ lab
    else:
        v1 = A @ lab
    v0 = deg - v1
    ll1 = -0.5 * np.log(var[1]) - 0.5 * ((x - mu[1]) ** 2) / var[1]
    ll0 = -0.5 * np.log(var[0]) - 0.5 * ((x - mu[0]) ** 2) / var[0]
    logit = np.clip((ll1 - ll0) + beta * (v1 - v0), -30, 30)
    return 1.0 / (1.0 + np.exp(-logit))


def correlation_modules(P: np.ndarray, min_size: int, min_corr: float) -> np.ndarray:
    """P is cells x genes. Returns module ids, -1 if below min_size."""
    n_genes = P.shape[1]
    labels = np.full(n_genes, -1, dtype=np.int32)
    if n_genes < 2:
        return labels
    sd = P.std(axis=0)
    keep = sd > 1e-6
    if int(keep.sum()) < 2:
        return labels
    idx = np.nonzero(keep)[0]
    C = np.corrcoef(P[:, idx], rowvar=False)
    C = np.nan_to_num(C, nan=0.0)
    np.fill_diagonal(C, 1.0)
    D = np.clip(1.0 - C, 0.0, 2.0)
    D = (D + D.T) / 2.0
    np.fill_diagonal(D, 0.0)
    condensed = squareform(D, checks=False)
    link = linkage(condensed, method="average")
    raw = fcluster(link, t=1.0 - min_corr, criterion="distance")
    # Drop modules smaller than min_size, then renumber.
    out = np.full(idx.size, -1, dtype=np.int32)
    next_id = 1
    for u in np.unique(raw):
        members = np.nonzero(raw == u)[0]
        if members.size >= min_size:
            out[members] = next_id
            next_id += 1
    labels[idx] = out
    return labels


def library_lognorm(counts: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=1).astype(np.float64)
    pos = totals[totals > 0]
    med = float(np.median(pos)) if pos.size else 1.0
    scaled = counts.astype(np.float64) * (med / np.maximum(totals, 1.0))[:, None]
    return np.log1p(scaled)


def mean_genes(X: np.ndarray, gmap: dict, names: list[str]) -> np.ndarray | None:
    idx = [gmap[g] for g in names if g in gmap]
    if not idx:
        return None
    return X[:, idx].mean(axis=1)


def spearman(a, b) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size < 10 or np.unique(a).size < 2 or np.unique(b).size < 2:
        return float("nan")
    rho = stats.spearmanr(a, b).statistic
    return float(rho)


def rank_resid(y, x):
    y = np.asarray(y, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    if np.unique(y).size < 2 or np.unique(x).size < 2:
        return None
    ry = stats.rankdata(y).astype(np.float64)
    rx = stats.rankdata(x).astype(np.float64)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = float(np.dot(rx, rx))
    if denom <= 0:
        return None
    return ry - (float(np.dot(rx, ry)) / denom) * rx


def ball_count(tree: cKDTree, query: np.ndarray, radius: float) -> np.ndarray:
    try:
        out = tree.query_ball_point(query, r=radius, return_length=True, workers=-1)
    except TypeError:
        out = tree.query_ball_point(query, r=radius, return_length=True)
    return np.asarray(out, dtype=np.float64)


def torus_shift(xy: np.ndarray, lo: np.ndarray, span: np.ndarray, shift: np.ndarray) -> np.ndarray:
    return lo + np.mod(xy - lo + shift, span)


def classify_members(member_names: set[str]) -> str:
    if "CLDN4" not in member_names:
        return "unassigned"
    n_ad = sum(g in member_names for g in ADHESION_WO)
    n_ke = sum(g in member_names for g in KERATIN)
    if n_ad >= 2 and n_ke >= 2:
        return "adhesion_and_keratin"
    if n_ad >= 2:
        return "adhesion_like"
    if n_ke >= 2:
        return "keratin_without_adhesion"
    return "other"


def module_name_set(mod_row: np.ndarray, genes: list[str], cldn4_code: int) -> set[str]:
    if cldn4_code <= 0:
        return set()
    return {genes[i] for i, v in enumerate(mod_row) if int(v) == cldn4_code}


def domain_hmrf(logx: np.ndarray, xy: np.ndarray):
    n = logx.shape[0]
    if n < 150:
        return None
    sd = logx.std(axis=0)
    keep = sd > 0
    if int(keep.sum()) < 5:
        return None
    Y = logx[:, keep]
    Y = (Y - Y.mean(axis=0)) / Y.std(axis=0)
    n_comp = int(min(10, n - 1, Y.shape[1]))
    pcs = PCA(n_components=n_comp, random_state=0).fit_transform(Y)
    lab = KMeans(n_clusters=5, n_init=10, random_state=0).fit_predict(pcs)
    A = binary_knn_adj(xy, k=HMRF_K)
    for _ in range(HMRF_ITERS):
        mu = np.zeros((5, n_comp))
        var = np.ones((5, n_comp))
        for s in range(5):
            m = lab == s
            if int(m.sum()) >= 5:
                mu[s] = pcs[m].mean(axis=0)
                var[s] = pcs[m].var(axis=0) + 1e-3
            else:
                mu[s] = pcs.mean(axis=0)
                var[s] = pcs.var(axis=0) + 1e-3
        ll = np.zeros((n, 5))
        for s in range(5):
            ll[:, s] = -0.5 * np.log(var[s]).sum() - 0.5 * np.sum((pcs - mu[s]) ** 2 / var[s], axis=1)
        votes = np.column_stack([A @ (lab == s).astype(np.float64) for s in range(5)])
        lab = np.argmax(ll + HMRF_BETA * votes, axis=1)
    return lab.astype(np.int32)


def naive_autocorr(X, neighbors, weights):
    n, k = neighbors.shape
    gstat = np.zeros(X.shape[1])
    for g in range(X.shape[1]):
        x = X[:, g]
        acc = 0.0
        for i in range(n):
            xi = x[i]
            for t in range(k):
                w = weights[i, t]
                if w == 0:
                    continue
                acc += xi * x[neighbors[i, t]] * w
        gstat[g] = acc
    return gstat / np.sqrt(np.sum(weights ** 2))


def naive_pair_lc(x, y, neighbors, weights):
    n, k = neighbors.shape
    acc = 0.0
    for i in range(n):
        for t in range(k):
            w = weights[i, t]
            if w == 0:
                continue
            j = neighbors[i, t]
            acc += w * (x[i] * y[j] + y[i] * x[j])
    return acc


def naive_eg2(x, neighbors, weights):
    n, k = neighbors.shape
    t1 = np.zeros(n)
    for i in range(n):
        for t in range(k):
            w = weights[i, t]
            if w == 0:
                continue
            j = int(neighbors[i, t])
            t1[i] += w * x[j]
            t1[j] += w * x[i]
    return float(np.sum(t1 ** 2))


def self_test() -> None:
    sys.setrecursionlimit(10000)
    rng = np.random.default_rng(0)
    xy = rng.normal(size=(40, 2))
    W, ind, weights = hotspot_graph(xy, k=5, neighborhood_factor=3)
    counts = rng.poisson(1.5, size=(40, 3)).astype(np.float64)
    X = danb_center(counts)
    z_vec, WX = autocorr_z(X, W)
    z_naive = naive_autocorr(X, ind, weights)
    if not np.allclose(z_vec, z_naive, atol=1e-6):
        raise AssertionError(f"autocorr mismatch {z_vec} vs {z_naive}")
    eg2, _ = eg2_from_W(X, W)
    Z = pairwise_from_products(X, WX, eg2)
    lc = naive_pair_lc(X[:, 0], X[:, 1], ind, weights)
    std0 = np.sqrt(max(naive_eg2(X[:, 0], ind, weights), 1e-12))
    std1 = np.sqrt(max(naive_eg2(X[:, 1], ind, weights), 1e-12))
    z01 = lc / std0
    z10 = lc / std1
    expect = z01 if abs(z01) < abs(z10) else z10
    if abs(Z[0, 1] - expect) > 1e-5:
        raise AssertionError(f"pair Z mismatch {Z[0, 1]} vs {expect}")

    # Spatial structure should beat a permutation of the same values.
    # Extra genes keep the DANB depth factor from forcing a two-gene anti-correlation.
    grid = np.stack(np.meshgrid(np.arange(20), np.arange(20), indexing="ij"), axis=-1).reshape(-1, 2).astype(float)
    signal = (grid[:, 0] > 10).astype(np.float64) * 30.0
    noise = rng.poisson(0.3, size=signal.shape[0]).astype(np.float64)
    bg_mat = rng.poisson(1.0, size=(signal.size, 20)).astype(np.float64)
    counts_s = np.column_stack([
        signal + noise,
        rng.permutation(signal) + noise,
        bg_mat,
    ])
    Xs = danb_center(counts_s)
    Ws, _, _ = hotspot_graph(grid, k=8)
    zs, _ = autocorr_z(Xs, Ws)
    if not zs[0] > zs[1] + 5:
        raise AssertionError(f"spatial autocorr did not separate signal {zs}")

    # Two tight blocks should not be merged by the core cutter.
    n_block = 25
    n_g = n_block * 2
    Zb = np.zeros((n_g, n_g))
    Zb[:n_block, :n_block] = 25
    Zb[n_block:, n_block:] = 25
    np.fill_diagonal(Zb, 0)
    mod = hotspot_modules(Zb, min_gene_threshold=20, fdr_threshold=0.05, core_only=True)
    if np.any(mod < 0):
        raise AssertionError(f"core modules left genes unassigned: {np.unique(mod, return_counts=True)}")
    if mod[0] == mod[n_block] or len(set(mod[:n_block])) != 1 or len(set(mod[n_block:])) != 1:
        raise AssertionError(f"block modules wrong: {mod[:3]} ... {mod[n_block:n_block+3]}")
    fine = hotspot_modules(Zb[:8, :8], min_gene_threshold=4, fdr_threshold=0.05, core_only=False)
    if len(set(fine.tolist())) != 1 or fine[0] <= 0:
        raise AssertionError(f"fine module cut failed: {fine}")

    # HMRF: left-high genes agree; a right-high gene does not.
    left = (grid[:, 0] < 10).astype(np.float64)
    right = 1.0 - left
    A = binary_knn_adj(grid, k=8)
    p_left = hmrf_high_prob(left * 5 + rng.normal(0, 0.05, left.size), A)
    p_left2 = hmrf_high_prob(left * 5 + rng.normal(0, 0.05, left.size), A)
    p_right = hmrf_high_prob(right * 5 + rng.normal(0, 0.05, right.size), A)
    c_same = float(np.corrcoef(p_left, p_left2)[0, 1])
    c_opp = float(np.corrcoef(p_left, p_right)[0, 1])
    if not (c_same > 0.8 and c_opp < 0):
        raise AssertionError(f"HMRF fields same={c_same} opp={c_opp}")
    print("self-test ok", flush=True)


def read_header(path: Path) -> list[str]:
    import csv
    with path.open(newline="") as handle:
        row = next(csv.reader(handle))
    return [c.replace("\ufeff", "").strip() for c in row]


def find_file(sample: str, kind: str) -> Path:
    hits = sorted(DATA.rglob(f"{sample}_{kind}_file.csv"))
    if not hits:
        raise FileNotFoundError(f"missing {sample} {kind} under {DATA}")
    return hits[0]


def panel_genes(header: list[str], fov_name: str, cell_name: str) -> list[str]:
    skip = {fov_name, cell_name}
    genes = []
    for col in header:
        if col in skip:
            continue
        if col.lower().startswith("neg"):
            continue
        genes.append(col)
    return genes


def load_sample(sample: str) -> dict:
    CACHE.mkdir(parents=True, exist_ok=True)
    cache = CACHE / f"{sample}.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=False)
        genes = z["genes"].astype(str).tolist()
        return {
            "genes": genes,
            "gmap": {g: i for i, g in enumerate(genes)},
            "fov": z["fov"].astype(np.int32),
            "cell_id": z["cell_id"].astype(np.int32),
            "xy": z["xy"].astype(np.float64),
            "counts": z["counts"].astype(np.float32),
        }
    t0 = time.time()
    expr_path = find_file(sample, "exprMat")
    meta_path = find_file(sample, "metadata")
    header = read_header(expr_path)
    fov_c = find_col(header, ["fov"])
    cell_c = find_col(header, ["cell_ID", "cell_id", "cellid"])
    if fov_c is None or cell_c is None:
        raise RuntimeError(f"{sample} expr header missing fov/cell_ID: {header[:12]}")
    genes = panel_genes(header, fov_c, cell_c)
    if "CLDN4" not in genes:
        raise RuntimeError(f"{sample} panel has no CLDN4")
    usecols = [fov_c, cell_c] + genes
    dtype = {fov_c: np.int32, cell_c: np.int32}
    dtype.update({g: np.float32 for g in genes})
    chunks = []
    for chunk in pd.read_csv(expr_path, usecols=usecols, dtype=dtype, chunksize=25000):
        chunk = chunk.rename(columns={fov_c: "fov", cell_c: "cell_ID"})
        chunk = chunk.loc[chunk["cell_ID"] != 0]
        chunks.append(chunk)
    expr = pd.concat(chunks, ignore_index=True)
    meta_head = read_header(meta_path)
    mf = find_col(meta_head, ["fov"])
    mc = find_col(meta_head, ["cell_ID", "cell_id", "cellid"])
    mx = find_col(meta_head, ["CenterX_local_px", "CenterX_local", "x_local_px", "CenterX_global_px"])
    my = find_col(meta_head, ["CenterY_local_px", "CenterY_local", "y_local_px", "CenterY_global_px"])
    if None in (mf, mc, mx, my):
        raise RuntimeError(f"{sample} metadata columns {meta_head[:40]}")
    meta = pd.read_csv(meta_path, usecols=[mf, mc, mx, my])
    meta = meta.rename(columns={mf: "fov", mc: "cell_ID", mx: "x", my: "y"})
    meta["fov"] = meta["fov"].astype(np.int32)
    meta["cell_ID"] = meta["cell_ID"].astype(np.int32)
    merged = expr.merge(meta, on=["fov", "cell_ID"], how="inner")
    merged = merged.drop_duplicates(["fov", "cell_ID"])
    counts = merged[genes].to_numpy(dtype=np.float32)
    ncount = counts.sum(axis=1)
    ngene = (counts > 0).sum(axis=1)
    keep = (ncount >= MIN_COUNTS) & (ngene >= MIN_GENES)
    counts = counts[keep]
    fov = merged.loc[keep, "fov"].to_numpy(np.int32)
    cell_id = merged.loc[keep, "cell_ID"].to_numpy(np.int32)
    xy = merged.loc[keep, ["x", "y"]].to_numpy(np.float64)
    np.savez(cache, genes=np.array(genes), fov=fov, cell_id=cell_id, xy=xy, counts=counts)
    print(
        f"[load] {sample} cells={counts.shape[0]} genes={counts.shape[1]} "
        f"fovs={len(np.unique(fov))} xy={mx},{my} in {time.time()-t0:.1f}s",
        flush=True,
    )
    return {
        "genes": genes,
        "gmap": {g: i for i, g in enumerate(genes)},
        "fov": fov,
        "cell_id": cell_id,
        "xy": xy,
        "counts": counts,
    }


def check_scale(xy: np.ndarray, sample: str, fov: int) -> float:
    if xy.shape[0] < 30:
        return float("nan")
    nn = NearestNeighbors(n_neighbors=2, algorithm="kd_tree").fit(xy)
    d_px = nn.kneighbors(xy)[0][:, 1]
    med_um = float(np.median(d_px) * PX_TO_UM)
    if not (2.0 <= med_um <= 40.0):
        span = np.ptp(xy, axis=0)
        raise RuntimeError(
            f"{sample} FOV {fov}: median NN {med_um:.3f} µm "
            f"(span_px={span.tolist()}). Refusing to treat 0.18 µm/px as valid."
        )
    return med_um


def process_sample(sample: str, n_perm: int, max_fovs: int, genes_ref: list[str]) -> None:
    ck = OUT / "checkpoints"
    ck.mkdir(parents=True, exist_ok=True)
    fov_path = ck / f"{sample}.fov.tsv"
    npz_path = ck / f"{sample}.npz"
    if fov_path.exists() and npz_path.exists():
        print(f"[resume] {sample}", flush=True)
        return
    pack = load_sample(sample)
    if pack["genes"] != genes_ref:
        raise RuntimeError(f"{sample} gene panel does not match the reference sample")
    gmap = pack["gmap"]
    logx_all = library_lognorm(pack["counts"]).astype(np.float32)
    epi_s = mean_genes(logx_all, gmap, EPI_MARKERS)
    imm_s = mean_genes(logx_all, gmap, IMM_MARKERS)
    str_s = mean_genes(logx_all, gmap, STR_MARKERS)
    S = np.vstack([epi_s, imm_s, str_s])
    comp = np.array(["epithelial", "immune", "stromal"])[S.argmax(axis=0)]
    comp[S.max(axis=0) <= 0] = "unassigned"
    is_epi = comp == "epithelial"
    cd8_raw = np.zeros(len(comp), dtype=np.float32)
    cd3_raw = np.zeros(len(comp), dtype=np.float32)
    for g in ("CD8A", "CD8B"):
        if g in gmap:
            cd8_raw += pack["counts"][:, gmap[g]]
    for g in ("CD3D", "CD3E", "CD3G"):
        if g in gmap:
            cd3_raw += pack["counts"][:, gmap[g]]
    is_cd8 = (cd8_raw > 0) & (cd3_raw > 0) & (~is_epi)
    cldn4_log = logx_all[:, gmap["CLDN4"]]
    adhesion = mean_genes(logx_all, gmap, ADHESION)
    adhesion_wo = mean_genes(logx_all, gmap, ADHESION_WO)
    keratin = mean_genes(logx_all, gmap, KERATIN)
    cd8_tx = mean_genes(logx_all, gmap, ["CD8A", "CD8B"])

    fov_ids = np.unique(pack["fov"])
    if max_fovs:
        fov_ids = fov_ids[:max_fovs]
    n_genes = len(genes_ref)
    heat_idx = [gmap[g] for g in HEATMAP if g in gmap]
    heat_names = [genes_ref[i] for i in heat_idx]
    rows = []
    z_rows = []
    hmrf_rows = []
    mods = {k: [] for k in ("hs_program", "hs_fine", "hmrf_program", "hmrf_fine")}
    pairs = []
    nulls = {k: [] for k in PERM_ENDPOINTS}
    r50 = RADIUS_UM[0] / PX_TO_UM
    r100 = RADIUS_UM[1] / PX_TO_UM
    t0 = time.time()
    scale_checked = False
    for fi, fov in enumerate(fov_ids):
        m = pack["fov"] == fov
        n_qc = int(m.sum())
        row = {
            "sample": sample,
            "patient": PATIENT[sample],
            "fov": int(fov),
            "n_qc": n_qc,
            "n_epi": int((m & is_epi).sum()),
            "n_cd8": int((m & is_cd8).sum()),
            "in_module_universe": False,
            "in_spatial_universe": False,
            "cldn4_ac_z": np.nan,
            "cldn4_ac_fdr": np.nan,
            "n_sig_genes": 0,
            "median_nn_um": np.nan,
            "rho_samecell_cldn4_cd8tx": np.nan,
            "delta_cd8_domain": np.nan,
            "domain_epi_purity": np.nan,
            "domain_max_frac": np.nan,
            "hs_fine_sign_agrees": np.nan,
            "hs_program_sign_agrees": np.nan,
        }
        for key in PERM_ENDPOINTS:
            row[f"rho_{key}"] = np.nan
        z_out = np.full(n_genes, np.nan, dtype=np.float32)
        hmrf_out = np.full(n_genes, np.nan, dtype=np.float32)
        mod_out = {k: np.zeros(n_genes, dtype=np.int16) for k in mods}
        pair = np.full((len(heat_idx), len(heat_idx)), np.nan, dtype=np.float32)
        null_row = {k: np.full(n_perm, np.nan, dtype=np.float32) for k in PERM_ENDPOINTS}
        if n_qc >= MIN_FOV_QC:
            idx = np.nonzero(m)[0]
            counts = pack["counts"][idx].astype(np.float64)
            xy = pack["xy"][idx]
            detected = ((counts > 0).sum(axis=0) >= 5)
            if counts[:, gmap["CLDN4"]].sum() > 0:
                detected[gmap["CLDN4"]] = True
            use = np.nonzero(detected & (counts.sum(axis=0) > 0))[0]
            if not scale_checked and n_qc >= 100:
                row["median_nn_um"] = check_scale(xy, sample, int(fov))
                scale_checked = True
            else:
                nn = NearestNeighbors(n_neighbors=2, algorithm="kd_tree").fit(xy)
                row["median_nn_um"] = float(np.median(nn.kneighbors(xy)[0][:, 1]) * PX_TO_UM)
            if use.size >= 10 and gmap["CLDN4"] in set(use.tolist()):
                row["in_module_universe"] = True
                X = danb_center(counts[:, use])
                W, _, _ = hotspot_graph(xy)
                ac, WX = autocorr_z(X, W)
                ac_p = stats.norm.sf(ac)
                ac_q = bh_fdr(ac_p)
                local_of = {int(g): j for j, g in enumerate(use)}
                c_local = local_of[gmap["CLDN4"]]
                row["cldn4_ac_z"] = float(ac[c_local])
                row["cldn4_ac_fdr"] = float(ac_q[c_local])
                eg2, _ = eg2_from_W(X, W)
                z_used = lc_z_against_one(X, WX, eg2, c_local)
                z_out[use] = z_used.astype(np.float32)
                sig_local = np.nonzero(ac_q < GENE_FDR)[0]
                row["n_sig_genes"] = int(sig_local.size)
                if sig_local.size >= 2:
                    Zsig = pairwise_from_products(X[:, sig_local], WX[:, sig_local], eg2[sig_local])
                    prog = hotspot_modules(Zsig, 20, PAIR_FDR, core_only=True)
                    fine = hotspot_modules(Zsig, 4, PAIR_FDR, core_only=False)
                    mod_out["hs_program"][use[sig_local]] = prog.astype(np.int16)
                    mod_out["hs_fine"][use[sig_local]] = fine.astype(np.int16)
                # HMRF on top autocorrelated genes plus watchlist.
                order = np.argsort(-ac)
                top = [int(use[j]) for j in order if ac_q[j] < GENE_FDR][:HMRF_TOP]
                forced = [gmap[g] for g in WATCHLIST if g in gmap and gmap[g] in local_of]
                hmrf_genes = list(dict.fromkeys(top + forced))
                if len(hmrf_genes) >= 4 and gmap["CLDN4"] in local_of:
                    A = binary_knn_adj(xy, k=HMRF_K)
                    # CLDN4 is included via the watchlist even if it is outside the top 150.
                    if gmap["CLDN4"] not in hmrf_genes:
                        hmrf_genes.append(gmap["CLDN4"])
                    P = np.column_stack([
                        hmrf_high_prob(X[:, local_of[g]], A) for g in hmrf_genes
                    ])
                    cpos = hmrf_genes.index(gmap["CLDN4"])
                    p_ref = P[:, cpos]
                    p_ref_c = p_ref - p_ref.mean()
                    denom_ref = float(np.dot(p_ref_c, p_ref_c))
                    if denom_ref > 0 and float(p_ref.std()) > 1e-6:
                        for j, g in enumerate(hmrf_genes):
                            if j == cpos:
                                continue
                            pj = P[:, j]
                            if float(pj.std()) <= 1e-6:
                                continue
                            pj_c = pj - pj.mean()
                            hmrf_out[g] = float(np.dot(p_ref_c, pj_c) / np.sqrt(denom_ref * np.dot(pj_c, pj_c)))
                    h_prog = correlation_modules(P, min_size=15, min_corr=HMRF_CORR)
                    h_fine = correlation_modules(P, min_size=4, min_corr=HMRF_CORR)
                    for j, g in enumerate(hmrf_genes):
                        mod_out["hmrf_program"][g] = np.int16(h_prog[j])
                        mod_out["hmrf_fine"][g] = np.int16(h_fine[j])
                if len(heat_idx) >= 2:
                    h_local = [local_of[g] for g in heat_idx if g in local_of]
                    h_global = [g for g in heat_idx if g in local_of]
                    if len(h_local) >= 2:
                        Zh = pairwise_from_products(X[:, h_local], WX[:, h_local], eg2[h_local])
                        for a, ga in enumerate(h_global):
                            for b, gb in enumerate(h_global):
                                pair[heat_idx.index(ga), heat_idx.index(gb)] = Zh[a, b]
                # Scores on all QC cells in the FOV, then subset.
                ref = X[:, c_local]
                def score_from(mod_vec, kind):
                    code = int(mod_vec[gmap["CLDN4"]])
                    if code <= 0:
                        return None, np.nan
                    members_global = np.nonzero(mod_vec == code)[0]
                    members_local = [local_of[int(g)] for g in members_global if int(g) in local_of]
                    if not members_local:
                        return None, np.nan
                    if kind == "hotspot":
                        sc, agrees = hotspot_module_score(X, W, np.array(members_local), c_local)
                    else:
                        sc = X[:, members_local].mean(axis=1)
                        sc, agrees = align_to_reference(sc, ref)
                    return sc, float(agrees)
                sc_hs_f, ag_f = score_from(mod_out["hs_fine"], "hotspot")
                sc_hs_p, ag_p = score_from(mod_out["hs_program"], "hotspot")
                sc_hm_f, _ = score_from(mod_out["hmrf_fine"], "hmrf")
                sc_hm_p, _ = score_from(mod_out["hmrf_program"], "hmrf")
                row["hs_fine_sign_agrees"] = ag_f
                row["hs_program_sign_agrees"] = ag_p
                epi_local = is_epi[idx]
                cd8_local = is_cd8[idx]
                same = spearman(cldn4_log[idx], cd8_tx[idx])
                row["rho_samecell_cldn4_cd8tx"] = same
                # Domain composition.
                lab = domain_hmrf(logx_all[idx].astype(np.float64), xy)
                if lab is not None and int(cd8_local.sum()) >= MIN_CD8:
                    means = []
                    for s in range(5):
                        msk = lab == s
                        if int(msk.sum()) >= 20:
                            means.append((float(cldn4_log[idx][msk].mean()), s, int(msk.sum())))
                    if means:
                        _, high, _ = max(means, key=lambda t: t[0])
                        high_m = lab == high
                        frac_high = float(cd8_local[high_m].mean()) if high_m.any() else np.nan
                        frac_other = float(cd8_local[~high_m].mean()) if (~high_m).any() else np.nan
                        row["delta_cd8_domain"] = frac_high - frac_other
                        row["domain_epi_purity"] = float(epi_local[high_m].mean())
                        sizes = np.bincount(lab, minlength=5).astype(float)
                        row["domain_max_frac"] = float(sizes.max() / sizes.sum())
                spatial = int(epi_local.sum()) >= MIN_EPI and int(cd8_local.sum()) >= MIN_CD8
                row["in_spatial_universe"] = bool(spatial)
                if spatial:
                    epi_xy = xy[epi_local]
                    cd8_xy = xy[cd8_local]
                    lo = xy.min(axis=0)
                    span = np.ptp(xy, axis=0)
                    span[span == 0] = 1.0
                    tree0 = cKDTree(cd8_xy)
                    c50 = ball_count(tree0, epi_xy, r50)
                    c100 = ball_count(tree0, epi_xy, r100)
                    cl = cldn4_log[idx][epi_local]
                    ker = None if keratin is None else keratin[idx][epi_local]
                    ad = None if adhesion is None else adhesion[idx][epi_local]
                    adw = None if adhesion_wo is None else adhesion_wo[idx][epi_local]
                    resid = None if ker is None else rank_resid(cl, ker)
                    scores = {
                        "cldn4_50": cl,
                        "cldn4_100": cl,
                        "keratin_50": ker,
                        "resid_50": resid,
                        "adhesion_50": ad,
                        "adhesion_wo_50": adw,
                        "hs_fine_50": None if sc_hs_f is None else sc_hs_f[epi_local],
                        "hs_program_50": None if sc_hs_p is None else sc_hs_p[epi_local],
                        "hmrf_fine_50": None if sc_hm_f is None else sc_hm_f[epi_local],
                        "hmrf_program_50": None if sc_hm_p is None else sc_hm_p[epi_local],
                    }
                    obs_counts = {"50": c50, "100": c100}
                    for key, vec in scores.items():
                        if vec is None:
                            continue
                        radius_key = "100" if key.endswith("100") else "50"
                        row[f"rho_{key}"] = spearman(vec, obs_counts[radius_key])
                    ss = np.random.SeedSequence([SEED, zlib.adler32(sample.encode()), int(fov)])
                    rng = np.random.default_rng(ss)
                    shifts = rng.random((n_perm, 2)) * span
                    for b in range(n_perm):
                        moved = torus_shift(cd8_xy, lo, span, shifts[b])
                        tree = cKDTree(moved)
                        cb50 = ball_count(tree, epi_xy, r50)
                        cb100 = ball_count(tree, epi_xy, r100)
                        for key, vec in scores.items():
                            if vec is None or not np.isfinite(row[f"rho_{key}"]):
                                continue
                            cnt = cb100 if key.endswith("100") else cb50
                            null_row[key][b] = spearman(vec, cnt)
        rows.append(row)
        z_rows.append(z_out)
        hmrf_rows.append(hmrf_out)
        for k in mods:
            mods[k].append(mod_out[k])
        pairs.append(pair)
        for k in PERM_ENDPOINTS:
            nulls[k].append(null_row[k])
        if (fi + 1) % 5 == 0 or fi + 1 == len(fov_ids):
            print(f"  {sample} FOV {fi+1}/{len(fov_ids)} ({time.time()-t0:.0f}s)", flush=True)

    fov_df = pd.DataFrame(rows)
    fov_df.to_csv(fov_path, sep="\t", index=False)
    np.savez_compressed(
        npz_path,
        fov=fov_df["fov"].to_numpy(np.int32),
        z=np.vstack(z_rows),
        hmrf_corr=np.vstack(hmrf_rows),
        pair=np.stack(pairs),
        heat_names=np.array(heat_names),
        **{f"mod_{k}": np.vstack(v) for k, v in mods.items()},
        **{f"null_{k}": np.vstack(v) for k, v in nulls.items()},
    )
    print(f"[done] {sample} fovs={len(fov_df)} in {time.time()-t0:.0f}s", flush=True)


def patient_mean(values: np.ndarray, patients: np.ndarray) -> float:
    means = []
    for patient in PATIENT_ORDER:
        v = values[patients == patient]
        v = v[np.isfinite(v)]
        if v.size == 0:
            return float("nan")
        means.append(float(v.mean()))
    return float(np.mean(means))


def perm_p_less(obs: float, null: np.ndarray) -> float:
    null = np.asarray(null, dtype=np.float64)
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size == 0:
        return float("nan")
    return float((1 + np.sum(null <= obs)) / (null.size + 1))


def class_counts(mod: np.ndarray, genes: list[str], mask: np.ndarray) -> dict:
    cldn = genes.index("CLDN4")
    counts = {
        "unassigned": 0,
        "adhesion_like": 0,
        "adhesion_and_keratin": 0,
        "keratin_without_adhesion": 0,
        "other": 0,
    }
    for i in np.nonzero(mask)[0]:
        code = int(mod[i, cldn])
        if code <= 0:
            counts["unassigned"] += 1
            continue
        names = module_name_set(mod[i], genes, code)
        counts[classify_members(names)] += 1
    counts["n"] = int(mask.sum())
    return counts


def co_membership(mod: np.ndarray, genes: list[str], mask: np.ndarray) -> pd.DataFrame:
    cldn = genes.index("CLDN4")
    rows = []
    assigned = mask & (mod[:, cldn] > 0)
    for j, gene in enumerate(genes):
        eligible = assigned & (mod[:, j] != 0)
        same = eligible & (mod[:, j] == mod[:, cldn])
        denom = int(eligible.sum())
        rows.append({
            "gene": gene,
            "n_eligible": denom,
            "n_comember": int(same.sum()),
            "frac_comember": float(same.sum() / denom) if denom else np.nan,
        })
    return pd.DataFrame(rows)


def aggregate(n_perm: int) -> None:
    ck = OUT / "checkpoints"
    frames = []
    packs = []
    for sample in SAMPLES:
        frames.append(pd.read_csv(ck / f"{sample}.fov.tsv", sep="\t"))
        packs.append(np.load(ck / f"{sample}.npz", allow_pickle=False))
    fov = pd.concat(frames, ignore_index=True)
    genes = json.loads((OUT / "genes.json").read_text())
    z = np.vstack([p["z"] for p in packs])
    hmrf_corr = np.vstack([p["hmrf_corr"] for p in packs])
    mod = {k: np.vstack([p[f"mod_{k}"] for p in packs]) for k in (
        "hs_program", "hs_fine", "hmrf_program", "hmrf_fine"
    )}
    null = {k: np.vstack([p[f"null_{k}"] for p in packs]) for k in PERM_ENDPOINTS}
    heat_names = packs[0]["heat_names"].astype(str).tolist()
    pair = np.concatenate([p["pair"] for p in packs], axis=0)
    if len(fov) != len(z):
        raise RuntimeError("checkpoint row mismatch")
    tab = OUT / "tables"
    fig = OUT / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    fov.to_csv(tab / "fov_metrics.tsv", sep="\t", index=False)

    module_mask = fov["in_module_universe"].to_numpy(bool)
    spatial_mask = fov["in_spatial_universe"].to_numpy(bool)
    patients = fov["patient"].to_numpy()
    samples = fov["sample"].to_numpy()

    # Coupling with CLDN4.
    cldn = genes.index("CLDN4")
    coupling_rows = []
    for j, gene in enumerate(genes):
        if gene == "CLDN4":
            continue
        col = z[:, j]
        use = module_mask & np.isfinite(col)
        vals = col[use]
        hcol = hmrf_corr[:, j]
        huse = module_mask & np.isfinite(hcol)
        hvals = hcol[huse]
        coupling_rows.append({
            "gene": gene,
            "n_fov": int(use.sum()),
            "median_lc_z": float(np.median(vals)) if vals.size else np.nan,
            "q25_lc_z": float(np.quantile(vals, 0.25)) if vals.size else np.nan,
            "q75_lc_z": float(np.quantile(vals, 0.75)) if vals.size else np.nan,
            "n_fov_hmrf": int(huse.sum()),
            "median_hmrf_corr": float(np.median(hvals)) if hvals.size else np.nan,
        })
    coupling = pd.DataFrame(coupling_rows).sort_values("median_lc_z", ascending=False)
    for key, label in (
        ("hs_program", "hotspot_program"),
        ("hs_fine", "hotspot_fine"),
        ("hmrf_program", "hmrf_program"),
        ("hmrf_fine", "hmrf_fine"),
    ):
        co = co_membership(mod[key], genes, module_mask).rename(columns={
            "n_eligible": f"n_eligible_{label}",
            "n_comember": f"n_comember_{label}",
            "frac_comember": f"frac_comember_{label}",
        })
        coupling = coupling.merge(co.drop(columns=["gene"]).assign(gene=co["gene"]), on="gene", how="left")
    coupling.to_csv(tab / "cldn4_spatial_coupling.tsv", sep="\t", index=False)

    # Membership long table for CLDN4's modules.
    mem_rows = []
    for key in mod:
        for i in np.nonzero(module_mask)[0]:
            code = int(mod[key][i, cldn])
            if code <= 0:
                continue
            for j, gene in enumerate(genes):
                if int(mod[key][i, j]) == code:
                    mem_rows.append({
                        "sample": fov.loc[i, "sample"],
                        "patient": fov.loc[i, "patient"],
                        "fov": int(fov.loc[i, "fov"]),
                        "method": key,
                        "gene": gene,
                    })
    pd.DataFrame(mem_rows).to_csv(tab / "cldn4_module_members.tsv", sep="\t", index=False)

    classes = {}
    for key in mod:
        classes[key] = {
            "module_fovs": class_counts(mod[key], genes, module_mask),
            "spatial_fovs": class_counts(mod[key], genes, spatial_mask),
        }

    # Pairwise median local-correlation heatmap.
    pair_med = np.full((len(heat_names), len(heat_names)), np.nan)
    for a in range(len(heat_names)):
        for b in range(len(heat_names)):
            vals = pair[:, a, b]
            vals = vals[np.isfinite(vals)]
            if vals.size:
                pair_med[a, b] = float(np.median(vals))
    pd.DataFrame(pair_med, index=heat_names, columns=heat_names).to_csv(
        tab / "median_local_correlation_z.tsv", sep="\t"
    )

    endpoint_stats = {}
    for key in PERM_ENDPOINTS:
        obs = fov[f"rho_{key}"].to_numpy(dtype=float)
        mask = np.isfinite(obs)
        T = patient_mean(obs, patients)
        null_T = np.array([
            patient_mean(np.where(mask, null[key][:, b], np.nan), patients)
            for b in range(null[key].shape[1])
        ])
        per_patient = {}
        per_sample = {}
        for patient in PATIENT_ORDER:
            v = obs[(patients == patient) & mask]
            per_patient[patient] = {
                "n_fov": int(v.size),
                "mean_rho": float(v.mean()) if v.size else None,
            }
        for sample in SAMPLES:
            v = obs[(samples == sample) & mask]
            per_sample[sample] = {
                "n_fov": int(v.size),
                "mean_rho": float(v.mean()) if v.size else None,
            }
        n_pat_neg = sum(
            1 for p in PATIENT_ORDER
            if per_patient[p]["mean_rho"] is not None and per_patient[p]["mean_rho"] < 0
        )
        n_pat = sum(1 for p in PATIENT_ORDER if per_patient[p]["mean_rho"] is not None)
        n_samp_neg = sum(
            1 for s in SAMPLES
            if per_sample[s]["mean_rho"] is not None and per_sample[s]["mean_rho"] < 0
        )
        n_samp = sum(1 for s in SAMPLES if per_sample[s]["mean_rho"] is not None)
        sign_p = None
        if n_pat:
            sign_p = float(stats.binomtest(n_pat_neg, n_pat, 0.5, alternative="greater").pvalue)
        patient_rhos = np.array([
            per_patient[p]["mean_rho"] if per_patient[p]["mean_rho"] is not None else np.nan
            for p in PATIENT_ORDER
        ], dtype=float)
        wilcox_p = None
        finite_patients = patient_rhos[np.isfinite(patient_rhos)]
        if finite_patients.size >= 5 and np.any(finite_patients != 0):
            wilcox_p = float(stats.wilcoxon(finite_patients, alternative="less").pvalue)
        endpoint_stats[key] = {
            "patient_mean_rho": T,
            "permutation_p_one_sided_less": perm_p_less(T, null_T),
            "n_perm": int(np.isfinite(null_T).sum()),
            "n_fov": int(mask.sum()),
            "median_fov_rho": float(np.median(obs[mask])) if mask.any() else None,
            "patients_rho_lt_0": f"{n_pat_neg}/{n_pat}",
            "samples_rho_lt_0": f"{n_samp_neg}/{n_samp}",
            "sign_test_p_patients": sign_p,
            "wilcoxon_p_patients_less": wilcox_p,
            "per_patient": per_patient,
            "per_sample": per_sample,
            "null_patient_mean": null_T.tolist(),
        }

    # Example FOV: closest to the median primary rho, with enough cells.
    obs = fov["rho_cldn4_50"].to_numpy(float)
    eligible = (
        spatial_mask
        & np.isfinite(obs)
        & (fov["n_epi"].to_numpy() >= 200)
        & (fov["n_cd8"].to_numpy() >= 20)
    )
    if not eligible.any():
        eligible = spatial_mask & np.isfinite(obs)
    med = float(np.median(obs[eligible]))
    cand = np.nonzero(eligible)[0]
    order = sorted(cand, key=lambda i: (abs(obs[i] - med), str(fov.loc[i, "sample"]), int(fov.loc[i, "fov"])))
    ex = int(order[0])
    example = {
        "sample": str(fov.loc[ex, "sample"]),
        "patient": str(fov.loc[ex, "patient"]),
        "fov": int(fov.loc[ex, "fov"]),
        "rho_cldn4_50": float(obs[ex]),
        "median_rho_cldn4_50": med,
        "n_epi": int(fov.loc[ex, "n_epi"]),
        "n_cd8": int(fov.loc[ex, "n_cd8"]),
    }

    panel = json.loads((OUT / "panel_status.json").read_text())
    summary = {
        "dataset": "He et al. 2022 CosMx SMI NSCLC 960-plex, NanoString public flat files",
        "n_sections": 8,
        "n_patients": 5,
        "px_to_um": PX_TO_UM,
        "n_fov_rows": int(len(fov)),
        "n_module_fovs": int(module_mask.sum()),
        "n_spatial_fovs": int(spatial_mask.sum()),
        "panel": panel,
        "module_classes": classes,
        "endpoints": {k: {kk: vv for kk, vv in st.items() if kk != "null_patient_mean"} for k, st in endpoint_stats.items()},
        "top_positive_partners": coupling.head(15)[["gene", "n_fov", "median_lc_z"]].to_dict(orient="records"),
        "top_negative_partners": coupling.dropna(subset=["median_lc_z"]).tail(10)[["gene", "n_fov", "median_lc_z"]].to_dict(orient="records"),
        "example_fov": example,
        "primary_endpoint": "cldn4_50",
        "module_endpoint": "hs_fine_50",
        "null": "toroidal shift of CD8 coordinates within the FOV bounding box",
        "n_perm_requested": n_perm,
    }
    (OUT / "stats.json").write_text(json.dumps(summary, indent=2))
    # Keep nulls for the figure without putting 199 x endpoints into the prose json twice.
    (tab / "patient_endpoint_means.tsv").write_text(patient_table(endpoint_stats))
    plot_all(fov, coupling, pair_med, heat_names, endpoint_stats, example, genes)
    print("[aggregate] wrote tables, figures, stats.json", flush=True)


def patient_table(endpoint_stats: dict) -> str:
    rows = []
    for key, st in endpoint_stats.items():
        for patient, rec in st["per_patient"].items():
            rows.append({
                "endpoint": key,
                "patient": patient,
                "n_fov": rec["n_fov"],
                "mean_rho": rec["mean_rho"],
                "patient_mean_of_means": st["patient_mean_rho"],
                "perm_p_less": st["permutation_p_one_sided_less"],
            })
    return pd.DataFrame(rows).to_csv(sep="\t", index=False)


def plot_all(fov, coupling, pair_med, heat_names, endpoint_stats, example, genes):
    figdir = OUT / "figures"
    plt.rcParams.update({
        "font.size": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.bbox": "tight",
    })
    # Forest of the primary gene test and the specificity scores.
    show = [
        ("cldn4_50", "CLDN4", "#E69F00"),
        ("hs_fine_50", "Hotspot fine module", "#009E73"),
        ("keratin_50", "KRT7/8/18/19", "#0072B2"),
        ("resid_50", "CLDN4 residual on keratin", "#CC79A7"),
    ]
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    for pi, patient in enumerate(PATIENT_ORDER):
        sub = fov[(fov["patient"] == patient) & fov["in_spatial_universe"]]
        for j, (key, label, color) in enumerate(show):
            vals = sub[f"rho_{key}"].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            jitter = (j - 1.5) * 0.12
            ax.scatter(vals, np.full(vals.size, pi) + jitter, s=12, color=color, alpha=0.35, linewidths=0, zorder=2)
            if vals.size:
                ax.scatter([vals.mean()], [pi + jitter], s=46, color=color, zorder=3, label=label if pi == 0 else None)
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(range(len(PATIENT_ORDER)))
    ax.set_yticklabels(PATIENT_ORDER)
    ax.set_xlabel("Spearman ρ, epithelial cells vs CD8 count within 50 µm")
    ax.set_title("Patient means (large) and FOVs (small)")
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    fig.savefig(figdir / "patient_rho_forest.png")
    fig.savefig(figdir / "patient_rho_forest.pdf")
    plt.close(fig)

    # Co-membership for the watchlist, Hotspot fine vs program.
    watch = [g for g in dict.fromkeys(WATCHLIST) if g != "CLDN4" and g in set(coupling["gene"])]
    sub = coupling.set_index("gene").loc[watch]
    fig, ax = plt.subplots(figsize=(8, 6.2))
    y = np.arange(len(watch))
    ax.barh(y - 0.18, sub["frac_comember_hotspot_fine"], height=0.36, color="#009E73", label="Hotspot fine (≥4 genes)")
    ax.barh(y + 0.18, sub["frac_comember_hotspot_program"], height=0.36, color="#E69F00", label="Hotspot program (≥20, core)")
    ax.set_yticks(y)
    ax.set_yticklabels(watch)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of module FOVs where the gene shares CLDN4's module")
    ax.invert_yaxis()
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("Who sits in CLDN4's Hotspot module")
    fig.savefig(figdir / "comembership.png")
    fig.savefig(figdir / "comembership.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    lim = np.nanmax(np.abs(pair_med))
    lim = 8 if not np.isfinite(lim) or lim == 0 else float(lim)
    im = ax.imshow(pair_med, cmap="RdBu_r", vmin=-lim, vmax=lim)
    ax.set_xticks(range(len(heat_names)))
    ax.set_yticks(range(len(heat_names)))
    ax.set_xticklabels(heat_names, rotation=90)
    ax.set_yticklabels(heat_names)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Median Hotspot local-correlation Z")
    ax.set_title("Expression co-variation in space\n(not the CD8-density exclusion test)")
    fig.savefig(figdir / "local_correlation_heatmap.png")
    fig.savefig(figdir / "local_correlation_heatmap.pdf")
    plt.close(fig)

    st = endpoint_stats["cldn4_50"]
    null = np.asarray(st["null_patient_mean"], dtype=float)
    null = null[np.isfinite(null)]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.hist(null, bins=30, color="#B0B0B0", edgecolor="white")
    ax.axvline(st["patient_mean_rho"], color="#E69F00", lw=2, label=f"observed {st['patient_mean_rho']:.3f}")
    ax.set_xlabel("Patient-mean Spearman ρ")
    ax.set_ylabel("Toroidal shifts")
    ax.set_title(f"CLDN4 vs CD8 count, 50 µm\none-sided p = {st['permutation_p_one_sided_less']:.3f}")
    ax.legend(frameon=False)
    fig.savefig(figdir / "primary_perm_null.png")
    fig.savefig(figdir / "primary_perm_null.pdf")
    plt.close(fig)

    plot_example_fov(example)
    # Silence unused in case genes is useful for debugging plots.
    _ = genes


def plot_example_fov(example: dict) -> None:
    sample = example["sample"]
    fov = int(example["fov"])
    pack = load_sample(sample)
    gmap = pack["gmap"]
    m = pack["fov"] == fov
    logx = library_lognorm(pack["counts"][m])
    epi_s = mean_genes(logx, gmap, EPI_MARKERS)
    imm_s = mean_genes(logx, gmap, IMM_MARKERS)
    str_s = mean_genes(logx, gmap, STR_MARKERS)
    S = np.vstack([epi_s, imm_s, str_s])
    comp = np.array(["epithelial", "immune", "stromal"])[S.argmax(axis=0)]
    comp[S.max(axis=0) <= 0] = "unassigned"
    is_epi = comp == "epithelial"
    cd8 = np.zeros(is_epi.size)
    cd3 = np.zeros(is_epi.size)
    for g in ("CD8A", "CD8B"):
        cd8 += pack["counts"][m][:, gmap[g]]
    for g in ("CD3D", "CD3E", "CD3G"):
        cd3 += pack["counts"][m][:, gmap[g]]
    is_cd8 = (cd8 > 0) & (cd3 > 0) & (~is_epi)
    xy = pack["xy"][m] * PX_TO_UM
    cldn = logx[:, gmap["CLDN4"]]
    fig, ax = plt.subplots(figsize=(6.2, 5.6))
    other = ~is_epi & ~is_cd8
    ax.scatter(xy[other, 0], xy[other, 1], s=3, c="#E6E6E6", linewidths=0, label="other")
    sc = ax.scatter(
        xy[is_epi, 0], xy[is_epi, 1], s=8, c=cldn[is_epi], cmap="cividis",
        linewidths=0, label="RNA-epithelial",
    )
    ax.scatter(xy[is_cd8, 0], xy[is_cd8, 1], s=16, facecolors="none", edgecolors="#0072B2", linewidths=0.6, label="CD8 T")
    fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label="CLDN4 log-norm")
    ax.set_aspect("equal")
    ax.set_xlabel("µm")
    ax.set_ylabel("µm")
    ax.set_title(
        f"{sample} FOV {fov}  (near median ρ)\n"
        f"CLDN4 vs CD8 50 µm ρ = {example['rho_cldn4_50']:.2f}"
    )
    ax.legend(frameon=False, markerscale=2, loc="best", fontsize=8)
    fig.savefig(OUT / "figures" / "example_fov.png")
    fig.savefig(OUT / "figures" / "example_fov.pdf")
    plt.close(fig)


def write_panel_status(genes: list[str]) -> None:
    present = set(genes)
    status = {
        "n_panel_genes": len(genes),
        "classical_tj_present": [g for g in CLASSICAL_TJ if g in present],
        "classical_tj_absent": [g for g in CLASSICAL_TJ if g not in present],
        "adhesion_present": [g for g in ADHESION if g in present],
        "adhesion_absent": [g for g in ADHESION if g not in present],
        "keratin_present": [g for g in KERATIN if g in present],
        "note": (
            "CLDN4 is the only claudin queried that is on this 960-plex. "
            "CDH1 is an adherens-junction gene. EPCAM and TACSTD2 are epithelial "
            "adhesion genes, not tight-junction genes. A multi-gene claudin/occludin/ZO "
            "module cannot be identified on this panel."
        ),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "panel_status.json").write_text(json.dumps(status, indent=2))
    (OUT / "genes.json").write_text(json.dumps(genes))
    pd.DataFrame({
        "gene": CLASSICAL_TJ + [g for g in ADHESION + KERATIN + ENDO_JUNC if g not in CLASSICAL_TJ],
        "group": (
            ["classical_tj_query"] * len(CLASSICAL_TJ)
            + ["adhesion"] * len([g for g in ADHESION if g not in CLASSICAL_TJ])
            + ["keratin"] * len(KERATIN)
            + ["endothelial_junction"] * len(ENDO_JUNC)
        ),
        "on_960_panel": [
            g in present
            for g in CLASSICAL_TJ + [g for g in ADHESION + KERATIN + ENDO_JUNC if g not in CLASSICAL_TJ]
        ],
    }).to_csv(OUT / "tables" / "panel_gene_status.tsv", sep="\t", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--perms", type=int, default=N_PERM_DEFAULT)
    ap.add_argument("--samples", nargs="*")
    ap.add_argument("--max-fovs", type=int, default=0)
    ap.add_argument("--skip-aggregate", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return
    sys.setrecursionlimit(10000)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    samples = args.samples or SAMPLES
    ref = load_sample(samples[0])
    write_panel_status(ref["genes"])
    # Reference sample was not checkpointed yet; process_sample reloads the cache.
    for sample in samples:
        process_sample(sample, args.perms, args.max_fovs, ref["genes"])
    if args.samples or args.max_fovs or args.skip_aggregate:
        if not args.skip_aggregate and not args.max_fovs and set(args.samples or []) == set(SAMPLES):
            aggregate(args.perms)
        else:
            print("[stop] partial run, aggregation skipped", flush=True)
        return
    aggregate(args.perms)


if __name__ == "__main__":
    main()
