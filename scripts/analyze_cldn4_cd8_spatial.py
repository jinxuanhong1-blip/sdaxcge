#!/usr/bin/env python3
"""CLDN4 vs CD8 spatial statistics on public GSE307534 Visium CytAssist sections.

CLDN4-only. No private gene lists. Statistics are computed from downloaded
filtered spot matrices; nothing is imputed or fabricated.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import deque
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread
from scipy.spatial import ConvexHull, QhullError, cKDTree

GENES = ("CLDN4", "CD8A", "NKG7", "GNLY", "EPCAM", "KRT8")
K_NBR = 6
GI_Z = 1.96
N_PERM_DEFAULT = 199
MIN_Q_SPOTS = 15
MIN_TUMOR_CC = 40
RNG_SEED = 307534


def load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_features(path: Path) -> list[str]:
    names: list[str] = []
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            names.append(parts[1] if len(parts) > 1 else parts[0])
    return names


def read_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as f:
        return [line.strip() for line in f if line.strip()]


def read_positions(spatial_dir: Path) -> pd.DataFrame:
    for name in ("tissue_positions.csv", "tissue_positions_list.csv"):
        p = spatial_dir / name
        if not p.is_file():
            continue
        df = pd.read_csv(p, header=None)
        first = str(df.iloc[0, 0])
        if first == "barcode" or first.lower().startswith("barcode"):
            df = pd.read_csv(p)
        else:
            df.columns = [
                "barcode",
                "in_tissue",
                "array_row",
                "array_col",
                "pxl_row_in_fullres",
                "pxl_col_in_fullres",
            ]
        return df
    raise FileNotFoundError(f"no tissue positions in {spatial_dir}")


def microns_per_pixel(spatial_dir: Path) -> float | None:
    p = spatial_dir / "scalefactors_json.json"
    if not p.is_file():
        return None
    sc = json.loads(p.read_text())
    d = sc.get("spot_diameter_fullres")
    if d and d > 0:
        return 55.0 / float(d)
    return None


def load_section(section_dir: Path) -> dict:
    feat_path = section_dir / "filtered_feature_bc_matrix" / "features.tsv.gz"
    bc_path = section_dir / "filtered_feature_bc_matrix" / "barcodes.tsv.gz"
    mtx_path = section_dir / "filtered_feature_bc_matrix" / "matrix.mtx.gz"
    genes = read_features(feat_path)
    barcodes = read_barcodes(bc_path)
    with gzip.open(mtx_path, "rb") as fh:
        mat = mmread(fh)
    if sparse.issparse(mat):
        mat = mat.tocsr()
    else:
        mat = sparse.csr_matrix(mat)
    if mat.shape[0] != len(genes) and mat.shape[1] == len(genes):
        mat = mat.T.tocsr()
    pos = read_positions(section_dir / "spatial")
    pos = pos.copy()
    pos["barcode"] = pos["barcode"].astype(str)
    pos = pos.set_index("barcode")
    bc_series = pd.Index(barcodes)
    keep_bc = [b for b in barcodes if b in pos.index]
    if len(keep_bc) != len(barcodes):
        idx = [i for i, b in enumerate(barcodes) if b in pos.index]
        mat = mat[:, idx]
        bc_series = pd.Index([barcodes[i] for i in idx])
    aligned = pos.loc[bc_series]
    in_tissue = aligned["in_tissue"].astype(int).to_numpy() == 1
    if in_tissue.any() and in_tissue.sum() < len(in_tissue):
        mat = mat[:, in_tissue]
        aligned = aligned.loc[in_tissue]
        bc_series = bc_series[in_tissue]
    um_per_px = microns_per_pixel(section_dir / "spatial")
    xy_px = np.column_stack(
        [
            aligned["pxl_col_in_fullres"].astype(float).to_numpy(),
            aligned["pxl_row_in_fullres"].astype(float).to_numpy(),
        ]
    )
    if um_per_px is None:
        tree = cKDTree(xy_px)
        nn = tree.query(xy_px, k=2)[0][:, 1]
        med = float(np.median(nn[nn > 0])) if np.any(nn > 0) else 1.0
        um_per_px = 100.0 / med if med > 0 else 1.0
    xy_um = xy_px * um_per_px
    return {
        "matrix": mat,
        "genes": genes,
        "barcodes": bc_series.to_numpy(),
        "xy_px": xy_px,
        "xy_um": xy_um,
        "array": np.column_stack(
            [
                aligned["array_row"].astype(float).to_numpy(),
                aligned["array_col"].astype(float).to_numpy(),
            ]
        ),
        "um_per_px": um_per_px,
    }


def gene_index_map(genes: list[str]) -> dict[str, list[int]]:
    idx: dict[str, list[int]] = {}
    for i, g in enumerate(genes):
        idx.setdefault(g, []).append(i)
    return idx


def extract_counts(mat: sparse.csr_matrix, genes: list[str], wanted: tuple[str, ...]) -> dict[str, np.ndarray]:
    gmap = gene_index_map(genes)
    out: dict[str, np.ndarray] = {}
    for g in wanted:
        ids = gmap.get(g, [])
        if not ids:
            continue
        if len(ids) == 1:
            out[g] = np.asarray(mat[ids[0], :].todense()).ravel().astype(float)
        else:
            out[g] = np.asarray(mat[ids, :].sum(axis=0)).ravel().astype(float)
    return out


def log_norm(counts: dict[str, np.ndarray], lib: np.ndarray) -> dict[str, np.ndarray]:
    scale = np.zeros_like(lib, dtype=float)
    nz = lib > 0
    scale[nz] = 1e4 / lib[nz]
    return {g: np.log1p(v * scale) for g, v in counts.items()}


def zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    sd = x.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return np.zeros_like(x)
    return (x - x.mean()) / sd


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    X = np.column_stack([np.ones(len(x)), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def knn_indices(xy: np.ndarray, k: int = K_NBR) -> np.ndarray:
    n = len(xy)
    kk = min(k + 1, n)
    tree = cKDTree(xy)
    _, idx = tree.query(xy, k=kk)
    if kk == 1:
        return np.zeros((n, 0), dtype=int)
    return np.asarray(idx[:, 1:], dtype=int)


def apply_W(values: np.ndarray, nbr: np.ndarray) -> np.ndarray:
    if nbr.size == 0:
        return np.zeros_like(values, dtype=float)
    return values[nbr].mean(axis=1)


def bivariate_moran(x: np.ndarray, y: np.ndarray, nbr: np.ndarray) -> float:
    zx = x - x.mean()
    zy = y - y.mean()
    den = float(zx @ zx)
    if den <= 0:
        return float("nan")
    return float(zx @ apply_W(zy, nbr) / den)


def lee_L(x: np.ndarray, y: np.ndarray, nbr: np.ndarray) -> float:
    zx = x - x.mean()
    zy = y - y.mean()
    n = len(x)
    den = float((zx @ zx) * (zy @ zy) / n)
    if den <= 0:
        return float("nan")
    return float(apply_W(zx, nbr) @ apply_W(zy, nbr) / den)


def moran_i(x: np.ndarray, nbr: np.ndarray) -> float:
    return bivariate_moran(x, x, nbr)


def gi_star(x: np.ndarray, nbr: np.ndarray) -> np.ndarray:
    n = len(x)
    if n < 4:
        return np.full(n, np.nan)
    incl = np.concatenate([np.arange(n)[:, None], nbr], axis=1)
    k = incl.shape[1]
    xbar = float(x.mean())
    s = float(x.std(ddof=1))
    if s == 0 or not np.isfinite(s):
        return np.zeros(n)
    sum_wx = x[incl].sum(axis=1)
    numer = sum_wx - xbar * k
    denom = s * math.sqrt((n * k - k * k) / (n - 1))
    if denom == 0:
        return np.zeros(n)
    return numer / denom


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 4 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return float("nan"), float("nan")
    r, p = stats.spearmanr(x, y)
    return float(r), float(p)


def empirical_p(obs: float, null: np.ndarray, alternative: str) -> float:
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or null.size == 0:
        return float("nan")
    n = null.size
    if alternative == "less":
        return float((1 + np.sum(null <= obs)) / (n + 1))
    if alternative == "greater":
        return float((1 + np.sum(null >= obs)) / (n + 1))
    p_hi = (1 + np.sum(null >= obs)) / (n + 1)
    p_lo = (1 + np.sum(null <= obs)) / (n + 1)
    return float(min(1.0, 2 * min(p_hi, p_lo)))


def nearest_mean_dist(query_xy: np.ndarray, target_xy: np.ndarray) -> float:
    if len(query_xy) == 0 or len(target_xy) == 0:
        return float("nan")
    tree = cKDTree(target_xy)
    d, _ = tree.query(query_xy, k=1)
    return float(np.mean(d))


def neighbor_sum(query_idx: np.ndarray, values: np.ndarray, nbr: np.ndarray) -> float:
    if query_idx.size == 0 or nbr.size == 0:
        return float("nan")
    return float(values[nbr[query_idx]].sum(axis=1).mean())


def largest_component(mask: np.ndarray, nbr: np.ndarray) -> np.ndarray:
    n = len(mask)
    seen = np.zeros(n, dtype=bool)
    best: list[int] = []
    for i in range(n):
        if not mask[i] or seen[i]:
            continue
        q = deque([i])
        seen[i] = True
        comp = [i]
        while q:
            u = q.popleft()
            for v in nbr[u]:
                if mask[v] and not seen[v]:
                    seen[v] = True
                    q.append(int(v))
                    comp.append(int(v))
        if len(comp) > len(best):
            best = comp
    out = np.zeros(n, dtype=bool)
    if best:
        out[np.array(best, dtype=int)] = True
    return out


def hop_depth(tumor: np.ndarray, nbr: np.ndarray) -> np.ndarray:
    n = len(tumor)
    depth = np.full(n, -1, dtype=int)
    q: deque[int] = deque()
    for i in range(n):
        if not tumor[i]:
            depth[i] = 0
            q.append(i)
    while q:
        i = q.popleft()
        for j in nbr[i]:
            if depth[j] == -1:
                depth[j] = depth[i] + 1
                q.append(int(j))
    return depth


def study_area_um2(xy: np.ndarray) -> float:
    if len(xy) < 3:
        if len(xy) == 0:
            return float("nan")
        return float(max(1.0, np.ptp(xy[:, 0]) * np.ptp(xy[:, 1])))
    try:
        hull = ConvexHull(xy)
        return float(hull.volume)
    except QhullError:
        return float(max(1.0, np.ptp(xy[:, 0]) * np.ptp(xy[:, 1])))


def cross_K(xy_a: np.ndarray, xy_b: np.ndarray, area: float, radii: np.ndarray) -> np.ndarray:
    n1, n2 = len(xy_a), len(xy_b)
    out = np.full(len(radii), np.nan)
    if n1 < 2 or n2 < 2 or not np.isfinite(area) or area <= 0:
        return out
    tree = cKDTree(xy_b)
    same = np.array_equal(xy_a, xy_b)
    for i, r in enumerate(radii):
        neigh = tree.query_ball_point(xy_a, r)
        tot = sum(len(ix) for ix in neigh)
        if same:
            tot -= n1
        out[i] = area * tot / (n1 * n2)
    return out


def L_of_K(K: np.ndarray, radii: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore"):
        return np.sqrt(np.maximum(K, 0.0) / math.pi) - radii


def high_mask(x: np.ndarray) -> np.ndarray:
    q = np.quantile(x, 0.75)
    if q == 0:
        return x > 0
    return x >= q


def analyze_section(meta: dict[str, str], section_dir: Path, n_perm: int, rng: np.random.Generator) -> dict:
    rec: dict = {
        "gsm": meta["gsm"],
        "label": meta["label"],
        "patient": meta["patient"],
        "histology": meta["histology"],
        "class": meta["class"],
        "status": "ok",
        "skip_reason": "",
    }
    data = load_section(section_dir)
    genes = data["genes"]
    present = {g: g in genes for g in GENES}
    for g, p in present.items():
        rec[f"has_{g}"] = bool(p)
    if not present["CLDN4"] and not present["CD8A"]:
        rec["status"] = "skipped"
        rec["skip_reason"] = "CLDN4 and CD8A absent"
        return rec
    if not present["CLDN4"] or not present["CD8A"]:
        rec["status"] = "skipped"
        rec["skip_reason"] = "CLDN4 or CD8A absent"
        return rec

    lib = np.asarray(data["matrix"].sum(axis=0)).ravel().astype(float)
    rec["n_spots"] = int(len(lib))
    rec["median_umi"] = float(np.median(lib))
    counts = extract_counts(data["matrix"], genes, GENES)
    expr = log_norm(counts, lib)
    cldn4 = expr["CLDN4"]
    cd8a = expr["CD8A"]
    nkg7 = expr.get("NKG7")
    krt8 = expr.get("KRT8")
    epcam = expr.get("EPCAM")
    cd8_nkg7 = cd8a + nkg7 if nkg7 is not None else None

    xy = data["xy_um"]
    nbr = knn_indices(xy, K_NBR)
    rec["median_nn_um"] = float(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1]))

    epi_parts = []
    if epcam is not None:
        epi_parts.append(zscore(epcam))
    if krt8 is not None:
        epi_parts.append(zscore(krt8))
    rec["epi_defined"] = bool(epi_parts)
    if epi_parts:
        epi_score = np.mean(np.vstack(epi_parts), axis=0)
        epi = epi_score >= np.quantile(epi_score, 0.75)
    else:
        epi_score = np.zeros(len(cldn4))
        epi = np.ones(len(cldn4), dtype=bool)
    rec["n_epi"] = int(epi.sum())
    rec["median_CLDN4_all"] = float(np.median(cldn4))
    rec["median_CLDN4_epi"] = float(np.median(cldn4[epi])) if epi.any() else float("nan")
    rec["median_CD8A_all"] = float(np.median(cd8a))
    rec["median_KRT8_all"] = float(np.median(krt8)) if krt8 is not None else float("nan")

    r, p = spearman(cldn4, cd8a)
    rec["rho_CLDN4_CD8A"] = r
    rec["p_CLDN4_CD8A"] = p
    if epi.sum() >= 10:
        r_e, p_e = spearman(cldn4[epi], cd8a[epi])
    else:
        r_e, p_e = float("nan"), float("nan")
    rec["rho_CLDN4_CD8A_epi"] = r_e
    rec["p_CLDN4_CD8A_epi"] = p_e
    if cd8_nkg7 is not None:
        r2, p2 = spearman(cldn4, cd8_nkg7)
        rec["rho_CLDN4_CD8A_NKG7"] = r2
        rec["p_CLDN4_CD8A_NKG7"] = p2
    else:
        rec["rho_CLDN4_CD8A_NKG7"] = float("nan")
        rec["p_CLDN4_CD8A_NKG7"] = float("nan")

    rec["I_CLDN4"] = moran_i(cldn4, nbr)
    rec["I_CD8A"] = moran_i(cd8a, nbr)
    rec["I_biv_CLDN4_CD8A"] = bivariate_moran(cldn4, cd8a, nbr)
    rec["LeeL_CLDN4_CD8A"] = lee_L(cldn4, cd8a, nbr)
    rec["rho_CLDN4_lagCD8A"] = spearman(cldn4, apply_W(cd8a, nbr))[0]

    if krt8 is not None:
        r_c = residualize(cldn4, krt8)
        r_8 = residualize(cd8a, krt8)
        rec["rho_resid_CLDN4_CD8A"] = spearman(r_c, r_8)[0]
        rec["p_resid_CLDN4_CD8A"] = spearman(r_c, r_8)[1]
        lag_r8 = apply_W(r_8, nbr)
        rec["rho_resid_CLDN4_lagCD8A"] = spearman(r_c, lag_r8)[0]
        rec["p_resid_CLDN4_lagCD8A"] = spearman(r_c, lag_r8)[1]
        rec["I_biv_resid"] = bivariate_moran(r_c, r_8, nbr)
        rec["LeeL_resid"] = lee_L(r_c, r_8, nbr)
    else:
        r_c = cldn4.copy()
        r_8 = cd8a.copy()
        rec["rho_resid_CLDN4_CD8A"] = float("nan")
        rec["p_resid_CLDN4_CD8A"] = float("nan")
        rec["rho_resid_CLDN4_lagCD8A"] = float("nan")
        rec["p_resid_CLDN4_lagCD8A"] = float("nan")
        rec["I_biv_resid"] = float("nan")
        rec["LeeL_resid"] = float("nan")

    gi_c = gi_star(cldn4, nbr)
    gi_8 = gi_star(cd8a, nbr)
    hot_c = gi_c > GI_Z
    cold_8 = gi_8 < -GI_Z
    overlap = hot_c & cold_8
    rec["n_Gi_CLDN4_hot"] = int(hot_c.sum())
    rec["n_Gi_CD8A_cold"] = int(cold_8.sum())
    rec["n_Gi_overlap_hotCLDN4_coldCD8A"] = int(overlap.sum())
    rec["frac_CLDN4hot_that_are_CD8cold"] = (
        float(overlap.sum() / hot_c.sum()) if hot_c.any() else float("nan")
    )
    union = hot_c | cold_8
    rec["jaccard_CLDN4hot_CD8cold"] = float(overlap.sum() / union.sum()) if union.any() else float("nan")
    rec["expected_overlap_indep"] = float(hot_c.mean() * cold_8.mean() * len(cldn4))

    epi_idx = np.flatnonzero(epi)
    if epi_idx.size >= MIN_Q_SPOTS * 2:
        q_cldn = cldn4[epi]
        q1 = q_cldn <= np.quantile(q_cldn, 0.25)
        q4 = q_cldn >= np.quantile(q_cldn, 0.75)
        q1_idx = epi_idx[q1]
        q4_idx = epi_idx[q4]
    else:
        q1_idx = np.array([], dtype=int)
        q4_idx = np.array([], dtype=int)
    rec["n_CLDN4_Q1_epi"] = int(q1_idx.size)
    rec["n_CLDN4_Q4_epi"] = int(q4_idx.size)

    cd8_hi = high_mask(cd8a)
    rec["n_CD8A_high"] = int(cd8_hi.sum())
    d_q4 = nearest_mean_dist(xy[q4_idx], xy[cd8_hi]) if q4_idx.size and cd8_hi.any() else float("nan")
    d_q1 = nearest_mean_dist(xy[q1_idx], xy[cd8_hi]) if q1_idx.size and cd8_hi.any() else float("nan")
    rec["mean_nn_um_Q4_to_CD8high"] = d_q4
    rec["mean_nn_um_Q1_to_CD8high"] = d_q1
    rec["delta_nn_um_Q4_minus_Q1"] = d_q4 - d_q1 if np.isfinite(d_q4) and np.isfinite(d_q1) else float("nan")
    rec["mean_kNN_CD8A_Q4"] = neighbor_sum(q4_idx, cd8a, nbr)
    rec["mean_kNN_CD8A_Q1"] = neighbor_sum(q1_idx, cd8a, nbr)
    rec["delta_kNN_CD8A_Q4_minus_Q1"] = (
        rec["mean_kNN_CD8A_Q4"] - rec["mean_kNN_CD8A_Q1"]
        if np.isfinite(rec["mean_kNN_CD8A_Q4"]) and np.isfinite(rec["mean_kNN_CD8A_Q1"])
        else float("nan")
    )

    area = study_area_um2(xy)
    spacing = rec["median_nn_um"]
    radii = np.array([1, 2, 3, 4, 5, 6], dtype=float) * spacing
    rec["K_radii_um"] = radii.tolist()
    if q4_idx.size >= 8 and cd8_hi.sum() >= 8:
        K = cross_K(xy[q4_idx], xy[cd8_hi], area, radii)
        L = L_of_K(K, radii)
        rec["crossK_Q4_CD8high"] = [float(v) for v in K]
        rec["crossL_Q4_CD8high"] = [float(v) for v in L]
        rec["crossL_at_2nn"] = float(L[1]) if len(L) > 1 else float("nan")
    else:
        rec["crossK_Q4_CD8high"] = []
        rec["crossL_Q4_CD8high"] = []
        rec["crossL_at_2nn"] = float("nan")

    tumor = largest_component(epi, nbr) if rec["epi_defined"] else np.zeros(len(cldn4), dtype=bool)
    rec["n_tumor_domain"] = int(tumor.sum())
    rec["tumor_domain_defined"] = bool(tumor.sum() >= MIN_TUMOR_CC)
    if rec["tumor_domain_defined"]:
        depth = hop_depth(tumor, nbr)
        t_idx = np.flatnonzero(tumor)
        t_depth = depth[t_idx].astype(float)
        rec["median_tumor_hop_depth"] = float(np.median(t_depth))
        rec["max_tumor_hop_depth"] = int(t_depth.max())
        rec["rho_CD8A_vs_hopdepth_tumor"] = spearman(cd8a[t_idx], t_depth)[0]
        rec["p_CD8A_vs_hopdepth_tumor"] = spearman(cd8a[t_idx], t_depth)[1]
        margin = tumor & (depth == 1)
        core_hop = max(2, int(np.quantile(t_depth, 0.75))) if t_depth.size else 2
        cldn4_t = np.zeros(len(cldn4), dtype=bool)
        cldn4_t[t_idx] = cldn4[t_idx] >= np.quantile(cldn4[t_idx], 0.75)
        core = tumor & (depth >= core_hop) & cldn4_t
        rec["n_margin"] = int(margin.sum())
        rec["n_CLDN4_high_core"] = int(core.sum())
        rec["mean_CD8A_margin"] = float(cd8a[margin].mean()) if margin.any() else float("nan")
        rec["mean_CD8A_CLDN4_core"] = float(cd8a[core].mean()) if core.any() else float("nan")
        rec["delta_CD8A_core_minus_margin"] = (
            rec["mean_CD8A_CLDN4_core"] - rec["mean_CD8A_margin"]
            if np.isfinite(rec["mean_CD8A_CLDN4_core"]) and np.isfinite(rec["mean_CD8A_margin"])
            else float("nan")
        )
        depth_bins = {}
        for d in sorted(set(t_depth.astype(int).tolist())):
            m = tumor & (depth == d)
            depth_bins[str(int(d))] = float(cd8a[m].mean()) if m.any() else float("nan")
        rec["CD8A_by_hop_depth"] = depth_bins
        if core.any() and cd8_hi.any():
            rec["mean_nn_um_CD8high_to_CLDN4core"] = nearest_mean_dist(xy[cd8_hi], xy[core])
        else:
            rec["mean_nn_um_CD8high_to_CLDN4core"] = float("nan")
    else:
        rec["median_tumor_hop_depth"] = float("nan")
        rec["max_tumor_hop_depth"] = 0
        rec["rho_CD8A_vs_hopdepth_tumor"] = float("nan")
        rec["p_CD8A_vs_hopdepth_tumor"] = float("nan")
        rec["n_margin"] = 0
        rec["n_CLDN4_high_core"] = 0
        rec["mean_CD8A_margin"] = float("nan")
        rec["mean_CD8A_CLDN4_core"] = float("nan")
        rec["delta_CD8A_core_minus_margin"] = float("nan")
        rec["CD8A_by_hop_depth"] = {}
        rec["mean_nn_um_CD8high_to_CLDN4core"] = float("nan")
        depth = None
        core = np.zeros(len(cldn4), dtype=bool)
        margin = np.zeros(len(cldn4), dtype=bool)

    # Permute CD8A among epithelial-like spots.
    perm_keys = [
        "rho_CLDN4_CD8A_epi",
        "I_biv_CLDN4_CD8A",
        "LeeL_CLDN4_CD8A",
        "rho_resid_CLDN4_lagCD8A",
        "delta_nn_um_Q4_minus_Q1",
        "delta_kNN_CD8A_Q4_minus_Q1",
        "n_Gi_overlap_hotCLDN4_coldCD8A",
        "crossL_at_2nn",
        "delta_CD8A_core_minus_margin",
        "rho_CD8A_vs_hopdepth_tumor",
    ]
    nulls = {k: [] for k in perm_keys}
    if epi_idx.size >= 20 and n_perm > 0:
        cd8_epi_vals = cd8a[epi_idx].copy()
        r_c_fixed = r_c
        krt8_ok = krt8 is not None
        for _ in range(n_perm):
            cd8_p = cd8a.copy()
            cd8_p[epi_idx] = rng.permutation(cd8_epi_vals)
            if epi_idx.size >= 10:
                nulls["rho_CLDN4_CD8A_epi"].append(spearman(cldn4[epi], cd8_p[epi])[0])
            nulls["I_biv_CLDN4_CD8A"].append(bivariate_moran(cldn4, cd8_p, nbr))
            nulls["LeeL_CLDN4_CD8A"].append(lee_L(cldn4, cd8_p, nbr))
            if krt8_ok:
                r8_p = residualize(cd8_p, krt8)
                nulls["rho_resid_CLDN4_lagCD8A"].append(spearman(r_c_fixed, apply_W(r8_p, nbr))[0])
            cd8_hi_p = high_mask(cd8_p)
            if q4_idx.size and q1_idx.size and cd8_hi_p.any():
                d4 = nearest_mean_dist(xy[q4_idx], xy[cd8_hi_p])
                d1 = nearest_mean_dist(xy[q1_idx], xy[cd8_hi_p])
                nulls["delta_nn_um_Q4_minus_Q1"].append(d4 - d1)
            if q4_idx.size and q1_idx.size:
                nulls["delta_kNN_CD8A_Q4_minus_Q1"].append(
                    neighbor_sum(q4_idx, cd8_p, nbr) - neighbor_sum(q1_idx, cd8_p, nbr)
                )
            gi8_p = gi_star(cd8_p, nbr)
            nulls["n_Gi_overlap_hotCLDN4_coldCD8A"].append(int((hot_c & (gi8_p < -GI_Z)).sum()))
            if q4_idx.size >= 8 and cd8_hi_p.sum() >= 8:
                Kp = cross_K(xy[q4_idx], xy[cd8_hi_p], area, radii)
                Lp = L_of_K(Kp, radii)
                nulls["crossL_at_2nn"].append(float(Lp[1]))
            if rec["tumor_domain_defined"] and depth is not None:
                t_idx = np.flatnonzero(tumor)
                nulls["rho_CD8A_vs_hopdepth_tumor"].append(spearman(cd8_p[t_idx], depth[t_idx].astype(float))[0])
                if core.any() and margin.any():
                    nulls["delta_CD8A_core_minus_margin"].append(cd8_p[core].mean() - cd8_p[margin].mean())

    alt = {
        "rho_CLDN4_CD8A_epi": "less",
        "I_biv_CLDN4_CD8A": "less",
        "LeeL_CLDN4_CD8A": "less",
        "rho_resid_CLDN4_lagCD8A": "less",
        "delta_nn_um_Q4_minus_Q1": "greater",
        "delta_kNN_CD8A_Q4_minus_Q1": "less",
        "n_Gi_overlap_hotCLDN4_coldCD8A": "greater",
        "crossL_at_2nn": "less",
        "delta_CD8A_core_minus_margin": "less",
        "rho_CD8A_vs_hopdepth_tumor": "less",
    }
    rec["n_perm"] = int(n_perm)
    rec["perm_nulls"] = {k: [float(v) for v in vs if np.isfinite(v)] for k, vs in nulls.items()}
    for k, how in alt.items():
        rec[f"p_perm_{k}"] = empirical_p(rec.get(k, float("nan")), np.asarray(nulls[k], float), how)
        rec[f"p_perm_two_{k}"] = empirical_p(rec.get(k, float("nan")), np.asarray(nulls[k], float), "two-sided")

    rec["_plot"] = {
        "xy": xy,
        "cldn4": cldn4,
        "cd8a": cd8a,
        "r_c": r_c,
        "lag_r8": apply_W(r_8, nbr),
        "gi_c": gi_c,
        "gi_8": gi_8,
        "hot_c": hot_c,
        "cold_8": cold_8,
        "overlap": overlap,
        "epi": epi,
        "q4_idx": q4_idx,
        "q1_idx": q1_idx,
        "cd8_hi": cd8_hi,
        "tumor": tumor,
        "core": core,
        "margin": margin,
        "radii": radii,
        "K": np.asarray(rec["crossK_Q4_CD8high"], float) if rec["crossK_Q4_CD8high"] else np.array([]),
        "L": np.asarray(rec["crossL_Q4_CD8high"], float) if rec["crossL_Q4_CD8high"] else np.array([]),
        "L_null": np.asarray(rec["perm_nulls"].get("crossL_at_2nn", []), float),
    }
    return rec


def drop_plot(rec: dict) -> dict:
    rec = dict(rec)
    rec.pop("_plot", None)
    rec.pop("perm_nulls", None)
    return rec


def fmt(x, nd=3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    return f"{float(x):.{nd}g}"


def wilcoxon_vs_zero(vals: np.ndarray) -> float:
    vals = vals[np.isfinite(vals)]
    if vals.size < 6 or np.allclose(vals, 0):
        return float("nan")
    try:
        return float(stats.wilcoxon(vals, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return float("nan")


def median_iqr(vals: np.ndarray) -> tuple[float, float, float]:
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan"), float("nan"), float("nan")
    return float(np.median(vals)), float(np.quantile(vals, 0.25)), float(np.quantile(vals, 0.75))


def stouffer_from_p(ps: np.ndarray, signs: np.ndarray) -> float:
    ps = np.clip(ps, 1e-12, 1 - 1e-12)
    z = -stats.norm.ppf(ps / 2.0) * np.sign(signs)
    z = z[np.isfinite(z)]
    if z.size == 0:
        return float("nan")
    zc = float(z.sum() / math.sqrt(len(z)))
    return float(2 * stats.norm.sf(abs(zc)))


def plot_expression_map(plot: dict, title: str, path: Path) -> None:
    xy, c, d = plot["xy"], plot["cldn4"], plot["cd8a"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    for ax, val, name, cmap in (
        (axes[0], c, "CLDN4 (log-norm)", "magma"),
        (axes[1], d, "CD8A (log-norm)", "viridis"),
    ):
        sc = ax.scatter(xy[:, 0], -xy[:, 1], c=val, s=6, cmap=cmap, linewidths=0)
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(name)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("x (µm)")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_gi_star(plot: dict, title: str, path: Path) -> None:
    xy = plot["xy"]
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.0))
    for ax, val, name, cmap, vmin, vmax in (
        (axes[0], plot["gi_c"], "Gi* CLDN4", "RdBu_r", -4, 4),
        (axes[1], plot["gi_8"], "Gi* CD8A", "RdBu_r", -4, 4),
    ):
        sc = ax.scatter(xy[:, 0], -xy[:, 1], c=val, s=6, cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0)
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(name)
        ax.set_aspect("equal")
        ax.axis("off")
    col = np.zeros(len(xy), dtype=int)
    col[plot["hot_c"]] = 1
    col[plot["cold_8"]] = 2
    col[plot["overlap"]] = 3
    cmap = matplotlib.colors.ListedColormap(["#d9d9d9", "#d62728", "#1f77b4", "#6a3d9a"])
    axes[2].scatter(xy[:, 0], -xy[:, 1], c=col, s=6, cmap=cmap, vmin=0, vmax=3, linewidths=0)
    axes[2].set_title("hot CLDN4 / cold CD8A / overlap")
    axes[2].set_aspect("equal")
    axes[2].axis("off")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_moran_scatter(plot: dict, title: str, path: Path) -> None:
    x, y = plot["r_c"], plot["lag_r8"]
    fig, ax = plt.subplots(figsize=(5.0, 4.6))
    ax.scatter(x, y, s=8, alpha=0.35, c="#222222", linewidths=0)
    if np.unique(x).size > 1:
        b, a = np.polyfit(x, y, 1)
        xs = np.linspace(x.min(), x.max(), 50)
        ax.plot(xs, a + b * xs, color="#d62728", lw=1.5)
    ax.axhline(0, color="#888", lw=0.6)
    ax.axvline(0, color="#888", lw=0.6)
    ax.set_xlabel("CLDN4 residual on KRT8")
    ax.set_ylabel("Spatial lag of CD8A residual on KRT8")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_kfunction(plot: dict, title: str, path: Path) -> None:
    radii, L = plot["radii"], plot["L"]
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    if L.size:
        ax.plot(radii, L, marker="o", color="#222", label="L₁₂(r) Q4 CLDN4 vs CD8A-high")
        ax.axhline(0, color="#888", lw=0.8, label="CSR reference L=0")
    else:
        ax.text(0.5, 0.5, "cross-K not computed", ha="center", va="center", transform=ax.transAxes)
    ax.set_xlabel("r (µm)")
    ax.set_ylabel("L₁₂(r) = √(K/π) − r")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def forest_plot(df: pd.DataFrame, col: str, title: str, path: Path, vline: float = 0.0) -> None:
    sub = df.dropna(subset=[col]).sort_values(["class", "label"])
    fig, ax = plt.subplots(figsize=(8.5, max(3.5, 0.22 * len(sub) + 1.2)))
    y = np.arange(len(sub))
    colors = {"invasive": "#d62728", "precursor": "#1f77b4", "normal": "#7f7f7f"}
    ax.scatter(sub[col], y, c=sub["class"].map(colors), s=18)
    ax.axvline(vline, color="#444", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(sub["label"].tolist(), fontsize=7)
    ax.set_xlabel(col)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def box_by_histology(df: pd.DataFrame, col: str, title: str, path: Path) -> None:
    order = [h for h in ("Normal", "AAH", "AIS", "MIA", "LUAD") if h in set(df["histology"])]
    data = [df.loc[df["histology"] == h, col].dropna().to_numpy() for h in order]
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.boxplot(data, tick_labels=order, showfliers=False)
    for i, vals in enumerate(data, start=1):
        ax.scatter(np.random.default_rng(1).normal(i, 0.06, size=len(vals)), vals, s=12, alpha=0.7, c="#333")
    ax.axhline(0, color="#888", lw=0.6)
    ax.set_title(title)
    ax.set_ylabel(col)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def summarize_block(df: pd.DataFrame, name: str) -> dict:
    out = {"set": name, "n_sections": int(len(df))}
    metrics = [
        "rho_CLDN4_CD8A",
        "rho_CLDN4_CD8A_epi",
        "rho_CLDN4_CD8A_NKG7",
        "rho_resid_CLDN4_CD8A",
        "rho_resid_CLDN4_lagCD8A",
        "I_biv_CLDN4_CD8A",
        "LeeL_CLDN4_CD8A",
        "I_biv_resid",
        "LeeL_resid",
        "delta_nn_um_Q4_minus_Q1",
        "delta_kNN_CD8A_Q4_minus_Q1",
        "n_Gi_overlap_hotCLDN4_coldCD8A",
        "gi_overlap_minus_expected",
        "frac_CLDN4hot_that_are_CD8cold",
        "crossL_at_2nn",
        "delta_CD8A_core_minus_margin",
        "rho_CD8A_vs_hopdepth_tumor",
    ]
    for m in metrics:
        med, q1, q3 = median_iqr(df[m].to_numpy(float) if m in df.columns else np.array([]))
        out[f"{m}_median"] = med
        out[f"{m}_q1"] = q1
        out[f"{m}_q3"] = q3
        out[f"{m}_wilcoxon_p"] = wilcoxon_vs_zero(df[m].to_numpy(float) if m in df.columns else np.array([]))
        if m in df.columns:
            vals = df[m].to_numpy(float)
            out[f"{m}_n"] = int(np.isfinite(vals).sum())
            out[f"{m}_n_neg"] = int(np.sum(vals[np.isfinite(vals)] < 0))
            out[f"{m}_n_pos"] = int(np.sum(vals[np.isfinite(vals)] > 0))
        pcol = f"p_perm_{m}"
        if pcol in df.columns:
            medp, _, _ = median_iqr(df[pcol].to_numpy(float))
            out[f"{pcol}_median"] = medp
            out[f"{pcol}_n_lt_0.05"] = int(np.sum(df[pcol].to_numpy(float) < 0.05))
    return out


def write_results_md(path: Path, ok: pd.DataFrame, meta_rows: list[dict], gene_note: dict) -> None:
    prim = ok[ok["class"] == "invasive"].copy()
    prec = ok[ok["class"] == "precursor"].copy()
    use_primary = prim if len(prim) else ok
    primary_name = "invasive LUAD" if len(prim) else "all analyzed sections"

    def block(df: pd.DataFrame, title: str) -> list[str]:
        if df.empty:
            return [f"### {title}", "", "No sections in this cut.", ""]
        s = summarize_block(df, title)
        lines = [
            f"### {title}",
            "",
            f"n sections = {s['n_sections']}.",
            "",
            "| statistic | n | median [Q1, Q3] | n<0 | n>0 | Wilcoxon p vs 0 | median epi-perm p (exclusion dir.) | n perm p<0.05 |",
            "|---|---:|---|---:|---:|---:|---:|---:|",
        ]
        rows = [
            ("Spearman CLDN4 vs CD8A (all spots)", "rho_CLDN4_CD8A"),
            ("Spearman CLDN4 vs CD8A (epi-like spots)", "rho_CLDN4_CD8A_epi"),
            ("Spearman CLDN4 vs CD8A+NKG7", "rho_CLDN4_CD8A_NKG7"),
            ("Spearman after KRT8 residual", "rho_resid_CLDN4_CD8A"),
            ("CLDN4 residual vs spatial lag of CD8A residual", "rho_resid_CLDN4_lagCD8A"),
            ("Bivariate Moran I(CLDN4, CD8A)", "I_biv_CLDN4_CD8A"),
            ("Lee's L(CLDN4, CD8A)", "LeeL_CLDN4_CD8A"),
            ("Δ nearest-µm (CLDN4 Q4 − Q1) to CD8A-high", "delta_nn_um_Q4_minus_Q1"),
            ("Δ kNN CD8A sum (Q4 − Q1)", "delta_kNN_CD8A_Q4_minus_Q1"),
            ("Gi* overlap n (CLDN4 hot ∩ CD8A cold)", "n_Gi_overlap_hotCLDN4_coldCD8A"),
            ("Gi* overlap minus independence expectation", "gi_overlap_minus_expected"),
            ("Cross-L at 2× median NN", "crossL_at_2nn"),
            ("Mean CD8A in CLDN4-high core − margin", "delta_CD8A_core_minus_margin"),
            ("Spearman CD8A vs hop-depth (tumor domain)", "rho_CD8A_vs_hopdepth_tumor"),
        ]
        for lab, key in rows:
            if f"{key}_n" not in s:
                continue
            med = fmt(s.get(f"{key}_median"))
            q1 = fmt(s.get(f"{key}_q1"))
            q3 = fmt(s.get(f"{key}_q3"))
            lines.append(
                f"| {lab} | {s.get(f'{key}_n', 0)} | {med} [{q1}, {q3}] | "
                f"{s.get(f'{key}_n_neg', 0)} | {s.get(f'{key}_n_pos', 0)} | "
                f"{fmt(s.get(f'{key}_wilcoxon_p'))} | {fmt(s.get(f'p_perm_{key}_median'))} | "
                f"{s.get(f'p_perm_{key}_n_lt_0.05', 'NA')} |"
            )
        lines.append("")
        return lines

    cldn4_by_h = (
        ok.groupby("histology")["median_CLDN4_epi"].median().reindex(["Normal", "AAH", "AIS", "MIA", "LUAD"])
    )
    lines = [
        "# GSE307534 Visium: CLDN4 vs CD8 spatial statistics",
        "",
        "Public dataset only (GEO [GSE307534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307534)): "
        "Visium CytAssist of normal lung, AAH, AIS, MIA, and invasive LUAD. "
        "Analysis is **CLDN4-only**. No private gene sets were used.",
        "",
        "## Question",
        "",
        "Do CLDN4-high tumor/epithelial spots spatially exclude CD8?",
        "",
        "## Data and gene check",
        "",
        f"- Sections with processed spots downloaded: {gene_note['n_downloaded']}.",
        f"- Sections analyzed (CLDN4 and CD8A present): {gene_note['n_ok']}.",
        f"- Sections skipped (both or one of CLDN4/CD8A absent, or load error): {gene_note['n_skip']}.",
        f"- CLDN4 present in {gene_note['n_CLDN4']} / {gene_note['n_downloaded']} downloaded sections.",
        f"- CD8A present in {gene_note['n_CD8A']} / {gene_note['n_downloaded']} downloaded sections.",
        f"- NKG7 present in {gene_note['n_NKG7']} / {gene_note['n_downloaded']}; "
        f"GNLY present in {gene_note['n_GNLY']} / {gene_note['n_downloaded']}.",
        f"- EPCAM present in {gene_note['n_EPCAM']} / {gene_note['n_downloaded']}; "
        f"KRT8 present in {gene_note['n_KRT8']} / {gene_note['n_downloaded']}.",
        "",
        "## Methods (computed, not assumed)",
        "",
        "- Filtered Space Ranger matrices and `tissue_positions.csv` from each GSM tar; H&E images were not used.",
        "- Expression: log1p(10⁴ × counts / spot UMI). Duplicate gene symbols were summed.",
        "- Epithelial-like spots: top quartile of the mean z-score of available EPCAM and KRT8.",
        "- CLDN4-high / low among epithelial-like spots: top vs bottom quartile (Q4 vs Q1).",
        "- CD8A-high: spots at or above the 75th percentile of CD8A (or CD8A>0 if that percentile is 0).",
        "- Spatial weights: row-standardized 6-nearest neighbors in micron coordinates "
        "(55 µm / `spot_diameter_fullres`).",
        "- Bivariate Moran I(x,y) = z_x' W z_y / z_x' z_x with mean-centered (not variance-standardized) x,y; "
        "|I| and |L| can exceed 1 when CLDN4 and CD8A have different variances. Lee's L as in Lee (2001) with row-standardized W.",
        "- Residual spatial lag: OLS-residualize CLDN4 and CD8A on KRT8, then Spearman of residual CLDN4 vs W·residual CD8A.",
        "- Getis-Ord Gi* with binary kNN including self; hotspot/coldspot |z|>1.96. Overlap = CLDN4 hot ∩ CD8A cold.",
        "- Cross-K / L₁₂(r) between CLDN4 Q4 epithelial-like points and CD8A-high points; L(r)=√(K/π)−r; "
        "summary radius = 2 × median nearest-neighbor distance.",
        "- Q4 vs Q1: mean distance (µm) to nearest CD8A-high spot, and mean CD8A sum in 6-NN.",
        "- Tumor domain (when defined): largest 6-NN connected component of epithelial-like spots with ≥40 spots. "
        "Hop-depth from non-domain spots; margin = depth 1; CLDN4-high core = CLDN4 Q4 among domain spots "
        "and depth ≥ max(2, domain-depth Q3).",
        f"- Empirical p: {int(ok['n_perm'].median()) if len(ok) else 'NA'} permutations of CD8A among epithelial-like spots "
        "(non-epithelial CD8A left in place). One-sided p is in the exclusion direction "
        "(negative association / larger Q4 distance / more Gi* overlap / lower core CD8A).",
        "- Primary cut is invasive LUAD if those sections are available; precursor vs invasive is a secondary cut.",
        "- No private 8-KL or other private signatures.",
        "",
        "## Histology check (whether precursors are normal-like for CLDN4)",
        "",
        "Median across sections of per-section median CLDN4 in epithelial-like spots:",
        "",
    ]
    for h in ("Normal", "AAH", "AIS", "MIA", "LUAD"):
        v = cldn4_by_h.get(h, np.nan)
        n_h = int((ok["histology"] == h).sum())
        lines.append(f"- {h}: {fmt(v)} (n={n_h} sections)")
    lines += [
        "",
        f"Primary reporting set: **{primary_name}** (n={len(use_primary)}).",
        "",
        "## Meta results",
        "",
    ]
    lines += block(use_primary, f"Primary: {primary_name}")
    lines += block(prec, "Secondary: precursor (AAH / AIS / MIA)")
    lines += block(ok, "All analyzed sections")
    if len(use_primary):
        s = summarize_block(use_primary, "primary")
        lines += [
            "## Observed direction on the primary set",
            "",
            f"- n={s['n_sections']} sections.",
            f"- Median Spearman CLDN4 vs CD8A = {fmt(s.get('rho_CLDN4_CD8A_median'))} "
            f"({s.get('rho_CLDN4_CD8A_n_neg', 0)} negative, {s.get('rho_CLDN4_CD8A_n_pos', 0)} positive; "
            f"Wilcoxon p={fmt(s.get('rho_CLDN4_CD8A_wilcoxon_p'))}).",
            f"- Median bivariate Moran I = {fmt(s.get('I_biv_CLDN4_CD8A_median'))}; "
            f"median Lee's L = {fmt(s.get('LeeL_CLDN4_CD8A_median'))}.",
            f"- Median residual CLDN4 vs lag residual CD8A ρ = {fmt(s.get('rho_resid_CLDN4_lagCD8A_median'))}.",
            f"- Median Δ nearest-µm (Q4 − Q1 to CD8A-high) = {fmt(s.get('delta_nn_um_Q4_minus_Q1_median'))} "
            f"(positive = Q4 farther from CD8A-high).",
            f"- Median Gi* overlap count (CLDN4 hot ∩ CD8A cold) = {fmt(s.get('n_Gi_overlap_hotCLDN4_coldCD8A_median'))}.",
            f"- Median cross-L at 2× NN = {fmt(s.get('crossL_at_2nn_median'))} "
            f"(negative = fewer CD8A-high near CLDN4 Q4 than CSR).",
            f"- Tumor-domain CD8A (core − margin) median = {fmt(s.get('delta_CD8A_core_minus_margin_median'))}.",
            "",
        ]
    lines += [
        "## How to read the exclusion-direction numbers",
        "",
        "- Negative Spearman / bivariate Moran / Lee's L / residual lag ρ: lower CD8A where CLDN4 is higher.",
        "- Positive Δ nearest-µm (Q4 − Q1): CLDN4-high epithelial-like spots are farther from CD8A-high spots than CLDN4-low.",
        "- Negative Δ kNN CD8A (Q4 − Q1): fewer CD8A UMIs in the neighborhood of CLDN4-high epithelial-like spots.",
        "- Positive Gi* overlap: CLDN4 hotspots coincide with CD8A coldspots.",
        "- Negative cross-L at 2× NN: fewer CD8A-high points near CLDN4 Q4 points than a CSR / labeling null.",
        "- Negative (core − margin) CD8A: lower CD8A in the CLDN4-high core than at the tumor-domain margin.",
        "",
        "These are observed statistics and permutation p-values. They are not a clinical claim.",
        "",
        "## Per-section table",
        "",
        "See `tables/section_stats.csv`. Figures:",
        "",
        "- `figures/summary_leeL_forest.png`, `figures/summary_biv_moran_forest.png`",
        "- `figures/summary_delta_nn_forest.png`, `figures/summary_gi_overlap.png`",
        "- `figures/summary_crossL_overlay.png`",
        "- `figures/summary_by_histology_leeL.png`, `figures/summary_by_histology_delta_nn.png`",
        "- `figures/maps/*_cldn4_cd8a.png` (H&E-free)",
        "- `figures/gi_star/*_gi_star.png`",
        "- `figures/moran/*_moran_lag.png`",
        "- `figures/kfunction/*_crossL.png`",
        "",
        "## Software",
        "",
        f"numpy {np.__version__}, pandas {pd.__version__}, scipy {__import__('scipy').__version__}, "
        f"matplotlib {matplotlib.__version__}.",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, default=Path("scripts/sample_manifest.tsv"))
    p.add_argument("--data-root", type=Path, default=Path("data/gse307534"))
    p.add_argument("--out-root", type=Path, default=Path("."))
    p.add_argument("--n-perm", type=int, default=N_PERM_DEFAULT)
    p.add_argument("--limit", type=int, default=0)
    args = p.parse_args()

    fig_root = args.out_root / "figures"
    for sub in ("maps", "gi_star", "moran", "kfunction"):
        (fig_root / sub).mkdir(parents=True, exist_ok=True)
    tab_root = args.out_root / "tables"
    tab_root.mkdir(parents=True, exist_ok=True)

    rows = load_manifest(args.manifest)
    if args.limit:
        rows = rows[: args.limit]
    rng = np.random.default_rng(RNG_SEED)
    records = []
    n_downloaded = 0
    for i, meta in enumerate(rows, start=1):
        section_dir = args.data_root / f"{meta['gsm']}_{meta['label']}"
        if not (section_dir / "filtered_feature_bc_matrix" / "matrix.mtx.gz").is_file():
            print(f"MISSING {meta['gsm']}_{meta['label']}", flush=True)
            records.append(
                {
                    "gsm": meta["gsm"],
                    "label": meta["label"],
                    "patient": meta["patient"],
                    "histology": meta["histology"],
                    "class": meta["class"],
                    "status": "skipped",
                    "skip_reason": "matrix not downloaded",
                }
            )
            continue
        n_downloaded += 1
        print(f"[{i}/{len(rows)}] {meta['gsm']}_{meta['label']}", flush=True)
        try:
            rec = analyze_section(meta, section_dir, args.n_perm, rng)
        except Exception as e:
            rec = {
                "gsm": meta["gsm"],
                "label": meta["label"],
                "patient": meta["patient"],
                "histology": meta["histology"],
                "class": meta["class"],
                "status": "skipped",
                "skip_reason": f"load/analysis error: {e}",
            }
            print(f"  ERROR {e}", flush=True)
            records.append(rec)
            continue
        plot = rec.get("_plot")
        if rec.get("status") == "ok" and plot is not None:
            stem = f"{meta['gsm']}_{meta['label']}"
            title = f"{meta['label']} ({meta['histology']})"
            plot_expression_map(plot, title, fig_root / "maps" / f"{stem}_cldn4_cd8a.png")
            plot_gi_star(plot, f"Gi* {title}", fig_root / "gi_star" / f"{stem}_gi_star.png")
            plot_moran_scatter(plot, f"Residual Moran scatter {title}", fig_root / "moran" / f"{stem}_moran_lag.png")
            plot_kfunction(plot, f"Cross-L {title}", fig_root / "kfunction" / f"{stem}_crossL.png")
        records.append(drop_plot(rec))
        print(
            f"  status={rec.get('status')} rho={rec.get('rho_CLDN4_CD8A')} "
            f"LeeL={rec.get('LeeL_CLDN4_CD8A')} dNN={rec.get('delta_nn_um_Q4_minus_Q1')}",
            flush=True,
        )

    slim = []
    for rec in records:
        row = {k: v for k, v in rec.items() if k not in {"crossK_Q4_CD8high", "crossL_Q4_CD8high", "K_radii_um", "CD8A_by_hop_depth"}}
        slim.append(row)
    df = pd.DataFrame(slim)
    df.to_csv(tab_root / "section_stats.csv", index=False)
    (tab_root / "section_stats.json").write_text(json.dumps(records, default=str, indent=2))

    ok = df[df["status"] == "ok"].copy()
    gene_note = {
        "n_downloaded": n_downloaded,
        "n_ok": int((df["status"] == "ok").sum()),
        "n_skip": int((df["status"] != "ok").sum()),
        "n_CLDN4": int(ok["has_CLDN4"].sum()) if len(ok) and "has_CLDN4" in ok.columns else 0,
        "n_CD8A": int(ok["has_CD8A"].sum()) if len(ok) and "has_CD8A" in ok.columns else 0,
        "n_NKG7": int(ok["has_NKG7"].sum()) if len(ok) and "has_NKG7" in ok.columns else 0,
        "n_GNLY": int(ok["has_GNLY"].sum()) if len(ok) and "has_GNLY" in ok.columns else 0,
        "n_EPCAM": int(ok["has_EPCAM"].sum()) if len(ok) and "has_EPCAM" in ok.columns else 0,
        "n_KRT8": int(ok["has_KRT8"].sum()) if len(ok) and "has_KRT8" in ok.columns else 0,
    }
    meta_rows = []
    if len(ok):
        for name, sub in (
            ("invasive", ok[ok["class"] == "invasive"]),
            ("precursor", ok[ok["class"] == "precursor"]),
            ("all", ok),
        ):
            meta_rows.append(summarize_block(sub, name))
        pd.DataFrame(meta_rows).to_csv(tab_root / "meta_stats.csv", index=False)
        forest_plot(ok, "LeeL_CLDN4_CD8A", "Lee's L(CLDN4, CD8A) per section", fig_root / "summary_leeL_forest.png")
        forest_plot(ok, "I_biv_CLDN4_CD8A", "Bivariate Moran I(CLDN4, CD8A)", fig_root / "summary_biv_moran_forest.png")
        forest_plot(
            ok,
            "delta_nn_um_Q4_minus_Q1",
            "Δ nearest distance (µm): CLDN4 Q4 − Q1 to CD8A-high",
            fig_root / "summary_delta_nn_forest.png",
        )
        forest_plot(
            ok,
            "n_Gi_overlap_hotCLDN4_coldCD8A",
            "Gi* overlap count: CLDN4 hot ∩ CD8A cold",
            fig_root / "summary_gi_overlap.png",
            vline=0.0,
        )
        box_by_histology(ok, "LeeL_CLDN4_CD8A", "Lee's L by histology", fig_root / "summary_by_histology_leeL.png")
        box_by_histology(
            ok,
            "delta_nn_um_Q4_minus_Q1",
            "Δ nearest-µm (Q4−Q1) by histology",
            fig_root / "summary_by_histology_delta_nn.png",
        )
        box_by_histology(ok, "median_CLDN4_epi", "Median CLDN4 in epi-like spots", fig_root / "summary_cldn4_epi_by_histology.png")
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        colors = {"invasive": "#d62728", "precursor": "#1f77b4", "normal": "#7f7f7f"}
        for rec in records:
            if rec.get("status") != "ok":
                continue
            L = rec.get("crossL_Q4_CD8high") or []
            radii = rec.get("K_radii_um") or []
            if len(L) and len(radii) == len(L):
                ax.plot(radii, L, color=colors.get(rec.get("class"), "#333"), alpha=0.35, lw=1.0)
        ax.axhline(0, color="#222", lw=0.8)
        ax.set_xlabel("r (µm)")
        ax.set_ylabel("L₁₂(r) = √(K/π) − r")
        ax.set_title("Cross-L: CLDN4 Q4 epi vs CD8A-high (red=LUAD, blue=precursor)")
        fig.tight_layout()
        fig.savefig(fig_root / "summary_crossL_overlay.png", dpi=140)
        plt.close(fig)
        if "rho_resid_CLDN4_lagCD8A" in ok.columns:
            forest_plot(
                ok,
                "rho_resid_CLDN4_lagCD8A",
                "Spearman: KRT8-residual CLDN4 vs lag of residual CD8A",
                fig_root / "summary_resid_lag_forest.png",
            )
    write_results_md(args.out_root / "RESULTS.md", ok, meta_rows, gene_note)
    print(f"wrote RESULTS.md n_ok={gene_note['n_ok']} n_skip={gene_note['n_skip']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
