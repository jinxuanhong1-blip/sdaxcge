#!/usr/bin/env python3
"""spaGCN-style spatial domains on CosMx NSCLC (figshare 25976224).

Histology-free sparse adaptation of SpaGCN (Hu et al., Nat Methods 2021):
Gaussian kNN graph convolution of log-normalized expression, then fixed-k
clustering. Labels are refined with a Giotto-style Potts HMRF (iterated
conditional modes; the model used by Giotto doHMRF).

CLDN4 is held out of the features that define domains. Domain CLDN4
enrichment is scored afterwards on author tumor cells and compared with
(i) Shannon diversity of compartment composition and (ii) a normalized
spatial mixing index (immune neighbors of tumor cells / sample immune
fraction). Cytotoxic neighbor counts at 50 µm and 100 µm link the domains
to the locked cell-level exclusion result.

Primary unit: 8 sections. Donor unit: 5 patients (Lung5 and Lung9 replicates
are averaged). No private cohort is read or written.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import read_h5ad
from scipy import sparse, stats
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_H5AD = ROOT / "data" / "cosmx" / "cosmx_human_nsclc_clustered.h5ad"
OUT = ROOT / "results" / "cosmx_spagcn_domains"

# CosMx global pixel coordinates in this object. Calibrated against the public
# Giotto spatial_locs (mm) on Lung5_Rep1 FOV 1: FOV span 0.97866 x 0.65142 mm
# vs 5437 x 3619 px → 0.180 µm/px. Median nearest-centroid distance is 9.0 µm.
UM_PER_PX = 0.18

K_SPATIAL = 10
N_HVG = 400
N_PCS = 30
K_DOMAINS = 8
BETA = 2.0
HMRF_ITERS = 25
MIN_AGREE = 0.55
BETA_IF_FRAGMENTED = 4.0
SEED = 20260921
MIN_DOMAIN_CELLS = 100
MIN_DOMAIN_TUMOR = 30
PURE_EXPECTED = 0.02

SAMPLE_CANON = {
    "LUAD-5 R1": "Lung5_Rep1",
    "LUAD-5 R2": "Lung5_Rep2",
    "LUAD-5 R3": "Lung5_Rep3",
    "LUSC-6": "Lung6",
    "LUAD-9 R1": "Lung9_Rep1",
    "LUAD-9 R2": "Lung9_Rep2",
    "LUAD-12": "Lung12",
    "LUAD-13": "Lung13",
}

IMMUNE = {
    "T CD4 naive",
    "T CD4 memory",
    "T CD8 naive",
    "T CD8 memory",
    "Treg",
    "NK",
    "B-cell",
    "plasmablast",
    "macrophage",
    "monocyte",
    "mDC",
    "pDC",
    "neutrophil",
    "mast",
}
CYTOTOXIC = {"T CD8 naive", "T CD8 memory", "NK"}
STROMAL = {"fibroblast", "endothelial"}
# CLDN4 is the held-out query. Epithelial controls are scored the same way
# and are not interpreted as a specificity claim.
CONTROL_GENES = ("EPCAM", "KRT8", "KRT18")
COMPARTMENTS = ("tumor", "immune", "stromal", "other")


def self_check() -> None:
    """Closed-form checks for Shannon and the normalized mixing ratio."""
    counts = np.array([50, 50, 0, 0], dtype=float)
    h = shannon_from_counts(counts)
    assert abs(h - 1.0) < 1e-9, h
    even = counts / counts.sum()
    # two equally likely classes among four bins: H = 1 bit
    assert abs(shannon_from_probs(even) - 1.0) < 1e-9
    # mixing ratio: observed 0.10 immune neighbors, sample immune fraction 0.20
    nm = normalized_mixing(0.10, 0.20)
    assert abs(nm - 0.5) < 1e-12
    assert np.isnan(normalized_mixing(0.10, 0.0))
    print("self_check ok", flush=True)


def shannon_from_probs(p: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if p.size == 0:
        return float("nan")
    return float(-(p * np.log2(p)).sum())


def shannon_from_counts(counts: np.ndarray) -> float:
    total = float(np.sum(counts))
    if total <= 0:
        return float("nan")
    return shannon_from_probs(np.asarray(counts, dtype=float) / total)


def normalized_mixing(obs_immune_neighbor_frac: float, sample_immune_frac: float) -> float:
    if not np.isfinite(sample_immune_frac) or sample_immune_frac <= 0:
        return float("nan")
    return float(obs_immune_neighbor_frac / sample_immune_frac)


def compartment_ids(cell_type: pd.Series) -> np.ndarray:
    ct = cell_type.astype(str)
    out = np.full(len(ct), 3, dtype=np.int8)  # other
    out[ct.str.startswith("tumor").to_numpy()] = 0
    out[ct.isin(IMMUNE).to_numpy()] = 1
    out[ct.isin(STROMAL).to_numpy()] = 2
    return out


def lognorm_selected(counts: sparse.spmatrix, gene_idx: np.ndarray) -> np.ndarray:
    lib = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    lib[lib <= 0] = 1.0
    x = counts[:, gene_idx].astype(np.float64)
    if sparse.issparse(x):
        x = x.toarray()
    x = np.log1p(x / lib[:, None] * 1e4)
    return x.astype(np.float32)


def select_genes(counts: sparse.spmatrix, var_names: np.ndarray, min_cells: int) -> np.ndarray:
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    keep = detected >= min_cells
    names = var_names.astype(str)
    keep &= names != "CLDN4"  # held out of domain features
    idx = np.flatnonzero(keep)
    if idx.size <= N_HVG:
        return idx
    # variance on a log1p(count) sketch, not the full normalized matrix
    sub = counts[:, idx]
    # mean/var of raw counts is enough to drop near-constant panel genes
    mean = np.asarray(sub.mean(axis=0)).ravel()
    mean_sq = np.asarray(sub.multiply(sub).mean(axis=0)).ravel()
    var = np.clip(mean_sq - mean**2, 0, None)
    order = np.argsort(var)[::-1][:N_HVG]
    return idx[order]


def spatial_graph(xy: np.ndarray, k: int = K_SPATIAL):
    n = xy.shape[0]
    kk = min(k + 1, n)
    nn = NearestNeighbors(n_neighbors=kk, algorithm="kd_tree").fit(xy)
    dist, idx = nn.kneighbors(xy)
    # column 0 is self (distance 0). Characteristic length = median k-th neighbor.
    l = float(np.median(dist[:, -1]))
    l = max(l, 1e-3)
    w = np.exp(-(dist**2) / (2.0 * l * l)).astype(np.float64)
    rows = np.repeat(np.arange(n), kk)
    cols = idx.ravel()
    vals = w.ravel()
    a = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
    a.sum_duplicates()
    rs = np.asarray(a.sum(axis=1)).ravel()
    rs[rs <= 0] = 1.0
    a = sparse.diags(1.0 / rs) @ a
    neigh = idx[:, 1:]  # drop self
    return a, neigh, l


def smooth_pca(x: np.ndarray, adj: sparse.spmatrix, n_pcs: int, seed: int) -> np.ndarray:
    xs = adj @ x
    xs = StandardScaler(with_mean=True, with_std=True).fit_transform(xs)
    n_comp = int(min(n_pcs, xs.shape[0] - 1, xs.shape[1]))
    pca = PCA(n_components=n_comp, svd_solver="randomized", random_state=seed)
    return pca.fit_transform(xs).astype(np.float64)


def kmeans_labels(emb: np.ndarray, k: int, seed: int) -> np.ndarray:
    k_use = int(min(k, emb.shape[0]))
    km = KMeans(n_clusters=k_use, n_init=10, random_state=seed)
    return km.fit_predict(emb).astype(np.int32)


def neighbor_agreement(labels: np.ndarray, neigh: np.ndarray) -> float:
    return float((labels[neigh] == labels[:, None]).mean())


def hmrf_refine(
    emb: np.ndarray,
    labels: np.ndarray,
    neigh: np.ndarray,
    beta: float,
    n_iter: int,
) -> tuple[np.ndarray, int, float]:
    """Jacobi Potts HMRF. Emission is squared distance in the smoothed PCA,
    scaled by the initial within-cluster median so beta is dimensionless
    relative to a 0–1 neighbor-vote fraction.
    """
    labels = labels.copy()
    k = int(labels.max()) + 1
    x2 = (emb**2).sum(axis=1, keepdims=True)

    def centroids(lab: np.ndarray) -> np.ndarray:
        mu = np.zeros((k, emb.shape[1]), dtype=np.float64)
        global_mu = emb.mean(axis=0)
        for c in range(k):
            m = lab == c
            mu[c] = emb[m].mean(axis=0) if m.any() else global_mu
        return mu

    mu0 = centroids(labels)
    d2_own = ((emb - mu0[labels]) ** 2).sum(axis=1)
    s2 = float(np.median(d2_own)) + 1e-6
    changed = 1.0
    it_done = 0
    for it in range(n_iter):
        it_done = it + 1
        mu = centroids(labels)
        m2 = (mu**2).sum(axis=1)
        d2 = x2 + m2[None, :] - 2.0 * emb @ mu.T
        neigh_lab = labels[neigh]
        votes = np.zeros((emb.shape[0], k), dtype=np.float64)
        # mean indicator per label; k is small (6–10)
        for c in range(k):
            votes[:, c] = (neigh_lab == c).mean(axis=1)
        energy = d2 / s2 - beta * votes
        new = energy.argmin(axis=1).astype(np.int32)
        changed = float(np.mean(new != labels))
        labels = new
        if changed < 0.001:
            break
    return labels, it_done, changed


def gene_vector(counts: sparse.spmatrix, var_names: np.ndarray, gene: str) -> np.ndarray | None:
    names = var_names.astype(str)
    hit = np.flatnonzero(names == gene)
    if hit.size == 0:
        return None
    col = counts[:, int(hit[0])]
    if sparse.issparse(col):
        return np.asarray(col.toarray()).ravel().astype(np.float64)
    return np.asarray(col).ravel().astype(np.float64)


def domain_table_for_sample(
    sample: str,
    patient: str,
    labels: np.ndarray,
    comp: np.ndarray,
    is_tumor: np.ndarray,
    is_immune: np.ndarray,
    is_cyt: np.ndarray,
    neigh: np.ndarray,
    cldn4: np.ndarray,
    cldn4_high: np.ndarray,
    cyt50: np.ndarray,
    cyt100: np.ndarray,
    controls: dict[str, np.ndarray],
    sample_immune_frac: float,
    k_requested: int,
    beta_used: float,
    length_scale_um: float,
    agree: float,
) -> list[dict]:
    rows = []
    for d in sorted(np.unique(labels)):
        m = labels == d
        n = int(m.sum())
        n_tumor = int((m & is_tumor).sum())
        n_immune = int((m & is_immune).sum())
        n_cyt = int((m & is_cyt).sum())
        n_stromal = int((m & (comp == 2)).sum())
        n_other = int((m & (comp == 3)).sum())
        counts = np.array([n_tumor, n_immune, n_stromal, n_other], dtype=float)
        h = shannon_from_counts(counts)
        h_even = h / np.log2(len(COMPARTMENTS)) if np.isfinite(h) else float("nan")
        tum = m & is_tumor
        if n_tumor > 0:
            obs_frac = float(is_immune[neigh[tum]].mean())
            mean_cl = float(cldn4[tum].mean())
            frac_hi = float(cldn4_high[tum].mean())
            mean_c50 = float(cyt50[tum].mean())
            mean_c100 = float(cyt100[tum].mean())
        else:
            obs_frac = float("nan")
            mean_cl = float("nan")
            frac_hi = float("nan")
            mean_c50 = float("nan")
            mean_c100 = float("nan")
        nm = normalized_mixing(obs_frac, sample_immune_frac)
        # composition-null spatial mixing (supplement): heterotypic kNN / Gini-Simpson
        neigh_c = comp[neigh[m]]
        obs_het = float((neigh_c != comp[m][:, None]).mean()) if n else float("nan")
        pi = counts / counts.sum() if n else np.zeros(4)
        exp_het = float(1.0 - np.sum(pi**2))
        nm_comp = obs_het / exp_het if exp_het >= PURE_EXPECTED else float("nan")
        row = {
            "sample": sample,
            "patient": patient,
            "domain": int(d),
            "k_domains": int(k_requested),
            "beta": beta_used,
            "length_scale_um": length_scale_um,
            "neighbor_agreement": agree,
            "n_cells": n,
            "n_tumor": n_tumor,
            "n_immune": n_immune,
            "n_cytotoxic": n_cyt,
            "n_stromal": n_stromal,
            "n_other": n_other,
            "frac_tumor": n_tumor / n if n else float("nan"),
            "frac_immune": n_immune / n if n else float("nan"),
            "frac_cytotoxic": n_cyt / n if n else float("nan"),
            "shannon": h,
            "shannon_evenness": h_even,
            "immune_neighbor_frac_tumor": obs_frac,
            "sample_immune_frac": sample_immune_frac,
            "normalized_mixing": nm,
            "heterotypic_frac": obs_het,
            "normalized_mixing_composition_null": nm_comp,
            "mean_cldn4_tumor": mean_cl,
            "frac_cldn4_high_tumor": frac_hi,
            "mean_cytotoxic_50um": mean_c50,
            "mean_cytotoxic_100um": mean_c100,
            "eligible": bool(n >= MIN_DOMAIN_CELLS and n_tumor >= MIN_DOMAIN_TUMOR),
        }
        for g, vec in controls.items():
            row[f"mean_{g}_tumor"] = float(vec[tum].mean()) if n_tumor else float("nan")
        rows.append(row)
    return rows


def assign_high_low(df: pd.DataFrame) -> pd.DataFrame:
    """Within each sample, median-split eligible domains on tumor-cell mean CLDN4."""
    out = df.copy()
    out["cldn4_arm"] = pd.NA
    for sample, idx in out.groupby("sample").groups.items():
        sub = out.loc[idx]
        elig = sub.index[sub["eligible"].to_numpy()]
        if len(elig) < 2:
            continue
        med = float(sub.loc[elig, "mean_cldn4_tumor"].median())
        hi = sub.loc[elig, "mean_cldn4_tumor"] >= med
        # if the median ties every domain, skip the split
        if hi.all() or (~hi).all():
            continue
        out.loc[elig, "cldn4_arm"] = np.where(hi, "high", "low")
    return out


def weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    v = values.to_numpy(dtype=float)
    w = weights.to_numpy(dtype=float)
    ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
    if not ok.any():
        return float("nan")
    return float(np.average(v[ok], weights=w[ok]))


def paired_summary(df: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    rows = []
    for sample, sub in df.groupby("sample"):
        patient = sub["patient"].iloc[0]
        row: dict = {"sample": sample, "patient": patient}
        for arm in ("high", "low"):
            part = sub[sub["cldn4_arm"] == arm]
            row[f"n_domains_{arm}"] = int(len(part))
            row[f"n_tumor_{arm}"] = int(part["n_tumor"].sum()) if len(part) else 0
            for m in metrics:
                row[f"{m}_{arm}"] = weighted_mean(part[m], part["n_tumor"]) if len(part) else float("nan")
        for m in metrics:
            row[f"delta_{m}"] = row[f"{m}_high"] - row[f"{m}_low"]
        rows.append(row)
    return pd.DataFrame(rows)


def donor_summary(sample_df: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    rows = []
    for patient, sub in sample_df.groupby("patient"):
        row: dict = {"patient": patient, "n_samples": int(len(sub))}
        for m in metrics:
            row[f"{m}_high"] = float(sub[f"{m}_high"].mean())
            row[f"{m}_low"] = float(sub[f"{m}_low"].mean())
            row[f"delta_{m}"] = row[f"{m}_high"] - row[f"{m}_low"]
        rows.append(row)
    return pd.DataFrame(rows)


def sign_test(deltas: np.ndarray, alternative: str = "less") -> dict:
    d = np.asarray(deltas, dtype=float)
    d = d[np.isfinite(d) & (d != 0)]
    n = int(d.size)
    if alternative == "less":
        n_hit = int(np.sum(d < 0))
    else:
        n_hit = int(np.sum(d > 0))
    if n == 0:
        return {"n": 0, "n_hit": 0, "p_one": float("nan"), "p_two": float("nan")}
    p_one = float(stats.binomtest(n_hit, n, 0.5, alternative="greater").pvalue)
    p_two = float(stats.binomtest(n_hit, n, 0.5, alternative="two-sided").pvalue)
    return {"n": n, "n_hit": n_hit, "p_one": p_one, "p_two": p_two}


def wilcoxon_less(deltas: np.ndarray) -> dict:
    d = np.asarray(deltas, dtype=float)
    d = d[np.isfinite(d)]
    if np.sum(d != 0) < 1:
        return {"n": int(d.size), "stat": float("nan"), "p_two": float("nan"), "p_less": float("nan")}
    # two-sided, then one-sided (high < low) from the signed rank direction
    res_two = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox")
    res_less = stats.wilcoxon(d, alternative="less", zero_method="wilcox")
    return {
        "n": int(d.size),
        "stat": float(res_two.statistic),
        "p_two": float(res_two.pvalue),
        "p_less": float(res_less.pvalue),
    }


def cldn4_residual_on_epcam(df: pd.DataFrame) -> list[dict]:
    """Within-section Spearman of CLDN4 after a rank-linear residual on EPCAM.

    Asks whether CLDN4 still tracks mixing once the shared epithelial axis is removed.
    """
    rows = []
    if "mean_EPCAM_tumor" not in df.columns:
        return rows
    for sample, sub in df.groupby("sample"):
        elig = sub[sub["eligible"]]
        if len(elig) < 4:
            continue
        x = elig["mean_EPCAM_tumor"].rank().to_numpy(dtype=float)
        y = elig["mean_cldn4_tumor"].rank().to_numpy(dtype=float)
        if not np.isfinite(x).all() or np.nanstd(x) == 0:
            continue
        slope, intercept = np.polyfit(x, y, 1)
        resid = y - (slope * x + intercept)
        for metric in ("normalized_mixing", "shannon", "mean_cytotoxic_50um"):
            b = elig[metric].to_numpy(dtype=float)
            ok = np.isfinite(resid) & np.isfinite(b)
            if ok.sum() < 4:
                continue
            rho, p = stats.spearmanr(resid[ok], b[ok])
            rows.append(
                {
                    "sample": sample,
                    "patient": elig["patient"].iloc[0],
                    "x": "cldn4_resid_epcam",
                    "y": metric,
                    "n": int(ok.sum()),
                    "rho": float(rho),
                    "p": float(p),
                }
            )
    return rows


def spearman_within(df: pd.DataFrame, x: str, y: str) -> list[dict]:
    rows = []
    for sample, sub in df.groupby("sample"):
        elig = sub[sub["eligible"]]
        if len(elig) < 4:
            continue
        a = elig[x].to_numpy(dtype=float)
        b = elig[y].to_numpy(dtype=float)
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum() < 4:
            continue
        rho, p = stats.spearmanr(a[ok], b[ok])
        rows.append({"sample": sample, "patient": elig["patient"].iloc[0], "x": x, "y": y, "n": int(ok.sum()), "rho": float(rho), "p": float(p)})
    return rows


def cell_exclusion_table(
    sample: str,
    patient: str,
    is_tumor: np.ndarray,
    cldn4: np.ndarray,
    cyt50: np.ndarray,
    cyt100: np.ndarray,
) -> dict:
    tum = is_tumor
    vals = cldn4[tum]
    q75 = float(np.quantile(vals, 0.75))
    q50 = float(np.quantile(vals, 0.50))
    high75 = tum & (cldn4 >= max(q75, 0)) & (cldn4 >= 1)
    # if q75 is 0, max(q75,0) & >=1 keeps the count threshold used previously
    low75 = tum & ~high75
    high50 = tum & (cldn4 >= q50) & (cldn4 >= 1)
    # median split among tumor: top half by rank, ties at the median go high if >= median and count>=1
    low50 = tum & ~high50 & (cldn4 < q50 if q50 > 0 else cldn4 < 1)

    def ratio(high_mask, low_mask, vec):
        hi = vec[high_mask]
        lo = vec[low_mask]
        mean_hi = float(hi.mean()) if hi.size else float("nan")
        mean_lo = float(lo.mean()) if lo.size else float("nan")
        r = mean_hi / mean_lo if mean_lo and np.isfinite(mean_lo) and mean_lo > 0 else float("nan")
        return mean_hi, mean_lo, r, int(high_mask.sum()), int(low_mask.sum())

    out = {
        "sample": sample,
        "patient": patient,
        "n_tumor": int(tum.sum()),
        "cldn4_q75": q75,
        "cldn4_q50": q50,
    }
    for name, hi, lo in (("q75", high75, low75), ("median", high50, low50)):
        for radius, vec in (("50", cyt50), ("100", cyt100)):
            mean_hi, mean_lo, r, n_hi, n_lo = ratio(hi, lo, vec)
            out[f"n_high_{name}"] = n_hi
            out[f"n_low_{name}"] = n_lo
            out[f"cyt{radius}_{name}_high"] = mean_hi
            out[f"cyt{radius}_{name}_low"] = mean_lo
            out[f"ratio_cyt{radius}_{name}"] = r
    return out


def cluster_domains(emb: np.ndarray, neigh: np.ndarray, k_domains: int) -> tuple[np.ndarray, float, float, int, float]:
    init = kmeans_labels(emb, k_domains, SEED)
    agree0 = neighbor_agreement(init, neigh)
    labels, n_iter, changed = hmrf_refine(emb, init, neigh, BETA, HMRF_ITERS)
    agree1 = neighbor_agreement(labels, neigh)
    beta_used = BETA
    if agree1 < MIN_AGREE:
        labels, n_iter, changed = hmrf_refine(emb, init, neigh, BETA_IF_FRAGMENTED, HMRF_ITERS)
        agree1 = neighbor_agreement(labels, neigh)
        beta_used = BETA_IF_FRAGMENTED
    return labels, beta_used, agree1, n_iter, agree0


def process_sample(adata_backed, sample_raw: str, genes: np.ndarray, k_list: list[int]) -> dict:
    obs_all = adata_backed.obs
    mask = (obs_all["sample"].astype(str) == sample_raw).to_numpy()
    print(f"{sample_raw}: loading {int(mask.sum())} cells", flush=True)
    sub = adata_backed[mask].to_memory()
    counts = sub.layers["counts"]
    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)
    else:
        counts = counts.tocsr()
    xy = np.asarray(sub.obsm["spatial"], dtype=np.float64) * UM_PER_PX
    # shift to local µm so the KD-tree is well conditioned
    xy = xy - xy.min(axis=0, keepdims=True)
    obs = sub.obs.copy()
    del sub

    canon = SAMPLE_CANON[sample_raw]
    patient = str(obs["patient"].iloc[0])
    comp = compartment_ids(obs["cell_type"])
    is_tumor = comp == 0
    is_immune = comp == 1
    is_cyt = obs["cell_type"].astype(str).isin(CYTOTOXIC).to_numpy()
    sample_immune_frac = float(is_immune.mean())

    min_cells = max(30, int(0.01 * counts.shape[0]))
    gene_idx = select_genes(counts, genes, min_cells)
    print(f"{canon}: {gene_idx.size} genes (CLDN4 held out), immune frac {sample_immune_frac:.3f}", flush=True)
    x = lognorm_selected(counts, gene_idx)
    cldn4 = gene_vector(counts, genes, "CLDN4")
    if cldn4 is None:
        raise RuntimeError("CLDN4 missing from the 960-gene panel")
    controls = {}
    for g in CONTROL_GENES:
        vec = gene_vector(counts, genes, g)
        if vec is not None:
            controls[g] = vec
    del counts

    adj, neigh, length_um = spatial_graph(xy, K_SPATIAL)
    emb = smooth_pca(x, adj, N_PCS, SEED)
    del x, adj

    # cytotoxic neighbor counts for every tumor cell (self is not in the tree)
    cyt50 = np.zeros(len(obs), dtype=np.float64)
    cyt100 = np.zeros(len(obs), dtype=np.float64)
    if is_cyt.any() and is_tumor.any():
        tree = cKDTree(xy[is_cyt])
        cyt50[is_tumor] = tree.query_ball_point(xy[is_tumor], r=50.0, return_length=True)
        cyt100[is_tumor] = tree.query_ball_point(xy[is_tumor], r=100.0, return_length=True)

    q75 = float(np.quantile(cldn4[is_tumor], 0.75)) if is_tumor.any() else float("nan")
    cldn4_high = is_tumor & (cldn4 >= q75) & (cldn4 >= 1)

    primary_k = k_list[0]
    rows = []
    labels_primary = None
    for k_domains in k_list:
        labels, beta_used, agree1, n_iter, agree0 = cluster_domains(emb, neigh, k_domains)
        print(
            f"{canon}: k={k_domains}, l={length_um:.1f} µm, agreement {agree0:.3f} -> {agree1:.3f}, "
            f"beta={beta_used}, hmrf iters={n_iter}",
            flush=True,
        )
        rows.extend(
            domain_table_for_sample(
                canon,
                patient,
                labels,
                comp,
                is_tumor,
                is_immune,
                is_cyt,
                neigh,
                cldn4,
                cldn4_high,
                cyt50,
                cyt100,
                controls,
                sample_immune_frac,
                k_domains,
                beta_used,
                length_um,
                agree1,
            )
        )
        if k_domains == primary_k:
            labels_primary = labels
    excl = cell_exclusion_table(canon, patient, is_tumor, cldn4, cyt50, cyt100)
    labels = labels_primary

    # author Giotto niches, same indices, descriptive
    niche = obs["niche"].astype(str).to_numpy()
    niche_rows = []
    for name in sorted(pd.unique(niche)):
        m = niche == name
        # reuse domain metric code by passing a single label
        fake = np.zeros(len(obs), dtype=np.int32)
        fake[m] = 0
        # only score the niche's own cells: restrict arrays
        # Build from the cells in the niche directly.
        n = int(m.sum())
        n_tumor = int((m & is_tumor).sum())
        n_immune = int((m & is_immune).sum())
        n_cyt = int((m & is_cyt).sum())
        n_stromal = int((m & (comp == 2)).sum())
        n_other = n - n_tumor - n_immune - n_stromal
        h = shannon_from_counts([n_tumor, n_immune, n_stromal, n_other])
        tum = m & is_tumor
        if n_tumor:
            obs_frac = float(is_immune[neigh[tum]].mean())
            mean_cl = float(cldn4[tum].mean())
            mean_c50 = float(cyt50[tum].mean())
            mean_c100 = float(cyt100[tum].mean())
        else:
            obs_frac = mean_cl = mean_c50 = mean_c100 = float("nan")
        niche_rows.append(
            {
                "sample": canon,
                "patient": patient,
                "niche": name,
                "n_cells": n,
                "n_tumor": n_tumor,
                "n_immune": n_immune,
                "n_cytotoxic": n_cyt,
                "shannon": h,
                "shannon_evenness": h / np.log2(4) if np.isfinite(h) else float("nan"),
                "normalized_mixing": normalized_mixing(obs_frac, sample_immune_frac),
                "mean_cldn4_tumor": mean_cl,
                "mean_cytotoxic_50um": mean_c50,
                "mean_cytotoxic_100um": mean_c100,
                "frac_immune": n_immune / n if n else float("nan"),
            }
        )

    # Crop the densest 700 µm radius so domain contiguity is visible.
    rng = np.random.default_rng(SEED)
    n_try = min(60, len(obs))
    centers = xy[rng.choice(len(obs), size=n_try, replace=False)]
    crop_tree = cKDTree(xy)
    n_near = np.asarray(crop_tree.query_ball_point(centers, r=700.0, return_length=True))
    center = centers[int(np.argmax(n_near))]
    take = np.asarray(crop_tree.query_ball_point(center, r=700.0), dtype=int)
    if take.size > 40000:
        take = rng.choice(take, size=40000, replace=False)
    plot_df = pd.DataFrame(
        {
            "x": xy[take, 0],
            "y": xy[take, 1],
            "domain": labels[take],
            "cldn4": cldn4[take],
            "is_tumor": is_tumor[take],
            "comp": comp[take],
        }
    )
    return {
        "domains": rows,
        "exclusion": excl,
        "niches": niche_rows,
        "plot": plot_df,
        "canon": canon,
    }


METRICS = [
    "shannon",
    "shannon_evenness",
    "normalized_mixing",
    "frac_immune",
    "mean_cytotoxic_50um",
    "mean_cytotoxic_100um",
    "mean_cldn4_tumor",
]


def fmt(x, digits=3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    ax = abs(float(x))
    if ax != 0 and (ax < 0.001 or ax >= 1000):
        return f"{float(x):.2e}"
    return f"{float(x):.{digits}f}"


def write_results(
    domain_df: pd.DataFrame,
    paired: pd.DataFrame,
    donors: pd.DataFrame,
    tests: dict,
    excl: pd.DataFrame,
    niches: pd.DataFrame,
    spearman_rows: list[dict],
    run_meta: dict,
) -> None:
    lines = []
    lines.append("# CosMx spatial domains: CLDN4 enrichment vs immune mixing")
    lines.append("")
    lines.append("Public CosMx SMI NSCLC, figshare 25976224 (`cosmx_human_nsclc_clustered.h5ad`): **765,771 cells, 960 genes, 8 sections, 5 donors** (He et al. 2022). CLDN4-only. No private cohort.")
    lines.append("")
    lines.append("This layer does **not** replace the locked cell-level exclusion result (CLDN4-high tumor has fewer nearby cytotoxic cells at 50/100 µm; 8/8 sections, 5/5 donors, one-sided sign P = 0.031). It asks whether that exclusion is visible as spaGCN-style spatial domains whose CLDN4 enrichment tracks Shannon diversity and a normalized mixing index.")
    lines.append("")
    lines.append("## Method")
    lines.append("")
    lines.append("Per section, histology-free SpaGCN-style aggregation (Hu et al. 2021) on a sparse graph, because a dense SpaGCN adjacency is not tractable at ~10^5 cells:")
    lines.append("")
    lines.append(f"- Coordinates: object pixels × {UM_PER_PX} µm/px (calibrated to the public Giotto locs; median nearest centroid = 9.0 µm on Lung5_Rep1 FOV 1).")
    lines.append(f"- Expression: raw `counts`, library-size normalize to 10^4, log1p. Up to {N_HVG} highest-variance genes detected in ≥1% of cells. **CLDN4 is removed from these features** before clustering.")
    lines.append(f"- Spatial graph: {K_SPATIAL} nearest neighbors, Gaussian weights with length scale l = median distance to the {K_SPATIAL}th neighbor, row-normalized, self included. One graph convolution, then PCA ({N_PCS} components).")
    lines.append(f"- Domains: KMeans k={K_DOMAINS} (fixed so every section has the same domain count; SpaGCN's original Louvain step is the piece replaced here), seed {SEED}.")
    lines.append(f"- Giotto-style Potts HMRF refinement: Jacobi updates, emission = squared PCA distance / initial within-cluster median, Potts term = beta × fraction of the {K_SPATIAL} neighbors sharing the label. Primary beta={BETA}. If neighbor agreement stays below {MIN_AGREE}, beta is raised to {BETA_IF_FRAGMENTED} (contiguity QC, not an outcome search).")
    lines.append("- CLDN4-high tumor cells, for the cell-level bridge only: within-section tumor CLDN4 count ≥ the tumor 75th percentile and ≥ 1. Domain split is separate: eligible domains (n≥100 and n_tumor≥30) are median-split on **mean CLDN4 of tumor cells in the domain**.")
    lines.append("- Shannon: H = −Σ p log2 p on four compartments (tumor, immune, stromal, other). Evenness = H / 2.")
    lines.append("- Normalized mixing: among tumor cells in the domain, the mean fraction of the 10 spatial neighbors that are immune, divided by the **section** immune fraction. Values < 1 mean tumor cells in that domain see fewer immune neighbors than the section base rate.")
    lines.append("- Exclusion link: mean number of author cytotoxic cells (T CD8 naive, T CD8 memory, NK) inside 50 µm and 100 µm of tumor cells. Sample tallies use tumor-cell-weighted means of domains. Donor tallies average the sample values (Lung5 = 3 sections, Lung9 = 2).")
    lines.append("- Tests: Wilcoxon signed-rank on 8 section deltas (high − low), and a sign test on 5 donor deltas. One-sided tests use the exclusion direction (high < low for Shannon, mixing, immune fraction, and cytotoxic counts). Two-sided p-values are reported next to them.")
    lines.append("- Epithelial controls (EPCAM, KRT8, KRT18), same tumor-cell mean vs immune-neighbor association, are descriptive. They are not a claim that CLDN4 is the strongest surface gene.")
    lines.append("")
    lines.append("Author `niche` labels already in the object (tumor interior, stroma, immune, …) are the published Giotto niches. Those niches were built from neighborhood composition, so immune fraction is partly by construction. Tumor-cell mean CLDN4 is also similar across tumor-interior and immune niches, because the immune niche still contains some CLDN4-high tumor cells. Niches are descriptive, not the primary test.")
    lines.append("")
    lines.append("## Cell-level bridge to exclusion")
    lines.append("")
    lines.append("Recomputed on this object (author cell types, cytotoxic = CD8 T + NK). Ratio = mean neighbor count around CLDN4-high tumor / mean around the remaining tumor. Q75 is the primary cut; the median cut is shown beside it.")
    lines.append("")
    lines.append("| Section | Donor | 50 µm Q75 high | 50 µm Q75 low | ratio 50 | ratio 100 | ratio 50 median |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for _, r in excl.sort_values("sample").iterrows():
        lines.append(
            f"| {r['sample']} | {r['patient']} | {fmt(r['cyt50_q75_high'])} | {fmt(r['cyt50_q75_low'])} | "
            f"{fmt(r['ratio_cyt50_q75'])} | {fmt(r['ratio_cyt100_q75'])} | {fmt(r['ratio_cyt50_median'])} |"
        )
    n50 = int(np.sum(excl["ratio_cyt50_q75"] < 1))
    n100 = int(np.sum(excl["ratio_cyt100_q75"] < 1))
    d50 = excl.groupby("patient")["ratio_cyt50_q75"].mean()
    d100 = excl.groupby("patient")["ratio_cyt100_q75"].mean()
    lines.append("")
    lines.append(
        f"Q75 count-ratio < 1 in **{n50}/8** sections at 50 µm and **{n100}/8** at 100 µm; "
        f"donor means < 1 in **{int(np.sum(d50 < 1))}/5** (50 µm) and **{int(np.sum(d100 < 1))}/5** (100 µm). "
        f"Median section ratios: 50 µm {fmt(float(excl['ratio_cyt50_q75'].median()))}, "
        f"100 µm {fmt(float(excl['ratio_cyt100_q75'].median()))}. "
        "The locked 0.36 / 0.52 figures are the prior cell-level result and are not overwritten by these recomputed ratios."
    )
    lines.append("")
    lines.append("## Domain result")
    lines.append("")
    lines.append(
        f"Eligible domains: {int(domain_df['eligible'].sum())} of {len(domain_df)}. "
        f"Median Gaussian length scale {fmt(float(domain_df.groupby('sample')['length_scale_um'].first().median()))} µm. "
        f"Median post-HMRF neighbor agreement {fmt(float(domain_df.groupby('sample')['neighbor_agreement'].first().median()))}."
    )
    lines.append("")
    lines.append("Section tallies (tumor-cell-weighted means of CLDN4-high vs CLDN4-low domains):")
    lines.append("")
    lines.append("| Section | Donor | Shannon high | Shannon low | Δ Shannon | Mixing high | Mixing low | Δ mixing | Cyt50 high | Cyt50 low |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for _, r in paired.sort_values("sample").iterrows():
        lines.append(
            f"| {r['sample']} | {r['patient']} | {fmt(r['shannon_high'])} | {fmt(r['shannon_low'])} | {fmt(r['delta_shannon'])} | "
            f"{fmt(r['normalized_mixing_high'])} | {fmt(r['normalized_mixing_low'])} | {fmt(r['delta_normalized_mixing'])} | "
            f"{fmt(r['mean_cytotoxic_50um_high'])} | {fmt(r['mean_cytotoxic_50um_low'])} |"
        )
    lines.append("")
    lines.append("Donor tallies (unweighted mean of that donor's sections):")
    lines.append("")
    lines.append("| Donor | Sections | Δ Shannon | Δ mixing | Δ immune frac | Δ cyt 50 µm | Δ cyt 100 µm |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for _, r in donors.sort_values("patient").iterrows():
        lines.append(
            f"| {r['patient']} | {int(r['n_samples'])} | {fmt(r['delta_shannon'])} | {fmt(r['delta_normalized_mixing'])} | "
            f"{fmt(r['delta_frac_immune'])} | {fmt(r['delta_mean_cytotoxic_50um'])} | {fmt(r['delta_mean_cytotoxic_100um'])} |"
        )
    lines.append("")
    lines.append("### Tests (high − low; exclusion direction is negative)")
    lines.append("")
    lines.append("| Metric | Sections <0 | Wilcoxon p two-sided | Wilcoxon p (high<low) | Donors <0 | Sign p one-sided | Sign p two-sided |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for key, label in (
        ("shannon", "Shannon"),
        ("normalized_mixing", "Normalized mixing"),
        ("frac_immune", "Immune fraction"),
        ("mean_cytotoxic_50um", "Cytotoxic count 50 µm"),
        ("mean_cytotoxic_100um", "Cytotoxic count 100 µm"),
    ):
        t = tests[key]
        lines.append(
            f"| {label} | {t['section_sign']['n_hit']}/{t['section_sign']['n']} | {fmt(t['wilcoxon']['p_two'], 4)} | "
            f"{fmt(t['wilcoxon']['p_less'], 4)} | {t['donor_sign']['n_hit']}/{t['donor_sign']['n']} | "
            f"{fmt(t['donor_sign']['p_one'], 4)} | {fmt(t['donor_sign']['p_two'], 4)} |"
        )
    lines.append("")
    lines.append("Within-section Spearman (eligible domains, tumor-mean CLDN4 vs metric):")
    lines.append("")
    if spearman_rows:
        sp = pd.DataFrame(spearman_rows)
        lines.append("| Metric | Sections with ρ<0 | Median ρ |")
        lines.append("|---|---:|---:|")
        for y, label in (
            ("shannon", "Shannon"),
            ("normalized_mixing", "Normalized mixing"),
            ("frac_immune", "Immune fraction"),
            ("mean_cytotoxic_50um", "Cytotoxic count 50 µm"),
        ):
            part = sp[(sp["y"] == y) & (sp["x"] == "mean_cldn4_tumor")]
            if part.empty:
                continue
            lines.append(
                f"| {label} | {int(np.sum(part['rho'] < 0))}/{len(part)} | {fmt(float(part['rho'].median()))} |"
            )
        resid = sp[sp["x"] == "cldn4_resid_epcam"]
        if not resid.empty:
            lines.append("")
            lines.append(
                "CLDN4 rank residualized on EPCAM within each section (epithelial-program control). "
                "A negative residual ρ would mean CLDN4 still tracks that metric after the shared epithelial axis is removed."
            )
            lines.append("")
            lines.append("| Residual CLDN4 vs | Sections with ρ<0 | Median ρ |")
            lines.append("|---|---:|---:|")
            for y, label in (
                ("shannon", "Shannon"),
                ("normalized_mixing", "Normalized mixing"),
                ("mean_cytotoxic_50um", "Cytotoxic count 50 µm"),
            ):
                part = resid[resid["y"] == y]
                if part.empty:
                    continue
                lines.append(
                    f"| {label} | {int(np.sum(part['rho'] < 0))}/{len(part)} | {fmt(float(part['rho'].median()))} |"
                )
        lines.append("")
        lines.append("Epithelial control genes, within-section Spearman of the tumor-cell mean against normalized mixing:")
        lines.append("")
        lines.append("| Gene | Sections with ρ<0 | Median ρ |")
        lines.append("|---|---:|---:|")
        for g in CONTROL_GENES:
            part = sp[(sp["x"] == f"mean_{g}_tumor") & (sp["y"] == "normalized_mixing")]
            if part.empty:
                continue
            lines.append(
                f"| {g} | {int(np.sum(part['rho'] < 0))}/{len(part)} | {fmt(float(part['rho'].median()))} |"
            )
    lines.append("")
    sens_path = OUT / "tables" / "sensitivity_k.csv"
    if sens_path.exists():
        sens = pd.read_csv(sens_path)
        lines.append("## Sensitivity to domain count")
        lines.append("")
        lines.append("Same smoothed embedding, KMeans k = 6, 8, and 10, each with its own HMRF pass. Deltas are still high − low CLDN4.")
        lines.append("")
        lines.append("| k | Metric | Sections <0 | Wilcoxon p (high<low) | Donors <0 | Sign p one-sided | Median section Δ |")
        lines.append("|---:|---|---:|---:|---:|---:|---:|")
        for _, r in sens.sort_values(["metric", "k_domains"]).iterrows():
            lines.append(
                f"| {int(r['k_domains'])} | {r['metric']} | {r['sections_lt0']} | {fmt(r['wilcoxon_p_less'], 4)} | "
                f"{r['donors_lt0']} | {fmt(r['sign_p_one'], 4)} | {fmt(r['median_section_delta'])} |"
            )
        lines.append("")
    lines.append("## How this sits with exclusion")
    lines.append("")
    n_q75_50 = int(np.sum(excl["ratio_cyt50_q75"] < 1))
    n_q75_100 = int(np.sum(excl["ratio_cyt100_q75"] < 1))
    n_sec = int(excl["sample"].nunique())

    def tally(key: str) -> str:
        sec = tests[key]["section_sign"]
        don = tests[key]["donor_sign"]
        return f"{sec['n_hit']}/{sec['n']} sections and {don['n_hit']}/{don['n']} donors"

    lines.append(
        "Two different aggregations are in this file. "
        f"Cell-level Q75 (each tumor cell scored by its own CLDN4 count) gives a cytotoxic-count ratio < 1 in {n_q75_50}/{n_sec} sections at 50 µm and {n_q75_100}/{n_sec} at 100 µm. "
        "That estimator is not the locked 0.36 / 0.52 result, and these recomputed ratios do not replace it. "
        "Domain-level aggregation puts tumor cells into spaGCN/HMRF domains and compares CLDN4-high vs CLDN4-low domains. "
        f"Shannon is lower in the CLDN4-high domains in {tally('shannon')}; "
        f"normalized mixing in {tally('normalized_mixing')}; "
        f"50 µm cytotoxic counts in {tally('mean_cytotoxic_50um')}; "
        f"100 µm cytotoxic counts in {tally('mean_cytotoxic_100um')}."
    )
    lines.append("")
    lines.append(
        "The domain contrast follows the epithelial program, not a CLDN4-only axis. "
        "EPCAM, KRT8, and KRT18 tumor-means correlate with lower normalized mixing at least as consistently as CLDN4, and the CLDN4 residual after EPCAM does not stay negative. "
        "Read the domain result as: CLDN4-enriched epithelial domains are the immune-poor domains (exclusion geography). "
        "Do not read it as public evidence that CLDN4 outranks other epithelial genes. "
        "Effector-cell state (GZMB/PRF1/NKG7/IFNG) is not retested here."
    )
    lines.append("")
    lines.append("## Author niches (descriptive)")
    lines.append("")
    if len(niches):
        g = (
            niches.groupby("niche", as_index=False)
            .agg(
                n_sections=("sample", "nunique"),
                mean_cldn4=("mean_cldn4_tumor", "mean"),
                mean_shannon=("shannon", "mean"),
                mean_mixing=("normalized_mixing", "mean"),
                mean_cyt50=("mean_cytotoxic_50um", "mean"),
                mean_frac_immune=("frac_immune", "mean"),
            )
            .sort_values("mean_cldn4", ascending=False)
        )
        lines.append("| Niche | Sections | Mean tumor CLDN4 | Shannon | Normalized mixing | Cytotoxic 50 µm | Immune fraction |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for _, r in g.iterrows():
            lines.append(
                f"| {r['niche']} | {int(r['n_sections'])} | {fmt(r['mean_cldn4'])} | {fmt(r['mean_shannon'])} | "
                f"{fmt(r['mean_mixing'])} | {fmt(r['mean_cyt50'])} | {fmt(r['mean_frac_immune'])} |"
            )
    lines.append("")
    lines.append("## Run")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/download_cosmx_figshare_25976224.py")
    lines.append("python3 scripts/cosmx_spagcn_domain_mixing.py")
    lines.append("```")
    lines.append("")
    lines.append(f"Seed {SEED}. k={run_meta['k']}. Cells used: {run_meta['n_cells']}.")
    lines.append("")
    text = "\n".join(lines) + "\n"
    (OUT / "RESULTS.md").write_text(text)
    (ROOT / "RESULTS.md").write_text(text)


def make_figures(domain_df: pd.DataFrame, paired: pd.DataFrame, plots: dict[str, pd.DataFrame], excl: pd.DataFrame) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # paired Shannon and mixing
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2))
    for ax, metric, title in (
        (axes[0], "shannon", "Shannon diversity"),
        (axes[1], "normalized_mixing", "Normalized mixing"),
    ):
        for _, r in paired.iterrows():
            ax.plot([0, 1], [r[f"{metric}_low"], r[f"{metric}_high"]], color="#4d4d4d", lw=0.8, zorder=1)
            ax.scatter([0], [r[f"{metric}_low"]], color="#2c7fb8", s=28, zorder=2)
            ax.scatter([1], [r[f"{metric}_high"]], color="#d95f0e", s=28, zorder=2)
        ax.set_xticks([0, 1], ["CLDN4-low\ndomains", "CLDN4-high\ndomains"])
        ax.set_title(title)
        ax.set_xlim(-0.3, 1.3)
    axes[0].set_ylabel("Tumor-cell-weighted domain mean")
    fig.tight_layout()
    fig.savefig(fig_dir / "section_paired_shannon_mixing.png", dpi=160)
    fig.savefig(fig_dir / "section_paired_shannon_mixing.pdf")
    plt.close(fig)

    elig = domain_df[domain_df["eligible"]].copy()
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.2))
    samples = sorted(elig["sample"].unique())
    cmap = plt.colormaps["tab10"].resampled(max(len(samples), 1))
    for ax, y, ylab in (
        (axes[0], "shannon", "Shannon"),
        (axes[1], "normalized_mixing", "Normalized mixing"),
    ):
        for i, s in enumerate(samples):
            part = elig[elig["sample"] == s]
            ax.scatter(part["mean_cldn4_tumor"], part[y], s=22, color=cmap(i), label=s, alpha=0.9)
        ax.set_xlabel("Mean CLDN4 in tumor cells (counts)")
        ax.set_ylabel(ylab)
        ax.axhline(1.0 if y == "normalized_mixing" else np.nan, color="#999999", lw=0.6, ls="--")
    axes[1].legend(fontsize=6, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(fig_dir / "domain_cldn4_vs_mixing.png", dpi=160)
    fig.savefig(fig_dir / "domain_cldn4_vs_mixing.pdf")
    plt.close(fig)

    # cell-level exclusion ratios
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    excl_s = excl.sort_values("sample")
    x = np.arange(len(excl_s))
    ax.bar(x - 0.15, excl_s["ratio_cyt50_q75"], width=0.3, label="50 µm", color="#2c7fb8")
    ax.bar(x + 0.15, excl_s["ratio_cyt100_q75"], width=0.3, label="100 µm", color="#7fcdbb")
    ax.axhline(1.0, color="#666666", lw=0.8, ls="--")
    ax.set_xticks(x, excl_s["sample"], rotation=40, ha="right")
    ax.set_ylabel("Cytotoxic neighbor-count ratio\n(CLDN4-high / CLDN4-low tumor)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "cell_exclusion_ratios.png", dpi=160)
    fig.savefig(fig_dir / "cell_exclusion_ratios.pdf")
    plt.close(fig)

    # spatial maps, pre-specified sections (densest 700 µm crop saved in `plots`)
    show = [s for s in ("Lung12", "Lung6") if s in plots]
    if show:
        fig, axes = plt.subplots(len(show), 2, figsize=(8.2, 4.0 * len(show)), squeeze=False)
        arm_lookup = {
            (r.sample, int(r.domain)): r.cldn4_arm for r in domain_df.itertuples(index=False)
        }
        for i, s in enumerate(show):
            p = plots[s]
            arm_color = []
            for d in p["domain"].to_numpy():
                arm = arm_lookup.get((s, int(d)))
                if arm == "high":
                    arm_color.append("#d95f0e")
                elif arm == "low":
                    arm_color.append("#2c7fb8")
                else:
                    arm_color.append("#bdbdbd")
            axes[i, 0].scatter(p["x"], p["y"], c=arm_color, s=2.0, linewidths=0, rasterized=True)
            axes[i, 0].set_title(f"{s}: CLDN4-high (orange) vs low (blue)")
            tum = p["is_tumor"].to_numpy()
            axes[i, 1].scatter(p.loc[~tum, "x"], p.loc[~tum, "y"], c="#d9d9d9", s=1.2, linewidths=0, rasterized=True)
            sc = axes[i, 1].scatter(
                p.loc[tum, "x"],
                p.loc[tum, "y"],
                c=np.log1p(p.loc[tum, "cldn4"]),
                s=2.0,
                cmap="viridis",
                linewidths=0,
                rasterized=True,
            )
            axes[i, 1].set_title(f"{s} tumor log1p(CLDN4)")
            fig.colorbar(sc, ax=axes[i, 1], fraction=0.046, pad=0.04)
            for ax in axes[i]:
                ax.set_aspect("equal")
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xlabel("µm")
        fig.tight_layout()
        fig.savefig(fig_dir / "spatial_lung12_lung6.png", dpi=160)
        fig.savefig(fig_dir / "spatial_lung12_lung6.pdf")
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--samples", nargs="*", default=None, help="Raw sample names in the h5ad (default: all 8)")
    parser.add_argument("--k", type=int, default=K_DOMAINS, help="Primary domain count. k=6 and k=10 are always run as sensitivity.")
    args = parser.parse_args()
    self_check()
    if not args.h5ad.exists():
        raise SystemExit(f"missing {args.h5ad}; run scripts/download_cosmx_figshare_25976224.py")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    adata = read_h5ad(args.h5ad, backed="r")
    if adata.n_obs != 765771:
        print(f"warning: expected 765771 cells, found {adata.n_obs}", flush=True)
    genes = np.asarray(adata.var_names)
    raw_samples = list(SAMPLE_CANON)
    if args.samples:
        raw_samples = args.samples

    domain_rows: list[dict] = []
    excl_rows: list[dict] = []
    niche_rows: list[dict] = []
    plots: dict[str, pd.DataFrame] = {}
    for raw in raw_samples:
        if raw not in SAMPLE_CANON:
            raise SystemExit(f"unknown sample {raw}")
        k_list = [args.k]
        for extra in (6, 10):
            if extra not in k_list:
                k_list.append(extra)
        result = process_sample(adata, raw, genes, k_list)
        domain_rows.extend(result["domains"])
        excl_rows.append(result["exclusion"])
        niche_rows.extend(result["niches"])
        plots[result["canon"]] = result["plot"]

    domain_all = pd.DataFrame(domain_rows)
    domain_df = assign_high_low(domain_all[domain_all["k_domains"] == args.k].copy())
    # enrichment relative to the section's tumor-cell mean, eligible domains only use raw mean for the split
    sec_rows = []
    for sample, g in domain_df[domain_df["eligible"]].groupby("sample"):
        w = g["n_tumor"].to_numpy(dtype=float)
        v = g["mean_cldn4_tumor"].to_numpy(dtype=float)
        sec_rows.append((sample, float(np.average(v, weights=w)) if w.sum() else float("nan")))
    sec_mean = dict(sec_rows)
    domain_df["cldn4_enrichment"] = domain_df["mean_cldn4_tumor"] / domain_df["sample"].map(sec_mean)

    paired = paired_summary(domain_df, METRICS)
    donors = donor_summary(paired, METRICS)
    excl = pd.DataFrame(excl_rows)
    niches = pd.DataFrame(niche_rows)

    tests = {}
    for m in ("shannon", "normalized_mixing", "frac_immune", "mean_cytotoxic_50um", "mean_cytotoxic_100um"):
        tests[m] = {
            "wilcoxon": wilcoxon_less(paired[f"delta_{m}"].to_numpy()),
            "section_sign": sign_test(paired[f"delta_{m}"].to_numpy(), "less"),
            "donor_sign": sign_test(donors[f"delta_{m}"].to_numpy(), "less"),
        }

    spearman_rows: list[dict] = []
    for y in ("shannon", "normalized_mixing", "frac_immune", "mean_cytotoxic_50um", "mean_cytotoxic_100um"):
        spearman_rows.extend(spearman_within(domain_df, "mean_cldn4_tumor", y))
    for g in CONTROL_GENES:
        col = f"mean_{g}_tumor"
        if col in domain_df.columns:
            spearman_rows.extend(spearman_within(domain_df, col, "normalized_mixing"))
            spearman_rows.extend(spearman_within(domain_df, col, "shannon"))
    spearman_rows.extend(cldn4_residual_on_epcam(domain_df))

    sens_rows = []
    for k_value, part in domain_all.groupby("k_domains"):
        part = assign_high_low(part.copy())
        paired_k = paired_summary(part, METRICS)
        donors_k = donor_summary(paired_k, METRICS)
        for m in ("shannon", "normalized_mixing", "mean_cytotoxic_50um"):
            sec = sign_test(paired_k[f"delta_{m}"].to_numpy(), "less")
            don = sign_test(donors_k[f"delta_{m}"].to_numpy(), "less")
            wx = wilcoxon_less(paired_k[f"delta_{m}"].to_numpy())
            sens_rows.append(
                {
                    "k_domains": int(k_value),
                    "metric": m,
                    "sections_lt0": f"{sec['n_hit']}/{sec['n']}",
                    "wilcoxon_p_less": wx["p_less"],
                    "donors_lt0": f"{don['n_hit']}/{don['n']}",
                    "sign_p_one": don["p_one"],
                    "median_section_delta": float(np.nanmedian(paired_k[f"delta_{m}"].to_numpy())),
                }
            )
    pd.DataFrame(sens_rows).to_csv(OUT / "tables" / "sensitivity_k.csv", index=False)

    domain_df.to_csv(OUT / "tables" / "domain_metrics.csv", index=False)
    paired.to_csv(OUT / "tables" / "sample_paired.csv", index=False)
    donors.to_csv(OUT / "tables" / "donor_paired.csv", index=False)
    excl.to_csv(OUT / "tables" / "cell_exclusion_ratios.csv", index=False)
    niches.to_csv(OUT / "tables" / "author_niche_metrics.csv", index=False)
    pd.DataFrame(spearman_rows).to_csv(OUT / "tables" / "within_sample_spearman.csv", index=False)

    summary = {
        "n_cells_object": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "um_per_px": UM_PER_PX,
        "k_spatial": K_SPATIAL,
        "k_domains": args.k,
        "beta": BETA,
        "seed": SEED,
        "cldn4_held_out": True,
        "tests": tests,
        "samples": paired.to_dict(orient="records"),
        "donors": donors.to_dict(orient="records"),
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    write_results(
        domain_df,
        paired,
        donors,
        tests,
        excl,
        niches,
        spearman_rows,
        {"k": args.k, "n_cells": int(adata.n_obs)},
    )
    make_figures(domain_df, paired, plots, excl)
    print(f"wrote {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
