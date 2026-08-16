#!/usr/bin/env python3
"""Harmony integration of public neoadjuvant lung scRNA (processed UMI only).

User A3 / GSE207422 CopyKAT is taken as given and is not re-audited.
Patient (or post-treatment sample) is the unit. No cell-level p-values.
"""
from __future__ import annotations

import gzip
import json
import sys
import warnings
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"
TAB = ROOT / "tables"
SCR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCR))
from gene_sets import (  # noqa: E402
    LINEAGE_SCORES,
    NORMAL_LUNG,
    TARGETS,
    core_panel,
)

FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(7)
N_PCS = 30
N_HVG = 800
MAX_PER_SAMPLE = 1200
MIN_MAL = 20
MIN_TNK = 20
N_NEIGHBORS = 15


# ---------------------------------------------------------------------------
# GEO helpers
# ---------------------------------------------------------------------------

def _feat_symbols(path: Path) -> list[str]:
    out = []
    with gzip.open(path, "rt") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2 and p[1] and not p[1].startswith("ENSG"):
                out.append(p[1])
            else:
                out.append(p[0])
    return out


def _gene_rows_1based(path: Path, wanted: set[str]) -> dict[str, int]:
    """Map symbol -> 1-based MTX row (first hit)."""
    hit = {}
    with gzip.open(path, "rt") as f:
        for i, ln in enumerate(f, start=1):
            p = ln.rstrip("\n").split("\t")
            cands = [p[0]]
            if len(p) >= 2:
                cands.append(p[1])
            for c in cands:
                if c in wanted and c not in hit:
                    hit[c] = i
    return hit


def stream_mtx(path: Path, keep_1based: dict[str, int], n_cells: int) -> tuple[np.ndarray, list[str]]:
    """Extract selected MTX rows (genes × cells, 1-based) into cells × genes float32."""
    genes = [g for g, _ in sorted(keep_1based.items(), key=lambda kv: kv[1])]
    row_to_j = {keep_1based[g]: j for j, g in enumerate(genes)}
    X = np.zeros((n_cells, len(genes)), dtype=np.float32)
    print(f"  stream MTX {path.name} keep={len(genes)} n_cells={n_cells}", flush=True)
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] != "%":
                break
        n_hit = 0
        for line in f:
            a, b, c = line.split()
            j = row_to_j.get(int(a))
            if j is None:
                continue
            X[int(b) - 1, j] = float(c)
            n_hit += 1
    print(f"  MTX nnz kept={n_hit}", flush=True)
    return X, genes


def gse207422_stats_and_extract(
    path: Path, core: list[str], n_hvg: int
) -> tuple[np.ndarray, list[str], np.ndarray, list[str]]:
    """One-and-a-half pass: means/vars + nUMI, then extract core+HVG on a second read.

    Returns X (cells × genes), genes, nUMI, barcodes.
    """
    print("GSE207422 pass1 mean/var + nUMI", flush=True)
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        n = len(barcodes)
        nUMI = np.zeros(n, dtype=np.float64)
        means = []
        m2s = []
        names = []
        for line in f:
            gene, rest = line.split("\t", 1)
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            nUMI += vals
            mu = float(vals.mean())
            var = float(vals.var())
            names.append(gene)
            means.append(mu)
            m2s.append(var)
    means = np.asarray(means)
    var = np.asarray(m2s)
    # dispersion among detected genes
    disp = var / (means + 1e-6)
    name_to_i = {g: i for i, g in enumerate(names)}
    must = [g for g in core if g in name_to_i]
    remaining = [g for g in names if g not in set(must)]
    order = sorted(remaining, key=lambda g: disp[name_to_i[g]], reverse=True)
    pick = must + [g for g in order if g][: max(0, n_hvg - len(must))]
    pick_set = set(pick)
    print(f"  n_cells={n} n_genes_full={len(names)} extract={len(pick)}", flush=True)

    print("GSE207422 pass2 extract", flush=True)
    X = np.zeros((n, len(pick)), dtype=np.float32)
    col = {g: j for j, g in enumerate(pick)}
    with gzip.open(path, "rt") as f:
        f.readline()
        for line in f:
            gene, rest = line.split("\t", 1)
            j = col.get(gene)
            if j is None:
                continue
            X[:, j] = np.fromstring(rest, sep="\t", dtype=np.float32)
    return X, pick, nUMI.astype(np.float32), barcodes


def map_mpr(x) -> str | float:
    s = str(x).strip()
    if s in {"non-MPR", "NMPR", "nonMPR", "Non-MPR", "non-mpr"}:
        return "NMPR"
    if s in {"MPR", "pCR", "MPR (pCR)", "pCR/MPR"}:
        return "MPR"
    if s in {"NE", "nan", "None", ""}:
        return np.nan
    return np.nan


def map_histo(x) -> str:
    s = str(x).strip().lower()
    if s in {"adeno", "luad", "adenocarcinoma"}:
        return "LUAD"
    if s in {"squamous", "lusc", "sqcc"}:
        return "LUSC"
    if s in {"asc", "adeno-squamous", "adenosquamous"}:
        return "ASC"
    return "unknown"


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------

def exact_mwu_p(x: np.ndarray, y: np.ndarray) -> float | None:
    """Two-sided exact Wilcoxon/MWU by enumerating assignments. None if too many."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n, k = len(x) + len(y), len(x)
    from math import comb

    if n == 0 or k == 0 or k == n or comb(n, k) > 20000:
        return None
    vals = np.concatenate([x, y])
    ranks = stats.rankdata(vals)
    obs = ranks[:k].sum()
    extreme = 0
    total = 0
    for idx in combinations(range(n), k):
        s = ranks[list(idx)].sum()
        total += 1
        if abs(s - ranks.mean() * k) >= abs(obs - ranks.mean() * k) - 1e-12:
            extreme += 1
    return extreme / total


def mwu_report(x, y) -> dict:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    out = {
        "n_x": int(len(x)),
        "n_y": int(len(y)),
        "median_x": float(np.median(x)) if len(x) else np.nan,
        "median_y": float(np.median(y)) if len(y) else np.nan,
        "mean_x": float(np.mean(x)) if len(x) else np.nan,
        "mean_y": float(np.mean(y)) if len(y) else np.nan,
        "delta_median_x_minus_y": float(np.median(x) - np.median(y)) if len(x) and len(y) else np.nan,
    }
    if len(x) >= 2 and len(y) >= 2:
        u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        out["U"] = float(u)
        out["p_mwu"] = float(p)
        pe = exact_mwu_p(x, y)
        if pe is not None:
            out["p_exact"] = float(pe)
            out["p"] = float(pe)
        else:
            out["p"] = float(p)
    else:
        out["p"] = np.nan
    return out


def spearman_report(a, b) -> dict:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    out = {"n": int(len(a))}
    if len(a) >= 4:
        rho, p = stats.spearmanr(a, b)
        out["rho"] = float(rho)
        out["p"] = float(p)
    else:
        out["rho"] = np.nan
        out["p"] = np.nan
    return out


def partial_spearman(a, b, cohort: pd.Series) -> dict:
    """Spearman on residuals after OLS on cohort dummies."""
    import statsmodels.api as sm

    a = np.asarray(a, float)
    b = np.asarray(b, float)
    dummies = pd.get_dummies(cohort, drop_first=True).astype(float)
    m = np.isfinite(a) & np.isfinite(b) & dummies.notna().all(axis=1).to_numpy()
    if m.sum() < 6 or dummies.shape[1] == 0:
        return {"n": int(m.sum()), "rho": np.nan, "p": np.nan, "note": "too few or one cohort"}
    X = sm.add_constant(dummies.to_numpy()[m])
    ra = sm.OLS(a[m], X).fit().resid
    rb = sm.OLS(b[m], X).fit().resid
    rho, p = stats.spearmanr(ra, rb)
    return {"n": int(m.sum()), "rho": float(rho), "p": float(p)}


def ols_mpr_cohort(score: pd.Series, mpr: pd.Series, cohort: pd.Series) -> dict:
    import statsmodels.formula.api as smf

    df = pd.DataFrame({"score": score, "mpr": mpr, "cohort": cohort}).dropna()
    if df["mpr"].nunique() < 2 or len(df) < 6:
        return {"n": int(len(df)), "note": "too few"}
    fit = smf.ols('score ~ C(mpr, Treatment("MPR")) + C(cohort)', data=df).fit()
    key = [k for k in fit.params.index if "mpr" in k.lower()][0]
    return {
        "n": int(len(df)),
        "coef_NMPR_vs_MPR": float(fit.params[key]),
        "se": float(fit.bse[key]),
        "p": float(fit.pvalues[key]),
        "r2": float(fit.rsquared),
        "formula": 'score ~ C(mpr, Treatment("MPR")) + C(cohort)',
    }


# ---------------------------------------------------------------------------
# Load each series onto a shared gene list
# ---------------------------------------------------------------------------

def inventory_gse205335() -> dict:
    p = DATA / "GSE205335" / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
    lin = df["lineage.total"].value_counts().to_dict()
    sub = df["lineage.sub"].value_counts().to_dict()
    return {
        "accession": "GSE205335",
        "n_cells_identity": int(len(df)),
        "n_patients_orig_ident": int(df["orig.ident"].nunique()),
        "lineage.total": lin,
        "lineage.sub": sub,
        "has_epithelium": bool(lin.get("Epithelial cells", 0) > 0 or sub.get("Malignant cells", 0) > 0),
        "included_in_harmony": False,
        "why_out": (
            "Author epithelium and malignant cells are present "
            f"(Epithelial cells={lin.get('Epithelial cells', 0)}, "
            f"Malignant cells={sub.get('Malignant cells', 0)}), but the series is "
            "palliative / advanced ICI with RECIST, not neoadjuvant MPR. "
            "UMI is an R dgCMatrix RDS (not MTX). Excluded to avoid an endpoint swap "
            "and an extra RDS parse. CellIdentity only was downloaded."
        ),
    }


def load_gse207422(shared_genes: list[str]) -> tuple[np.ndarray, pd.DataFrame, list[str]]:
    xlsx = pd.read_excel(DATA / "GSE207422" / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    xlsx = xlsx.dropna(subset=["Sample"]).copy()
    xlsx = xlsx[xlsx["Sample"].astype(str).str.startswith("BD_")]
    meta_s = xlsx.set_index("Sample")
    path = DATA / "GSE207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
    core = [g for g in core_panel() if g in set(shared_genes)]
    X, genes, nUMI, barcodes = gse207422_stats_and_extract(path, core, N_HVG)
    # restrict to shared (already from this matrix); caller will reindex later
    rows = []
    for bc in barcodes:
        sample = "_".join(bc.split("_")[:2])
        r = meta_s.loc[sample]
        rows.append(
            {
                "barcode": bc,
                "sample_id": sample,
                "patient_id": str(r["Patient"]),
                "dataset": "GSE207422",
                "cohort": "GSE207422",
                "mpr": map_mpr(r["Pathologic Response"]),
                "pathologic_response": str(r["Pathologic Response"]),
                "histology": map_histo(r["Pathology"]),
                "is_post": "Post" in str(r["Resource"]),
                "resource": str(r["Resource"]),
                "nUMI": np.nan,  # filled below
                "author_lineage": "",
            }
        )
    obs = pd.DataFrame(rows)
    obs["nUMI"] = nUMI
    return X, obs, genes


def load_gse241934(tag: str, shared: set[str], gene_order: list[str] | None) -> tuple[np.ndarray, pd.DataFrame, list[str]]:
    if tag == "IIT":
        feat = DATA / "GSE241934" / "GSE241934_IIT_features.tsv.gz"
        bc_p = DATA / "GSE241934" / "GSE241934_IIT_barcodes.tsv.gz"
        meta_p = DATA / "GSE241934" / "GSE241934_IIT_Meta.txt.gz"
        mtx = DATA / "GSE241934" / "GSE241934_IIT_Matrix.mtx.gz"
        dataset = "GSE241934_IIT"
    else:
        feat = DATA / "GSE241934" / "GSE241934_RWC_features.tsv.gz"
        bc_p = DATA / "GSE241934" / "GSE241934_RWC_barcodes.tsv.gz"
        meta_p = DATA / "GSE241934" / "GSE241934_Real_Meta.txt.gz"
        mtx = DATA / "GSE241934" / "GSE241934_Real_Matrix.mtx.gz"
        dataset = "GSE241934_REAL"
    with gzip.open(bc_p, "rt") as f:
        barcodes = [ln.rstrip("\n") for ln in f]
    md = pd.read_csv(meta_p, sep="\t", dtype=str, low_memory=False)
    md = md.set_index("cellID", drop=False)
    md = md.reindex(barcodes)
    wanted = set(gene_order) if gene_order is not None else (shared & set(core_panel()))
    keep = _gene_rows_1based(feat, wanted)
    X, genes = stream_mtx(mtx, keep, len(barcodes))
    obs = pd.DataFrame(
        {
            "barcode": barcodes,
            "sample_id": md["sampleID"].astype(str).to_numpy(),
            "patient_id": md["sampleID"].astype(str).to_numpy(),
            "dataset": dataset,
            "cohort": "GSE241934",
            "mpr": [map_mpr(x) for x in md["Pathological Response"]],
            "pathologic_response": md["Pathological Response"].astype(str).to_numpy(),
            "histology": [map_histo(x) for x in md["Histology"]],
            "is_post": True,
            "resource": "post-tx tumor",
            "nUMI": pd.to_numeric(md["nCount_RNA"], errors="coerce").to_numpy(),
            "author_lineage": md["major.cell.type"].astype(str).to_numpy(),
        }
    )
    return X, obs, genes


def load_gse291670(shared: set[str], gene_order: list[str]) -> tuple[np.ndarray, pd.DataFrame, list[str]]:
    raw = DATA / "GSE291670" / "raw"
    blocks_X = []
    blocks_obs = []
    genes_ref = None
    for mtx in sorted(raw.glob("*_matrix.mtx.gz")):
        stem = mtx.name.replace("_matrix.mtx.gz", "")
        # GSM8839599_MPR-1
        gsm, title = stem.split("_", 1)
        feat = raw / f"{stem}_features.tsv.gz"
        bc_p = raw / f"{stem}_barcodes.tsv.gz"
        with gzip.open(bc_p, "rt") as f:
            barcodes = [ln.rstrip("\n") for ln in f]
        keep = _gene_rows_1based(feat, set(gene_order))
        X, genes = stream_mtx(mtx, keep, len(barcodes))
        # nUMI from full MTX (small)
        import scipy.io

        full = scipy.io.mmread(gzip.open(mtx, "rb")).tocsc()
        nUMI = np.asarray(full.sum(axis=0)).ravel().astype(np.float32)
        del full
        mpr = "MPR" if title.upper().startswith("MPR") else "NMPR"
        obs = pd.DataFrame(
            {
                "barcode": [f"{title}_{b}" for b in barcodes],
                "sample_id": title,
                "patient_id": title,
                "dataset": "GSE291670",
                "cohort": "GSE291670",
                "mpr": mpr,
                "pathologic_response": title,
                "histology": "unknown",
                "is_post": True,
                "resource": "post-tx tumor (snRNA)",
                "nUMI": nUMI,
                "author_lineage": "",
            }
        )
        if genes_ref is None:
            genes_ref = genes
        else:
            # align columns
            idx = {g: i for i, g in enumerate(genes)}
            X = X[:, [idx[g] for g in genes_ref]]
        blocks_X.append(X)
        blocks_obs.append(obs)
        print(f"  GSE291670 {title} cells={len(barcodes)} mpr={mpr}", flush=True)
    return np.vstack(blocks_X), pd.concat(blocks_obs, ignore_index=True), genes_ref


def reindex_X(X: np.ndarray, have: list[str], want: list[str]) -> np.ndarray:
    idx = {g: i for i, g in enumerate(have)}
    out = np.zeros((X.shape[0], len(want)), dtype=np.float32)
    for j, g in enumerate(want):
        i = idx.get(g)
        if i is not None:
            out[:, j] = X[:, i]
    return out


# ---------------------------------------------------------------------------
# Lineage / Harmony
# ---------------------------------------------------------------------------

def log_cp10k(X: np.ndarray, nUMI: np.ndarray) -> np.ndarray:
    lib = np.maximum(nUMI.astype(np.float64), 1.0)
    return np.log1p((X.astype(np.float64) / lib[:, None]) * 1e4).astype(np.float32)


def gene_index(genes: list[str]) -> dict[str, int]:
    return {g: i for i, g in enumerate(genes)}


def mean_genes(L: np.ndarray, genes: list[str], names: list[str]) -> np.ndarray:
    gi = gene_index(genes)
    cols = [gi[g] for g in names if g in gi]
    if not cols:
        return np.zeros(L.shape[0], dtype=np.float32)
    return L[:, cols].mean(axis=1)


def assign_marker_lineage(L: np.ndarray, genes: list[str]) -> np.ndarray:
    scores = {k: mean_genes(L, genes, v) for k, v in LINEAGE_SCORES.items()}
    names = list(scores)
    M = np.vstack([scores[k] for k in names])
    arg = M.argmax(axis=0)
    # require the winner to beat immune/epi sanity a bit
    lab = np.array(names, dtype=object)[arg]
    ptprc = mean_genes(L, genes, ["PTPRC"])
    epi = scores["epithelial"]
    # if PTPRC high and epi not winning by much, keep immune winner
    lab[(epi > ptprc + 0.15) & (lab != "epithelial") & (epi >= np.quantile(epi, 0.80))] = "epithelial"
    return lab.astype(str)


def malignant_like(L: np.ndarray, genes: list[str], lineage: np.ndarray) -> np.ndarray:
    nl = mean_genes(L, genes, NORMAL_LUNG)
    epi = lineage == "epithelial"
    # low normal-lung among epithelial: below 75th percentile of epithelial normal-lung
    out = np.zeros(len(lineage), dtype=bool)
    if epi.sum() >= 20:
        thr = np.quantile(nl[epi], 0.75)
        out[epi & (nl <= thr)] = True
    else:
        out[epi] = True
    return out


def run_harmony(Z: np.ndarray, batch: list[str]) -> np.ndarray:
    import harmonypy as hm

    meta = pd.DataFrame({"batch": batch})
    # harmonypy wants cells × PCs
    ho = hm.run_harmony(Z, meta, ["batch"], max_iter_harmony=20)
    Zc = np.ascontiguousarray(np.array(ho.Z_corr).T)
    if Zc.shape != Z.shape:
        # some versions return PCs × cells already
        if Zc.T.shape == Z.shape:
            Zc = np.ascontiguousarray(np.array(ho.Z_corr))
        else:
            raise RuntimeError(f"Harmony shape {Zc.shape} vs PCA {Z.shape}")
    return Zc.astype(np.float32)


def leiden_umap(Z: np.ndarray, n_neighbors: int = 15, res: float = 0.8):
    import igraph as ig
    import leidenalg
    import umap
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
    nn.fit(Z)
    dist, ind = nn.kneighbors()
    edges = []
    w = []
    n = Z.shape[0]
    for i in range(n):
        for d, j in zip(dist[i, 1:], ind[i, 1:]):
            edges.append((i, int(j)))
            w.append(float(1.0 / (d + 1e-6)))
    g = ig.Graph(n=n, edges=edges, directed=False)
    g.es["weight"] = w
    part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, weights="weight", resolution_parameter=res, seed=7)
    cl = np.array(part.membership)
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.3, metric="euclidean", random_state=7)
    uv = reducer.fit_transform(Z)
    return cl, uv


def label_clusters(cl: np.ndarray, marker_lin: np.ndarray) -> dict[int, str]:
    out = {}
    for c in np.unique(cl):
        sub = marker_lin[cl == c]
        if len(sub) == 0:
            out[int(c)] = "other"
            continue
        vals, cnt = np.unique(sub, return_counts=True)
        out[int(c)] = str(vals[cnt.argmax()])
    return out


# ---------------------------------------------------------------------------
# Patient table + figures
# ---------------------------------------------------------------------------

def patient_table(obs: pd.DataFrame, L: np.ndarray, genes: list[str]) -> pd.DataFrame:
    gi = gene_index(genes)
    rows = []
    for (cohort, sid), sub in obs.groupby(["cohort", "sample_id"], sort=False):
        idx = sub.index.to_numpy()
        n = len(sub)
        mal = sub["is_malignant"].to_numpy()
        tnk = sub["is_tnk"].to_numpy()
        rec = {
            "cohort": cohort,
            "dataset": sub["dataset"].iloc[0],
            "sample_id": sid,
            "patient_id": sub["patient_id"].iloc[0],
            "mpr": sub["mpr"].iloc[0],
            "pathologic_response": sub["pathologic_response"].iloc[0],
            "histology": sub["histology"].iloc[0],
            "is_post": bool(sub["is_post"].iloc[0]),
            "n_cells": n,
            "n_malignant": int(mal.sum()),
            "n_tnk": int(tnk.sum()),
            "frac_tnk": float(tnk.mean()) if n else np.nan,
            "frac_malignant": float(mal.mean()) if n else np.nan,
        }
        for gene in TARGETS:
            j = gi.get(gene)
            if j is None:
                rec[f"{gene}_mean_log1p_cp10k"] = np.nan
                rec[f"{gene}_pct_pos"] = np.nan
                continue
            if mal.sum() == 0:
                rec[f"{gene}_mean_log1p_cp10k"] = np.nan
                rec[f"{gene}_pct_pos"] = np.nan
            else:
                rec[f"{gene}_mean_log1p_cp10k"] = float(L[idx[mal], j].mean())
                # pct pos from raw-ish: log1p(CP10k)>0 iff UMI>0
                rec[f"{gene}_pct_pos"] = float((L[idx[mal], j] > 0).mean())
        rows.append(rec)
    return pd.DataFrame(rows)


def eligible(per: pd.DataFrame, gene: str) -> pd.DataFrame:
    a = per[
        per["is_post"]
        & per["mpr"].isin(["MPR", "NMPR"])
        & (per["n_malignant"] >= MIN_MAL)
        & per[f"{gene}_mean_log1p_cp10k"].notna()
    ].copy()
    return a


def make_figures(obs_s: pd.DataFrame, umap_xy: np.ndarray, per: pd.DataFrame, stats_d: dict) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 9, "figure.dpi": 140, "savefig.bbox": "tight"})

    # 1. UMAP cohort + lineage
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 4.0))
    cohorts = sorted(obs_s["dataset"].unique())
    cmap_c = plt.get_cmap("tab10")
    for i, c in enumerate(cohorts):
        m = obs_s["dataset"].to_numpy() == c
        ax[0].scatter(umap_xy[m, 0], umap_xy[m, 1], s=2, alpha=0.35, label=c, color=cmap_c(i), linewidths=0)
    ax[0].legend(markerscale=4, fontsize=7, loc="best")
    ax[0].set_title("Harmony UMAP — dataset")
    ax[0].set_xlabel("UMAP1")
    ax[0].set_ylabel("UMAP2")
    lins = ["epithelial", "t", "nk", "b", "myeloid", "endothelial", "fibroblast", "mast"]
    cmap_l = plt.get_cmap("tab20")
    lin = obs_s["harmony_lineage"].astype(str).to_numpy()
    for i, c in enumerate(lins):
        m = lin == c
        if m.any():
            ax[1].scatter(umap_xy[m, 0], umap_xy[m, 1], s=2, alpha=0.35, label=c, color=cmap_l(i), linewidths=0)
    ax[1].legend(markerscale=4, fontsize=7, loc="best")
    ax[1].set_title("Harmony UMAP — transferred lineage")
    ax[1].set_xlabel("UMAP1")
    fig.savefig(FIG / "fig1_umap_harmony_dataset_lineage.png")
    plt.close()

    # 2. MPR boxplots
    fig, ax = plt.subplots(1, 2, figsize=(8.8, 4.2))
    colors = {"MPR": "#4C78A8", "NMPR": "#F58518"}
    for j, gene in enumerate(TARGETS):
        a = eligible(per, gene)
        data, labs, cols = [], [], []
        for grp in ["MPR", "NMPR"]:
            v = a.loc[a.mpr == grp, f"{gene}_mean_log1p_cp10k"].to_numpy()
            data.append(v)
            labs.append(f"{grp}\nn={len(v)}")
            cols.append(colors[grp])
        bp = ax[j].boxplot(data, tick_labels=labs, patch_artist=True, widths=0.55)
        for patch, c in zip(bp["boxes"], cols):
            patch.set_facecolor(c)
            patch.set_alpha(0.7)
        # jitter
        for i, v in enumerate(data, start=1):
            ax[j].scatter(np.full(len(v), i) + RNG.uniform(-0.08, 0.08, len(v)), v, s=14, c="k", alpha=0.7, zorder=3)
        st = stats_d.get(f"{gene}_mpr", {})
        ax[j].set_title(f"{gene} malignant-like vs MPR\nΔ={st.get('delta_median_x_minus_y', float('nan')):.3f} p={st.get('p', float('nan')):.3g}")
        ax[j].set_ylabel("mean log1p(CP10k) in malignant-like")
    fig.savefig(FIG / "fig2_malignant_tacstd2_cldn4_by_mpr.png")
    plt.close()

    # 3. vs T/NK
    fig, ax = plt.subplots(1, 2, figsize=(8.8, 4.2))
    cmap = {c: plt.get_cmap("tab10")(i) for i, c in enumerate(sorted(per["cohort"].unique()))}
    for j, gene in enumerate(TARGETS):
        a = eligible(per, gene)
        a = a[a["n_tnk"] >= MIN_TNK]
        for c, sub in a.groupby("cohort"):
            ax[j].scatter(
                sub["frac_tnk"],
                sub[f"{gene}_mean_log1p_cp10k"],
                s=28,
                alpha=0.85,
                label=c,
                color=cmap[c],
            )
        st = stats_d.get(f"{gene}_tnk", {})
        ax[j].set_title(f"{gene} vs T/NK\nρ={st.get('rho', float('nan')):.2f} p={st.get('p', float('nan')):.3g} n={st.get('n', 0)}")
        ax[j].set_xlabel("T/NK fraction (joint embedding)")
        ax[j].set_ylabel(f"malignant-like {gene} mean log1p(CP10k)")
        ax[j].legend(fontsize=7)
    fig.savefig(FIG / "fig3_malignant_score_vs_tnk.png")
    plt.close()

    # 4. histology
    fig, ax = plt.subplots(1, 2, figsize=(8.8, 4.2))
    for j, gene in enumerate(TARGETS):
        a = per[per["is_post"] & (per["n_malignant"] >= MIN_MAL) & per["histology"].isin(["LUAD", "LUSC"])].copy()
        data, labs = [], []
        for h in ["LUAD", "LUSC"]:
            v = a.loc[a.histology == h, f"{gene}_mean_log1p_cp10k"].to_numpy()
            data.append(v)
            labs.append(f"{h}\nn={len(v)}")
        if all(len(v) for v in data):
            bp = ax[j].boxplot(data, tick_labels=labs, patch_artist=True, widths=0.55)
            for patch, c in zip(bp["boxes"], ["#54A24B", "#B279A2"]):
                patch.set_facecolor(c)
                patch.set_alpha(0.7)
            for i, v in enumerate(data, start=1):
                ax[j].scatter(np.full(len(v), i) + RNG.uniform(-0.08, 0.08, len(v)), v, s=14, c="k", alpha=0.7)
        st = stats_d.get(f"{gene}_histo", {})
        ax[j].set_title(f"{gene} LUAD vs LUSC\np={st.get('p', float('nan')):.3g}")
        ax[j].set_ylabel("mean log1p(CP10k) in malignant-like")
    fig.savefig(FIG / "fig4_luad_vs_lusc.png")
    plt.close()


def write_finding(inv: dict, per: pd.DataFrame, stats_d: dict, n_cells: int, n_embed: int) -> None:
    n_pat = int(per["patient_id"].nunique())
    n_samp = int(len(per))
    post = per[per["is_post"] & per["mpr"].isin(["MPR", "NMPR"])]
    lines = []
    lines.append("# Harmony neoadjuvant lung scRNA — FINDING")
    lines.append("")
    lines.append("**Public data only. Additive to User A3. The GSE207422 CopyKAT slide was not re-run.**")
    lines.append("")
    lines.append("## 一句话结论 / TL;DR")
    t = stats_d.get("TACSTD2_mpr", {})
    c = stats_d.get("CLDN4_mpr", {})
    tb = stats_d.get("TACSTD2_tnk", {})
    cb = stats_d.get("CLDN4_tnk", {})
    lines.append(
        f"Joint Harmony object: **{n_cells:,} cells** from **{n_pat} patients** "
        f"({n_samp} samples) across GSE207422 + GSE241934 + GSE291670 "
        f"(embedding used {n_embed:,} stratified cells). "
        f"Malignant-like TACSTD2 NMPR vs MPR: median {t.get('median_x', float('nan')):.3f} vs "
        f"{t.get('median_y', float('nan')):.3f} (n={t.get('n_x', 0)} vs {t.get('n_y', 0)}, "
        f"p={t.get('p', float('nan')):.3g}). "
        f"CLDN4: {c.get('median_x', float('nan')):.3f} vs {c.get('median_y', float('nan')):.3f}, "
        f"p={c.get('p', float('nan')):.3g}. "
        f"TACSTD2 vs T/NK ρ={tb.get('rho', float('nan')):.2f} (n={tb.get('n', 0)}, p={tb.get('p', float('nan')):.3g}); "
        f"CLDN4 ρ={cb.get('rho', float('nan')):.2f} (p={cb.get('p', float('nan')):.3g}). "
        "This does not retract or replace User A3."
    )
    lines.append("")
    lines.append("## n")
    lines.append("| | Count |")
    lines.append("|---|---|")
    lines.append(f"| Cells in deposited matrices (included series) | {n_cells:,} |")
    lines.append(f"| Cells in Harmony embedding (stratified cap {MAX_PER_SAMPLE}/sample) | {n_embed:,} |")
    lines.append(f"| Patients | {n_pat} |")
    lines.append(f"| Post-tx samples with MPR label | {len(post)} (MPR {(post.mpr=='MPR').sum()}, NMPR {(post.mpr=='NMPR').sum()}) |")
    for coh, sub in per.groupby("cohort"):
        lines.append(f"| {coh} samples / cells | {len(sub)} / {int(sub.n_cells.sum()):,} |")
    lines.append("")
    lines.append("## Left out")
    g = inv["GSE205335"]
    lines.append(f"- **GSE205335** — {g['why_out']}")
    lines.append("")
    lines.append("## Primary (pre-specified)")
    lines.append("| Test | n | Result | p |")
    lines.append("|---|---|---|---|")
    for key, lab in [
        ("TACSTD2_mpr", "TACSTD2 NMPR vs MPR (median log1p CP10k)"),
        ("CLDN4_mpr", "CLDN4 NMPR vs MPR"),
        ("TACSTD2_mpr_ols", "TACSTD2 OLS NMPR vs MPR + cohort"),
        ("CLDN4_mpr_ols", "CLDN4 OLS NMPR vs MPR + cohort"),
        ("TACSTD2_tnk", "TACSTD2 vs T/NK Spearman"),
        ("CLDN4_tnk", "CLDN4 vs T/NK Spearman"),
        ("TACSTD2_tnk_partial", "TACSTD2 vs T/NK partial Spearman (cohort)"),
        ("CLDN4_tnk_partial", "CLDN4 vs T/NK partial Spearman (cohort)"),
        ("TACSTD2_histo", "TACSTD2 LUAD vs LUSC"),
        ("CLDN4_histo", "CLDN4 LUAD vs LUSC"),
    ]:
        st = stats_d.get(key, {})
        if "rho" in st:
            res = f"ρ={st.get('rho', float('nan')):.3f}"
            p = st.get("p", float("nan"))
            n = st.get("n", "")
        elif "coef_NMPR_vs_MPR" in st:
            res = f"β={st.get('coef_NMPR_vs_MPR', float('nan')):.3f}"
            p = st.get("p", float("nan"))
            n = st.get("n", "")
        else:
            res = f"med NMPR {st.get('median_x', float('nan')):.3f} vs MPR {st.get('median_y', float('nan')):.3f}"
            p = st.get("p", float("nan"))
            n = f"{st.get('n_x', '')} vs {st.get('n_y', '')}"
        lines.append(f"| {lab} | {n} | {res} | {p if p==p else 'NA'} |")
    lines.append("")
    lines.append("## Read this as extra, not as a copy of the given slide")
    lines.append("- Cell-level p-values are not reported (pseudoreplication).")
    lines.append("- GSE207422 CopyKAT author IDs are still not on GEO; malignant-like here is Harmony epithelium minus a high normal-lung score (plus author Epi as sensitivity on GSE241934).")
    lines.append("- GSE291670 is nuclear RNA; Harmony batch = dataset.")
    lines.append("- GSE241934 has no LUSC; LUAD vs LUSC is powered mainly by GSE207422 Squamous vs Adeno.")
    lines.append("")
    lines.append("## Files")
    lines.append("- `figures/fig1_umap_harmony_dataset_lineage.png`")
    lines.append("- `figures/fig2_malignant_tacstd2_cldn4_by_mpr.png`")
    lines.append("- `figures/fig3_malignant_score_vs_tnk.png`")
    lines.append("- `figures/fig4_luad_vs_lusc.png`")
    lines.append("- `tables/per_sample.tsv`, `tables/stats.json`, `tables/inventory.json`")
    (ROOT / "FINDING.md").write_text("\n".join(lines) + "\n")


def _dir_phrase(delta: float, p: float, higher: str, lower: str) -> str:
    if not np.isfinite(p):
        return "was not testable"
    if p >= 0.05:
        if np.isfinite(delta) and delta > 0:
            return f"was higher in {higher} than {lower} but not significant"
        if np.isfinite(delta) and delta < 0:
            return f"was higher in {lower} than {higher} but not significant"
        return "did not differ significantly"
    if np.isfinite(delta) and delta > 0:
        return f"was higher in {higher} than {lower}"
    if np.isfinite(delta) and delta < 0:
        return f"was higher in {lower} than {higher}"
    return "differed"


def write_snippet(stats_d: dict, n_cells: int, n_pat: int, n_embed: int) -> None:
    t = stats_d.get("TACSTD2_mpr", {})
    c = stats_d.get("CLDN4_mpr", {})
    tb = stats_d.get("TACSTD2_tnk", {})
    cb = stats_d.get("CLDN4_tnk", {})
    to = stats_d.get("TACSTD2_mpr_ols", {})
    tp = stats_d.get("TACSTD2_tnk_partial", {})
    th = stats_d.get("TACSTD2_histo", {})
    ch = stats_d.get("CLDN4_histo", {})
    t_phrase = _dir_phrase(t.get("delta_median_x_minus_y", float("nan")), t.get("p", float("nan")), "NMPR", "MPR")
    c_phrase = _dir_phrase(c.get("delta_median_x_minus_y", float("nan")), c.get("p", float("nan")), "NMPR", "MPR")
    txt = f"""# Paper snippet — joint Harmony neoadjuvant scRNA (additive)

*Placement: extra panel after User A3. Do **not** replace the GSE207422 CopyKAT slide.*

Public processed UMI matrices from three neoadjuvant NSCLC series with epithelium and immune cells (GSE207422, GSE241934, GSE291670) were intersected on gene symbols and integrated with Harmony on a shared space (batch = dataset). GSE205335 has author malignant epithelium but was left out: it is palliative ICI with RECIST, not neoadjuvant MPR. The joint object contains **{n_cells:,} cells** from **{n_pat} patients** (Harmony/UMAP used a sample-stratified subset of {n_embed:,} cells). Malignant-like cells were Harmony epithelial clusters after dropping a high normal-lung program. Patient-level malignant-like *TACSTD2* {t_phrase} (median {t.get('median_x', float('nan')):.2f} vs {t.get('median_y', float('nan')):.2f}; n={t.get('n_x', 0)} vs {t.get('n_y', 0)}; p={t.get('p', float('nan')):.3g}). *CLDN4* {c_phrase} (median {c.get('median_x', float('nan')):.2f} vs {c.get('median_y', float('nan')):.2f}; p={c.get('p', float('nan')):.3g}). Cohort-adjusted OLS for *TACSTD2* (NMPR vs MPR + cohort) gave β={to.get('coef_NMPR_vs_MPR', float('nan')):.3f} (p={to.get('p', float('nan')):.3g}). Per-patient malignant-like *TACSTD2* vs joint-embedding T/NK fraction was ρ={tb.get('rho', float('nan')):.2f} (n={tb.get('n', 0)}, p={tb.get('p', float('nan')):.3g}); partial Spearman after cohort residualization ρ={tp.get('rho', float('nan')):.2f} (p={tp.get('p', float('nan')):.3g}). *CLDN4* vs T/NK: ρ={cb.get('rho', float('nan')):.2f} (p={cb.get('p', float('nan')):.3g}). Where histology was labeled, malignant-like *TACSTD2* LUSC vs LUAD p={th.get('p', float('nan')):.3g} (n_LUSC={th.get('n_x', 0)}, n_LUAD={th.get('n_y', 0)}); *CLDN4* p={ch.get('p', float('nan')):.3g}. These public joint-embedding tests are additive and do not replace User A3.

**Caption (Fig. Harmony-neoadj).** Harmony UMAP of public neoadjuvant NSCLC scRNA (GSE207422, GSE241934 IIT+real-world, GSE291670). Malignant-like *TACSTD2*/*CLDN4* mean log1p(CP10k) by MPR/NMPR (patient unit; pCR counted as MPR) and versus T/NK fraction in the joint embedding. Cohort is a covariate in the OLS / partial Spearman. GSE205335 omitted (RECIST, not MPR).
"""
    (ROOT / "paper_snippet.md").write_text(txt)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    warnings.filterwarnings("ignore", category=UserWarning)
    inv205 = inventory_gse205335()
    print("GSE205335", json.dumps({k: inv205[k] for k in ("n_cells_identity", "has_epithelium", "included_in_harmony")}, indent=2), flush=True)

    # gene universes
    def tsv_genes(path):
        with gzip.open(path, "rt") as f:
            f.readline()
            return [ln.split("\t", 1)[0] for ln in f]

    g207 = set(tsv_genes(DATA / "GSE207422" / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"))
    g241 = set(_feat_symbols(DATA / "GSE241934" / "GSE241934_IIT_features.tsv.gz"))
    g291 = set(_feat_symbols(next((DATA / "GSE291670" / "raw").glob("*MPR-1_features.tsv.gz"))))
    shared = sorted(g207 & g241 & g291)
    print(f"shared symbols={len(shared)} TACSTD2={('TACSTD2' in shared)} CLDN4={('CLDN4' in shared)}", flush=True)
    if "TACSTD2" not in shared or "CLDN4" not in shared:
        raise SystemExit("TACSTD2/CLDN4 missing from intersection")

    # Load GSE207422 first (defines HVG list)
    X207, obs207, genes207 = load_gse207422(shared)
    gene_order = [g for g in genes207 if g in set(shared)]
    if "TACSTD2" not in gene_order:
        gene_order = ["TACSTD2"] + gene_order
    if "CLDN4" not in gene_order:
        gene_order = ["CLDN4"] + gene_order
    # unique preserve
    seen = set()
    go = []
    for g in gene_order:
        if g not in seen:
            seen.add(g)
            go.append(g)
    gene_order = go
    print(f"gene_order n={len(gene_order)}", flush=True)
    X207 = reindex_X(X207, genes207, gene_order)

    X_iit, obs_iit, g_iit = load_gse241934("IIT", set(shared), gene_order)
    X_iit = reindex_X(X_iit, g_iit, gene_order)
    X_real, obs_real, g_real = load_gse241934("REAL", set(shared), gene_order)
    X_real = reindex_X(X_real, g_real, gene_order)
    X291, obs291, g291e = load_gse291670(set(shared), gene_order)
    X291 = reindex_X(X291, g291e, gene_order)

    X = np.vstack([X207, X_iit, X_real, X291])
    obs = pd.concat([obs207, obs_iit, obs_real, obs291], ignore_index=True)
    del X207, X_iit, X_real, X291
    print(f"concat cells={len(obs)} genes={len(gene_order)}", flush=True)

    nUMI = obs["nUMI"].to_numpy(dtype=np.float64)
    nUMI = np.where(np.isfinite(nUMI) & (nUMI > 0), nUMI, X.sum(axis=1) + 1)
    L = log_cp10k(X, nUMI)
    del X

    marker_lin = assign_marker_lineage(L, gene_order)
    obs["marker_lineage"] = marker_lin
    # author override where present
    al = obs["author_lineage"].astype(str)
    map_auth = {"Epi": "epithelial", "T": "t", "NK": "nk", "B": "b", "Myeloid": "myeloid", "Endo": "endothelial", "Fibro": "fibroblast", "Mast": "mast"}
    auth_lin = al.map(map_auth)
    obs["native_lineage"] = np.where(auth_lin.notna() & (al != "") & (al != "nan"), auth_lin, marker_lin)

    # subsample for Harmony
    take = []
    for sid, sub in obs.groupby("sample_id"):
        k = min(MAX_PER_SAMPLE, len(sub))
        take.extend(sub.sample(k, random_state=7).index.to_list())
    take = np.array(sorted(take))
    print(f"Harmony subsample {len(take)} / {len(obs)}", flush=True)

    Ls = L[take]
    # z-score genes on subsample, apply to all
    mu = Ls.mean(axis=0)
    sd = Ls.std(axis=0)
    sd[sd < 1e-6] = 1.0
    Zall = (L - mu) / sd
    from sklearn.decomposition import PCA

    pca = PCA(n_components=min(N_PCS, Zall.shape[1] - 1), random_state=7)
    pcs_s = pca.fit_transform(Zall[take])
    print("PCA done", pcs_s.shape, flush=True)
    Z_h = run_harmony(pcs_s, obs.loc[take, "dataset"].astype(str).tolist())
    print("Harmony done", Z_h.shape, flush=True)
    cl, uv = leiden_umap(Z_h)
    cl_map = label_clusters(cl, obs.loc[take, "native_lineage"].to_numpy())
    harm_lin_s = np.array([cl_map[int(c)] for c in cl])
    obs_s = obs.loc[take].copy()
    obs_s["harmony_lineage"] = harm_lin_s
    obs_s["leiden"] = cl

    # transfer to all cells: kNN in PCA space of subsample
    from sklearn.neighbors import NearestNeighbors

    pcs_all = pca.transform(Zall)
    nn = NearestNeighbors(n_neighbors=1)
    nn.fit(pcs_s)
    _, ind = nn.kneighbors(pcs_all)
    obs["harmony_lineage"] = harm_lin_s[ind[:, 0]]
    obs["leiden"] = cl[ind[:, 0]]

    # malignant-like / T/NK from joint embedding labels
    # refine epithelial with normal-lung on ALL cells
    nl = mean_genes(L, gene_order, NORMAL_LUNG)
    is_epi = obs["harmony_lineage"].to_numpy() == "epithelial"
    is_mal = np.zeros(len(obs), dtype=bool)
    if is_epi.sum() >= 20:
        thr = np.quantile(nl[is_epi], 0.75)
        is_mal[is_epi & (nl <= thr)] = True
    else:
        is_mal[is_epi] = True
    # GSE241934 sensitivity already in native; primary = joint
    obs["is_malignant"] = is_mal
    obs["is_tnk"] = obs["harmony_lineage"].isin(["t", "nk"]).to_numpy()

    print(
        "lineage counts",
        obs["harmony_lineage"].value_counts().to_dict(),
        "n_mal",
        int(is_mal.sum()),
        "n_tnk",
        int(obs["is_tnk"].sum()),
        flush=True,
    )

    per = patient_table(obs, L, gene_order)
    per.to_csv(TAB / "per_sample.tsv", sep="\t", index=False)

    stats_d = {}
    for gene in TARGETS:
        a = eligible(per, gene)
        nmpr = a.loc[a.mpr == "NMPR", f"{gene}_mean_log1p_cp10k"].to_numpy()
        mpr = a.loc[a.mpr == "MPR", f"{gene}_mean_log1p_cp10k"].to_numpy()
        stats_d[f"{gene}_mpr"] = mwu_report(nmpr, mpr)
        stats_d[f"{gene}_mpr_ols"] = ols_mpr_cohort(
            a[f"{gene}_mean_log1p_cp10k"], a["mpr"], a["cohort"]
        )
        b = a[a["n_tnk"] >= MIN_TNK]
        stats_d[f"{gene}_tnk"] = spearman_report(b[f"{gene}_mean_log1p_cp10k"], b["frac_tnk"])
        stats_d[f"{gene}_tnk_partial"] = partial_spearman(b[f"{gene}_mean_log1p_cp10k"], b["frac_tnk"], b["cohort"])
        h = per[per["is_post"] & (per["n_malignant"] >= MIN_MAL) & per["histology"].isin(["LUAD", "LUSC"])]
        # x = LUSC, y = LUAD in mwu_report call order? use LUSC vs LUAD as x vs y for delta LUSC-LUAD
        stats_d[f"{gene}_histo"] = mwu_report(
            h.loc[h.histology == "LUSC", f"{gene}_mean_log1p_cp10k"].to_numpy(),
            h.loc[h.histology == "LUAD", f"{gene}_mean_log1p_cp10k"].to_numpy(),
        )
        # per-cohort MWU
        for coh, sub in a.groupby("cohort"):
            stats_d[f"{gene}_mpr_{coh}"] = mwu_report(
                sub.loc[sub.mpr == "NMPR", f"{gene}_mean_log1p_cp10k"].to_numpy(),
                sub.loc[sub.mpr == "MPR", f"{gene}_mean_log1p_cp10k"].to_numpy(),
            )

    n_cells = int(len(obs))
    n_embed = int(len(take))
    n_pat = int(per["patient_id"].nunique())

    inventory = {
        "GSE205335": inv205,
        "shared_n_genes_intersection": len(shared),
        "harmony_n_genes": len(gene_order),
        "harmony_genes": gene_order,
        "n_cells_included": n_cells,
        "n_cells_embedding": n_embed,
        "n_patients": n_pat,
        "n_samples": int(len(per)),
        "by_dataset": obs.groupby("dataset").size().to_dict(),
        "by_cohort_patients": per.groupby("cohort")["patient_id"].nunique().to_dict(),
        "included": ["GSE207422", "GSE241934", "GSE291670"],
        "excluded": ["GSE205335"],
        "max_per_sample_embed": MAX_PER_SAMPLE,
        "min_malignant": MIN_MAL,
        "min_tnk": MIN_TNK,
    }
    (TAB / "inventory.json").write_text(json.dumps(inventory, indent=2))
    (TAB / "stats.json").write_text(json.dumps(stats_d, indent=2))
    obs.groupby(["dataset", "harmony_lineage"]).size().unstack(fill_value=0).to_csv(TAB / "lineage_by_dataset.tsv", sep="\t")

    make_figures(obs_s.reset_index(drop=True), uv, per, stats_d)
    write_finding(inventory, per, stats_d, n_cells, n_embed)
    write_snippet(stats_d, n_cells, n_pat, n_embed)

    # dropouts
    drop = per[per["is_post"] & per["mpr"].isin(["MPR", "NMPR"]) & (per["n_malignant"] < MIN_MAL)]
    drop.to_csv(TAB / "dropouts_lt20_malignant.tsv", sep="\t", index=False)
    print("STATS", json.dumps(stats_d, indent=2)[:2000], flush=True)
    print("wrote figures + FINDING", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
