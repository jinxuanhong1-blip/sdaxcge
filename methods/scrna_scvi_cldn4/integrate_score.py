#!/usr/bin/env python3
"""scVI (else Harmony) latent on GSE207422, scored by malignant CLDN4.

Additive. TACSTD2 is a companion gene and is never a gate. The given A3
TACSTD2 slice and the marker-only CLDN4 dualhigh slice are not re-argued.

Unit of every test is the 10x sample. Primary clinical unit is the 12
post-treatment patients. Cell-level p-values are not reported.
"""
from __future__ import annotations

import argparse
import json
import traceback
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.sparse import issparse

MIN_UMI = 200
MIN_GENES = 200
MIN_MAL = 10
MIN_TNK = 20
N_HVG = 2000
N_LATENT = 10
MAX_EPOCHS = 40

LINEAGE = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT7", "ELF3", "CDH1"],
    "T/NK": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "NKG7", "GNLY", "KLRD1"],
    "B/Plasma": ["CD79A", "CD79B", "MS4A1", "JCHAIN", "MZB1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA", "FCN1"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
}
NORMAL_LUNG = [
    "SFTPA1",
    "SFTPA2",
    "SFTPB",
    "SFTPC",
    "AGER",
    "SCGB1A1",
    "SCGB3A2",
    "TPPP3",
    "FOXJ1",
]
A3_NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]
CYTO_GENES = ["GZMB", "GZMA", "PRF1", "IFNG", "NKG7"]


def present(genes, columns) -> list[str]:
    return [g for g in genes if g in columns]


def gene_umi(adata: ad.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.zeros(adata.n_obs, dtype=np.float64)
    x = adata[:, gene].X
    if issparse(x):
        x = np.asarray(x.todense()).ravel()
    else:
        x = np.asarray(x).ravel()
    return np.asarray(x, dtype=np.float64).ravel()


def gene_log1p_cp10k(adata: ad.AnnData, gene: str) -> np.ndarray:
    x = gene_umi(adata, gene)
    tot = np.asarray(adata.obs["n_umi"], dtype=np.float64)
    tot[tot <= 0] = np.nan
    return np.log1p(x / tot * 1e4)


def mean_score(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    genes = present(genes, adata.var_names)
    if not genes:
        return np.zeros(adata.n_obs, dtype=np.float64)
    mats = [gene_log1p_cp10k(adata, g) for g in genes]
    return np.mean(np.vstack(mats), axis=0)


def assign_marker_lineage(adata: ad.AnnData) -> None:
    scores = {k: mean_score(adata, gs) for k, gs in LINEAGE.items()}
    score_df = pd.DataFrame(scores, index=adata.obs_names)
    for col in score_df.columns:
        adata.obs[f"score_{col}"] = score_df[col].to_numpy()
    adata.obs["lineage"] = score_df.idxmax(axis=1)
    adata.obs["normal_lung"] = mean_score(adata, NORMAL_LUNG)
    epi = adata.obs["lineage"] == "Epithelial"
    cut = float(adata.obs.loc[epi, "normal_lung"].quantile(0.75)) if epi.any() else 0.0
    adata.uns["normal_lung_cut"] = cut
    adata.obs["malignant_like"] = epi & (adata.obs["normal_lung"] <= cut)
    adata.obs["is_tnk"] = adata.obs["lineage"] == "T/NK"
    a3 = np.zeros(adata.n_obs, dtype=np.int32)
    for g in A3_NORMAL_LUNG:
        a3 += (gene_umi(adata, g) > 0).astype(np.int32)
    adata.obs["a3_malignant"] = epi.to_numpy() & (a3 == 0)
    adata.obs["CLDN4_log1p"] = gene_log1p_cp10k(adata, "CLDN4")
    adata.obs["CLDN4_umi"] = gene_umi(adata, "CLDN4")
    adata.obs["TACSTD2_log1p"] = gene_log1p_cp10k(adata, "TACSTD2")
    adata.obs["TACSTD2_umi"] = gene_umi(adata, "TACSTD2")
    adata.obs["CXCL13_umi"] = gene_umi(adata, "CXCL13")
    adata.obs["is_cxcl13_tnk"] = adata.obs["is_tnk"] & (adata.obs["CXCL13_umi"] >= 1)
    cyto = mean_score(adata, CYTO_GENES)
    adata.obs["cyto_score"] = cyto
    tnk = adata.obs["is_tnk"].to_numpy()
    if tnk.any():
        cyto_cut = float(np.quantile(cyto[tnk], 0.75))
    else:
        cyto_cut = 0.0
    adata.uns["cyto_cut"] = cyto_cut
    adata.obs["is_cyto_high"] = adata.obs["is_tnk"] & (adata.obs["cyto_score"] >= cyto_cut)


def annotate_leiden_clusters(adata: ad.AnnData) -> None:
    rows = []
    for cl, g in adata.obs.groupby("leiden", observed=True):
        rec = {"leiden": str(cl), "n": int(len(g))}
        for lin in LINEAGE:
            rec[f"score_{lin}"] = float(g[f"score_{lin}"].mean())
        rec["score_normal_lung"] = float(g["normal_lung"].mean())
        rec["frac_marker_epi"] = float((g["lineage"] == "Epithelial").mean())
        rec["frac_marker_tnk"] = float((g["lineage"] == "T/NK").mean())
        rec["mean_CLDN4"] = float(g["CLDN4_log1p"].mean())
        rows.append(rec)
    tab = pd.DataFrame(rows)
    score_cols = [c for c in tab.columns if c.startswith("score_") and c != "score_normal_lung"]
    tab["cluster_lineage"] = tab[score_cols].idxmax(axis=1).str.replace("score_", "", regex=False)
    epi_cl = tab["cluster_lineage"] == "Epithelial"
    cl_cut = float(tab.loc[epi_cl, "score_normal_lung"].quantile(0.75)) if epi_cl.any() else 0.0
    tab["cluster_malignant"] = (tab["cluster_lineage"] == "Epithelial") & (
        tab["score_normal_lung"] <= cl_cut
    )
    mapping = tab.set_index("leiden")
    adata.obs["cluster_lineage"] = adata.obs["leiden"].astype(str).map(mapping["cluster_lineage"])
    adata.obs["cluster_malignant"] = (
        adata.obs["leiden"].astype(str).map(mapping["cluster_malignant"]).astype(bool)
    )
    adata.obs["cluster_tnk"] = adata.obs["cluster_lineage"] == "T/NK"
    adata.uns["cluster_annotation"] = tab.to_dict(orient="records")
    adata.uns["cluster_nl_cut"] = cl_cut


def spearman(x, y) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4:
        return float("nan"), float("nan"), n
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), n


def exact_mwu(a, b) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    n_a, n_b = int(len(a)), int(len(b))
    if n_a < 2 or n_b < 2:
        return {
            "n_a": n_a,
            "n_b": n_b,
            "mean_a": float(np.mean(a)) if n_a else np.nan,
            "mean_b": float(np.mean(b)) if n_b else np.nan,
            "delta": np.nan,
            "exact_p": np.nan,
            "note": "too_few_samples",
        }
    method = "exact" if (n_a + n_b) <= 20 else "asymptotic"
    _u, p = stats.mannwhitneyu(a, b, alternative="two-sided", method=method)
    return {
        "n_a": n_a,
        "n_b": n_b,
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "delta": float(np.mean(a) - np.mean(b)),
        "exact_p": float(p),
        "note": "",
    }


def sample_table(adata: ad.AnnData, mal_col: str, tnk_col: str) -> pd.DataFrame:
    rows = []
    for sample, g in adata.obs.groupby("sample", observed=True):
        n = len(g)
        n_mal = int(g[mal_col].sum())
        n_tnk = int(g[tnk_col].sum())
        n_cx = int(g["is_cxcl13_tnk"].sum())
        n_cy = int(g["is_cyto_high"].sum())
        rec = {
            "sample": sample,
            "patient": g["patient"].iloc[0] if "patient" in g else sample,
            "timing": g["timing"].iloc[0] if "timing" in g else "",
            "mpr": g["mpr"].iloc[0] if "mpr" in g else "",
            "histology": g["histology"].iloc[0] if "histology" in g else "",
            "n_cells": n,
            "n_malignant": n_mal,
            "n_tnk": n_tnk,
            "n_cxcl13": n_cx,
            "n_cyto_high": n_cy,
            "frac_tnk": n_tnk / n if n else np.nan,
            "frac_cxcl13": n_cx / n if n else np.nan,
            "frac_cyto_high": n_cy / n if n else np.nan,
        }
        if n_mal >= MIN_MAL:
            mal = g.loc[g[mal_col]]
            rec["mal_CLDN4_mean"] = float(mal["CLDN4_log1p"].mean())
            rec["mal_CLDN4_pct_pos"] = float((mal["CLDN4_umi"] > 0).mean() * 100)
            rec["mal_TACSTD2_mean"] = float(mal["TACSTD2_log1p"].mean())
            rec["mal_TACSTD2_pct_pos"] = float((mal["TACSTD2_umi"] > 0).mean() * 100)
        else:
            rec["mal_CLDN4_mean"] = np.nan
            rec["mal_CLDN4_pct_pos"] = np.nan
            rec["mal_TACSTD2_mean"] = np.nan
            rec["mal_TACSTD2_pct_pos"] = np.nan
        rec["eligible"] = bool(
            n_mal >= MIN_MAL and n_tnk >= MIN_TNK and np.isfinite(rec["mal_CLDN4_mean"])
        )
        rows.append(rec)
    return pd.DataFrame(rows)


def score_block(samp: pd.DataFrame, label: str, xcol: str, ycol: str) -> dict:
    ok = samp.loc[samp["eligible"]].copy()
    r, p, n = spearman(ok[xcol], ok[ycol])
    return {
        "contrast": label,
        "n": n,
        "n_samples_total": int(len(samp)),
        "n_eligible": int(ok["eligible"].sum()) if len(samp) else 0,
        "dropped": samp.loc[~samp["eligible"], "sample"].astype(str).tolist(),
        "rho": r,
        "p": p,
        "negative": bool(np.isfinite(r) and r < 0),
        "x": xcol,
        "y": ycol,
    }


def try_scvi(adata: ad.AnnData, batch_key: str) -> str:
    import scvi
    import torch

    torch.set_num_threads(4)
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
    model = scvi.model.SCVI(
        adata,
        n_hidden=64,
        n_layers=1,
        n_latent=N_LATENT,
        gene_likelihood="nb",
    )
    model.train(
        max_epochs=MAX_EPOCHS,
        early_stopping=True,
        early_stopping_patience=8,
        batch_size=256,
        accelerator="cpu",
        devices=1,
        enable_progress_bar=True,
    )
    adata.obsm["X_int"] = model.get_latent_representation()
    adata.uns["integrator"] = {
        "method": "scVI",
        "scvi_version": getattr(scvi, "__version__", "unknown"),
        "n_latent": N_LATENT,
        "max_epochs": MAX_EPOCHS,
        "n_hvg": int(adata.n_vars),
        "n_cells": int(adata.n_obs),
        "batch_key": batch_key,
    }
    return "scVI"


def try_harmony(adata: ad.AnnData, batch_key: str, note: str) -> str:
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30, svd_solver="arpack")
    try:
        sc.external.pp.harmony_integrate(
            adata, key=batch_key, basis="X_pca", adjusted_basis="X_int"
        )
    except Exception:
        import harmonypy as hm

        ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, batch_key)
        adata.obsm["X_int"] = np.array(ho.Z_corr).T
    adata.uns["integrator"] = {
        "method": "Harmony",
        "n_pcs": 30,
        "n_hvg": int(adata.n_vars),
        "n_cells": int(adata.n_obs),
        "batch_key": batch_key,
        "note": note,
    }
    return "Harmony"


def integrate(adata_full: ad.AnnData, batch_key: str, force_harmony: bool) -> tuple[ad.AnnData, str]:
    if "n_umi" not in adata_full.obs:
        adata_full.obs["n_umi"] = np.asarray(adata_full.X.sum(axis=1)).ravel()
    sc.pp.highly_variable_genes(
        adata_full,
        n_top_genes=N_HVG,
        flavor="seurat_v3",
        batch_key=batch_key,
        layer="counts",
        subset=False,
    )
    hvg = adata_full[:, adata_full.var["highly_variable"]].copy()
    hvg.layers["counts"] = hvg.X.copy()
    method = None
    err = None
    if not force_harmony:
        try:
            method = try_scvi(hvg, batch_key)
        except Exception as exc:  # noqa: BLE001
            err = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()}"
            print(f"scVI failed, falling back to Harmony:\n{err}", flush=True)
    if method is None:
        method = try_harmony(
            hvg,
            batch_key,
            note="Harmony fallback (scVI failed or --force-harmony)" if err or force_harmony else "Harmony",
        )
        if err:
            hvg.uns["integrator"]["scvi_error"] = err.splitlines()[0]
    adata_full.obsm["X_int"] = np.zeros((adata_full.n_obs, hvg.obsm["X_int"].shape[1]), dtype=np.float32)
    adata_full.obsm["X_int"] = hvg.obsm["X_int"]
    adata_full.uns["integrator"] = hvg.uns["integrator"]
    adata_full.var["highly_variable"] = adata_full.var["highly_variable"]
    return adata_full, method


def neighbors_umap_leiden(adata: ad.AnnData) -> None:
    sc.pp.neighbors(adata, use_rep="X_int", n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, resolution=0.6, key_added="leiden")


def plot_umap(adata: ad.AnnData, out: Path, color: str, fname: str) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    sc.pl.umap(adata, color=color, ax=ax, show=False, frameon=False, legend_loc="on data")
    ax.set_title(color)
    fig.tight_layout()
    fig.savefig(out / fname, dpi=140)
    plt.close(fig)


def plot_scatter(samp: pd.DataFrame, block: dict, out: Path, fname: str, title: str) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 4.6))
    ok = samp.copy()
    colors = {"MPR": "#d1495b", "NMPR": "#2c6eaf", "NE": "#6b6b6b"}
    for mpr, g in ok.groupby("mpr"):
        ax.scatter(
            g["mal_CLDN4_mean"],
            g["frac_tnk"],
            c=colors.get(str(mpr), "#333"),
            s=48,
            label=str(mpr),
            edgecolors="white",
            linewidths=0.4,
        )
        for _, row in g.iterrows():
            ax.annotate(str(row["patient"]), (row["mal_CLDN4_mean"], row["frac_tnk"]), fontsize=7, alpha=0.8)
    rho, p, n = block.get("rho"), block.get("p"), block.get("n")
    rho_s = "NA" if rho is None or not np.isfinite(rho) else f"{rho:.3f}"
    p_s = "NA" if p is None or not np.isfinite(p) else f"{p:.3g}"
    ax.set_xlabel("Malignant CLDN4 mean log1p(CP10k)")
    ax.set_ylabel("T/NK fraction")
    ax.set_title(f"{title}\nn={n}  ρ={rho_s}  p={p_s}")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / fname, dpi=140)
    fig.savefig(out / fname.replace(".png", ".pdf"))
    plt.close(fig)


def fmt_rho(block: dict) -> str:
    if block["n"] < 4 or not np.isfinite(block.get("rho", np.nan)):
        return f"n={block['n']} (too few)"
    return f"ρ={block['rho']:.2f}, p={block['p']:.2g}, n={block['n']}"


def fmt_mwu(row: dict) -> str:
    if row.get("note") == "too_few_samples":
        return f"n={row['n_a']} vs {row['n_b']} (too few)"
    return (
        f"mean {row['mean_a']:.3f} vs {row['mean_b']:.3f} "
        f"(Δ={row['delta']:+.3f}); exact p={row['exact_p']:.2g}; "
        f"n={row['n_a']} vs {row['n_b']}"
    )


def nmpr_mpr(samp: pd.DataFrame, col: str) -> dict:
    post = samp.loc[samp["timing"] == "post"].copy()
    a = post.loc[post["mpr"] == "NMPR", col]
    b = post.loc[post["mpr"] == "MPR", col]
    out = exact_mwu(a, b)
    out["contrast"] = f"NMPR vs MPR {col}"
    out["dropped"] = post.loc[~np.isfinite(post[col]), "sample"].astype(str).tolist()
    return out


def write_finding(outdir: Path, summary: dict) -> None:
    method = summary["method"]
    integ = summary.get("integrator", {})
    scores = summary["scores"]
    mpr_rows = summary["nmpr_vs_mpr"]
    primary = scores["cluster_post_cldn4_tnk"]
    marker = scores["marker_post_cldn4_tnk"]
    a3 = scores["a3_post_cldn4_tnk"]
    lines = [
        "# GSE207422 scVI/Harmony latent — malignant CLDN4 (not TACSTD2 redo)",
        "",
        "Additive methods only. Public GEO UMI. **CLDN4 is the score.** TACSTD2 is a companion gene and is **never a gate**. The given A3 TACSTD2 analysis and the marker-only CLDN4 dualhigh slice are not re-argued.",
        "",
        f"**Verdict (honest n):** after a real {method} latent on GSE207422 (batch = 10x sample), Leiden-cluster malignant CLDN4 vs T/NK in the 12 post-treatment patients is {fmt_rho(primary)}. Marker-lineage malignant CLDN4 vs T/NK is {fmt_rho(marker)}. A3-malignant CLDN4 vs T/NK is {fmt_rho(a3)}. None of these is a significant anti-correlation. NMPR vs MPR on cluster-malignant CLDN4 is {fmt_mwu(mpr_rows['cluster_cldn4'])}.",
        "",
        "## What this is (and is not)",
        "",
        "| Existing slice | This slice |",
        "|---|---|",
        "| scVI combo GSE207422+GSE241934 IIT scored **TACSTD2** (n=26) | GSE207422 only; scored **CLDN4** |",
        "| Marker-only GSE207422 malignant CLDN4 (no latent) | Same public UMI, but malignant calls also from the **latent Leiden** |",
        "| A3 TACSTD2 NMPR/MPR / T/NK | Taken as given; not re-run as the primary |",
        "",
        "## Data and n",
        "",
        f"- Public GEO UMI only: **{summary['n_cells_raw']}** cells × **{summary['n_genes_raw']}** genes. After QC (≥{MIN_UMI} UMI, ≥{MIN_GENES} genes): **{summary['n_cells']}** cells. Raw GSA-Human HRA001033 was not used. Author CopyKAT barcodes are not on GEO.",
        f"- **Header n = 12 post-treatment patients** (MPR n=4 including pCR P06; NMPR n=8). The three pre-treatment biopsies (P01/P05/P08) are excluded from primary tests. All 15 samples are in the tables.",
        f"- Eligible for a Spearman: ≥{MIN_MAL} malignant-like cells and ≥{MIN_TNK} T/NK cells and a finite CLDN4 mean. Cell-level p-values are not reported.",
        f"- Marker lineage after QC: {json.dumps(summary['lineage_counts'], sort_keys=True)}.",
        f"- Leiden cluster lineage after {method}: {json.dumps(summary['cluster_lineage_counts'], sort_keys=True)}.",
        f"- Cluster-malignant empty/ineligible post samples: {primary.get('dropped') or 'none'}.",
        f"- Marker-malignant empty/ineligible post samples: {marker.get('dropped') or 'none'}.",
        f"- A3-malignant empty/ineligible post samples: {a3.get('dropped') or 'none'}.",
        "",
        "## Integration",
        "",
        f"- Method: **{method}**",
        f"- Batch key: `{integ.get('batch_key', 'sample')}` (the 15 10x libraries)",
        f"- HVGs: {summary.get('n_hvg')} (`seurat_v3`, batch-aware); latent dim {integ.get('n_latent', integ.get('n_pcs', 'NA'))}",
        f"- Cells in the latent: {integ.get('n_cells', summary['n_cells'])}",
        f"- scVI version / note: {integ.get('scvi_version', integ.get('note', ''))}",
        "",
        "CLDN4 and T/NK are scored on raw log1p(CP10k) / lineage fractions, not on the latent coordinates. The latent is used to build a joint embedding and a Leiden malignant call. That is the additive piece versus the marker-only CLDN4 slice.",
        "",
        "## Definitions",
        "",
        "| Item | Rule |",
        "|---|---|",
        "| Marker malignant-like | argmax lineage = Epithelial AND normal-lung score ≤ epithelial 75th percentile |",
        "| Cluster malignant | Leiden cluster on `X_int` whose mean marker score is Epithelial and whose mean normal-lung score is ≤ the epithelial-cluster 75th percentile |",
        "| A3-malignant | epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (same rule as the given A3 TACSTD2 slice; **not** CopyKAT). Sensitivity only. |",
        "| CLDN4 mean | sample mean `log1p(CP10k)` inside the named malignant call |",
        "| T/NK | matching lineage fraction of all QC cells |",
        "| CXCL13+ | (T/NK) AND CXCL13 UMI ≥ 1 / all QC cells |",
        "| cyto-high | T/NK AND mean log1p(GZMB, GZMA, PRF1, IFNG, NKG7) ≥ pooled T/NK 75th percentile / all QC cells |",
        "| TACSTD2 | companion only; not a gate |",
        "",
        "## Primary — post n=12, malignant CLDN4 vs immune",
        "",
        "| Malignant call | vs T/NK | vs CXCL13+ | vs cyto-high |",
        "|---|---|---|---|",
        f"| Cluster (latent Leiden) | {fmt_rho(scores['cluster_post_cldn4_tnk'])} | {fmt_rho(scores['cluster_post_cldn4_cxcl13'])} | {fmt_rho(scores['cluster_post_cldn4_cyto'])} |",
        f"| Marker lineage | {fmt_rho(scores['marker_post_cldn4_tnk'])} | {fmt_rho(scores['marker_post_cldn4_cxcl13'])} | {fmt_rho(scores['marker_post_cldn4_cyto'])} |",
        f"| A3-malignant (sensitivity) | {fmt_rho(scores['a3_post_cldn4_tnk'])} | {fmt_rho(scores['a3_post_cldn4_cxcl13'])} | {fmt_rho(scores['a3_post_cldn4_cyto'])} |",
        "",
        "## Primary — NMPR vs MPR (post only)",
        "",
        "| Score | Result |",
        "|---|---|",
        f"| Cluster-malignant CLDN4 mean | {fmt_mwu(mpr_rows['cluster_cldn4'])} |",
        f"| Marker-malignant CLDN4 mean | {fmt_mwu(mpr_rows['marker_cldn4'])} |",
        f"| A3-malignant CLDN4 mean | {fmt_mwu(mpr_rows['a3_cldn4'])} |",
        f"| Cluster T/NK fraction | {fmt_mwu(mpr_rows['cluster_tnk'])} |",
        f"| Marker T/NK fraction | {fmt_mwu(mpr_rows['marker_tnk'])} |",
        "",
        "## Sensitivity — all 15 samples (3 pre + 12 post)",
        "",
        "Pre-treatment biopsies are not a second clinical unit. They are shown so the 15-library latent is not silently subset.",
        "",
        "| Malignant call | vs T/NK | vs CXCL13+ | vs cyto-high |",
        "|---|---|---|---|",
        f"| Cluster (latent Leiden) | {fmt_rho(scores['cluster_all_cldn4_tnk'])} | {fmt_rho(scores['cluster_all_cldn4_cxcl13'])} | {fmt_rho(scores['cluster_all_cldn4_cyto'])} |",
        f"| Marker lineage | {fmt_rho(scores['marker_all_cldn4_tnk'])} | {fmt_rho(scores['marker_all_cldn4_cxcl13'])} | {fmt_rho(scores['marker_all_cldn4_cyto'])} |",
        "",
        "## Companion TACSTD2 (not a gate; not a redo)",
        "",
        "Reported so CLDN4 can be read next to TACSTD2 on the **same latent malignant call**. This does not replace the given A3 TACSTD2 write-up or the scVI combo TACSTD2 score.",
        "",
        f"- Cluster-malignant CLDN4 vs TACSTD2 (post eligible): {fmt_rho(scores['cluster_post_cldn4_vs_tacstd2'])}",
        f"- Cluster-malignant TACSTD2 vs T/NK (post eligible): {fmt_rho(scores['cluster_post_tacstd2_tnk'])}",
        f"- Marker-malignant TACSTD2 vs T/NK (post eligible): {fmt_rho(scores['marker_post_tacstd2_tnk'])}",
        "",
        "## Honest limits",
        "",
        "1. Header n=12 (4 MPR vs 8 NMPR) is the cohort. Spearman |ρ|≈0.45 has two-sided p≈0.14 at n=12. Eligible n can be lower when a malignant call is empty.",
        "2. A3-malignant-like empties some MPR residual tumors (normal-lung program). That is reported, not patched.",
        "3. Leiden-cluster malignant is a latent annotation, not Hu et al. CopyKAT. Residual unmarked epithelium can leak.",
        "4. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is not a gate.",
        "5. Expression is log1p(CP10k) from the public UMI. The latent does not change the CLDN4 numbers; it changes who is called malignant.",
        "6. This is GSE207422 only. It is not the GSE207422+GSE241934 scVI combo.",
        "",
        "## Files",
        "",
        "- `sample_table_cluster.tsv` / `sample_table_marker.tsv` / `sample_table_a3.tsv`",
        "- `cluster_annotation.tsv` / `spearman.tsv` / `nmpr_vs_mpr.tsv` / `summary.json`",
        "- Figures: cluster and marker CLDN4 vs T/NK; UMAPs (sample, marker lineage, cluster lineage)",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/scrna_scvi_cldn4/download.py",
        "python3 methods/scrna_scvi_cldn4/prepare.py",
        "python3 methods/scrna_scvi_cldn4/integrate_score.py",
        "```",
        "",
        "scVI is preferred. If training fails, the script falls back to Harmony and records the error.",
        "",
    ]
    (outdir / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, default=Path("/tmp/scrna_scvi_cldn4/h5ad/GSE207422.h5ad"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/scrna_scvi_cldn4"))
    ap.add_argument("--force-harmony", action="store_true")
    args = ap.parse_args()
    figdir = args.outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)
    n_cells_raw = int(adata.n_obs)
    n_genes_raw = int(adata.n_vars)
    if issparse(adata.X):
        adata.obs["n_umi"] = np.asarray(adata.X.sum(axis=1)).ravel()
        adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    else:
        adata.obs["n_umi"] = np.asarray(adata.X.sum(axis=1)).ravel()
        adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    keep = (adata.obs["n_umi"] >= MIN_UMI) & (adata.obs["n_genes"] >= MIN_GENES)
    adata = adata[keep].copy()
    adata.layers["counts"] = adata.X.copy()
    print(f"QC {n_cells_raw} -> {adata.n_obs}", flush=True)

    assign_marker_lineage(adata)
    adata, method = integrate(adata, batch_key="sample", force_harmony=args.force_harmony)
    neighbors_umap_leiden(adata)
    annotate_leiden_clusters(adata)

    marker = sample_table(adata, "malignant_like", "is_tnk")
    cluster = sample_table(adata, "cluster_malignant", "cluster_tnk")
    a3 = sample_table(adata, "a3_malignant", "is_tnk")
    marker.to_csv(args.outdir / "sample_table_marker.tsv", sep="\t", index=False)
    cluster.to_csv(args.outdir / "sample_table_cluster.tsv", sep="\t", index=False)
    a3.to_csv(args.outdir / "sample_table_a3.tsv", sep="\t", index=False)
    pd.DataFrame(adata.uns.get("cluster_annotation", [])).to_csv(
        args.outdir / "cluster_annotation.tsv", sep="\t", index=False
    )

    post_m = marker["timing"] == "post"
    post_c = cluster["timing"] == "post"
    post_a = a3["timing"] == "post"

    scores = {
        "cluster_post_cldn4_tnk": score_block(cluster.loc[post_c], "post cluster CLDN4 vs T/NK", "mal_CLDN4_mean", "frac_tnk"),
        "cluster_post_cldn4_cxcl13": score_block(cluster.loc[post_c], "post cluster CLDN4 vs CXCL13+", "mal_CLDN4_mean", "frac_cxcl13"),
        "cluster_post_cldn4_cyto": score_block(cluster.loc[post_c], "post cluster CLDN4 vs cyto-high", "mal_CLDN4_mean", "frac_cyto_high"),
        "marker_post_cldn4_tnk": score_block(marker.loc[post_m], "post marker CLDN4 vs T/NK", "mal_CLDN4_mean", "frac_tnk"),
        "marker_post_cldn4_cxcl13": score_block(marker.loc[post_m], "post marker CLDN4 vs CXCL13+", "mal_CLDN4_mean", "frac_cxcl13"),
        "marker_post_cldn4_cyto": score_block(marker.loc[post_m], "post marker CLDN4 vs cyto-high", "mal_CLDN4_mean", "frac_cyto_high"),
        "a3_post_cldn4_tnk": score_block(a3.loc[post_a], "post A3 CLDN4 vs T/NK", "mal_CLDN4_mean", "frac_tnk"),
        "a3_post_cldn4_cxcl13": score_block(a3.loc[post_a], "post A3 CLDN4 vs CXCL13+", "mal_CLDN4_mean", "frac_cxcl13"),
        "a3_post_cldn4_cyto": score_block(a3.loc[post_a], "post A3 CLDN4 vs cyto-high", "mal_CLDN4_mean", "frac_cyto_high"),
        "cluster_all_cldn4_tnk": score_block(cluster, "all15 cluster CLDN4 vs T/NK", "mal_CLDN4_mean", "frac_tnk"),
        "cluster_all_cldn4_cxcl13": score_block(cluster, "all15 cluster CLDN4 vs CXCL13+", "mal_CLDN4_mean", "frac_cxcl13"),
        "cluster_all_cldn4_cyto": score_block(cluster, "all15 cluster CLDN4 vs cyto-high", "mal_CLDN4_mean", "frac_cyto_high"),
        "marker_all_cldn4_tnk": score_block(marker, "all15 marker CLDN4 vs T/NK", "mal_CLDN4_mean", "frac_tnk"),
        "marker_all_cldn4_cxcl13": score_block(marker, "all15 marker CLDN4 vs CXCL13+", "mal_CLDN4_mean", "frac_cxcl13"),
        "marker_all_cldn4_cyto": score_block(marker, "all15 marker CLDN4 vs cyto-high", "mal_CLDN4_mean", "frac_cyto_high"),
        "cluster_post_cldn4_vs_tacstd2": score_block(cluster.loc[post_c], "post cluster CLDN4 vs TACSTD2", "mal_CLDN4_mean", "mal_TACSTD2_mean"),
        "cluster_post_tacstd2_tnk": score_block(cluster.loc[post_c], "post cluster TACSTD2 vs T/NK (companion)", "mal_TACSTD2_mean", "frac_tnk"),
        "marker_post_tacstd2_tnk": score_block(marker.loc[post_m], "post marker TACSTD2 vs T/NK (companion)", "mal_TACSTD2_mean", "frac_tnk"),
    }
    pd.DataFrame(list(scores.values())).to_csv(args.outdir / "spearman.tsv", sep="\t", index=False)

    mpr_rows = {
        "cluster_cldn4": nmpr_mpr(cluster, "mal_CLDN4_mean"),
        "marker_cldn4": nmpr_mpr(marker, "mal_CLDN4_mean"),
        "a3_cldn4": nmpr_mpr(a3, "mal_CLDN4_mean"),
        "cluster_tnk": nmpr_mpr(cluster, "frac_tnk"),
        "marker_tnk": nmpr_mpr(marker, "frac_tnk"),
    }
    pd.DataFrame(list(mpr_rows.values())).to_csv(args.outdir / "nmpr_vs_mpr.tsv", sep="\t", index=False)

    plot_umap(adata, figdir, "sample", "umap_sample.png")
    plot_umap(adata, figdir, "lineage", "umap_marker_lineage.png")
    plot_umap(adata, figdir, "cluster_lineage", "umap_cluster_lineage.png")
    plot_scatter(
        cluster.loc[post_c],
        scores["cluster_post_cldn4_tnk"],
        figdir,
        "cldn4_vs_tnk_cluster.png",
        f"Post cluster-malignant CLDN4 vs T/NK ({method})",
    )
    plot_scatter(
        marker.loc[post_m],
        scores["marker_post_cldn4_tnk"],
        figdir,
        "cldn4_vs_tnk_marker.png",
        "Post marker-malignant CLDN4 vs T/NK",
    )

    summary = {
        "dataset": "GSE207422",
        "task": "scVI or Harmony latent; malignant CLDN4 (not TACSTD2 redo)",
        "tacstd2_gate": False,
        "dual_high": False,
        "n_cells_raw": n_cells_raw,
        "n_genes_raw": n_genes_raw,
        "n_cells": int(adata.n_obs),
        "n_hvg": int(adata.var["highly_variable"].sum()),
        "method": method,
        "integrator": adata.uns.get("integrator", {}),
        "normal_lung_cut": float(adata.uns.get("normal_lung_cut", np.nan)),
        "cluster_nl_cut": float(adata.uns.get("cluster_nl_cut", np.nan)),
        "lineage_counts": adata.obs["lineage"].value_counts().to_dict(),
        "cluster_lineage_counts": adata.obs["cluster_lineage"].value_counts().to_dict(),
        "n_marker_malignant": int(adata.obs["malignant_like"].sum()),
        "n_cluster_malignant": int(adata.obs["cluster_malignant"].sum()),
        "n_a3_malignant": int(adata.obs["a3_malignant"].sum()),
        "scores": scores,
        "nmpr_vs_mpr": mpr_rows,
        "CLDN4_present": "CLDN4" in adata.var_names,
        "TACSTD2_present": "TACSTD2" in adata.var_names,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    write_finding(args.outdir, summary)
    print(json.dumps({"method": method, "primary": scores["cluster_post_cldn4_tnk"]}, indent=2, default=str), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
