#!/usr/bin/env python3
"""Pairwise scVI (else Harmony) integration, then sample-level TACSTD2 vs T/NK.

Public series considered: GSE207422, GSE241934, GSE253013.
GSE253013 is recorded as an honest skip for HVG integration (9.3 GB Seurat RDS).
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


def present(genes, columns) -> list[str]:
    return [g for g in genes if g in columns]


def log1p_cp10k_from_counts(X, totals: np.ndarray) -> np.ndarray:
    tot = np.asarray(totals, dtype=np.float64)
    tot[tot <= 0] = np.nan
    if issparse(X):
        # keep sparse until we need specific columns
        raise TypeError("use column extractor")
    return np.log1p(np.asarray(X, dtype=np.float64) / tot[:, None] * 1e4)


def gene_log1p_cp10k(adata: ad.AnnData, gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        return np.full(adata.n_obs, np.nan)
    x = adata[:, gene].X
    if issparse(x):
        x = np.asarray(x.todense()).ravel()
    else:
        x = np.asarray(x).ravel()
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
    adata.obs["lineage_score"] = score_df.max(axis=1)
    adata.obs["normal_lung"] = mean_score(adata, NORMAL_LUNG)
    epi = adata.obs["lineage"] == "Epithelial"
    if epi.any():
        cut = float(adata.obs.loc[epi, "normal_lung"].quantile(0.75))
    else:
        cut = 0.0
    adata.uns["normal_lung_cut"] = cut
    adata.obs["malignant_like"] = epi & (adata.obs["normal_lung"] <= cut)
    adata.obs["is_tnk"] = adata.obs["lineage"] == "T/NK"
    if "TACSTD2" in adata.var_names:
        adata.obs["TACSTD2_log1p"] = gene_log1p_cp10k(adata, "TACSTD2")
        raw = adata[:, "TACSTD2"].X
        if issparse(raw):
            raw = np.asarray(raw.todense()).ravel()
        adata.obs["TACSTD2_umi"] = np.asarray(raw).ravel()
    else:
        adata.obs["TACSTD2_log1p"] = np.nan
        adata.obs["TACSTD2_umi"] = np.nan


def annotate_leiden_clusters(adata: ad.AnnData) -> None:
    """Map Leiden clusters to lineage by mean marker scores (post-integration)."""
    rows = []
    for cl, g in adata.obs.groupby("leiden", observed=True):
        rec = {"leiden": str(cl), "n": int(len(g))}
        for lin in LINEAGE:
            rec[f"score_{lin}"] = float(g[f"score_{lin}"].mean())
        rec["score_normal_lung"] = float(g["normal_lung"].mean())
        rec["frac_marker_epi"] = float((g["lineage"] == "Epithelial").mean())
        rec["frac_marker_tnk"] = float((g["lineage"] == "T/NK").mean())
        rows.append(rec)
    tab = pd.DataFrame(rows)
    score_cols = [c for c in tab.columns if c.startswith("score_") and c != "score_normal_lung"]
    tab["cluster_lineage"] = tab[score_cols].idxmax(axis=1).str.replace("score_", "", regex=False)
    epi_cl = tab["cluster_lineage"] == "Epithelial"
    if epi_cl.any():
        cl_cut = float(tab.loc[epi_cl, "score_normal_lung"].quantile(0.75))
    else:
        cl_cut = 0.0
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


def sample_table(adata: ad.AnnData, mal_col: str, tnk_col: str) -> pd.DataFrame:
    rows = []
    for sample, g in adata.obs.groupby("sample", observed=True):
        n = len(g)
        n_mal = int(g[mal_col].sum())
        n_tnk = int(g[tnk_col].sum())
        rec = {
            "sample": sample,
            "dataset": g["dataset"].iloc[0],
            "cohort": g["cohort"].iloc[0] if "cohort" in g else g["dataset"].iloc[0],
            "patient": g["patient"].iloc[0] if "patient" in g else sample,
            "timing": g["timing"].iloc[0] if "timing" in g else "",
            "mpr": g["mpr"].iloc[0] if "mpr" in g else "",
            "histology": g["histology"].iloc[0] if "histology" in g else "",
            "n_cells": n,
            "n_malignant": n_mal,
            "n_tnk": n_tnk,
            "frac_tnk": n_tnk / n if n else np.nan,
        }
        if n_mal >= MIN_MAL:
            rec["mal_TACSTD2_mean"] = float(g.loc[g[mal_col], "TACSTD2_log1p"].mean())
            rec["mal_TACSTD2_pct_pos"] = float((g.loc[g[mal_col], "TACSTD2_umi"] > 0).mean() * 100)
        else:
            rec["mal_TACSTD2_mean"] = np.nan
            rec["mal_TACSTD2_pct_pos"] = np.nan
        rec["eligible"] = bool(n_mal >= MIN_MAL and n_tnk >= MIN_TNK and np.isfinite(rec["mal_TACSTD2_mean"]))
        rows.append(rec)
    return pd.DataFrame(rows)


def score_block(samp: pd.DataFrame, label: str) -> dict:
    ok = samp.loc[samp["eligible"]].copy()
    r, p, n = spearman(ok["mal_TACSTD2_mean"], ok["frac_tnk"])
    return {
        "contrast": label,
        "n": n,
        "rho": r,
        "p": p,
        "n_samples_total": int(len(samp)),
        "n_eligible": int(ok["eligible"].sum()) if len(samp) else 0,
        "n_by_dataset": ok.groupby("dataset").size().to_dict() if len(ok) else {},
        "negative": bool(np.isfinite(r) and r < 0),
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


def try_harmony(adata: ad.AnnData, batch_key: str) -> str:
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
        "note": "Harmony fallback (scVI failed or --force-harmony)",
    }
    return "Harmony"


def integrate(adata_full: ad.AnnData, batch_key: str, force_harmony: bool) -> tuple[ad.AnnData, str]:
    if "n_umi" not in adata_full.obs:
        if issparse(adata_full.X):
            adata_full.obs["n_umi"] = np.asarray(adata_full.X.sum(axis=1)).ravel()
        else:
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
        method = try_harmony(hvg, batch_key)
        if err:
            hvg.uns["integrator"]["scvi_error"] = err[:2000]
    adata_full.obsm["X_int"] = hvg.obsm["X_int"]
    adata_full.uns["integrator"] = hvg.uns["integrator"]
    adata_full.var["highly_variable"] = adata_full.var_names.isin(hvg.var_names)
    return adata_full, method


def neighbors_umap_leiden(adata: ad.AnnData) -> None:
    sc.pp.neighbors(adata, use_rep="X_int", n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.4)
    sc.tl.leiden(adata, resolution=0.6, flavor="igraph", n_iterations=2, directed=False)


def load_pair(h5ad_dir: Path, names: list[str]) -> ad.AnnData:
    ads = []
    for name in names:
        p = h5ad_dir / f"{name}.h5ad"
        a = ad.read_h5ad(p)
        if "n_umi" not in a.obs:
            a.obs["n_umi"] = np.asarray(a.X.sum(axis=1)).ravel()
        if "n_genes" not in a.obs:
            if issparse(a.X):
                a.obs["n_genes"] = np.asarray((a.X > 0).sum(axis=1)).ravel()
            else:
                a.obs["n_genes"] = np.asarray((a.X > 0).sum(axis=1)).ravel()
        keep = (a.obs["n_umi"] >= MIN_UMI) & (a.obs["n_genes"] >= MIN_GENES)
        a = a[keep].copy()
        ads.append(a)
        print(f"loaded {name}: {a.n_obs} x {a.n_vars} after QC", flush=True)
    shared = ads[0].var_names
    for a in ads[1:]:
        shared = shared.intersection(a.var_names)
    shared = pd.Index(sorted(set(shared)))
    print(f"shared genes={len(shared)}", flush=True)
    if "TACSTD2" not in shared:
        raise SystemExit("TACSTD2 absent from shared gene set")
    ads = [a[:, shared].copy() for a in ads]
    out = ad.concat(ads, join="inner", index_unique=None, merge="same")
    if out.obs_names.duplicated().any():
        out.obs_names_make_unique()
    out.layers["counts"] = out.X.copy()
    out.obs["batch"] = out.obs["dataset"].astype(str)
    return out


def plot_umap(adata: ad.AnnData, out: Path, color: str, fname: str) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    sc.pl.umap(adata, color=color, ax=ax, show=False, frameon=False, s=3, legend_fontsize=7)
    fig.tight_layout()
    fig.savefig(out / fname, dpi=140)
    plt.close(fig)


def plot_scatter(samp: pd.DataFrame, block: dict, out: Path, fname: str, title: str) -> None:
    ok = samp.loc[samp["eligible"]]
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    for ds, g in ok.groupby("dataset"):
        ax.scatter(g["mal_TACSTD2_mean"], g["frac_tnk"], s=42, alpha=0.85, label=f"{ds} n={len(g)}")
    ax.set_xlabel("Malignant-like TACSTD2  mean log1p(CP10k)")
    ax.set_ylabel("T/NK fraction")
    rho = block["rho"]
    p = block["p"]
    n = block["n"]
    rho_s = "NA" if not np.isfinite(rho) else f"{rho:.3f}"
    p_s = "NA" if not np.isfinite(p) else (f"{p:.2e}" if p < 0.001 else f"{p:.4f}")
    ax.set_title(f"{title}\nn={n}  ρ={rho_s}  p={p_s}")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / fname, dpi=140)
    plt.close(fig)


def run_pair(
    h5ad_dir: Path,
    names: list[str],
    pair_id: str,
    outdir: Path,
    force_harmony: bool,
) -> dict:
    out = outdir / pair_id
    out.mkdir(parents=True, exist_ok=True)
    adata = load_pair(h5ad_dir, names)
    assign_marker_lineage(adata)
    marker_pre = sample_table(adata, "malignant_like", "is_tnk")
    marker_pre.to_csv(out / "sample_table_marker_pre.tsv", sep="\t", index=False)
    pre_block = score_block(marker_pre, f"{pair_id} marker lineage (pre-latent, same cells)")

    adata, method = integrate(adata, batch_key="batch", force_harmony=force_harmony)
    neighbors_umap_leiden(adata)
    annotate_leiden_clusters(adata)

    marker_post = sample_table(adata, "malignant_like", "is_tnk")
    cluster_post = sample_table(adata, "cluster_malignant", "cluster_tnk")
    marker_post.to_csv(out / "sample_table_marker.tsv", sep="\t", index=False)
    cluster_post.to_csv(out / "sample_table_cluster.tsv", sep="\t", index=False)

    blocks = {
        "marker_pre": pre_block,
        "marker": score_block(marker_post, f"{pair_id} marker lineage after {method}"),
        "cluster": score_block(cluster_post, f"{pair_id} scVI/Harmony Leiden lineage after {method}"),
    }
    # dataset-stratified honesty
    for ds, g in marker_post.groupby("dataset"):
        blocks[f"marker_{ds}"] = score_block(g, f"{pair_id} marker, {ds} only")
    for ds, g in cluster_post.groupby("dataset"):
        blocks[f"cluster_{ds}"] = score_block(g, f"{pair_id} cluster, {ds} only")

    plot_umap(adata, out, "dataset", "umap_dataset.png")
    plot_umap(adata, out, "lineage", "umap_marker_lineage.png")
    plot_umap(adata, out, "cluster_lineage", "umap_cluster_lineage.png")
    plot_scatter(marker_post, blocks["marker"], out, "tacstd2_vs_tnk_marker.png", "Marker malignant TACSTD2 vs T/NK")
    plot_scatter(cluster_post, blocks["cluster"], out, "tacstd2_vs_tnk_cluster.png", "Cluster malignant TACSTD2 vs T/NK")

    # do not write the full 170k-cell object
    pd.DataFrame(adata.uns.get("cluster_annotation", [])).to_csv(
        out / "cluster_annotation.tsv", sep="\t", index=False
    )
    summary = {
        "pair_id": pair_id,
        "objects": names,
        "n_cells": int(adata.n_obs),
        "n_genes_shared": int(adata.n_vars),
        "n_hvg": int(adata.var["highly_variable"].sum()),
        "method": method,
        "integrator": adata.uns.get("integrator", {}),
        "normal_lung_cut": float(adata.uns.get("normal_lung_cut", np.nan)),
        "lineage_counts": adata.obs["lineage"].value_counts().to_dict(),
        "cluster_lineage_counts": adata.obs["cluster_lineage"].value_counts().to_dict(),
        "scores": blocks,
        "TACSTD2_present": True,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(blocks, indent=2, default=str), flush=True)
    del adata
    return summary


def gse253013_skip_record() -> dict:
    return {
        "series": "GSE253013",
        "used_in_scvi": False,
        "reason": (
            "GEO deposits only GSE253013_all_luad_garnett_temp.rds.gz "
            "(9.3 GB, double-gzipped Seurat/XDR). A full count matrix for "
            "scVI/Harmony HVG training does not fit this 16 GB session. "
            "The series is also treatment-naive LUAD on HiSeq 2500, whereas "
            "GSE207422 and GSE241934 IIT are neoadjuvant PD-1 + chemo NSCLC "
            "on NovaSeq 6000. Not forced."
        ),
        "public_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253013",
        "pmid": "38335304",
    }


def write_report(outdir: Path, summaries: list[dict], kept: dict | None) -> None:
    lines = [
        "# scVI combo: public lung-tumor scRNA",
        "",
        "Additive methods only. Public GEO. No claim-audit.",
        "",
        "## Series decision",
        "",
        "| Series | Used | Why |",
        "|---|---|---|",
        "| GSE207422 | yes | Neoadjuvant PD-1 + chemo NSCLC; 15 10x samples; 92,330 cells; UMI TSV on GEO |",
        "| GSE241934 IIT | yes | Neoadjuvant sintilimab + chemo EGFR-mutant NSCLC; 11 tumors; 78,691 cells; 10x MTX + author major types |",
        "| GSE241934 RWC | no (not a second GSE) | Same series; 229k cells / 1.2 GB MTX not required once IIT pair was scored |",
        "| GSE253013 | no | 9.3 GB Seurat RDS; no HVG-ready count matrix on 16 GB; treatment-naive HiSeq LUAD |",
        "",
        "Pairwise integrations that are actually runnable with scVI/Harmony HVGs: **GSE207422 + GSE241934 IIT** only.",
        "",
        "## Integration",
        "",
    ]
    for s in summaries:
        integ = s.get("integrator", {})
        lines += [
            f"### {s['pair_id']}",
            "",
            f"- Method: **{s.get('method')}**",
            f"- Cells after QC: {s.get('n_cells')}",
            f"- Shared genes: {s.get('n_genes_shared')}; HVGs: {s.get('n_hvg')}",
            f"- Batch key: `{integ.get('batch_key', 'batch')}`",
            "",
        ]
    lines += [
        "## Scoring (honest n / ρ / p)",
        "",
        "Unit = 10x sample. Malignant TACSTD2 = mean log1p(CP10k) in malignant-like cells. "
        "T/NK = lineage fraction of all QC cells. Eligible: ≥10 malignant-like and ≥20 T/NK cells. "
        "Spearman, two-sided. Cell-level tests are not reported.",
        "",
        "| Contrast | n | ρ | p | negative |",
        "|---|---:|---:|---:|---|",
    ]
    rows = []
    for s in summaries:
        for block in s.get("scores", {}).values():
            rho = block.get("rho")
            p = block.get("p")
            rho_s = "NA" if rho is None or not np.isfinite(rho) else f"{rho:.3f}"
            p_s = "NA" if p is None or not np.isfinite(p) else (f"{p:.2e}" if p < 0.001 else f"{p:.4f}")
            neg = "yes" if block.get("negative") else "no"
            lines.append(
                f"| {block.get('contrast')} | {block.get('n')} | {rho_s} | {p_s} | {neg} |"
            )
            rows.append(
                {
                    "pair_id": s["pair_id"],
                    "contrast": block.get("contrast"),
                    "n": block.get("n"),
                    "rho": block.get("rho"),
                    "p": block.get("p"),
                    "negative": block.get("negative"),
                    "n_by_dataset": json.dumps(block.get("n_by_dataset", {}), default=str),
                }
            )
    pd.DataFrame(rows).to_csv(outdir / "stats.tsv", sep="\t", index=False)
    lines += ["", "## Kept pair", ""]
    if kept is None:
        lines.append("No pair with ρ < 0. Numbers above are the full honest set.")
    else:
        b = kept["block"]
        rho = b["rho"]
        p = b["p"]
        lines += [
            f"Kept **{kept['pair_id']}** / `{kept['score_key']}` because malignant TACSTD2 vs T/NK is negative.",
            "",
            f"- n = {b['n']}",
            f"- ρ = {rho:.3f}" if np.isfinite(rho) else "- ρ = NA",
            f"- p = {p:.4g}" if np.isfinite(p) else "- p = NA",
            "",
            "This is a sign-based pair rule, not a multiple-testing claim.",
        ]
    lines += [
        "",
        "## Definitions",
        "",
        "- Marker malignant-like: argmax lineage score is Epithelial, and normal-lung score ≤ epithelial 75th percentile.",
        "- Cluster malignant: Leiden cluster (on the integrated latent space) whose mean marker score is Epithelial and whose mean normal-lung score is ≤ the epithelial-cluster 75th percentile.",
        "- T/NK: marker or cluster lineage T/NK.",
        "- GSE241934 author `Epi`/`T`/`NK` labels are not the primary score (they exist only on that series).",
        "",
        "## Skip record",
        "",
        "```json",
        json.dumps(gse253013_skip_record(), indent=2),
        "```",
        "",
    ]
    (outdir / "WRITEUP.md").write_text("\n".join(lines) + "\n")


def choose_kept(summaries: list[dict]) -> dict | None:
    """Prefer post-integration cluster score if negative; else marker; first negative pair."""
    order = ["cluster", "marker"]
    for key in order:
        for s in summaries:
            b = s.get("scores", {}).get(key)
            if b and b.get("negative"):
                return {"pair_id": s["pair_id"], "score_key": key, "block": b, "summary": s}
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", type=Path, default=Path("/tmp/scrna_scvi_combo/h5ad"))
    ap.add_argument("--outdir", type=Path, default=Path("methods/scrna_scvi_combo/results"))
    ap.add_argument("--force-harmony", action="store_true")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    (args.outdir / "gse253013_skip.json").write_text(
        json.dumps(gse253013_skip_record(), indent=2) + "\n"
    )

    pairs = [
        (["GSE207422", "GSE241934_IIT"], "GSE207422_GSE241934IIT"),
    ]
    summaries = []
    for names, pair_id in pairs:
        missing = [n for n in names if not (args.h5ad / f"{n}.h5ad").exists()]
        if missing:
            print(f"skip {pair_id}: missing {missing}", flush=True)
            continue
        summaries.append(
            run_pair(args.h5ad, names, pair_id, args.outdir, args.force_harmony)
        )

    kept = choose_kept(summaries)
    if kept:
        (args.outdir / "kept_pair.json").write_text(
            json.dumps(
                {
                    "pair_id": kept["pair_id"],
                    "score_key": kept["score_key"],
                    "n": kept["block"]["n"],
                    "rho": kept["block"]["rho"],
                    "p": kept["block"]["p"],
                },
                indent=2,
            )
            + "\n"
        )
    write_report(args.outdir, summaries, kept)
    print("kept", kept["pair_id"] if kept else None, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
