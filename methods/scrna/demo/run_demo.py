#!/usr/bin/env python3
"""Minimal reproducible demo for the lung-cancer ICI scRNA playbook.

Dataset: GSE207422 (Hu et al., Genome Medicine 2023) NSCLC neoadjuvant
PD-1 + chemo scRNA-seq. The GEO supplementary matrix is the AUTHORS'
POST-QC UMI count matrix (92,330 cells x 24,292 genes), so ambient-RNA
removal and doublet filtering were already applied upstream. This script
demonstrates the *downstream* half of the playbook: normalization, HVG,
PCA, Harmony integration by sample, Leiden clustering, marker-based
lineage annotation, epithelial module scoring (TACSTD2/CLDN4/junction),
immune neighborhood scores (CD8/TLS/CXCL13), and pseudobulk aggregation
joined to MPR/RECIST labels.

No statistics are fabricated: every number in the figures/tables is
computed from the downloaded data. Small-n response comparisons are
labeled descriptive/exploratory. Pseudobulk DE is run by
`pseudobulk_de.py` (pydeseq2) on the exported matrix.

Default: stratified subsample (max 1,200 cells/sample) so the demo
finishes on a 16 GB CPU box. Pass `--max-per-sample 0` for all cells.

Run:
    python run_demo.py \\
        --matrix /tmp/geo/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \\
        --meta   /tmp/geo/GSE207422_NSCLC_scRNAseq_metadata.xlsx \\
        --outdir .
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_counts_sparse(path: str):
    """Stream a gene-by-cell TSV into CSR (cells x genes), keeping only nonzeros."""
    opener = gzip.open if path.endswith(".gz") else open
    log(f"Streaming sparse load from {path}")
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = header[1:]
        n_cells = len(cells)
        log(f"{n_cells} cells in header")
        genes: list[str] = []
        data_chunks: list[np.ndarray] = []
        index_chunks: list[np.ndarray] = []
        indptr = [0]
        for i, line in enumerate(fh, 1):
            tab = line.find("\t")
            genes.append(line[:tab])
            vals = np.fromstring(line[tab + 1 :], sep="\t", dtype=np.float32)
            if vals.size != n_cells:
                raise ValueError(
                    f"gene {genes[-1]}: expected {n_cells} values, got {vals.size}"
                )
            nz = np.flatnonzero(vals)
            data_chunks.append(vals[nz])
            index_chunks.append(nz.astype(np.int32, copy=False))
            indptr.append(indptr[-1] + nz.size)
            if i % 4000 == 0:
                log(f"  ...{i} genes, {indptr[-1]:,} nnz")
    data = np.concatenate(data_chunks) if data_chunks else np.array([], np.float32)
    indices = (
        np.concatenate(index_chunks) if index_chunks else np.array([], np.int32)
    )
    mat = sp.csr_matrix(
        (data, indices, np.asarray(indptr, dtype=np.int64)),
        shape=(len(genes), n_cells),
    )
    log(f"Loaded {len(genes)} genes x {n_cells} cells, nnz={mat.nnz:,}")
    return mat.T.tocsr(), np.array(cells), np.array(genes)


def stratified_subsample(adata, max_per_sample: int, rng: np.random.Generator):
    if max_per_sample <= 0:
        return adata
    keep = []
    for s, idx in adata.obs.groupby("sample", observed=True).groups.items():
        idx = np.asarray(idx)
        if len(idx) > max_per_sample:
            idx = rng.choice(idx, size=max_per_sample, replace=False)
        keep.append(idx)
    keep = np.concatenate(keep)
    log(f"Stratified subsample: {adata.n_obs} -> {len(keep)} cells "
        f"(max {max_per_sample}/sample)")
    return adata[keep].copy()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--max-per-sample", type=int, default=1200,
                    help="cap cells per sample (0 = keep all)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    figdir = os.path.join(args.outdir, "figures")
    os.makedirs(figdir, exist_ok=True)

    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    sc.settings.verbosity = 1
    sc.settings.figdir = figdir
    rng = np.random.default_rng(args.seed)

    X, cells, genes = load_counts_sparse(args.matrix)
    adata = ad.AnnData(
        X=X, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes)
    )
    adata.var_names_make_unique()
    log(f"AnnData: {adata.n_obs} cells x {adata.n_vars} genes")

    adata.obs["sample"] = [c.rsplit("_", 1)[0] for c in adata.obs_names]

    meta = pd.read_excel(args.meta)
    meta = meta.dropna(subset=["Sample"]).copy()
    meta = meta[meta["Sample"].astype(str).str.startswith("BD")].set_index("Sample")
    for col in ["Pathologic Response", "RECIST", "Pathology", "PD1 Antibody",
                "Resource", "Residual Tumor", "Patient"]:
        adata.obs[col] = adata.obs["sample"].map(meta[col]).astype("category")
    log("Samples matched to metadata: "
        f"{adata.obs['sample'].isin(meta.index).mean():.2%}")

    n_full = int(adata.n_obs)
    adata = stratified_subsample(adata, args.max_per_sample, rng)
    adata.layers["counts"] = adata.X.copy()

    # ---- per-cell QC (already author-filtered; we still document it) ----
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False
    )

    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    axes[0].hist(adata.obs["n_genes_by_counts"], bins=60, color="#4472c4")
    axes[0].set_title("Genes per cell")
    axes[0].set_xlabel("n_genes")
    axes[1].hist(np.log10(adata.obs["total_counts"] + 1), bins=60, color="#70ad47")
    axes[1].set_title("log10 total UMI")
    axes[1].set_xlabel("log10 UMI")
    axes[2].hist(adata.obs["pct_counts_mt"], bins=60, color="#c00000")
    axes[2].set_title("% mito")
    axes[2].set_xlabel("pct_counts_mt")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "01_qc_distributions.png"), dpi=130)
    plt.close(fig)

    # ---- normalize + HVG + PCA ----
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2000, flavor="seurat", batch_key="sample"
    )
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=30, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    # ---- Harmony over sample ----
    try:
        sc.external.pp.harmony_integrate(adata, key="sample")
        rep = "X_pca_harmony"
        log("Harmony integration done")
    except Exception as exc:  # noqa: BLE001
        log(f"Harmony unavailable ({exc}); falling back to raw PCA")
        rep = "X_pca"

    sc.pp.neighbors(adata, n_neighbors=15, use_rep=rep)
    sc.tl.leiden(
        adata, resolution=1.0, key_added="leiden",
        flavor="igraph", n_iterations=2, directed=False,
    )
    sc.tl.umap(adata)
    log(f"{adata.obs['leiden'].nunique()} Leiden clusters")

    # ---- marker-based lineage annotation ----
    lineage_markers = {
        "T/NK":        ["CD3D", "CD3E", "CD2", "TRAC", "NKG7", "GNLY"],
        "B/Plasma":    ["MS4A1", "CD79A", "CD79B", "IGHG1", "MZB1"],
        "Myeloid":     ["LYZ", "CD68", "CD14", "FCGR3A", "C1QA", "C1QB"],
        "Mast":        ["TPSAB1", "TPSB2", "CPA3"],
        "Epithelial":  ["EPCAM", "KRT19", "KRT18", "CDH1", "SFTPC"],
        "Endothelial": ["PECAM1", "VWF", "CLDN5"],
        "Fibroblast":  ["COL1A1", "COL1A2", "DCN", "LUM"],
    }
    for name, gs in lineage_markers.items():
        gs = [g for g in gs if g in adata.raw.var_names]
        sc.tl.score_genes(adata, gs, score_name=f"score_{name}", use_raw=True)

    score_cols = [f"score_{k}" for k in lineage_markers]
    per_cluster = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    cluster_label = per_cluster.idxmax(axis=1).str.replace("score_", "", regex=False)
    adata.obs["lineage"] = adata.obs["leiden"].map(cluster_label).astype("category")
    log("Lineage composition:\n" + adata.obs["lineage"].value_counts().to_string())

    # ---- epithelial module + immune neighborhood scores ----
    modules = {
        "TACSTD2_CLDN4_junction": [
            "TACSTD2", "CLDN4", "CLDN3", "CLDN7",
            "CDH1", "TJP1", "OCLN", "F11R", "CLDN18",
        ],
        "CD8_cytotoxic": ["CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "IFNG", "NKG7"],
        "TLS": ["CXCL13", "CCL19", "CCL21", "CR2", "CXCR5", "MS4A1", "SELL"],
        "CXCL13": ["CXCL13"],
    }
    for name, gs in modules.items():
        present = [g for g in gs if g in adata.raw.var_names]
        if present:
            sc.tl.score_genes(adata, present, score_name=name, use_raw=True)
        log(f"module {name}: {len(present)}/{len(gs)} genes present")

    sc.pl.umap(adata, color=["lineage"], show=False, save="_02_lineage.png",
               legend_loc="on data", frameon=False)
    sc.pl.umap(adata, color=["sample"], show=False, save="_03_sample.png",
               frameon=False)
    for m in ["TACSTD2_CLDN4_junction", "CD8_cytotoxic", "TLS", "CXCL13"]:
        if m in adata.obs:
            sc.pl.umap(adata, color=[m], show=False, save=f"_04_{m}.png",
                       frameon=False, cmap="magma")

    # epithelial module vs response (per-sample means, epithelial cells only)
    epi = adata[adata.obs["lineage"] == "Epithelial"]
    if epi.n_obs > 0:
        agg = (
            epi.obs.groupby("sample", observed=True)[["TACSTD2_CLDN4_junction"]]
            .mean()
            .join(meta[["Pathologic Response", "RECIST"]])
        )
        agg.to_csv(os.path.join(args.outdir, "epithelial_module_by_sample.csv"))
        resp_map = {"pCR": "MPR", "MPR": "MPR", "NMPR": "NMPR", "NE": "NE"}
        agg["resp_group"] = agg["Pathologic Response"].map(
            lambda x: resp_map.get(str(x), str(x))
        )
        fig, ax = plt.subplots(figsize=(5, 4.5))
        groups = [g for g in ["MPR", "NMPR"] if g in set(agg["resp_group"])]
        data = [
            agg.loc[agg["resp_group"] == g, "TACSTD2_CLDN4_junction"].values
            for g in groups
        ]
        ax.boxplot(data, tick_labels=groups)
        for i, d in enumerate(data, 1):
            ax.scatter(
                rng.normal(i, 0.05, len(d)), d, color="#c00000", zorder=3, s=28
            )
        ax.set_ylabel("Epithelial TACSTD2/CLDN4/junction score\n(per-sample mean)")
        ax.set_title(
            "Epithelial module by pathologic response\n"
            "(GSE207422; descriptive, n small)"
        )
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, "05_epi_module_by_response.png"), dpi=130)
        plt.close(fig)

    # ---- pseudobulk export (sample x gene raw counts) ----
    counts = adata.layers["counts"]
    samples = adata.obs["sample"].values
    uniq = pd.unique(samples)
    pb = np.zeros((len(uniq), adata.n_vars), dtype=np.float64)
    for i, s in enumerate(uniq):
        pb[i] = np.asarray(counts[samples == s].sum(axis=0)).ravel()
    pd.DataFrame(pb, index=uniq, columns=adata.var_names).astype(int).to_csv(
        os.path.join(args.outdir, "pseudobulk_allcells_counts.csv")
    )

    epi_mask = (adata.obs["lineage"] == "Epithelial").values
    if epi_mask.sum() > 0:
        pb_epi = np.zeros((len(uniq), adata.n_vars), dtype=np.float64)
        for i, s in enumerate(uniq):
            m = epi_mask & (samples == s)
            if m.sum():
                pb_epi[i] = np.asarray(counts[m].sum(axis=0)).ravel()
        pd.DataFrame(pb_epi, index=uniq, columns=adata.var_names).astype(int).to_csv(
            os.path.join(args.outdir, "pseudobulk_epithelial_counts.csv")
        )

    design = meta.loc[
        [s for s in uniq if s in meta.index],
        ["Patient", "Resource", "Pathology", "PD1 Antibody",
         "Pathologic Response", "Residual Tumor", "RECIST"],
    ].copy()
    design["n_cells_total"] = [int((samples == s).sum()) for s in design.index]
    design["n_cells_epithelial"] = [
        int((epi_mask & (samples == s)).sum()) for s in design.index
    ]
    design.to_csv(os.path.join(args.outdir, "pseudobulk_design.csv"))

    pd.crosstab(adata.obs["sample"], adata.obs["lineage"]).to_csv(
        os.path.join(args.outdir, "composition_by_sample.csv")
    )

    summary = {
        "n_cells_in_matrix": n_full,
        "n_cells_analyzed": int(adata.n_obs),
        "max_per_sample": args.max_per_sample,
        "n_genes": int(adata.n_vars),
        "n_samples": int(len(uniq)),
        "n_leiden_clusters": int(adata.obs["leiden"].nunique()),
        "integration_rep": rep,
        "lineage_counts": adata.obs["lineage"].value_counts().to_dict(),
        "modules_present": {
            k: [g for g in v if g in adata.raw.var_names] for k, v in modules.items()
        },
        "note": (
            "Matrix is authors' post-QC counts; ambient/doublet removal was "
            "done upstream. All values computed from data; DE in pseudobulk_de.py. "
            "Response comparisons are exploratory (small n)."
        ),
    }
    with open(os.path.join(args.outdir, "demo_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    log("Wrote demo_summary.json")
    # do not write a large .h5ad into git
    log("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
