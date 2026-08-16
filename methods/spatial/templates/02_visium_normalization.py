"""02 · Visium normalization + embedding. Visium 归一化 + 降维聚类。

Log-normalize (or Pearson residuals), select HVGs, PCA/neighbors/UMAP, Leiden
clustering, and per-spot marker-set scores (tumor TACSTD2/CLDN4, T, B, TLS).
对数归一化（或 Pearson 残差），选取高变基因，PCA/邻接/UMAP，Leiden 聚类，并对
每个 spot 计算标记基因集评分（肿瘤 TACSTD2/CLDN4、T、B、TLS）。

Tools: scanpy >= 1.10.  No results are fabricated.

Usage:
    python 02_visium_normalization.py --input <sample_qc.h5ad> --sample <id>
"""
from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc

from _utils import load_config, resolve, ensure_dir, genes_in


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="output of 01_visium_qc.py")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    n = cfg["normalization"]["visium"]
    out = ensure_dir(resolve(cfg["paths"]["output_dir"]) / "visium" / args.sample)

    adata = sc.read_h5ad(args.input)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    if n["method"] == "pearson_residuals":
        # Analytic Pearson residuals (recommended for count-based clustering).
        sc.experimental.pp.normalize_pearson_residuals(adata)
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n["n_hvg"], flavor="pearson_residuals",
            layer="counts",
        )
    else:  # lognorm (default)
        sc.pp.normalize_total(adata, target_sum=n["target_sum"])
        sc.pp.log1p(adata)
        adata.layers["lognorm"] = adata.X.copy()
        sc.pp.highly_variable_genes(adata, n_top_genes=n["n_hvg"])

    sc.pp.pca(adata, n_comps=30, use_highly_variable=True)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden")

    # Marker-set scores. 标记基因集评分。
    # These feed the niche template's spot labelling (tumor vs T/B/TLS).
    for label, genes in cfg["markers"].items():
        present = genes_in(adata.var_names, genes)
        if present:
            sc.tl.score_genes(adata, present, score_name=f"score_{label}")

    adata.write(out / f"{args.sample}_norm.h5ad")
    print(
        f"[02_norm] {args.sample}: {adata.n_obs} spots, method={n['method']}, "
        f"{adata.obs['leiden'].nunique()} leiden clusters "
        f"-> {out / (args.sample + '_norm.h5ad')}"
    )


if __name__ == "__main__":
    main()
