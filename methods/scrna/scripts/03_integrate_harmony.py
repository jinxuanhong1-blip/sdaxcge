"""Integration with Harmony (fast, CPU-friendly). Default when scVI compute is limited.
Harmony 在 PCA 空间快速校正批次；算力有限时的稳妥默认（本仓库 demo 亦用此法）。

Usage: python 03_integrate_harmony.py concat.h5ad out.h5ad [batch_key]
"""
import sys
import scanpy as sc

def main(in_h5ad, out_h5ad, batch_key="sample"):
    adata = sc.read_h5ad(in_h5ad)
    if "log1p" not in adata.uns_keys():
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=batch_key)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=30, svd_solver="arpack")
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]

    # Call harmonypy directly. scanpy 1.12 sc.external.pp.harmony_integrate can
    # fail to write X_pca_harmony (non-contiguous Z_corr.T fails the obsm check).
    import numpy as np
    import harmonypy as hm
    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, batch_key)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[1] == adata.n_obs:
        Z = Z.T
    adata.obsm["X_pca_harmony"] = np.ascontiguousarray(Z)
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata)
    adata.write(out_h5ad)
    print(f"Harmony-integrated -> {out_h5ad}  ({adata.obs['leiden'].nunique()} clusters)")

if __name__ == "__main__":
    main(*sys.argv[1:])
