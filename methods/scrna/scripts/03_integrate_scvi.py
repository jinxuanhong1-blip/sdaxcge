"""Integration with scVI (scvi-tools) on raw counts. GPU recommended.

Concatenate per-sample QC'd objects first (all with a `sample` and `batch` obs).
先把各样本质控后的对象合并，保证 obs 含 `sample`/`batch`。GPU 推荐。

Usage: python 03_integrate_scvi.py concat.h5ad out.h5ad [batch_key]
"""
import sys
import scanpy as sc
import scvi

def main(in_h5ad, out_h5ad, batch_key="sample"):
    scvi.settings.seed = 0
    adata = sc.read_h5ad(in_h5ad)

    # keep raw counts in .layers['counts'] for the model
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",
                                layer="counts", batch_key=batch_key, subset=True)

    scvi.model.SCVI.setup_anndata(
        adata, layer="counts", batch_key=batch_key,
        # add nuisance covariates if present:
        continuous_covariate_keys=[c for c in ["pct_counts_mt"] if c in adata.obs],
    )
    model = scvi.model.SCVI(adata, n_latent=30, n_layers=2, gene_likelihood="nb")
    model.train(max_epochs=200, early_stopping=True)

    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.layers["scvi_normalized"] = model.get_normalized_expression()

    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
    sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata)

    model.save(out_h5ad.replace(".h5ad", "_scvi_model"), overwrite=True)
    adata.write(out_h5ad)
    print(f"scVI-integrated -> {out_h5ad}  ({adata.obs['leiden'].nunique()} clusters)")

if __name__ == "__main__":
    main(*sys.argv[1:])
