"""QC, normalize, integrate (Harmony), cluster (Leiden), and annotate compartments.

Also assigns response groups and writes an annotated .h5ad plus QC/marker figures.
"""
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 1
sc.settings.figdir = "results/align_gse207422"
RAW = "data/gse207422/gse207422_raw.h5ad"
OUT = "data/gse207422/gse207422_annotated.h5ad"

adata = sc.read_h5ad(RAW)
print("loaded", adata.shape)

# ---- response group / treatment status ----
resp = adata.obs["Pathologic Response"].astype(str)
group = np.where(resp.isin(["MPR", "pCR"]), "MPR",
         np.where(resp == "NMPR", "NMPR", "other"))
adata.obs["response_group"] = pd.Categorical(group)
adata.obs["treatment"] = np.where(
    adata.obs["Resource"].str.contains("Pre-treatment"), "Pre", "Post")
print(adata.obs.groupby(["Sample", "treatment", "response_group"], observed=True).size())

# ---- QC ----
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "MT."))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
print("median genes/cell:", np.median(adata.obs["n_genes_by_counts"]))
print("median pct_mt:", np.median(adata.obs["pct_counts_mt"]))

n0 = adata.n_obs
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
sc.pp.filter_genes(adata, min_cells=3)
print(f"QC: {n0} -> {adata.n_obs} cells, {adata.n_vars} genes")

# ---- normalize ----
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# ---- HVG / PCA / Harmony / clustering ----
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", batch_key="Sample")
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50)
adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
import harmonypy
ho = harmonypy.run_harmony(adata_hvg.obsm["X_pca"], adata.obs, ["Sample"], max_iter_harmony=20)
Z = np.asarray(ho.Z_corr)
if Z.shape[0] != adata.n_obs:
    Z = Z.T
assert Z.shape[0] == adata.n_obs, Z.shape
adata.obsm["X_pca_harmony"] = Z
del adata_hvg

sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=50)
sc.tl.leiden(adata, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata)
print("leiden clusters:", adata.obs["leiden"].nunique())

# ---- compartment marker scoring ----
markers = {
    "Epithelial": ["EPCAM", "KRT19", "KRT18", "KRT8", "KRT7", "CDH1", "SFN", "ELF3", "KRT17"],
    "T": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "IL7R"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "KLRB1"],
    "B_Plasma": ["CD79A", "MS4A1", "CD19", "MZB1", "IGHG1", "JCHAIN"],
    "Myeloid": ["LYZ", "CD68", "CD14", "FCGR3A", "C1QA", "AIF1", "CD163"],
    "Mast": ["TPSAB1", "CPA3", "MS4A2", "KIT"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "COL3A1"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CLEC14A", "RAMP2"],
}
for name, genes in markers.items():
    g = [x for x in genes if x in adata.var_names]
    sc.tl.score_genes(adata, g, score_name=f"score_{name}")

score_cols = [f"score_{k}" for k in markers]
cl_scores = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
cl_assign = cl_scores.idxmax(axis=1).str.replace("score_", "", regex=False)
print("cluster -> compartment:\n", cl_assign)
adata.obs["compartment"] = adata.obs["leiden"].map(cl_assign).astype("category")
print(adata.obs["compartment"].value_counts())

adata.write(OUT)
print("wrote", OUT)

# ---- figures ----
sc.pl.umap(adata, color=["compartment", "leiden"], save="_compartment.png", show=False)
sc.pl.umap(adata, color=["EPCAM", "CD3D", "NKG7", "LYZ", "TACSTD2"], save="_markers.png", show=False, ncols=3)
dot_genes = ["EPCAM", "KRT19", "TACSTD2", "CD3D", "CD3E", "NKG7", "GNLY",
             "CD79A", "MS4A1", "LYZ", "CD68", "COL1A1", "DCN", "PECAM1", "VWF", "TPSAB1"]
dot_genes = [g for g in dot_genes if g in adata.var_names]
sc.pl.dotplot(adata, dot_genes, groupby="compartment", save="_compartment_markers.png", show=False)
print("figures written")
