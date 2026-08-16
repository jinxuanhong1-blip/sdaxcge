"""Annotation: marker-based coarse lineage + CellTypist automated calls.

Reconcile automated labels with markers and cluster DE before trusting them.
自动注释需与 marker、簇DE互证后方可采信。

Usage: python 04_annotate.py integrated.h5ad out.h5ad
"""
import sys
import scanpy as sc
from gene_sets import LINEAGE_MARKERS

def main(in_h5ad, out_h5ad):
    adata = sc.read_h5ad(in_h5ad)

    # 1) marker-based per-cluster scoring -> assign lineage by argmax
    use_raw = adata.raw is not None
    for name, gs in LINEAGE_MARKERS.items():
        present = [g for g in gs if g in (adata.raw.var_names if use_raw
                                          else adata.var_names)]
        sc.tl.score_genes(adata, present, score_name=f"score_{name}", use_raw=use_raw)
    score_cols = [f"score_{k}" for k in LINEAGE_MARKERS]
    per_cluster = adata.obs.groupby("leiden")[score_cols].mean()
    lab = per_cluster.idxmax(axis=1).str.replace("score_", "", regex=False)
    adata.obs["lineage"] = adata.obs["leiden"].map(lab).astype("category")

    # 2) CellTypist (optional automated prior for immune)
    try:
        import celltypist
        from celltypist import models
        models.download_models(model=["Immune_All_Low.pkl"])
        pred = celltypist.annotate(adata, model="Immune_All_Low.pkl",
                                   majority_voting=True)
        adata.obs["celltypist"] = pred.predicted_labels["majority_voting"].values
    except Exception as exc:  # noqa: BLE001
        print(f"CellTypist skipped: {exc}")

    # 3) cluster marker DE for the record (Wilcoxon here is fine: marker discovery)
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

    adata.write(out_h5ad)
    print("Lineage counts:\n", adata.obs["lineage"].value_counts())

if __name__ == "__main__":
    main(*sys.argv[1:])
