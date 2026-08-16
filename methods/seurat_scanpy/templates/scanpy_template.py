#!/usr/bin/env python
# =============================================================================
# Scanpy template — lung ICI scRNA-seq (methods only, no data included)
# -----------------------------------------------------------------------------
# Covers: AnnData layers (.X vs .layers["counts"]), per-sample QC + scrublet,
#         normalization, batch-aware HVGs, PCA/neighbors/UMAP/Leiden (igraph
#         flavor), Harmony integration, and TACSTD2 / CLDN4 module scoring with
#         sc.tl.score_genes.
#
# Every `# TODO` must be supplied for YOUR cohort. No dataset, gene list,
# threshold, or result is included. APIs reflect scanpy >= 1.10 (2025-2026).
# scVI/scANVI integration + CellTypist annotation live in
# celltypist_scvi_template.py.
# =============================================================================
from __future__ import annotations

import numpy as np
import scanpy as sc
import scanpy.external as sce
import anndata as ad

sc.settings.verbosity = 1
RANDOM_STATE = 0

# ------------------------------------------------------------------------------
# 1. Ingest each sample; keep raw counts in a named layer from the start
# ------------------------------------------------------------------------------
# sample_sheet: list of dicts you provide, e.g.
# sample_sheet = [{"sample_id": "S1", "path": ".../filtered_feature_bc_matrix.h5",
#                  "patient": "P1", "timepoint": "pre", "response": "R"}, ...]  # TODO
sample_sheet: list[dict] = []  # TODO
adatas = []
for s in sample_sheet:
    a = sc.read_10x_h5(s["path"])
    a.var_names_make_unique()
    a.obs["sample_id"] = s["sample_id"]
    for k, v in s.items():
        if k not in ("path",):
            a.obs[k] = v  # TODO covariates: patient, timepoint, response, site, chemistry
    adatas.append(a)
adata = ad.concat(adatas, join="outer", label="batch", index_unique="-")
adata.layers["counts"] = adata.X.copy()   # <-- immutable raw counts; scVI etc. read this

# ------------------------------------------------------------------------------
# 2. QC (derive thresholds from YOUR per-sample distributions)
# ------------------------------------------------------------------------------
adata.var["mt"]   = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"]   = adata.var_names.str.contains(r"^HB[^P]")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=None)

# Inspect first (sc.pl.violin per sample), then set explicit cutoffs:
MIN_GENES = None   # TODO
MAX_GENES = None   # TODO
MAX_PCT_MT = None  # TODO (tumor/epithelial MT often elevated; set per cohort)
assert None not in (MIN_GENES, MAX_GENES, MAX_PCT_MT), "Set QC thresholds from your data."
adata = adata[
    (adata.obs["n_genes_by_counts"] > MIN_GENES)
    & (adata.obs["n_genes_by_counts"] < MAX_GENES)
    & (adata.obs["pct_counts_mt"] < MAX_PCT_MT)
].copy()

# Doublets: run per sample, never on the merged object.
for sid in adata.obs["sample_id"].unique():
    m = adata.obs["sample_id"] == sid
    sub = adata[m].copy()
    sc.pp.scrublet(sub, random_state=RANDOM_STATE)
    adata.obs.loc[m, "predicted_doublet"] = sub.obs["predicted_doublet"].values
    adata.obs.loc[m, "doublet_score"] = sub.obs["doublet_score"].values
adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()

# ------------------------------------------------------------------------------
# 3. Normalize (log1p CPM-10k) -> batch-aware HVGs -> PCA
# ------------------------------------------------------------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)                       # adata.X now log1p CPM-10k (CellTypist-compatible)
adata.raw = adata                        # freeze full log-normalized state (optional)

sc.pp.highly_variable_genes(
    adata, n_top_genes=2000, flavor="seurat_v3",
    layer="counts", batch_key="sample_id",   # batch-aware selection
)  # TODO n_top_genes
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.pp.pca(adata, n_comps=50, random_state=RANDOM_STATE)   # TODO n_comps

# ------------------------------------------------------------------------------
# 4. Batch integration (pick ONE; integrate over sample, not response)
# ------------------------------------------------------------------------------
sce.pp.harmony_integrate(adata, key="sample_id")   # -> adata.obsm["X_pca_harmony"]
use_rep = "X_pca_harmony"
# Alternatives:
# sce.pp.bbknn(adata, batch_key="sample_id")                 # modifies neighbor graph
# sce.pp.scanorama_integrate(adata, key="sample_id"); use_rep = "X_scanorama"
# scVI/scANVI: see celltypist_scvi_template.py -> use_rep = "X_scVI"

# ------------------------------------------------------------------------------
# 5. Neighbors + UMAP + Leiden (igraph flavor, recommended in current scanpy)
# ------------------------------------------------------------------------------
sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15)   # TODO n_neighbors
sc.tl.umap(adata, random_state=RANDOM_STATE)
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False,
             random_state=RANDOM_STATE)   # TODO sweep resolution

# ------------------------------------------------------------------------------
# 6. Marker-based confirmation (annotation via CellTypist -> other template)
# ------------------------------------------------------------------------------
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
# Reconcile with canonical lung TME markers (you supply the panels) using dotplots.

# ------------------------------------------------------------------------------
# 7. Module scoring: TACSTD2 / CLDN4 (YOU supply curated gene sets; none invented)
# ------------------------------------------------------------------------------
# Curate + document provenance in ../signatures/tacstd2_cldn4_modules.md.
TACSTD2_MODULE: list[str] = []  # TODO curated list, must include "TACSTD2"
CLDN4_MODULE:   list[str] = []  # TODO curated list, must include "CLDN4"
assert "TACSTD2" in TACSTD2_MODULE and "CLDN4" in CLDN4_MODULE, "Include anchor genes."

# Score on the full log-normalized matrix (use .raw if HVG-subset above):
score_adata = adata.raw.to_adata() if adata.raw is not None else adata
tacstd2_genes = [g for g in TACSTD2_MODULE if g in score_adata.var_names]
cldn4_genes   = [g for g in CLDN4_MODULE   if g in score_adata.var_names]
sc.tl.score_genes(score_adata, tacstd2_genes, score_name="TACSTD2_module",
                  ctrl_size=100, random_state=RANDOM_STATE)
sc.tl.score_genes(score_adata, cldn4_genes, score_name="CLDN4_module",
                  ctrl_size=100, random_state=RANDOM_STATE)
adata.obs["TACSTD2_module"] = score_adata.obs["TACSTD2_module"].values
adata.obs["CLDN4_module"]   = score_adata.obs["CLDN4_module"].values
# Sanity checks (playbook 7.4): localization to epithelial/malignant cells,
# correlation with the anchor gene, per-sample stability, independence from total counts.

# ------------------------------------------------------------------------------
# 8. Persist + record versions
# ------------------------------------------------------------------------------
adata.write_h5ad("lung_ici_scanpy.h5ad")   # TODO path
sc.logging.print_versions()
