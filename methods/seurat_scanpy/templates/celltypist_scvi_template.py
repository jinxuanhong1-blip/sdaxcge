#!/usr/bin/env python
# =============================================================================
# CellTypist + scVI/scANVI template — lung ICI scRNA-seq (methods only, no data)
# -----------------------------------------------------------------------------
# Covers: scVI integration (probabilistic latent space), scANVI label transfer,
#         and CellTypist automated annotation (immune + Human_Lung_Atlas models).
#
# Assumes you already have an AnnData with raw counts in adata.layers["counts"]
# and log1p CPM-10k in adata.X (see scanpy_template.py steps 1-3). Every `# TODO`
# is cohort-specific. No dataset, gene list, threshold, or result is included.
# APIs reflect scvi-tools + celltypist >= 1.7 (2025-2026).
# =============================================================================
from __future__ import annotations

import scanpy as sc
import scvi
import celltypist
from celltypist import models

RANDOM_STATE = 0
scvi.settings.seed = RANDOM_STATE

adata = sc.read_h5ad("lung_ici_scanpy.h5ad")   # TODO from scanpy_template.py
assert "counts" in adata.layers, "scVI needs raw counts in adata.layers['counts']."

# ------------------------------------------------------------------------------
# 1. scVI integration — batch = sample_id, optional covariates
# ------------------------------------------------------------------------------
# Train on HVGs computed from counts; restrict adata to HVGs before setup if desired.
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="sample_id",                       # technical batch (NOT response)
    # continuous_covariate_keys=["pct_counts_mt"],   # TODO optional nuisance covariates
    # categorical_covariate_keys=["chemistry"],      # TODO optional
)
vae = scvi.model.SCVI(adata, n_latent=30)          # TODO n_latent
vae.train(max_epochs=None)                          # TODO epochs / early stopping
adata.obsm["X_scVI"] = vae.get_latent_representation()

# Use X_scVI as the integrated embedding downstream:
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata, random_state=RANDOM_STATE)
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False, random_state=RANDOM_STATE)

# ------------------------------------------------------------------------------
# 2. (Optional) scANVI — semi-supervised label transfer from partial labels
# ------------------------------------------------------------------------------
# If you have trusted labels for some cells (e.g. from Azimuth/CellTypist/manual),
# put them in adata.obs["seed_labels"] with "Unknown" for the rest, then:
# lvae = scvi.model.SCANVI.from_scvi_model(vae, adata=adata,
#                                          labels_key="seed_labels", unlabeled_category="Unknown")
# lvae.train(max_epochs=None)
# adata.obs["scANVI_labels"] = lvae.predict()
# adata.obsm["X_scANVI"] = lvae.get_latent_representation()

# ------------------------------------------------------------------------------
# 3. CellTypist annotation (verified API + model names)
# ------------------------------------------------------------------------------
# Input scale requirement: log1p normalized to 10,000 counts/cell (adata.X already is).
models.download_models(model=["Immune_All_Low.pkl", "Human_Lung_Atlas.pkl"])  # TODO choose
# models.models_description()   # list all available models

pred_immune = celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True)
adata.obs["celltypist_immune"] = pred_immune.predicted_labels["majority_voting"].values
adata.obs["celltypist_immune_conf"] = pred_immune.predicted_labels["conf_score"].values

pred_lung = celltypist.annotate(adata, model="Human_Lung_Atlas.pkl", majority_voting=True)
adata.obs["celltypist_lung"] = pred_lung.predicted_labels["majority_voting"].values
adata.obs["celltypist_lung_conf"] = pred_lung.predicted_labels["conf_score"].values

# Reconcile CellTypist (immune + lung), scANVI, and marker curation into a final
# consensus label; send low-confidence / disagreeing clusters to manual review.  # TODO

# ------------------------------------------------------------------------------
# 4. Persist + record versions
# ------------------------------------------------------------------------------
adata.write_h5ad("lung_ici_scvi_celltypist.h5ad")   # TODO path
sc.logging.print_versions()
