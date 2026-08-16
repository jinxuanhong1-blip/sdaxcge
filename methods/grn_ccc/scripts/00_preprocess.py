# %% [markdown]
# # 00 - Preprocess & define the TACSTD2/CLDN4-high epithelial state
#
# **Goal.** Turn a raw 10x lung ICI object into an analysis-ready AnnData with
# (1) ambient-RNA-aware counts, (2) doublets removed, (3) clean cell types,
# (4) malignant-vs-normal epithelium separated, and (5) a reproducible
# `epi_state` label (`TACSTD2_CLDN4_high` / `_low`).
#
# Everything downstream (LIANA, CellChat, NicheNet, SCENIC) consumes this
# object. **The quality of this step decides whether a "reduced T-cell
# recruitment signal" is biology or an artifact** (see playbook §"When these
# tools overclaim"). This is a *template*: read the NOTE comments and adapt.
#
# Run as a script (`python 00_preprocess.py`) or open as cells in Jupyter/VSCode.

# %%
from pathlib import Path
import yaml
import numpy as np
import scanpy as sc

HERE = Path(__file__).resolve().parent
CFG = yaml.safe_load((HERE.parent / "config" / "gene_sets.yaml").read_text())
P = CFG["params"]
sc.settings.verbosity = 1

# ---- EDIT THESE PATHS -------------------------------------------------------
RAW_H5AD = "data/lung_ici_raw.h5ad"          # cells x genes, RAW counts in .X
OUT_H5AD = "data/lung_ici_preprocessed.h5ad"
# -----------------------------------------------------------------------------

# %% [markdown]
# ## 1. Ambient RNA removal (do this FIRST, on raw counts)
#
# Dissociated solid lung leaks highly-expressed epithelial transcripts
# (surfactant, keratins, and TACSTD2/CLDN4 themselves) into every droplet.
# Uncorrected, this manufactures "epithelial ligand" signal in non-epithelial
# cells AND inflates the epithelial compartment. Correct BEFORE clustering.
#
# Preferred: run **CellBender** (raw+filtered matrices) or **SoupX** upstream,
# or **decontX** (below) if you only have a filtered matrix. Then treat the
# corrected matrix as the new raw counts.

# %%
try:
    import celltypist  # noqa: F401  (placeholder; not required)
except Exception:
    pass

adata = sc.read_h5ad(RAW_H5AD)
adata.layers["counts_raw"] = adata.X.copy()

# decontX via scanpy-external if available; otherwise assume upstream CellBender/SoupX.
try:
    import scanpy.external as sce
    sce.pp.scrublet  # noqa: B018  (import probe)
except Exception:
    pass

# NOTE: decontX lives in the R 'celda' package or the python 'decontx' wrapper.
# If you cannot run it here, run CellBender/SoupX upstream and load the cleaned
# matrix as RAW_H5AD. Do not skip ambient correction silently -- record it.
adata.uns["ambient_correction"] = "ASSUMED_UPSTREAM (CellBender/SoupX). EDIT ME."

# %% [markdown]
# ## 2. QC + doublet removal
# Epithelial-T doublets create chimeric expression -> phantom LR pairs.

# %%
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
adata = adata[
    (adata.obs.n_genes_by_counts > 200)
    & (adata.obs.n_genes_by_counts < 8000)
    & (adata.obs.pct_counts_mt < 20)
].copy()

try:
    sc.pp.scrublet(adata, random_state=P["random_seed"])
    adata = adata[~adata.obs["predicted_doublet"]].copy()
except Exception as e:  # scrublet optional
    print(f"[warn] scrublet skipped ({e}); use scDblFinder upstream.")

# %% [markdown]
# ## 3. Normalize (keep raw counts around for pseudobulk & SCENIC)

# %%
adata.layers["counts"] = adata.X.copy()          # post-QC integer counts
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata                                # log-norm snapshot for CCC tools

# %% [markdown]
# ## 4. Cluster & annotate
# Use your preferred reference/marker workflow. IMPORTANT: run clustering on a
# batch-corrected embedding, but keep `adata.X` / `adata.raw` as *uncorrected*
# log-norm values -- CCC and DE tools must not see Harmony/scVI-imputed values.

# %%
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key=P.get("sample_key"))
sc.pp.pca(adata, n_comps=30, mask_var="highly_variable")
# Batch integration for the *embedding only* (e.g. Harmony); optional:
try:
    sc.external.pp.harmony_integrate(adata, key=P["sample_key"])
    rep = "X_pca_harmony"
except Exception:
    rep = "X_pca"
sc.pp.neighbors(adata, use_rep=rep)
sc.tl.leiden(adata, resolution=1.0, random_state=P["random_seed"])

# NOTE: map leiden clusters -> cell_type using lineage_markers from the config.
# Fill this dict after inspecting sc.tl.rank_genes_groups / dotplots.
leiden_to_celltype = {}  # e.g. {"0": "T_cell", "1": "Epithelial", ...}
if leiden_to_celltype:
    adata.obs[P["celltype_key"]] = (
        adata.obs["leiden"].map(leiden_to_celltype).astype("category")
    )
else:
    print("[action needed] fill leiden_to_celltype before continuing.")

# %% [markdown]
# ## 5. Separate malignant vs normal epithelium
# TACSTD2/CLDN4 are high in *malignant* epithelium but also in reactive normal
# epithelium. Run inferCNV / CopyKAT (R) or numbat to call malignant cells so
# the "state" is not confounded by normal alveolar/airway epithelium.

# %%
# Placeholder: expect adata.obs['malignant'] in {True, False} from inferCNV/CopyKAT.
if "malignant" not in adata.obs:
    adata.obs["malignant"] = np.nan
    print("[action needed] add adata.obs['malignant'] from inferCNV/CopyKAT.")

# %% [markdown]
# ## 6. Define the TACSTD2/CLDN4-high state (reproducibly)
# Score epithelium on the positive markers, then split by quantile (default),
# Otsu, or a dedicated cluster. Record the method + threshold in `.uns`.

# %%
def define_epi_state(ad, celltype_key):
    epi_mask = ad.obs[celltype_key].astype(str).str.contains("pithel")
    epi = ad[epi_mask].copy()
    pos = [g for g in CFG["state_markers"]["positive"] if g in epi.var_names]
    sc.tl.score_genes(epi, pos, score_name="tacstd2_cldn4_score",
                      random_state=P["random_seed"])
    method = P["state_score_method"]
    s = epi.obs["tacstd2_cldn4_score"].to_numpy()
    if method == "quantile":
        thr = np.quantile(s, P["state_high_quantile"])
    elif method == "otsu":
        from skimage.filters import threshold_otsu
        thr = threshold_otsu(s)
    else:  # cluster: expect a precomputed sub-cluster label
        thr = None
    label = np.where(s >= thr, "TACSTD2_CLDN4_high", "TACSTD2_CLDN4_low")
    ad.obs["epi_state"] = "non_epithelial"
    ad.obs.loc[epi.obs_names, "epi_state"] = label
    ad.obs["epi_state"] = ad.obs["epi_state"].astype("category")
    ad.uns["epi_state_def"] = {"method": method, "threshold": float(thr),
                               "markers": pos}
    return ad

if P["celltype_key"] in adata.obs:
    adata = define_epi_state(adata, P["celltype_key"])
    print(adata.obs["epi_state"].value_counts())

# %% [markdown]
# ## 7. Save
# Downstream scripts read OUT_H5AD. Also emit a small metadata report so the
# provenance of `epi_state` (and whether ambient correction ran) is auditable.

# %%
Path(OUT_H5AD).parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(OUT_H5AD)
print(f"[done] wrote {OUT_H5AD}")
print("ambient_correction:", adata.uns.get("ambient_correction"))
print("epi_state_def:", adata.uns.get("epi_state_def"))
