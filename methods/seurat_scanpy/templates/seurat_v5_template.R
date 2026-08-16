#!/usr/bin/env Rscript
# =============================================================================
# Seurat v5 template — lung ICI scRNA-seq (methods only, no data included)
# -----------------------------------------------------------------------------
# Covers: BPCells on-disk counts, layers/assays, QC, normalization, HVG/PCA,
#         SketchData + ProjectData, IntegrateLayers (5 methods), Azimuth (lungref),
#         and TACSTD2 / CLDN4 module scoring.
#
# Everything marked `# TODO` must be supplied for YOUR cohort. This script does
# not ship any dataset, gene list, threshold, or result. APIs reflect Seurat v5
# (2025-2026). Run section by section; do not treat it as a turnkey pipeline.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)        # v5.x
  library(SeuratObject)  # v5.x
  library(BPCells)       # on-disk matrices
  library(dplyr)
})

set.seed(42)
options(future.globals.maxSize = 8 * 1024^3)  # TODO: tune for your machine

# ------------------------------------------------------------------------------
# 1. Ingest each sample into BPCells (on-disk) and build a Seurat object
# ------------------------------------------------------------------------------
# sample_sheet: a data.frame you provide with at least columns:
#   sample_id, h5_path (10x CellRanger filtered_feature_bc_matrix.h5), + covariates
sample_sheet <- read.csv("sample_sheet.csv", stringsAsFactors = FALSE)  # TODO
bpcells_root <- "bpcells"                                               # TODO writable dir
dir.create(bpcells_root, showWarnings = FALSE)

counts_list <- list()
for (i in seq_len(nrow(sample_sheet))) {
  sid  <- sample_sheet$sample_id[i]
  dest <- file.path(bpcells_root, sid)
  if (!dir.exists(dest)) {
    mat <- open_matrix_10x_hdf5(path = sample_sheet$h5_path[i])  # streams from disk
    write_matrix_dir(mat = mat, dir = dest)                       # one-time conversion
  }
  counts_list[[sid]] <- open_matrix_dir(dir = dest)               # disk-backed pointer
}

# Merged, disk-backed counts -> single Seurat object; then split into per-sample layers.
obj <- CreateSeuratObject(counts = counts_list, project = "lung_ici")  # counts stay on disk
# attach per-cell covariates from the sample sheet (join on orig.ident == sample_id)
meta <- sample_sheet
rownames(meta) <- meta$sample_id
obj$sample_id <- obj$orig.ident
for (col in setdiff(colnames(meta), "sample_id")) {
  obj[[col]] <- meta[obj$sample_id, col]                # TODO: e.g. patient, timepoint, response
}

# ------------------------------------------------------------------------------
# 2. QC (derive thresholds from YOUR distributions; nothing is hard-coded here)
# ------------------------------------------------------------------------------
obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
obj[["percent.hb"]]   <- PercentageFeatureSet(obj, pattern = "^HB[^(P)]")

# Inspect BEFORE filtering (write plots, look at per-sample distributions):
# VlnPlot(obj, c("nFeature_RNA","nCount_RNA","percent.mt"), group.by="sample_id")
# Prefer MAD-based, per-sample cutoffs over universal constants:
#   e.g. via scater::isOutlier on nCount/nFeature/percent.mt (see OSCA template).
min_features <- NA_real_   # TODO
max_features <- NA_real_   # TODO
max_mt       <- NA_real_   # TODO (cohort-specific; tumor/epithelial MT is often high)
stopifnot(!is.na(min_features), !is.na(max_features), !is.na(max_mt))  # force explicit choice
obj <- subset(obj,
              subset = nFeature_RNA > min_features &
                       nFeature_RNA < max_features &
                       percent.mt   < max_mt)

# Doublets: run per sample. scDblFinder operates on an SCE; convert per sample,
# call scDblFinder::scDblFinder(sce, samples=...), then carry the calls back.
# (Left as a documented step; see OSCA template for the SCE-side code.)

# ------------------------------------------------------------------------------
# 3. Split into per-sample layers -> normalize -> HVG -> PCA
# ------------------------------------------------------------------------------
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample_id)   # creates counts.<sample> layers
obj <- NormalizeData(obj)                                # LogNormalize (streams over disk)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)  # TODO nfeatures
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)                            # TODO npcs

# ------------------------------------------------------------------------------
# 4. (Large cohorts) Sketch: leverage-score subsample -> analyze in memory
# ------------------------------------------------------------------------------
obj <- SketchData(obj, ncells = 5000, method = "LeverageScore", sketched.assay = "sketch")  # TODO ncells
DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)

# ------------------------------------------------------------------------------
# 5. Integration on the (sketch) layers: pick ONE method (verified Seurat v5 API)
# ------------------------------------------------------------------------------
# Integrate over the technical batch (usually sample/patient), NOT the response var.
obj <- IntegrateLayers(obj, method = RPCAIntegration,
                       orig.reduction = "pca", new.reduction = "integrated.rpca", verbose = FALSE)
# Alternatives (choose one; see playbook 5.1):
# obj <- IntegrateLayers(obj, method = CCAIntegration,    orig.reduction="pca", new.reduction="integrated.cca")
# obj <- IntegrateLayers(obj, method = HarmonyIntegration, orig.reduction="pca", new.reduction="harmony")
# obj <- IntegrateLayers(obj, method = FastMNNIntegration, orig.reduction="pca", new.reduction="integrated.mnn")
# obj <- IntegrateLayers(obj, method = scVIIntegration,   new.reduction="integrated.scvi",
#                        conda_env = "/path/to/scvi-env")   # TODO scvi-tools env
integrated_reduction <- "integrated.rpca"   # TODO keep in sync with the method above

# ------------------------------------------------------------------------------
# 6. Cluster + UMAP on the integrated reduction (sketch), then project to full
# ------------------------------------------------------------------------------
obj <- FindNeighbors(obj, reduction = integrated_reduction, dims = 1:30)
obj <- FindClusters(obj, resolution = 1.0, cluster.name = "sketch_clusters")  # TODO sweep resolution
obj <- RunUMAP(obj, reduction = integrated_reduction, dims = 1:30, return.model = TRUE)

# Extend sketch results to the full on-disk dataset:
obj <- ProjectData(obj, assay = "RNA", sketched.assay = "sketch",
                   sketched.reduction = integrated_reduction,
                   full.reduction = paste0(integrated_reduction, ".full"),
                   umap.model = "umap", dims = 1:30,
                   refdata = list(cluster_full = "sketch_clusters"))
DefaultAssay(obj) <- "RNA"   # back to full, on-disk data

# ------------------------------------------------------------------------------
# 7. Annotation — Azimuth lung reference (verified: lungref v2.0.0 == HLCA v2)
# ------------------------------------------------------------------------------
# Requires: library(Azimuth); library(SeuratData); SeuratData::InstallData("lungref")
# obj <- Azimuth::RunAzimuth(obj, reference = "lungref")
#   -> adds predicted.ann_level_*, predicted.ann_finest_level, mapping.score, prediction.score.*
# Cross-check with SingleR (OSCA template) and manual marker curation before trusting labels.

# ------------------------------------------------------------------------------
# 8. Module scoring: TACSTD2 / CLDN4 (YOU supply curated gene sets; none invented)
# ------------------------------------------------------------------------------
# Curate + document provenance in ../signatures/tacstd2_cldn4_modules.md.
# Each vector must include its anchor gene and only genes with a citable/derivable source.
tacstd2_genes <- character(0)  # TODO: curated vector, must include "TACSTD2"
cldn4_genes   <- character(0)  # TODO: curated vector, must include "CLDN4"
stopifnot(length(tacstd2_genes) > 0, "TACSTD2" %in% tacstd2_genes,
          length(cldn4_genes)   > 0, "CLDN4"   %in% cldn4_genes)

obj <- JoinLayers(obj)                       # collapse per-sample layers before scoring/DE
modules <- list(TACSTD2_module = intersect(tacstd2_genes, rownames(obj)),
                CLDN4_module   = intersect(cldn4_genes,   rownames(obj)))
obj <- AddModuleScore(obj, features = modules, name = names(modules), ctrl = 100, seed = 42)
# -> meta.data: TACSTD2_module1, CLDN4_module2 (Seurat appends a positional index)

# Recommended rank-based alternative (robust to depth/normalization):
# library(UCell)
# obj <- AddModuleScore_UCell(obj, features = modules)   # -> *_UCell columns

# Sanity checks (playbook 7.4): score should localize to epithelial/malignant cells,
# correlate with the anchor gene, be stable per-sample, and not merely track nCount_RNA.

# ------------------------------------------------------------------------------
# 9. Persist + record versions
# ------------------------------------------------------------------------------
saveRDS(obj, file = "lung_ici_seurat_v5.rds")  # TODO path (BPCells matrices stay on disk)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
