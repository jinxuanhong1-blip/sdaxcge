#!/usr/bin/env Rscript
# =============================================================================
# OSCA / Bioconductor template — lung ICI scRNA-seq (methods only, no data)
# -----------------------------------------------------------------------------
# Covers: SingleCellExperiment, emptyDrops, MAD-based QC, scDblFinder doublets,
#         scran normalization + HVGs + denoisePCA, batchelor::fastMNN integration,
#         graph clustering, SingleR annotation, and AUCell module scoring for
#         the TACSTD2 / CLDN4 programs.
#
# `# TODO` marks every cohort-specific value. No dataset, gene list, threshold,
# or result is included. APIs reflect the current Bioconductor release (2025-2026).
# =============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DropletUtils)   # emptyDrops
  library(scater)         # QC metrics, isOutlier, runPCA/UMAP
  library(scran)          # normalization, modelGeneVar, denoisePCA, graph clustering
  library(batchelor)      # fastMNN integration
  library(scDblFinder)    # doublets
  library(SingleR)        # reference annotation
  library(AUCell)         # rank-based module scoring
  library(bluster)        # clusterCells
})

set.seed(42)

# ------------------------------------------------------------------------------
# 1. Build a SingleCellExperiment per sample; call cells with emptyDrops
# ------------------------------------------------------------------------------
sample_sheet <- read.csv("sample_sheet.csv", stringsAsFactors = FALSE)  # TODO
sce_list <- lapply(seq_len(nrow(sample_sheet)), function(i) {
  raw <- DropletUtils::read10xCounts(sample_sheet$raw_path[i])          # TODO raw (unfiltered) matrix
  ed  <- DropletUtils::emptyDrops(counts(raw))                          # ambient vs. real cells
  keep <- !is.na(ed$FDR) & ed$FDR < 0.001                               # TODO FDR
  sce <- raw[, keep]
  sce$sample_id <- sample_sheet$sample_id[i]
  sce
})
sce <- do.call(cbind, sce_list)   # requires matching rowData across samples

# ------------------------------------------------------------------------------
# 2. QC via MAD outliers (no universal constants) + per-sample doublets
# ------------------------------------------------------------------------------
is_mito <- grepl("^MT-", rowData(sce)$Symbol %||% rownames(sce))
sce <- scater::addPerCellQC(sce, subsets = list(Mito = is_mito))
# MAD-based outlier flags, computed per sample (batch = sample_id):
low_lib   <- scater::isOutlier(sce$sum,               log = TRUE, type = "lower", batch = sce$sample_id)
low_genes <- scater::isOutlier(sce$detected,          log = TRUE, type = "lower", batch = sce$sample_id)
high_mito <- scater::isOutlier(sce$subsets_Mito_percent, type = "higher",         batch = sce$sample_id)
sce <- sce[, !(low_lib | low_genes | high_mito)]

sce <- scDblFinder::scDblFinder(sce, samples = "sample_id")   # per-sample doublet calls
sce <- sce[, sce$scDblFinder.class == "singlet"]

# ------------------------------------------------------------------------------
# 3. Normalization (deconvolution size factors) + HVGs + denoised PCA
# ------------------------------------------------------------------------------
clust <- scran::quickCluster(sce)
sce   <- scran::computeSumFactors(sce, cluster = clust)
sce   <- scater::logNormCounts(sce)                       # populates logcounts assay

dec  <- scran::modelGeneVar(sce, block = sce$sample_id)   # batch-aware variance modelling
hvg  <- scran::getTopHVGs(dec, n = 2000)                  # TODO n
sce  <- scran::denoisePCA(sce, technical = dec, subset.row = hvg)  # PCs by biological variance

# ------------------------------------------------------------------------------
# 4. Batch integration with fastMNN (over sample_id, not the response variable)
# ------------------------------------------------------------------------------
mnn <- batchelor::fastMNN(sce, batch = sce$sample_id, subset.row = hvg)
reducedDim(sce, "corrected") <- reducedDim(mnn, "corrected")
# Alternative lightweight options: batchelor::rescaleBatches / regressBatches,
# or harmony::RunHarmony(sce, "sample_id", reduction = "PCA").

# ------------------------------------------------------------------------------
# 5. Clustering + UMAP on the corrected embedding
# ------------------------------------------------------------------------------
g <- scran::buildSNNGraph(sce, use.dimred = "corrected")
sce$cluster <- factor(igraph::cluster_louvain(g)$membership)   # TODO try Leiden / sweep k
sce <- scater::runUMAP(sce, dimred = "corrected")

# ------------------------------------------------------------------------------
# 6. Annotation with SingleR (cross-check against CellTypist / Azimuth)
# ------------------------------------------------------------------------------
# ref <- celldex::HumanPrimaryCellAtlasData()   # TODO or a custom lung reference SCE
# pred <- SingleR::SingleR(test = sce, ref = ref, labels = ref$label.main,
#                          assay.type.test = "logcounts")
# sce$SingleR.labels <- pred$labels
# Confirm with canonical lung TME markers via scran::scoreMarkers(sce, sce$cluster).

# ------------------------------------------------------------------------------
# 7. Module scoring: TACSTD2 / CLDN4 with AUCell (rank-based, YOU supply genes)
# ------------------------------------------------------------------------------
# Curate + document provenance in ../signatures/tacstd2_cldn4_modules.md.
tacstd2_genes <- character(0)   # TODO curated vector incl. "TACSTD2"
cldn4_genes   <- character(0)   # TODO curated vector incl. "CLDN4"
stopifnot(length(tacstd2_genes) > 0, "TACSTD2" %in% tacstd2_genes,
          length(cldn4_genes)   > 0, "CLDN4"   %in% cldn4_genes)

gene_sets <- list(TACSTD2_module = intersect(tacstd2_genes, rownames(sce)),
                  CLDN4_module   = intersect(cldn4_genes,   rownames(sce)))
rankings  <- AUCell::AUCell_buildRankings(as.matrix(logcounts(sce)), plotStats = FALSE)
auc       <- AUCell::AUCell_calcAUC(gene_sets, rankings)
sce$TACSTD2_module <- AUCell::getAUC(auc)["TACSTD2_module", ]
sce$CLDN4_module   <- AUCell::getAUC(auc)["CLDN4_module",   ]
# Sanity checks per playbook 7.4 (localization, anchor correlation, per-sample stability).

# ------------------------------------------------------------------------------
# 8. Persist + record versions
# ------------------------------------------------------------------------------
saveRDS(sce, file = "lung_ici_sce.rds")   # TODO path
writeLines(capture.output(sessionInfo()), "sessionInfo_osca.txt")
