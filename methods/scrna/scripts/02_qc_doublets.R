## Per-sample doublet detection (scDblFinder) + MAD-based QC.
## Run BEFORE integration, on each sample separately. 整合前、每样本单独运行。
suppressPackageStartupMessages({
  library(scDblFinder); library(SingleCellExperiment); library(scater); library(scran)
})
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
in_rds  <- args[[1]]     # SingleCellExperiment or Seurat counts for ONE sample
out_rds <- args[[2]]

sce <- readRDS(in_rds)
if (inherits(sce, "Seurat")) sce <- Seurat::as.SingleCellExperiment(sce)

## --- MAD-based per-sample QC (log scale) ---
sce <- addPerCellQC(sce, subsets = list(
  mt = grep("^MT-", rownames(sce), ignore.case = TRUE)))
qc_lib   <- isOutlier(sce$sum,        log = TRUE,  nmads = 3, type = "lower")
qc_ngene <- isOutlier(sce$detected,   log = TRUE,  nmads = 3, type = "lower")
qc_mt    <- isOutlier(sce$subsets_mt_percent,      nmads = 3, type = "higher")
discard  <- qc_lib | qc_ngene | qc_mt
message(sprintf("QC discards: %d / %d cells", sum(discard), ncol(sce)))
sce <- sce[, !discard]

## --- doublets ---
sce <- scDblFinder(sce)                      # adds sce$scDblFinder.class / .score
message("Doublets called: ",
        sum(sce$scDblFinder.class == "doublet"), " / ", ncol(sce))
sce <- sce[, sce$scDblFinder.class == "singlet"]

saveRDS(sce, out_rds)
message("Clean singlets -> ", out_rds)
