#!/usr/bin/env Rscript
# 08 · GeoMx cell-type deconvolution with SpatialDecon. GeoMx 细胞类型去卷积。
# ---------------------------------------------------------------------------
# GeoMx WTA AOIs are multi-cellular. SpatialDecon estimates cell-type abundances
# per AOI from a cell-profile (signature) matrix, using the negative-probe
# background as the noise model. Pair with a lung/TME reference profile matrix
# (e.g. a SafeTME-style or study-matched scRNA-derived signature).
# GeoMx WTA 的 AOI 为多细胞。SpatialDecon 基于细胞特征（signature）矩阵，以阴性
# 探针背景为噪声模型，估计每个 AOI 的细胞类型丰度。需配套肺/TME 参考特征矩阵。
#
# Tools: SpatialDecon, GeomxTools (Bioconductor 3.18+, 2024+).
# No results are fabricated; abundances come from the model.
# 工具：SpatialDecon、GeomxTools。不伪造结果；丰度由模型估计。
#
# Usage:
#   Rscript 08_geomx_deconvolution_spatialdecon.R --rds <geomx_target_qnorm.rds> \
#       --profile <cell_profile_matrix.csv> --out methods/spatial/demo/out/geomx
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(GeomxTools)
  library(SpatialDecon)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--rds", type = "character", help = "output of template 06"),
  make_option("--profile", type = "character",
              help = "genes x cell-types signature matrix (CSV); supply a lung/TME profile"),
  make_option("--out", type = "character", default = "methods/spatial/demo/out/geomx")
)))

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
target_geomx <- readRDS(opt$rds)

# Background model from negative probes. 由阴性探针估计背景。
bg <- derive_GeoMx_background(
  norm = assayDataElement(target_geomx, "q_norm"),
  probepool = fData(target_geomx)$Module,
  negnames = fData(target_geomx)$TargetName[fData(target_geomx)$Negative]
)

profile <- as.matrix(read.csv(opt$profile, row.names = 1, check.names = FALSE))

res <- runspatialdecon(
  object = target_geomx,
  norm_elt = "q_norm",
  raw_elt = "exprs",
  X = profile,
  align_genes = TRUE
)

# Cell abundances (beta) per AOI. 每个 AOI 的细胞丰度（beta）。
beta <- t(res$beta)
write.csv(beta, file.path(opt$out, "spatialdecon_abundance.csv"))
saveRDS(res, file.path(opt$out, "spatialdecon_result.rds"))
message(sprintf("[08_geomx] deconvolved %d AOIs into %d cell types -> %s",
                nrow(beta), ncol(beta), opt$out))
