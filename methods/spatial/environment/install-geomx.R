#!/usr/bin/env Rscript
# R/Bioconductor environment for the GeoMx track (templates 06-09).
# GeoMx 流程（模板 06-09）的 R/Bioconductor 环境。
# Target: R >= 4.3 with Bioconductor >= 3.18 (2024) or newer.
# 目标：R >= 4.3，Bioconductor >= 3.18（2024）或更新。
#
#   Rscript methods/spatial/environment/install-geomx.R

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c(
  "GeomxTools",       # GeoMxSet, QC, aggregation, Q3 norm, mixedModelDE
  "GeoMxWorkflows",   # end-to-end GeoMx helpers
  "NanoStringNCTools",
  "SpatialDecon",     # GeoMx cell-type deconvolution
  "standR",           # SpatialExperiment bridge + limma-voom workflow
  "SpatialExperiment",
  "limma", "edgeR"
)
cran_pkgs <- c("optparse", "survival", "survminer", "dplyr", "R.utils", "ggplot2")

BiocManager::install(bioc_pkgs, update = FALSE, ask = FALSE)
install.packages(setdiff(cran_pkgs, rownames(installed.packages())))

message("[install-geomx] done. Verify with: library(GeomxTools); library(SpatialDecon)")
