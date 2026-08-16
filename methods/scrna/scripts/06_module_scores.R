## Module scoring with UCell (rank-based, robust to depth/normalization) + Seurat.
## Recommended over AddModuleScore for cross-sample comparison. See playbook.md 6.2.
## UCell 基于排序，对深度/归一化更稳健，跨样本比较优于 AddModuleScore。
suppressPackageStartupMessages({library(Seurat); library(UCell)})
source(file.path(dirname(sys.frame(1)$ofile %||% "."), "gene_sets.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE)
in_rds  <- args[[1]]
out_rds <- args[[2]]

if (!exists("gene_sets")) source("gene_sets.R")
seu <- readRDS(in_rds)

## UCell works on the RNA counts/data assay; scores are per-cell [0,1].
seu <- AddModuleScore_UCell(seu, features = gene_sets, name = "_UCell")

## Compare at SAMPLE level (e.g. epithelial cells only for the tumor module):
## epi <- subset(seu, subset = lineage == "Epithelial")
## agg <- aggregate(TACSTD2_CLDN4_junction_UCell ~ sample, epi@meta.data, mean)
saveRDS(seu, out_rds)
message("UCell module scores -> ", out_rds)
