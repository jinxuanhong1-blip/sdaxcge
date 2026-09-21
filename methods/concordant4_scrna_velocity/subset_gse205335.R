#!/usr/bin/env Rscript
# Subset the public GSE205335 UMI RDS to epithelial cells and write MatrixMarket.
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
rds <- "/tmp/c4data/GSE205335.rds"
barcodes <- "/tmp/c4data/GSE205335_epi_barcodes.txt"
message("readRDS")
m <- readRDS(rds)
if (!inherits(m, "dgCMatrix")) m <- as(m, "CsparseMatrix")
message("dim ", nrow(m), " x ", ncol(m))
want <- readLines(barcodes)
want <- want[nzchar(want)]
hit <- intersect(want, colnames(m))
message("epithelial barcodes requested ", length(want), " found ", length(hit))
sub <- m[, hit, drop = FALSE]
out <- "/tmp/c4data/GSE205335_epi.mtx"
writeMM(sub, out)
writeLines(rownames(sub), "/tmp/c4data/GSE205335_epi_genes.txt")
writeLines(colnames(sub), "/tmp/c4data/GSE205335_epi_cells.txt")
message("wrote ", out, " ", nrow(sub), " x ", ncol(sub))
