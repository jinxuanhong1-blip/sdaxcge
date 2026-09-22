#!/usr/bin/env Rscript
# Peel double-gzipped GSE205335 UMI RDS and write a gene-subset Matrix Market
# plus full-matrix library sizes. The RDS is a dgCMatrix, not a Seurat object.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: export_gse205335.R <rds.gz> <genes.txt> <out_dir>")
rds_gz <- args[[1]]
genes_file <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

is_gzip <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 2)
  length(magic) == 2 && magic[[1]] == as.raw(0x1f) && magic[[2]] == as.raw(0x8b)
}
cur <- rds_gz
temps <- character()
on.exit(unlink(temps[file.exists(temps)]), add = TRUE)
for (i in seq_len(4)) {
  if (!is_gzip(cur)) break
  dest <- tempfile(pattern = paste0("gse205335_", i, "_"), fileext = ".rds")
  temps <- c(temps, dest)
  message("gzip -dc layer ", i)
  st <- system2("gzip", c("-dc", cur), stdout = dest)
  if (!is.null(st) && !identical(st, 0L)) stop("gzip failed on layer ", i)
  cur <- dest
}
message("readRDS ", cur, " bytes ", file.info(cur)$size)
mat <- readRDS(cur)
if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "CsparseMatrix")
message("full ", nrow(mat), " x ", ncol(mat))
wanted <- unique(toupper(scan(genes_file, what = character(), quiet = TRUE)))
rn <- toupper(rownames(mat))
keep <- rn %in% wanted
if (!any(keep)) stop("none of the requested genes are in the RDS")
lib <- as.numeric(Matrix::colSums(mat))
cells <- colnames(mat)
sub <- mat[keep, , drop = FALSE]
rownames(sub) <- rn[keep]
rm(mat)
gc()
if (anyDuplicated(rownames(sub))) sub <- sub[!duplicated(rownames(sub)), , drop = FALSE]
Matrix::writeMM(sub, file.path(out_dir, "matrix.mtx"))
writeLines(rownames(sub), file.path(out_dir, "genes.txt"))
writeLines(cells, file.path(out_dir, "cells.txt"))
utils::write.table(
  data.frame(cell = cells, libsize = lib),
  file.path(out_dir, "libsize.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("done ", nrow(sub), " x ", ncol(sub))
