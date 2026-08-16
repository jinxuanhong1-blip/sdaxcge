#!/usr/bin/env Rscript
# Extract per-cell TACSTD2 / EPCAM / PTPRC UMI and total UMI from the
# author-processed GSE205335 dgCMatrix. The GEO .rds.gz is double-gzipped:
# gunzip once, then readRDS() consumes the remaining gzip layer.

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
rds_path <- if (length(args) >= 1) args[[1]] else "/workspace/data/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds"
out_path <- if (length(args) >= 2) args[[2]] else "/workspace/data/GSE205335/percell_tacstd2.csv.gz"

cat("Reading RDS:", rds_path, "\n")
obj <- readRDS(rds_path)
cat("Class:", paste(class(obj), collapse = ","), "\n")

if (inherits(obj, "list")) {
  cat("List names:", paste(names(obj), collapse = ", "), "\n")
  for (el in obj) if (inherits(el, "sparseMatrix")) { obj <- el; break }
}

stopifnot(inherits(obj, "sparseMatrix") || is.matrix(obj))
cat("Dims:", nrow(obj), "genes x", ncol(obj), "cells\n")

genes <- rownames(obj)
need <- c("TACSTD2", "EPCAM", "PTPRC")
idx <- match(need, genes)
if (any(is.na(idx))) {
  cat("Missing genes:\n")
  print(need[is.na(idx)])
  cat("TACSTD candidates:\n")
  print(grep("TACSTD", genes, value = TRUE))
  quit(status = 1)
}

tacstd2 <- as.numeric(obj[idx[1], ])
epcam <- as.numeric(obj[idx[2], ])
ptprc <- as.numeric(obj[idx[3], ])
total <- as.numeric(Matrix::colSums(obj))

df <- data.frame(
  barcode = colnames(obj),
  tacstd2_umi = tacstd2,
  epcam_umi = epcam,
  ptprc_umi = ptprc,
  total_umi = total
)
con <- gzfile(out_path, "w")
write.csv(df, con, row.names = FALSE)
close(con)
cat("Wrote", out_path, "with", nrow(df), "cells\n")
cat("TACSTD2 UMI>0 fraction:", mean(tacstd2 > 0), "\n")
