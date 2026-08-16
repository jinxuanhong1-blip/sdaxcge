#!/usr/bin/env Rscript
# Extract per-cell TACSTD2 UMI counts and total UMI counts from the
# GSE205335 UMI matrix (dgCMatrix expected).

suppressPackageStartupMessages(library(Matrix))

path <- "/workspace/data/GSE205335_Lung_IO_UMI_matrix.rds"
cat("Reading RDS...\n")
obj <- readRDS(path)
cat("Class:", class(obj), "\n")

if (inherits(obj, "list")) {
  cat("List names:", names(obj), "\n")
  # take the first sparse matrix element
  for (el in obj) if (inherits(el, "sparseMatrix")) { obj <- el; break }
}

stopifnot(inherits(obj, "sparseMatrix") || is.matrix(obj))
cat("Dims:", dim(obj)[1], "genes x", dim(obj)[2], "cells\n")

genes <- rownames(obj)
hit <- grep("^TACSTD2$", genes)
cat("TACSTD2 row index:", hit, "\n")
if (length(hit) == 0) {
  cat("Exact TACSTD2 not found; candidates:\n")
  print(grep("TACSTD", genes, value = TRUE))
  quit(status = 1)
}

tacstd2 <- as.numeric(obj[hit, ])
total   <- Matrix::colSums(obj)
# EPCAM as a sanity-check epithelial marker
epcam_idx <- grep("^EPCAM$", genes)
epcam <- if (length(epcam_idx) == 1) as.numeric(obj[epcam_idx, ]) else NA
# PTPRC (CD45) as immune sanity check
ptprc_idx <- grep("^PTPRC$", genes)
ptprc <- if (length(ptprc_idx) == 1) as.numeric(obj[ptprc_idx, ]) else NA

df <- data.frame(
  barcode = colnames(obj),
  tacstd2_umi = tacstd2,
  epcam_umi = epcam,
  ptprc_umi = ptprc,
  total_umi = as.numeric(total)
)
out <- "/workspace/data/percell_tacstd2.csv.gz"
con <- gzfile(out, "w")
write.csv(df, con, row.names = FALSE)
close(con)
cat("Wrote", out, "with", nrow(df), "cells\n")
