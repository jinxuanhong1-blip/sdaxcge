#!/usr/bin/env Rscript
# Slice the public GSE205335 dgCMatrix for CNV cells + malignant CLDN4.
# The GEO file is double-gzipped. Identity parsing is done in Python.

args <- commandArgs(trailingOnly = TRUE)
raw <- if (length(args) >= 1) args[[1]] else "/tmp/concordant4_raw"
cache <- if (length(args) >= 2) args[[2]] else "/tmp/concordant4_cnv_cache/v1/GSE205335_mats"
flags_path <- if (length(args) >= 3) args[[3]] else "/tmp/concordant4_cnv_cache/v1/gse205335_flags.tsv"

suppressPackageStartupMessages(library(Matrix))
dir.create(cache, recursive = TRUE, showWarnings = FALSE)

rds <- file.path(raw, "GSE205335_umi.rds")
if (!file.exists(rds) || file.info(rds)$size < 1e8) {
  src <- file.path(raw, "GSE205335_umi.rds.gz")
  message("double-gunzip ", src)
  system(sprintf("gzip -dc %s | gzip -dc > %s", shQuote(src), shQuote(rds)))
}

message("readRDS")
mat <- readRDS(rds)
if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
message("dim ", nrow(mat), " x ", ncol(mat))
flags <- read.delim(flags_path, stringsAsFactors = FALSE, check.names = FALSE)
message("flags ", nrow(flags), " colnames example ", colnames(mat)[1])

# Identity uses the same barcode string as the matrix (orig.ident_cell-1).
in_mat <- flags$barcode %in% colnames(mat)
message("barcode overlap ", sum(in_mat), " / ", nrow(flags))
if (sum(in_mat) < 1000) stop("GSE205335 barcode overlap too small")
flags <- flags[in_mat, , drop = FALSE]

mal <- flags[flags$mal == 1, , drop = FALSE]
cl <- rep(NA_real_, nrow(mal))
if (nrow(mal) && "CLDN4" %in% rownames(mat)) {
  cl <- as.numeric(mat["CLDN4", mal$barcode])
}
utils::write.table(
  data.frame(barcode = mal$barcode, patient = mal$patient, cldn4 = cl, stringsAsFactors = FALSE),
  file.path(cache, "cldn4_mal.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cnv <- flags[flags$epi == 1 | flags$ref == 1, , drop = FALSE]
for (pat in unique(cnv$patient)) {
  sub <- cnv[cnv$patient == pat, , drop = FALSE]
  # Stable order: epithelial rows first, then reference, as written by Python.
  m <- mat[, sub$barcode, drop = FALSE]
  keep <- Matrix::rowSums(m) > 0
  m <- m[keep, , drop = FALSE]
  Matrix::writeMM(m, file.path(cache, paste0(pat, ".mtx")))
  writeLines(rownames(m), file.path(cache, paste0(pat, ".genes")))
  writeLines(sub$barcode, file.path(cache, paste0(pat, ".barcodes")))
  utils::write.table(
    sub[, c("barcode", "patient", "epi", "ref", "mal")],
    file.path(cache, paste0(pat, ".flags.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  message("wrote ", pat, " genes ", nrow(m), " cells ", ncol(m))
}
message("export done")
