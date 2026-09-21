#!/usr/bin/env Rscript
# Subset the GSE205335 dgCMatrix RDS to malignant barcodes x requested genes.
# The GEO file is gzip-wrapped twice. Python passes an already-inflated RDS.

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: subset_gse205335.R <rds> <barcode_patient.tsv> <genes.txt> <out_dir>")
}
rds <- args[[1]]
bc_file <- args[[2]]
gene_file <- args[[3]]
out_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wanted_bc <- read.delim(bc_file, stringsAsFactors = FALSE, sep = "\t")
stopifnot(all(c("barcode", "patient") %in% colnames(wanted_bc)))
genes <- unique(toupper(readLines(gene_file)))
genes <- genes[nzchar(genes)]

message("reading RDS")
mat <- readRDS(rds)
if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
message("full dim ", nrow(mat), " x ", ncol(mat))

colnames(mat) <- make.unique(colnames(mat))
keep <- wanted_bc$barcode %in% colnames(mat)
message("barcodes requested ", nrow(wanted_bc), " matched ", sum(keep))
if (sum(keep) < 1000) stop("too few GSE205335 barcodes matched")
wanted_bc <- wanted_bc[keep, , drop = FALSE]
# preserve requested order
wanted_bc <- wanted_bc[!duplicated(wanted_bc$barcode), , drop = FALSE]
sub <- mat[, wanted_bc$barcode, drop = FALSE]
rm(mat)
gc(verbose = FALSE)

n_count <- as.numeric(Matrix::colSums(sub))
n_feat <- as.numeric(Matrix::colSums(sub > 0))
mt <- grepl("^MT-", toupper(rownames(sub)))
n_mt <- if (any(mt)) as.numeric(Matrix::colSums(sub[mt, , drop = FALSE])) else rep(0, ncol(sub))

rn <- toupper(rownames(sub))
want <- rn %in% genes
sub_g <- sub[want, , drop = FALSE]
rn <- rn[want]
if (anyDuplicated(rn)) {
  dup <- unique(rn[duplicated(rn)])
  drop_rows <- integer(0)
  for (g in dup) {
    rows <- which(rn == g)
    summed <- Matrix::colSums(sub_g[rows, , drop = FALSE])
    sub_g[rows[[1]], ] <- summed
    drop_rows <- c(drop_rows, rows[-1])
  }
  sub_g <- sub_g[-drop_rows, , drop = FALSE]
  rn <- rn[-drop_rows]
}
message("gene subset ", nrow(sub_g), " x ", ncol(sub_g))

write.table(
  data.frame(
    barcode = wanted_bc$barcode,
    patient = wanted_bc$patient,
    n_count = n_count,
    n_feat = n_feat,
    n_mt = n_mt,
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "cell_qc.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(rn, file.path(out_dir, "genes.txt"))
Matrix::writeMM(sub_g, file.path(out_dir, "matrix.mtx"))
message("wrote ", out_dir)
