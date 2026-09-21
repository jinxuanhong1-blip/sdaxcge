#!/usr/bin/env Rscript
# Peel is done by the caller. Read the dgCMatrix RDS, keep MEBOCOST genes,
# and write a Matrix Market plus the full-transcriptome library size.
args <- commandArgs(trailingOnly = TRUE)
rds <- args[[1]]
genes_file <- args[[2]]
out <- args[[3]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
wanted <- unique(readLines(genes_file))
wanted <- wanted[nzchar(wanted)]
suppressPackageStartupMessages(library(Matrix))
message(format(Sys.time()), " readRDS ", rds)
mat <- readRDS(rds)
message(format(Sys.time()), " class=", paste(class(mat), collapse = ","),
        " dim=", nrow(mat), "x", ncol(mat))
if (is.null(rownames(mat)) || is.null(colnames(mat))) stop("matrix missing dimnames")
up <- toupper(rownames(mat))
keep <- up %in% toupper(wanted)
message("genes kept ", sum(keep), " / wanted ", length(wanted))
lib <- Matrix::colSums(mat)
sub <- mat[keep, , drop = FALSE]
rn <- toupper(rownames(sub))
if (any(duplicated(rn))) sub <- sub[!duplicated(rn), , drop = FALSE]
rownames(sub) <- toupper(rownames(sub))
message(format(Sys.time()), " write subset ", nrow(sub), "x", ncol(sub))
Matrix::writeMM(sub, file.path(out, "matrix.mtx"))
writeLines(rownames(sub), file.path(out, "genes.txt"))
writeLines(colnames(sub), file.path(out, "barcodes.txt"))
con <- gzfile(file.path(out, "lib.tsv.gz"), "wt")
writeLines("barcode\tlib", con)
bc <- colnames(mat)
n <- length(bc)
step <- 20000L
for (i in seq(1L, n, by = step)) {
  j <- min(n, i + step - 1L)
  writeLines(paste(bc[i:j], sprintf("%.6f", lib[i:j]), sep = "\t"), con)
}
close(con)
message(format(Sys.time()), " done nnz=", length(sub@x))
