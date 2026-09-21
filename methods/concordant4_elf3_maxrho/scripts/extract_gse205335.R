#!/usr/bin/env Rscript
# Malignant-cell panel counts for the locked GSE205335 patients.
# Identity map is built in Python; this step only subsets the dgCMatrix.

suppressPackageStartupMessages(library(Matrix))
options(warn = 1)

say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

rds_path <- Sys.getenv("RDS")
map_path <- Sys.getenv("MAP")
out_path <- Sys.getenv("OUT")
genes_path <- Sys.getenv("GENES")
if (!nzchar(rds_path) || !nzchar(map_path) || !nzchar(out_path) || !nzchar(genes_path)) {
  stop("Set RDS, MAP, OUT, and GENES")
}

say("readRDS", rds_path)
mat <- readRDS(rds_path)
if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
rn <- toupper(rownames(mat))
if (any(duplicated(rn) | !nzchar(rn))) {
  keep_gene <- !duplicated(rn) & nzchar(rn)
  mat <- mat[keep_gene, , drop = FALSE]
  rn <- rn[keep_gene]
}
rownames(mat) <- rn
if (!"CLDN4" %in% rn) stop("CLDN4 missing from rownames; example ", paste(head(rn, 6), collapse = ","))
cn <- colnames(mat)
say("matrix", nrow(mat), "genes", ncol(mat), "cells; example", cn[1])

map <- read.delim(map_path, stringsAsFactors = FALSE)
map$barcode <- as.character(map$barcode)
stopifnot(all(c("barcode", "unit_id", "is_mal", "is_tnk", "is_normal") %in% names(map)))

n_hit <- function(labels) sum(labels %in% map$barcode)
cands <- list(
  exact = cn,
  dash_to_dot = chartr("-", ".", cn),
  us_to_dash = chartr("_", "-", cn),
  dash_to_us = chartr("-", "_", cn),
  tail = sub("^.*_", "", cn)
)
hits <- vapply(cands, n_hit, numeric(1))
say("barcode overlaps", paste(names(hits), hits, sep = "=", collapse = " "))
best <- names(which.max(hits))
if (hits[[best]] < 1000) {
  stop("barcode mismatch. ident example ", paste(head(map$barcode, 3), collapse = " | "),
       " matrix example ", cn[1])
}
labels <- cands[[best]]
m <- match(labels, map$barcode)
ok <- which(!is.na(m))
say("using", best, "matched", length(ok))
submap <- map[m[ok], , drop = FALSE]
lib_all <- as.numeric(Matrix::colSums(mat)[ok])
unit <- submap$unit_id
is_mal <- submap$is_mal %in% c(TRUE, "TRUE", "True", 1, "1")
is_normal <- submap$is_normal %in% c(TRUE, "TRUE", "True", 1, "1")
tumor <- !is_normal
mal <- tumor & is_mal

wanted <- unique(toupper(readLines(genes_path)))
wanted <- wanted[nzchar(wanted)]
present <- intersect(wanted, rownames(mat))
missing <- setdiff(wanted, rownames(mat))
say("panel present", length(present), "missing", paste(missing, collapse = ","))

mal_cols <- ok[mal]
piece <- mat[present, mal_cols, drop = FALSE]
# dense panel is small (genes x malignant cells)
dense <- as.matrix(piece)
lib <- lib_all[mal]
unit_mal <- unit[mal]
say("malignant cells", length(mal_cols), "panel", nrow(dense))

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
con <- gzfile(out_path, "wt")
writeLines(paste(c("unit_id", "lib", present), collapse = "\t"), con)
# write in column chunks to limit the character buffer
block <- 2000L
n <- ncol(dense)
for (start in seq.int(1L, n, by = block)) {
  end <- min(n, start + block - 1L)
  sl <- start:end
  chunk <- data.frame(unit_id = unit_mal[sl], lib = lib[sl], t(dense[, sl, drop = FALSE]), check.names = FALSE)
  write.table(chunk, con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  if (end %% 10000L < block) say("  wrote", end)
}
close(con)
say("DONE", out_path)
