#!/usr/bin/env Rscript
# GSE205335 is a double-gzipped dgCMatrix. Keep one copy in memory and
# write malignant summaries for analyze.py.

suppressPackageStartupMessages(library(Matrix))
options(warn = 1)

say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

rds_path <- Sys.getenv("RDS")
map_path <- Sys.getenv("MAP")
out_dir <- Sys.getenv("OUT")
wanted_path <- Sys.getenv("WANTED")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!nzchar(rds_path) || !nzchar(map_path) || !nzchar(out_dir) || !nzchar(wanted_path)) {
  stop("Set RDS, MAP, OUT, and WANTED")
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
cn <- colnames(mat)
say("matrix", nrow(mat), "genes", ncol(mat), "cells; example", cn[1])

map <- read.delim(map_path, stringsAsFactors = FALSE)
map$barcode <- as.character(map$barcode)
stopifnot(all(c("barcode", "unit_id", "is_mal", "is_tnk") %in% names(map)))

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
  stop("barcode mismatch. ident example ", paste(head(map$barcode, 3), collapse = " | "))
}
labels <- cands[[best]]
m <- match(labels, map$barcode)
ok <- which(!is.na(m))
say("using", best, "matched", length(ok))
submap <- map[m[ok], , drop = FALSE]
lib <- as.numeric(Matrix::colSums(mat)[ok])
unit <- submap$unit_id
is_mal <- submap$is_mal %in% c(TRUE, "TRUE", "True", 1, "1")
is_tnk <- submap$is_tnk %in% c(TRUE, "TRUE", "True", 1, "1")

units <- sort(unique(unit))
count_tab <- do.call(rbind, lapply(units, function(u) {
  ix <- unit == u
  data.frame(
    unit_id = u,
    n_cells = sum(ix),
    n_malignant = sum(ix & is_mal),
    n_tnk = sum(ix & is_tnk),
    stringsAsFactors = FALSE
  )
}))
write.table(count_tab, file.path(out_dir, "counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
say("counts", nrow(count_tab), "units", sum(count_tab$n_malignant), "malignant")

mal_cols <- ok[is_mal]
mal_unit <- unit[is_mal]
pats <- count_tab$unit_id[count_tab$n_malignant > 0]

sum_log1p_cols <- function(cols) {
  piece <- mat[, cols, drop = FALSE]
  out <- numeric(nrow(mat))
  if (length(piece@x)) {
    rs <- rowsum(log1p(piece@x), group = piece@i + 1L, reorder = FALSE)
    out[as.integer(rownames(rs))] <- rs[, 1]
  }
  out / length(cols)
}

say("mean log1p over", length(pats), "patients")
sum_log <- matrix(0, nrow = nrow(mat), ncol = length(pats),
                  dimnames = list(rownames(mat), pats))
for (i in seq_along(pats)) {
  cols <- mal_cols[mal_unit == pats[i]]
  sum_log[, i] <- sum_log1p_cols(cols)
  if (i %% 5 == 0) say("  patients", i)
}
keep_rows <- rowSums(sum_log) > 0
sum_log <- sum_log[keep_rows, , drop = FALSE]
write.table(sum_log, file.path(out_dir, "all_mean_log1p.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
say("all_mean genes", nrow(sum_log))
rm(sum_log)

focus <- intersect(c("ELF3", "GRHL1", "GRHL2", "GRHL3", "TACSTD2", "CLDN4", "CLDN7"), rownames(mat))
dense <- as.matrix(mat[focus, mal_cols, drop = FALSE])
cell <- data.frame(unit_id = mal_unit, t(dense), check.names = FALSE)
write.table(cell, file.path(out_dir, "cell_focus.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
say("cell_focus", nrow(cell), "genes", paste(focus, collapse = ","))
rm(cell, dense)

wanted <- unique(intersect(toupper(readLines(wanted_path)), rownames(mat)))
say("cp10k genes", length(wanted))
cp_sum <- matrix(0, nrow = length(wanted), ncol = length(pats),
                 dimnames = list(wanted, pats))
for (i in seq_along(pats)) {
  cols <- mal_cols[mal_unit == pats[i]]
  piece <- mat[wanted, cols, drop = FALSE]
  libs <- pmax(lib[is_mal][mal_unit == pats[i]], 1)
  if (length(piece@x)) {
    j <- rep.int(seq_len(ncol(piece)), diff(piece@p))
    vals <- log1p(piece@x / libs[j] * 1e4)
    rs <- rowsum(vals, group = piece@i + 1L, reorder = FALSE)
    cp_sum[as.integer(rownames(rs)), i] <- rs[, 1]
  }
  cp_sum[, i] <- cp_sum[, i] / length(cols)
}
write.table(cp_sum, file.path(out_dir, "target_cp10k.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
say("DONE")
