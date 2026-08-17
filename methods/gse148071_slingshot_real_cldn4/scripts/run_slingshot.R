#!/usr/bin/env Rscript
# REAL Slingshot (Street et al. 2018) on exported PCA + Leiden.
# Root cluster is chosen in Python (Alveolar / AT2; never CLDN4-high).
suppressPackageStartupMessages({
  library(slingshot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: run_slingshot.R <work_dir>")
}
work <- args[[1]]
pca_path <- file.path(work, "sling_pca.tsv")
clu_path <- file.path(work, "sling_clusters.tsv")
start_path <- file.path(work, "sling_start_cluster.txt")
out_pt <- file.path(work, "sling_pseudotime.tsv")
out_w <- file.path(work, "sling_weights.tsv")
out_lin <- file.path(work, "sling_lineages.json")
out_log <- file.path(work, "sling_run.log")

sink(out_log, split = TRUE)
on.exit(sink(), add = TRUE)

cat("slingshot", as.character(packageVersion("slingshot")), "\n")
cat("work", work, "\n")

rd <- as.matrix(read.delim(pca_path, row.names = 1, check.names = FALSE))
clu_df <- read.delim(clu_path, row.names = 1, check.names = FALSE)
if (!all(rownames(rd) == rownames(clu_df))) {
  stop("PCA and cluster row names differ")
}
cl <- as.character(clu_df$leiden)
names(cl) <- rownames(clu_df)
start <- trimws(readLines(start_path, warn = FALSE)[[1]])
cat("n_cells", nrow(rd), "n_pcs", ncol(rd), "n_clusters", length(unique(cl)), "start", start, "\n")
if (!start %in% unique(cl)) {
  stop(sprintf("start cluster %s not in Leiden labels", start))
}

set.seed(0)
sds <- slingshot(
  rd,
  clusterLabels = cl,
  start.clus = start,
  stretch = 0
)
pt <- as.data.frame(slingPseudotime(sds))
wt <- as.data.frame(slingCurveWeights(sds))
pt$cell <- rownames(rd)
wt$cell <- rownames(rd)
# put cell first
pt <- pt[, c("cell", setdiff(colnames(pt), "cell"))]
wt <- wt[, c("cell", setdiff(colnames(wt), "cell"))]

lineages <- slingLineages(sds)
lin_list <- lapply(names(lineages), function(nm) {
  list(name = nm, clusters = as.character(lineages[[nm]]))
})
n_on <- vapply(seq_len(ncol(pt) - 1), function(i) sum(is.finite(pt[[i + 1]])), integer(1))
cat("lineages", length(lineages), "cells_on", paste(n_on, collapse = ","), "\n")

write.table(pt, out_pt, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wt, out_w, sep = "\t", quote = FALSE, row.names = FALSE)

# write lineages as JSON without extra packages
esc <- function(x) gsub("\"", "\\\"", x)
parts <- vapply(lin_list, function(L) {
  clus <- paste(sprintf("\"%s\"", esc(L$clusters)), collapse = ",")
  sprintf("{\"name\":\"%s\",\"clusters\":[%s]}", esc(L$name), clus)
}, character(1))
writeLines(sprintf("[%s]", paste(parts, collapse = ",")), out_lin)
cat("wrote", out_pt, out_w, out_lin, "\n")
cat("SLINGSHOT_DONE\n")
