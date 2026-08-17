#!/usr/bin/env Rscript
# Real Slingshot (Street et al. 2018) on a reduced-dimension matrix.
# Args: embedding.csv clusters.csv start_cluster out_pseudotime.csv
# embedding.csv: cells x dims, row names = cell ids, no header required beyond colnames.
# clusters.csv: two columns cell,cluster (header).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: run_slingshot.R embedding.csv clusters.csv start_cluster out.csv")
}
emb_path <- args[[1]]
clu_path <- args[[2]]
start_cluster <- args[[3]]
out_path <- args[[4]]

suppressPackageStartupMessages({
  ok <- requireNamespace("slingshot", quietly = TRUE)
  if (!ok) stop("slingshot is not installed")
  library(slingshot)
})

emb <- as.matrix(read.csv(emb_path, row.names = 1, check.names = FALSE))
clu <- read.csv(clu_path, stringsAsFactors = FALSE)
if (!all(c("cell", "cluster") %in% colnames(clu))) {
  stop("clusters.csv must have columns cell,cluster")
}
rownames(clu) <- clu$cell
common <- intersect(rownames(emb), rownames(clu))
if (length(common) < 20) {
  stop(sprintf("too few shared cells: %d", length(common)))
}
emb <- emb[common, , drop = FALSE]
clusters <- as.character(clu[common, "cluster"])
if (!(start_cluster %in% clusters)) {
  # fall back to the most frequent cluster among the first 50 rows tagged start-like
  stop(sprintf("start cluster %s not in cluster vector", start_cluster))
}

sds <- slingshot::slingshot(
  data = emb,
  clusterLabels = clusters,
  start.clus = start_cluster
)
pt <- slingshot::slingPseudotime(sds)
avg <- rowMeans(pt, na.rm = TRUE)
out <- data.frame(
  cell = rownames(emb),
  slingshot_pseudotime = as.numeric(avg),
  n_lineages = rowSums(is.finite(pt)),
  stringsAsFactors = FALSE
)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_path, row.names = FALSE)
cat(sprintf("wrote %s n=%d start=%s\n", out_path, nrow(out), start_cluster))
