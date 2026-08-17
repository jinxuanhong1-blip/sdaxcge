#!/usr/bin/env Rscript
# Real Slingshot (Street et al. 2018) on a reduced-dimension matrix.
# For n_cells >> 10k, lineages are fit on a stratified subsample
# (max 80 cells / cluster, approx_points=150) and written with the
# subsample embedding so Python can kNN-project the rest.
# Args: embedding.csv clusters.csv start_cluster out.csv

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
  stop(sprintf("start cluster %s not in cluster vector", start_cluster))
}

set.seed(1)
max_per <- 80
keep <- unlist(lapply(split(seq_along(clusters), clusters), function(ix) {
  if (length(ix) <= max_per) ix else sample(ix, max_per)
}), use.names = FALSE)
cat(sprintf("slingshot subsample %d / %d cells; start=%s\n", length(keep), nrow(emb), start_cluster))

sds <- slingshot::slingshot(
  data = emb[keep, , drop = FALSE],
  clusterLabels = clusters[keep],
  start.clus = start_cluster,
  approx_points = 150
)
pt <- slingshot::slingPseudotime(sds)
avg <- rowMeans(pt, na.rm = TRUE)
out <- data.frame(
  cell = rownames(emb)[keep],
  slingshot_pseudotime = as.numeric(avg),
  n_lineages = rowSums(is.finite(pt)),
  stringsAsFactors = FALSE
)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(out, out_path, row.names = FALSE)
# subsample embedding for kNN projection
write.csv(
  as.data.frame(emb[keep, , drop = FALSE]),
  file.path(dirname(out_path), "subsample_embedding.csv")
)
cat(sprintf("wrote %s n=%d start=%s\n", out_path, nrow(out), start_cluster))
