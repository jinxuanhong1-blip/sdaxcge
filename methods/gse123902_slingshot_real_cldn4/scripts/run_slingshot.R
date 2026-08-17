#!/usr/bin/env Rscript
# Real Slingshot (Street et al. 2018) on PCA + Leiden.
# Writes lineage membership and per-cell slingPseudotime.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: run_slingshot.R <indir> <outdir>")
}
indir <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(slingshot)
})

rd <- as.matrix(read.csv(file.path(indir, "pca.csv"), row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(indir, "clusters.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(rd) == nrow(meta))
cl <- as.character(meta$leiden)
start <- as.character(readLines(file.path(indir, "start_cluster.txt"), warn = FALSE)[[1]])

sds <- slingshot(rd, clusterLabels = cl, start.clus = start, stretch = 2)
pt <- slingPseudotime(sds)
w <- slingCurveWeights(sds)
lin <- slingLineages(sds)

pt_df <- data.frame(cell = rownames(rd), cluster = cl, pt, check.names = FALSE)
w_df <- data.frame(cell = rownames(rd), w, check.names = FALSE)
write.csv(pt_df, file.path(outdir, "sling_pseudotime.csv"), row.names = FALSE)
write.csv(w_df, file.path(outdir, "sling_weights.csv"), row.names = FALSE)

rows <- list()
for (i in seq_along(lin)) {
  path <- lin[[i]]
  rows[[i]] <- data.frame(
    lineage_id = names(lin)[[i]],
    start_cluster = path[[1]],
    end_cluster = path[[length(path)]],
    path = paste(path, collapse = ">"),
    n_clusters = length(path),
    stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, rows), file.path(outdir, "sling_lineages.csv"), row.names = FALSE)
cat("ok lineages=", length(lin), " cells=", nrow(rd), "\n", sep = "")
