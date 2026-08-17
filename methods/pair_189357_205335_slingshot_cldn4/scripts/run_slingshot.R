#!/usr/bin/env Rscript
# REAL Slingshot (Street et al. 2018) on Harmony/PCA + Leiden clusters.
# Start cluster is supplied by Python and is never the CLDN4-high cluster.

suppressPackageStartupMessages({
  if (nzchar(Sys.getenv("R_LIBS_USER"))) {
    .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  }
  library(slingshot)
})

args <- commandArgs(trailingOnly = TRUE)
kv <- list()
if (length(args) %% 2L != 0L) {
  stop("expected --key value pairs")
}
if (length(args)) {
  keys <- args[seq(1, length(args), by = 2)]
  vals <- args[seq(2, length(args), by = 2)]
  kv <- setNames(as.list(vals), sub("^--", "", keys))
}

emb_path <- kv$embedding
cl_path <- kv$clusters
start <- kv$start
out_prefix <- kv$out
if (is.null(emb_path) || is.null(cl_path) || is.null(start) || is.null(out_prefix)) {
  stop("usage: run_slingshot.R --embedding emb.csv --clusters cl.csv --start START --out prefix")
}

emb <- as.matrix(read.csv(emb_path, row.names = 1, check.names = FALSE))
cl_df <- read.csv(cl_path, row.names = 1, check.names = FALSE)
cl <- as.character(cl_df[[1]])
names(cl) <- rownames(cl_df)
if (!identical(rownames(emb), names(cl))) {
  common <- intersect(rownames(emb), names(cl))
  if (length(common) < 50) {
    stop("embedding and cluster barcodes do not align")
  }
  emb <- emb[common, , drop = FALSE]
  cl <- cl[common]
}

if (!(start %in% unique(cl))) {
  stop(sprintf("start cluster %s not in cluster labels", start))
}

sds <- slingshot(emb, clusterLabels = cl, start.clus = start, stretch = 2)
pt <- as.data.frame(slingPseudotime(sds))
pt$barcode <- rownames(pt)
w <- as.data.frame(slingCurveWeights(sds))
w$barcode <- rownames(w)
lineages <- slingLineages(sds)

lin_rows <- lapply(names(lineages), function(nm) {
  path <- lineages[[nm]]
  data.frame(
    lineage_id = nm,
    cluster_path = paste(path, collapse = "->"),
    n_clusters = length(path),
    start_cluster = path[[1]],
    end_cluster = path[[length(path)]],
    stringsAsFactors = FALSE
  )
})
lin_tab <- do.call(rbind, lin_rows)

write.csv(pt, paste0(out_prefix, "_pseudotime.csv"), row.names = FALSE)
write.csv(w, paste0(out_prefix, "_weights.csv"), row.names = FALSE)
write.csv(lin_tab, paste0(out_prefix, "_lineages.csv"), row.names = FALSE)
writeLines(sprintf("slingshot_ok lineages=%d cells=%d start=%s", nrow(lin_tab), nrow(pt), start))
