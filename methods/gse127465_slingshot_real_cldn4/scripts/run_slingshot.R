#!/usr/bin/env Rscript
# REAL Slingshot (Street et al. 2018) on precomputed PCA + Leiden.
# Start cluster is supplied by Python and must not be the CLDN4-high cluster.

suppressPackageStartupMessages({
  if (nzchar(Sys.getenv("R_LIBS_USER"))) {
    .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  }
  library(slingshot)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  args[[i + 1]]
}

pca_path <- get_arg("--pca")
cl_path <- get_arg("--clusters")
start <- get_arg("--start")
out_dir <- get_arg("--out")
if (is.null(pca_path) || is.null(cl_path) || is.null(start) || is.null(out_dir)) {
  stop("Usage: run_slingshot.R --pca PCA.csv --clusters clusters.csv --start START --out DIR")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pca <- as.matrix(read.csv(pca_path, row.names = 1, check.names = FALSE))
cl_df <- read.csv(cl_path, row.names = 1, check.names = FALSE)
if (!all(rownames(pca) %in% rownames(cl_df))) {
  stop("cluster rownames do not cover PCA rownames")
}
cl <- as.character(cl_df[rownames(pca), 1])
names(cl) <- rownames(pca)

if (!(start %in% unique(cl))) {
  stop(sprintf("start cluster %s not in cluster labels: %s",
               start, paste(sort(unique(cl)), collapse = ",")))
}

sds <- slingshot(
  pca,
  clusterLabels = cl,
  start.clus = start,
  stretch = 0
)

pt <- as.data.frame(slingPseudotime(sds))
pt$cell <- rownames(pca)
write.csv(pt, file.path(out_dir, "slingshot_pseudotime.csv"), row.names = FALSE)

w <- tryCatch(as.data.frame(slingCurveWeights(sds)), error = function(e) NULL)
if (!is.null(w)) {
  w$cell <- rownames(pca)
  write.csv(w, file.path(out_dir, "slingshot_weights.csv"), row.names = FALSE)
}

lins <- slingLineages(sds)
lin_rows <- lapply(names(lins), function(nm) {
  path <- as.character(lins[[nm]])
  data.frame(
    lineage_id = nm,
    start_cluster = path[[1]],
    end_cluster = path[[length(path)]],
    n_clusters = length(path),
    cluster_path = paste(path, collapse = "->"),
    stringsAsFactors = FALSE
  )
})
lin_df <- do.call(rbind, lin_rows)
write.csv(lin_df, file.path(out_dir, "slingshot_lineages_raw.csv"), row.names = FALSE)

sink(file.path(out_dir, "slingshot_session.txt"))
cat("slingshot", as.character(packageVersion("slingshot")), "\n")
cat("start.clus", start, "\n")
print(lins)
sink()
cat("OK lineages=", nrow(lin_df), " start=", start, "\n", sep = "")
