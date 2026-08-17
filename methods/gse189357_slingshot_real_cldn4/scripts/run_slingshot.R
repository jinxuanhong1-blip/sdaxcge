#!/usr/bin/env Rscript
# REAL Slingshot (Street et al. 2018) on GSE189357 tumor epithelium.
# Root / start cluster is supplied by Python and is never the CLDN4-high cluster.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("usage: run_slingshot.R <workdir>")
}
wd <- args[[1]]
r_libs <- Sys.getenv("R_LIBS_USER", "/tmp/r-libs")
dir.create(r_libs, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(r_libs, .libPaths()))

suppressPackageStartupMessages({
  library(slingshot)
})

rd <- as.matrix(read.csv(file.path(wd, "reduced_dim.csv"), row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(wd, "clusters.csv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(rownames(rd) == rownames(meta))) {
  stop("row names of reduced_dim.csv and clusters.csv do not match")
}
cl <- factor(as.character(meta$leiden))
start <- readLines(file.path(wd, "start_cluster.txt"), warn = FALSE)[[1]]
start <- gsub("^\\s+|\\s+$", "", start)
if (!start %in% levels(cl)) {
  stop(sprintf("start cluster %s not in Leiden labels", start))
}

set.seed(0)
sds <- slingshot(rd, clusterLabels = cl, start.clus = start, stretch = 2)

pt <- as.data.frame(slingPseudotime(sds))
pt$cell <- rownames(rd)
weights <- as.data.frame(slingCurveWeights(sds))
weights$cell <- rownames(rd)

lineages <- slingLineages(sds)
lin_rows <- list()
for (i in seq_along(lineages)) {
  path <- lineages[[i]]
  col <- colnames(slingPseudotime(sds))[i]
  vals <- slingPseudotime(sds)[, i]
  keep <- is.finite(vals)
  lin_rows[[i]] <- data.frame(
    lineage = col,
    lineage_index = i,
    clusters = paste(path, collapse = "->"),
    n_clusters = length(path),
    start_cluster = path[[1]],
    end_cluster = path[[length(path)]],
    n_cells = sum(keep),
    mean_pseudotime = mean(vals[keep]),
    stringsAsFactors = FALSE
  )
}
lin_df <- do.call(rbind, lin_rows)

mst <- tryCatch(slingMST(sds), error = function(e) NULL)
if (!is.null(mst) && requireNamespace("igraph", quietly = TRUE) && inherits(mst, "igraph")) {
  mst_mat <- as.matrix(igraph::as_adjacency_matrix(mst, sparse = FALSE))
  write.csv(mst_mat, file.path(wd, "slingshot_mst.csv"), quote = TRUE)
}

write.csv(pt, file.path(wd, "slingshot_pseudotime.csv"), row.names = FALSE)
write.csv(weights, file.path(wd, "slingshot_weights.csv"), row.names = FALSE)
write.csv(lin_df, file.path(wd, "slingshot_lineages.csv"), row.names = FALSE)
writeLines(as.character(packageVersion("slingshot")), file.path(wd, "slingshot_version.txt"))
cat("wrote", nrow(lin_df), "lineages; start=", start, "\n")
