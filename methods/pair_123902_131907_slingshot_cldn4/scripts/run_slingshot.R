#!/usr/bin/env Rscript
# REAL Slingshot (Street et al. 2018) on Harmony/PCA + Leiden.
# Root cluster is supplied by Python and is never the CLDN4-high cluster.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("usage: run_slingshot.R --embed PATH --clusters PATH --start CLUSTER --outdir DIR")
}

kv <- list()
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[[i]])
  if (i == length(args)) stop(paste("missing value for", args[[i]]))
  kv[[key]] <- args[[i + 1]]
  i <- i + 2
}
for (need in c("embed", "clusters", "start", "outdir")) {
  if (is.null(kv[[need]])) stop(paste("missing --", need, sep = ""))
}

user_lib <- Sys.getenv("R_LIBS_USER", unset = file.path(Sys.getenv("HOME"), "R", "library"))
if (nzchar(user_lib) && dir.exists(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

suppressPackageStartupMessages({
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("slingshot is not installed")
  }
  library(slingshot)
})

embed_path <- kv[["embed"]]
cl_path <- kv[["clusters"]]
start_clus <- kv[["start"]]
outdir <- kv[["outdir"]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

embed <- as.matrix(read.delim(embed_path, row.names = 1, check.names = FALSE))
cl <- read.delim(cl_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("cell", "cluster") %in% colnames(cl))) {
  stop("clusters table must have columns cell, cluster")
}
if (!identical(rownames(embed), cl$cell)) {
  map <- match(rownames(embed), cl$cell)
  if (any(is.na(map))) stop("cluster table missing embedding cells")
  cl <- cl[map, , drop = FALSE]
}
clusters <- as.character(cl$cluster)
if (!(start_clus %in% clusters)) {
  stop(paste("start cluster", start_clus, "not present"))
}

# Use the first 10 Harmony PCs for the curve fit (Street 2018 used low-D
# embeddings). Neighbors / Leiden stay on 30 PCs in Python.
if (ncol(embed) > 10) {
  embed <- embed[, seq_len(10), drop = FALSE]
}

cat(sprintf(
  "slingshot n_cells=%d n_dim=%d n_clusters=%d start=%s version=%s\n",
  nrow(embed), ncol(embed), length(unique(clusters)), start_clus,
  as.character(packageVersion("slingshot"))
))
flush.console()

cat("getLineages...\n")
flush.console()
sds <- getLineages(
  embed,
  clusterLabels = clusters,
  start.clus = start_clus
)
cat(sprintf("lineages=%s\n", paste(names(slingLineages(sds)), collapse = ",")))
flush.console()

cat("getCurves (approx_points=150)...\n")
flush.console()
# approx_points is the official large-n option in slingshot >=2.0
# (Street 2018 algorithm; not a DPT fallback).
sds <- getCurves(
  sds,
  stretch = 2,
  approx_points = 150
)
cat("curves done\n")
flush.console()

pt <- as.data.frame(slingPseudotime(sds))
pt$cell <- rownames(embed)
w <- as.data.frame(slingCurveWeights(sds))
w$cell <- rownames(embed)

lineages <- slingLineages(sds)
lin_rows <- list()
for (nm in names(lineages)) {
  path <- lineages[[nm]]
  col <- nm
  if (!(col %in% colnames(pt))) {
    # slingshot sometimes names columns Lineage1 / curve1
    alt <- grep(paste0("^", nm, "$|^curve", sub("Lineage", "", nm), "$"), colnames(pt), value = TRUE)
    col <- if (length(alt)) alt[[1]] else colnames(pt)[1]
  }
  vals <- pt[[col]]
  lin_rows[[length(lin_rows) + 1]] <- data.frame(
    lineage = nm,
    n_clusters = length(path),
    path = paste(path, collapse = ">"),
    start_cluster = path[[1]],
    end_cluster = path[[length(path)]],
    n_cells_on_lineage = sum(is.finite(vals)),
    pt_column = col,
    stringsAsFactors = FALSE
  )
}
lin_df <- do.call(rbind, lin_rows)

write.table(pt, file.path(outdir, "slingshot_pseudotime.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(w, file.path(outdir, "slingshot_weights.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(lin_df, file.path(outdir, "slingshot_lineages.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

info <- list(
  n_cells = nrow(embed),
  n_lineages = length(lineages),
  start_cluster = start_clus,
  version = as.character(packageVersion("slingshot")),
  lineages = names(lineages)
)
writeLines(
  sprintf(
    "n_cells\t%d\nn_lineages\t%d\nstart_cluster\t%s\nversion\t%s\nlineages\t%s\n",
    info$n_cells, info$n_lineages, info$start_cluster, info$version,
    paste(info$lineages, collapse = ",")
  ),
  file.path(outdir, "slingshot_info.tsv")
)
cat(sprintf("wrote %d lineages to %s\n", info$n_lineages, outdir))
