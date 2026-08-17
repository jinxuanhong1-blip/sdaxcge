#!/usr/bin/env Rscript
# Real Slingshot (Street et al. 2018) on malignant PCA + Leiden.
# Root cluster is supplied and must not be the CLDN4-high cluster.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: run_slingshot.R <input_dir> <output_dir>")
}
indir <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(slingshot)
  library(SingleCellExperiment)
})

pca <- as.matrix(read.csv(file.path(indir, "pca.csv"), row.names = 1, check.names = FALSE))
umap <- as.matrix(read.csv(file.path(indir, "umap.csv"), row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(indir, "cell_meta.csv"), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
start_clus <- trimws(readLines(file.path(indir, "start_cluster.txt"), warn = FALSE)[[1]])

if (!all(rownames(pca) == rownames(meta))) {
  stop("pca / meta rownames differ")
}
clusters <- as.character(meta$leiden)
if (!(start_clus %in% clusters)) {
  stop(sprintf("start cluster %s not in Leiden labels", start_clus))
}

sce <- SingleCellExperiment(assays = list(logcounts = t(pca)))
reducedDims(sce) <- list(PCA = pca, UMAP = umap)
colData(sce)$leiden <- clusters

set.seed(1)
sds <- slingshot(
  sce,
  clusterLabels = "leiden",
  reducedDim = "PCA",
  start.clus = start_clus
)

pt <- slingPseudotime(sds)
avg <- slingAvgPseudotime(sds)
out_pt <- data.frame(
  barcode = rownames(pca),
  sling_avg_pseudotime = as.numeric(avg),
  stringsAsFactors = FALSE
)
for (j in seq_len(ncol(pt))) {
  out_pt[[colnames(pt)[j]]] <- as.numeric(pt[, j])
}
write.csv(out_pt, file.path(outdir, "slingshot_pseudotime.csv"), row.names = FALSE)

lineages <- slingLineages(sds)
lin_rows <- list()
for (nm in names(lineages)) {
  lin_rows[[length(lin_rows) + 1]] <- data.frame(
    lineage = nm,
    order = seq_along(lineages[[nm]]),
    leiden = as.character(lineages[[nm]]),
    stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, lin_rows), file.path(outdir, "slingshot_lineages.csv"), row.names = FALSE)

# Principal curves in UMAP for the lineage plot (geometry only).
sds_umap <- slingshot(
  sce,
  clusterLabels = "leiden",
  reducedDim = "UMAP",
  start.clus = start_clus
)
curves <- slingCurves(sds_umap)
curve_rows <- list()
for (nm in names(curves)) {
  crv <- curves[[nm]]$s
  if (is.null(crv) || !nrow(crv)) next
  curve_rows[[length(curve_rows) + 1]] <- data.frame(
    lineage = nm,
    step = seq_len(nrow(crv)),
    UMAP1 = crv[, 1],
    UMAP2 = crv[, 2],
    stringsAsFactors = FALSE
  )
}
if (length(curve_rows)) {
  write.csv(do.call(rbind, curve_rows), file.path(outdir, "slingshot_umap_curves.csv"), row.names = FALSE)
} else {
  write.csv(
    data.frame(lineage = character(), step = integer(), UMAP1 = numeric(), UMAP2 = numeric()),
    file.path(outdir, "slingshot_umap_curves.csv"),
    row.names = FALSE
  )
}

info <- list(
  start_cluster = start_clus,
  n_cells = nrow(pca),
  n_lineages = length(lineages),
  lineages = paste(vapply(names(lineages), function(nm) {
    paste0(nm, ":", paste(lineages[[nm]], collapse = ">"))
  }, character(1)), collapse = "; "),
  slingshot_version = as.character(packageVersion("slingshot"))
)
writeLines(
  sprintf(
    "start_cluster=%s\nn_cells=%s\nn_lineages=%s\nlineages=%s\nslingshot_version=%s\n",
    info$start_cluster, info$n_cells, info$n_lineages, info$lineages, info$slingshot_version
  ),
  file.path(outdir, "slingshot_info.txt")
)
cat("SLINGSHOT_DONE\n")
cat(info$lineages, "\n")
