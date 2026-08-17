#!/usr/bin/env Rscript
# Real Bioconductor slingshot (Street et al. 2018) on Harmony/PCA + Leiden.
# Root cluster is supplied from Python and is never the CLDN4-high cluster.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("usage: run_slingshot.R --pca FILE --cluster FILE --start CLUSTER --outdir DIR")
}

kv <- list()
i <- 1
while (i <= length(args)) {
  key <- sub("^--", "", args[[i]])
  if (i == length(args)) stop("missing value for ", key)
  kv[[key]] <- args[[i + 1]]
  i <- i + 2
}
for (need in c("pca", "cluster", "start", "outdir")) {
  if (is.null(kv[[need]])) stop("missing --", need)
}

suppressPackageStartupMessages({
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("slingshot is not installed. Run scripts/install_r_slingshot.sh")
  }
  library(slingshot)
})

pca <- as.matrix(read.delim(kv$pca, header = TRUE, row.names = 1, check.names = FALSE))
cl <- read.delim(kv$cluster, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("cell_id", "cluster") %in% names(cl))) {
  stop("cluster table must have cell_id and cluster")
}
if (nrow(pca) != nrow(cl) || !all(rownames(pca) == cl$cell_id)) {
  stop("PCA rownames must match cluster$cell_id in the same order")
}

start <- as.character(kv$start)
clusters <- as.character(cl$cluster)
if (!start %in% clusters) {
  stop("start cluster ", start, " is not among Leiden labels")
}

message("slingshot n_cells=", nrow(pca), " n_clusters=", length(unique(clusters)),
        " start=", start)

sds <- slingshot(pca, clusterLabels = clusters, start.clus = start, stretch = 2)
lin <- slingLineages(sds)
pt <- as.data.frame(slingPseudotime(sds))
pt$cell_id <- rownames(pca)
wts <- as.data.frame(slingCurveWeights(sds))
wts$cell_id <- rownames(pca)

outdir <- kv$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

lineage_rows <- list()
for (name in names(lin)) {
  path <- as.character(lin[[name]])
  ptv <- pt[[name]]
  keep <- is.finite(ptv)
  lineage_rows[[name]] <- data.frame(
    lineage_id = name,
    start_cluster = path[[1]],
    end_cluster = path[[length(path)]],
    path = paste(path, collapse = "->"),
    n_clusters = length(path),
    n_cells_finite_pt = sum(keep),
    pt_min = if (any(keep)) min(ptv[keep]) else NA_real_,
    pt_max = if (any(keep)) max(ptv[keep]) else NA_real_,
    stringsAsFactors = FALSE
  )
}
lineage_df <- do.call(rbind, lineage_rows)
write.table(lineage_df, file.path(outdir, "slingshot_lineages.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pt, file.path(outdir, "slingshot_pseudotime.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(wts, file.path(outdir, "slingshot_weights.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

mst_ok <- FALSE
try({
  mst <- slingMST(sds)
  if (!is.null(mst) && requireNamespace("igraph", quietly = TRUE)) {
    ed <- igraph::as_data_frame(mst, what = "edges")
    write.table(ed, file.path(outdir, "slingshot_mst_edges.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    mst_ok <- TRUE
  }
}, silent = TRUE)
if (!mst_ok) message("MST edges not written (optional)")

info <- data.frame(
  package = "slingshot",
  version = as.character(packageVersion("slingshot")),
  start_cluster = start,
  n_lineages = length(lin),
  n_cells = nrow(pca),
  n_clusters = length(unique(clusters)),
  stringsAsFactors = FALSE
)
write.table(info, file.path(outdir, "slingshot_run_info.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", file.path(outdir, "slingshot_lineages.tsv"))
