#!/usr/bin/env Rscript
# Real Slingshot (Street et al. 2018) on Harmony/PCA + Leiden.
# Start cluster is the AT2-rooted cluster, never a CLDN4-high cluster.

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args)) stop("missing value for ", args[[i]])
    out[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

opt <- parse_args(args)
need <- c("rd", "clusters", "start", "outdir")
miss <- setdiff(need, names(opt))
if (length(miss)) stop("missing args: ", paste(miss, collapse = ", "))

user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib)) {
  .libPaths(c(user_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(slingshot)
})

rd_path <- opt$rd
cl_path <- opt$clusters
start_clus <- as.character(opt$start)
outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rd <- as.matrix(read.delim(rd_path, header = TRUE, row.names = 1, check.names = FALSE))
cl <- read.delim(cl_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("cell_id", "cluster") %in% colnames(cl))) {
  stop("clusters TSV must have cell_id and cluster")
}
if (nrow(rd) != nrow(cl)) {
  stop(sprintf("rd rows %d != cluster rows %d", nrow(rd), nrow(cl)))
}
if (!all(rownames(rd) == cl$cell_id)) {
  idx <- match(rownames(rd), cl$cell_id)
  if (any(is.na(idx))) stop("cell_id mismatch between rd and clusters")
  cl <- cl[idx, , drop = FALSE]
}

clusters <- as.character(cl$cluster)
if (!(start_clus %in% clusters)) {
  stop(sprintf("start cluster %s not in cluster labels", start_clus))
}

cat(sprintf(
  "slingshot n_cells=%d n_pcs=%d n_clusters=%d start=%s\n",
  nrow(rd), ncol(rd), length(unique(clusters)), start_clus
), flush = TRUE)

sds <- slingshot(
  rd,
  clusterLabels = clusters,
  start.clus = start_clus,
  stretch = 0
)

lin <- slingLineages(sds)
pt <- as.data.frame(slingPseudotime(sds))
wt <- as.data.frame(slingCurveWeights(sds))
pt$cell_id <- rownames(rd)
wt$cell_id <- rownames(rd)
pt$cluster <- clusters
wt$cluster <- clusters

lineage_rows <- list()
for (i in seq_along(lin)) {
  path <- as.character(lin[[i]])
  col <- colnames(pt)[i]
  vals <- pt[[col]]
  assigned <- is.finite(vals)
  lineage_rows[[i]] <- data.frame(
    lineage_id = sprintf("Lineage%s", i),
    slingshot_column = col,
    start_cluster = path[1],
    end_cluster = path[length(path)],
    path = paste(path, collapse = "->"),
    n_clusters = length(path),
    n_cells_finite_pt = sum(assigned),
    stringsAsFactors = FALSE
  )
}
lineage_df <- do.call(rbind, lineage_rows)

write.table(
  lineage_df,
  file = file.path(outdir, "slingshot_lineages.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  pt[, c("cell_id", "cluster", setdiff(colnames(pt), c("cell_id", "cluster")))],
  file = file.path(outdir, "slingshot_cell_pseudotime.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  wt[, c("cell_id", "cluster", setdiff(colnames(wt), c("cell_id", "cluster")))],
  file = file.path(outdir, "slingshot_curve_weights.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

info <- list(
  n_cells = nrow(rd),
  n_lineages = length(lin),
  start_cluster = start_clus,
  slingshot_version = as.character(packageVersion("slingshot")),
  lineages = lapply(lin, as.character)
)
if (requireNamespace("jsonlite", quietly = TRUE)) {
  writeLines(
    jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE),
    file.path(outdir, "slingshot_info.json")
  )
} else {
  lin_txt <- paste(
    vapply(seq_along(info$lineages), function(i) {
      sprintf(
        "    \"Lineage%d\": [\"%s\"]",
        i,
        paste(info$lineages[[i]], collapse = "\", \"")
      )
    }, character(1)),
    collapse = ",\n"
  )
  writeLines(
    sprintf(
      "{\n  \"n_cells\": %d,\n  \"n_lineages\": %d,\n  \"start_cluster\": \"%s\",\n  \"slingshot_version\": \"%s\",\n  \"lineages\": {\n%s\n  }\n}\n",
      info$n_cells,
      info$n_lineages,
      info$start_cluster,
      info$slingshot_version,
      lin_txt
    ),
    file.path(outdir, "slingshot_info.json")
  )
}

cat(
  sprintf("wrote %d lineages to %s\n", nrow(lineage_df), outdir),
  flush = TRUE
)
