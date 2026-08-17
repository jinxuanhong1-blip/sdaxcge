#!/usr/bin/env Rscript
# Official Bioconductor slingshot (Street et al. 2018) if installed.
# Reads TSV: cell_id, cluster, dim1..dimK  and a start cluster.
# Writes lineage pseudotime TSV.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: run_slingshot.R embedding.tsv start_cluster out.tsv")
}
emb_path <- args[[1]]
start_clus <- args[[2]]
out_path <- args[[3]]

if (!requireNamespace("slingshot", quietly = TRUE)) {
  stop("slingshot is not installed")
}

suppressPackageStartupMessages({
  library(slingshot)
})

tab <- read.delim(emb_path, check.names = FALSE, stringsAsFactors = FALSE)
rd_cols <- grep("^dim", colnames(tab), value = TRUE)
if (length(rd_cols) < 2) {
  stop("need dim1.. columns")
}
rd <- as.matrix(tab[, rd_cols, drop = FALSE])
storage.mode(rd) <- "double"
clus <- as.character(tab$cluster)
fit <- slingshot(rd, clusterLabels = clus, start.clus = start_clus)
pt <- slingPseudotime(fit)
w <- slingCurveWeights(fit)
out <- data.frame(cell_id = tab$cell_id, cluster = clus, pt, check.names = FALSE)
if (!is.null(w)) {
  colnames(w) <- paste0("weight_", colnames(w))
  out <- cbind(out, w)
}
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s lineages=%s cells=%s\n", out_path, ncol(pt), nrow(out)))
