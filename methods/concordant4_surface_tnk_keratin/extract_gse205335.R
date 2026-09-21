#!/usr/bin/env Rscript
# Read the double-gzipped GSE205335 dgCMatrix and write the surface/keratin rows.
# R is required only because the public UMI file is an RDS S4 matrix.

args <- commandArgs(trailingOnly = TRUE)
geo <- if (length(args) >= 1) args[[1]] else "/tmp/geo_c4"
gz <- file.path(geo, "GSE205335_Lung_IO_UMI_matrix.rds.gz")
rds <- file.path(geo, "GSE205335_Lung_IO_UMI_matrix.rds")
out <- file.path(geo, "GSE205335_surface_genes.tsv.gz")
if (!file.exists(rds)) {
  if (!file.exists(gz)) stop("missing ", gz)
  message("double-gunzip ", gz)
  system2("bash", c("-c", sprintf("gzip -dc %s | gzip -dc > %s", shQuote(gz), shQuote(rds))))
}

suppressPackageStartupMessages(library(Matrix))
m <- readRDS(rds)
wanted <- c("CLDN1", "CLDN4", "CLDN7", "EPCAM", "KRT8", "KRT18", "KRT19", "MUC1")
rn <- toupper(rownames(m))
idx <- integer(0)
used <- character(0)
for (g in wanted) {
  hit <- which(rn == g)
  if (!length(hit)) stop("GSE205335 missing gene ", g)
  if (length(hit) > 1) message("duplicate ", g, " n=", length(hit), " using first row")
  idx <- c(idx, hit[[1]])
  used <- c(used, g)
}
sub <- as.matrix(m[idx, , drop = FALSE])
rownames(sub) <- used
message("export ", paste(dim(sub), collapse = " x "))
con <- gzfile(out, "wt")
writeLines(paste(c("barcode", used), collapse = "\t"), con)
bc <- colnames(m)
chunk <- 5000L
for (i in seq(1, length(bc), by = chunk)) {
  j <- min(i + chunk - 1L, length(bc))
  block <- t(sub[, i:j, drop = FALSE])
  lines <- apply(cbind(bc[i:j], format(block, trim = TRUE, scientific = FALSE)), 1, paste, collapse = "\t")
  writeLines(lines, con)
}
close(con)
message("wrote ", out)
