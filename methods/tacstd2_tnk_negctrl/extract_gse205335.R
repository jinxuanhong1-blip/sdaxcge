#!/usr/bin/env Rscript
# Read the double-gzipped GSE205335 dgCMatrix and write selected gene rows.
# R is required only because the public UMI file is an RDS S4 matrix.
# TACSTD2 is exported as a score. It is not a malignant gate.

args <- commandArgs(trailingOnly = TRUE)
geo <- if (length(args) >= 1) args[[1]] else "/tmp/geo_c4"
gz <- file.path(geo, "GSE205335_Lung_IO_UMI_matrix.rds.gz")
rds <- file.path(geo, "GSE205335_Lung_IO_UMI_matrix.rds")
out <- file.path(geo, "GSE205335_negctrl_genes.tsv.gz")
if (!file.exists(rds) || file.info(rds)$size < 1000) {
  if (file.exists(rds)) unlink(rds)
  if (!file.exists(gz)) stop("missing ", gz)
  message("gunzip ", gz)
  # system() (one string) is required so the shell sees the redirect.
  # The public file is a single gzip. If the result is still gzip, peel once more.
  status <- system(sprintf("gzip -dc %s > %s", shQuote(gz), shQuote(rds)))
  if (!identical(status, 0L)) stop("gunzip failed with status ", status)
  magic <- readBin(rds, "raw", 2L)
  if (length(magic) == 2L && identical(as.integer(magic), c(0x1fL, 0x8bL))) {
    message("inner gzip; peeling once more")
    tmp <- paste0(rds, ".inner")
    status <- system(sprintf("gzip -dc %s > %s", shQuote(rds), shQuote(tmp)))
    if (!identical(status, 0L)) stop("inner gunzip failed with status ", status)
    file.rename(tmp, rds)
  }
}

suppressPackageStartupMessages(library(Matrix))
message("readRDS")
m <- readRDS(rds)
if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
wanted <- c(
  "CLDN3", "CLDN4", "CLDN7", "EPCAM", "MUC1", "KRT19", "TACSTD2",
  "KRT8", "KRT18"
)
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
rm(m)
gc()
message("export ", paste(dim(sub), collapse = " x "))
con <- gzfile(out, "wt")
writeLines(paste(c("barcode", used), collapse = "\t"), con)
bc <- colnames(sub)
chunk <- 5000L
for (i in seq(1, length(bc), by = chunk)) {
  j <- min(i + chunk - 1L, length(bc))
  block <- t(sub[, i:j, drop = FALSE])
  lines <- apply(
    cbind(bc[i:j], format(block, trim = TRUE, scientific = FALSE)),
    1, paste, collapse = "\t"
  )
  writeLines(lines, con)
  if (i %% 20000L == 1L) message("  wrote through ", j)
}
close(con)
message("wrote ", out)
