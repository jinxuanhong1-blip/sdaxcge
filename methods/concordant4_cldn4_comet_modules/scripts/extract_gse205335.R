#!/usr/bin/env Rscript
# Export author-malignant, non-normal cells from the GEO RDS (often double-gzipped).
# Writes a MatrixMarket counts matrix plus meta. Library size is the full-gene colSum.

args <- commandArgs(trailingOnly = TRUE)
RAW <- if (length(args) >= 1) args[[1]] else "/tmp/concordant4_raw"
OUT <- if (length(args) >= 2) args[[2]] else "/tmp/concordant4_malig/GSE205335"
HERE <- Sys.getenv("COMET_HERE", unset = "")
if (!nzchar(HERE)) {
  HERE <- normalizePath(file.path(getwd(), "methods/concordant4_cldn4_comet_modules"), mustWork = FALSE)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(Matrix))

logmsg <- function(...) {
  cat(paste0("[gse205335] ", paste(..., collapse = " "), "\n"))
  flush.console()
}

is_gzip <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 2)
  length(magic) == 2 && magic[[1]] == as.raw(0x1f) && magic[[2]] == as.raw(0x8b)
}

read_geo_rds <- function(path) {
  cur <- path
  temps <- character()
  on.exit(unlink(temps[file.exists(temps)]), add = TRUE)
  for (i in seq_len(4)) {
    if (!is_gzip(cur)) break
    dest <- tempfile(pattern = paste0("geo_rds_", i, "_"), tmpdir = tempdir())
    temps <- c(temps, dest)
    logmsg("gzip -dc layer", i)
    st <- system2("gzip", c("-dc", cur), stdout = dest)
    if (!identical(st, 0L)) stop("gzip -dc failed status=", st)
    cur <- dest
  }
  logmsg("readRDS bytes", file.info(cur)$size)
  readRDS(cur)
}

ident_path <- file.path(RAW, "GSE205335", "GSE205335_Lung_IO_CellIdentity.txt.gz")
rds_path <- file.path(RAW, "GSE205335", "GSE205335_Lung_IO_UMI_matrix.rds.gz")
map_path <- file.path(HERE, "data", "GSE205335_gsm_map.tsv")
locked_path <- file.path(HERE, "data", "locked_units.tsv")

ident <- read.delim(ident_path, stringsAsFactors = FALSE, check.names = FALSE)
gsm <- read.delim(map_path, stringsAsFactors = FALSE, check.names = FALSE)
locked <- read.delim(locked_path, stringsAsFactors = FALSE, check.names = FALSE)
locked_pt <- locked$unit_id[locked$dataset == "GSE205335"]

ident <- merge(ident, gsm[, c("orig.ident", "patient", "tissue")], by = "orig.ident",
               all.x = TRUE, sort = FALSE)
logmsg("identity rows", nrow(ident), "unmapped", sum(is.na(ident$patient) | ident$patient == ""))

mat <- read_geo_rds(rds_path)
if (is.data.frame(mat)) mat <- as.matrix(mat)
if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
logmsg("matrix", nrow(mat), "x", ncol(mat))

bc <- ident$barcode
common <- intersect(colnames(mat), bc)
if (length(common) < 1000) stop("barcode overlap too small: ", length(common))
mat <- mat[, common, drop = FALSE]
ident <- ident[match(common, ident$barcode), , drop = FALSE]
is_normal <- grepl("^Normal", ident$tissue)
is_normal[is.na(is_normal)] <- FALSE
keep <- !is.na(ident$lineage.sub) & ident$lineage.sub == "Malignant cells" &
  !is_normal & ident$patient %in% locked_pt
logmsg("malignant non-normal in locked patients", sum(keep))
sub <- mat[, keep, drop = FALSE]
rm(mat)
gc()
ident <- ident[keep, , drop = FALSE]
lib <- as.numeric(Matrix::colSums(sub))

Matrix::writeMM(sub, file.path(OUT, "counts.mtx"))
writeLines(colnames(sub), file.path(OUT, "barcodes.tsv"))
writeLines(rownames(sub), file.path(OUT, "features.tsv"))
meta <- data.frame(
  barcode = colnames(sub),
  unit_id = ident$patient,
  dataset = "GSE205335",
  stringsAsFactors = FALSE
)
write.table(meta, file.path(OUT, "meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(libsize = lib), file.path(OUT, "libsize.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
tab <- as.data.frame(table(ident$patient), stringsAsFactors = FALSE)
names(tab) <- c("unit_id", "n_malignant_export")
write.table(tab, file.path(OUT, "unit_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
logmsg("wrote", OUT, "cells", ncol(sub), "genes", nrow(sub))
