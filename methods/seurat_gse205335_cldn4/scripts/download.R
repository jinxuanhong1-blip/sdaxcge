#!/usr/bin/env Rscript
# Public GEO processed files for GSE205335 only. No EGA raw. No GSE148071.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "/tmp/gse205335"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

files <- c(
  "GSE205335_Lung_IO_CellIdentity.txt.gz" =
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
  "GSE205335_Lung_IO_UMI_matrix.rds.gz" =
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
  "GSE205335_family.soft.gz" =
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/soft/GSE205335_family.soft.gz"
)

for (name in names(files)) {
  dest <- file.path(out, name)
  if (file.exists(dest) && file.info(dest)$size > 1000) {
    message("have ", dest, " (", file.info(dest)$size, " bytes)")
    next
  }
  message("GET ", files[[name]])
  utils::download.file(files[[name]], destfile = dest, mode = "wb", quiet = FALSE)
  message("wrote ", dest, " (", file.info(dest)$size, " bytes)")
}
