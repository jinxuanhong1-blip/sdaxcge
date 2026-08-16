# Dump counts, phenotype, and feature tables from IMvigor210CoreBiologies
# cds (CountDataSet / DESeq) without installing Bioconductor.
# S4 slots are stored as attributes; assayData is an environment.
#
# Usage: Rscript extract_cds.R <cds.RData> <out_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript extract_cds.R <cds.RData> <out_dir>")
}
rdata <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(rdata)
if (!exists("cds")) stop("cds object not found in ", rdata)

ad <- attr(cds, "assayData")
counts <- get("counts", envir = ad)
pd <- attr(attr(cds, "phenoData"), "data")
fd <- attr(attr(cds, "featureData"), "data")

write.csv(pd, file.path(out_dir, "pData.csv"))
write.csv(fd, file.path(out_dir, "fData.csv"))
gz <- gzfile(file.path(out_dir, "counts.csv.gz"), "w")
write.csv(counts, gz)
close(gz)

cat(sprintf(
  "Wrote pData %s, fData %s, counts %s x %s\n",
  nrow(pd), nrow(fd), nrow(counts), ncol(counts)
))
