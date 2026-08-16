# Extract IMvigor210 pheno + target-gene raw counts from the public
# IMvigor210CoreBiologies package (Mariathasan et al., Nature 2018).
# The CountDataSet S4 class is not required; slots are read as attributes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript extract_imvigor.R <cds.RData> <outdir>")
}
cds_path <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

e <- new.env()
nm <- load(cds_path, envir = e)
if (!("cds" %in% nm)) {
  stop("cds object not found in ", cds_path, "; objects: ", paste(nm, collapse = ", "))
}
cds <- e$cds
cnt <- get("counts", envir = attr(cds, "assayData"))
pd <- attr(attr(cds, "phenoData"), "data")
fd <- attr(attr(cds, "featureData"), "data")

genes <- c(
  "TACSTD2", "CLDN4", "CD274", "CD8A", "GZMA", "PRF1",
  "EPCAM", "KRT5", "KRT20", "GATA3", "FOXA1", "UPK2",
  "CDH1", "CXCL9", "CXCL13", "TGFB1"
)
keep <- fd$symbol %in% genes
if (sum(keep) != length(genes)) {
  missing <- setdiff(genes, fd$symbol[keep])
  stop("Missing genes in featureData: ", paste(missing, collapse = ", "))
}

lib <- colSums(cnt)
sub <- as.data.frame(t(cnt[keep, , drop = FALSE]))
colnames(sub) <- fd$symbol[keep]
sub$sample_id <- rownames(sub)
sub$libsize <- as.numeric(lib[rownames(sub)])

pd_out <- pd
pd_out$sample_id <- rownames(pd)

write.csv(pd_out, file.path(outdir, "imvigor210_pdata.csv"), row.names = FALSE)
write.csv(sub, file.path(outdir, "imvigor210_gene_counts.csv"), row.names = FALSE)
write.csv(fd[keep, ], file.path(outdir, "imvigor210_gene_annot.csv"), row.names = TRUE)
cat("Wrote IMvigor210 extracts to ", outdir, "\n", sep = "")
cat("samples=", nrow(pd), " genes=", sum(keep), "\n", sep = "")
