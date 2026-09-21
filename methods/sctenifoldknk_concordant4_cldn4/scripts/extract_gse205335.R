#!/usr/bin/env Rscript
# GSE205335 is a double-gzipped dgCMatrix. Extract tumor malignant cells
# for the locked patients. Normal-tissue samples are excluded.
# Pseudobulk = all author-malignant tumor cells. GRN = QC malignant, <=100/patient.

.libPaths(c("/home/ubuntu/R/library", .libPaths()))
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
geo <- if (length(args) >= 1) args[[1]] else "/tmp/geo_c4"
out <- if (length(args) >= 2) args[[2]] else "/tmp/c4_work/GSE205335"
root <- Sys.getenv("C4_ROOT", unset = "/workspace/methods/sctenifoldknk_concordant4_cldn4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

locked <- read.delim(file.path(root, "data/locked_patient_units.tsv"), stringsAsFactors = FALSE)
keep <- locked$unit_id[locked$dataset == "GSE205335"]
gsm <- read.delim(file.path(root, "data/GSE205335_gsm_map.tsv"), stringsAsFactors = FALSE)
if (any(duplicated(gsm$orig.ident))) stop("duplicated orig.ident in gsm map")

ident <- read.delim(gzfile(file.path(geo, "GSE205335_Lung_IO_CellIdentity.txt.gz")), stringsAsFactors = FALSE)
ident <- merge(ident, gsm[, c("orig.ident", "patient", "tissue")], by = "orig.ident", all.x = TRUE)
ident$tissue[is.na(ident$tissue)] <- ""
ident$is_normal <- grepl("^Normal", ident$tissue)
ident$author_malignant <- ident$lineage.sub == "Malignant cells"
ident$author_tnk <- ident$lineage.total == "T/NK cells"

plain <- file.path(dirname(out), "GSE205335_Lung_IO_UMI_matrix.rds")
gz <- file.path(geo, "GSE205335_Lung_IO_UMI_matrix.rds.gz")
if (!file.exists(plain)) {
  message("double-gunzip RDS")
  status <- system(sprintf("gzip -dc %s | gzip -dc > %s", shQuote(gz), shQuote(plain)))
  if (status != 0) stop("double gunzip failed")
}
message("readRDS")
mat <- readRDS(plain)
if (!inherits(mat, "dgCMatrix")) mat <- as(as(mat, "CsparseMatrix"), "dgCMatrix")
message("matrix ", nrow(mat), " x ", ncol(mat), " example col ", colnames(mat)[1])

bc <- ident$barcode
common <- intersect(colnames(mat), bc)
if (length(common) < 1000) {
  alt <- gsub("_", "-", colnames(mat))
  names(alt) <- colnames(mat)
  common <- colnames(mat)[alt %in% bc]
  if (length(common) < 1000) {
    stop("barcode mismatch matrix ", colnames(mat)[1], " ident ", bc[1])
  }
  rownames(ident) <- ident$barcode
  ident2 <- ident[alt[common], ]
  ident2$mat_barcode <- common
} else {
  rownames(ident) <- ident$barcode
  ident2 <- ident[common, ]
  ident2$mat_barcode <- common
}
message("matched cells ", nrow(ident2))

# Uppercase symbols. Collapse duplicates without make.names (keeps hyphens).
rn <- toupper(rownames(mat))
if (any(duplicated(rn))) {
  message("collapsing ", sum(duplicated(rn)), " duplicate symbols")
  fac <- factor(rn, levels = unique(rn))
  design <- sparseMatrix(
    i = seq_along(rn), j = as.integer(fac), x = 1,
    dims = c(length(rn), nlevels(fac))
  )
  mat <- as(Matrix::t(design) %*% mat, "dgCMatrix")
  rownames(mat) <- levels(fac)
} else {
  rownames(mat) <- rn
}
if (!"CLDN4" %in% rownames(mat)) stop("CLDN4 missing after symbol collapse")

tumor <- ident2[ident2$patient %in% keep & !ident2$is_normal, ]
tumor <- tumor[!is.na(tumor$patient), ]
message("tumor cells in locked patients ", nrow(tumor))
message("lineage.sub malignant ", sum(tumor$author_malignant))

n_cells <- as.integer(table(factor(tumor$patient, levels = keep)))
names(n_cells) <- keep
n_tnk <- as.integer(table(factor(tumor$patient[tumor$author_tnk], levels = keep)))
names(n_tnk) <- keep

mal <- tumor[tumor$author_malignant, ]
mal <- mal[mal$mat_barcode %in% colnames(mat), ]
message("malignant barcodes ", nrow(mal))
mat_mal <- mat[, mal$mat_barcode, drop = FALSE]
rm(mat)
gc(verbose = FALSE)

fac <- factor(mal$patient, levels = keep)
group <- Matrix::sparse.model.matrix(~ 0 + fac)
colnames(group) <- keep
pb <- as.matrix(mat_mal %*% group)
rownames(pb) <- rownames(mat_mal)

cl <- as.numeric(mat_mal["CLDN4", ])
pct <- tapply(cl > 0, fac, function(z) 100 * mean(z))
meanlog <- tapply(log1p(cl), fac, mean)
n_mal <- as.integer(table(fac))

ncount <- Matrix::colSums(mat_mal)
nfeat <- Matrix::colSums(mat_mal > 0)
mt <- grep("^MT-", rownames(mat_mal), value = TRUE)
mt_pct <- if (length(mt)) 100 * Matrix::colSums(mat_mal[mt, , drop = FALSE]) / pmax(ncount, 1) else rep(0, ncol(mat_mal))
qc <- (nfeat >= 200) & (ncount >= 500) & (mt_pct < 20)
n_qc <- as.integer(tapply(qc, fac, sum))

set.seed(1)
keep_cells <- integer(0)
keep_unit <- character(0)
for (pat in keep) {
  if (!is.finite(n_mal[[pat]]) || n_mal[[pat]] < 30) next
  ix <- which(fac == pat & qc)
  if (!length(ix)) next
  if (length(ix) > 100) ix <- sample(ix, 100)
  keep_cells <- c(keep_cells, ix)
  keep_unit <- c(keep_unit, rep(pat, length(ix)))
}
grn <- mat_mal[, keep_cells, drop = FALSE]
message("GRN cells ", ncol(grn), " genes ", nrow(grn))

units <- data.frame(
  dataset = "GSE205335",
  unit_id = keep,
  n_cells = as.integer(n_cells[keep]),
  n_malignant = as.integer(n_mal[keep]),
  n_tnk = as.integer(n_tnk[keep]),
  n_malignant_qc = as.integer(n_qc[keep]),
  n_grn = as.integer(table(factor(keep_unit, levels = keep))),
  mal_CLDN4_pct = as.numeric(pct[keep]),
  mal_CLDN4_mean = as.numeric(meanlog[keep]),
  stringsAsFactors = FALSE
)
write.table(units, file.path(out, "units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pb_out <- data.frame(gene = rownames(pb), pb, check.names = FALSE)
colnames(pb_out) <- c("gene", keep)
gz <- gzfile(file.path(out, "pseudobulk.tsv.gz"), "wt")
write.table(pb_out, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

writeMM(grn, file.path(out, "grn_matrix.mtx"))
writeLines(rownames(grn), file.path(out, "grn_genes.tsv"))
write.table(
  data.frame(unit_id = keep_unit),
  file.path(out, "grn_cells.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("wrote ", out)
unlink(plain)
message("removed uncompressed RDS")
