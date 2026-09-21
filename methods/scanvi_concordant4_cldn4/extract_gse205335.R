#!/usr/bin/env Rscript
# GSE205335: author labels, patient is the locked unit (libraries pooled, Normal* out).
# Writes full-unit CLDN4 / T/NK table and a QC-capped subsample matrix for scVI.

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
raw <- if (length(args) >= 1) args[[1]] else "/tmp/geo_c4"
out <- if (length(args) >= 2) args[[2]] else "/workspace/methods/scanvi_concordant4_cldn4/results/cache/GSE205335"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

MAL_CAP <- 220L
TNK_CAP <- 140L
OTHER_CAP <- 40L
QC_MIN_GENES <- 200
QC_MIN_UMI <- 500
QC_MAX_MT <- 20

patients_path <- "/workspace/methods/scanvi_concordant4_cldn4/data/GSE205335_patients.tsv"
keep_pats <- read.delim(patients_path, stringsAsFactors = FALSE)$patient

ident <- read.delim(gzfile(file.path(raw, "GSE205335_Lung_IO_CellIdentity.txt.gz")), stringsAsFactors = FALSE)
soft_lines <- readLines(gzfile(file.path(raw, "GSE205335_family.soft.gz")))
soft <- list()
cur <- list()
for (line in soft_lines) {
  line <- sub("\r$", "", line)
  if (startsWith(line, "^SAMPLE")) {
    if (!is.null(cur$title)) soft[[length(soft) + 1]] <- cur
    cur <- list()
  } else if (startsWith(line, "!Sample_title = ")) {
    cur$title <- sub("^!Sample_title = ", "", line)
  } else if (startsWith(line, "!Sample_characteristics_ch1 = ")) {
    val <- sub("^!Sample_characteristics_ch1 = ", "", line)
    if (grepl(": ", val, fixed = TRUE)) {
      kv <- strsplit(val, ": ", fixed = TRUE)[[1]]
      cur[[kv[[1]]]] <- paste(kv[-1], collapse = ": ")
    }
  }
}
if (!is.null(cur$title)) soft[[length(soft) + 1]] <- cur

norm_code <- function(x) {
  x <- toupper(gsub("-", "_", x, fixed = TRUE))
  sub("_[35]P$", "", x)
}

soft_map <- do.call(rbind, lapply(soft, function(s) {
  code <- if (!is.null(s$title) && grepl(" ", s$title, fixed = TRUE)) sub("^\\S+\\s+", "", s$title) else NA_character_
  data.frame(
    patient = if (is.null(s$patient)) NA_character_ else s$patient,
    tissue = if (is.null(s$tissue)) NA_character_ else s$tissue,
    code = toupper(gsub("-", "_", code, fixed = TRUE)),
    stringsAsFactors = FALSE
  )
}))
soft_map <- soft_map[!is.na(soft_map$code) & !duplicated(soft_map$code), ]

ident$orig_code <- norm_code(ident$orig.ident)
ident <- merge(ident, soft_map, by.x = "orig_code", by.y = "code", all.x = TRUE, sort = FALSE)
ident$tissue[is.na(ident$tissue)] <- ""
ident$is_normal <- grepl("^Normal", ident$tissue)
ident$author_malignant <- ident$lineage.sub == "Malignant cells"
ident$author_malignant[is.na(ident$author_malignant)] <- FALSE
ident$author_tnk <- ident$lineage.total == "T/NK cells"
ident$author_tnk[is.na(ident$author_tnk)] <- FALSE
say("identity rows", nrow(ident), "with patient", sum(!is.na(ident$patient)), "kept patients target", length(keep_pats))

rds_gz <- file.path(raw, "GSE205335_Lung_IO_UMI_matrix.rds.gz")
say("readRDS via gzcon(gzfile) — double gzip")
mat <- tryCatch(
  readRDS(gzcon(gzfile(rds_gz, "rb"))),
  error = function(e) {
    say("gzcon failed:", conditionMessage(e), "-- expanding to plain RDS")
    plain <- file.path(raw, "GSE205335_plain.rds")
    if (!file.exists(plain)) {
      status <- system2("bash", c("-c", sprintf("gzip -dc %s | gzip -dc > %s", shQuote(rds_gz), shQuote(plain))))
      if (status != 0) stop("double gunzip failed")
    }
    readRDS(plain)
  }
)
if (!inherits(mat, "dgCMatrix")) mat <- as(as(mat, "CsparseMatrix"), "dgCMatrix")
say("matrix", nrow(mat), "x", ncol(mat), "example col", colnames(mat)[1])

rn <- toupper(rownames(mat))
keep_g <- !duplicated(rn)
mat <- mat[keep_g, , drop = FALSE]
rownames(mat) <- rn[keep_g]
if (!"CLDN4" %in% rownames(mat)) stop("CLDN4 missing from GSE205335")
if (!"TACSTD2" %in% rownames(mat)) stop("TACSTD2 missing from GSE205335")

bc_ident <- ident$barcode
common <- intersect(colnames(mat), bc_ident)
if (length(common) < 1000) {
  alt <- gsub("_", "-", colnames(mat), fixed = TRUE)
  names(alt) <- colnames(mat)
  hit <- alt %in% bc_ident
  if (sum(hit) < 1000) {
    alt2 <- gsub("-", "_", colnames(mat), fixed = TRUE)
    names(alt2) <- colnames(mat)
    hit <- alt2 %in% bc_ident
    alt <- alt2
  }
  if (sum(hit) < 1000) {
    stop("barcode mismatch matrix ", colnames(mat)[1], " ident ", bc_ident[1], " n_common ", length(common))
  }
  map <- ident
  rownames(map) <- map$barcode
  ident2 <- map[alt[hit], , drop = FALSE]
  ident2$mat_barcode <- colnames(mat)[hit]
} else {
  map <- ident
  rownames(map) <- map$barcode
  ident2 <- map[common, , drop = FALSE]
  ident2$mat_barcode <- common
}
say("matched cells", nrow(ident2), "patients", length(unique(ident2$patient)))

qc_pass <- function(sub) {
  ncount <- Matrix::colSums(sub)
  nfeat <- Matrix::colSums(sub > 0)
  mt <- grep("^MT-", rownames(sub), value = TRUE)
  mt_pct <- if (length(mt)) 100 * Matrix::colSums(sub[mt, , drop = FALSE]) / pmax(ncount, 1) else rep(0, ncol(sub))
  nfeat >= QC_MIN_GENES & ncount >= QC_MIN_UMI & mt_pct < QC_MAX_MT
}

stable_keep <- function(keys) {
  if (!length(keys)) return(logical())
  script <- "/workspace/methods/scanvi_concordant4_cldn4/stable_keep.py"
  tf <- tempfile()
  on.exit(unlink(tf), add = TRUE)
  writeLines(keys, tf)
  out <- system(sprintf("python3 %s < %s", shQuote(script), shQuote(tf)), intern = TRUE)
  if (length(out) != length(keys)) stop("stable_keep length mismatch")
  out == "1"
}
unit_rows <- list()
meta_rows <- list()
keep_barcodes <- character()

for (pat in keep_pats) {
  rows <- ident2[ident2$patient == pat & !is.na(ident2$patient), , drop = FALSE]
  if (!nrow(rows)) {
    say("  missing patient", pat)
    next
  }
  tumor <- rows[!rows$is_normal, , drop = FALSE]
  if (!nrow(tumor)) tumor <- rows
  n_cells <- nrow(tumor)
  n_mal <- sum(tumor$author_malignant)
  n_tnk <- sum(tumor$author_tnk)
  mal_bc <- tumor$mat_barcode[tumor$author_malignant]
  cldn4 <- if (length(mal_bc)) as.numeric(mat["CLDN4", mal_bc]) else numeric()
  tacstd2 <- if (length(mal_bc)) as.numeric(mat["TACSTD2", mal_bc]) else numeric()
  unit_rows[[pat]] <- data.frame(
    dataset = "GSE205335",
    unit_id = pat,
    patient_id = pat,
    unit_type = "patient",
    tissue = paste(sort(unique(tumor$tissue)), collapse = ","),
    n_cells = n_cells,
    n_malignant = n_mal,
    n_tnk = n_tnk,
    frac_tnk = if (n_cells) n_tnk / n_cells else NA_real_,
    mal_CLDN4_pct = if (length(cldn4)) 100 * mean(cldn4 > 0) else NA_real_,
    mal_CLDN4_mean_log1p = if (length(cldn4)) mean(log1p(cldn4)) else NA_real_,
    mal_TACSTD2_pct = if (length(tacstd2)) 100 * mean(tacstd2 > 0) else NA_real_,
    mal_TACSTD2_mean_log1p = if (length(tacstd2)) mean(log1p(tacstd2)) else NA_real_,
    malig_def = "author_malig",
    stringsAsFactors = FALSE
  )
  bcs <- tumor$mat_barcode
  ok <- qc_pass(mat[, bcs, drop = FALSE])
  bcs_qc <- bcs[ok]
  tumor_qc <- tumor[match(bcs_qc, tumor$mat_barcode), , drop = FALSE]
  chosen_bc <- bcs_qc[stable_keep(paste0(pat, "|", bcs_qc))]
  if (!length(chosen_bc)) next
  chosen <- tumor_qc[match(chosen_bc, tumor_qc$mat_barcode), , drop = FALSE]
  chosen$cell_id <- paste0("GSE205335|", pat, "|", chosen$mat_barcode)
  meta_rows[[pat]] <- data.frame(
    cell_id = chosen$cell_id,
    dataset = "GSE205335",
    unit_id = pat,
    patient_id = pat,
    unit_type = "patient",
    barcode = chosen$mat_barcode,
    seed_class = ifelse(chosen$author_malignant, "malignant", ifelse(chosen$author_tnk, "T/NK", "other")),
    stringsAsFactors = FALSE
  )
  keep_barcodes <- c(keep_barcodes, chosen$mat_barcode)
  say("  ", pat, "cells", n_cells, "mal", n_mal, "tnk", n_tnk, "sub", nrow(chosen),
      "cldn4", round(unit_rows[[pat]]$mal_CLDN4_pct, 2),
      "tacstd2", round(unit_rows[[pat]]$mal_TACSTD2_pct, 2))
}

units <- do.call(rbind, unit_rows)
meta <- do.call(rbind, meta_rows)
rownames(meta) <- NULL
# unique cells if a barcode was selected once
meta <- meta[!duplicated(meta$barcode), , drop = FALSE]
sub <- mat[, meta$barcode, drop = FALSE]
colnames(sub) <- meta$cell_id
rm(mat)
gc(verbose = FALSE)

write.table(units, file.path(out, "units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta, file.path(out, "meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
Matrix::writeMM(sub, file.path(out, "matrix.mtx"))
writeLines(rownames(sub), file.path(out, "features.tsv"))
writeLines(colnames(sub), file.path(out, "barcodes.tsv"))
say("wrote", out, "units", nrow(units), "subcells", ncol(sub), "genes", nrow(sub))
