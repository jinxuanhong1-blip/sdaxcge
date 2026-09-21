#!/usr/bin/env Rscript
# Slim export of GSE205335 raw counts for the state panel.
# Author labels. Normal tissue dropped. Locked patients only.
# Does not score. extract_scores.py applies the same log1p(CP10k) math as the other cohorts.

suppressPackageStartupMessages(library(Matrix))
options(warn = 1, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
raw <- if (length(args) >= 1) args[[1]] else "/tmp/concordant4_raw"
out <- if (length(args) >= 2) args[[2]] else "/tmp/concordant4_state/gse205335_cells.tsv.gz"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

panel <- c(
  "CLDN4", "EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC",
  "CD3D", "CD3E", "CD8A", "CD8B", "CD4",
  "KLRD1", "NCR1", "KLRF1", "NCAM1", "GNLY", "NKG7",
  "GZMB", "PRF1", "IFNG", "GZMA", "GZMH",
  "PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CTLA4", "ENTPD1", "CXCL13",
  "IL7R", "TCF7", "CCR7", "SELL", "LEF1",
  "CX3CR1", "FGFBP2"
)

locked <- c(
  "P0031", "P1006", "P1015", "P1016", "P1017", "P1018", "P1025", "P1027",
  "P1030", "P1037", "P1056", "P1062", "P1063", "P1072", "P1076", "P1079",
  "P1084", "P1089", "P1090", "P1115", "P1119", "P4001"
)

say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

norm_code <- function(x) {
  x <- toupper(gsub("-", "_", x, fixed = TRUE))
  sub("_[35]P$", "", x)
}

say("parse SOFT")
soft_lines <- readLines(file.path(raw, "GSE205335", "GSE205335_family.soft.gz"))
soft <- list()
cur <- list()
for (line in soft_lines) {
  if (startsWith(line, "^SAMPLE")) {
    if (!is.null(cur$title)) soft[[length(soft) + 1L]] <- cur
    cur <- list()
  } else if (startsWith(line, "!Sample_title")) {
    cur$title <- sub("^!Sample_title = ", "", line)
  } else if (startsWith(line, "!Sample_characteristics_ch1")) {
    val <- sub("^!Sample_characteristics_ch1 = ", "", line)
    if (grepl(": ", val, fixed = TRUE)) {
      kv <- strsplit(val, ": ", fixed = TRUE)[[1]]
      cur[[kv[[1]]]] <- paste(kv[-1], collapse = ": ")
    }
  }
}
if (!is.null(cur$title)) soft[[length(soft) + 1L]] <- cur

soft_map <- do.call(rbind, lapply(soft, function(s) {
  code <- sub("^\\S+\\s+", "", s$title)
  data.frame(
    title = s$title,
    patient = s$patient,
    tissue = s$tissue,
    code = toupper(gsub("-", "_", code, fixed = TRUE)),
    stringsAsFactors = FALSE
  )
}))
if (anyDuplicated(soft_map$code)) stop("duplicated SOFT sample codes")

say("identity")
ident <- read.delim(
  file.path(raw, "GSE205335", "GSE205335_Lung_IO_CellIdentity.txt.gz"),
  stringsAsFactors = FALSE
)
ident$orig_code <- norm_code(ident$orig.ident)
ident <- merge(ident, soft_map, by.x = "orig_code", by.y = "code", all.x = TRUE)
if (anyNA(ident$patient)) {
  say("WARN unmapped identity rows", sum(is.na(ident$patient)))
}
ident$is_normal <- grepl("^Normal", ident$tissue)
ident <- ident[ident$patient %in% locked & !ident$is_normal & !is.na(ident$patient), ]
say("identity rows after patient/normal filter", nrow(ident))

rds <- file.path(raw, "GSE205335", "GSE205335_Lung_IO_UMI_matrix.rds")
if (!file.exists(rds)) stop("missing plain RDS: ", rds)
say("readRDS")
mat <- readRDS(rds)
if (!inherits(mat, "dgCMatrix")) mat <- as(as(mat, "CsparseMatrix"), "dgCMatrix")
say("matrix", nrow(mat), "x", ncol(mat), "class", paste(class(mat), collapse = ","))

rn <- toupper(rownames(mat))
keep_rn <- !duplicated(rn)
if (!all(keep_rn)) {
  say("drop duplicated rownames", sum(!keep_rn))
  mat <- mat[keep_rn, , drop = FALSE]
  rn <- rn[keep_rn]
}
rownames(mat) <- rn
present <- intersect(panel, rownames(mat))
missing <- setdiff(panel, present)
say("panel present", length(present), "missing", paste(missing, collapse = ","))
if (!"CLDN4" %in% present) stop("CLDN4 absent")
if (!"GZMB" %in% present || !"PDCD1" %in% present) stop("effector/exhaustion gene absent")

say("library size")
lib_all <- as.numeric(Matrix::colSums(mat))
names(lib_all) <- colnames(mat)

bc_mat <- colnames(mat)
m <- match(ident$barcode, bc_mat)
say("barcode match", sum(!is.na(m)), "/", nrow(ident))
if (sum(!is.na(m)) < 0.9 * nrow(ident)) {
  say("example ident", ident$barcode[1], "example matrix", bc_mat[1])
  stop("barcode match below 90%")
}
ident <- ident[!is.na(m), ]
ident$mat_col <- bc_mat[m[!is.na(m)]]

comp <- rep(NA_character_, nrow(ident))
is_amb <- ident$celltype == "AMB cells"
comp[ident$lineage.sub == "Malignant cells"] <- "malignant"
comp[ident$lineage.sub == "CD8+ T cells" & !is_amb] <- "cd8"
comp[ident$lineage.sub == "NK cells" & !is_amb] <- "nk"
say("AMB dropped", sum(is_amb & ident$lineage.sub %in% c("CD8+ T cells", "NK cells")))
ident$compartment <- comp
ident <- ident[!is.na(ident$compartment), ]
say("kept cells", nrow(ident), "by compartment:")
print(table(ident$compartment))

say("subset panel counts")
sub <- mat[present, ident$mat_col, drop = FALSE]
rm(mat)
gc(verbose = FALSE)
dense <- as.matrix(sub)
rm(sub)
gc(verbose = FALSE)

slice <- ident$celltype
slice[is.na(slice) | slice == ""] <- "unassigned"
meta <- data.frame(
  dataset = "GSE205335",
  unit_id = ident$patient,
  compartment = ident$compartment,
  author_subset = slice,
  libsize = lib_all[ident$mat_col],
  stringsAsFactors = FALSE
)
# dense is genes x cells
stopifnot(ncol(dense) == nrow(meta))
out_df <- cbind(meta, t(dense))
colnames(out_df) <- c(colnames(meta), present)

say("write", out, "rows", nrow(out_df))
con <- gzfile(out, "wt")
write.table(out_df, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
say("done", file.info(out)$size, "bytes")
