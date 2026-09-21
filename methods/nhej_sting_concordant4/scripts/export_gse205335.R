#!/usr/bin/env Rscript
# GSE205335 is shipped as a double-gzipped RDS. Export malignant UMI sums only.
# Patient map matches the concordant-4 Seurat script: SOFT patient id, drop
# tissue starting with "Normal", author lineage.sub == "Malignant cells".

suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
geo <- if (length(args) >= 1) args[[1]] else "/tmp/geo_nhej"
out_dir <- if (length(args) >= 2) args[[2]] else "/tmp/geo_nhej/gse205335_export"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

inner <- file.path(geo, "GSE205335_inner.rds")
if (!file.exists(inner)) {
  message("double-gunzip RDS")
  system(sprintf("gzip -dc %s > %s", file.path(geo, "GSE205335_Lung_IO_UMI_matrix.rds.gz"), inner))
}

message("readRDS")
mat <- readRDS(gzfile(inner))
message("class ", paste(class(mat), collapse = "/"), " dim ", paste(dim(mat), collapse = " x "))
if (inherits(mat, "Seurat")) stop("unexpected Seurat object")
if (!inherits(mat, "dgCMatrix")) mat <- as(as(mat, "CsparseMatrix"), "dgCMatrix")
rn <- toupper(rownames(mat))
keep <- !duplicated(rn)
mat <- mat[keep, , drop = FALSE]
rownames(mat) <- rn[keep]

norm_code <- function(x) {
  x <- toupper(gsub("-", "_", x, fixed = TRUE))
  sub("_[35]P$", "", x)
}

soft_lines <- readLines(file.path(geo, "GSE205335_family.soft.gz"))
soft <- list()
cur <- list()
for (line in soft_lines) {
  if (startsWith(line, "^SAMPLE")) {
    if (!is.null(cur$title)) soft[[length(soft) + 1]] <- cur
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
if (!is.null(cur$title)) soft[[length(soft) + 1]] <- cur

soft_map <- do.call(rbind, lapply(soft, function(s) {
  code <- if (!is.null(s$title) && grepl(" ", s$title)) sub("^\\S+\\s+", "", s$title) else NA
  data.frame(
    patient = s$patient,
    tissue = if (is.null(s$tissue)) "" else s$tissue,
    code = toupper(gsub("-", "_", code, fixed = TRUE)),
    stringsAsFactors = FALSE
  )
}))

ident <- read.delim(file.path(geo, "GSE205335_Lung_IO_CellIdentity.txt.gz"), stringsAsFactors = FALSE, check.names = FALSE)
ident$orig_code <- norm_code(ident$orig.ident)
ident <- merge(ident, soft_map, by.x = "orig_code", by.y = "code", all.x = TRUE)
ident$tissue[is.na(ident$tissue)] <- ""
ident$is_normal <- grepl("^Normal", ident$tissue)
ident$author_malignant <- ident$lineage.sub == "Malignant cells"
ident$author_tnk <- ident$lineage.total == "T/NK cells"

message("example matrix bc ", paste(head(colnames(mat), 2), collapse = " | "))
message("example ident bc ", paste(head(ident$barcode, 2), collapse = " | "))

bc_ident <- ident$barcode
common <- intersect(colnames(mat), bc_ident)
if (length(common) < 1000) {
  alt <- gsub("_", "-", colnames(mat), fixed = TRUE)
  names(alt) <- colnames(mat)
  hit <- alt %in% bc_ident
  if (sum(hit) < 1000) stop("barcode mismatch")
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
message("matched cells ", nrow(ident2))

locked_path <- Sys.getenv(
  "LOCKED_UNITS",
  unset = "/workspace/methods/nhej_sting_concordant4/data/locked_patient_units_pr539.tsv"
)
lk <- read.delim(locked_path, stringsAsFactors = FALSE)
keep_pats <- lk$unit_id[lk$dataset == "GSE205335"]

wanted_path <- Sys.getenv(
  "WANTED_GENES",
  unset = "/workspace/methods/nhej_sting_concordant4/data/wanted_genes.txt"
)
wanted <- toupper(readLines(wanted_path))
wanted <- wanted[nzchar(wanted)]
present <- intersect(unique(wanted), rownames(mat))
message("genes present ", length(present), " / requested ", length(unique(wanted)))

meta_rows <- list()
gene_rows <- list()
for (pat in keep_pats) {
  rows <- ident2[ident2$patient == pat & !ident2$is_normal, , drop = FALSE]
  n_cells <- nrow(rows)
  mal <- rows[rows$author_malignant %in% c(TRUE, "TRUE"), , drop = FALSE]
  n_mal <- nrow(mal)
  n_tnk <- sum(rows$author_tnk %in% c(TRUE, "TRUE"))
  bcs <- intersect(mal$mat_barcode, colnames(mat))
  if (!length(bcs) || n_mal == 0) {
    message("skip empty ", pat)
    next
  }
  sub <- mat[, bcs, drop = FALSE]
  lib <- sum(Matrix::colSums(sub))
  cldn4 <- if ("CLDN4" %in% rownames(sub)) as.numeric(sub["CLDN4", ]) else rep(0, ncol(sub))
  meta_rows[[pat]] <- data.frame(
    dataset = "GSE205335",
    unit_id = pat,
    n_cells = n_cells,
    n_malignant = n_mal,
    n_tnk = n_tnk,
    mal_CLDN4_pct = 100 * mean(cldn4 > 0),
    mal_CLDN4_mean = mean(log1p(cldn4)),
    lib_umi = as.numeric(lib),
    stringsAsFactors = FALSE
  )
  sums <- Matrix::rowSums(sub[present, , drop = FALSE])
  gene_rows[[pat]] <- data.frame(
    dataset = "GSE205335",
    unit_id = pat,
    gene = names(sums),
    umi = as.numeric(sums),
    stringsAsFactors = FALSE
  )
  rm(sub)
}
meta <- do.call(rbind, meta_rows)
genes <- do.call(rbind, gene_rows)
write.table(meta, file.path(out_dir, "unit_meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes, file.path(out_dir, "gene_umi.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote ", out_dir, " units ", nrow(meta), " gene rows ", nrow(genes))
