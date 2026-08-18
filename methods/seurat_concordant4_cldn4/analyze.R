#!/usr/bin/env Rscript
# Seurat / Harmony primary for the concordant-4 CLDN4-only pool.
# Honest unit = patient / donor / sample. Do not cite cell n as the test n.

suppressPackageStartupMessages({
  .libPaths(c("/home/ubuntu/R/library", .libPaths()))
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

options(warn = 1)
set.seed(1)

ROOT <- "/workspace/methods/seurat_concordant4_cldn4"
GEO <- "/tmp/geo_seurat"
OUT_FIG <- file.path(ROOT, "results", "figures")
OUT_TAB <- file.path(ROOT, "results", "tables")
OUT_OBJ <- file.path(ROOT, "results", "objects")
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_OBJ, recursive = TRUE, showWarnings = FALSE)

CAP <- 350L
SEED <- 1L
QC_MIN_GENES <- 200L
QC_MIN_UMI <- 500L
QC_MAX_MT <- 20
PCA_DIMS <- 1:30

LOCKED_123902 <- file.path(ROOT, "data", "GSE123902_marker_units.tsv")
LOCKED_205335 <- file.path(ROOT, "data", "GSE205335_patients.tsv")
LOCKED_189357 <- file.path(ROOT, "data", "GSE189357_marker_units.tsv")
SETS_JSON <- file.path(ROOT, "data", "a8_sets.json")

say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

`%||%` <- function(a, b) if (!is.null(a)) a else b

upper_rownames <- function(mat) {
  rn <- toupper(rownames(mat))
  # keep first occurrence of duplicated symbols after uppercasing
  keep <- !duplicated(rn)
  mat <- mat[keep, , drop = FALSE]
  rownames(mat) <- rn[keep]
  mat
}

gene_or_zero <- function(mat, gene) {
  gene <- toupper(gene)
  if (gene %in% rownames(mat)) as.numeric(mat[gene, ]) else rep(0, ncol(mat))
}

marker_classes <- function(mat) {
  epcam <- gene_or_zero(mat, "EPCAM")
  krt8 <- gene_or_zero(mat, "KRT8")
  krt18 <- gene_or_zero(mat, "KRT18")
  krt19 <- gene_or_zero(mat, "KRT19")
  ptprc <- gene_or_zero(mat, "PTPRC")
  cd3d <- gene_or_zero(mat, "CD3D")
  cd3e <- gene_or_zero(mat, "CD3E")
  cd8a <- gene_or_zero(mat, "CD8A")
  nkg7 <- gene_or_zero(mat, "NKG7")
  gnly <- gene_or_zero(mat, "GNLY")
  klrd1 <- gene_or_zero(mat, "KLRD1")
  mal <- (epcam > 0 | krt8 > 0 | krt18 > 0 | krt19 > 0) & ptprc == 0
  tnk <- (cd3d > 0 | cd3e > 0 | cd8a > 0 | nkg7 > 0 | gnly > 0 | klrd1 > 0) & !mal
  list(malignant = mal, tnk = tnk)
}

qc_pass <- function(mat) {
  nfeat <- Matrix::colSums(mat > 0)
  ncount <- Matrix::colSums(mat)
  mt <- grep("^MT-", rownames(mat), value = TRUE)
  mt_pct <- if (length(mt)) 100 * Matrix::colSums(mat[mt, , drop = FALSE]) / pmax(ncount, 1) else rep(0, ncol(mat))
  nfeat >= QC_MIN_GENES & ncount >= QC_MIN_UMI & mt_pct < QC_MAX_MT
}

cap_cells <- function(cells, cap = CAP, seed = SEED) {
  if (length(cells) <= cap) return(cells)
  set.seed(seed)
  sample(cells, cap)
}

as_dgC <- function(mat) {
  if (!inherits(mat, "dgCMatrix")) mat <- as(as(mat, "CsparseMatrix"), "dgCMatrix")
  mat
}

seurat_from_mat <- function(mat, meta, project) {
  mat <- as_dgC(mat)
  colnames(mat) <- make.unique(colnames(mat))
  rownames(meta) <- colnames(mat)
  CreateSeuratObject(counts = mat, meta.data = meta, project = project, min.cells = 0, min.features = 0)
}

read_sets <- function(path) {
  js <- fromJSON(path, simplifyVector = TRUE)
  ifn <- unique(toupper(c(
    js$sets$HALLMARK_INTERFERON_ALPHA_RESPONSE,
    js$sets$HALLMARK_INTERFERON_GAMMA_RESPONSE
  )))
  mhc <- unique(toupper(js$sets$CUSTOM_MHC_I_ANTIGEN_PRESENTATION))
  tj <- unique(toupper(c(
    js$sets$KEGG_TIGHT_JUNCTION,
    js$sets$GOBP_TIGHT_JUNCTION_ORGANIZATION
  )))
  extras <- setdiff(toupper(c("CDH1", "VIM", "ZEB1")), c("CLDN4", "KRT5", "KRT7", "KRT8", "KRT17", "KRT18", "KRT19"))
  tj <- setdiff(unique(c(tj, extras)), "CLDN4")
  list(IFN = ifn, MHC = mhc, TJ = tj)
}

tmm_factors <- function(counts) {
  counts <- as.matrix(counts)
  lib <- colSums(counts)
  lib[lib == 0] <- 1
  rel <- sweep(counts, 2, lib, "/")
  q75 <- apply(rel, 2, function(x) {
    x <- x[x > 0]
    if (!length(x)) 0 else as.numeric(stats::quantile(x, 0.75, names = FALSE))
  })
  ref <- which.min(abs(q75 - mean(q75)))
  sf <- vapply(seq_len(ncol(counts)), function(j) {
    if (j == ref) return(1)
    keep <- counts[, j] > 0 & counts[, ref] > 0
    if (sum(keep) < 20) return(1)
    m <- log2((counts[keep, j] / lib[j]) / (counts[keep, ref] / lib[ref]))
    a <- 0.5 * log2((counts[keep, j] / lib[j]) * (counts[keep, ref] / lib[ref]))
    n <- length(m)
    keep_m <- order(m)[ceiling(0.3 * n):max(ceiling(0.3 * n), floor(0.7 * n))]
    keep_a <- order(a)[ceiling(0.05 * n):max(ceiling(0.05 * n), floor(0.95 * n))]
    keep2 <- intersect(keep_m, keep_a)
    if (length(keep2) < 10) keep2 <- seq_len(n)
    w <- 1 / (1 / pmax(counts[keep, j][keep2], 1) + 1 / pmax(counts[keep, ref][keep2], 1))
    as.numeric(2^stats::weighted.mean(m[keep2], w))
  }, numeric(1))
  sf / mean(sf)
}

dl_spearman <- function(rhos, ns) {
  ok <- is.finite(rhos) & is.finite(ns) & ns > 3
  rhos <- rhos[ok]
  ns <- ns[ok]
  z <- atanh(pmin(pmax(rhos, -0.999999), 0.999999))
  var_z <- 1 / (ns - 3)
  w <- 1 / var_z
  zbar <- sum(w * z) / sum(w)
  q <- sum(w * (z - zbar)^2)
  k <- length(rhos)
  dfree <- k - 1
  cdenom <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (dfree > 0 && cdenom > 0) max(0, (q - dfree) / cdenom) else 0
  wstar <- 1 / (var_z + tau2)
  zre <- sum(wstar * z) / sum(wstar)
  se <- sqrt(1 / sum(wstar))
  p <- 2 * stats::pnorm(-abs(zre / se))
  i2 <- if (q > 0) max(0, (q - dfree) / q) else 0
  ci <- tanh(zre + c(-1, 1) * 1.96 * se)
  list(rho = tanh(zre), p = p, I2 = i2, ci_lo = ci[1], ci_hi = ci[2], k = k, N = sum(ns))
}

within_quartile <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  qs <- tryCatch(
    as.character(cut(r, breaks = stats::quantile(r, probs = seq(0, 1, 0.25), na.rm = TRUE, type = 7),
                     include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))),
    error = function(e) rep(NA_character_, length(x))
  )
  qs
}

rank_biserial <- function(x_q4, x_q1) {
  wt <- stats::wilcox.test(x_q4, x_q1, alternative = "two.sided", exact = FALSE)
  n4 <- length(x_q4)
  n1 <- length(x_q1)
  u <- as.numeric(wt$statistic)
  r <- 2 * u / (n4 * n1) - 1
  list(r = r, p = wt$p.value, n_q1 = n1, n_q4 = n4)
}

family_present <- function(genes, universe) intersect(toupper(genes), toupper(universe))

# ---------------------------------------------------------------------------
# Gene sets
# ---------------------------------------------------------------------------
sets <- read_sets(SETS_JSON)
say("sets IFN", length(sets$IFN), "MHC", length(sets$MHC), "TJ", length(sets$TJ))

# ---------------------------------------------------------------------------
# GSE123902 (Laughney) — donor, marker gate, drop NORMAL
# ---------------------------------------------------------------------------
say("GSE123902 extract")
d123 <- file.path(GEO, "gse123902")
dir.create(d123, showWarnings = FALSE)
if (length(list.files(d123, pattern = "dense\\.csv\\.gz$")) < 17) {
  system(sprintf("tar -C %s -xf %s", d123, file.path(GEO, "GSE123902_RAW.tar")))
}
locked_123 <- read.delim(LOCKED_123902, stringsAsFactors = FALSE)
# one tumor row per donor; drop NORMAL; eligible n_mal>=20 & n_tnk>=20
tumor_123 <- locked_123[locked_123$tissue %in% c("PRIMARY", "METASTASIS"), ]
tumor_123 <- tumor_123[order(tumor_123$patient, tumor_123$tissue), ]
tumor_123 <- tumor_123[!duplicated(tumor_123$patient), ]
tumor_123 <- tumor_123[tumor_123$n_malignant >= 20 & tumor_123$n_tnk >= 20, ]

units_123 <- list()
pb_123 <- list()
objs <- list()

for (i in seq_len(nrow(tumor_123))) {
  fn <- tumor_123$file[i]
  donor <- tumor_123$patient[i]
  say("  GSE123902", donor, fn)
  path <- file.path(d123, fn)
  raw <- as.matrix(utils::read.csv(path, row.names = 1, check.names = FALSE))
  mat <- t(raw)
  storage.mode(mat) <- "double"
  mat <- as_dgC(mat)
  mat <- upper_rownames(mat)
  cls <- marker_classes(mat)
  cldn4 <- gene_or_zero(mat, "CLDN4")
  mal_idx <- which(cls$malignant)
  units_123[[donor]] <- data.frame(
    dataset = "GSE123902",
    unit_id = donor,
    unit_type = "donor",
    tissue = tumor_123$tissue[i],
    n_cells = ncol(mat),
    n_malignant = sum(cls$malignant),
    n_tnk = sum(cls$tnk),
    frac_tnk = mean(cls$tnk),
    mal_CLDN4_pct = if (length(mal_idx)) 100 * mean(cldn4[mal_idx] > 0) else NA_real_,
    mal_CLDN4_mean = if (length(mal_idx)) mean(log1p(cldn4[mal_idx])) else NA_real_,
    eligible = TRUE,
    malig_def = "marker_malig",
    stringsAsFactors = FALSE
  )
  if (length(mal_idx)) {
    pb_123[[donor]] <- Matrix::rowSums(mat[, mal_idx, drop = FALSE])
  }
  keep <- colnames(mat)[qc_pass(mat)]
  keep <- cap_cells(keep)
  meta <- data.frame(
    dataset = "GSE123902",
    unit_id = donor,
    unit_type = "donor",
    cell_class_src = ifelse(cls$malignant[keep], "malignant",
                     ifelse(cls$tnk[keep], "T/NK", "other")),
    author_malignant = FALSE,
    stringsAsFactors = FALSE,
    row.names = paste0("GSE123902_", donor, "_", keep)
  )
  sub <- mat[, keep, drop = FALSE]
  colnames(sub) <- rownames(meta)
  objs[[paste0("GSE123902_", donor)]] <- seurat_from_mat(sub, meta, "GSE123902")
  rm(raw, mat, sub)
  gc(verbose = FALSE)
}
units_123 <- do.call(rbind, units_123)

# ---------------------------------------------------------------------------
# GSE189357 — patient, marker gate, 10x MTX
# ---------------------------------------------------------------------------
say("GSE189357 extract")
d189 <- file.path(GEO, "gse189357")
dir.create(d189, showWarnings = FALSE)
if (!length(list.files(d189, pattern = "_matrix.mtx.gz"))) {
  system(sprintf("tar -C %s -xf %s", d189, file.path(GEO, "GSE189357_RAW.tar")))
}
locked_189 <- read.delim(LOCKED_189357, stringsAsFactors = FALSE)
units_189 <- list()
pb_189 <- list()
for (i in seq_len(nrow(locked_189))) {
  pat <- locked_189$patient[i]
  say("  GSE189357", pat)
  mtx <- list.files(d189, pattern = paste0("_", pat, "_matrix.mtx.gz$"), full.names = TRUE)
  feat <- list.files(d189, pattern = paste0("_", pat, "_features.tsv.gz$"), full.names = TRUE)
  bc <- list.files(d189, pattern = paste0("_", pat, "_barcodes.tsv.gz$"), full.names = TRUE)
  if (!length(mtx) || !length(feat) || !length(bc)) stop("missing 10x files for ", pat)
  td <- file.path(d189, paste0("td_", pat))
  dir.create(td, showWarnings = FALSE)
  file.copy(mtx, file.path(td, "matrix.mtx.gz"), overwrite = TRUE)
  file.copy(feat, file.path(td, "features.tsv.gz"), overwrite = TRUE)
  file.copy(bc, file.path(td, "barcodes.tsv.gz"), overwrite = TRUE)
  mat <- Read10X(td, gene.column = 2)
  if (is.list(mat)) mat <- mat[[1]]
  mat <- upper_rownames(as_dgC(mat))
  cls <- marker_classes(mat)
  cldn4 <- gene_or_zero(mat, "CLDN4")
  mal_idx <- which(cls$malignant)
  units_189[[pat]] <- data.frame(
    dataset = "GSE189357",
    unit_id = pat,
    unit_type = "patient",
    tissue = "TUMOR",
    n_cells = ncol(mat),
    n_malignant = sum(cls$malignant),
    n_tnk = sum(cls$tnk),
    frac_tnk = mean(cls$tnk),
    mal_CLDN4_pct = if (length(mal_idx)) 100 * mean(cldn4[mal_idx] > 0) else NA_real_,
    mal_CLDN4_mean = if (length(mal_idx)) mean(log1p(cldn4[mal_idx])) else NA_real_,
    eligible = TRUE,
    malig_def = "marker_malig",
    stringsAsFactors = FALSE
  )
  if (length(mal_idx)) pb_189[[pat]] <- Matrix::rowSums(mat[, mal_idx, drop = FALSE])
  keep <- colnames(mat)[qc_pass(mat)]
  keep <- cap_cells(keep)
  meta <- data.frame(
    dataset = "GSE189357",
    unit_id = pat,
    unit_type = "patient",
    cell_class_src = ifelse(cls$malignant[keep], "malignant",
                     ifelse(cls$tnk[keep], "T/NK", "other")),
    author_malignant = FALSE,
    stringsAsFactors = FALSE,
    row.names = paste0("GSE189357_", pat, "_", gsub("-", ".", keep))
  )
  sub <- mat[, keep, drop = FALSE]
  colnames(sub) <- rownames(meta)
  objs[[paste0("GSE189357_", pat)]] <- seurat_from_mat(sub, meta, "GSE189357")
  rm(mat, sub)
  gc(verbose = FALSE)
}
units_189 <- do.call(rbind, units_189)

# ---------------------------------------------------------------------------
# GSE205335 — patient, author labels; drop Normal * from T/NK denom
# ---------------------------------------------------------------------------
say("GSE205335 load RDS")
ident <- read.delim(file.path(GEO, "GSE205335_Lung_IO_CellIdentity.txt.gz"), stringsAsFactors = FALSE)
soft_lines <- readLines(file.path(GEO, "GSE205335_family.soft.gz"))
soft <- list()
cur <- list()
for (line in soft_lines) {
  if (startsWith(line, "^SAMPLE")) {
    if (!is.null(cur$title)) soft[[cur$title]] <- cur
    cur <- list()
  } else if (startsWith(line, "!Sample_title")) {
    cur$title <- sub("^!Sample_title = ", "", line)
  } else if (startsWith(line, "!Sample_characteristics_ch1")) {
    val <- sub("^!Sample_characteristics_ch1 = ", "", line)
    if (grepl(": ", val, fixed = TRUE)) {
      kv <- strsplit(val, ": ", fixed = TRUE)[[1]]
      cur[[kv[1]]] <- paste(kv[-1], collapse = ": ")
    }
  }
}
if (!is.null(cur$title)) soft[[cur$title]] <- cur

norm_code <- function(x) {
  x <- toupper(gsub("-", "_", x))
  sub("_[35]P$", "", x)
}
soft_map <- do.call(rbind, lapply(soft, function(s) {
  code <- if (!is.null(s$title) && grepl(" ", s$title)) sub("^\\S+\\s+", "", s$title) else NA
  data.frame(
    title = s$title %||% NA,
    patient = s$patient %||% NA,
    tissue = s$tissue %||% NA,
    recist = s$recist %||% NA,
    code = toupper(gsub("-", "_", code)),
    stringsAsFactors = FALSE
  )
}))
ident$orig_code <- norm_code(ident$orig.ident)
ident <- merge(ident, soft_map, by.x = "orig_code", by.y = "code", all.x = TRUE)
ident$tissue[is.na(ident$tissue)] <- ""
ident$is_normal <- grepl("^Normal", ident$tissue)
ident$author_malignant <- ident$lineage.sub == "Malignant cells"
ident$author_tnk <- ident$lineage.total == "T/NK cells"

locked_205 <- read.delim(LOCKED_205335, stringsAsFactors = FALSE)
keep_pats <- locked_205$patient
ident <- ident[ident$patient %in% keep_pats, ]

rds_plain <- file.path(GEO, "GSE205335_Lung_IO_UMI_matrix.rds")
if (!file.exists(rds_plain)) {
  say("  double-gunzip RDS")
  system(sprintf("gunzip -c %s | gunzip > %s", file.path(GEO, "GSE205335_Lung_IO_UMI_matrix.rds.gz"), rds_plain))
}
mat205 <- readRDS(rds_plain)
if (!inherits(mat205, "dgCMatrix")) mat205 <- as_dgC(mat205)
mat205 <- upper_rownames(mat205)
# barcode match
bc_ident <- ident$barcode
common_bc <- intersect(colnames(mat205), bc_ident)
if (length(common_bc) < 1000) {
  # try without suffix / with underscore swaps
  alt <- gsub("_", "-", colnames(mat205))
  names(alt) <- colnames(mat205)
  common_bc <- colnames(mat205)[alt %in% bc_ident]
  if (length(common_bc)) {
    map <- ident
    rownames(map) <- map$barcode
    ident2 <- map[alt[common_bc], ]
    ident2$mat_barcode <- common_bc
  } else {
    stop("GSE205335 barcode mismatch: ident ", length(bc_ident), " matrix ", ncol(mat205),
         " example matrix ", colnames(mat205)[1], " example ident ", bc_ident[1])
  }
} else {
  ident2 <- ident
  rownames(ident2) <- ident2$barcode
  ident2 <- ident2[common_bc, ]
  ident2$mat_barcode <- common_bc
}
say("  GSE205335 matched cells", nrow(ident2))

units_205 <- list()
pb_205 <- list()
for (pat in keep_pats) {
  rows <- ident2[ident2$patient == pat, , drop = FALSE]
  if (!nrow(rows)) next
  tumor_rows <- rows[!rows$is_normal, , drop = FALSE]
  if (!nrow(tumor_rows)) tumor_rows <- rows
  n_cells <- nrow(tumor_rows)
  n_mal <- sum(tumor_rows$author_malignant %in% c(TRUE, "TRUE"), na.rm = TRUE)
  n_tnk <- sum(tumor_rows$author_tnk %in% c(TRUE, "TRUE"), na.rm = TRUE)
  mal_bc <- tumor_rows$mat_barcode[tumor_rows$author_malignant %in% c(TRUE, "TRUE")]
  mal_bc <- intersect(mal_bc, colnames(mat205))
  cldn4 <- if (length(mal_bc) && "CLDN4" %in% rownames(mat205)) as.numeric(mat205["CLDN4", mal_bc]) else numeric()
  units_205[[pat]] <- data.frame(
    dataset = "GSE205335",
    unit_id = pat,
    unit_type = "patient",
    tissue = paste(sort(unique(tumor_rows$tissue)), collapse = ","),
    n_cells = n_cells,
    n_malignant = n_mal,
    n_tnk = n_tnk,
    frac_tnk = if (n_cells) n_tnk / n_cells else NA_real_,
    mal_CLDN4_pct = if (length(cldn4)) 100 * mean(cldn4 > 0) else NA_real_,
    mal_CLDN4_mean = if (length(cldn4)) mean(log1p(cldn4)) else NA_real_,
    eligible = TRUE,
    malig_def = "author_malig",
    stringsAsFactors = FALSE
  )
  if (length(mal_bc)) pb_205[[pat]] <- Matrix::rowSums(mat205[, mal_bc, drop = FALSE])
  keep_bc <- tumor_rows$mat_barcode
  keep_bc <- keep_bc[qc_pass(mat205[, keep_bc, drop = FALSE])]
  keep_bc <- cap_cells(keep_bc)
  sub_meta <- tumor_rows
  rownames(sub_meta) <- sub_meta$mat_barcode
  sub_meta <- sub_meta[keep_bc, , drop = FALSE]
  meta <- data.frame(
    dataset = "GSE205335",
    unit_id = pat,
    unit_type = "patient",
    cell_class_src = ifelse(sub_meta$author_malignant, "malignant",
                     ifelse(sub_meta$author_tnk, "T/NK", "other")),
    author_malignant = sub_meta$author_malignant,
    stringsAsFactors = FALSE,
    row.names = paste0("GSE205335_", pat, "_", keep_bc)
  )
  sub <- mat205[, keep_bc, drop = FALSE]
  colnames(sub) <- rownames(meta)
  objs[[paste0("GSE205335_", pat)]] <- seurat_from_mat(sub, meta, "GSE205335")
}
rm(mat205)
gc(verbose = FALSE)
units_205 <- do.call(rbind, units_205)

# ---------------------------------------------------------------------------
# GSE131907 — sample unit; Python stream for atlas + malignant sums
# ---------------------------------------------------------------------------
ext <- file.path(GEO, "gse131907_extract")
if (!file.exists(file.path(ext, "atlas_10x", "matrix.mtx.gz"))) {
  say("GSE131907 python stream (I/O only)")
  rc <- system2("python3", c(file.path(ROOT, "prep_gse131907.py")), stdout = "", stderr = "")
  if (rc != 0) stop("prep_gse131907.py failed")
}
say("GSE131907 read extract")
u131 <- read.delim(file.path(ext, "units.tsv"), stringsAsFactors = FALSE)
u131 <- u131[u131$eligible %in% c("True", "TRUE", TRUE), ]
units_131 <- data.frame(
  dataset = "GSE131907",
  unit_id = u131$unit_id,
  unit_type = "sample",
  tissue = u131$origin,
  n_cells = u131$n_cells,
  n_malignant = u131$n_malignant,
  n_tnk = u131$n_tnk,
  frac_tnk = u131$frac_tnk,
  mal_CLDN4_pct = as.numeric(u131$mal_CLDN4_pct),
  mal_CLDN4_mean = as.numeric(u131$mal_CLDN4_mean),
  eligible = TRUE,
  malig_def = "author_malig",
  stringsAsFactors = FALSE
)
pb131_mat <- as.matrix(read.delim(file.path(ext, "malig_pseudobulk.tsv.gz"), row.names = 1, check.names = FALSE))
rownames(pb131_mat) <- toupper(rownames(pb131_mat))
pb_131 <- lapply(colnames(pb131_mat), function(s) {
  v <- pb131_mat[, s]
  names(v) <- rownames(pb131_mat)
  v
})
names(pb_131) <- colnames(pb131_mat)

mat131 <- Read10X(file.path(ext, "atlas_10x"), gene.column = 2)
if (is.list(mat131)) mat131 <- mat131[[1]]
mat131 <- upper_rownames(as_dgC(mat131))
meta131 <- read.delim(file.path(ext, "atlas_meta.tsv"), stringsAsFactors = FALSE)
stopifnot(ncol(mat131) == nrow(meta131))
# Read10X may strip / alter barcodes; match by order
colnames(mat131) <- paste0("GSE131907_", meta131$unit_id, "_", meta131$barcode)
meta <- data.frame(
  dataset = "GSE131907",
  unit_id = meta131$unit_id,
  unit_type = "sample",
  cell_class_src = ifelse(meta131$author_malignant, "malignant",
                   ifelse(meta131$author_tnk, "T/NK", "other")),
  author_malignant = meta131$author_malignant,
  stringsAsFactors = FALSE,
  row.names = colnames(mat131)
)
objs[["GSE131907"]] <- seurat_from_mat(mat131, meta, "GSE131907")
rm(mat131, pb131_mat)
gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# Patient / donor / sample table + tests
# ---------------------------------------------------------------------------
patients <- rbind(units_123, units_131, units_205, units_189)
patients$cldn4_quartile <- NA_character_
for (ds in unique(patients$dataset)) {
  ix <- patients$dataset == ds
  patients$cldn4_quartile[ix] <- within_quartile(patients$mal_CLDN4_pct[ix])
}
write.table(patients, file.path(OUT_TAB, "patient_units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Spearman per cohort + DL meta
cohorts <- c("GSE123902", "GSE131907", "GSE205335", "GSE189357")
tnk_single <- do.call(rbind, lapply(cohorts, function(ds) {
  d <- patients[patients$dataset == ds, ]
  ct <- suppressWarnings(stats::cor.test(d$mal_CLDN4_pct, d$frac_tnk, method = "spearman", exact = FALSE))
  data.frame(dataset = ds, n = nrow(d), rho = unname(ct$estimate), p = ct$p.value, stringsAsFactors = FALSE)
}))
dl <- dl_spearman(tnk_single$rho, tnk_single$n)
say("T/NK DL rho", signif(dl$rho, 3), "p", signif(dl$p, 3), "N", dl$N)

q4 <- patients$cldn4_quartile == "Q4"
q1 <- patients$cldn4_quartile == "Q1"
rb <- rank_biserial(patients$frac_tnk[q4], patients$frac_tnk[q1])
say("T/NK stacked Q4vsQ1 r", signif(rb$r, 3), "p", signif(rb$p, 3), "n", rb$n_q1, rb$n_q4)

# Malignant family scores from full-data UMI-sum (not the 350-cell cap)
align_pb <- function(pb_list) {
  genes <- Reduce(union, lapply(pb_list, names))
  mat <- matrix(0, nrow = length(genes), ncol = length(pb_list), dimnames = list(genes, names(pb_list)))
  for (nm in names(pb_list)) {
    v <- pb_list[[nm]]
    names(v) <- toupper(names(v))
    common <- intersect(names(v), genes)
    mat[common, nm] <- as.numeric(v[common])
  }
  mat
}
pb_all <- c(
  setNames(pb_123, paste0("GSE123902|", names(pb_123))),
  setNames(pb_131, paste0("GSE131907|", names(pb_131))),
  setNames(pb_205, paste0("GSE205335|", names(pb_205))),
  setNames(pb_189, paste0("GSE189357|", names(pb_189)))
)
pb_mat <- align_pb(pb_all)
# drop sparse noisy units from DE (n_mal < 30), matching prior P4001 exclusion
de_units <- patients[is.finite(patients$n_malignant) & patients$n_malignant >= 30, ]
de_ids <- paste(de_units$dataset, de_units$unit_id, sep = "|")
pb_mat <- pb_mat[, intersect(colnames(pb_mat), de_ids), drop = FALSE]
de_units <- de_units[match(colnames(pb_mat), paste(de_units$dataset, de_units$unit_id, sep = "|")), ]
sf <- tmm_factors(pb_mat)
lib <- colSums(pb_mat)
cpm <- sweep(pb_mat, 2, (lib / sf) / 1e6, "/")
logcpm <- log2(cpm + 1)
score_fun <- function(genes) {
  g <- family_present(genes, rownames(logcpm))
  if (!length(g)) return(rep(NA_real_, ncol(logcpm)))
  colMeans(logcpm[g, , drop = FALSE])
}
de_units$ifn_score <- score_fun(sets$IFN)
de_units$mhc_score <- score_fun(sets$MHC)
de_units$tj_score <- score_fun(sets$TJ)
# Keep the within-cohort quartiles from the full unit table. Do not recut.
# attach scores onto the full patient table where possible
patients$ifn_score <- de_units$ifn_score[match(paste(patients$dataset, patients$unit_id, sep = "|"),
                                              paste(de_units$dataset, de_units$unit_id, sep = "|"))]
patients$mhc_score <- de_units$mhc_score[match(paste(patients$dataset, patients$unit_id, sep = "|"),
                                              paste(de_units$dataset, de_units$unit_id, sep = "|"))]
patients$tj_score <- de_units$tj_score[match(paste(patients$dataset, patients$unit_id, sep = "|"),
                                            paste(de_units$dataset, de_units$unit_id, sep = "|"))]
write.table(patients, file.path(OUT_TAB, "patient_units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(de_units, file.path(OUT_TAB, "malignant_pseudobulk_units.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

family_ols <- function(score, d) {
  d2 <- d[d$cldn4_quartile %in% c("Q1", "Q4") & is.finite(score), ]
  sc <- score[d$cldn4_quartile %in% c("Q1", "Q4") & is.finite(score)]
  d2$CLDN4_Q4 <- as.integer(d2$cldn4_quartile == "Q4")
  d2$cohort <- factor(d2$dataset, levels = cohorts)
  fit <- stats::lm(sc ~ cohort + CLDN4_Q4, data = d2)
  sm <- summary(fit)
  coef <- sm$coefficients
  data.frame(
    n_q1 = sum(d2$cldn4_quartile == "Q1"),
    n_q4 = sum(d2$cldn4_quartile == "Q4"),
    logFC = unname(coef["CLDN4_Q4", "Estimate"]),
    p = unname(coef["CLDN4_Q4", "Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
}
fam_tab <- rbind(
  cbind(family = "IFN", family_ols(de_units$ifn_score, de_units)),
  cbind(family = "MHC-I/APM", family_ols(de_units$mhc_score, de_units)),
  cbind(family = "TJ", family_ols(de_units$tj_score, de_units))
)
fam_tab$n_genes <- c(
  length(family_present(sets$IFN, rownames(logcpm))),
  length(family_present(sets$MHC, rownames(logcpm))),
  length(family_present(sets$TJ, rownames(logcpm)))
)
write.table(fam_tab, file.path(OUT_TAB, "family_q4q1_ols.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tnk_single, file.path(OUT_TAB, "tnk_spearman_singles.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# Seurat merge + Harmony + UMAP
# ---------------------------------------------------------------------------
say("merge", length(objs), "unit objects")
# restrict to shared genes before merge to cut RAM
common <- Reduce(intersect, lapply(objs, rownames))
say("shared genes", length(common))
objs <- lapply(objs, function(o) o[common, ])
merged <- merge(x = objs[[1]], y = objs[-1])
rm(objs)
gc(verbose = FALSE)
if (packageVersion("Seurat") >= "5.0.0") {
  merged <- JoinLayers(merged)
}
say("merged cells", ncol(merged), "genes", nrow(merged))
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
say("RunHarmony dataset")
merged <- RunHarmony(merged, group.by.vars = "dataset", reduction.use = "pca",
                     dims.use = PCA_DIMS, reduction.save = "harmony", verbose = TRUE)
merged <- RunUMAP(merged, reduction = "harmony", dims = PCA_DIMS, verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "harmony", dims = PCA_DIMS, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.6, verbose = FALSE)

# cluster-level cell class
DefaultAssay(merged) <- "RNA"
markers <- list(
  malignant = c("EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4"),
  T = c("CD3D", "CD3E"),
  NK = c("NKG7", "GNLY", "KLRD1"),
  myeloid = c("CD68", "LYZ"),
  B = c("MS4A1", "CD79A"),
  endothelial = c("PECAM1", "VWF"),
  fibroblast = c("COL1A1")
)
Idents(merged) <- "seurat_clusters"
clust <- sort(unique(as.character(merged$seurat_clusters)))
ann <- do.call(rbind, lapply(clust, function(cl) {
  cells <- WhichCells(merged, idents = cl)
  md <- merged[[]][cells, , drop = FALSE]
  expr <- GetAssayData(merged, layer = "data")[, cells, drop = FALSE]
  scores <- sapply(markers, function(gs) {
    gs <- intersect(gs, rownames(expr))
    if (!length(gs)) return(0)
    mean(Matrix::colMeans(expr[gs, , drop = FALSE]))
  })
  author_frac <- mean(md$author_malignant %in% c(TRUE, "TRUE"))
  ord <- order(scores, decreasing = TRUE)
  top <- names(scores)[ord[1]]
  second <- names(scores)[ord[2]]
  rule <- "marker"
  label <- top
  if (is.finite(author_frac) && author_frac >= 0.5) {
    label <- "malignant"
    rule <- "author_malignant_frac>=0.50"
  } else if (scores[top] < 0.15) {
    label <- "other"
    rule <- "top<0.15"
  } else if (scores[second] >= 0.15 && abs(scores[top] - scores[second]) < 0.05) {
    label <- "other"
    rule <- "mixed"
  }
  coarse <- if (label %in% c("endothelial", "fibroblast")) "other" else label
  if (label == "malignant") coarse <- "malignant"
  data.frame(
    cluster = cl, n = length(cells), label = coarse, fine = label,
    author_malignant_frac = author_frac, top_score = unname(scores[top]),
    rule = rule, stringsAsFactors = FALSE
  )
}))
write.table(ann, file.path(OUT_TAB, "cluster_annotation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
map <- setNames(ann$label, ann$cluster)
merged$cell_class <- unname(map[as.character(merged$seurat_clusters)])

# slim object: drop scale.data to keep RDS smaller; keep counts/data + reductions
if ("scale.data" %in% Layers(merged[["RNA"]])) {
  merged[["RNA"]]$scale.data <- NULL
}
saveRDS(merged, file.path(OUT_OBJ, "seurat_harmony_integrated.rds"))
saveRDS(
  list(
    umap = Embeddings(merged, "umap"),
    harmony = Embeddings(merged, "harmony"),
    meta = merged[[]]
  ),
  file.path(OUT_OBJ, "seurat_harmony_embeddings.rds")
)
say("saved Seurat object", ncol(merged), "cells")

# DimPlots
theme_set(theme_bw(base_size = 11))
p_ds <- DimPlot(merged, reduction = "umap", group.by = "dataset", pt.size = 0.2) +
  ggtitle("Harmony UMAP by dataset") +
  theme(legend.position = "right")
p_cl <- DimPlot(merged, reduction = "umap", group.by = "cell_class", pt.size = 0.2) +
  ggtitle("Harmony UMAP by cell class") +
  theme(legend.position = "right")
ggsave(file.path(OUT_FIG, "DimPlot_dataset.png"), p_ds, width = 7.2, height = 5.6, dpi = 300)
ggsave(file.path(OUT_FIG, "DimPlot_dataset.pdf"), p_ds, width = 7.2, height = 5.6)
ggsave(file.path(OUT_FIG, "DimPlot_cell_class.png"), p_cl, width = 7.2, height = 5.6, dpi = 300)
ggsave(file.path(OUT_FIG, "DimPlot_cell_class.pdf"), p_cl, width = 7.2, height = 5.6)
p_both <- p_ds + p_cl
ggsave(file.path(OUT_FIG, "DimPlot_dataset_and_cell_class.png"), p_both, width = 12.5, height = 5.6, dpi = 300)
ggsave(file.path(OUT_FIG, "DimPlot_dataset_and_cell_class.pdf"), p_both, width = 12.5, height = 5.6)

inv <- as.data.frame(table(merged$dataset), stringsAsFactors = FALSE)
names(inv) <- c("dataset", "n_cells_used")
inv$n_units <- as.integer(tapply(merged$unit_id, merged$dataset, function(x) length(unique(x)))[inv$dataset])
write.table(inv, file.path(OUT_TAB, "inventory_used_cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# FINDING.md
# ---------------------------------------------------------------------------
fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 1e-4) format(p, digits = 3, scientific = TRUE) else sprintf("%.4g", p)
}
n_used <- ncol(merged)
n_units <- nrow(patients)
class_tab <- table(merged$cell_class)
finding <- paste0(
  "# Seurat concordant-4 CLDN4-only\n\n",
  "ADDITIVE. **CLDN4-only.** No TACSTD2∩CLDN4 dual-high. Not a mega-merge.\n",
  "The four sets that already point the same way: **GSE123902 + GSE131907 +\n",
  "GSE205335 + GSE189357**. Not GSE148071, GSE127465, GSE207422, GSE154826,\n",
  "or CD45+/T-only extracts.\n\n",
  "Thesis (already correct; not re-derived): CLDN4-high malignant cells have\n",
  "lower own IFN/MHC-I, and patients have lower T/NK. CLDN4-low/KD opens IFN/MHC.\n",
  "Matching extras here are IFN/MHC **DOWN** in CLDN4-high.\n\n",
  "Primary stack is **R + Seurat ", as.character(packageVersion("Seurat")),
  " + Harmony ", as.character(packageVersion("harmony")),
  "** (`RunHarmony` on the Seurat object, `group.by.vars = dataset`).\n",
  "Python was used only to stream the GSE131907 genes×cells text into a sparse\n",
  "atlas subset. Patient-level numbers use **full-unit** counts (not the UMAP cap).\n\n",
  "Patient / donor / sample is the unit. Do not quote cell counts as n.\n",
  "GSE123902 = donor. GSE131907 = sample (tumor-bearing). GSE205335 = patient\n",
  "(RECIST not required). GSE189357 = patient. p-values are descriptive.\n\n",
  "## Honest n\n\n",
  "- **n_units = ", n_units, "** (",
  sum(patients$dataset == "GSE123902"), " donors + ",
  sum(patients$dataset == "GSE131907"), " samples + ",
  sum(patients$dataset == "GSE205335"), " patients + ",
  sum(patients$dataset == "GSE189357"), " patients).\n",
  "- **n_cells in the Seurat object (after QC + ≤", CAP, "/unit cap) = ", n_used, "**.\n",
  "- Do not replace the unit n with this cell count.\n\n",
  paste0(apply(inv, 1, function(r) paste0("- ", r[["dataset"]], ": ", r[["n_units"]],
         " units, ", r[["n_cells_used"]], " cells used")), collapse = "\n"), "\n\n",
  "UMAP cell class (used cells, not the test n): ",
  paste(sprintf("%s=%s", names(class_tab), as.integer(class_tab)), collapse = "; "), ".\n\n",
  "## 1. Patient-level malignant CLDN4 vs T/NK fraction\n\n",
  "Primary score = malignant CLDN4 **%pos** (full unit, not the cap).\n",
  "Pooling = DerSimonian–Laird on Fisher-z of the four cohort Spearmans.\n",
  "Q4 vs Q1 = **within-cohort quartiles stacked**, then Mann–Whitney on T/NK\n",
  "(rank-biserial r). Marker gate (GSE123902, GSE189357): ",
  "(EPCAM|KRT8|KRT18|KRT19)>0 AND PTPRC==0 vs (CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1)>0.\n",
  "Author labels (GSE131907, GSE205335). `frac_tnk = n_tnk / n_cells`.\n\n",
  "| score | k | N | ρ (p, I², 95% CI) | stacked Q4 vs Q1 r (n_Q1/n_Q4, p) |\n",
  "|---|---:|---:|---|---|\n",
  sprintf("| %%pos | %d | %d | %.3f (%s, I²=%.1f%%, %.3f to %.3f) | %.3f (%d/%d, %s) |\n",
          dl$k, dl$N, dl$rho, fmt_p(dl$p), 100 * dl$I2, dl$ci_lo, dl$ci_hi,
          rb$r, rb$n_q1, rb$n_q4, fmt_p(rb$p)),
  "\n### Singles (context; not the new number)\n\n",
  "| cohort | unit | n | ρ | p |\n|---|---|---:|---:|---:|\n",
  paste0(apply(tnk_single, 1, function(r) {
    sprintf("| %s | %s | %s | %.3f | %s |",
            r[["dataset"]],
            c(GSE123902 = "donor", GSE131907 = "sample", GSE205335 = "patient", GSE189357 = "patient")[r[["dataset"]]],
            r[["n"]], as.numeric(r[["rho"]]), fmt_p(as.numeric(r[["p"]])))
  }), collapse = "\n"), "\n\n",
  "## 2. Malignant Q4 vs Q1 IFN / MHC / TJ\n\n",
  "Patient-pseudobulk OLS on log2(TMM-CPM+1) family means, `~ cohort + CLDN4_Q4`.\n",
  "Within-cohort %pos quartiles. Units with n_malignant < 30 are out of this DE\n",
  "(P4001-style noise gate). Positive logFC = higher in CLDN4-high.\n",
  "IFN = Hallmark IFNα ∪ IFNγ. MHC = custom MHC-I/APM. TJ = KEGG ∪ GOBP organization\n",
  "plus CDH1/VIM/ZEB1; **CLDN4 held out**.\n\n",
  "| family | n_Q1 / n_Q4 | n_genes | logFC | p |\n",
  "|---|---|---:|---:|---:|\n",
  paste0(apply(fam_tab, 1, function(r) {
    sprintf("| %s | %s / %s | %s | %.3f | %s |",
            r[["family"]], r[["n_q1"]], r[["n_q4"]], r[["n_genes"]],
            as.numeric(r[["logFC"]]), fmt_p(as.numeric(r[["p"]])))
  }), collapse = "\n"), "\n\n",
  "Expected under the thesis: IFN down, MHC-I/APM down, TJ up or held.\n",
  "Do not quote N=", n_units, " as the DE n. Do not quote ", n_used, " cells as n.\n\n",
  "## 3. Seurat / Harmony UMAP\n\n",
  "- QC: n_genes ≥ ", QC_MIN_GENES, ", n_UMI ≥ ", QC_MIN_UMI, ", mitochondrial % < ", QC_MAX_MT, ".\n",
  "- Cap: ≤", CAP, " cells / unit (memory).\n",
  "- Inner gene join → NormalizeData → VST 2000 HVG → ScaleData → PCA 30.\n",
  "- Harmony: `group.by.vars = dataset` only (sample was not a second key).\n",
  "- Neighbors / UMAP / FindClusters resolution 0.6 on the Harmony embedding.\n",
  "- Cluster labels: author-malignant fraction ≥ 0.50 → malignant; else highest\n",
  "  mean log-normalized marker score (EPCAM/KRT/CLDN4, CD3, NK, myeloid, B).\n",
  "  Endothelium / fibroblast / mixed / low-score → **other**.\n\n",
  "DimPlots: `results/figures/DimPlot_dataset.png`, `DimPlot_cell_class.png`.\n",
  "Seurat object: `results/objects/seurat_harmony_integrated.rds`.\n",
  "Patient table: `results/tables/patient_units.tsv`.\n\n",
  "## Reproduce\n\n",
  "```\nRscript methods/seurat_concordant4_cldn4/analyze.R\n```\n"
)
writeLines(finding, file.path(ROOT, "FINDING.md"))
say("wrote FINDING.md")
say("DONE units", n_units, "cells", n_used)
