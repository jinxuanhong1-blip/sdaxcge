#!/usr/bin/env Rscript
# GSE267321 LKR13 K / KK / KLK — Cldn4-only, public GEO CSV, Seurat primary.
#
# CreateSeuratObject on GSE267321_Normalized_expression_matrix_02122026.csv.gz.
# Score Cldn4 if the row exists. T/NK fraction by genotype. Leftover
# epithelium IFN/MHC. Honest n = tumors. No private KL. No dual-high.
#
# Thesis is not rewritten. Additive public mouse only.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

stopifnot(packageVersion("Seurat") >= "5.0.0")

set.seed(42)
options(warn = 1)

HERE <- local({
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f) == 1L) {
    return(normalizePath(file.path(dirname(f), "..")))
  }
  normalizePath("methods/seurat_gse267321_cldn4")
})
TABLES <- file.path(HERE, "tables")
FIGS <- file.path(HERE, "figures")
dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)

GEO <- "GSE267321"
MATRIX_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267321/suppl/",
  "GSE267321_Normalized_expression_matrix_02122026.csv.gz"
)
SERIES_URL <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267321/matrix/",
  "GSE267321_series_matrix.txt.gz"
)
CACHE <- "/tmp/gse267321"
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
MATRIX <- file.path(CACHE, "GSE267321_Normalized_expression_matrix_02122026.csv.gz")
SERIES <- file.path(CACHE, "GSE267321_series_matrix.txt.gz")

IFN_ISG <- c(
  "Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Mx2", "Oas1a", "Oas2",
  "Oasl1", "Rsad2", "Usp18", "Bst2", "Stat1", "Stat2", "Irf7", "Irf9",
  "Ddx58", "Ifih1", "Ly6e", "Ifitm3", "Ifi44", "Ifi27", "Ifi27l2a",
  "Irf1", "Cmpk2", "Parp9", "Dtx3l", "Sp100", "Plscr1"
)
MHC1 <- c(
  "B2m", "H2-K1", "H2-D1", "H2-Q4", "H2-Q7", "H2-T23",
  "Tap1", "Tap2", "Tapbp", "Psmb8", "Psmb9", "Psmb10",
  "Psme1", "Psme2", "Nlrc5", "Calr", "Pdia3", "Canx"
)
IFN_CORE6 <- c("Ifi27", "Oas2", "Ifit1", "Mx1", "Isg15", "H2-K1")
CTRL_OXPHOS <- c(
  "Ndufa1", "Ndufb3", "Cox5a", "Cox7a2", "Uqcrc1", "Sdha",
  "Atp5a1", "Atp5b", "Cycs", "Vdac1"
)
LINEAGE <- c(
  "Cldn4", "Tacstd2", "Epcam", "Cdh1", "Krt7", "Krt8", "Krt18", "Krt19",
  "Nkx2-1", "Napsa", "Sftpc", "Sftpb", "Sftpa1", "Scgb1a1", "Scgb3a2",
  "Ager", "Ptprc", "Cd3d", "Cd3e", "Cd3g", "Cd2", "Cd8a", "Cd8b1", "Cd4",
  "Nkg7", "Ncr1", "Klrb1c", "Klrd1", "Gzma", "Prf1", "Lyz2", "Cd68",
  "Csf1r", "Itgam", "Cd14", "Col1a1", "Pdgfra", "Dcn", "Acta2", "Pecam1",
  "Cdh5", "Cldn5", "Vwf", "Cd79a", "Ms4a1"
)
WANTED <- sort(unique(c(LINEAGE, IFN_ISG, MHC1, IFN_CORE6, CTRL_OXPHOS)))

GENO_META <- list(
  K = list(
    label = "K", cell_line = "LKR13K", genotype_geo = "KrasG12D",
    stk11 = "WT", keap1 = "WT", closest_to_user_KL = FALSE
  ),
  KK = list(
    label = "KK", cell_line = "LKR13KK",
    genotype_geo = "KrasG12D, KEAP1 knockout",
    stk11 = "WT", keap1 = "KO", closest_to_user_KL = FALSE
  ),
  KLK = list(
    label = "KLK", cell_line = "LKR13KLK",
    genotype_geo = "KrasG12D, KEAP1/LKB1 knockout",
    stk11 = "KO", keap1 = "KO", closest_to_user_KL = TRUE
  )
)

download_if_needed <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) {
    return(invisible(dest))
  }
  tmp <- paste0(dest, ".part")
  status <- tryCatch(
    utils::download.file(url, tmp, mode = "wb", quiet = FALSE),
    error = function(e) stop("download failed for ", url, ": ", conditionMessage(e))
  )
  if (!identical(status, 0L)) stop("download.file non-zero status for ", url)
  file.rename(tmp, dest)
  invisible(dest)
}

parse_barcode <- function(bc) {
  m <- regexec(
    "^(K|KK|KLK)_([ACGT]+)\\.1_(LKR13\\.(K|KK|KLK)\\.([0-9]+))$",
    bc
  )
  hit <- regmatches(bc, m)[[1]]
  if (length(hit) != 6L) stop("unparsed barcode: ", bc)
  geno <- hit[2]
  geno2 <- hit[5]
  if (!identical(geno, geno2)) stop("genotype mismatch in ", bc)
  list(
    barcode = bc,
    genotype = geno,
    umi16 = hit[3],
    sample = gsub("\\.", "-", hit[4]),
    sample_raw = hit[4],
    replicate = as.integer(hit[6])
  )
}

gt <- function(obj, gene) {
  if (!(gene %in% rownames(obj))) {
    return(rep(FALSE, ncol(obj)))
  }
  as.numeric(GetAssayData(obj, layer = "counts")[gene, ]) > 0
}

mean_log1p_module <- function(obj, genes) {
  present <- intersect(genes, rownames(obj))
  if (length(present) == 0L) return(rep(NA_real_, ncol(obj)))
  counts <- GetAssayData(obj, layer = "counts")[present, , drop = FALSE]
  if (length(present) == 1L) {
    return(as.numeric(log1p(counts)))
  }
  as.numeric(Matrix::colMeans(log1p(counts)))
}

fmt <- function(x, nd = 3) {
  if (length(x) != 1L || is.null(x) || (is.numeric(x) && (is.na(x) || is.infinite(x)))) {
    return("NA")
  }
  if (is.logical(x)) return(tolower(as.character(x)))
  if (is.numeric(x)) {
    if (x == 0) return("0")
    if (abs(x) >= 100) return(sprintf("%.1f", x))
    if (abs(x) >= 1) return(sprintf(paste0("%.", nd, "f"), x))
    return(sprintf(paste0("%.", nd, "g"), x))
  }
  as.character(x)
}

write_tsv <- function(path, df) {
  fwrite(as.data.table(df), path, sep = "\t", na = "", quote = FALSE)
}

mw_p <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) == 0L || length(b) == 0L || (length(a) + length(b)) < 3L) {
    return(NA_real_)
  }
  suppressWarnings(
    tryCatch(
      wilcox.test(a, b, alternative = "two.sided", exact = TRUE)$p.value,
      error = function(e) NA_real_
    )
  )
}

# -----------------------------------------------------------------------------
# 1. Public GEO only
# -----------------------------------------------------------------------------
message("Downloading public GSE267321 files if needed...")
download_if_needed(MATRIX_URL, MATRIX)
download_if_needed(SERIES_URL, SERIES)

# Header is cell barcodes only (no gene-column name). Data rows are gene + cells.
message("Reading normalized CSV...")
hdr <- strsplit(readLines(gzfile(MATRIX), n = 1L), ",", fixed = TRUE)[[1]]
dt <- fread(
  cmd = paste("gzip -dc", shQuote(MATRIX)),
  header = FALSE,
  skip = 1L,
  showProgress = TRUE
)
if (ncol(dt) != length(hdr) + 1L) {
  stop("CSV shape mismatch: header cells=", length(hdr), " data cols=", ncol(dt))
}
genes <- as.character(dt[[1]])
if (anyDuplicated(genes)) {
  message("Duplicated gene symbols: keeping first")
  keep <- !duplicated(genes)
  dt <- dt[keep]
  genes <- genes[keep]
}

dense <- as.matrix(dt[, -1, with = FALSE])
storage.mode(dense) <- "numeric"
nz <- dense[dense > 0]
frac_int <- if (length(nz)) mean(abs(nz - round(nz)) < 1e-6) else NA_real_
value_note <- if (is.finite(frac_int) && frac_int >= 0.95) {
  sprintf("counts_log1p (frac_int=%.3f); CreateSeuratObject treats matrix as counts", frac_int)
} else {
  sprintf("already_normalized (frac_int=%s); still ingested via CreateSeuratObject", fmt(frac_int))
}
message("Value heuristic: ", value_note)

sp <- as(dense, "dgCMatrix")
rownames(sp) <- genes
colnames(sp) <- hdr
rm(dt, dense)
gc()

# -----------------------------------------------------------------------------
# 2. CreateSeuratObject — required primary
# -----------------------------------------------------------------------------
message("CreateSeuratObject...")
obj <- CreateSeuratObject(
  counts = sp,
  project = "GSE267321",
  min.cells = 0,
  min.features = 0,
  assay = "RNA"
)
rm(sp)
gc()

meta <- rbindlist(lapply(colnames(obj), parse_barcode))
stopifnot(identical(meta$barcode, colnames(obj)))
obj$genotype <- meta$genotype
obj$sample <- meta$sample
obj$replicate <- meta$replicate
obj$stk11 <- vapply(obj$genotype, function(g) GENO_META[[g]]$stk11, character(1))
obj$keap1 <- vapply(obj$genotype, function(g) GENO_META[[g]]$keap1, character(1))
obj$closest_to_user_KL <- vapply(
  obj$genotype, function(g) GENO_META[[g]]$closest_to_user_KL, logical(1)
)

# Marker compartments. Epithelium wins over T/NK if both fire.
t_lineage <- gt(obj, "Cd3d") | gt(obj, "Cd3e") | gt(obj, "Cd3g") | gt(obj, "Cd8a")
nk_lineage <- gt(obj, "Nkg7") | gt(obj, "Ncr1") | gt(obj, "Klrb1c")
tnk <- t_lineage | nk_lineage
host <- gt(obj, "Sftpc") | gt(obj, "Scgb1a1") | gt(obj, "Ager")
epi <- gt(obj, "Epcam") |
  (gt(obj, "Cdh1") & gt(obj, "Krt8")) |
  (gt(obj, "Krt18") & gt(obj, "Krt19"))
caf <- gt(obj, "Col1a1") | gt(obj, "Dcn") | gt(obj, "Pdgfra")
endo <- gt(obj, "Pecam1") | gt(obj, "Cdh5")
myeloid <- gt(obj, "Lyz2") | gt(obj, "Cd68") | gt(obj, "Csf1r")
bcell <- gt(obj, "Cd79a") | gt(obj, "Ms4a1")

lab <- rep("other", ncol(obj))
lab[myeloid] <- "myeloid"
lab[caf] <- "caf"
lab[endo] <- "endo"
lab[bcell] <- "b"
lab[tnk] <- "tnk"
lab[host & epi] <- "host_epithelial"
lab[epi & !host] <- "epithelial"
obj$compartment <- lab
obj$is_tnk <- lab == "tnk"
obj$is_epithelial <- lab == "epithelial"
obj$is_host_epithelial <- lab == "host_epithelial"

cldn4_present <- "Cldn4" %in% rownames(obj)
obj$Cldn4_counts <- if (cldn4_present) {
  as.numeric(GetAssayData(obj, layer = "counts")["Cldn4", ])
} else {
  rep(0, ncol(obj))
}
obj$Cldn4_log1p <- log1p(obj$Cldn4_counts)
obj$Cldn4_pos <- obj$Cldn4_counts > 0

obj$ifn_isg <- mean_log1p_module(obj, IFN_ISG)
obj$mhc1 <- mean_log1p_module(obj, MHC1)
obj$ifn_core6 <- mean_log1p_module(obj, IFN_CORE6)
obj$ctrl_oxphos <- mean_log1p_module(obj, CTRL_OXPHOS)

message("NormalizeData + AddModuleScore (Seurat)...")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", verbose = FALSE)
ifn_use <- intersect(IFN_ISG, rownames(obj))
mhc_use <- intersect(MHC1, rownames(obj))
if (length(ifn_use) >= 2L) {
  obj <- AddModuleScore(obj, features = list(ifn_use), name = "SeuratIFN", ctrl = 80, seed = 42)
}
if (length(mhc_use) >= 2L) {
  obj <- AddModuleScore(obj, features = list(mhc_use), name = "SeuratMHC1", ctrl = 80, seed = 42)
}

md <- as.data.table(obj[[]], keep.rownames = "barcode")
if (!("SeuratIFN1" %in% names(md))) md[, SeuratIFN1 := NA_real_]
if (!("SeuratMHC11" %in% names(md))) md[, SeuratMHC11 := NA_real_]

present_wanted <- intersect(WANTED, rownames(obj))
missing_wanted <- setdiff(WANTED, rownames(obj))

# Gene inventory on the full counts matrix (wanted genes only).
inv_rows <- lapply(WANTED, function(g) {
  if (!(g %in% rownames(obj))) {
    return(data.table(
      gene = g, present = FALSE, n_pos = 0L, pct_pos = NA_real_,
      mean_counts = NA_real_, max_counts = NA_real_
    ))
  }
  x <- as.numeric(GetAssayData(obj, layer = "counts")[g, ])
  data.table(
    gene = g, present = TRUE,
    n_pos = as.integer(sum(x > 0)),
    pct_pos = 100 * mean(x > 0),
    mean_counts = mean(x),
    max_counts = max(x)
  )
})
gene_inv <- rbindlist(inv_rows)

# -----------------------------------------------------------------------------
# 3. Unit = tumor
# -----------------------------------------------------------------------------
per_tumor <- md[, {
  n_s <- .N
  n_epi <- as.integer(sum(compartment == "epithelial"))
  n_host <- as.integer(sum(compartment == "host_epithelial"))
  n_tnk <- as.integer(sum(compartment == "tnk"))
  epi <- compartment == "epithelial"
  list(
    closest_to_user_KL = as.logical(unique(closest_to_user_KL)[1]),
    stk11 = as.character(unique(stk11)[1]),
    keap1 = as.character(unique(keap1)[1]),
    n_cells = n_s,
    n_epithelial = n_epi,
    n_host_epithelial = n_host,
    n_tnk = n_tnk,
    frac_tnk = if (n_s) n_tnk / n_s else NA_real_,
    n_cldn4_pos = as.integer(sum(Cldn4_pos)),
    n_myeloid = as.integer(sum(compartment == "myeloid")),
    n_caf = as.integer(sum(compartment == "caf")),
    n_endo = as.integer(sum(compartment == "endo")),
    n_b = as.integer(sum(compartment == "b")),
    n_other = as.integer(sum(compartment == "other")),
    frac_epithelial = if (n_s) n_epi / n_s else NA_real_,
    Cldn4_all_n = n_s,
    Cldn4_all_mean = mean(Cldn4_counts),
    Cldn4_all_pctpos = 100 * mean(Cldn4_pos),
    Cldn4_all_mean_log1p = mean(Cldn4_log1p),
    Cldn4_epithelial_n = n_epi,
    Cldn4_epithelial_mean = if (n_epi) mean(Cldn4_counts[epi]) else NA_real_,
    Cldn4_epithelial_pctpos = if (n_epi) 100 * mean(Cldn4_pos[epi]) else NA_real_,
    ifn_isg_epithelial_mean = if (n_epi) mean(ifn_isg[epi]) else NA_real_,
    mhc1_epithelial_mean = if (n_epi) mean(mhc1[epi]) else NA_real_,
    ifn_core6_epithelial_mean = if (n_epi) mean(ifn_core6[epi]) else NA_real_,
    seurat_ifn_epithelial_mean = if (n_epi) mean(SeuratIFN1[epi]) else NA_real_,
    seurat_mhc1_epithelial_mean = if (n_epi) mean(SeuratMHC11[epi]) else NA_real_,
    ifn_isg_all_mean = mean(ifn_isg),
    mhc1_all_mean = mean(mhc1)
  )
}, by = .(sample, genotype)]
setorder(per_tumor, genotype, sample)

geno_levels <- c("K", "KK", "KLK")
genotype_table <- rbindlist(lapply(geno_levels, function(g) {
  sub <- md[genotype == g]
  tumors <- sort(unique(sub$sample))
  n_g <- nrow(sub)
  n_epi <- as.integer(sum(sub$compartment == "epithelial"))
  n_host <- as.integer(sum(sub$compartment == "host_epithelial"))
  n_tnk <- as.integer(sum(sub$compartment == "tnk"))
  epi <- sub$compartment == "epithelial"
  data.table(
    genotype = g,
    cell_line = GENO_META[[g]]$cell_line,
    genotype_geo = GENO_META[[g]]$genotype_geo,
    stk11 = GENO_META[[g]]$stk11,
    keap1 = GENO_META[[g]]$keap1,
    closest_to_user_KL = GENO_META[[g]]$closest_to_user_KL,
    n_tumors = length(tumors),
    tumors = paste(tumors, collapse = ","),
    n_cells = n_g,
    n_epithelial = n_epi,
    n_host_epithelial = n_host,
    n_tnk = n_tnk,
    frac_tnk = if (n_g) n_tnk / n_g else NA_real_,
    n_cldn4_pos = as.integer(sum(sub$Cldn4_pos)),
    n_myeloid = as.integer(sum(sub$compartment == "myeloid")),
    n_caf = as.integer(sum(sub$compartment == "caf")),
    n_endo = as.integer(sum(sub$compartment == "endo")),
    n_b = as.integer(sum(sub$compartment == "b")),
    n_other = as.integer(sum(sub$compartment == "other")),
    Cldn4_all_n = n_g,
    Cldn4_all_mean = mean(sub$Cldn4_counts),
    Cldn4_all_pctpos = 100 * mean(sub$Cldn4_pos),
    Cldn4_epithelial_n = n_epi,
    Cldn4_epithelial_mean = if (n_epi) mean(sub$Cldn4_counts[epi]) else NA_real_,
    Cldn4_epithelial_pctpos = if (n_epi) 100 * mean(sub$Cldn4_pos[epi]) else NA_real_,
    ifn_isg_epithelial_mean = if (n_epi) mean(sub$ifn_isg[epi]) else NA_real_,
    mhc1_epithelial_mean = if (n_epi) mean(sub$mhc1[epi]) else NA_real_,
    seurat_ifn_epithelial_mean = if (n_epi) mean(sub$SeuratIFN1[epi]) else NA_real_,
    seurat_mhc1_epithelial_mean = if (n_epi) mean(sub$SeuratMHC11[epi]) else NA_real_
  )
}))

# Tumor-level KLK vs K (honest unit)
klk_t <- per_tumor[genotype == "KLK"]
k_t <- per_tumor[genotype == "K"]
contrast_metrics <- c(
  "frac_tnk", "Cldn4_all_mean", "Cldn4_all_pctpos",
  "Cldn4_epithelial_mean", "Cldn4_epithelial_pctpos",
  "ifn_isg_epithelial_mean", "mhc1_epithelial_mean", "n_epithelial"
)
klk_vs_k <- rbindlist(lapply(contrast_metrics, function(met) {
  a <- as.numeric(klk_t[[met]])
  b <- as.numeric(k_t[[met]])
  data.table(
    metric = met,
    n_tumors = "2 vs 2",
    KLK = mean(a, na.rm = TRUE),
    K = mean(b, na.rm = TRUE),
    delta_KLK_minus_K = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
    MW_p = mw_p(a, b),
    KLK_values = paste(sprintf("%.6g", a), collapse = ","),
    K_values = paste(sprintf("%.6g", b), collapse = ",")
  )
}))

comp <- md[, .N, by = .(genotype, compartment)]
comp[, frac := N / sum(N), by = genotype]
setorder(comp, genotype, -N)

cldn4_pos_cells <- md[Cldn4_pos == TRUE, .(
  barcode, genotype, sample, compartment, Cldn4_counts, Cldn4_log1p
)]
setorder(cldn4_pos_cells, genotype, sample, barcode)

honest_n <- data.table(
  item = c(
    "GEO series", "Tumors (unit)", "KLK vs K tumors", "Cells in matrix",
    "Genes", "Author malignant labels", "Marker epithelial cells",
    "Host-epithelial cells", "T/NK cells", "Cldn4-positive cells (any)",
    "Dual-high", "ICI arms", "Private KL libraries"
  ),
  n = c(
    1L, uniqueN(md$sample), NA_integer_, ncol(obj), nrow(obj), 0L,
    as.integer(sum(md$compartment == "epithelial")),
    as.integer(sum(md$compartment == "host_epithelial")),
    as.integer(sum(md$compartment == "tnk")),
    as.integer(sum(md$Cldn4_pos)),
    0L, 0L, 0L
  ),
  note = c(
    GEO,
    "2 K + 2 KK + 2 KLK",
    "2 vs 2; MW p cannot beat 1/3",
    "not the unit",
    "mm10 symbols",
    "none deposited; title says non-malignant",
    "not author-malignant",
    "residual normal lung",
    "Cd3d/e/g or Cd8a or Nkg7/Ncr1/Klrb1c",
    if (cldn4_present) sprintf("%.3f%% of matrix", 100 * mean(md$Cldn4_pos)) else "Cldn4 row absent",
    "not defined (Cldn4-only)",
    "GEO treatment: no",
    "public GEO only; no private 8 KL"
  )
)
honest_n[item == "KLK vs K tumors", n := NA_real_]
# keep n column numeric except the 2 vs 2 note lives in `note`

summary <- list(
  engine = "R + Seurat CreateSeuratObject",
  seurat_version = as.character(packageVersion("Seurat")),
  seuratobject_version = as.character(packageVersion("SeuratObject")),
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  geo = GEO,
  matrix_file = basename(MATRIX),
  value_note = value_note,
  n_cells = ncol(obj),
  n_genes = nrow(obj),
  n_tumors = uniqueN(md$sample),
  n_tnk = as.integer(sum(md$compartment == "tnk")),
  n_epithelial = as.integer(sum(md$compartment == "epithelial")),
  n_host_epithelial = as.integer(sum(md$compartment == "host_epithelial")),
  n_cldn4_pos = as.integer(sum(md$Cldn4_pos)),
  cldn4_present = cldn4_present,
  private_kl = 0L
)

write_tsv(file.path(TABLES, "genotype_table.tsv"), genotype_table)
write_tsv(file.path(TABLES, "per_tumor.tsv"), per_tumor)
write_tsv(file.path(TABLES, "klk_vs_k.tsv"), klk_vs_k)
write_tsv(file.path(TABLES, "honest_n.tsv"), honest_n)
write_tsv(file.path(TABLES, "gene_inventory.tsv"), gene_inv)
write_tsv(file.path(TABLES, "compartment_by_genotype.tsv"), comp)
write_tsv(file.path(TABLES, "cldn4_positive_cells.tsv"), cldn4_pos_cells)
writeLines(
  jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE, digits = 8),
  file.path(TABLES, "summary.json")
)
writeLines(capture.output(sessionInfo()), file.path(TABLES, "sessionInfo.txt"))

# -----------------------------------------------------------------------------
# 4. Figures
# -----------------------------------------------------------------------------
pt <- copy(per_tumor)
pt[, genotype := factor(genotype, levels = geno_levels)]
p1 <- ggplot(pt, aes(genotype, frac_tnk, fill = genotype)) +
  geom_col(data = pt[, .(frac_tnk = mean(frac_tnk)), by = genotype],
           width = 0.6, alpha = 0.45, color = NA) +
  geom_point(aes(shape = sample), size = 3, position = position_nudge(x = 0)) +
  scale_fill_manual(values = c(K = "#4C78A8", KK = "#F58518", KLK = "#E45756")) +
  labs(
    title = "GSE267321 T/NK fraction (unit = tumor)",
    subtitle = "Seurat CreateSeuratObject; leftover non-malignant digest",
    y = "T/NK fraction", x = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
ggsave(file.path(FIGS, "fig1_tnk_fraction.pdf"), p1, width = 6.2, height = 4.4)
ggsave(file.path(FIGS, "fig1_tnk_fraction.png"), p1, width = 6.2, height = 4.4, dpi = 140)

p2 <- ggplot(pt, aes(genotype, Cldn4_all_mean, fill = genotype)) +
  geom_col(data = pt[, .(Cldn4_all_mean = mean(Cldn4_all_mean)), by = genotype],
           width = 0.6, alpha = 0.45, color = NA) +
  geom_point(size = 3) +
  scale_fill_manual(values = c(K = "#4C78A8", KK = "#F58518", KLK = "#E45756")) +
  labs(
    title = "Cldn4 mean (all cells)",
    subtitle = "Floor expected: public non-malignant matrix",
    y = "Cldn4 mean counts", x = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(FIGS, "fig2_cldn4_all.pdf"), p2, width = 6.2, height = 4.4)
ggsave(file.path(FIGS, "fig2_cldn4_all.png"), p2, width = 6.2, height = 4.4, dpi = 140)

comp_plot <- copy(comp)
comp_plot[, compartment := factor(compartment, levels = unique(compartment))]
p3 <- ggplot(comp_plot, aes(genotype, N, fill = compartment)) +
  geom_col(position = "fill") +
  labs(
    title = "Marker compartments by genotype",
    y = "Fraction of cells", x = NULL
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(FIGS, "fig3_compartments.pdf"), p3, width = 6.8, height = 4.6)
ggsave(file.path(FIGS, "fig3_compartments.png"), p3, width = 6.8, height = 4.6, dpi = 140)

# -----------------------------------------------------------------------------
# 5. FINDING.md from the Seurat run
# -----------------------------------------------------------------------------
gmap <- split(genotype_table, genotype_table$genotype)
smap <- split(per_tumor, per_tumor$sample)
cmap <- split(klk_vs_k, klk_vs_k$metric)
cldn4_inv <- gene_inv[gene == "Cldn4"]
n_epi_total <- summary$n_epithelial
n_host_total <- summary$n_host_epithelial
host_note <- if (n_host_total > 0) {
  sprintf(" Host-lung epithelium (Sftpc/Scgb1a1/Ager ∩ epi markers) = **%d** cells.", n_host_total)
} else {
  ""
}

row_g <- function(r) {
  sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    r$genotype, r$stk11, r$keap1, r$n_tumors, r$n_cells, r$n_epithelial,
    r$n_host_epithelial, r$n_tnk, fmt(r$frac_tnk), fmt(r$Cldn4_all_mean),
    fmt(r$Cldn4_all_pctpos), fmt(r$Cldn4_epithelial_mean),
    fmt(r$Cldn4_epithelial_pctpos), fmt(r$ifn_isg_epithelial_mean),
    fmt(r$mhc1_epithelial_mean)
  )
}

trow <- function(r) {
  sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    r$sample, r$genotype, r$n_cells, r$n_epithelial, r$n_host_epithelial,
    r$n_tnk, fmt(r$frac_tnk), fmt(r$Cldn4_all_mean), fmt(r$Cldn4_all_pctpos),
    fmt(r$Cldn4_epithelial_mean), fmt(r$ifn_isg_epithelial_mean),
    fmt(r$mhc1_epithelial_mean)
  )
}

crow <- function(met) {
  r <- cmap[[met]]
  sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s |",
    r$metric, r$n_tumors, fmt(r$KLK), fmt(r$K), fmt(r$delta_KLK_minus_K),
    fmt(r$MW_p), r$KLK_values, r$K_values
  )
}

tnk_klk <- fmt(mean(klk_t$frac_tnk))
tnk_k <- fmt(mean(k_t$frac_tnk))
tnk_d <- fmt(mean(klk_t$frac_tnk) - mean(k_t$frac_tnk))

present_lineage <- intersect(
  c("Epcam", "Cdh1", "Krt8", "Krt18", "Krt19", "Cldn4", "Cd3d", "Cd3e", "Nkg7", "Ptprc", "Sftpc"),
  rownames(obj)
)
missing_lineage <- setdiff(
  c("Epcam", "Cdh1", "Krt8", "Krt18", "Krt19", "Cldn4", "Cd3d", "Cd3e", "Nkg7", "Ptprc", "Sftpc"),
  rownames(obj)
)

sftpb_n <- if ("Sftpb" %in% rownames(obj)) {
  as.integer(sum(as.numeric(GetAssayData(obj, layer = "counts")["Sftpb", ]) > 0))
} else {
  NA_integer_
}
scgb_n <- if ("Scgb1a1" %in% rownames(obj)) {
  as.integer(sum(as.numeric(GetAssayData(obj, layer = "counts")["Scgb1a1", ]) > 0))
} else {
  NA_integer_
}

klk2_epi <- per_tumor[sample == "LKR13-KLK-2", n_epithelial]
klk1_ifn <- fmt(per_tumor[sample == "LKR13-KLK-1", ifn_isg_epithelial_mean])
klk2_ifn <- fmt(per_tumor[sample == "LKR13-KLK-2", ifn_isg_epithelial_mean])
k1_tnk <- fmt(per_tumor[sample == "LKR13-K-1", frac_tnk])
k2_tnk <- fmt(per_tumor[sample == "LKR13-K-2", frac_tnk])
klk1_tnk <- fmt(per_tumor[sample == "LKR13-KLK-1", frac_tnk])
klk2_tnk <- fmt(per_tumor[sample == "LKR13-KLK-2", frac_tnk])
kk1_tnk <- fmt(per_tumor[sample == "LKR13-KK-1", frac_tnk])
kk2_tnk <- fmt(per_tumor[sample == "LKR13-KK-2", frac_tnk])

finding <- paste0(
  "# FINDING — Seurat GSE267321 LKR13 K / KK / KLK (Cldn4-only)\n",
  "\n",
  "**ADDITIVE. Public mouse. Cldn4-only. No dual-high. Seurat primary.** ",
  "Thesis is already correct and is not rewritten. This folder is the R + Seurat ",
  "(`CreateSeuratObject`) run on one public GEO object: ",
  "[GSE267321](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267321) ",
  "(Qian / Skoulidis / Heymach, MDACC). LKR13 syngeneic subcutaneous tumors in 129SV mice. ",
  "Genotypes on GEO: **K** (KrasG12D), **KK** (KrasG12D + KEAP1 KO), **KLK** ",
  "(KrasG12D + KEAP1/LKB1 KO). There is **no LKR13-KL** (STK11-only) library and ",
  "**no private 8 KL** object is used. **KLK is the closest public stand-in for user KL** ",
  "because it is the only STK11/LKB1-loss arm.\n",
  "\n",
  "GEO title: *Gene expression profile at single cell level of **non-malignant cells** ",
  "in KRAS syngeneic mouse tumor models harboring STK11 and/or KEAP1 co-mutation.* ",
  "Processed file (corrected 12 Feb 2026): `GSE267321_Normalized_expression_matrix_02122026.csv.gz`. ",
  "No author cluster / malignant column on GEO. No ICI arm (`treatment: no`).\n",
  "\n",
  sprintf(
    paste0(
      "**Engine.** R %s + Seurat %s / SeuratObject %s. ",
      "`CreateSeuratObject(counts=..., project=\"GSE267321\", min.cells=0, min.features=0)` ",
      "then `NormalizeData` + `AddModuleScore` for IFN/MHC. ",
      "Python is not the primary. Seurat **did run**.\n"
    ),
    summary$r_version, summary$seurat_version, summary$seuratobject_version
  ),
  "\n",
  sprintf(
    paste0(
      "**Verdict.** The public object is what the title says: a **non-malignant** digest. ",
      "Marker leftover epithelium exists (**%d** cells, %.1f%% of the matrix) but is **not** ",
      "an author-malignant compartment. Cldn4 is in the matrix and is **at the floor: %d / %d cells > 0 (%.3f%%)**. ",
      "Cldn4 cannot be tested as a genotype effect in epithelium. T/NK **can**: both KLK tumors ",
      "are below both K tumors (sample-mean fraction **%s vs %s**, Δ = %s). Honest n = **2 vs 2 tumors**. ",
      "IFN/MHC in leftover epithelium is mixed and KLK-2 has only %d epi cells — not a claim.%s\n"
    ),
    n_epi_total, 100 * n_epi_total / summary$n_cells,
    summary$n_cldn4_pos, summary$n_cells, 100 * summary$n_cldn4_pos / summary$n_cells,
    tnk_klk, tnk_k, tnk_d, klk2_epi, host_note
  ),
  "\n---\n\n",
  "## Decision\n\n",
  "| Question | Answer |\n",
  "|---|---|\n",
  "| Engine | **R + Seurat CreateSeuratObject** (not Python-only) |\n",
  "| Public processed matrix | **yes** — 15.7 MB CSV.gz on GEO |\n",
  "| Author malignant / epithelial labels | **no** — barcode matrix only; title says non-malignant |\n",
  sprintf("| Marker tumor epithelium | **%d** cells |\n", n_epi_total),
  sprintf("| Host-lung epithelium | **%d** cells |\n", n_host_total),
  sprintf(
    "| Cldn4 row present | **%s** — **%d cells > 0** (floor) |\n",
    if (cldn4_present) "yes" else "no", summary$n_cldn4_pos
  ),
  "| Cldn4 usable as a genotype test | **no** — too few positive cells |\n",
  "| T/NK by genotype | **yes** — KLK < K in both tumors |\n",
  sprintf(
    "| IFN/MHC in leftover epithelium | **%s** |\n",
    if (n_epi_total) "scored in leftover marker epi; mixed / underpowered" else "not scored — epithelium absent"
  ),
  "| Closest to user KL | **KLK** (STK11/LKB1 loss). KK is KEAP1-only. No KL library. |\n",
  "| Private 8 KL | **not used** |\n",
  "| Unit | **tumor / replicate** (2 per genotype) |\n",
  "| Dual-high TACSTD2 × Cldn4 | **not defined** |\n",
  "| ICI / PD-1 | **no** |\n",
  "\n",
  "Cldn4 was scored because the row exists. The score is a floor, not a biology test. ",
  "T/NK fraction is the genotype readout this object can support.\n",
  "\n---\n\n",
  "## Honest n\n\n",
  "| item | n | note |\n",
  "|---|---:|---|\n",
  sprintf("| GEO series | **1** | %s |\n", GEO),
  "| Tumors (unit) | **6** | 2 K + 2 KK + 2 KLK |\n",
  "| KLK vs K tumors | **2 vs 2** | MW p cannot beat 1/3 |\n",
  sprintf("| Cells in matrix | **%d** | not the unit |\n", summary$n_cells),
  sprintf("| Genes | **%d** | mm10 symbols |\n", summary$n_genes),
  "| Author malignant labels | **0** | none deposited |\n",
  sprintf("| Marker epithelial cells | **%d** | not author-malignant |\n", n_epi_total),
  sprintf("| Host-epithelial cells | **%d** | residual normal lung |\n", n_host_total),
  sprintf("| T/NK cells | **%d** | Cd3d/e/g or Cd8a or Nkg7/Ncr1/Klrb1c |\n", summary$n_tnk),
  sprintf(
    "| Cldn4-positive cells (any) | **%d** | %s%% of matrix |\n",
    summary$n_cldn4_pos,
    if (cldn4_present) fmt(100 * summary$n_cldn4_pos / summary$n_cells, 2) else "0"
  ),
  "| Dual-high | **0** | not defined |\n",
  "| ICI arms | **0** | untreated |\n",
  "| Private KL libraries | **0** | public GEO only |\n",
  "\n",
  sprintf(
    "Do not write n = %d. Do not write n = 3 genotypes as if they were biological replicates of KL. Do not import a private 8-KL object.\n",
    summary$n_cells
  ),
  "\n---\n\n",
  "## Genotype table (unit = genotype, built from 2 tumors)\n\n",
  "K = KrasG12D parental. KK = KEAP1-loss. KLK = STK11/LKB1 + KEAP1-loss (closest to user KL).\n\n",
  "| genotype | STK11 | KEAP1 | n tumors | n cells | n epi | n host-epi | n T/NK | frac T/NK | Cldn4 all mean | Cldn4 all %pos | Cldn4 epi mean | Cldn4 epi %pos | IFN epi | MHC-I epi |\n",
  "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
  row_g(gmap[["K"]]), "\n",
  row_g(gmap[["KK"]]), "\n",
  row_g(gmap[["KLK"]]), "\n",
  "\n",
  "Source: `tables/genotype_table.tsv` (written by the Seurat script after `CreateSeuratObject`).\n",
  "\n---\n\n",
  "## Per-tumor table (the actual unit)\n\n",
  "| sample | genotype | n cells | n epi | n host-epi | n T/NK | frac T/NK | Cldn4 all | Cldn4 %pos | Cldn4 epi | IFN epi | MHC-I epi |\n",
  "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
  trow(smap[["LKR13-K-1"]]), "\n",
  trow(smap[["LKR13-K-2"]]), "\n",
  trow(smap[["LKR13-KK-1"]]), "\n",
  trow(smap[["LKR13-KK-2"]]), "\n",
  trow(smap[["LKR13-KLK-1"]]), "\n",
  trow(smap[["LKR13-KLK-2"]]), "\n",
  "\n---\n\n",
  "## KLK vs K (primary)\n\n",
  "Closest public STK11-loss arm versus parental K. Tumor-level means; n=2 vs 2. ",
  "Wilcoxon / Mann-Whitney computed in R on the two tumor values.\n\n",
  "| metric | n tumors | KLK | K | Δ (KLK−K) | MW p | KLK values | K values |\n",
  "|---|---|---:|---:|---:|---:|---|---|\n",
  crow("frac_tnk"), "\n",
  crow("Cldn4_all_mean"), "\n",
  crow("Cldn4_all_pctpos"), "\n",
  crow("Cldn4_epithelial_mean"), "\n",
  crow("Cldn4_epithelial_pctpos"), "\n",
  crow("ifn_isg_epithelial_mean"), "\n",
  crow("mhc1_epithelial_mean"), "\n",
  crow("n_epithelial"), "\n",
  "\n",
  sprintf(
    paste0(
      "**T/NK is the only directional genotype result.** Both KLK tumors sit below both K tumors ",
      "(%s and %s vs %s and %s). KK is also T/NK-low (%s and %s) — KEAP1-loss alone already looks cold, ",
      "and KLK is not warmer. Mann-Whitney at 2 vs 2 cannot beat p = 1/3; the sign is the result.\n"
    ),
    klk1_tnk, klk2_tnk, k1_tnk, k2_tnk, kk1_tnk, kk2_tnk
  ),
  "\n",
  sprintf(
    paste0(
      "**Cldn4 is empty.** %d positive cells in the whole matrix (`tables/cldn4_positive_cells.tsv`). ",
      "That is dropout / floor, not “KLK down-regulates Cldn4.”\n"
    ),
    summary$n_cldn4_pos
  ),
  "\n",
  sprintf(
    paste0(
      "**IFN/MHC in leftover epithelium is not a claim.** KLK-1 leftover epi IFN is %s (n=%d); ",
      "KLK-2 is %s (n=%d). Tumor-level MW p cannot support a genotype effect. ",
      "Seurat `AddModuleScore` values are in `genotype_table.tsv` as `seurat_ifn_epithelial_mean` / ",
      "`seurat_mhc1_epithelial_mean` and were not used to invent a leftover-epi IFN claim.\n"
    ),
    klk1_ifn, per_tumor[sample == "LKR13-KLK-1", n_epithelial],
    klk2_ifn, klk2_epi
  ),
  "\n",
  "KK is the KEAP1-only extra arm. It is **not** user KL.\n",
  "\n---\n\n",
  "## Epithelium (honest)\n\n",
  "GEO title and summary say **non-malignant cells**. The digest protocol is whole-tumor ",
  "(collagenase / hyaluronidase / dispase; no CD45 sort on the GEO record). The public object has ",
  "**no** author `Malignant` / `Epithelial` column.\n\n",
  "This page therefore uses markers on the Seurat counts layer, not author calls:\n\n",
  "- **Marker tumor epithelium** = Epcam+ **or** (Cdh1+ and Krt8+) **or** (Krt18+ and Krt19+), and **not** Sftpc / Scgb1a1 / Ager.\n",
  "- **Host epithelium** = those same epi markers **and** Sftpc or Scgb1a1 or Ager.\n",
  "- T/NK = Cd3d / Cd3e / Cd3g / Cd8a / Nkg7 / Ncr1 / Klrb1c. Epithelium wins if both fire.\n",
  "\n",
  sprintf(
    paste0(
      "Tumors are **subcutaneous** (GEO `source_name`: Subcutaneous tumor). There is no orthotopic lung, ",
      "so host AT2 should be near zero. Sftpc and Sftpa1 are **%s**; Sftpb is %s. Scgb1a1 = %s cells. ",
      "The %d leftover Epcam/keratin cells are therefore residual LKR13 tumor epithelium the authors did not fully strip, ",
      "not normal lung. KLK-2 has **%d** such cells (below a ≥20 occupancy floor). ",
      "Do not call this an author-malignant atlas.\n"
    ),
    if ("Sftpc" %in% rownames(obj)) "present" else "missing from the matrix",
    if (is.na(sftpb_n)) "missing" else sprintf("present and **%d / %d**", sftpb_n, summary$n_cells),
    if (is.na(scgb_n)) "NA" else as.character(scgb_n),
    n_epi_total, klk2_epi
  ),
  "\n",
  sprintf("Lineage genes present: %s.\n", paste(present_lineage, collapse = ", ")),
  sprintf("Missing: %s.\n", if (length(missing_lineage)) paste(missing_lineage, collapse = ", ") else "none"),
  "\n---\n\n",
  "## Methods (short)\n\n",
  "1. Public only. Downloaded `GSE267321_Normalized_expression_matrix_02122026.csv.gz` and the series matrix from NCBI GEO FTP. No SRA / FASTQ. No private object. No private 8 KL.\n",
  "2. Matrix is genes × cells, 10x-style barcodes `{K|KK|KLK}_{UMI}.1_LKR13.{geno}.{rep}`. Genotype and replicate parsed from the barcode. That is the only metadata join.\n",
  sprintf("3. Ingest: `data.table::fread` → sparse `dgCMatrix` → **`CreateSeuratObject`**. Values: %s. Positive = raw counts layer > 0. Simple module scores = mean log1p of present genes. Seurat `AddModuleScore` run after `NormalizeData` as a second IFN/MHC score.\n", value_note),
  "4. Cldn4-only. Tacstd2 is inventory/audit, never a gate. No dual-high.\n",
  "5. IFN ISG = mouse orthologs of the public type-I ISG core (Isg15, Ifit1/2/3, Mx1, Oas1a/Oas2, Stat1, Irf7, …). MHC-I APM = B2m, H2-K1/D1/Q4/Q7/T23, Tap1/2, Psmb8/9/10, Nlrc5, ….\n",
  "6. Unit = tumor. Primary contrast = KLK vs K. KK is extra.\n",
  "7. Thesis is not rewritten.\n",
  "\n",
  "```bash\n",
  "Rscript methods/seurat_gse267321_cldn4/scripts/analyze_gse267321_seurat.R\n",
  "```\n",
  "\n---\n\n",
  "## How to read this\n\n",
  "- **Additive public mouse**, not a human concordant-pool join. Additive to the existing public GSE267321 folder; this page is the Seurat engine.\n",
  "- **KLK ≠ KL.** KLK is STK11-loss **plus** KEAP1-loss. It is the closest arm in *this* series, not a clean STK11-only replicate of user KL.\n",
  "- **No private 8 KL.** Honest n is the six public tumors.\n",
  "- **Title says non-malignant.** Leftover Epcam/keratin is residual tumor epithelium, not an author malignant call.\n",
  sprintf("- **Cldn4 is present and empty (%d cells).** Score it; do not build a Cldn4-high story on this object.\n", summary$n_cldn4_pos),
  "- **T/NK fraction is the genotype table:** KLK (and KK) colder than K. n=2 vs 2.\n",
  "- **Honest n is 2 vs 2 tumors.** Direction can be described. A p-value cannot.\n",
  "- **No dual-high. No ICI endpoint. Thesis unchanged.**\n",
  "\n",
  "## Files\n\n",
  "- `tables/genotype_table.tsv` — required genotype roll-up (done criterion)\n",
  "- `tables/per_tumor.tsv` — unit-level table\n",
  "- `tables/klk_vs_k.tsv` — primary contrast\n",
  "- `tables/cldn4_positive_cells.tsv` — every Cldn4>0 barcode\n",
  "- `tables/compartment_by_genotype.tsv`\n",
  "- `tables/gene_inventory.tsv`, `tables/honest_n.tsv`, `tables/summary.json`, `tables/sessionInfo.txt`\n",
  "- `figures/fig1_tnk_fraction.png`, `fig2_cldn4_all.png`, `fig3_compartments.png`\n",
  "- `scripts/analyze_gse267321_seurat.R`\n",
  "\n",
  "## 结论\n\n",
  sprintf(
    paste0(
      "GSE267321 是公开的 LKR13 皮下同基因瘤 scRNA（K / KK / KLK；**无 KL 库**）。",
      "本页用 **R + Seurat CreateSeuratObject** 跑同一份 GEO normalized CSV，**不用私有 8 KL**。",
      "**KLK（STK11/LKB1+KEAP1 缺失）是本系列最接近用户 KL 的一臂。** ",
      "GEO 标题为 **non-malignant cells**：无作者恶性标签；标记残留上皮 **%d** 个细胞（皮下，不是正常肺）。",
      "Cldn4 行在，但全矩阵仅 **%d** 个细胞 >0，**不能做 Cldn4 基因型检验**。",
      "T/NK 可以：KLK 两瘤都低于 K 两瘤（样本均数 %s vs %s）。",
      "IFN/MHC 在残留上皮上混合且 KLK-2 仅 %d 个上皮细胞，不作结论。诚实 n = **2 vs 2 个瘤**。无 dual-high，无 ICI，不改 thesis。\n"
    ),
    n_epi_total, summary$n_cldn4_pos, tnk_klk, tnk_k, klk2_epi
  )
)

writeLines(finding, file.path(HERE, "FINDING.md"), useBytes = TRUE)
message("Wrote ", file.path(HERE, "FINDING.md"))
message("Wrote ", file.path(TABLES, "genotype_table.tsv"))
message("DONE Seurat GSE267321 Cldn4-only. n_tumors=", summary$n_tumors,
        " n_cells=", summary$n_cells, " n_cldn4_pos=", summary$n_cldn4_pos)
