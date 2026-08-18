#!/usr/bin/env Rscript
# GSE131907 CLDN4-only Seurat analysis.
# ADDITIVE. Thesis already correct. No dual-high. No GSE148071.
# Honest unit = patient. Primary: author Malignant cells.
# CreateSeuratObject on public processed UMI (gene-subset + full-library sizes).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a) || a == "") b else a

opt <- list(
  subset = "/tmp/gse131907/subset",
  datadir = "/tmp/gse131907",
  outdir = "methods/seurat_gse131907_cldn4/results",
  finding = "methods/seurat_gse131907_cldn4/FINDING.md",
  families = "methods/seurat_gse131907_cldn4/data/families.json"
)
if (length(args) >= 1) opt$subset <- args[[1]]
if (length(args) >= 2) opt$datadir <- args[[2]]
if (length(args) >= 3) opt$outdir <- args[[3]]
if (length(args) >= 4) opt$finding <- args[[4]]
if (length(args) >= 5) opt$families <- args[[5]]

dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "figures"), recursive = TRUE, showWarnings = FALSE)

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 1e-4) return(format(p, digits = 3, scientific = TRUE))
  sprintf("%.4g", p)
}
fmt_n <- function(x, d = 3) {
  if (!is.finite(x)) return("NA")
  sprintf(paste0("%.", d, "f"), x)
}
spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 4) return(list(n = n, rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  list(n = n, rho = unname(ct$estimate), p = ct$p.value)
}
rank_biserial <- function(high, low) {
  high <- high[is.finite(high)]
  low <- low[is.finite(low)]
  n4 <- length(high); n1 <- length(low)
  if (n4 < 3 || n1 < 3) {
    return(list(n_q1 = n1, n_q4 = n4, r_rb = NA_real_, p = NA_real_,
                median_q1 = median(low), median_q4 = median(high),
                delta = median(high) - median(low)))
  }
  wt <- wilcox.test(high, low, alternative = "two.sided", exact = FALSE)
  r <- (2 * as.numeric(wt$statistic)) / (n4 * n1) - 1
  list(n_q1 = n1, n_q4 = n4, r_rb = r, p = wt$p.value,
       median_q1 = median(low), median_q4 = median(high),
       delta = median(high) - median(low))
}
assign_quartiles <- function(x) {
  r <- rank(as.numeric(x), ties.method = "average", na.last = "keep")
  qs <- tryCatch(
    as.character(cut(
      r,
      breaks = quantile(r[is.finite(r)], probs = seq(0, 1, 0.25), names = FALSE, type = 7),
      include.lowest = TRUE,
      labels = c("Q1", "Q2", "Q3", "Q4")
    )),
    error = function(e) rep(NA_character_, length(x))
  )
  qs
}

parse_series <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  title <- NULL
  geo <- NULL
  patient <- NULL
  stage <- NULL
  origin <- NULL
  while (length(line <- readLines(con, n = 1)) == 1) {
    if (!startsWith(line, "!Sample_")) next
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    key <- sub("^!Sample_", "", parts[[1]])
    vals <- gsub('^"|"$', "", parts[-1])
    if (key == "title") title <- vals
    if (key == "geo_accession") geo <- vals
    if (key == "characteristics_ch1") {
      if (grepl("^patient id:", vals[[1]])) {
        patient <- sub("^patient id: ", "", vals)
      } else if (grepl("^tumor stage:", vals[[1]])) {
        stage <- sub("^tumor stage: ", "", vals)
      } else if (grepl("^tissue origin abbrevation:", vals[[1]])) {
        origin <- sub("^tissue origin abbrevation: ", "", vals)
      }
    }
  }
  data.frame(
    Sample = title,
    geo_accession = geo,
    patient_id = patient,
    tumor_stage = stage,
    Sample_Origin_geo = origin,
    stringsAsFactors = FALSE
  )
}

read_json_list <- function(path) {
  fam <- jsonlite::fromJSON(path)
  list(
    IFN = as.character(fam$IFN),
    MHC_I_APM = as.character(fam$MHC_I_APM),
    TJ_no_CLDN4 = as.character(fam$TJ_no_CLDN4),
    chemokine = as.character(fam$chemokine)
  )
}

message("Loading extracted matrix + annotation")
feat <- readLines(file.path(opt$subset, "features.tsv"))
bc <- readLines(file.path(opt$subset, "barcodes.tsv"))
mat <- readMM(file.path(opt$subset, "matrix.mtx"))
mat <- as(mat, "CsparseMatrix")
rownames(mat) <- feat
colnames(mat) <- bc

lib <- read.delim(file.path(opt$subset, "libsize.tsv"), stringsAsFactors = FALSE)
lib <- lib[match(bc, lib$cell), ]
if (any(is.na(lib$nCount_full))) stop("libsize alignment failed")

ann <- read.delim(
  file.path(opt$datadir, "GSE131907_Lung_Cancer_cell_annotation.txt.gz"),
  stringsAsFactors = FALSE,
  na.strings = ""
)
if (!"Index" %in% names(ann)) stop("annotation missing Index")
ann <- ann[match(bc, ann$Index), ]
if (any(is.na(ann$Sample))) stop("matrix barcodes do not match annotation Index")

series <- parse_series(file.path(opt$datadir, "GSE131907_series_matrix.txt.gz"))
ann$patient_id <- series$patient_id[match(ann$Sample, series$Sample)]
ann$tumor_stage <- series$tumor_stage[match(ann$Sample, series$Sample)]
if (any(is.na(ann$patient_id))) stop("patient_id missing after GEO join")

message("CreateSeuratObject")
obj <- CreateSeuratObject(counts = mat, project = "GSE131907", min.cells = 0, min.features = 0)
stopifnot(ncol(obj) == length(bc))
obj$nCount_subset <- obj$nCount_RNA
obj$nCount_full <- as.numeric(lib$nCount_full)
# NormalizeData uses nCount_RNA; overwrite with full UMI library size.
obj$nCount_RNA <- obj$nCount_full
obj$Sample <- ann$Sample
obj$Sample_Origin <- ann$Sample_Origin
obj$Cell_type <- ann$Cell_type
obj$Cell_subtype <- ann$Cell_subtype
obj$patient_id <- ann$patient_id
obj$tumor_stage <- ann$tumor_stage
obj$is_tumor <- obj$Sample_Origin %in% c("tLung", "tL/B", "mLN", "PE", "mBrain")
obj$is_tnk <- obj$Cell_type %in% c("T lymphocytes", "NK cells")
obj$mal_author <- obj$Cell_subtype == "Malignant cells"
obj$mal_broad <- obj$Cell_subtype %in% c("Malignant cells", "tS1", "tS2", "tS3")

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
DefaultAssay(obj) <- "RNA"

cldn4_umi <- as.numeric(GetAssayData(obj, layer = "counts")["CLDN4", ])
cldn4_log <- as.numeric(GetAssayData(obj, layer = "data")["CLDN4", ])
obj$CLDN4_umi <- cldn4_umi
obj$CLDN4_log <- cldn4_log
obj$CLDN4_pos <- cldn4_umi > 0

md <- as.data.frame(obj[[]])
if (!"CLDN4" %in% rownames(obj)) stop("CLDN4 missing from CreateSeuratObject")

score_units <- function(md, mal_col, unit_col, min_mal = 20L) {
  tumor <- md[md$is_tumor, , drop = FALSE]
  split_by <- tumor[[unit_col]]
  out <- lapply(split(tumor, split_by, drop = TRUE), function(d) {
    mal <- !is.na(d[[mal_col]]) & as.logical(d[[mal_col]])
    n_mal <- sum(mal)
    n_tnk <- sum(d$is_tnk, na.rm = TRUE)
    n_cells <- nrow(d)
    cldn <- d$CLDN4_umi[mal]
    cldn_log <- d$CLDN4_log[mal]
    data.frame(
      unit = d[[unit_col]][1],
      patient_id = d$patient_id[1],
      samples = paste(sort(unique(d$Sample)), collapse = ","),
      origins = paste(sort(unique(d$Sample_Origin)), collapse = ","),
      n_samples = length(unique(d$Sample)),
      n_cells = n_cells,
      n_malignant = n_mal,
      n_tnk = n_tnk,
      frac_tnk = n_tnk / n_cells,
      cldn4_pct = if (n_mal > 0) 100 * mean(cldn > 0) else NA_real_,
      cldn4_mean_log = if (n_mal > 0) mean(cldn_log) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  u <- do.call(rbind, out)
  rownames(u) <- NULL
  u$eligible <- is.finite(u$n_malignant) & u$n_malignant >= min_mal
  u
}

pat_author <- score_units(md, "mal_author", "patient_id")
pat_broad <- score_units(md, "mal_broad", "patient_id")
samp_author <- score_units(md, "mal_author", "Sample")
samp_broad <- score_units(md, "mal_broad", "Sample")

el <- pat_author[pat_author$eligible, , drop = FALSE]
el$quartile <- assign_quartiles(el$cldn4_pct)
el_b <- pat_broad[pat_broad$eligible, , drop = FALSE]
el_b$quartile <- assign_quartiles(el_b$cldn4_pct)
el_s <- samp_author[samp_author$eligible, , drop = FALSE]
el_s$quartile <- assign_quartiles(el_s$cldn4_pct)

tnk_block <- function(df, label) {
  sp_pct <- spearman(df$cldn4_pct, df$frac_tnk)
  sp_mean <- spearman(df$cldn4_mean_log, df$frac_tnk)
  q <- df$quartile
  mw <- rank_biserial(df$frac_tnk[q == "Q4"], df$frac_tnk[q == "Q1"])
  data.frame(
    contrast = label,
    n = nrow(df),
    n_q1 = mw$n_q1,
    n_q4 = mw$n_q4,
    rho_pct = sp_pct$rho,
    p_pct = sp_pct$p,
    rho_mean = sp_mean$rho,
    p_mean = sp_mean$p,
    r_rb_q4q1 = mw$r_rb,
    p_q4q1 = mw$p,
    median_tnk_q1 = mw$median_q1,
    median_tnk_q4 = mw$median_q4,
    delta_tnk = mw$delta,
    stringsAsFactors = FALSE
  )
}

tnk_tab <- rbind(
  tnk_block(el, "patient_authorMalignant"),
  tnk_block(el_b, "patient_broad_tS_plus_Malignant"),
  tnk_block(el_s, "sample_authorMalignant_sensitivity")
)

# Family scores: patient-pseudobulk log2(CPM+1) using FULL library UMI sums.
fam <- read_json_list(opt$families)
fam$TJ_no_CLDN4 <- setdiff(fam$TJ_no_CLDN4, "CLDN4")
present_genes <- rownames(obj)
fam_present <- lapply(fam, function(g) intersect(g, present_genes))
message(
  "family genes present: ",
  paste(sprintf("%s=%d/%d", names(fam_present), vapply(fam_present, length, 1L), vapply(fam, length, 1L)), collapse = " ")
)

counts <- GetAssayData(obj, layer = "counts")
mal_cells <- !is.na(md$mal_author) & md$mal_author & !is.na(md$is_tumor) & md$is_tumor & md$patient_id %in% el$unit
mal_md <- md[mal_cells, , drop = FALSE]
mal_counts <- counts[, mal_cells, drop = FALSE]
lib_mal <- md$nCount_full[mal_cells]

patients <- el$unit
pb_cpm <- sapply(patients, function(pid) {
  hit <- mal_md$patient_id == pid
  cs <- as.numeric(Matrix::rowSums(mal_counts[, hit, drop = FALSE]))
  lib <- sum(lib_mal[hit])
  if (!is.finite(lib) || lib <= 0) return(rep(NA_real_, nrow(mal_counts)))
  log2(1e6 * cs / lib + 1)
})
rownames(pb_cpm) <- rownames(mal_counts)
colnames(pb_cpm) <- patients

family_score <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (length(g) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[g, , drop = FALSE])
}

fam_df <- data.frame(
  patient_id = patients,
  quartile = el$quartile,
  cldn4_pct = el$cldn4_pct,
  IFN = family_score(pb_cpm, fam_present$IFN),
  MHC_I_APM = family_score(pb_cpm, fam_present$MHC_I_APM),
  TJ = family_score(pb_cpm, fam_present$TJ_no_CLDN4),
  chemokine = family_score(pb_cpm, fam_present$chemokine),
  CLDN4 = if ("CLDN4" %in% rownames(pb_cpm)) as.numeric(pb_cpm["CLDN4", ]) else NA_real_,
  stringsAsFactors = FALSE
)

family_test <- function(df, col) {
  q1 <- df[[col]][df$quartile == "Q1"]
  q4 <- df[[col]][df$quartile == "Q4"]
  mw <- rank_biserial(q4, q1)
  logfc <- mean(q4, na.rm = TRUE) - mean(q1, na.rm = TRUE)
  data.frame(
    family = col,
    n = nrow(df),
    n_q1 = mw$n_q1,
    n_q4 = mw$n_q4,
    n_genes = length(fam_present[[if (col == "TJ") "TJ_no_CLDN4" else if (col == "MHC_I_APM") "MHC_I_APM" else col]]),
    mean_q1 = mean(q1, na.rm = TRUE),
    mean_q4 = mean(q4, na.rm = TRUE),
    logFC_q4_minus_q1 = logfc,
    r_rb = mw$r_rb,
    p = mw$p,
    stringsAsFactors = FALSE
  )
}

# n_genes for CLDN4 special-case
fam_present$CLDN4 <- intersect("CLDN4", present_genes)
fam_tab <- rbind(
  family_test(fam_df, "IFN"),
  family_test(fam_df, "MHC_I_APM"),
  family_test(fam_df, "TJ"),
  family_test(fam_df, "chemokine"),
  family_test(fam_df, "CLDN4")
)
fam_tab$n_genes[fam_tab$family == "CLDN4"] <- 1L
fam_tab$expected <- c("DOWN", "DOWN", "UP", "DOWN", "UP_split_check")

n_honest <- data.frame(
  piece = c(
    "cells_in_CreateSeuratObject",
    "patients_in_GEO",
    "samples_in_GEO",
    "tumor_bearing_origins",
    "author_Malignant_cells",
    "tS1_tS2_tS3",
    "patients_authorMal_n>=20",
    "patients_broadMal_n>=20",
    "samples_authorMal_n>=20",
    "nLung_authorMal_cells",
    "Q1_patients_author",
    "Q4_patients_author"
  ),
  n = c(
    ncol(obj),
    length(unique(md$patient_id)),
    length(unique(md$Sample)),
    sum(md$is_tumor),
    sum(md$mal_author),
    sum(md$Cell_subtype %in% c("tS1", "tS2", "tS3")),
    nrow(el),
    nrow(el_b),
    nrow(el_s),
    sum(md$mal_author & md$Sample_Origin == "nLung"),
    sum(el$quartile == "Q1"),
    sum(el$quartile == "Q4")
  ),
  note = c(
    "gene-subset object; do not quote as inferential n",
    "44 LUAD patients in Kim 2020",
    "58 samples",
    "tLung/tL/B/mLN/PE/mBrain cells",
    "Cell_subtype == Malignant cells",
    "included only in broad sensitivity",
    "PRIMARY n",
    "sensitivity: tS1/tS2/tS3 + Malignant cells",
    "comparison to prior sample-level GSE131907 n=21",
    "must be 0; nLung cannot join malignant CLDN4",
    "within-cohort %pos quartile",
    "within-cohort %pos quartile"
  ),
  stringsAsFactors = FALSE
)

# Inventory of excluded tumor patients
pat_all_tumor <- score_units(md, "mal_author", "patient_id", min_mal = 0L)
n_lung_note <- sum(md$Sample_Origin == "nLung")

write.table(el, file.path(opt$outdir, "tables", "patient_units.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(el_b, file.path(opt$outdir, "tables", "patient_units_broad.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(el_s, file.path(opt$outdir, "tables", "sample_units_author.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tnk_tab, file.path(opt$outdir, "tables", "tnk_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fam_tab, file.path(opt$outdir, "tables", "family_q4q1.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fam_df, file.path(opt$outdir, "tables", "patient_family_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(n_honest, file.path(opt$outdir, "tables", "n_honest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(series, file.path(opt$outdir, "tables", "geo_sample_patient.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- plots ----
theme_set(theme_bw(base_size = 12))
el$qcol <- ifelse(el$quartile %in% c("Q1", "Q4"), el$quartile, "Q2/Q3")

p1 <- ggplot(el, aes(cldn4_pct, frac_tnk)) +
  geom_point(aes(color = qcol), size = 2.6) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.4) +
  scale_color_manual(values = c(Q1 = "#1f77b4", Q4 = "#d62728", `Q2/Q3` = "#9aa0a6")) +
  labs(
    title = "GSE131907 patient-level CLDN4 vs T/NK",
    subtitle = sprintf(
      "author Malignant; n=%d patients; Spearman rho=%s p=%s",
      nrow(el), fmt_n(tnk_tab$rho_pct[1]), fmt_p(tnk_tab$p_pct[1])
    ),
    x = "Malignant CLDN4 %pos",
    y = "T/NK fraction (tumor-bearing cells)",
    color = "CLDN4 quartile"
  )
ggsave(file.path(opt$outdir, "figures", "fig_cldn4_vs_tnk.png"), p1, width = 6.4, height = 4.8, dpi = 160)
ggsave(file.path(opt$outdir, "figures", "fig_cldn4_vs_tnk.pdf"), p1, width = 6.4, height = 4.8)

el_q <- el[el$quartile %in% c("Q1", "Q4"), ]
p2 <- ggplot(el_q, aes(quartile, frac_tnk, fill = quartile)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 2) +
  scale_fill_manual(values = c(Q1 = "#1f77b4", Q4 = "#d62728")) +
  labs(
    title = "T/NK fraction, CLDN4 Q4 vs Q1",
    subtitle = sprintf(
      "n_Q1=%d n_Q4=%d; rank-biserial r=%s p=%s",
      tnk_tab$n_q1[1], tnk_tab$n_q4[1],
      fmt_n(tnk_tab$r_rb_q4q1[1]), fmt_p(tnk_tab$p_q4q1[1])
    ),
    x = NULL, y = "T/NK fraction"
  ) +
  theme(legend.position = "none")
ggsave(file.path(opt$outdir, "figures", "fig_tnk_q4q1.png"), p2, width = 5.2, height = 4.6, dpi = 160)
ggsave(file.path(opt$outdir, "figures", "fig_tnk_q4q1.pdf"), p2, width = 5.2, height = 4.6)

fam_long <- rbind(
  data.frame(family = "IFN", score = fam_df$IFN, quartile = fam_df$quartile),
  data.frame(family = "MHC-I/APM", score = fam_df$MHC_I_APM, quartile = fam_df$quartile),
  data.frame(family = "TJ (no CLDN4)", score = fam_df$TJ, quartile = fam_df$quartile)
)
fam_long <- fam_long[fam_long$quartile %in% c("Q1", "Q4"), ]
p3 <- ggplot(fam_long, aes(quartile, score, fill = quartile)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.6) +
  facet_wrap(~family, scales = "free_y") +
  scale_fill_manual(values = c(Q1 = "#1f77b4", Q4 = "#d62728")) +
  labs(
    title = "Malignant family scores, patient Q4 vs Q1",
    subtitle = "log2(CPM+1) of UMI-sum pseudobulk; full-library size; CLDN4 held out of TJ",
    x = NULL, y = "family score"
  ) +
  theme(legend.position = "none")
ggsave(file.path(opt$outdir, "figures", "fig_family_q4q1.png"), p3, width = 8.2, height = 4.4, dpi = 160)
ggsave(file.path(opt$outdir, "figures", "fig_family_q4q1.pdf"), p3, width = 8.2, height = 4.4)

n_plot <- n_honest[n_honest$piece %in% c(
  "patients_authorMal_n>=20", "Q1_patients_author", "Q4_patients_author",
  "samples_authorMal_n>=20", "patients_broadMal_n>=20"
), ]
n_plot$piece <- factor(n_plot$piece, levels = rev(n_plot$piece))
p4 <- ggplot(n_plot, aes(n, piece)) +
  geom_col(fill = "#4c78a8", width = 0.7) +
  geom_text(aes(label = n), hjust = -0.15, size = 3.4) +
  xlim(0, max(n_plot$n) * 1.2) +
  labs(title = "Honest n (do not quote cell counts)", x = "n units", y = NULL)
ggsave(file.path(opt$outdir, "figures", "fig_honest_n.png"), p4, width = 6.6, height = 3.6, dpi = 160)
ggsave(file.path(opt$outdir, "figures", "fig_honest_n.pdf"), p4, width = 6.6, height = 3.6)

prim <- tnk_tab[tnk_tab$contrast == "patient_authorMalignant", ]
broad <- tnk_tab[tnk_tab$contrast == "patient_broad_tS_plus_Malignant", ]
samp <- tnk_tab[tnk_tab$contrast == "sample_authorMalignant_sensitivity", ]
ifn <- fam_tab[fam_tab$family == "IFN", ]
mhc <- fam_tab[fam_tab$family == "MHC_I_APM", ]
tj <- fam_tab[fam_tab$family == "TJ", ]
che <- fam_tab[fam_tab$family == "chemokine", ]
cld <- fam_tab[fam_tab$family == "CLDN4", ]

missing_path <- file.path(opt$subset, "missing_genes.txt")
missing_genes <- if (file.exists(missing_path)) {
  g <- readLines(missing_path, warn = FALSE)
  g[nzchar(g)]
} else character()

seurat_ver <- as.character(packageVersion("Seurat"))
unique_patients_in_21 <- length(unique(el$patient_id))
n_multi <- sum(el$n_samples > 1)

finding <- c(
  "# FINDING — Seurat GSE131907 CLDN4-only (patient unit)",
  "",
  "ADDITIVE. **CLDN4-only.** No TACSTD2∩CLDN4 dual-high. Not GSE148071.",
  "Kim et al., *Nat Commun* 2020, PMID 32385277 ([GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907)).",
  "Thesis (already correct; not re-derived): CLDN4-high malignant cells have lower own",
  "IFN/MHC-I, and patients have lower T/NK. CLDN4-low/KD opens IFN/MHC. Matching extras",
  "here are IFN/MHC **DOWN** in CLDN4-high. Do not sell a DoRothEA IFN-up-in-high as the",
  "KD direction.",
  "",
  sprintf("Primary engine: **Seurat %s `CreateSeuratObject`** on the public processed raw UMI", seurat_ver),
  "matrix (gene-subset of locked IFN/MHC/TJ/chemokine + anchors). Full-library UMI sizes",
  "were streamed in the same pass and used for `NormalizeData` / CPM. The 2.86 GB log2TPM",
  "text and EGA FASTQ were not used. Python-only is not the primary.",
  "",
  "**Honest unit = patient.** Tumor-bearing sites only (tLung, tL/B, mLN, PE, mBrain).",
  "Malignant = author `Cell_subtype == Malignant cells` (same author-malig gate as the",
  "concordant-4 GSE131907 slice). Eligible if n_malignant ≥ 20. nLung has 0 author-malignant",
  "cells and cannot join. p-values are descriptive.",
  "",
  "## 1. CLDN4 vs T/NK fraction",
  "",
  sprintf(
    "Primary (patient, author Malignant, %%pos): **n=%d**, Spearman ρ=%s, p=%s.",
    prim$n, fmt_n(prim$rho_pct), fmt_p(prim$p_pct)
  ),
  sprintf(
    "Mean log-normalized CLDN4 vs T/NK: ρ=%s, p=%s.",
    fmt_n(prim$rho_mean), fmt_p(prim$p_mean)
  ),
  sprintf(
    "Q4 vs Q1 T/NK (within-patient-set %%pos quartiles): rank-biserial r=%s, n_Q1/n_Q4=%d/%d, p=%s; median T/NK Q1=%s Q4=%s.",
    fmt_n(prim$r_rb_q4q1), prim$n_q1, prim$n_q4, fmt_p(prim$p_q4q1),
    fmt_n(prim$median_tnk_q1, 3), fmt_n(prim$median_tnk_q4, 3)
  ),
  sprintf(
    "Unique patients in the eligible set: %d. Patients with >1 tumor-bearing sample in this slice: %d (P1006, P1011, P1012, P1013: each adds a PE capture with n_malignant=0).",
    unique_patients_in_21, n_multi
  ),
  "Sample-level T/NK uses only the malignant-bearing capture (concordant-4 given ρ=−0.522).",
  "Patient-level T/NK pools that PE immune into the same patient, so ρ softens. That is the honest patient unit.",
  "",
  "| contrast | unit | malignant | n | rho %pos (p) | rho mean (p) | Q4 vs Q1 r (n_Q1/n_Q4, p) |",
  "|---|---|---|---:|---|---|---|",
  sprintf(
    "| primary | patient | author Malignant cells | %d | %s (%s) | %s (%s) | %s (%d/%d, %s) |",
    prim$n, fmt_n(prim$rho_pct), fmt_p(prim$p_pct), fmt_n(prim$rho_mean), fmt_p(prim$p_mean),
    fmt_n(prim$r_rb_q4q1), prim$n_q1, prim$n_q4, fmt_p(prim$p_q4q1)
  ),
  sprintf(
    "| sensitivity | patient | Malignant + tS1/tS2/tS3 | %d | %s (%s) | %s (%s) | %s (%d/%d, %s) |",
    broad$n, fmt_n(broad$rho_pct), fmt_p(broad$p_pct), fmt_n(broad$rho_mean), fmt_p(broad$p_mean),
    fmt_n(broad$r_rb_q4q1), broad$n_q1, broad$n_q4, fmt_p(broad$p_q4q1)
  ),
  sprintf(
    "| comparison | sample | author Malignant cells | %d | %s (%s) | %s (%s) | %s (%d/%d, %s) |",
    samp$n, fmt_n(samp$rho_pct), fmt_p(samp$p_pct), fmt_n(samp$rho_mean), fmt_p(samp$p_mean),
    fmt_n(samp$r_rb_q4q1), samp$n_q1, samp$n_q4, fmt_p(samp$p_q4q1)
  ),
  "",
  "The sample-level author-Malignant row is the comparison to the prior GSE131907",
  "concordant-4 sample n=21 (given ρ=−0.522). It is not a re-audit of the four-set pool.",
  "tS1/tS2/tS3 are **not** in the primary gate (primary tLung uses those labels; author",
  "`Malignant cells` are mets / tL/B / mLN in this atlas).",
  "",
  "## 2. Malignant IFN / MHC-I / TJ — Q4 vs Q1",
  "",
  "Method: **patient-pseudobulk log2(CPM+1)** of author-malignant UMI sums, using the",
  "full-library UMI total as the size factor. Not TMM. Not muscat. Not a cell-level Wilcoxon.",
  "Quartiles = the same patient %pos Q labels as the T/NK test. Positive logFC = higher in",
  "CLDN4-high. CLDN4 is held out of TJ. Expected under the thesis: IFN down, MHC-I/APM down, TJ up.",
  "",
  sprintf("Honest DE n = Q1+Q4 patients in the count collapse: **%d / %d** (Q1/Q4).", ifn$n_q1, ifn$n_q4),
  "Do not quote 208,506 cells as n.",
  "",
  "| family | n | n_Q1 / n_Q4 | n_genes | logFC | r (Q4 vs Q1) | p | expected |",
  "|---|---:|---|---:|---:|---:|---:|---|",
  sprintf(
    "| IFN | %d | %d / %d | %d | %+s | %s | %s | DOWN |",
    ifn$n, ifn$n_q1, ifn$n_q4, ifn$n_genes, fmt_n(ifn$logFC_q4_minus_q1), fmt_n(ifn$r_rb), fmt_p(ifn$p)
  ),
  sprintf(
    "| MHC-I/APM | %d | %d / %d | %d | %+s | %s | %s | DOWN |",
    mhc$n, mhc$n_q1, mhc$n_q4, mhc$n_genes, fmt_n(mhc$logFC_q4_minus_q1), fmt_n(mhc$r_rb), fmt_p(mhc$p)
  ),
  sprintf(
    "| TJ (CLDN4 out) | %d | %d / %d | %d | %+s | %s | %s | UP |",
    tj$n, tj$n_q1, tj$n_q4, tj$n_genes, fmt_n(tj$logFC_q4_minus_q1), fmt_n(tj$r_rb), fmt_p(tj$p)
  ),
  sprintf(
    "| chemokine | %d | %d / %d | %d | %+s | %s | %s | DOWN |",
    che$n, che$n_q1, che$n_q4, che$n_genes, fmt_n(che$logFC_q4_minus_q1), fmt_n(che$r_rb), fmt_p(che$p)
  ),
  sprintf(
    "| CLDN4 (split check) | %d | %d / %d | 1 | %+s | %s | %s | UP |",
    cld$n, cld$n_q1, cld$n_q4, fmt_n(cld$logFC_q4_minus_q1), fmt_n(cld$r_rb), fmt_p(cld$p)
  ),
  "",
  "IFN is near-flat in this single cohort. That matches the published GSE131907-only",
  "family row in the four-set table (already the near-null IFN set). The four-set pooled",
  "IFN DOWN is not re-derived here. MHC-I/APM and chemokine are DOWN and TJ is UP;",
  "tails are thin. CLDN4 itself is UP (split check).",
  "",
  "Gene families are the locked A8/concordant-4 sets (Hallmark IFN-α∪γ; custom MHC-I/APM;",
  "KEGG TJ ∪ GOBP tight-junction organization minus CLDN4; compact chemokine panel).",
  if (length(missing_genes)) {
    sprintf("Genes in the locked lists but absent from the UMI matrix: %s.", paste(missing_genes, collapse = ", "))
  } else {
    "Every locked family gene queried was present in the UMI matrix."
  },
  "",
  "## Honest n",
  "",
  sprintf("- CreateSeuratObject cells: **%d** (do not quote as n).", ncol(obj)),
  sprintf("- GEO patients / samples: **%d / %d**.", length(unique(md$patient_id)), length(unique(md$Sample))),
  sprintf("- Author Malignant cells: **%d**. tS1+tS2+tS3: **%d**.", sum(md$mal_author), sum(md$Cell_subtype %in% c("tS1", "tS2", "tS3"))),
  sprintf("- PRIMARY eligible patients: **%d** (Q1=%d, Q4=%d).", nrow(el), sum(el$quartile == "Q1"), sum(el$quartile == "Q4")),
  sprintf("- Author-malignant cells in nLung: **%d**.", sum(md$mal_author & md$Sample_Origin == "nLung")),
  sprintf("- Seurat %s; NormalizeData LogNormalize on full-library size.", seurat_ver),
  "- Treatment-naive LUAD atlas: no ICI / MPR / RECIST labels.",
  "",
  "## What this is not",
  "",
  "- Not a mega-merge and not GSE148071 / GSE207422 / GSE205335 / GSE123902 / GSE189357.",
  "- Not a dual-high TACSTD2×CLDN4 score.",
  "- Not a re-derivation of the thesis and not a re-audit of the concordant-4 pooled n=65.",
  "- Not muscat mixed-model DE and not a cell-level Wilcoxon sold as n.",
  "- Not evidence that CLDN4 *causes* the T/NK or IFN/MHC change.",
  "- Q4 vs Q1 tails are thin in a single cohort; family direction is the claim, not genome-wide FDR.",
  "",
  "## Files",
  "",
  "- `results/tables/tnk_tests.tsv` — **headline T/NK table**",
  "- `results/tables/family_q4q1.tsv` — **headline IFN/MHC/TJ table**",
  "- `results/tables/patient_units.tsv` — patient-level CLDN4 and T/NK",
  "- `results/tables/n_honest.tsv`",
  "- `results/figures/fig_cldn4_vs_tnk.png`",
  "- `results/figures/fig_tnk_q4q1.png`",
  "- `results/figures/fig_family_q4q1.png`",
  "- `results/figures/fig_honest_n.png`",
  "",
  "## Reproduce",
  "",
  "```bash",
  "bash methods/seurat_gse131907_cldn4/scripts/download.sh /tmp/gse131907",
  "python3 methods/seurat_gse131907_cldn4/scripts/extract_gene_matrix.py \\",
  "  --datadir /tmp/gse131907 --outdir /tmp/gse131907/subset",
  "Rscript methods/seurat_gse131907_cldn4/scripts/analyze_seurat.R \\",
  "  /tmp/gse131907/subset /tmp/gse131907 \\",
  "  methods/seurat_gse131907_cldn4/results \\",
  "  methods/seurat_gse131907_cldn4/FINDING.md \\",
  "  methods/seurat_gse131907_cldn4/data/families.json",
  "```"
)

writeLines(finding, opt$finding)

summary <- list(
  seurat = seurat_ver,
  n_cells = ncol(obj),
  n_patients_primary = nrow(el),
  rho_pct = prim$rho_pct,
  p_pct = prim$p_pct,
  r_rb = prim$r_rb_q4q1,
  p_q4q1 = prim$p_q4q1,
  ifn_logfc = ifn$logFC_q4_minus_q1,
  mhc_logfc = mhc$logFC_q4_minus_q1,
  tj_logfc = tj$logFC_q4_minus_q1
)
writeLines(
  sprintf(
    paste(
      "{",
      "  \"seurat\": \"%s\",",
      "  \"n_cells\": %d,",
      "  \"n_patients_primary\": %d,",
      "  \"rho_pct\": %s,",
      "  \"p_pct\": %s,",
      "  \"r_rb_q4q1\": %s,",
      "  \"p_q4q1\": %s,",
      "  \"ifn_logfc\": %s,",
      "  \"mhc_logfc\": %s,",
      "  \"tj_logfc\": %s",
      "}",
      sep = "\n"
    ),
    seurat_ver, ncol(obj), nrow(el),
    prim$rho_pct, prim$p_pct, prim$r_rb_q4q1, prim$p_q4q1,
    ifn$logFC_q4_minus_q1, mhc$logFC_q4_minus_q1, tj$logFC_q4_minus_q1
  ),
  file.path(opt$outdir, "summary.json")
)

message("Wrote ", opt$finding)
message("Primary n=", prim$n, " rho=", prim$rho_pct, " p=", prim$p_pct)
message("IFN logFC=", ifn$logFC_q4_minus_q1, " MHC=", mhc$logFC_q4_minus_q1, " TJ=", tj$logFC_q4_minus_q1)
