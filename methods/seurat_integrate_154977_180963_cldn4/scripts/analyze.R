#!/usr/bin/env Rscript
# ADDITIVE pair: GSE154977 + GSE180963 only. R + Seurat + Harmony.
# Cldn4-only. Honest unit = mouse. Dataset covariate + leave-one-out.
# Not the triple. No GSE267321. No human. No private 8-KL.
# No FACS-epi-only primary (154977 is FACS CD45− by design; T/NK is no-go there).
# No Python-only primary.
#
# Usage:
#   Rscript methods/seurat_integrate_154977_180963_cldn4/scripts/analyze.R \
#     --data /tmp/pair_154977_180963 --out methods/seurat_integrate_154977_180963_cldn4

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(hdf5r)
  library(Matrix)
  library(ggplot2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- "/tmp/pair_154977_180963"
out_dir <- "methods/seurat_integrate_154977_180963_cldn4"
i <- 1
while (i <= length(args)) {
  if (args[[i]] == "--data" && i < length(args)) {
    data_dir <- args[[i + 1]]
    i <- i + 2
  } else if (args[[i]] == "--out" && i < length(args)) {
    out_dir <- args[[i + 1]]
    i <- i + 2
  } else {
    stop("unknown arg: ", args[[i]])
  }
}

tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
obj_dir <- file.path(out_dir, "objects")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1)

# Locked compact mouse sets (same family as GSE180963 / public_kl_vs_kp).
# One set for both datasets so pair scores are comparable.
# Cldn4 is the target, never an epithelial caller. Tacstd2 is inventory only.
CLDN4 <- "Cldn4"
T_NK_CALL <- c("Cd3d", "Cd3e", "Cd3g", "Cd8a", "Nkg7", "Ncr1", "Klrb1c")
T_NK_SCORE <- c(
  "Cd3d", "Cd3e", "Cd3g", "Cd2", "Cd8a", "Cd8b1", "Cd4", "Nkg7",
  "Gzma", "Gzmb", "Prf1", "Klrb1c", "Ncr1", "Klrd1", "Klrc1", "Ifng"
)
IFN <- c(
  "Stat1", "Stat2", "Irf1", "Irf7", "Irf9", "Isg15",
  "Ifit1", "Ifit2", "Ifit3", "Mx1", "Oasl2", "Rsad2", "Ifih1", "Ddx58", "Ifnb1"
)
MHC <- c(
  "B2m", "H2-K1", "H2-D1", "H2-Q4", "H2-Q6", "H2-Q7",
  "H2-Aa", "H2-Ab1", "H2-Eb1", "Tap1", "Tap2", "Psmb8", "Psmb9", "Nlrc5", "Ciita"
)
EPI_CORE <- c("Epcam", "Cdh1", "Krt8", "Krt18", "Krt19")
HOST_LUNG <- c("Sftpc", "Scgb1a1", "Ager")
AUDIT <- c("Ptprc", "Stk11", "Tacstd2", "Cldn18", "Nkx2-1")

MIN_N_SPEARMAN <- 4L
MIN_N_Q4 <- 8L

present_in <- function(genes, universe) genes[genes %in% universe]

pos_any <- function(mat, genes) {
  genes <- present_in(genes, rownames(mat))
  if (length(genes) == 0) return(rep(FALSE, ncol(mat)))
  as.numeric(Matrix::colSums(mat[genes, , drop = FALSE] > 0)) > 0
}

mean_mod <- function(mat, genes) {
  genes <- present_in(genes, rownames(mat))
  if (length(genes) == 0) return(rep(NA_real_, ncol(mat)))
  as.numeric(Matrix::colMeans(mat[genes, , drop = FALSE]))
}

spearman_safe <- function(x, y, min_n = MIN_N_SPEARMAN, exact = TRUE) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  out <- list(n = n, rho = NA_real_, p = NA_real_, p_exact = NA_real_,
              p_asymp = NA_real_, usable = FALSE)
  if (n < min_n) return(out)
  ct_a <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  out$rho <- unname(ct_a$estimate)
  out$p_asymp <- unname(ct_a$p.value)
  if (exact && n <= 10) {
    ct_e <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = TRUE))
    out$p_exact <- unname(ct_e$p.value)
    out$p <- out$p_exact
  } else {
    out$p <- out$p_asymp
  }
  out$usable <- TRUE
  out
}

fmt <- function(x, d = 3) {
  if (length(x) == 0 || is.null(x) || (is.numeric(x) && !is.finite(x[1]))) return("—")
  formatC(as.numeric(x[1]), format = "f", digits = d)
}

read_coo_h5 <- function(path, genes, cells) {
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  ii <- as.numeric(h5[["i"]][])
  jj <- as.numeric(h5[["j"]][])
  vv <- as.numeric(h5[["v"]][])
  if (min(ii) == 0 || min(jj) == 0) {
    ii <- ii + 1
    jj <- jj + 1
  }
  n_g <- length(genes)
  n_c <- length(cells)
  if (max(ii) > n_g || max(jj) > n_c) {
    stop(sprintf("COO index out of range: max i=%s (n_genes=%s) max j=%s (n_cells=%s)",
                 max(ii), n_g, max(jj), n_c))
  }
  sparseMatrix(
    i = as.integer(ii),
    j = as.integer(jj),
    x = vv,
    dims = c(n_g, n_c),
    dimnames = list(genes, cells)
  )
}

# ---------------------------------------------------------------------------
# Load GSE154977 (KP 30w 10x; FACS CD45− tumor). Not FACS-epi-only as a pair.
# ---------------------------------------------------------------------------
d154 <- file.path(data_dir, "GSE154977")
need154 <- file.path(d154, c(
  "GSE154977_mmLung10x_cis_dSp_rawCount.h5",
  "GSE154977_mmLung10x_cis_smpTable.csv.gz",
  "GSE154977_mmLung10x_cis_geneTable.csv.gz"
))
if (any(!file.exists(need154))) {
  stop("Missing GSE154977 files. Run scripts/download.sh ", data_dir)
}

smp <- read.csv(file.path(d154, "GSE154977_mmLung10x_cis_smpTable.csv.gz"),
                stringsAsFactors = FALSE)
gene_tbl <- read.csv(file.path(d154, "GSE154977_mmLung10x_cis_geneTable.csv.gz"),
                     stringsAsFactors = FALSE)
if (anyDuplicated(gene_tbl$geneID)) {
  gene_tbl$geneID <- make.unique(gene_tbl$geneID)
}
genes154 <- gene_tbl$geneID
cells154 <- smp$sampleID
stopifnot(length(unique(cells154)) == length(cells154))

message("Reading GSE154977 COO h5 ...")
counts154 <- read_coo_h5(
  file.path(d154, "GSE154977_mmLung10x_cis_dSp_rawCount.h5"),
  genes154, cells154
)

meta154 <- data.frame(row.names = cells154, sampleID = cells154, stringsAsFactors = FALSE)
meta154$barcode <- smp$barcode[match(cells154, smp$sampleID)]
meta154$library <- sub("_id-.*$", "", meta154$sampleID)
meta154$mouse <- sub("_PT$", "", meta154$library)
meta154$treatment <- ifelse(grepl("Cis72", meta154$library), "Cis72", "ND")
meta154$genotype <- "KP"
meta154$arm <- "KP"
meta154$dataset <- "GSE154977"
meta154$design <- "FACS_CD45neg_tumor"
meta154$tnk_honest <- FALSE
meta154$gsm <- NA_character_
meta154$gsm[meta154$library == "KP_30w_ND_m3_PT"] <- "GSM4685281"
meta154$gsm[meta154$library == "KP_30w_ND_m4_PT"] <- "GSM4685282"
meta154$gsm[meta154$library == "KP_30w_Cis72_m5_PT"] <- "GSM4685283"
meta154$gsm[meta154$library == "KP_30w_Cis72_m6_PT"] <- "GSM4685284"

message("CreateSeuratObject GSE154977")
obj154 <- CreateSeuratObject(
  counts = counts154,
  meta.data = meta154,
  project = "GSE154977",
  min.cells = 0,
  min.features = 0
)
rm(counts154)
gc()

# ---------------------------------------------------------------------------
# Load GSE180963 (K vs KL whole-tumor 10x MTX).
# ---------------------------------------------------------------------------
d180 <- file.path(data_dir, "GSE180963")
need180 <- file.path(d180, c(
  "K/matrix.mtx", "K/genes.tsv", "K/barcodes.tsv",
  "KL/matrix.mtx", "KL/genes.tsv", "KL/barcodes.tsv"
))
if (any(!file.exists(need180))) {
  stop("Missing GSE180963 MTX. Run scripts/download.sh ", data_dir)
}

read_180 <- function(label, gsm, genotype, arm) {
  counts <- Read10X(data.dir = file.path(d180, label), gene.column = 2, unique.features = TRUE)
  obj <- CreateSeuratObject(counts = counts, project = label, min.cells = 0, min.features = 0)
  obj$mouse <- label
  obj$gsm <- gsm
  obj$genotype <- genotype
  obj$arm <- arm
  obj$dataset <- "GSE180963"
  obj$design <- "whole_tumor_mixed_library"
  obj$treatment <- "ND"
  obj$library <- label
  obj$tnk_honest <- TRUE
  obj$orig.ident <- label
  obj
}

message("CreateSeuratObject GSE180963 K")
obj_k <- read_180("K", "GSM5481386", "KrasG12D/+", "K")
message("CreateSeuratObject GSE180963 KL")
obj_kl <- read_180("KL", "GSM5481387", "KrasG12D/+;Lkb1fl/fl", "KL")
obj180 <- merge(obj_k, y = obj_kl, add.cell.ids = c("K", "KL"), project = "GSE180963")
rm(obj_k, obj_kl)
gc()
if (packageVersion("SeuratObject") >= "5.0.0") {
  obj180 <- JoinLayers(obj180)
}

# ---------------------------------------------------------------------------
# Pair merge on shared genes. Harmony on dataset.
# ---------------------------------------------------------------------------
shared <- intersect(rownames(obj154), rownames(obj180))
message(sprintf("shared genes: %s (154977=%s, 180963=%s)",
                length(shared), nrow(obj154), nrow(obj180)))
if (!(CLDN4 %in% shared)) {
  stop("Cldn4 missing from shared gene universe — cannot score Cldn4-only pair.")
}

obj154 <- subset(obj154, features = shared)
obj180 <- subset(obj180, features = shared)

# Prefix cell barcodes so merge is unique.
obj154 <- RenameCells(obj154, add.cell.id = "GSE154977")

message("merge pair")
obj <- merge(obj154, y = obj180, add.cell.ids = NULL, project = "pair_154977_180963")
rm(obj154, obj180)
gc()
if (packageVersion("SeuratObject") >= "5.0.0") {
  obj <- JoinLayers(obj)
}

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
# Author matrices are already QC'd. Light pair filter only.
qc_keep <- obj$nFeature_RNA >= 200 & obj$percent.mt < 25
message(sprintf("QC keep %s / %s", sum(qc_keep), ncol(obj)))
obj$qc_pass <- qc_keep
obj <- subset(obj, subset = qc_pass)

message("Normalize / HVG / PCA / Harmony(dataset)")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# Harmony 2.x Seurat method. Dataset is the integration batch.
obj <- RunHarmony(
  obj,
  group.by.vars = "dataset",
  reduction.use = "pca",
  dims.use = 1:20,
  reduction.save = "harmony",
  verbose = TRUE
)

harm_name <- if ("harmony" %in% Reductions(obj)) "harmony" else {
  cand <- grep("harmony", Reductions(obj), value = TRUE, ignore.case = TRUE)
  if (!length(cand)) stop("RunHarmony did not write a harmony reduction. Reductions: ",
                          paste(Reductions(obj), collapse = ", "))
  cand[[1]]
}
message("harmony reduction: ", harm_name)

obj <- FindNeighbors(obj, reduction = harm_name, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
obj <- RunUMAP(obj, reduction = harm_name, dims = 1:20, verbose = FALSE)

# ---------------------------------------------------------------------------
# Marker compartments. Cldn4 is never a caller.
# ---------------------------------------------------------------------------
universe <- rownames(obj)
inv_rows <- lapply(unique(c(CLDN4, T_NK_SCORE, IFN, MHC, EPI_CORE, HOST_LUNG, AUDIT)), function(g) {
  data.frame(gene = g, present = g %in% universe, stringsAsFactors = FALSE)
})
inv <- do.call(rbind, inv_rows)
write.table(inv, file.path(tab_dir, "gene_inventory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
logn <- GetAssayData(obj, assay = "RNA", layer = "data")

epcam <- if ("Epcam" %in% universe) as.numeric(counts["Epcam", ] > 0) else 0
cdh1 <- if ("Cdh1" %in% universe) as.numeric(counts["Cdh1", ] > 0) else 0
krt8 <- if ("Krt8" %in% universe) as.numeric(counts["Krt8", ] > 0) else 0
krt18 <- if ("Krt18" %in% universe) as.numeric(counts["Krt18", ] > 0) else 0
krt19 <- if ("Krt19" %in% universe) as.numeric(counts["Krt19", ] > 0) else 0
cldn18 <- if ("Cldn18" %in% universe) as.numeric(counts["Cldn18", ] > 0) else 0
ptprc <- if ("Ptprc" %in% universe) as.numeric(counts["Ptprc", ] > 0) else 0
sftpc <- if ("Sftpc" %in% universe) as.numeric(counts["Sftpc", ] > 0) else 0
struct <- (cdh1 > 0) | (krt8 > 0) | (krt18 > 0) | (krt19 > 0) | (cldn18 > 0)
epi_tight <- (epcam > 0) & (struct > 0) & (ptprc == 0)
tnk_mark <- pos_any(counts, T_NK_CALL)
is_epi <- epi_tight
is_tnk <- tnk_mark & !is_epi

obj$is_epi <- is_epi
obj$is_tnk <- is_tnk
obj$sftpc_pos <- sftpc > 0
obj$ptprc_pos <- ptprc > 0
obj$compartment <- ifelse(is_epi, "epithelial",
                   ifelse(is_tnk, "T/NK", "other"))

obj$Cldn4_count <- as.numeric(counts[CLDN4, ])
obj$Cldn4_logn <- as.numeric(logn[CLDN4, ])
obj$Cldn4_pos <- obj$Cldn4_count > 0
obj$ifn_score <- mean_mod(logn, IFN)
obj$mhc_score <- mean_mod(logn, MHC)
obj$tnk_score <- mean_mod(logn, T_NK_SCORE)

ifn_use <- present_in(IFN, universe)
mhc_use <- present_in(MHC, universe)
if (length(ifn_use) >= 2 && length(mhc_use) >= 2) {
  obj <- AddModuleScore(
    obj,
    features = list(IFN = ifn_use, MHC = mhc_use),
    name = c("IFN_seurat", "MHC_seurat"),
    ctrl = min(50, max(10, length(universe) %/% 50)),
    seed = 1
  )
}

# ---------------------------------------------------------------------------
# Mouse-level table (the actual unit). n = 6 mice, not cells.
# ---------------------------------------------------------------------------
md <- slot(obj, "meta.data")
mice <- unique(md$mouse)
mouse_rows <- lapply(mice, function(lab) {
  w <- md$mouse == lab
  we <- w & md$is_epi
  wn <- w & md$is_tnk
  ds <- unique(as.character(md$dataset[w]))
  tnk_ok <- all(as.logical(md$tnk_honest[w]))
  data.frame(
    mouse = lab,
    dataset = ds,
    gsm = unique(as.character(md$gsm[w])),
    genotype = unique(as.character(md$genotype[w])),
    arm = unique(as.character(md$arm[w])),
    treatment = unique(as.character(md$treatment[w])),
    design = unique(as.character(md$design[w])),
    tnk_honest = tnk_ok,
    n_cells = sum(w),
    n_epi = sum(we),
    n_tnk = sum(wn),
    frac_tnk = if (tnk_ok) mean(md$is_tnk[w]) else NA_real_,
    frac_epi = mean(md$is_epi[w]),
    frac_ptprc = mean(md$ptprc_pos[w]),
    Cldn4_all_mean = mean(md$Cldn4_logn[w]),
    Cldn4_all_pctpos = 100 * mean(md$Cldn4_pos[w]),
    Cldn4_epi_mean = if (any(we)) mean(md$Cldn4_logn[we]) else NA_real_,
    Cldn4_epi_pctpos = if (any(we)) 100 * mean(md$Cldn4_pos[we]) else NA_real_,
    IFN_epi_mean = if (any(we)) mean(md$ifn_score[we]) else NA_real_,
    MHC_epi_mean = if (any(we)) mean(md$mhc_score[we]) else NA_real_,
    percent_mt_mean = mean(md$percent.mt[w]),
    stringsAsFactors = FALSE
  )
})
per_mouse <- do.call(rbind, mouse_rows)
per_mouse <- per_mouse[order(per_mouse$dataset, per_mouse$mouse), ]
write.table(per_mouse, file.path(tab_dir, "per_mouse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# Cldn4 vs T/NK (mouse). Only GSE180963 is honest. Pair n = 2 → no-go.
# ---------------------------------------------------------------------------
tnk_mice <- per_mouse[per_mouse$tnk_honest, ]
sp_tnk <- spearman_safe(tnk_mice$Cldn4_epi_mean, tnk_mice$frac_tnk)
tnk_test <- data.frame(
  contrast = "Cldn4_epi vs T/NK fraction",
  n_mice_honest = nrow(tnk_mice),
  n_mice_pair = nrow(per_mouse),
  rho = sp_tnk$rho,
  p_exact = sp_tnk$p_exact,
  p_asymp = sp_tnk$p_asymp,
  usable = sp_tnk$usable,
  note = "T/NK honest only in GSE180963 (whole tumor). GSE154977 is FACS CD45− tumor — design no-go, not used. Pair T/NK n=2 < lock (4).",
  stringsAsFactors = FALSE
)
write.table(tnk_test, file.path(tab_dir, "cldn4_vs_tnk.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# Cldn4 vs epithelial IFN/MHC (mouse). n = 6. Dataset covariate + LOO.
# ---------------------------------------------------------------------------
epi_ok <- is.finite(per_mouse$Cldn4_epi_mean) &
  is.finite(per_mouse$IFN_epi_mean) &
  is.finite(per_mouse$MHC_epi_mean) &
  per_mouse$n_epi >= 10
pm <- per_mouse[epi_ok, ]
stopifnot(nrow(pm) == nrow(per_mouse))

sp_ifn <- spearman_safe(pm$Cldn4_epi_mean, pm$IFN_epi_mean)
sp_mhc <- spearman_safe(pm$Cldn4_epi_mean, pm$MHC_epi_mean)

# Residualize on dataset (factor), then Spearman — the covariate-adjusted rank test.
resid_on_dataset <- function(y, dataset) {
  d <- data.frame(y = y, dataset = factor(dataset))
  as.numeric(residuals(lm(y ~ dataset, data = d)))
}
r_cldn <- resid_on_dataset(pm$Cldn4_epi_mean, pm$dataset)
r_ifn <- resid_on_dataset(pm$IFN_epi_mean, pm$dataset)
r_mhc <- resid_on_dataset(pm$MHC_epi_mean, pm$dataset)
sp_ifn_adj <- spearman_safe(r_cldn, r_ifn)
sp_mhc_adj <- spearman_safe(r_cldn, r_mhc)

# OLS with dataset covariate (descriptive; n=6).
fit_ifn <- lm(IFN_epi_mean ~ Cldn4_epi_mean + dataset, data = pm)
fit_mhc <- lm(MHC_epi_mean ~ Cldn4_epi_mean + dataset, data = pm)
s_ifn <- summary(fit_ifn)
s_mhc <- summary(fit_mhc)

lm_tab <- data.frame(
  outcome = c("IFN_epi", "MHC_epi"),
  n = c(nrow(pm), nrow(pm)),
  beta_Cldn4 = c(coef(fit_ifn)[["Cldn4_epi_mean"]], coef(fit_mhc)[["Cldn4_epi_mean"]]),
  se_Cldn4 = c(s_ifn$coefficients["Cldn4_epi_mean", "Std. Error"],
               s_mhc$coefficients["Cldn4_epi_mean", "Std. Error"]),
  p_Cldn4 = c(s_ifn$coefficients["Cldn4_epi_mean", "Pr(>|t|)"],
              s_mhc$coefficients["Cldn4_epi_mean", "Pr(>|t|)"]),
  beta_dataset180963 = c(coef(fit_ifn)[["datasetGSE180963"]], coef(fit_mhc)[["datasetGSE180963"]]),
  p_dataset = c(s_ifn$coefficients["datasetGSE180963", "Pr(>|t|)"],
                s_mhc$coefficients["datasetGSE180963", "Pr(>|t|)"]),
  r2 = c(s_ifn$r.squared, s_mhc$r.squared),
  note = "OLS at mouse unit; dataset factor; n=6 is thin; not a powered claim",
  stringsAsFactors = FALSE
)
write.table(lm_tab, file.path(tab_dir, "dataset_covariate_lm.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Within-dataset (154977 n=4 usable; 180963 n=2 locked off).
within_rows <- list()
for (ds in c("GSE154977", "GSE180963")) {
  sub <- pm[pm$dataset == ds, ]
  si <- spearman_safe(sub$Cldn4_epi_mean, sub$IFN_epi_mean)
  sm <- spearman_safe(sub$Cldn4_epi_mean, sub$MHC_epi_mean)
  within_rows[[length(within_rows) + 1]] <- data.frame(
    dataset = ds, family = "IFN", n = si$n, rho = si$rho,
    p_exact = si$p_exact, p_asymp = si$p_asymp, usable = si$usable,
    stringsAsFactors = FALSE
  )
  within_rows[[length(within_rows) + 1]] <- data.frame(
    dataset = ds, family = "MHC", n = sm$n, rho = sm$rho,
    p_exact = sm$p_exact, p_asymp = sm$p_asymp, usable = sm$usable,
    stringsAsFactors = FALSE
  )
}
within_tab <- do.call(rbind, within_rows)
write.table(within_tab, file.path(tab_dir, "within_dataset_spearman.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Leave-one-mouse-out.
loo_rows <- list()
for (drop in pm$mouse) {
  sub <- pm[pm$mouse != drop, ]
  si <- spearman_safe(sub$Cldn4_epi_mean, sub$IFN_epi_mean)
  sm <- spearman_safe(sub$Cldn4_epi_mean, sub$MHC_epi_mean)
  si_a <- spearman_safe(resid_on_dataset(sub$Cldn4_epi_mean, sub$dataset),
                        resid_on_dataset(sub$IFN_epi_mean, sub$dataset))
  sm_a <- spearman_safe(resid_on_dataset(sub$Cldn4_epi_mean, sub$dataset),
                        resid_on_dataset(sub$MHC_epi_mean, sub$dataset))
  loo_rows[[length(loo_rows) + 1]] <- data.frame(
    left_out = drop,
    left_out_dataset = pm$dataset[pm$mouse == drop],
    n = nrow(sub),
    IFN_rho = si$rho, IFN_p_exact = si$p_exact, IFN_usable = si$usable,
    MHC_rho = sm$rho, MHC_p_exact = sm$p_exact, MHC_usable = sm$usable,
    IFN_rho_dataset_adj = si_a$rho, IFN_adj_usable = si_a$usable,
    MHC_rho_dataset_adj = sm_a$rho, MHC_adj_usable = sm_a$usable,
    stringsAsFactors = FALSE
  )
}
loo_tab <- do.call(rbind, loo_rows)
write.table(loo_tab, file.path(tab_dir, "loo_mouse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Leave-one-dataset-out (still mouse unit).
lod_rows <- list()
for (ds in c("GSE154977", "GSE180963")) {
  sub <- pm[pm$dataset != ds, ]
  si <- spearman_safe(sub$Cldn4_epi_mean, sub$IFN_epi_mean)
  sm <- spearman_safe(sub$Cldn4_epi_mean, sub$MHC_epi_mean)
  lod_rows[[length(lod_rows) + 1]] <- data.frame(
    left_out_dataset = ds,
    kept_dataset = unique(sub$dataset),
    n = nrow(sub),
    IFN_rho = si$rho, IFN_p_exact = si$p_exact, IFN_usable = si$usable,
    MHC_rho = sm$rho, MHC_p_exact = sm$p_exact, MHC_usable = sm$usable,
    note = if (nrow(sub) < MIN_N_SPEARMAN) "locked: n < 4" else "usable",
    stringsAsFactors = FALSE
  )
}
lod_tab <- do.call(rbind, lod_rows)
write.table(lod_tab, file.path(tab_dir, "loo_dataset.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ND-only sensitivity (drop Cis72). n = 2 ND KP + 2 GSE180963 = 4.
nd <- pm[pm$treatment == "ND", ]
sp_ifn_nd <- spearman_safe(nd$Cldn4_epi_mean, nd$IFN_epi_mean)
sp_mhc_nd <- spearman_safe(nd$Cldn4_epi_mean, nd$MHC_epi_mean)
sp_ifn_nd_adj <- if (nrow(nd) >= MIN_N_SPEARMAN && length(unique(nd$dataset)) > 1) {
  spearman_safe(resid_on_dataset(nd$Cldn4_epi_mean, nd$dataset),
                resid_on_dataset(nd$IFN_epi_mean, nd$dataset))
} else {
  list(n = nrow(nd), rho = NA_real_, p_exact = NA_real_, usable = FALSE)
}
sp_mhc_nd_adj <- if (nrow(nd) >= MIN_N_SPEARMAN && length(unique(nd$dataset)) > 1) {
  spearman_safe(resid_on_dataset(nd$Cldn4_epi_mean, nd$dataset),
                resid_on_dataset(nd$MHC_epi_mean, nd$dataset))
} else {
  list(n = nrow(nd), rho = NA_real_, p_exact = NA_real_, usable = FALSE)
}

family_tab <- rbind(
  data.frame(
    cohort = "pair_6_mice", family = "IFN", adjustment = "none",
    n = sp_ifn$n, rho = sp_ifn$rho, p_exact = sp_ifn$p_exact, p_asymp = sp_ifn$p_asymp,
    usable = sp_ifn$usable, stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "pair_6_mice", family = "MHC", adjustment = "none",
    n = sp_mhc$n, rho = sp_mhc$rho, p_exact = sp_mhc$p_exact, p_asymp = sp_mhc$p_asymp,
    usable = sp_mhc$usable, stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "pair_6_mice", family = "IFN", adjustment = "dataset_residual",
    n = sp_ifn_adj$n, rho = sp_ifn_adj$rho, p_exact = sp_ifn_adj$p_exact, p_asymp = sp_ifn_adj$p_asymp,
    usable = sp_ifn_adj$usable, stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "pair_6_mice", family = "MHC", adjustment = "dataset_residual",
    n = sp_mhc_adj$n, rho = sp_mhc_adj$rho, p_exact = sp_mhc_adj$p_exact, p_asymp = sp_mhc_adj$p_asymp,
    usable = sp_mhc_adj$usable, stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "ND_only_4_mice", family = "IFN", adjustment = "none",
    n = sp_ifn_nd$n, rho = sp_ifn_nd$rho, p_exact = sp_ifn_nd$p_exact, p_asymp = sp_ifn_nd$p_asymp,
    usable = sp_ifn_nd$usable, stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "ND_only_4_mice", family = "MHC", adjustment = "none",
    n = sp_mhc_nd$n, rho = sp_mhc_nd$rho, p_exact = sp_mhc_nd$p_exact, p_asymp = sp_mhc_nd$p_asymp,
    usable = sp_mhc_nd$usable, stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "ND_only_4_mice", family = "IFN", adjustment = "dataset_residual",
    n = sp_ifn_nd_adj$n, rho = sp_ifn_nd_adj$rho, p_exact = sp_ifn_nd_adj$p_exact,
    p_asymp = NA_real_, usable = isTRUE(sp_ifn_nd_adj$usable), stringsAsFactors = FALSE
  ),
  data.frame(
    cohort = "ND_only_4_mice", family = "MHC", adjustment = "dataset_residual",
    n = sp_mhc_nd_adj$n, rho = sp_mhc_nd_adj$rho, p_exact = sp_mhc_nd_adj$p_exact,
    p_asymp = NA_real_, usable = isTRUE(sp_mhc_nd_adj$usable), stringsAsFactors = FALSE
  )
)
write.table(family_tab, file.path(tab_dir, "family_spearman.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Q4 vs Q1 locked: n=6 < 8.
q_lock <- data.frame(
  test = "Q4 vs Q1 Cldn4",
  n_mice = nrow(pm),
  usable = FALSE,
  note = "locked: needs n_units>=8",
  stringsAsFactors = FALSE
)
write.table(q_lock, file.path(tab_dir, "q4q1_locked.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Cell-level Spearman inside epithelium is exploratory / pseudoreplication.
cell_sp <- do.call(rbind, lapply(unique(md$dataset), function(ds) {
  we <- md$dataset == ds & md$is_epi
  si <- spearman_safe(md$Cldn4_logn[we], md$ifn_score[we], min_n = 10, exact = FALSE)
  sm <- spearman_safe(md$Cldn4_logn[we], md$mhc_score[we], min_n = 10, exact = FALSE)
  data.frame(
    subset = paste0(ds, "_epi_cells"),
    n_cells = si$n,
    IFN_rho = si$rho, IFN_p = si$p_asymp,
    MHC_rho = sm$rho, MHC_p = sm$p_asymp,
    note = "cell-level; pseudoreplication; not the unit",
    stringsAsFactors = FALSE
  )
}))
write.table(cell_sp, file.path(tab_dir, "cldn4_vs_ifn_mhc_cells_exploratory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

comp <- as.data.frame(table(dataset = md$dataset, mouse = md$mouse, compartment = md$compartment),
                      stringsAsFactors = FALSE)
comp <- comp[comp$Freq > 0, ]
write.table(comp, file.path(tab_dir, "compartment_by_mouse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

honest <- data.frame(
  item = c(
    "GEO series (pair, not triple)",
    "mice (unit)",
    "GSE154977 KP mice",
    "GSE180963 mice",
    "T/NK-honest mice",
    "IFN/MHC mice (tight epi, n_epi>=10)",
    "cells after pair QC (not the unit)",
    "shared genes",
    "Cldn4-positive cells",
    "tight epithelial cells",
    "T/NK cells (marker)",
    "Q4 vs Q1",
    "GSE267321",
    "human cohorts",
    "private 8-KL mice",
    "FACS-epi-only as pair primary",
    "Harmony batch"
  ),
  n = c(
    2,
    nrow(per_mouse),
    sum(per_mouse$dataset == "GSE154977"),
    sum(per_mouse$dataset == "GSE180963"),
    sum(per_mouse$tnk_honest),
    nrow(pm),
    ncol(obj),
    length(shared),
    sum(md$Cldn4_pos),
    sum(md$is_epi),
    sum(md$is_tnk),
    0,
    0, 0, 0, 0, 1
  ),
  note = c(
    "GSE154977 + GSE180963 only",
    "4 KP + 1 K + 1 KL",
    "FACS CD45− 30w 10x; 2 ND + 2 Cis72",
    "1 K + 1 KL; mixed 10x library",
    "GSE180963 only; 154977 T/NK is design no-go",
    "all 6 pass",
    "nFeature>=200 & percent.mt<25",
    "intersected symbols before merge",
    "count > 0",
    "Epcam+ structural+ Ptprc−; Cldn4 not a caller",
    "Cd3d/e/g or Cd8a or Nkg7/Ncr1/Klrb1c; tight epi wins",
    "locked: needs n_units>=8",
    "excluded (not the triple)",
    "excluded",
    "excluded",
    "excluded; 154977 design is noted, not sold as T/NK",
    "group.by.vars = dataset"
  ),
  stringsAsFactors = FALSE
)
write.table(honest, file.path(tab_dir, "honest_n.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
message("Plots ...")
p_ds <- DimPlot(obj, group.by = "dataset", reduction = "umap") +
  ggtitle("Harmony UMAP — dataset (pair, not triple)")
p_mouse <- DimPlot(obj, group.by = "mouse", reduction = "umap") +
  ggtitle("Harmony UMAP — mouse (honest unit)")
p_comp <- DimPlot(obj, group.by = "compartment", reduction = "umap") +
  ggtitle("Harmony UMAP — marker compartment")
p_arm <- DimPlot(obj, group.by = "arm", reduction = "umap") +
  ggtitle("Harmony UMAP — arm (KP / K / KL)")
ggsave(file.path(fig_dir, "fig_umap_dataset.png"), p_ds, width = 7.5, height = 5.5, dpi = 140)
ggsave(file.path(fig_dir, "fig_umap_mouse.png"), p_mouse, width = 8, height = 5.5, dpi = 140)
ggsave(file.path(fig_dir, "fig_umap_compartment.png"), p_comp, width = 7.5, height = 5.5, dpi = 140)
ggsave(file.path(fig_dir, "fig_umap_arm.png"), p_arm, width = 7.5, height = 5.5, dpi = 140)

p_feat <- FeaturePlot(obj, features = intersect(c("Cldn4", "Epcam", "Ptprc"), rownames(obj)), ncol = 3)
ggsave(file.path(fig_dir, "fig_feature_cldn4_epcam_ptprc.png"), p_feat, width = 12, height = 4, dpi = 140)

p_n <- ggplot(per_mouse, aes(x = mouse, y = n_cells, fill = dataset)) +
  geom_col() + coord_flip() + theme_bw() +
  labs(title = "Honest n = 6 mice (not cells)", y = "cells after pair QC", x = NULL)
ggsave(file.path(fig_dir, "fig_honest_n.png"), p_n, width = 7.2, height = 3.4, dpi = 140)

p_ifn <- ggplot(pm, aes(x = Cldn4_epi_mean, y = IFN_epi_mean, color = dataset, shape = treatment)) +
  geom_point(size = 3.2) +
  geom_text(aes(label = mouse), vjust = -0.8, size = 3, show.legend = FALSE) +
  theme_bw() +
  labs(title = "Mouse-level Cldn4 vs epithelial IFN (n=6)",
       x = "mean log-norm Cldn4 (tight epi)", y = "mean IFN (tight epi)")
p_mhc <- ggplot(pm, aes(x = Cldn4_epi_mean, y = MHC_epi_mean, color = dataset, shape = treatment)) +
  geom_point(size = 3.2) +
  geom_text(aes(label = mouse), vjust = -0.8, size = 3, show.legend = FALSE) +
  theme_bw() +
  labs(title = "Mouse-level Cldn4 vs epithelial MHC/APM (n=6)",
       x = "mean log-norm Cldn4 (tight epi)", y = "mean MHC/APM (tight epi)")
ggsave(file.path(fig_dir, "fig_cldn4_vs_ifn.png"), p_ifn, width = 7.2, height = 5, dpi = 140)
ggsave(file.path(fig_dir, "fig_cldn4_vs_mhc.png"), p_mhc, width = 7.2, height = 5, dpi = 140)

# T/NK: only two honest points.
p_tnk <- ggplot(per_mouse, aes(x = Cldn4_epi_mean, y = frac_tnk, color = dataset, label = mouse)) +
  geom_point(size = 3.2, na.rm = TRUE) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE, na.rm = TRUE) +
  theme_bw() +
  labs(title = "Cldn4 vs T/NK: honest n = 2 (GSE180963 only)",
       subtitle = "GSE154977 FACS CD45− — T/NK not plotted",
       x = "mean log-norm Cldn4 (tight epi)", y = "T/NK fraction")
ggsave(file.path(fig_dir, "fig_cldn4_vs_tnk.png"), p_tnk, width = 7.2, height = 5, dpi = 140)

p_loo <- ggplot(loo_tab, aes(x = left_out, y = IFN_rho, fill = left_out_dataset)) +
  geom_col() + coord_flip() + theme_bw() +
  geom_hline(yintercept = sp_ifn$rho, linetype = 2) +
  labs(title = "LOO mouse: Spearman Cldn4 vs IFN",
       subtitle = sprintf("full-pair rho=%s (dashed)", fmt(sp_ifn$rho)),
       x = "left-out mouse", y = "Spearman rho on remaining mice")
ggsave(file.path(fig_dir, "fig_loo_ifn.png"), p_loo, width = 7.2, height = 3.6, dpi = 140)

p_loo2 <- ggplot(loo_tab, aes(x = left_out, y = MHC_rho, fill = left_out_dataset)) +
  geom_col() + coord_flip() + theme_bw() +
  geom_hline(yintercept = sp_mhc$rho, linetype = 2) +
  labs(title = "LOO mouse: Spearman Cldn4 vs MHC",
       subtitle = sprintf("full-pair rho=%s (dashed)", fmt(sp_mhc$rho)),
       x = "left-out mouse", y = "Spearman rho on remaining mice")
ggsave(file.path(fig_dir, "fig_loo_mhc.png"), p_loo2, width = 7.2, height = 3.6, dpi = 140)

# Do not save a 25k-cell RDS into git. Write a tiny pointer instead.
writeLines(
  c("Seurat object not committed (pair ~25k cells).",
    "Re-run scripts/analyze.R to rebuild."),
  file.path(obj_dir, "README.md")
)

summary <- list(
  pair = c("GSE154977", "GSE180963"),
  not_included = c("GSE267321", "human", "private_8_KL", "FACS_epi_only_primary", "triple"),
  seurat_version = as.character(packageVersion("Seurat")),
  harmony_version = as.character(packageVersion("harmony")),
  r_version = as.character(getRversion()),
  harmony_reduction = harm_name,
  n_cells = ncol(obj),
  n_shared_genes = length(shared),
  n_mice = nrow(per_mouse),
  n_tnk_honest_mice = sum(per_mouse$tnk_honest),
  mice = per_mouse$mouse,
  tnk = list(n = sp_tnk$n, rho = sp_tnk$rho, usable = sp_tnk$usable),
  IFN = list(n = sp_ifn$n, rho = sp_ifn$rho, p_exact = sp_ifn$p_exact, p_asymp = sp_ifn$p_asymp),
  MHC = list(n = sp_mhc$n, rho = sp_mhc$rho, p_exact = sp_mhc$p_exact, p_asymp = sp_mhc$p_asymp),
  IFN_dataset_adj = list(n = sp_ifn_adj$n, rho = sp_ifn_adj$rho, p_exact = sp_ifn_adj$p_exact),
  MHC_dataset_adj = list(n = sp_mhc_adj$n, rho = sp_mhc_adj$rho, p_exact = sp_mhc_adj$p_exact),
  lm_IFN_beta_Cldn4 = unname(coef(fit_ifn)[["Cldn4_epi_mean"]]),
  lm_IFN_p_Cldn4 = unname(s_ifn$coefficients["Cldn4_epi_mean", "Pr(>|t|)"]),
  lm_MHC_beta_Cldn4 = unname(coef(fit_mhc)[["Cldn4_epi_mean"]]),
  lm_MHC_p_Cldn4 = unname(s_mhc$coefficients["Cldn4_epi_mean", "Pr(>|t|)"]),
  loo_IFN_rho_range = range(loo_tab$IFN_rho, na.rm = TRUE),
  loo_MHC_rho_range = range(loo_tab$MHC_rho, na.rm = TRUE),
  q4q1_usable = FALSE,
  gene_sets = list(IFN = ifn_use, MHC = mhc_use)
)
write_json(summary, file.path(tab_dir, "summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = 6)

sink(file.path(tab_dir, "sessionInfo.txt"))
cat("Seurat", as.character(packageVersion("Seurat")), "\n")
cat("harmony", as.character(packageVersion("harmony")), "\n")
print(sessionInfo())
sink()

message("DONE")
message("mice: ", paste(pm$mouse, collapse = ", "))
message(sprintf("IFN rho=%.3f p_exact=%s; MHC rho=%.3f p_exact=%s",
                sp_ifn$rho, fmt(sp_ifn$p_exact, 3),
                sp_mhc$rho, fmt(sp_mhc$p_exact, 3)))
message(sprintf("dataset-adj IFN rho=%.3f; MHC rho=%.3f",
                sp_ifn_adj$rho, sp_mhc_adj$rho))
message(sprintf("T/NK usable=%s n=%s", sp_tnk$usable, sp_tnk$n))
