#!/usr/bin/env Rscript
# INTEGRATION: GSE165641 + GSE180963 public KL/K GEMM 10x.
# CreateSeuratObject per dataset, then Harmony via Seurat::IntegrateLayers.
# Cldn4-only. Honest unit = mouse after merge. No private 8-KL. No Python primary.
#
# Usage:
#   Rscript methods/seurat_integrate_165641_180963_cldn4/scripts/analyze.R \
#     --data /tmp/geo/work --out methods/seurat_integrate_165641_180963_cldn4

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("STOP: Seurat is not installed.")
}
if (!requireNamespace("harmony", quietly = TRUE)) {
  stop("STOP: harmony is not installed; Seurat IntegrateLayers(HarmonyIntegration) cannot run.")
}

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- "/tmp/geo/work"
out_dir <- "methods/seurat_integrate_165641_180963_cldn4"
i <- 1
while (i <= length(args)) {
  if (args[[i]] == "--data" && i < length(args)) {
    data_dir <- args[[i + 1]]; i <- i + 2
  } else if (args[[i]] == "--out" && i < length(args)) {
    out_dir <- args[[i + 1]]; i <- i + 2
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

# Locked Cldn4-only mouse sets (same family as seurat_gse180963_cldn4 / public_kl_vs_kp).
# Cldn4 is the target, never an epithelial caller. Tacstd2 is inventory only.
# TJ excludes Cldn4 so the module is not circular with the ranking gene.
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
TJ <- c(
  "Cldn1", "Cldn3", "Cldn7", "Cldn18",
  "Tjp1", "Tjp2", "Tjp3", "Ocln", "F11r",
  "Marveld2", "Marveld3", "Cgn", "Crb3"
)
EPI_CORE <- c("Epcam", "Cdh1", "Krt8", "Krt18", "Krt19")
HOST_LUNG <- c("Sftpc", "Scgb1a1", "Ager")
AUDIT <- c("Ptprc", "Stk11", "Tacstd2", "Cldn18", "Nkx2-1")

present_in <- function(genes, universe) genes[genes %in% universe]

fmt <- function(x, d = 3) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  formatC(as.numeric(x), format = "f", digits = d)
}

spearman_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4) {
    return(list(n = sum(ok), rho = NA_real_, p = NA_real_, note = "n<4; Spearman not defined"))
  }
  s <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(n = sum(ok), rho = unname(s$estimate), p = s$p.value, note = "mouse-level Spearman")
}

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

read_10x_dir <- function(d) {
  if (!dir.exists(d)) stop("missing 10x dir: ", d)
  counts <- Read10X(data.dir = d, gene.column = 2, unique.features = TRUE)
  if (is.list(counts)) counts <- counts[[1]]
  counts
}

# ---- 1. CreateSeuratObject PER DATASET (mice kept as metadata) ----
p165_kl1 <- file.path(data_dir, "GSE165641/KL1_count/filtered_feature_bc_matrix")
p165_kl2 <- file.path(data_dir, "GSE165641/KL2_count/filtered_feature_bc_matrix")
p180_k <- file.path(data_dir, "GSE180963/K")
p180_kl <- file.path(data_dir, "GSE180963/KL")

have_165 <- dir.exists(p165_kl1) && dir.exists(p165_kl2)
have_180 <- file.exists(file.path(p180_k, "matrix.mtx")) && file.exists(file.path(p180_kl, "matrix.mtx"))

series_status <- data.frame(
  series = c("GSE165641", "GSE180963"),
  matrix_found = c(have_165, have_180),
  n_libraries = c(if (have_165) 2L else 0L, if (have_180) 2L else 0L),
  note = c(
    if (have_165) "Cell Ranger filtered MTX (KL1, KL2)" else "NO-GO: no processed 10x MTX",
    if (have_180) "author-QC 10x MTX (K, KL)" else "NO-GO: no processed 10x MTX"
  ),
  stringsAsFactors = FALSE
)
write.table(series_status, file.path(tab_dir, "series_status.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (!have_165 && !have_180) {
  writeLines(
    paste0(
      "# FINDING — GSE165641 + GSE180963 integration (Cldn4-only)\n\n",
      "Both series are **no-go**: no processed 10x matrix found. Stopped. Thesis unchanged.\n"
    ),
    file.path(out_dir, "FINDING.md")
  )
  stop("Both series no-go: no matrices.")
}

objs_ds <- list()

if (have_165) {
  message("CreateSeuratObject: GSE165641 KL1 + KL2")
  m1 <- CreateSeuratObject(read_10x_dir(p165_kl1), project = "GSE165641", min.cells = 0, min.features = 0)
  m1$mouse <- "GSE165641_KL1"
  m1$gsm <- "GSM5047302"
  m1$genotype <- "KL"
  m1$arm <- "KL"
  m1$strain <- "C57BL/6"
  m1$dataset <- "GSE165641"
  m1$orig.ident <- "GSE165641_KL1"
  m2 <- CreateSeuratObject(read_10x_dir(p165_kl2), project = "GSE165641", min.cells = 0, min.features = 0)
  m2$mouse <- "GSE165641_KL2"
  m2$gsm <- "GSM5047303"
  m2$genotype <- "KL"
  m2$arm <- "KL"
  m2$strain <- "C57BL/6"
  m2$dataset <- "GSE165641"
  m2$orig.ident <- "GSE165641_KL2"
  ds165 <- merge(m1, y = m2, add.cell.ids = c("KL1", "KL2"), project = "GSE165641")
  if (packageVersion("SeuratObject") >= "5.0.0") ds165 <- JoinLayers(ds165)
  ds165[["percent.mt"]] <- PercentageFeatureSet(ds165, pattern = "^mt-")
  n_before_165 <- ncol(ds165)
  ds165 <- subset(ds165, subset = nFeature_RNA >= 200 & nCount_RNA >= 500 & percent.mt < 25)
  ds165$n_pre_qc <- n_before_165
  objs_ds[["GSE165641"]] <- ds165
  rm(m1, m2, ds165); gc()
}

if (have_180) {
  message("CreateSeuratObject: GSE180963 K + KL")
  m1 <- CreateSeuratObject(read_10x_dir(p180_k), project = "GSE180963", min.cells = 0, min.features = 0)
  m1$mouse <- "GSE180963_K"
  m1$gsm <- "GSM5481386"
  m1$genotype <- "K"
  m1$arm <- "K"
  m1$strain <- "FVB"
  m1$dataset <- "GSE180963"
  m1$orig.ident <- "GSE180963_K"
  m2 <- CreateSeuratObject(read_10x_dir(p180_kl), project = "GSE180963", min.cells = 0, min.features = 0)
  m2$mouse <- "GSE180963_KL"
  m2$gsm <- "GSM5481387"
  m2$genotype <- "KL"
  m2$arm <- "KL"
  m2$strain <- "FVB"
  m2$dataset <- "GSE180963"
  m2$orig.ident <- "GSE180963_KL"
  ds180 <- merge(m1, y = m2, add.cell.ids = c("K", "KL"), project = "GSE180963")
  if (packageVersion("SeuratObject") >= "5.0.0") ds180 <- JoinLayers(ds180)
  ds180[["percent.mt"]] <- PercentageFeatureSet(ds180, pattern = "^mt-")
  n_before_180 <- ncol(ds180)
  ds180 <- subset(ds180, subset = nFeature_RNA >= 200 & nCount_RNA >= 500 & percent.mt < 25)
  ds180$n_pre_qc <- n_before_180
  objs_ds[["GSE180963"]] <- ds180
  rm(m1, m2, ds180); gc()
}

n_ds <- length(objs_ds)
integration_method <- NA_character_
harmony_used <- FALSE

if (n_ds == 1) {
  message("Only one series has a matrix — no cross-series integration; keep that object.")
  obj <- objs_ds[[1]]
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  integration_method <- "single_series_no_merge"
} else {
  message("Merge datasets then Harmony IntegrateLayers")
  obj <- merge(
    objs_ds[[1]],
    y = objs_ds[[2]],
    add.cell.ids = names(objs_ds),
    project = "GSE165641_GSE180963"
  )
  rm(objs_ds); gc()
  # Split RNA layers by dataset so IntegrateLayers sees two batches.
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$dataset)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  integ_ok <- FALSE
  tryCatch({
    obj <- IntegrateLayers(
      object = obj,
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "harmony",
      verbose = TRUE
    )
    integ_ok <- TRUE
    harmony_used <- TRUE
    integration_method <- "Seurat_IntegrateLayers_HarmonyIntegration"
  }, error = function(e) {
    message("HarmonyIntegration failed: ", conditionMessage(e), " — falling back to RPCA")
  })
  if (!integ_ok) {
    obj <- IntegrateLayers(
      object = obj,
      method = RPCAIntegration,
      orig.reduction = "pca",
      new.reduction = "integrated.rpca",
      verbose = TRUE
    )
    integration_method <- "Seurat_IntegrateLayers_RPCAIntegration"
    red <- "integrated.rpca"
  } else {
    red <- "harmony"
  }
  obj <- FindNeighbors(obj, reduction = red, dims = 1:20, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = red, dims = 1:20, verbose = FALSE)
  obj <- JoinLayers(obj)
}

# ---- 2. Marker compartments on the integrated object ----
universe <- rownames(obj)
inv_rows <- lapply(unique(c(CLDN4, T_NK_SCORE, IFN, MHC, TJ, EPI_CORE, HOST_LUNG, AUDIT)), function(g) {
  data.frame(gene = g, present = g %in% universe, stringsAsFactors = FALSE)
})
inv <- do.call(rbind, inv_rows)
inv$set <- ifelse(inv$gene == CLDN4, "Cldn4",
           ifelse(inv$gene %in% T_NK_SCORE, "T_NK",
           ifelse(inv$gene %in% IFN, "IFN",
           ifelse(inv$gene %in% MHC, "MHC",
           ifelse(inv$gene %in% TJ, "TJ", "lineage")))))
write.table(inv, file.path(tab_dir, "gene_inventory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (!(CLDN4 %in% universe)) {
  stop("Cldn4 row missing after merge — cannot score Cldn4-only.")
}

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
epi_loose <- (epcam > 0) | (cdh1 > 0 & krt8 > 0) | (krt18 > 0 & krt19 > 0)
epi_tight <- (epcam > 0) & (struct > 0) & (ptprc == 0)
tnk_mark <- pos_any(counts, T_NK_CALL)
is_epi <- epi_tight
is_tnk <- tnk_mark & !is_epi

obj$is_epi <- is_epi
obj$is_epi_loose <- epi_loose
obj$is_tnk <- is_tnk
obj$sftpc_pos <- sftpc > 0
obj$ptprc_pos <- ptprc > 0
obj$cell_class <- ifelse(is_epi, "epithelial",
                  ifelse(is_tnk, "T/NK", "other"))

obj$Cldn4_count <- as.numeric(counts[CLDN4, ])
obj$Cldn4_logn <- as.numeric(logn[CLDN4, ])
obj$Cldn4_pos <- obj$Cldn4_count > 0
obj$ifn_score <- mean_mod(logn, IFN)
obj$mhc_score <- mean_mod(logn, MHC)
obj$tj_score <- mean_mod(logn, TJ)
obj$tnk_score <- mean_mod(logn, T_NK_SCORE)

ifn_use <- present_in(IFN, universe)
mhc_use <- present_in(MHC, universe)
tj_use <- present_in(TJ, universe)
if (length(ifn_use) >= 2 && length(mhc_use) >= 2 && length(tj_use) >= 2) {
  obj <- AddModuleScore(
    obj,
    features = list(IFN = ifn_use, MHC = mhc_use, TJ = tj_use),
    name = c("IFN_seurat", "MHC_seurat", "TJ_seurat"),
    ctrl = min(50, max(10, length(universe) %/% 50)),
    seed = 1
  )
}

# ---- 3. Mouse-level table (the unit) ----
md <- slot(obj, "meta.data")
mice <- sort(unique(md$mouse))
n_mice <- length(mice)

mouse_rows <- lapply(mice, function(lab) {
  w <- md$mouse == lab
  we <- w & md$is_epi
  wl <- w & md$is_epi_loose
  wn <- w & md$is_tnk
  data.frame(
    mouse = lab,
    dataset = unique(md$dataset[w]),
    gsm = unique(md$gsm[w]),
    genotype = unique(md$genotype[w]),
    arm = unique(md$arm[w]),
    strain = unique(md$strain[w]),
    n_cells = sum(w),
    n_epi = sum(we),
    n_epi_loose = sum(wl),
    n_tnk = sum(wn),
    frac_tnk = mean(md$is_tnk[w]),
    frac_epi = mean(md$is_epi[w]),
    frac_sftpc = mean(md$sftpc_pos[w]),
    Cldn4_all_mean = mean(md$Cldn4_logn[w]),
    Cldn4_all_pctpos = 100 * mean(md$Cldn4_pos[w]),
    Cldn4_epi_mean = if (any(we)) mean(md$Cldn4_logn[we]) else NA_real_,
    Cldn4_epi_pctpos = if (any(we)) 100 * mean(md$Cldn4_pos[we]) else NA_real_,
    IFN_epi_mean = if (any(we)) mean(md$ifn_score[we]) else NA_real_,
    MHC_epi_mean = if (any(we)) mean(md$mhc_score[we]) else NA_real_,
    TJ_epi_mean = if (any(we)) mean(md$tj_score[we]) else NA_real_,
    stringsAsFactors = FALSE
  )
})
per_mouse <- do.call(rbind, mouse_rows)
write.table(per_mouse, file.path(tab_dir, "per_mouse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(per_mouse, file.path(tab_dir, "mouse_level_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

comp <- as.data.frame(table(mouse = md$mouse, cell_class = md$cell_class),
                      stringsAsFactors = FALSE)
write.table(comp, file.path(tab_dir, "compartment_by_mouse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cldn4_pos <- md[md$Cldn4_pos, c(
  "mouse", "dataset", "gsm", "arm", "cell_class", "is_epi", "is_epi_loose",
  "is_tnk", "sftpc_pos", "ptprc_pos", "Cldn4_count", "Cldn4_logn",
  "ifn_score", "mhc_score", "tj_score"
)]
write.table(cldn4_pos, file.path(tab_dir, "cldn4_positive_cells.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

n_cells <- nrow(md)
n_epi <- sum(md$is_epi)
n_tnk <- sum(md$is_tnk)
n_cldn4 <- sum(md$Cldn4_pos)
n_k <- sum(per_mouse$arm == "K")
n_kl <- sum(per_mouse$arm == "KL")

sp_tnk <- spearman_safe(per_mouse$Cldn4_epi_mean, per_mouse$frac_tnk)
sp_tnk_pct <- spearman_safe(per_mouse$Cldn4_epi_pctpos, per_mouse$frac_tnk)
sp_ifn <- spearman_safe(per_mouse$Cldn4_epi_mean, per_mouse$IFN_epi_mean)
sp_mhc <- spearman_safe(per_mouse$Cldn4_epi_mean, per_mouse$MHC_epi_mean)
sp_tj <- spearman_safe(per_mouse$Cldn4_epi_mean, per_mouse$TJ_epi_mean)

q4q1_allowed <- n_mice >= 6
q4q1 <- data.frame(
  contrast = c("epithelial IFN Q4 vs Q1", "epithelial MHC Q4 vs Q1", "epithelial TJ Q4 vs Q1"),
  n_mice = n_mice,
  allowed = q4q1_allowed,
  note = if (q4q1_allowed) "would split on epithelial Cldn4" else "NO-GO: n_mice<6; Q4 vs Q1 not run",
  stringsAsFactors = FALSE
)
if (q4q1_allowed) {
  q <- as.integer(cut(per_mouse$Cldn4_epi_mean,
                      breaks = quantile(per_mouse$Cldn4_epi_mean, probs = c(0, 0.25, 0.75, 1), na.rm = TRUE),
                      include.lowest = TRUE, labels = FALSE))
  # 1 = Q1, 3 = Q4 after 3-bin cut on 0/25/75/1 is messy; use rank quartiles
  rq <- dplyr::ntile(per_mouse$Cldn4_epi_mean, 4)
  q1 <- rq == 1
  q4 <- rq == 4
  wil <- function(a, b) {
    if (sum(is.finite(a)) < 2 || sum(is.finite(b)) < 2) return(NA_real_)
    suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value)
  }
  q4q1 <- rbind(
    data.frame(contrast = "epithelial IFN Q4 vs Q1", n_mice = n_mice, allowed = TRUE,
               n_Q4 = sum(q4), n_Q1 = sum(q1),
               mean_Q4 = mean(per_mouse$IFN_epi_mean[q4], na.rm = TRUE),
               mean_Q1 = mean(per_mouse$IFN_epi_mean[q1], na.rm = TRUE),
               MW_p = wil(per_mouse$IFN_epi_mean[q4], per_mouse$IFN_epi_mean[q1]),
               note = "mouse-level; split on epithelial Cldn4 mean",
               stringsAsFactors = FALSE),
    data.frame(contrast = "epithelial MHC Q4 vs Q1", n_mice = n_mice, allowed = TRUE,
               n_Q4 = sum(q4), n_Q1 = sum(q1),
               mean_Q4 = mean(per_mouse$MHC_epi_mean[q4], na.rm = TRUE),
               mean_Q1 = mean(per_mouse$MHC_epi_mean[q1], na.rm = TRUE),
               MW_p = wil(per_mouse$MHC_epi_mean[q4], per_mouse$MHC_epi_mean[q1]),
               note = "mouse-level; split on epithelial Cldn4 mean",
               stringsAsFactors = FALSE),
    data.frame(contrast = "epithelial TJ Q4 vs Q1", n_mice = n_mice, allowed = TRUE,
               n_Q4 = sum(q4), n_Q1 = sum(q1),
               mean_Q4 = mean(per_mouse$TJ_epi_mean[q4], na.rm = TRUE),
               mean_Q1 = mean(per_mouse$TJ_epi_mean[q1], na.rm = TRUE),
               MW_p = wil(per_mouse$TJ_epi_mean[q4], per_mouse$TJ_epi_mean[q1]),
               note = "mouse-level; split on epithelial Cldn4 mean",
               stringsAsFactors = FALSE)
  )
}
write.table(q4q1, file.path(tab_dir, "q4q1_ifn_mhc_tj.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

assoc <- data.frame(
  contrast = c(
    "Cldn4 epi mean vs T/NK fraction",
    "Cldn4 epi %pos vs T/NK fraction",
    "Cldn4 epi mean vs IFN epi",
    "Cldn4 epi mean vs MHC epi",
    "Cldn4 epi mean vs TJ epi"
  ),
  n_mice = c(sp_tnk$n, sp_tnk_pct$n, sp_ifn$n, sp_mhc$n, sp_tj$n),
  spearman_rho = c(sp_tnk$rho, sp_tnk_pct$rho, sp_ifn$rho, sp_mhc$rho, sp_tj$rho),
  spearman_p = c(sp_tnk$p, sp_tnk_pct$p, sp_ifn$p, sp_mhc$p, sp_tj$p),
  note = c(sp_tnk$note, sp_tnk_pct$note, sp_ifn$note, sp_mhc$note, sp_tj$note),
  stringsAsFactors = FALSE
)
write.table(assoc, file.path(tab_dir, "mouse_level_associations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

honest <- data.frame(
  item = c(
    "GEO series with matrix", "mice (unit, after merge)", "K mice", "KL mice",
    "datasets integrated", "cells (not the unit)", "genes after merge",
    "marker epithelial cells (tight)", "T/NK cells", "Cldn4-positive cells",
    "dual-high TACSTD2 x Cldn4", "private 8-KL mice used", "ICI / PD-1 arms",
    "Q4 vs Q1 IFN/MHC/TJ", "Spearman Cldn4 vs T/NK"
  ),
  n = c(
    n_ds, n_mice, n_k, n_kl,
    n_ds, n_cells, nrow(obj),
    n_epi, n_tnk, n_cldn4,
    0, 0, 0,
    if (q4q1_allowed) n_mice else 0,
    if (is.finite(sp_tnk$rho)) n_mice else 0
  ),
  note = c(
    paste(series_status$series[series_status$matrix_found], collapse = "+"),
    "one row per GEO library / mouse",
    paste(per_mouse$mouse[per_mouse$arm == "K"], collapse = ","),
    paste(per_mouse$mouse[per_mouse$arm == "KL"], collapse = ","),
    integration_method,
    "QC: nFeature>=200, nCount>=500, percent.mt<25",
    "symbol union; Cldn4 required",
    "Epcam+ AND structural+ AND Ptprc-; Cldn4 not a caller",
    "Cd3d/e/g or Cd8a or Nkg7/Ncr1/Klrb1c; tight epi wins",
    "count > 0",
    "not defined",
    "public GEO only",
    "untreated GEMM",
    if (q4q1_allowed) "run" else "NO-GO n_mice<6",
    if (is.finite(sp_tnk$rho)) sprintf("rho=%.3f p=%.3g", sp_tnk$rho, sp_tnk$p) else "n<4"
  ),
  stringsAsFactors = FALSE
)
write.table(honest, file.path(tab_dir, "honest_n.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 4. Figures ----
theme_set(theme_bw(base_size = 11))
pm <- per_mouse
pm$mouse <- factor(pm$mouse, levels = pm$mouse)

p_n <- ggplot(
  honest[honest$item %in% c("mice (unit, after merge)", "K mice", "KL mice", "cells (not the unit)"), ],
  aes(x = item, y = as.numeric(n))
) +
  geom_col(fill = "#4C72B0") +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  labs(title = "Integrated honest n = mice",
       subtitle = paste0("GSE165641 + GSE180963. Do not write n = ", n_cells, " cells."),
       x = NULL, y = "n") +
  ylim(0, max(as.numeric(honest$n[honest$item == "cells (not the unit)"]), 10) * 1.15)
ggsave(file.path(fig_dir, "fig_honest_n.png"), p_n, width = 7.4, height = 4.2, dpi = 140)
ggsave(file.path(fig_dir, "fig_honest_n.pdf"), p_n, width = 7.4, height = 4.2)

p_xy <- ggplot(pm, aes(x = Cldn4_epi_mean, y = frac_tnk, color = dataset, shape = arm, label = mouse)) +
  geom_point(size = 3.6) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 20) +
  scale_color_manual(values = c(GSE165641 = "#4C72B0", GSE180963 = "#C44E52")) +
  labs(
    title = "Integrated mouse-level Cldn4 (epithelium) vs T/NK fraction",
    subtitle = sprintf(
      "n = %d mice. Spearman %s. Q4 vs Q1 IFN/MHC/TJ %s.",
      n_mice,
      if (is.finite(sp_tnk$rho)) sprintf("rho=%.2f p=%.3g", sp_tnk$rho, sp_tnk$p) else "not defined",
      if (q4q1_allowed) "run" else "no-go (n<6)"
    ),
    x = "Cldn4 mean lognorm in marker epithelium",
    y = "T/NK fraction"
  )
ggsave(file.path(fig_dir, "fig_cldn4_vs_tnk.png"), p_xy, width = 6.8, height = 5.0, dpi = 140)
ggsave(file.path(fig_dir, "fig_cldn4_vs_tnk.pdf"), p_xy, width = 6.8, height = 5.0)

p_c4 <- ggplot(pm, aes(x = mouse, y = Cldn4_epi_pctpos, fill = dataset)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%\nn_epi=%d", Cldn4_epi_pctpos, n_epi)),
            vjust = -0.1, size = 2.8) +
  scale_fill_manual(values = c(GSE165641 = "#4C72B0", GSE180963 = "#C44E52")) +
  labs(title = "Cldn4 % positive in marker epithelium",
       subtitle = "Integrated object. Cldn4-only. Not a caller.",
       y = "Cldn4+ % of epithelium", x = NULL) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  ylim(0, max(pm$Cldn4_epi_pctpos, na.rm = TRUE) * 1.4)
ggsave(file.path(fig_dir, "fig_cldn4_epi.png"), p_c4, width = 7.2, height = 4.6, dpi = 140)
ggsave(file.path(fig_dir, "fig_cldn4_epi.pdf"), p_c4, width = 7.2, height = 4.6)

mod_long <- rbind(
  data.frame(mouse = pm$mouse, dataset = pm$dataset, module = "IFN", score = pm$IFN_epi_mean),
  data.frame(mouse = pm$mouse, dataset = pm$dataset, module = "MHC", score = pm$MHC_epi_mean),
  data.frame(mouse = pm$mouse, dataset = pm$dataset, module = "TJ", score = pm$TJ_epi_mean)
)
p_mod <- ggplot(mod_long, aes(x = mouse, y = score, fill = dataset)) +
  geom_col(width = 0.65) +
  facet_wrap(~module, scales = "free_y") +
  scale_fill_manual(values = c(GSE165641 = "#4C72B0", GSE180963 = "#C44E52")) +
  labs(title = "Epithelial IFN / MHC / TJ (mean lognorm)",
       subtitle = if (q4q1_allowed) "Q4 vs Q1 tabulated." else "n_mice<6: Q4 vs Q1 not run.",
       y = "mean lognorm", x = NULL) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 8))
ggsave(file.path(fig_dir, "fig_epi_ifn_mhc_tj.png"), p_mod, width = 9.0, height = 4.6, dpi = 140)
ggsave(file.path(fig_dir, "fig_epi_ifn_mhc_tj.pdf"), p_mod, width = 9.0, height = 4.6)

p_umap_ds <- DimPlot(obj, group.by = "dataset", reduction = "umap") +
  ggtitle(sprintf("UMAP by dataset (%s)", integration_method))
ggsave(file.path(fig_dir, "fig_umap_dataset.png"), p_umap_ds, width = 6.6, height = 5.2, dpi = 140)
ggsave(file.path(fig_dir, "fig_umap_dataset.pdf"), p_umap_ds, width = 6.6, height = 5.2)

p_umap_cc <- DimPlot(obj, group.by = "cell_class", reduction = "umap") +
  ggtitle("UMAP by cell class (marker epithelial / T/NK / other)")
ggsave(file.path(fig_dir, "fig_umap_cellclass.png"), p_umap_cc, width = 6.6, height = 5.2, dpi = 140)
ggsave(file.path(fig_dir, "fig_umap_cellclass.pdf"), p_umap_cc, width = 6.6, height = 5.2)

p_umap_m <- DimPlot(obj, group.by = "mouse", reduction = "umap") +
  ggtitle("UMAP by mouse (unit)")
ggsave(file.path(fig_dir, "fig_umap_mouse.png"), p_umap_m, width = 6.8, height = 5.2, dpi = 140)

p_umap_4 <- FeaturePlot(obj, features = "Cldn4", reduction = "umap", order = TRUE) +
  ggtitle("UMAP Cldn4 (lognorm)")
ggsave(file.path(fig_dir, "fig_umap_cldn4.png"), p_umap_4, width = 6.6, height = 5.2, dpi = 140)

# Slim integrated object (counts + data + reductions). Large RDS is gitignored.
keep_red <- intersect(c("pca", "umap", "harmony", "integrated.rpca"), names(obj@reductions))
slim <- tryCatch({
  DietSeurat(obj, assays = "RNA", dimreducs = keep_red, layers = c("counts", "data"))
}, error = function(e) obj)
saveRDS(slim, file.path(obj_dir, "integrated_gse165641_gse180963.rds"))
writeLines(
  paste0(
    "Integrated Seurat object written by analyze.R.\n",
    "cells=", n_cells, " mice=", n_mice, " method=", integration_method, "\n",
    "Seurat ", as.character(packageVersion("Seurat")),
    " harmony ", as.character(packageVersion("harmony")), "\n"
  ),
  file.path(obj_dir, "OBJECT.md")
)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(
    list(
      series = c("GSE165641", "GSE180963")[c(have_165, have_180)],
      engine = "R_Seurat",
      seurat_version = as.character(packageVersion("Seurat")),
      harmony_version = as.character(packageVersion("harmony")),
      create_seurat_object_per_dataset = TRUE,
      integration_method = integration_method,
      harmony_used = harmony_used,
      python_primary = FALSE,
      private_8_kl = FALSE,
      cldn4_only = TRUE,
      dual_high = FALSE,
      unit = "mouse",
      n_mice = n_mice,
      n_k = n_k,
      n_kl = n_kl,
      n_cells = n_cells,
      n_genes = nrow(obj),
      n_epi = n_epi,
      n_tnk = n_tnk,
      n_cldn4_pos = n_cldn4,
      q4q1_allowed = q4q1_allowed,
      associations = assoc,
      per_mouse = per_mouse
    ),
    file.path(tab_dir, "summary.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = 8
  )
}

# ---- 5. FINDING.md (live numbers) ----
mouse_md_rows <- apply(per_mouse, 1, function(r) {
  paste0(
    "| ", r[["mouse"]], " | ", r[["dataset"]], " | ", r[["gsm"]], " | ",
    r[["genotype"]], " | ", r[["strain"]], " | ",
    r[["n_cells"]], " | ", r[["n_epi"]], " | ", r[["n_tnk"]], " | ",
    fmt(as.numeric(r[["frac_tnk"]])), " | ",
    fmt(as.numeric(r[["Cldn4_all_mean"]]), 4), " | ",
    fmt(as.numeric(r[["Cldn4_all_pctpos"]]), 2), " | ",
    fmt(as.numeric(r[["Cldn4_epi_mean"]]), 4), " | ",
    fmt(as.numeric(r[["Cldn4_epi_pctpos"]]), 2), " | ",
    fmt(as.numeric(r[["IFN_epi_mean"]]), 3), " | ",
    fmt(as.numeric(r[["MHC_epi_mean"]]), 3), " | ",
    fmt(as.numeric(r[["TJ_epi_mean"]]), 3), " |"
  )
})

rho_txt <- function(sp) {
  if (is.na(sp$rho)) return(sprintf("n=%d, ρ undefined", sp$n))
  sprintf("n=%d mice, ρ=%.3f, p=%.3g", sp$n, sp$rho, sp$p)
}

lineage_present <- present_in(c(EPI_CORE, HOST_LUNG, T_NK_CALL, CLDN4, AUDIT, TJ), universe)
lineage_missing <- setdiff(c(EPI_CORE, HOST_LUNG, T_NK_CALL, CLDN4, AUDIT, TJ), universe)

finding <- paste0(
  "# FINDING — integrate GSE165641 + GSE180963 (Cldn4-only, Seurat + Harmony)\n\n",
  "**ADDITIVE. Public mouse. Cldn4-only. No dual-high.** Thesis is already correct and is not rewritten. ",
  "This folder **integrates** two public KL/K GEMM 10x series — ",
  "[GSE165641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165641) ",
  "(Wang / Zhong, *Adv Sci* 2021, [PMID 34369094](https://pubmed.ncbi.nlm.nih.gov/34369094/)) and ",
  "[GSE180963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963) ",
  "(Bai / Guo / Zhang / Long / Dong, [CAN-22-1740](https://doi.org/10.1158/0008-5472.can-22-1740)) — ",
  "into **one** Seurat object. The point is mouse-level n **after merge**, not two separate catalogs. ",
  "Primary engine is **R + Seurat ", as.character(packageVersion("Seurat")),
  "** (`CreateSeuratObject` per dataset → `IntegrateLayers` / ", integration_method, "). ",
  "No Python-only primary. Private 8-KL matrices were not opened. ",
  "FACS-only epithelium (GSE179502/GSE154989), CD45-only (GSE127465), subcutaneous Cldn4-floor (GSE267321), ",
  "and injury KO (GSE50927) were **not** merged.\n\n",
  "**Verdict.** Both series have processed 10x MTX. An integrated object **was built**. ",
  "Honest n after merge = **", n_mice, " mice** (", n_kl, " KL + ", n_k, " K). ",
  "Cldn4 vs T/NK is scored at that unit. ",
  "Epithelial IFN/MHC/TJ Q4 vs Q1 is ",
  if (q4q1_allowed) "**run**" else "**no-go (n_mice<6)**",
  ". Do not write n = ", n_cells, " cells. Thesis unchanged.\n\n",
  "---\n\n",
  "## Decision\n\n",
  "| Question | Answer |\n",
  "|---|---|\n",
  "| GSE165641 matrix | ", if (have_165) "**yes** — Cell Ranger filtered MTX, 2 KL mice" else "**no-go — dropped**", " |\n",
  "| GSE180963 matrix | ", if (have_180) "**yes** — author-QC MTX, 1 K + 1 KL" else "**no-go — dropped**", " |\n",
  "| CreateSeuratObject per dataset | **yes** — Seurat ", as.character(packageVersion("Seurat")), " |\n",
  "| Integration | **", integration_method, "** (Harmony package ", as.character(packageVersion("harmony")), ") |\n",
  "| Mice after merge (unit) | **", n_mice, "** |\n",
  "| Marker epithelium (tight) | **", n_epi, "** cells |\n",
  "| Cldn4 row present | **yes** — **", n_cldn4, "** cells > 0 |\n",
  "| Cldn4 vs T/NK at mouse unit | **scored** — ", rho_txt(sp_tnk), " |\n",
  "| Epithelial IFN/MHC/TJ Q4 vs Q1 | ",
  if (q4q1_allowed) "**run**" else "**no-go — n_mice<6**", " |\n",
  "| Unit | **mouse** |\n",
  "| Dual-high TACSTD2 × Cldn4 | **not defined** |\n",
  "| Private 8 KL | **not used** |\n",
  "| ICI / PD-1 | **no** |\n\n",
  "---\n\n",
  "## Honest n\n\n",
  "| item | n | note |\n",
  "|---|---:|---|\n",
  "| GEO series with matrix | **", n_ds, "** | ", paste(series_status$series[series_status$matrix_found], collapse = " + "), " |\n",
  "| Mice (unit, after merge) | **", n_mice, "** | not cells |\n",
  "| K mice | **", n_k, "** | GSE180963 only |\n",
  "| KL mice | **", n_kl, "** | GSE165641 KL1/KL2 + GSE180963 KL |\n",
  "| Cells after QC | **", n_cells, "** | not the unit |\n",
  "| Genes after merge | **", nrow(obj), "** | symbol union |\n",
  "| Marker epithelial cells | **", n_epi, "** | Epcam+ structural+ Ptprc- |\n",
  "| T/NK cells | **", n_tnk, "** | marker call |\n",
  "| Cldn4-positive cells | **", n_cldn4, "** | count > 0 |\n",
  "| Dual-high | **0** | not defined |\n",
  "| Private 8-KL mice | **0** | public GEO only |\n",
  "| ICI arms | **0** | untreated GEMM |\n\n",
  "Do not write n = ", n_cells, ". Do not treat two GEO series as two mice.\n\n",
  "---\n\n",
  "## Per-mouse table (the actual unit)\n\n",
  "| mouse | dataset | GSM | genotype | strain | n cells | n epi | n T/NK | frac T/NK | Cldn4 all | Cldn4 %pos | Cldn4 epi | Cldn4 epi %pos | IFN epi | MHC epi | TJ epi |\n",
  "|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
  paste(mouse_md_rows, collapse = "\n"), "\n\n",
  "---\n\n",
  "## Cldn4 vs T/NK (mouse unit)\n\n",
  rho_txt(sp_tnk), " for epithelial Cldn4 **mean** vs T/NK fraction. ",
  rho_txt(sp_tnk_pct), " for epithelial Cldn4 **%pos** vs T/NK fraction. ",
  "n = ", n_mice, " is the floor for Spearman and is **not powered**. ",
  "Strain / dataset is a confounder (C57BL/6 vs FVB). Harmony corrects the embedding, not the mouse-level scores. ",
  "Do not write a Cldn4–T/NK law from this object.\n\n",
  "Cldn4 vs epithelial IFN / MHC / TJ (descriptive Spearman, same n): IFN ", rho_txt(sp_ifn),
  "; MHC ", rho_txt(sp_mhc), "; TJ ", rho_txt(sp_tj), ".\n\n",
  "## Epithelial IFN / MHC / TJ Q4 vs Q1\n\n",
  if (q4q1_allowed) {
    paste0("n_mice≥6, contrast was run. See `tables/q4q1_ifn_mhc_tj.tsv`.\n")
  } else {
    paste0(
      "**No-go.** User lock: Q4 vs Q1 only if n_mice≥6. Here n_mice = **", n_mice,
      "**. Quartile split on 4 mice is 1 vs 1 and is not a test. ",
      "Mouse-level IFN / MHC / TJ means are in the table above; they are not a high-vs-low claim.\n"
    )
  },
  "\n",
  "---\n\n",
  "## Epithelium / cell class (honest)\n\n",
  "GEO deposits **no** author `Malignant` / `Epithelial` column on either series. Calling is markers on the **integrated** object:\n\n",
  "- **Tight epithelium (primary)** = Epcam+ **and** (Cdh1 or Krt8/18/19 or Cldn18)+ **and** Ptprc−.\n",
  "- **T/NK** = Cd3d / Cd3e / Cd3g / Cd8a / Nkg7 / Ncr1 / Klrb1c. Tight epithelium wins if both fire.\n",
  "- **Cldn4 is not a caller.** Tacstd2 is inventory only. No dual-high gate.\n",
  "- GSE180963 has heavy Sftpc ambient in the sibling digest; Sftpc is an ambient flag, not host AT2.\n\n",
  "Lineage genes present: ", paste(lineage_present, collapse = ", "), ".\n",
  "Missing: ", if (length(lineage_missing) == 0) "none of the core set" else paste(lineage_missing, collapse = ", "), ".\n\n",
  "---\n\n",
  "## Methods (short)\n\n",
  "1. Public only. Downloaded per-GSM 10x MTX from NCBI GEO. No SRA / FASTQ. No private 8-KL object.\n",
  "2. **GSE165641:** Cell Ranger filtered MTX, 2 KL mice (C57BL/6), Ad-Cre 10 weeks. ",
  "GSM5047303 SOFT `Lkb2fl/fl` treated as a typo for Lkb1 (series title + PMID).\n",
  "3. **GSE180963:** author-QC MTX (Seurat 3.1.5: 500–6000 features, mito < 20%). 1 K + 1 KL (FVB). ",
  "Sibling digest: the two samples were mixed in one 10x library and demultiplexed by label.\n",
  "4. `CreateSeuratObject` **per dataset** (mice as metadata), then `merge`. ",
  "QC: nFeature≥200, nCount≥500, percent.mt<25.\n",
  "5. Integration: split RNA layers by `dataset`, NormalizeData / VST / ScaleData / PCA, then ",
  "`Seurat::IntegrateLayers(method = HarmonyIntegration)` (harmony ", as.character(packageVersion("harmony")),
  "). UMAP on the Harmony reduction. Scores stay on log-normalized counts (not Harmony-corrected expression).\n",
  "6. Cldn4-only. Tacstd2 is inventory. No dual-high. TJ module **excludes** Cldn4.\n",
  "7. Positive = raw count > 0. Module scores = mean lognorm of present genes.\n",
  "8. Unit = mouse after merge. Q4 vs Q1 only if n_mice≥6.\n",
  "9. Thesis is not rewritten.\n\n",
  "```bash\n",
  "bash methods/seurat_integrate_165641_180963_cldn4/scripts/install_r.sh\n",
  "bash methods/seurat_integrate_165641_180963_cldn4/scripts/download.sh /tmp/geo/work\n",
  "Rscript methods/seurat_integrate_165641_180963_cldn4/scripts/analyze.R --data /tmp/geo/work --out methods/seurat_integrate_165641_180963_cldn4\n",
  "```\n\n",
  "---\n\n",
  "## How to read this\n\n",
  "- **This is an integration job.** Separate GSE165641 and GSE180963 catalogs already exist; they are not repeated here as the result.\n",
  "- **Honest n is ", n_mice, " mice after merge.** Not ", n_cells, " cells. Not 2 series.\n",
  "- **Cldn4 vs T/NK is four (or fewer) points.** Report the table. Do not write a law.\n",
  "- **Q4 vs Q1 IFN/MHC/TJ did not fire** unless n_mice≥6.\n",
  "- **Strain / lab batch remains.** Harmony is for the UMAP, not a license to ignore C57BL/6 vs FVB.\n",
  "- **No dual-high. No private 8 KL. No FACS-only / CD45-only / subQ / injury-KO merge. Thesis unchanged.**\n\n",
  "## Files\n\n",
  "- `tables/per_mouse.tsv` / `tables/mouse_level_scores.tsv` — unit-level table\n",
  "- `tables/mouse_level_associations.tsv` — Cldn4 vs T/NK and IFN/MHC/TJ Spearman\n",
  "- `tables/q4q1_ifn_mhc_tj.tsv` — no-go or results\n",
  "- `tables/cldn4_positive_cells.tsv`, `tables/compartment_by_mouse.tsv`\n",
  "- `tables/gene_inventory.tsv`, `tables/honest_n.tsv`, `tables/series_status.tsv`, `tables/summary.json`\n",
  "- `figures/fig_umap_dataset.png`, `fig_umap_cellclass.png` — required DimPlots\n",
  "- `figures/fig_cldn4_vs_tnk.png`, `fig_cldn4_epi.png`, `fig_epi_ifn_mhc_tj.png`, `fig_honest_n.png`\n",
  "- `objects/integrated_gse165641_gse180963.rds` — integrated Seurat object (gitignored if large)\n",
  "- `scripts/analyze.R` — R + Seurat primary\n\n",
  "## 结论\n\n",
  "GSE165641 与 GSE180963 都有公开 10x MTX，已用 **R + Seurat `CreateSeuratObject`（每个 dataset）+ Harmony `IntegrateLayers`** 做成 **一个** 整合对象。",
  "诚实 n = **", n_mice, " 只鼠**（合并后），不是 ", n_cells, " 个细胞。",
  "Cldn4 对 T/NK 在鼠单位上已计分；上皮 IFN/MHC/TJ 的 Q4 vs Q1 ",
  if (q4q1_allowed) "已做" else "因 n<6 未做",
  "。未打开私有 8 只 KL。无 dual-high，不改 thesis。\n"
)

writeLines(finding, file.path(out_dir, "FINDING.md"))

readme <- paste0(
  "# seurat_integrate_165641_180963_cldn4\n\n",
  "Integrate public [GSE165641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165641) + ",
  "[GSE180963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963) KL/K GEMM 10x. ",
  "**R + Seurat** (`CreateSeuratObject` per dataset) then **Harmony** via `IntegrateLayers`. ",
  "Cldn4-only. Honest n = mice after merge. See `FINDING.md`.\n\n",
  "```bash\n",
  "bash methods/seurat_integrate_165641_180963_cldn4/scripts/install_r.sh\n",
  "bash methods/seurat_integrate_165641_180963_cldn4/scripts/download.sh /tmp/geo/work\n",
  "Rscript methods/seurat_integrate_165641_180963_cldn4/scripts/analyze.R --data /tmp/geo/work --out methods/seurat_integrate_165641_180963_cldn4\n",
  "```\n"
)
writeLines(readme, file.path(out_dir, "README.md"))

message("wrote FINDING.md and tables under ", out_dir)
message("n_mice=", n_mice, " n_cells=", n_cells, " n_epi=", n_epi, " n_tnk=", n_tnk,
        " n_cldn4=", n_cldn4, " method=", integration_method)
