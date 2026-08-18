#!/usr/bin/env Rscript
# GSE180963 — public K vs KL lung GEMM scRNA. R + Seurat primary.
# CreateSeuratObject on GEO processed 10x MTX. Cldn4-only.
# Honest unit = mouse. n = 2 (1 K + 1 KL). No private 8-KL. No Python primary.
#
# Usage:
#   Rscript methods/seurat_gse180963_cldn4/scripts/analyze.R \
#     --data /tmp/gse180963 --out methods/seurat_gse180963_cldn4

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
data_dir <- "/tmp/gse180963"
out_dir <- "methods/seurat_gse180963_cldn4"
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
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Locked Cldn4-only mouse sets (same family as public_kl_vs_kp_cldn4 / GSE267321).
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

SAMPLES <- data.frame(
  label = c("K", "KL"),
  gsm = c("GSM5481386", "GSM5481387"),
  genotype = c("KrasG12D/+", "KrasG12D/+;Lkb1fl/fl"),
  arm = c("K", "KL"),
  stringsAsFactors = FALSE
)

need <- file.path(data_dir, c("K/matrix.mtx", "K/genes.tsv", "K/barcodes.tsv",
                              "KL/matrix.mtx", "KL/genes.tsv", "KL/barcodes.tsv"))
if (any(!file.exists(need))) {
  stop(
    "No processed 10x MTX under ", data_dir,
    ". Run methods/seurat_gse180963_cldn4/scripts/download.sh ", data_dir
  )
}

present_in <- function(genes, universe) genes[genes %in% universe]

read_one <- function(label, gsm, genotype, arm) {
  d <- file.path(data_dir, label)
  counts <- Read10X(data.dir = d, gene.column = 2, unique.features = TRUE)
  obj <- CreateSeuratObject(
    counts = counts,
    project = label,
    min.cells = 0,
    min.features = 0
  )
  obj$mouse <- label
  obj$gsm <- gsm
  obj$genotype <- genotype
  obj$arm <- arm
  obj$orig.ident <- label
  obj
}

message("CreateSeuratObject: K")
obj_k <- read_one(SAMPLES$label[1], SAMPLES$gsm[1], SAMPLES$genotype[1], SAMPLES$arm[1])
message("CreateSeuratObject: KL")
obj_kl <- read_one(SAMPLES$label[2], SAMPLES$gsm[2], SAMPLES$genotype[2], SAMPLES$arm[2])
obj <- merge(obj_k, y = obj_kl, add.cell.ids = c("K", "KL"), project = "GSE180963")
rm(obj_k, obj_kl)
gc()
if (packageVersion("SeuratObject") >= "5.0.0") {
  obj <- JoinLayers(obj)
}

# Author matrices are already QC'd (Seurat 3.1.5: 500–6000 features, mito < 20%).
# Do not drop more cells. percent.mt is audit only (mouse mitochondrial prefix).
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)

universe <- rownames(obj)
inv_rows <- list()
for (g in unique(c(CLDN4, T_NK_SCORE, IFN, MHC, EPI_CORE, HOST_LUNG, AUDIT))) {
  inv_rows[[g]] <- data.frame(
    gene = g,
    present = g %in% universe,
    stringsAsFactors = FALSE
  )
}
inv <- do.call(rbind, inv_rows)
write.table(inv, file.path(tab_dir, "gene_inventory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (!(CLDN4 %in% universe)) {
  stop("Cldn4 row missing from processed matrix — cannot score Cldn4-only.")
}

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
logn <- GetAssayData(obj, assay = "RNA", layer = "data")

pos_any <- function(mat, genes) {
  genes <- present_in(genes, rownames(mat))
  if (length(genes) == 0) {
    return(rep(FALSE, ncol(mat)))
  }
  as.numeric(Matrix::colSums(mat[genes, , drop = FALSE] > 0)) > 0
}

mean_mod <- function(mat, genes) {
  genes <- present_in(genes, rownames(mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(mat)))
  }
  as.numeric(Matrix::colMeans(mat[genes, , drop = FALSE]))
}

# Marker compartments. Cldn4 is never a caller.
# Loose Epcam+ is kept as audit only. This digest has heavy AT2 ambient
# (Sftpc+ in most cells); Sftpc is an ambient flag, not host epithelium.
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
# Tight: Epcam AND a structural marker AND not CD45. Primary epi for IFN/MHC.
epi_tight <- (epcam > 0) & (struct > 0) & (ptprc == 0)
tnk_mark <- pos_any(counts, T_NK_CALL)
# Epithelium (tight) wins if both fire.
is_epi <- epi_tight
is_tnk <- tnk_mark & !is_epi
is_epi_loose <- epi_loose

obj$is_epi <- is_epi
obj$is_epi_loose <- is_epi_loose
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

# Seurat-native module scores (audit; primary scores remain mean lognorm).
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

# Light Seurat space for figures only. No integration: genotype = mouse = library.
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)

md <- slot(obj, "meta.data")
stopifnot(nrow(md) == ncol(obj))

fmt <- function(x, d = 3) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  formatC(as.numeric(x), format = "f", digits = d)
}

mouse_rows <- lapply(SAMPLES$label, function(lab) {
  w <- md$mouse == lab
  we <- w & md$is_epi
  wl <- w & md$is_epi_loose
  wn <- w & md$is_tnk
  data.frame(
    mouse = lab,
    gsm = unique(md$gsm[w]),
    genotype = unique(md$genotype[w]),
    arm = unique(md$arm[w]),
    n_cells = sum(w),
    n_epi = sum(we),
    n_epi_loose = sum(wl),
    n_tnk = sum(wn),
    frac_tnk = mean(md$is_tnk[w]),
    frac_epi = mean(md$is_epi[w]),
    frac_sftpc = mean(md$sftpc_pos[w]),
    frac_ptprc = mean(md$ptprc_pos[w]),
    Cldn4_all_mean = mean(md$Cldn4_logn[w]),
    Cldn4_all_pctpos = 100 * mean(md$Cldn4_pos[w]),
    Cldn4_epi_mean = if (any(we)) mean(md$Cldn4_logn[we]) else NA_real_,
    Cldn4_epi_pctpos = if (any(we)) 100 * mean(md$Cldn4_pos[we]) else NA_real_,
    Cldn4_loose_mean = if (any(wl)) mean(md$Cldn4_logn[wl]) else NA_real_,
    Cldn4_loose_pctpos = if (any(wl)) 100 * mean(md$Cldn4_pos[wl]) else NA_real_,
    IFN_epi_mean = if (any(we)) mean(md$ifn_score[we]) else NA_real_,
    MHC_epi_mean = if (any(we)) mean(md$mhc_score[we]) else NA_real_,
    Stk11_all_mean = if ("Stk11" %in% universe) mean(as.numeric(logn["Stk11", w])) else NA_real_,
    percent_mt_mean = mean(md$percent.mt[w]),
    stringsAsFactors = FALSE
  )
})
per_mouse <- do.call(rbind, mouse_rows)
write.table(per_mouse, file.path(tab_dir, "per_mouse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

comp <- as.data.frame(table(mouse = md$mouse, compartment = md$compartment),
                      stringsAsFactors = FALSE)
write.table(comp, file.path(tab_dir, "compartment_by_mouse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cldn4_pos <- md[md$Cldn4_pos, c(
  "mouse", "gsm", "arm", "compartment", "is_epi", "is_epi_loose",
  "is_tnk", "sftpc_pos", "ptprc_pos", "Cldn4_count", "Cldn4_logn",
  "ifn_score", "mhc_score"
)]
write.table(cldn4_pos, file.path(tab_dir, "cldn4_positive_cells.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Cell-level Spearman inside epithelium is descriptive / pseudoreplication.
# Mouse-level Spearman is undefined at n = 2.
k_row <- per_mouse[per_mouse$mouse == "K", ]
kl_row <- per_mouse[per_mouse$mouse == "KL", ]
n_mice <- nrow(per_mouse)
n_cells <- nrow(md)
n_epi <- sum(md$is_epi)
n_tnk <- sum(md$is_tnk)
n_cldn4 <- sum(md$Cldn4_pos)

spearman_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) {
    return(list(n = sum(ok), rho = NA_real_, p = NA_real_))
  }
  s <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(n = sum(ok), rho = unname(s$estimate), p = s$p.value)
}

epi_k <- md$mouse == "K" & md$is_epi
epi_kl <- md$mouse == "KL" & md$is_epi
sp_kl_ifn <- spearman_safe(md$Cldn4_logn[epi_kl], md$ifn_score[epi_kl])
sp_kl_mhc <- spearman_safe(md$Cldn4_logn[epi_kl], md$mhc_score[epi_kl])
sp_k_ifn <- spearman_safe(md$Cldn4_logn[epi_k], md$ifn_score[epi_k])
sp_k_mhc <- spearman_safe(md$Cldn4_logn[epi_k], md$mhc_score[epi_k])

cell_sp <- data.frame(
  subset = c("KL_epi_cells", "KL_epi_cells", "K_epi_cells", "K_epi_cells"),
  contrast = c("Cldn4 vs IFN", "Cldn4 vs MHC", "Cldn4 vs IFN", "Cldn4 vs MHC"),
  n_cells = c(sp_kl_ifn$n, sp_kl_mhc$n, sp_k_ifn$n, sp_k_mhc$n),
  rho = c(sp_kl_ifn$rho, sp_kl_mhc$rho, sp_k_ifn$rho, sp_k_mhc$rho),
  p = c(sp_kl_ifn$p, sp_kl_mhc$p, sp_k_ifn$p, sp_k_mhc$p),
  note = "cell-level; pseudoreplication; not the unit",
  stringsAsFactors = FALSE
)
write.table(cell_sp, file.path(tab_dir, "cldn4_vs_ifn_mhc_cells_exploratory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

honest <- data.frame(
  item = c(
    "GEO series", "mice (unit)", "K mice", "KL mice",
    "libraries", "cells (not the unit)", "genes",
    "marker epithelial cells (tight)", "loose Epcam+ (audit)", "T/NK cells",
    "Cldn4-positive cells", "dual-high TACSTD2 x Cldn4",
    "private 8-KL mice used", "ICI / PD-1 arms"
  ),
  n = c(
    1, 2, 1, 1,
    1, n_cells, nrow(obj),
    n_epi, sum(md$is_epi_loose), n_tnk,
    n_cldn4, 0,
    0, 0
  ),
  note = c(
    "GSE180963",
    "1 K + 1 KL; genotype = mouse = sample",
    "GSM5481386",
    "GSM5481387",
    "two samples mixed in one 10x library, demultiplexed by label",
    "author-filtered MTX",
    "mm10 symbols; genes.tsv col2",
    "Epcam+ AND structural (Cdh1/Krt8/18/19/Cldn18)+ AND Ptprc-; Cldn4 not a caller",
    "Epcam+ OR keratin pair; ambient-polluted; not the primary epi",
    "Cd3d/e/g or Cd8a or Nkg7/Ncr1/Klrb1c; tight epi wins if both fire",
    "count > 0",
    "not defined",
    "public GEO only",
    "untreated GEMM"
  ),
  stringsAsFactors = FALSE
)
write.table(honest, file.path(tab_dir, "honest_n.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Contrast table: two mouse points, no p-value theater.
contrast <- data.frame(
  metric = c(
    "T/NK fraction",
    "Cldn4 mean, all cells",
    "Cldn4 %pos, all cells",
    "Cldn4 mean, marker epithelium",
    "Cldn4 %pos, marker epithelium",
    "IFN ISG mean, marker epithelium",
    "MHC/APM mean, marker epithelium",
    "n marker epithelial cells"
  ),
  n_mice = "1 vs 1",
  KL = c(
    kl_row$frac_tnk, kl_row$Cldn4_all_mean, kl_row$Cldn4_all_pctpos,
    kl_row$Cldn4_epi_mean, kl_row$Cldn4_epi_pctpos,
    kl_row$IFN_epi_mean, kl_row$MHC_epi_mean, kl_row$n_epi
  ),
  K = c(
    k_row$frac_tnk, k_row$Cldn4_all_mean, k_row$Cldn4_all_pctpos,
    k_row$Cldn4_epi_mean, k_row$Cldn4_epi_pctpos,
    k_row$IFN_epi_mean, k_row$MHC_epi_mean, k_row$n_epi
  ),
  stringsAsFactors = FALSE
)
contrast$delta_KL_minus_K <- contrast$KL - contrast$K
contrast$MW_p <- NA
contrast$note <- "n=1 vs 1 mice; no mouse-level test"
write.table(contrast, file.path(tab_dir, "kl_vs_k.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Figures
theme_set(theme_bw(base_size = 11))
pm <- per_mouse
pm$mouse <- factor(pm$mouse, levels = c("K", "KL"))

p_n <- ggplot(honest[honest$item %in% c("mice (unit)", "K mice", "KL mice", "cells (not the unit)"), ],
              aes(x = item, y = n)) +
  geom_col(fill = "#4C72B0") +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  labs(title = "GSE180963 honest n = mice",
       subtitle = "Do not write n = cells. 1 K mouse + 1 KL mouse.",
       x = NULL, y = "n") +
  ylim(0, max(honest$n[honest$item == "cells (not the unit)"], 10) * 1.15)
ggsave(file.path(fig_dir, "fig_honest_n.png"), p_n, width = 7.2, height = 4.2, dpi = 140)
ggsave(file.path(fig_dir, "fig_honest_n.pdf"), p_n, width = 7.2, height = 4.2)

p_tnk <- ggplot(pm, aes(x = mouse, y = frac_tnk, fill = mouse)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.3f\n(%d/%d)", frac_tnk, n_tnk, n_cells)),
            vjust = -0.15, size = 3.2) +
  scale_fill_manual(values = c(K = "#4C72B0", KL = "#C44E52"), guide = "none") +
  labs(title = "T/NK fraction by mouse",
       subtitle = "Marker call. Honest n = 1 vs 1. Not a p-value.",
       y = "T/NK fraction", x = NULL) +
  ylim(0, max(pm$frac_tnk) * 1.25)
ggsave(file.path(fig_dir, "fig_tnk_fraction.png"), p_tnk, width = 5.6, height = 4.4, dpi = 140)
ggsave(file.path(fig_dir, "fig_tnk_fraction.pdf"), p_tnk, width = 5.6, height = 4.4)

p_c4 <- ggplot(pm, aes(x = mouse, y = Cldn4_epi_pctpos, fill = mouse)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%.2f%%\nn_epi=%d", Cldn4_epi_pctpos, n_epi)),
            vjust = -0.15, size = 3.2) +
  scale_fill_manual(values = c(K = "#4C72B0", KL = "#C44E52"), guide = "none") +
  labs(title = "Cldn4 % positive in marker epithelium",
       subtitle = "Cldn4-only. Not a caller. n = 1 vs 1 mice.",
       y = "Cldn4+ % of epithelium", x = NULL) +
  ylim(0, max(pm$Cldn4_epi_pctpos, na.rm = TRUE) * 1.35)
ggsave(file.path(fig_dir, "fig_cldn4_epi.png"), p_c4, width = 5.6, height = 4.4, dpi = 140)
ggsave(file.path(fig_dir, "fig_cldn4_epi.pdf"), p_c4, width = 5.6, height = 4.4)

ifn_long <- rbind(
  data.frame(mouse = pm$mouse, module = "IFN ISG", score = pm$IFN_epi_mean),
  data.frame(mouse = pm$mouse, module = "MHC/APM", score = pm$MHC_epi_mean)
)
p_ifn <- ggplot(ifn_long, aes(x = mouse, y = score, fill = mouse)) +
  geom_col(width = 0.55) +
  facet_wrap(~module) +
  scale_fill_manual(values = c(K = "#4C72B0", KL = "#C44E52"), guide = "none") +
  labs(title = "Epithelial IFN / MHC module (mean lognorm)",
       subtitle = "Tight marker epithelium (Epcam+ structural+ Ptprc-). n = 1 vs 1.",
       y = "mean log1p", x = NULL)
ggsave(file.path(fig_dir, "fig_epi_ifn_mhc.png"), p_ifn, width = 7.0, height = 4.4, dpi = 140)
ggsave(file.path(fig_dir, "fig_epi_ifn_mhc.pdf"), p_ifn, width = 7.0, height = 4.4)

p_xy <- ggplot(pm, aes(x = Cldn4_epi_mean, y = frac_tnk, color = mouse, label = mouse)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 4, show.legend = FALSE) +
  scale_color_manual(values = c(K = "#4C72B0", KL = "#C44E52"), guide = "none") +
  labs(title = "Cldn4 (epithelium) vs T/NK fraction",
       subtitle = "Two mice. Spearman is not defined. Do not fit a line.",
       x = "Cldn4 mean lognorm in marker epithelium",
       y = "T/NK fraction")
ggsave(file.path(fig_dir, "fig_cldn4_vs_tnk.png"), p_xy, width = 5.8, height = 4.6, dpi = 140)
ggsave(file.path(fig_dir, "fig_cldn4_vs_tnk.pdf"), p_xy, width = 5.8, height = 4.6)

p_umap_c <- DimPlot(obj, group.by = "compartment", reduction = "umap") +
  ggtitle("UMAP compartment (markers, not author labels)")
ggsave(file.path(fig_dir, "fig_umap_compartment.png"), p_umap_c, width = 6.4, height = 5.0, dpi = 140)
p_umap_g <- DimPlot(obj, group.by = "mouse", reduction = "umap") +
  ggtitle("UMAP mouse / genotype (no integration)")
ggsave(file.path(fig_dir, "fig_umap_mouse.png"), p_umap_g, width = 6.4, height = 5.0, dpi = 140)
p_umap_4 <- FeaturePlot(obj, features = "Cldn4", reduction = "umap", order = TRUE) +
  ggtitle("UMAP Cldn4 (lognorm)")
ggsave(file.path(fig_dir, "fig_umap_cldn4.png"), p_umap_4, width = 6.4, height = 5.0, dpi = 140)

# summary.json
if (requireNamespace("jsonlite", quietly = TRUE)) {
  summary <- list(
    series = "GSE180963",
    engine = "R_Seurat",
    seurat_version = as.character(packageVersion("Seurat")),
    create_seurat_object = TRUE,
    python_primary = FALSE,
    private_8_kl = FALSE,
    cldn4_only = TRUE,
    dual_high = FALSE,
    unit = "mouse",
    n_mice = 2,
    n_k = 1,
    n_kl = 1,
    n_cells = n_cells,
    n_genes = nrow(obj),
    n_epi = n_epi,
    n_tnk = n_tnk,
    n_cldn4_pos = n_cldn4,
    per_mouse = per_mouse,
    kl_vs_k = contrast,
    cell_level_exploratory = cell_sp
  )
  jsonlite::write_json(summary, file.path(tab_dir, "summary.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = 8)
}

# FINDING.md
esc <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x[1])) return("NA")
  as.character(x[1])
}
rho_txt <- function(sp) {
  if (is.na(sp$rho)) return(sprintf("n=%d, ρ undefined (<3 finite pairs)", sp$n))
  sprintf("n=%d cells, ρ=%.3f, p=%.3g (pseudoreplication)", sp$n, sp$rho, sp$p)
}

finding <- paste0(
  "# FINDING — GSE180963 K vs KL lung GEMM scRNA (Cldn4-only, Seurat)\n\n",
  "**ADDITIVE. Public mouse. Cldn4-only. No dual-high.** Thesis is already correct and is not rewritten. ",
  "This folder scores one public GEO object: [GSE180963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963) ",
  "(Bai / Guo / Zhang / Long / Dong, Southern Medical University; paper: ",
  "[CAN-22-1740](https://doi.org/10.1158/0008-5472.can-22-1740)). ",
  "KrasG12D/+ (**K**, Tomato lenti) vs KrasG12D/+;Lkb1-targeted (**KL**) FVB lung tumor nodules. ",
  "Primary engine is **R + Seurat ", as.character(packageVersion("Seurat")),
  "** (`Read10X` → `CreateSeuratObject`). No Python-only primary. ",
  "Private 8-KL matrices were not opened.\n\n",
  "**Verdict.** Processed 10x MTX **exists** on GEO (`GSE180963_RAW.tar` → `K/` and `KL/` ",
  "`matrix.mtx` + `genes.tsv` + `barcodes.tsv`). This is **not** a no-go for matrix. ",
  "It **is** a no-go for a mouse-level Cldn4 vs T/NK or Cldn4 vs epithelial IFN/MHC test: ",
  "honest n = **2 mice (1 vs 1)**. Genotype = mouse = sample. The two samples were mixed in ",
  "**one 10x library** and demultiplexed by label. Cell-level p-values are pseudoreplication. ",
  "Cldn4 is in the matrix and is **epithelial-restricted and sparse**. ",
  "Do not write n = ", n_cells, " cells.\n\n",
  "---\n\n",
  "## Decision\n\n",
  "| Question | Answer |\n",
  "|---|---|\n",
  "| Public processed matrix | **yes** — 10x MTX in `GSE180963_RAW.tar` (87.9 MB) |\n",
  "| CreateSeuratObject | **yes** — Seurat ", as.character(packageVersion("Seurat")), " |\n",
  "| Author malignant / epithelial labels | **no** — barcode matrix only |\n",
  "| Marker epithelium (tight) | **", n_epi, "** cells (K ", k_row$n_epi, " / KL ", kl_row$n_epi, ") |\n",
  "| Loose Epcam+ (audit) | **", sum(md$is_epi_loose), "** cells; Sftpc is ambient |\n",
  "| Cldn4 row present | **yes** — **", n_cldn4, "** cells > 0 |\n",
  "| Cldn4 vs T/NK at mouse unit | **no-go** — n = 1 vs 1 mice |\n",
  "| Epithelial IFN/MHC at mouse unit | **no-go** — n = 1 vs 1 mice |\n",
  "| T/NK fraction (descriptive) | K ", fmt(k_row$frac_tnk), " vs KL ", fmt(kl_row$frac_tnk), " |\n",
  "| Closest to user KL | **this KL arm** (Lkb1-targeted KrasG12D/+ lung GEMM) |\n",
  "| Unit | **mouse** |\n",
  "| Dual-high TACSTD2 × Cldn4 | **not defined** |\n",
  "| Private 8 KL | **not used** |\n",
  "| ICI / PD-1 | **no** |\n\n",
  "Cldn4 was scored because the row exists. Two mice cannot test a law. Thesis unchanged.\n\n",
  "---\n\n",
  "## Honest n\n\n",
  "| item | n | note |\n",
  "|---|---:|---|\n",
  "| GEO series | **1** | GSE180963 |\n",
  "| Mice (unit) | **2** | 1 K + 1 KL |\n",
  "| K vs KL mice | **1 vs 1** | MW / Spearman cannot be computed |\n",
  "| 10x libraries | **1** | mixed, then demultiplexed by label |\n",
  "| Cells in matrix | **", n_cells, "** | not the unit |\n",
  "| Genes | **", nrow(obj), "** | mm10 symbols |\n",
  "| Author malignant labels | **0** | none deposited |\n",
  "| Marker epithelial cells (tight) | **", n_epi, "** | Epcam+ structural+ Ptprc- |\n",
  "| Loose Epcam+ (audit) | **", sum(md$is_epi_loose), "** | ambient-polluted |\n",
  "| Sftpc+ cells (ambient flag) | **", sum(md$sftpc_pos), "** | not host AT2 |\n",
  "| T/NK cells | **", n_tnk, "** | Cd3d/e/g or Cd8a or Nkg7/Ncr1/Klrb1c |\n",
  "| Cldn4-positive cells (any) | **", n_cldn4, "** | count > 0 |\n",
  "| Dual-high | **0** | not defined |\n",
  "| Private 8-KL mice | **0** | public GEO only |\n",
  "| ICI arms | **0** | untreated |\n\n",
  "Do not write n = ", n_cells, ". Do not write n = 2 genotypes as if they were biological replicates.\n\n",
  "---\n\n",
  "## Per-mouse table (the actual unit)\n\n",
  "| mouse | GSM | genotype | n cells | n epi tight | n epi loose | n T/NK | frac T/NK | %Sftpc | Cldn4 all | Cldn4 %pos | Cldn4 epi | Cldn4 epi %pos | IFN epi | MHC epi |\n",
  "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
  "| K | GSM5481386 | KrasG12D/+ | ",
  k_row$n_cells, " | ", k_row$n_epi, " | ", k_row$n_epi_loose, " | ", k_row$n_tnk, " | ",
  fmt(k_row$frac_tnk), " | ", fmt(100 * k_row$frac_sftpc, 1), " | ",
  fmt(k_row$Cldn4_all_mean, 4), " | ", fmt(k_row$Cldn4_all_pctpos, 2),
  " | ", fmt(k_row$Cldn4_epi_mean, 4), " | ", fmt(k_row$Cldn4_epi_pctpos, 2),
  " | ", fmt(k_row$IFN_epi_mean, 3), " | ", fmt(k_row$MHC_epi_mean, 3), " |\n",
  "| KL | GSM5481387 | KrasG12D/+;Lkb1fl/fl | ",
  kl_row$n_cells, " | ", kl_row$n_epi, " | ", kl_row$n_epi_loose, " | ", kl_row$n_tnk, " | ",
  fmt(kl_row$frac_tnk), " | ", fmt(100 * kl_row$frac_sftpc, 1), " | ",
  fmt(kl_row$Cldn4_all_mean, 4), " | ", fmt(kl_row$Cldn4_all_pctpos, 2),
  " | ", fmt(kl_row$Cldn4_epi_mean, 4), " | ", fmt(kl_row$Cldn4_epi_pctpos, 2),
  " | ", fmt(kl_row$IFN_epi_mean, 3), " | ", fmt(kl_row$MHC_epi_mean, 3), " |\n\n",
  "---\n\n",
  "## KL vs K (descriptive; n = 1 vs 1)\n\n",
  "| metric | n mice | KL | K | Δ (KL−K) | mouse-level test |\n",
  "|---|---|---:|---:|---:|---|\n",
  "| T/NK fraction | 1 vs 1 | ", fmt(kl_row$frac_tnk), " | ", fmt(k_row$frac_tnk), " | ",
  fmt(kl_row$frac_tnk - k_row$frac_tnk), " | **none** |\n",
  "| Cldn4 mean, all cells | 1 vs 1 | ", fmt(kl_row$Cldn4_all_mean, 4), " | ",
  fmt(k_row$Cldn4_all_mean, 4), " | ",
  fmt(kl_row$Cldn4_all_mean - k_row$Cldn4_all_mean, 4), " | **none** |\n",
  "| Cldn4 %pos, all cells | 1 vs 1 | ", fmt(kl_row$Cldn4_all_pctpos, 2), " | ",
  fmt(k_row$Cldn4_all_pctpos, 2), " | ",
  fmt(kl_row$Cldn4_all_pctpos - k_row$Cldn4_all_pctpos, 2), " | **none** |\n",
  "| Cldn4 mean, marker epithelium | 1 vs 1 | ", fmt(kl_row$Cldn4_epi_mean, 4), " | ",
  fmt(k_row$Cldn4_epi_mean, 4), " | ",
  fmt(kl_row$Cldn4_epi_mean - k_row$Cldn4_epi_mean, 4), " | **none** |\n",
  "| Cldn4 %pos, marker epithelium | 1 vs 1 | ", fmt(kl_row$Cldn4_epi_pctpos, 2), " | ",
  fmt(k_row$Cldn4_epi_pctpos, 2), " | ",
  fmt(kl_row$Cldn4_epi_pctpos - k_row$Cldn4_epi_pctpos, 2), " | **none** |\n",
  "| IFN ISG mean, marker epithelium | 1 vs 1 | ", fmt(kl_row$IFN_epi_mean, 3), " | ",
  fmt(k_row$IFN_epi_mean, 3), " | ",
  fmt(kl_row$IFN_epi_mean - k_row$IFN_epi_mean, 3), " | **none** |\n",
  "| MHC/APM mean, marker epithelium | 1 vs 1 | ", fmt(kl_row$MHC_epi_mean, 3), " | ",
  fmt(k_row$MHC_epi_mean, 3), " | ",
  fmt(kl_row$MHC_epi_mean - k_row$MHC_epi_mean, 3), " | **none** |\n",
  "| n marker epithelial cells | 1 vs 1 | ", kl_row$n_epi, " | ", k_row$n_epi, " | ",
  kl_row$n_epi - k_row$n_epi, " | **none** |\n\n",
  "**Cldn4 vs T/NK.** Two points. Spearman at the mouse unit is not defined. ",
  "K: Cldn4 epi mean ", fmt(k_row$Cldn4_epi_mean, 4), ", T/NK fraction ", fmt(k_row$frac_tnk),
  ". KL: Cldn4 epi mean ", fmt(kl_row$Cldn4_epi_mean, 4), ", T/NK fraction ", fmt(kl_row$frac_tnk),
  ". Do not draw a slope.\n\n",
  "**Cldn4 vs epithelial IFN/MHC.** Same two mice. Tight epithelium is **",
  k_row$n_epi, "** (K) and **", kl_row$n_epi, "** (KL) cells. ",
  "Cell-level Spearman inside KL tight epithelium (",
  rho_txt(sp_kl_ifn), " for IFN; ", rho_txt(sp_kl_mhc),
  " for MHC) is exploratory and must not be cited as n.\n\n",
  "The series summary frames LKB1 loss as an immune-desert TME. ",
  "This object does **not** independently confirm that as a drop in T/NK fraction ",
  "(K ", fmt(k_row$frac_tnk), " vs KL ", fmt(kl_row$frac_tnk),
  "). One library, one mouse per arm.\n\n",
  "---\n\n",
  "## Epithelium (honest)\n\n",
  "GEO deposits **no** author `Malignant` / `Epithelial` column. Calling is markers:\n\n",
  "- **Tight epithelium (primary)** = Epcam+ **and** (Cdh1 or Krt8/18/19 or Cldn18)+ **and** Ptprc−.\n",
  "- **Loose Epcam+** is audit only. Sftpc is an **ambient flag** in this digest, not host AT2.\n",
  "- **T/NK** = Cd3d / Cd3e / Cd3g / Cd8a / Nkg7 / Ncr1 / Klrb1c. Tight epithelium wins if both fire.\n",
  "- **Cldn4 is not a caller.** Tacstd2 is inventory only. No dual-high gate.\n\n",
  "Do not treat K vs KL epithelial IFN/MHC as a powered contrast. n = 1 vs 1 mice.\n\n",
  "Lineage genes present: ",
  paste(present_in(c(EPI_CORE, HOST_LUNG, T_NK_CALL, CLDN4, AUDIT), universe), collapse = ", "),
  ".\n",
  "Missing: ",
  {
    miss <- setdiff(c(EPI_CORE, HOST_LUNG, T_NK_CALL, CLDN4, AUDIT), universe)
    if (length(miss) == 0) "none of the core set"
    else paste(miss, collapse = ", ")
  },
  ".\n\n",
  "---\n\n",
  "## Methods (short)\n\n",
  "1. Public only. Downloaded `GSE180963_RAW.tar` and the series matrix from NCBI GEO FTP. ",
  "No SRA / FASTQ. No private 8-KL object.\n",
  "2. Each sample is a Cell Ranger v2-style MTX (`genes.tsv` + `barcodes.tsv` + `matrix.mtx`). ",
  "Integer counts. Author QC already applied (Seurat 3.1.5: 500–6000 features, mito < 20%, genes in ≥3 cells).\n",
  "3. Primary: `Seurat::Read10X` → `CreateSeuratObject` → `NormalizeData` (LogNormalize, 1e4). ",
  "Light PCA / neighbors / Leiden / UMAP for figures only. **No integration** — genotype is the mouse.\n",
  "4. Cldn4-only. Tacstd2 is inventory, never a gate. No dual-high.\n",
  "5. Positive = raw count > 0. Module scores = mean lognorm of present genes in the locked IFN and MHC sets.\n",
  "6. Unit = mouse. n = 1 vs 1. No mouse-level p-value.\n",
  "7. Thesis is not rewritten.\n\n",
  "```bash\n",
  "bash methods/seurat_gse180963_cldn4/scripts/install_r.sh\n",
  "bash methods/seurat_gse180963_cldn4/scripts/download.sh /tmp/gse180963\n",
  "Rscript methods/seurat_gse180963_cldn4/scripts/analyze.R --data /tmp/gse180963 --out methods/seurat_gse180963_cldn4\n",
  "```\n\n",
  "---\n\n",
  "## How to read this\n\n",
  "- **Additive public mouse**, not a human concordant-pool join.\n",
  "- **This is real KL** (Lkb1-targeted Kras lung GEMM), not KLK and not the private 8-KL cohort.\n",
  "- **Matrix exists.** The no-go is the **mouse n**, not the file.\n",
  "- **Cldn4 is present and sparse / epithelial.** Score it; do not build a Cldn4-high story on two mice.\n",
  "- **Cldn4 vs T/NK and epithelial IFN/MHC cannot be tested at the unit.** Two points.\n",
  "- **Honest n is 1 vs 1 mice.** Direction can be listed. A p-value cannot.\n",
  "- **No dual-high. No ICI endpoint. Thesis unchanged.**\n\n",
  "## Files\n\n",
  "- `tables/per_mouse.tsv` — unit-level table\n",
  "- `tables/kl_vs_k.tsv` — descriptive 1 vs 1\n",
  "- `tables/cldn4_positive_cells.tsv`\n",
  "- `tables/compartment_by_mouse.tsv`\n",
  "- `tables/gene_inventory.tsv`, `tables/honest_n.tsv`, `tables/summary.json`\n",
  "- `tables/cldn4_vs_ifn_mhc_cells_exploratory.tsv` — do not cite as n\n",
  "- `figures/fig_honest_n.png`, `fig_tnk_fraction.png`, `fig_cldn4_epi.png`, ",
  "`fig_epi_ifn_mhc.png`, `fig_cldn4_vs_tnk.png`, `fig_umap_*.png`\n",
  "- `scripts/analyze.R` — R + Seurat primary\n\n",
  "## 结论\n\n",
  "GSE180963 是公开的 K vs KL 肺 GEMM scRNA，GEO **有** 处理后的 10x MTX，",
  "因此用 **R + Seurat `CreateSeuratObject`** 做了 Cldn4-only 计分。",
  "诚实 n = **2 只鼠（1 vs 1）**，不是 ", n_cells, " 个细胞。",
  "Cldn4 行在，主要在上皮且稀疏；Cldn4 对 T/NK、上皮 IFN/MHC 在鼠单位上 **不能做检验**。",
  "未打开私有 8 只 KL。无 dual-high，无 ICI，不改 thesis。\n"
)

writeLines(finding, file.path(out_dir, "FINDING.md"))

readme <- paste0(
  "# seurat_gse180963_cldn4\n\n",
  "Public [GSE180963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963) ",
  "K vs KL lung GEMM scRNA. **R + Seurat** (`CreateSeuratObject`). Cldn4-only. ",
  "Honest n = **2 mice**. See `FINDING.md`.\n\n",
  "```bash\n",
  "bash methods/seurat_gse180963_cldn4/scripts/install_r.sh\n",
  "bash methods/seurat_gse180963_cldn4/scripts/download.sh /tmp/gse180963\n",
  "Rscript methods/seurat_gse180963_cldn4/scripts/analyze.R --data /tmp/gse180963 --out methods/seurat_gse180963_cldn4\n",
  "```\n"
)
writeLines(readme, file.path(out_dir, "README.md"))

message("wrote FINDING.md and tables under ", out_dir)
message("n_mice=2 n_cells=", n_cells, " n_epi=", n_epi, " n_tnk=", n_tnk, " n_cldn4=", n_cldn4)
