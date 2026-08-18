#!/usr/bin/env Rscript
# GSE165641 KL GEMM scRNA — Cldn4-only, Seurat CreateSeuratObject, mouse-level.
# Public Cell Ranger filtered matrices only. No private 8-KL. No Python primary.
# Honest n = mice (n=2). Do not overclaim.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
data_root <- if (length(args) >= 1) args[[1]] else "/tmp/geo_gse165641"
out_dir   <- if (length(args) >= 2) args[[2]] else "methods/seurat_gse165641_cldn4"
fig_dir   <- file.path(out_dir, "figures")
tab_dir   <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ---- locked gene sets (Cldn4-only; mouse symbols) ----
CLDN4 <- "Cldn4"

T_NK <- c(
  "Cd3d", "Cd3e", "Cd3g", "Cd2", "Cd8a", "Cd8b1", "Cd4",
  "Nkg7", "Gzma", "Gzmb", "Prf1", "Klrb1c", "Ncr1", "Klrd1", "Klrc1", "Ifng"
)

IFN <- c(
  "Stat1", "Stat2", "Irf1", "Irf7", "Irf9", "Isg15",
  "Ifit1", "Ifit2", "Ifit3", "Mx1", "Oasl2", "Rsad2", "Ifih1", "Ddx58", "Ifnb1"
)

MHC <- c(
  "B2m", "H2-K1", "H2-D1", "H2-Q4", "H2-Q6", "H2-Q7",
  "H2-Aa", "H2-Ab1", "H2-Eb1", "Tap1", "Tap2", "Psmb8", "Psmb9", "Nlrc5", "Ciita"
)

IFN_MHC <- c(IFN, MHC)

# Compartment markers (annotation only; not endpoints)
EPI_MARKERS <- c("Epcam", "Krt8", "Krt18", "Krt19", "Cdh1", "Nkx2-1", "Sftpc", "Scgb1a1")
MYE_MARKERS <- c("Lyz2", "Csf1r", "Cd68", "Itgam", "Fcgr1", "Adgre1", "C1qa")
B_MARKERS   <- c("Ms4a1", "Cd79a", "Cd79b", "Cd19")
NEU_MARKERS <- c("S100a8", "S100a9", "Ly6g", "Csf3r")
END_MARKERS <- c("Pecam1", "Cldn5", "Cdh5", "Kdr")
FIB_MARKERS <- c("Col1a1", "Dcn", "Pdgfra", "Col3a1")

present <- function(genes, features) intersect(genes, features)

mean_expr <- function(obj, genes, layer = "data") {
  genes <- present(genes, rownames(obj))
  if (length(genes) == 0) return(rep(NA_real_, ncol(obj)))
  mat <- GetAssayData(obj, layer = layer)[genes, , drop = FALSE]
  if (length(genes) == 1) as.numeric(mat) else Matrix::colMeans(mat)
}

any_pos <- function(obj, genes, layer = "counts") {
  genes <- present(genes, rownames(obj))
  if (length(genes) == 0) return(rep(FALSE, ncol(obj)))
  mat <- GetAssayData(obj, layer = layer)[genes, , drop = FALSE]
  Matrix::colSums(mat > 0) > 0
}

# ---- 1. Read10X + CreateSeuratObject per mouse ----
kl1_dir <- file.path(data_root, "KL1_count", "filtered_feature_bc_matrix")
kl2_dir <- file.path(data_root, "KL2_count", "filtered_feature_bc_matrix")
stopifnot(dir.exists(kl1_dir), dir.exists(kl2_dir))

kl1_counts <- Read10X(kl1_dir)
kl2_counts <- Read10X(kl2_dir)
# Ensembl + symbol; use gene symbols (second column already used by Read10X)
if (is.list(kl1_counts)) kl1_counts <- kl1_counts[[1]]
if (is.list(kl2_counts)) kl2_counts <- kl2_counts[[1]]

kl1 <- CreateSeuratObject(
  counts = kl1_counts,
  project = "GSE165641",
  min.cells = 0,
  min.features = 0
)
kl1$mouse    <- "KL1"
kl1$gsm      <- "GSM5047302"
kl1$genotype <- "KL"
kl1$orig.ident <- "KL1"

kl2 <- CreateSeuratObject(
  counts = kl2_counts,
  project = "GSE165641",
  min.cells = 0,
  min.features = 0
)
kl2$mouse    <- "KL2"
kl2$gsm      <- "GSM5047303"
kl2$genotype <- "KL"
kl2$orig.ident <- "KL2"

obj <- merge(kl1, y = kl2, add.cell.ids = c("KL1", "KL2"))
obj <- JoinLayers(obj)

# Mouse mito prefix
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

# Mild QC: drop empty / extreme debris. Thresholds chosen after seeing Cell Ranger
# filtered matrices (already cell-called: 4400 + 2780). Keep liberal.
n_before <- ncol(obj)
obj <- subset(obj, subset = nFeature_RNA >= 200 & nCount_RNA >= 500 & percent.mt < 25)
n_after <- ncol(obj)

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)

features <- rownames(obj)
stopifnot(CLDN4 %in% features)

# ---- 2. Compartment scores + assignment ----
obj <- AddModuleScore(obj, features = list(present(EPI_MARKERS, features)), name = "epi_mod", ctrl = 20)
obj <- AddModuleScore(obj, features = list(present(T_NK, features)),         name = "tnk_mod", ctrl = 20)
obj <- AddModuleScore(obj, features = list(present(MYE_MARKERS, features)),  name = "mye_mod", ctrl = 20)
obj <- AddModuleScore(obj, features = list(present(B_MARKERS, features)),    name = "b_mod",   ctrl = 20)
obj <- AddModuleScore(obj, features = list(present(NEU_MARKERS, features)),  name = "neu_mod", ctrl = 20)
obj <- AddModuleScore(obj, features = list(present(END_MARKERS, features)),  name = "end_mod", ctrl = 20)
obj <- AddModuleScore(obj, features = list(present(FIB_MARKERS, features)),  name = "fib_mod", ctrl = 20)

score_cols <- c(
  epithelial = "epi_mod1",
  T_NK       = "tnk_mod1",
  myeloid    = "mye_mod1",
  B          = "b_mod1",
  neutrophil = "neu_mod1",
  endothelial = "end_mod1",
  fibroblast = "fib_mod1"
)
score_mat <- as.matrix(obj[[unname(score_cols)]])
colnames(score_mat) <- names(score_cols)
best <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]
best_score <- apply(score_mat, 1, max)
# Require a positive module score; otherwise "other"
obj$compartment <- ifelse(best_score > 0.1, best, "other")

# Marker overrides for high-confidence T/NK and epithelium
obj$epcam_pos <- any_pos(obj, "Epcam")
obj$ptprc_pos <- any_pos(obj, "Ptprc")
obj$tnk_marker_pos <- any_pos(obj, c("Cd3d", "Cd3e", "Cd3g", "Nkg7", "Ncr1", "Klrb1c", "Cd8a"))
# If a cell is Epcam+ and not Ptprc+, force epithelium (Cldn4 lives here)
obj$compartment[obj$epcam_pos & !obj$ptprc_pos] <- "epithelial"
# If a cell is T/NK-marker+ and Ptprc+ and not Epcam+, force T/NK
obj$compartment[obj$tnk_marker_pos & obj$ptprc_pos & !obj$epcam_pos] <- "T_NK"

obj$cldn4 <- as.numeric(GetAssayData(obj, layer = "data")[CLDN4, ])
obj$cldn4_counts <- as.numeric(GetAssayData(obj, layer = "counts")[CLDN4, ])
obj$cldn4_pos <- obj$cldn4_counts > 0
obj$ifn_score <- mean_expr(obj, IFN)
obj$mhc_score <- mean_expr(obj, MHC)
obj$ifn_mhc_score <- mean_expr(obj, IFN_MHC)
obj$tnk_gene_score <- mean_expr(obj, T_NK)

# ---- 3. Mouse-level table (honest n = 2) ----
# Seurat [[ ]] returns the cell metadata data.frame (avoid @ slot syntax).
md <- obj[[]]
mice <- sort(unique(md$mouse))

mouse_rows <- lapply(mice, function(m) {
  d <- md[md$mouse == m, ]
  epi <- d[d$compartment == "epithelial", ]
  tnk <- d[d$compartment == "T_NK", ]
  data.frame(
    mouse = m,
    gsm = unique(d$gsm),
    genotype = "KL",
    n_cells_qc = nrow(d),
    n_epithelial = nrow(epi),
    n_T_NK = nrow(tnk),
    frac_T_NK = if (nrow(d) > 0) nrow(tnk) / nrow(d) else NA_real_,
    frac_epithelial = if (nrow(d) > 0) nrow(epi) / nrow(d) else NA_real_,
    cldn4_mean_all = mean(d$cldn4),
    cldn4_pct_pos_all = mean(d$cldn4_pos),
    cldn4_mean_epithelial = if (nrow(epi) > 0) mean(epi$cldn4) else NA_real_,
    cldn4_pct_pos_epithelial = if (nrow(epi) > 0) mean(epi$cldn4_pos) else NA_real_,
    ifn_mean_epithelial = if (nrow(epi) > 0) mean(epi$ifn_score) else NA_real_,
    mhc_mean_epithelial = if (nrow(epi) > 0) mean(epi$mhc_score) else NA_real_,
    ifn_mhc_mean_epithelial = if (nrow(epi) > 0) mean(epi$ifn_mhc_score) else NA_real_,
    tnk_gene_score_all = mean(d$tnk_gene_score),
    stringsAsFactors = FALSE
  )
})
mouse_tab <- do.call(rbind, mouse_rows)

# Spearman on n=2 is always |r|=1 or undefined. Report the two values, not a p.
if (nrow(mouse_tab) == 2 && all(is.finite(mouse_tab$cldn4_mean_epithelial))) {
  d_cldn4 <- diff(mouse_tab$cldn4_mean_epithelial)
  d_tnk   <- diff(mouse_tab$frac_T_NK)
  d_ifn   <- diff(mouse_tab$ifn_mhc_mean_epithelial)
  mouse_tab$note <- "n=2; Spearman not reported (always |r|=1)"
} else {
  d_cldn4 <- d_tnk <- d_ifn <- NA_real_
}

# Cluster audit
cluster_tab <- do.call(rbind, lapply(sort(unique(md$seurat_clusters)), function(cl) {
  d <- md[md$seurat_clusters == cl, ]
  data.frame(
    cluster = as.character(cl),
    n = nrow(d),
    n_KL1 = sum(d$mouse == "KL1"),
    n_KL2 = sum(d$mouse == "KL2"),
    pct_Epcam = 100 * mean(d$epcam_pos),
    pct_Ptprc = 100 * mean(d$ptprc_pos),
    pct_Cldn4 = 100 * mean(d$cldn4_pos),
    pct_TNK_marker = 100 * mean(d$tnk_marker_pos),
    mean_Cldn4 = mean(d$cldn4),
    mean_IFN_MHC = mean(d$ifn_mhc_score),
    majority_compartment = names(sort(table(d$compartment), decreasing = TRUE))[1],
    stringsAsFactors = FALSE
  )
}))

# Genes used
feat <- rownames(obj)
genes_used <- data.frame(
  set = c(
    rep("Cldn4", 1),
    rep("T_NK", length(T_NK)),
    rep("IFN", length(IFN)),
    rep("MHC", length(MHC))
  ),
  gene = c(CLDN4, T_NK, IFN, MHC),
  present = c(CLDN4, T_NK, IFN, MHC) %in% feat,
  stringsAsFactors = FALSE
)

compartment_tab <- as.data.frame(table(md$mouse, md$compartment), stringsAsFactors = FALSE)
colnames(compartment_tab) <- c("mouse", "compartment", "n")

# Session
sess <- list(
  seurat = as.character(packageVersion("Seurat")),
  seuratobject = as.character(packageVersion("SeuratObject")),
  r = paste(R.version$major, R.version$minor, sep = "."),
  n_cells_cellranger = n_before,
  n_cells_qc = n_after,
  n_mice = length(mice),
  create_seurat_object = TRUE,
  python_primary = FALSE,
  private_8kl = FALSE,
  accession = "GSE165641",
  processed_matrix = "GSM5047302_KL1_count.tar.gz + GSM5047303_KL2_count.tar.gz (Cell Ranger filtered MTX) ; also GSE165641_processed_normalized_matrix_data.Rdata.gz (31,053 x 7,180 gene x cell data.frame, not used as primary counts)"
)

# ---- write tables ----
write.table(mouse_tab, file.path(tab_dir, "mouse_level_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cluster_tab, file.path(tab_dir, "cluster_audit.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes_used, file.path(tab_dir, "genes_used.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(compartment_tab, file.path(tab_dir, "compartment_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(
  item = c("GEO_QC_cells_cellranger", "cells_after_QC", "biological_mice",
           "KL_mice", "Spearman_allowed", "high_vs_low_allowed",
           "epithelial_compartment", "T_NK_compartment"),
  n_or_flag = c(
    n_before, n_after, 2, 2,
    "NO — n=2 forces |r|=1",
    "NO — 1 vs 1 is not a split",
    as.integer(any(md$compartment == "epithelial")),
    as.integer(any(md$compartment == "T_NK"))
  ),
  stringsAsFactors = FALSE
), file.path(tab_dir, "honest_n.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

jsonlite_ok <- requireNamespace("jsonlite", quietly = TRUE)
if (jsonlite_ok) {
  jsonlite::write_json(sess, file.path(tab_dir, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
} else {
  writeLines(paste(names(sess), unlist(sess), sep = "\t"), file.path(tab_dir, "summary.json"))
}

# ---- figures ----
pdf(file.path(fig_dir, "umap_compartment_cldn4.pdf"), width = 10, height = 8)
p1 <- DimPlot(obj, group.by = "compartment", split.by = "mouse", label = TRUE) +
  ggplot2::ggtitle("GSE165641 KL GEMM — compartment (n=2 mice)")
p2 <- FeaturePlot(obj, features = "Cldn4", split.by = "mouse", order = TRUE) +
  ggplot2::ggtitle("Cldn4 (log-norm)")
print(p1)
print(p2)
print(DimPlot(obj, group.by = "seurat_clusters", split.by = "mouse", label = TRUE))
print(VlnPlot(obj, features = c("Cldn4", "Epcam", "Cd3d", "Nkg7", "Ptprc"),
              group.by = "compartment", pt.size = 0, ncol = 3))
dev.off()

png(file.path(fig_dir, "umap_compartment.png"), width = 1400, height = 600, res = 120)
print(p1)
dev.off()
png(file.path(fig_dir, "umap_cldn4.png"), width = 1400, height = 600, res = 120)
print(p2)
dev.off()

# Mouse-level two-point plot (not a correlation claim)
png(file.path(fig_dir, "mouse_two_point.png"), width = 900, height = 420, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(mouse_tab$cldn4_mean_epithelial, mouse_tab$frac_T_NK,
     pch = 19, cex = 1.6, col = c("#1f77b4", "#d62728"),
     xlab = "epithelial Cldn4 mean (log-norm)",
     ylab = "T/NK fraction of QC cells",
     main = "n=2 mice — not a Spearman")
text(mouse_tab$cldn4_mean_epithelial, mouse_tab$frac_T_NK,
     labels = mouse_tab$mouse, pos = 3, cex = 0.9)
plot(mouse_tab$cldn4_mean_epithelial, mouse_tab$ifn_mhc_mean_epithelial,
     pch = 19, cex = 1.6, col = c("#1f77b4", "#d62728"),
     xlab = "epithelial Cldn4 mean (log-norm)",
     ylab = "epithelial IFN/MHC mean (log-norm)",
     main = "n=2 mice — not a Spearman")
text(mouse_tab$cldn4_mean_epithelial, mouse_tab$ifn_mhc_mean_epithelial,
     labels = mouse_tab$mouse, pos = 3, cex = 0.9)
par(op)
dev.off()

# Per-cell metadata (small enough)
write.table(
  md[, c("mouse", "gsm", "seurat_clusters", "compartment",
         "nFeature_RNA", "nCount_RNA", "percent.mt",
         "cldn4", "cldn4_pos", "ifn_score", "mhc_score", "ifn_mhc_score",
         "tnk_gene_score", "epcam_pos", "ptprc_pos", "tnk_marker_pos")],
  file.path(tab_dir, "cell_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE
)

cat("WROTE tables to", tab_dir, "\n")
print(mouse_tab)
cat("cells before/after QC:", n_before, n_after, "\n")
cat("compartments:\n")
print(table(md$mouse, md$compartment))
cat("Seurat", as.character(packageVersion("Seurat")),
    "CreateSeuratObject used. Honest n=2 mice.\n")
