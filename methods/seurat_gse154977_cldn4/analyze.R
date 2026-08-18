#!/usr/bin/env Rscript
# ADDITIVE public mouse: GSE154977 Marjanovic KP 30w 10x, Cldn4-only.
# Primary is R + Seurat. T/NK is expected design no-go (FACS CD45−).
# Thesis is not rewritten. No dual-high. No private 8 KL. No Python-only primary.

suppressPackageStartupMessages({
  library(hdf5r)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

HERE <- {
  args_all <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
  if (length(f) == 1 && nzchar(f)) {
    dirname(normalizePath(f))
  } else if (dir.exists(file.path(getwd(), "methods/seurat_gse154977_cldn4"))) {
    normalizePath(file.path(getwd(), "methods/seurat_gse154977_cldn4"))
  } else {
    normalizePath(getwd())
  }
}

DATA <- file.path(HERE, "data")
TABLES <- file.path(HERE, "tables")
FIGS <- file.path(HERE, "figures")
OBJS <- file.path(HERE, "objects")
dir.create(DATA, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLES, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGS, showWarnings = FALSE, recursive = TRUE)
dir.create(OBJS, showWarnings = FALSE, recursive = TRUE)

GEO_BASE <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154977/suppl"
FILES <- list(
  raw = "GSE154977_mmLung10x_cis_dSp_rawCount.h5",
  smp = "GSE154977_mmLung10x_cis_smpTable.csv.gz",
  gene = "GSE154977_mmLung10x_cis_geneTable.csv.gz",
  annot = "GSE154977_mmLung10x_cis_dZ_annot_annot_smpTable.csv.gz",
  qc = "GSE154977_mmLung10x_cis_dZ_QCstat_QCstat_smpTable.csv.gz"
)

MIN_CELLS_PRIMARY <- 20L
MIN_N_SPEARMAN <- 4L
MIN_N_Q4 <- 8L

MOUSE_MHC_I_APM <- c(
  "H2-K1", "H2-D1", "H2-Q1", "H2-Q2", "H2-Q4", "H2-Q6", "H2-Q7",
  "H2-Q10", "H2-T23", "H2-M3", "B2m", "Tap1", "Tap2", "Tapbp",
  "Tapbpl", "Nlrc5", "Psmb8", "Psmb9", "Psmb10", "Erap1", "Calr",
  "Canx", "Pdia3", "Irf1"
)

HUMAN_TO_MOUSE <- c(
  WARS1 = "Wars", C1R = "C1ra", C1S = "C1s1", FCGR1A = "Fcgr1",
  `HLA-A` = "H2-K1", `HLA-B` = "H2-D1", `HLA-C` = "H2-Q4",
  `HLA-E` = "H2-T23", `HLA-F` = "H2-Q10", `HLA-G` = "H2-Q6",
  `HLA-DMA` = "H2-DMa", `HLA-DQA1` = "H2-Aa", `HLA-DRB1` = "H2-Eb1",
  SECTM1 = "Sectm1a", TENT5A = "Tent5a", MARCHF1 = "Marchf1"
)

TNK_AUDIT <- c("Ptprc", "Cd3d", "Cd3e", "Cd3g", "Cd2", "Cd8a", "Cd8b1", "Nkg7", "Klrd1", "Gzma")
EPI_AUDIT <- c("Cldn4", "Epcam", "Cdh1", "Krt8", "Krt18", "Sftpc", "Nkx2-1")

download_if_needed <- function(name) {
  dest <- file.path(DATA, name)
  if (file.exists(dest) && file.info(dest)$size > 1000) return(dest)
  url <- paste0(GEO_BASE, "/", name)
  message("DOWNLOAD ", url)
  utils::download.file(url, dest, mode = "wb", quiet = FALSE)
  dest
}

human_to_mouse <- function(sym) {
  if (sym %in% names(HUMAN_TO_MOUSE)) return(unname(HUMAN_TO_MOUSE[[sym]]))
  if (startsWith(sym, "HLA-")) return("")
  clean <- gsub("[-_]", "", sym)
  if (identical(clean, toupper(clean))) {
    parts <- unlist(strsplit(sym, "(?<=-)|(?=-)", perl = TRUE))
    out <- vapply(parts, function(p) {
      if (p == "-") return(p)
      paste0(substr(p, 1, 1), tolower(substr(p, 2, nchar(p))))
    }, character(1))
    return(paste(out, collapse = ""))
  }
  sym
}

load_families <- function(available) {
  a8 <- jsonlite::fromJSON(file.path(DATA, "a8_families.json"), simplifyVector = TRUE)
  sets <- a8$sets
  ifn_h <- unique(c(sets$HALLMARK_INTERFERON_GAMMA_RESPONSE, sets$HALLMARK_INTERFERON_ALPHA_RESPONSE))
  ifn <- character()
  seen <- character()
  for (g in sort(ifn_h)) {
    m <- human_to_mouse(g)
    if (nzchar(m) && m %in% available && !(m %in% seen)) {
      ifn <- c(ifn, m)
      seen <- c(seen, m)
    }
  }
  mhc <- MOUSE_MHC_I_APM[MOUSE_MHC_I_APM %in% available]
  tj_h <- unique(c(sets$KEGG_TIGHT_JUNCTION, sets$GOBP_TIGHT_JUNCTION_ORGANIZATION, a8$focal_genes))
  krt <- sets$KRT_EPITHELIAL
  tj_h <- setdiff(tj_h, c("CLDN4", krt))
  tj <- character()
  seen <- character()
  for (g in sort(tj_h)) {
    m <- human_to_mouse(g)
    if (identical(m, "Cldn4")) next
    if (nzchar(m) && m %in% available && !(m %in% seen)) {
      tj <- c(tj, m)
      seen <- c(seen, m)
    }
  }
  list(IFN = ifn, `MHC-I/APM` = mhc, TJ = tj)
}

read_coo_h5 <- function(path, genes, cells) {
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  read_vec <- function(name) {
    d <- h5[[name]]
    dims <- d$dims
    if (length(dims) == 1) {
      as.numeric(d[])
    } else if (length(dims) == 2) {
      as.numeric(d[, ])
    } else {
      stop("unexpected rank for ", name, ": ", paste(dims, collapse = "x"))
    }
  }
  i <- read_vec("i")
  j <- read_vec("j")
  v <- read_vec("v")
  if (min(i) == 0 || min(j) == 0) {
    i <- i + 1
    j <- j + 1
  }
  n_g <- length(genes)
  n_c <- length(cells)
  if (max(i) > n_g || max(j) > n_c) {
    stop(sprintf("COO index out of range: max i=%s (n_genes=%s) max j=%s (n_cells=%s)",
                 max(i), n_g, max(j), n_c))
  }
  sparseMatrix(
    i = as.integer(i),
    j = as.integer(j),
    x = v,
    dims = c(n_g, n_c),
    dimnames = list(genes, cells)
  )
}

spearman_safe <- function(x, y, min_n = MIN_N_SPEARMAN) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  out <- list(n = n, rho = NA_real_, p = NA_real_, usable = FALSE)
  if (n < min_n) return(out)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = n <= 10))
  out$rho <- unname(ct$estimate)
  out$p <- unname(ct$p.value)
  out$usable <- TRUE
  out
}

q4_vs_q1 <- function(cldn4, y, min_n = MIN_N_Q4) {
  ok <- is.finite(cldn4) & is.finite(y)
  n <- sum(ok)
  out <- list(
    n = n, n_q1 = NA_integer_, n_q4 = NA_integer_,
    median_q1 = NA_real_, median_q4 = NA_real_,
    delta_median = NA_real_, p = NA_real_, usable = FALSE
  )
  if (n < min_n) return(out)
  qs <- quantile(cldn4[ok], probs = c(0.25, 0.75), names = FALSE, type = 7)
  q1 <- cldn4[ok] <= qs[1]
  q4 <- cldn4[ok] >= qs[2]
  out$n_q1 <- sum(q1)
  out$n_q4 <- sum(q4)
  if (out$n_q1 < 2 || out$n_q4 < 2) return(out)
  y1 <- y[ok][q1]
  y4 <- y[ok][q4]
  out$median_q1 <- median(y1)
  out$median_q4 <- median(y4)
  out$delta_median <- out$median_q4 - out$median_q1
  out$p <- wilcox.test(y4, y1, exact = FALSE)$p.value
  out$usable <- TRUE
  out
}

fmt <- function(x, nd = 3) {
  if (length(x) == 0 || is.null(x) || (is.numeric(x) && !is.finite(x))) return("—")
  sprintf(paste0("%.", nd, "f"), x)
}

message("=== GSE154977 Seurat Cldn4-only ===")
invisible(lapply(FILES, download_if_needed))

smp <- read.csv(file.path(DATA, FILES$smp), stringsAsFactors = FALSE)
gene_tbl <- read.csv(file.path(DATA, FILES$gene), stringsAsFactors = FALSE)
annot <- read.csv(file.path(DATA, FILES$annot), stringsAsFactors = FALSE)
qcstat <- read.csv(file.path(DATA, FILES$qc), stringsAsFactors = FALSE)

if (anyDuplicated(gene_tbl$geneID)) {
  gene_tbl$geneID <- make.unique(gene_tbl$geneID)
}
genes <- gene_tbl$geneID
cells <- smp$sampleID
stopifnot(length(unique(cells)) == length(cells))

message("Reading rawCount COO h5 ...")
counts <- read_coo_h5(file.path(DATA, FILES$raw), genes, cells)
message(sprintf("matrix %s genes x %s cells; nnz=%s", nrow(counts), ncol(counts), nnzero(counts)))

meta <- data.frame(row.names = cells, sampleID = cells, stringsAsFactors = FALSE)
meta$barcode <- smp$barcode[match(cells, smp$sampleID)]
meta$library <- sub("_id-.*$", "", meta$sampleID)
meta$mouse <- sub("_PT$", "", meta$library)
meta$treatment <- ifelse(grepl("Cis72", meta$library), "Cis72", "ND")
meta$genotype <- "KP"
meta$week <- "30w"
meta$gsm <- NA_character_
meta$gsm[meta$library == "KP_30w_ND_m3_PT"] <- "GSM4685281"
meta$gsm[meta$library == "KP_30w_ND_m4_PT"] <- "GSM4685282"
meta$gsm[meta$library == "KP_30w_Cis72_m5_PT"] <- "GSM4685283"
meta$gsm[meta$library == "KP_30w_Cis72_m6_PT"] <- "GSM4685284"

annot <- annot[match(cells, annot$sampleID), ]
qcstat <- qcstat[match(cells, qcstat$sampleID), ]
meta$author_cluster <- as.character(annot$timecourse_pred_cluster)
meta$author_tSNE_X <- annot$tSNE_X
meta$author_tSNE_Y <- annot$tSNE_Y
meta$author_mitoPct <- qcstat$mitoPct

message("CreateSeuratObject ...")
obj <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  project = "GSE154977_KP30w",
  min.cells = 0,
  min.features = 0
)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
obj$nCount_RNA <- obj$nCount_RNA
obj$nFeature_RNA <- obj$nFeature_RNA

# Basic QC (author already dropped empty droplets; still apply a light Seurat filter)
qc_keep <- obj$nFeature_RNA >= 200 & obj$percent.mt < 20
message(sprintf(
  "QC keep %s / %s (nFeature>=200 & percent.mt<20); author mitoPct median=%.3f",
  sum(qc_keep), ncol(obj), median(obj$author_mitoPct, na.rm = TRUE)
))
obj$qc_pass <- qc_keep
obj <- subset(obj, subset = qc_pass)

message("Normalize / PCA / UMAP / clusters ...")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)

present <- rownames(obj)
families <- load_families(present)
write.table(
  data.frame(
    family = names(families),
    n_requested_or_mapped = vapply(families, length, integer(1)),
    stringsAsFactors = FALSE
  ),
  file.path(TABLES, "geneset_coverage.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

mod_input <- list(IFN = families$IFN, MHC = families$`MHC-I/APM`, TJ = families$TJ)
obj <- AddModuleScore(obj, features = mod_input, name = "mod_", search = FALSE)
mdn <- colnames(obj[[]])
mod_cols <- grep("^mod_", mdn, value = TRUE)
if (length(mod_cols) < 3) stop("AddModuleScore did not return 3 scores: ", paste(mod_cols, collapse = ","))
obj$score_IFN <- obj[[mod_cols[1]]][, 1]
obj$score_MHC <- obj[[mod_cols[2]]][, 1]
obj$score_TJ <- obj[[mod_cols[3]]][, 1]

mean_present <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  if (!length(g)) return(rep(NA_real_, ncol(mat)))
  if (length(g) == 1) return(as.numeric(mat[g, ]))
  Matrix::colMeans(mat[g, , drop = FALSE])
}

data_mat <- GetAssayData(obj, layer = "data")
obj$mean_IFN <- mean_present(data_mat, families$IFN)
obj$mean_MHC <- mean_present(data_mat, families$`MHC-I/APM`)
obj$mean_TJ <- mean_present(data_mat, families$TJ)
obj$Cldn4_log <- if ("Cldn4" %in% rownames(obj)) as.numeric(data_mat["Cldn4", ]) else NA_real_
obj$Cldn4_counts <- if ("Cldn4" %in% rownames(obj)) as.numeric(GetAssayData(obj, layer = "counts")["Cldn4", ]) else 0

# Marker audit (counts > 0)
count_mat <- GetAssayData(obj, layer = "counts")
audit_genes <- unique(c(TNK_AUDIT, EPI_AUDIT))
audit_n <- sapply(audit_genes, function(g) {
  if (!g %in% rownames(count_mat)) return(NA_integer_)
  as.integer(sum(count_mat[g, ] > 0))
})
audit_df <- data.frame(
  gene = audit_genes,
  n_pos = as.integer(audit_n),
  n_cells = ncol(obj),
  frac = as.numeric(audit_n) / ncol(obj),
  present_in_matrix = audit_genes %in% rownames(count_mat),
  stringsAsFactors = FALSE
)
write.table(audit_df, file.path(TABLES, "marker_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Author cluster audit
cl <- unique(obj$author_cluster)
cluster_audit <- do.call(rbind, lapply(sort(cl), function(k) {
  ii <- obj$author_cluster == k
  getf <- function(g) {
    if (!g %in% rownames(count_mat)) return(NA_real_)
    mean(count_mat[g, ii] > 0)
  }
  data.frame(
    author_cluster = k,
    n_cells = sum(ii),
    frac_Cldn4 = getf("Cldn4"),
    frac_Epcam = getf("Epcam"),
    frac_Krt8 = getf("Krt8"),
    frac_Ptprc = getf("Ptprc"),
    frac_Cd3d = getf("Cd3d"),
    frac_Nkg7 = getf("Nkg7"),
    stringsAsFactors = FALSE
  )
}))
write.table(cluster_audit, file.path(TABLES, "cluster_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Compartment call: author clusters 1–12 are tumor epithelial states (paper).
# Marker leftover T/NK: Cd3d|Cd3e|Nkg7 plus Ptprc.
tnk_any <- if (all(c("Cd3d", "Nkg7", "Ptprc") %in% rownames(count_mat))) {
  (count_mat["Ptprc", ] > 0) & ((count_mat["Cd3d", ] > 0) | (if ("Cd3e" %in% rownames(count_mat)) count_mat["Cd3e", ] > 0 else 0) | (count_mat["Nkg7", ] > 0))
} else {
  rep(FALSE, ncol(obj))
}
epi_any <- if ("Epcam" %in% rownames(count_mat)) count_mat["Epcam", ] > 0 else rep(TRUE, ncol(obj))
obj$tnk_leak <- as.logical(tnk_any)
obj$epcam_pos <- as.logical(epi_any)
obj$compartment <- "epithelial_tumor"
obj$compartment[obj$tnk_leak & !obj$epcam_pos] <- "TNK_leak"
# still not a real T/NK compartment if rate is tiny

tnk_n <- sum(obj$tnk_leak)
tnk_frac <- tnk_n / ncol(obj)
tnk_absent <- tnk_frac < 0.02
message(sprintf("T/NK leak (Ptprc & (Cd3d|Cd3e|Nkg7)): %s / %s (%.3f). absent_rule=%s",
                tnk_n, ncol(obj), tnk_frac, tnk_absent))

# Mouse-level table (PRIMARY unit)
md <- obj[[]]
mouse_tab <- do.call(rbind, lapply(split(seq_len(nrow(md)), md$mouse), function(ii) {
  d <- md[ii, ]
  data.frame(
    mouse = d$mouse[1],
    library = d$library[1],
    gsm = d$gsm[1],
    genotype = "KP",
    week = "30w",
    treatment = d$treatment[1],
    n_cells = nrow(d),
    n_Cldn4_pos = sum(d$Cldn4_counts > 0, na.rm = TRUE),
    frac_Cldn4_pos = mean(d$Cldn4_counts > 0, na.rm = TRUE),
    mean_Cldn4 = mean(d$Cldn4_log, na.rm = TRUE),
    mean_IFN = mean(d$mean_IFN, na.rm = TRUE),
    mean_MHC = mean(d$mean_MHC, na.rm = TRUE),
    mean_TJ = mean(d$mean_TJ, na.rm = TRUE),
    mod_IFN = mean(d$score_IFN, na.rm = TRUE),
    mod_MHC = mean(d$score_MHC, na.rm = TRUE),
    mod_TJ = mean(d$score_TJ, na.rm = TRUE),
    frac_Epcam = mean(d$epcam_pos),
    frac_TNK_leak = mean(d$tnk_leak),
    median_nFeature = median(d$nFeature_RNA),
    median_nCount = median(d$nCount_RNA),
    median_percent.mt = median(d$percent.mt),
    stringsAsFactors = FALSE
  )
}))
mouse_tab <- mouse_tab[order(mouse_tab$mouse), ]
write.table(mouse_tab, file.path(TABLES, "mouse_units.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

primary <- mouse_tab[mouse_tab$n_cells >= MIN_CELLS_PRIMARY, ]
write.table(primary, file.path(TABLES, "primary_KP_mice.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

family_rows <- list()
for (fam in c("IFN", "MHC", "TJ")) {
  ycol <- paste0("mean_", fam)
  sp <- spearman_safe(primary$mean_Cldn4, primary[[ycol]])
  q <- q4_vs_q1(primary$mean_Cldn4, primary[[ycol]])
  family_rows[[length(family_rows) + 1]] <- data.frame(
    cohort = "all_4_mice",
    family = fam,
    n = sp$n,
    rho = sp$rho,
    p = sp$p,
    spearman_usable = sp$usable,
    n_q1 = q$n_q1,
    n_q4 = q$n_q4,
    delta_median = q$delta_median,
    q_p = q$p,
    q_usable = q$usable,
    stringsAsFactors = FALSE
  )
}
# ND-only descriptive (n=2): Spearman locked off
nd <- primary[primary$treatment == "ND", ]
if (nrow(nd) >= 1) {
  for (fam in c("IFN", "MHC", "TJ")) {
    ycol <- paste0("mean_", fam)
    sp <- spearman_safe(nd$mean_Cldn4, nd[[ycol]])
    family_rows[[length(family_rows) + 1]] <- data.frame(
      cohort = "ND_only",
      family = fam,
      n = nrow(nd),
      rho = sp$rho,
      p = sp$p,
      spearman_usable = sp$usable,
      n_q1 = NA, n_q4 = NA, delta_median = NA, q_p = NA, q_usable = FALSE,
      stringsAsFactors = FALSE
    )
  }
}
family_tab <- do.call(rbind, family_rows)
write.table(family_tab, file.path(TABLES, "family_spearman.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

lib_qc <- do.call(rbind, lapply(split(seq_len(nrow(md)), md$library), function(ii) {
  d <- md[ii, ]
  data.frame(
    library = d$library[1],
    mouse = d$mouse[1],
    n_cells = nrow(d),
    median_nFeature = median(d$nFeature_RNA),
    median_nCount = median(d$nCount_RNA),
    median_percent.mt = median(d$percent.mt),
    stringsAsFactors = FALSE
  )
}))
write.table(lib_qc, file.path(TABLES, "library_qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

honest <- data.frame(
  item = c(
    "GEO QC cells (depositor)",
    "cells after Seurat QC",
    "libraries / GSM",
    "biological mice",
    "ND mice",
    "Cis72 mice",
    "PRIMARY mice (n_cells>=20)",
    "T/NK-fraction mice",
    "Q4 vs Q1 usable?"
  ),
  n = c(
    11017,
    ncol(obj),
    length(unique(md$library)),
    length(unique(md$mouse)),
    sum(primary$treatment == "ND"),
    sum(primary$treatment == "Cis72"),
    nrow(primary),
    0,
    0
  ),
  note = c(
    "do not quote as analysis n",
    "nFeature>=200 & percent.mt<20",
    "one library per mouse",
    "honest unit",
    "no-drug 30w PT",
    "cisplatin 72h 30w PT",
    "all 4 pass",
    "FACS CD45− tumor-only; design no-go",
    "locked off: needs n_units>=8"
  ),
  stringsAsFactors = FALSE
)
write.table(honest, file.path(TABLES, "honest_n.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# Plots — must be Seurat DimPlot / VlnPlot
message("Plots ...")
p_dim_lib <- DimPlot(obj, group.by = "library", reduction = "umap") + ggtitle("GSE154977 KP 30w 10x - library (mouse)")
p_dim_cl <- DimPlot(obj, group.by = "seurat_clusters", reduction = "umap", label = TRUE) + ggtitle("Seurat clusters")
p_dim_auth <- DimPlot(obj, group.by = "author_cluster", reduction = "umap", label = TRUE) + ggtitle("Author timecourse_pred_cluster")
p_dim_trt <- DimPlot(obj, group.by = "treatment", reduction = "umap") + ggtitle("Treatment")
ggsave(file.path(FIGS, "DimPlot_library.png"), p_dim_lib, width = 8, height = 6, dpi = 150)
ggsave(file.path(FIGS, "DimPlot_library.pdf"), p_dim_lib, width = 8, height = 6)
ggsave(file.path(FIGS, "DimPlot_seurat_clusters.png"), p_dim_cl, width = 8, height = 6, dpi = 150)
ggsave(file.path(FIGS, "DimPlot_seurat_clusters.pdf"), p_dim_cl, width = 8, height = 6)
ggsave(file.path(FIGS, "DimPlot_author_cluster.png"), p_dim_auth, width = 8, height = 6, dpi = 150)
ggsave(file.path(FIGS, "DimPlot_author_cluster.pdf"), p_dim_auth, width = 8, height = 6)
ggsave(file.path(FIGS, "DimPlot_treatment.png"), p_dim_trt, width = 7, height = 5.5, dpi = 150)
ggsave(file.path(FIGS, "DimPlot_treatment.pdf"), p_dim_trt, width = 7, height = 5.5)

feat <- c("Cldn4", "Epcam", "Ptprc", "Cd3d", "Nkg7")
feat <- feat[feat %in% rownames(obj)]
vln_feats <- paste0("rna_", feat)
p_vln_lib <- VlnPlot(obj, features = vln_feats, group.by = "library", pt.size = 0, ncol = 2) +
  plot_annotation(title = "GSE154977 - Cldn4 / epithelium / T-NK audit by mouse library")
p_vln_cl <- VlnPlot(
  obj,
  features = paste0("rna_", intersect(c("Cldn4", "Epcam", "Ptprc", "Cd3d"), rownames(obj))),
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 2
)
ggsave(file.path(FIGS, "VlnPlot_markers_by_library.png"), p_vln_lib, width = 10, height = 8, dpi = 150)
ggsave(file.path(FIGS, "VlnPlot_markers_by_library.pdf"), p_vln_lib, width = 10, height = 8)
ggsave(file.path(FIGS, "VlnPlot_markers_by_cluster.png"), p_vln_cl, width = 10, height = 7, dpi = 150)
ggsave(file.path(FIGS, "VlnPlot_markers_by_cluster.pdf"), p_vln_cl, width = 10, height = 7)

p_feat <- FeaturePlot(
  obj,
  features = paste0("rna_", intersect(c("Cldn4", "Epcam", "Ptprc"), rownames(obj))),
  ncol = 3
)
ggsave(file.path(FIGS, "FeaturePlot_Cldn4_Epcam_Ptprc.png"), p_feat, width = 12, height = 4, dpi = 150)
ggsave(file.path(FIGS, "FeaturePlot_Cldn4_Epcam_Ptprc.pdf"), p_feat, width = 12, height = 4)

# Mouse-level scatter
sc <- primary
sc$treat <- sc$treatment
p_sc <- ggplot(sc, aes(x = mean_Cldn4, y = mean_IFN, color = treat, label = mouse)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(title = "Mouse-level Cldn4 vs IFN (n=4 mice)", x = "mean log-norm Cldn4", y = "mean IFN family") +
  theme_bw()
p_sc2 <- ggplot(sc, aes(x = mean_Cldn4, y = mean_MHC, color = treat, label = mouse)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(title = "Mouse-level Cldn4 vs MHC-I/APM (n=4 mice)", x = "mean log-norm Cldn4", y = "mean MHC-I/APM") +
  theme_bw()
p_sc3 <- ggplot(sc, aes(x = mean_Cldn4, y = mean_TJ, color = treat, label = mouse)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(title = "Mouse-level Cldn4 vs TJ hold-out (n=4 mice)", x = "mean log-norm Cldn4", y = "mean TJ (Cldn4 held out)") +
  theme_bw()
p_mice <- p_sc + p_sc2 + p_sc3 + plot_layout(ncol = 3, guides = "collect")
ggsave(file.path(FIGS, "mouse_cldn4_vs_ifn_mhc_tj.png"), p_mice, width = 14, height = 4.5, dpi = 150)
ggsave(file.path(FIGS, "mouse_cldn4_vs_ifn_mhc_tj.pdf"), p_mice, width = 14, height = 4.5)

p_n <- ggplot(primary, aes(x = mouse, y = n_cells, fill = treatment)) +
  geom_col() +
  coord_flip() +
  labs(title = "Honest n = 4 KP mice (not 11,017 cells)", y = "cells after Seurat QC", x = NULL) +
  theme_bw()
ggsave(file.path(FIGS, "honest_n.png"), p_n, width = 7, height = 3.2, dpi = 150)
ggsave(file.path(FIGS, "honest_n.pdf"), p_n, width = 7, height = 3.2)

message("Saving Seurat object ...")
obj_path <- file.path(OBJS, "gse154977_kp30w_seurat.rds")
saveRDS(obj, obj_path)
message("Wrote ", obj_path, " size=", file.info(obj_path)$size)

sp_ifn <- spearman_safe(primary$mean_Cldn4, primary$mean_IFN)
sp_mhc <- spearman_safe(primary$mean_Cldn4, primary$mean_MHC)
sp_tj <- spearman_safe(primary$mean_Cldn4, primary$mean_TJ)

summary <- list(
  accession = "GSE154977",
  paper = "Marjanovic et al. Cancer Cell 2020 PMID 32707077",
  seurat_version = as.character(packageVersion("Seurat")),
  r_version = as.character(getRversion()),
  n_cells_depositor = 11017,
  n_cells_seurat_qc = ncol(obj),
  n_mice = nrow(primary),
  mice = primary$mouse,
  treatments = as.list(setNames(primary$treatment, primary$mouse)),
  tnk_design_nogo = TRUE,
  tnk_leak_n = unname(tnk_n),
  tnk_leak_frac = unname(tnk_frac),
  epcam_n = unname(sum(obj$epcam_pos)),
  cldn4_n = unname(sum(obj$Cldn4_counts > 0)),
  families_n = lapply(families, length),
  spearman = list(
    IFN = list(n = sp_ifn$n, rho = sp_ifn$rho, p = sp_ifn$p),
    MHC = list(n = sp_mhc$n, rho = sp_mhc$rho, p = sp_mhc$p),
    TJ = list(n = sp_tj$n, rho = sp_tj$rho, p = sp_tj$p)
  ),
  q4q1_usable = FALSE,
  object = obj_path,
  mouse_table = file.path(TABLES, "mouse_units.tsv")
)
write_json(summary, file.path(TABLES, "summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = 6)
write.table(data.frame(gene_set = names(families), n = vapply(families, length, integer(1))),
            file.path(TABLES, "geneset_n.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

message("DONE")
message("mice: ", paste(primary$mouse, collapse = ", "))
message(sprintf("IFN rho=%.3f p=%.3f; MHC rho=%.3f p=%.3f; TJ rho=%.3f p=%.3f",
                sp_ifn$rho, sp_ifn$p, sp_mhc$rho, sp_mhc$p, sp_tj$rho, sp_tj$p))
