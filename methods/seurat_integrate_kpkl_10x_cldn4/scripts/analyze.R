#!/usr/bin/env Rscript
# Integrate public mouse 10x lung GEMM: GSE154977 (KP) + GSE180963 (K/KL) + GSE165641 (KL)
# if each has a processed matrix. R + Seurat / Harmony. Cldn4-only.
# Honest unit = mouse. Dataset is a covariate. Within-genotype + leave-one-dataset-out.
# No GSE179502 / GSE154989 / GSE267321 / GSE127465 / GSE50927 / human / private 8-KL.
# No dual-high. No Python-only primary. Drop any series with no matrix.
#
# Usage:
#   Rscript methods/seurat_integrate_kpkl_10x_cldn4/scripts/analyze.R \
#     --data /tmp/kpkl_10x --out methods/seurat_integrate_kpkl_10x_cldn4

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(hdf5r)
})

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
data_root <- "/tmp/kpkl_10x"
out_dir <- "methods/seurat_integrate_kpkl_10x_cldn4"
i <- 1
while (i <= length(args)) {
  if (args[[i]] == "--data" && i < length(args)) {
    data_root <- args[[i + 1]]
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

# Locked Cldn4-only mouse sets (same family as GSE180963 / GSE165641).
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
MIN_EPI_CELLS <- 10L
MIN_TNK_CELLS <- 20L
MIN_TNK_FRAC <- 0.02

has_matrix <- function(gse) {
  f <- file.path(data_root, gse, "HAS_MATRIX")
  file.exists(f) && identical(trimws(readLines(f, warn = FALSE)[1]), "yes")
}

present_in <- function(genes, universe) genes[genes %in% universe]

assay_mat <- function(obj, layer) {
  tryCatch(
    GetAssayData(obj, assay = "RNA", layer = layer),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = layer)
  )
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

read_h5_col <- function(ds) {
  # MATLAB-style COO dumps as N x 1, not a 1-d vector.
  if (length(ds$dims) >= 2 && ds$dims[2] == 1) {
    as.numeric(ds[, 1])
  } else {
    as.numeric(ds[])
  }
}

read_coo_h5 <- function(path, genes, cells) {
  h5 <- H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  i <- read_h5_col(h5[["i"]])
  j <- read_h5_col(h5[["j"]])
  v <- read_h5_col(h5[["v"]])
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

tag_meta <- function(obj, dataset, mouse, gsm, genotype, treatment, digest, library_id) {
  obj$dataset <- dataset
  obj$mouse <- mouse
  obj$gsm <- gsm
  obj$genotype <- genotype
  obj$treatment <- treatment
  obj$digest <- digest
  obj$library_id <- library_id
  obj$orig.ident <- mouse
  obj
}

load_gse154977 <- function() {
  d <- file.path(data_root, "GSE154977")
  smp <- read.csv(file.path(d, "GSE154977_mmLung10x_cis_smpTable.csv.gz"),
                  stringsAsFactors = FALSE)
  gene_tbl <- read.csv(file.path(d, "GSE154977_mmLung10x_cis_geneTable.csv.gz"),
                       stringsAsFactors = FALSE)
  if (anyDuplicated(gene_tbl$geneID)) {
    gene_tbl$geneID <- make.unique(gene_tbl$geneID)
  }
  genes <- gene_tbl$geneID
  cells <- smp$sampleID
  counts <- read_coo_h5(file.path(d, "GSE154977_mmLung10x_cis_dSp_rawCount.h5"),
                        genes, cells)
  meta <- data.frame(row.names = cells, stringsAsFactors = FALSE)
  meta$library_raw <- sub("_id-.*$", "", rownames(meta))
  meta$mouse_raw <- sub("_PT$", "", meta$library_raw)
  gsm_map <- c(
    KP_30w_ND_m3_PT = "GSM4685281",
    KP_30w_ND_m4_PT = "GSM4685282",
    KP_30w_Cis72_m5_PT = "GSM4685283",
    KP_30w_Cis72_m6_PT = "GSM4685284"
  )
  obj <- CreateSeuratObject(counts = counts, project = "GSE154977",
                            min.cells = 0, min.features = 0)
  obj <- tag_meta(
    obj,
    dataset = "GSE154977",
    mouse = paste0("GSE154977_", meta$mouse_raw),
    gsm = unname(gsm_map[meta$library_raw]),
    genotype = "KP",
    treatment = ifelse(grepl("Cis72", meta$library_raw), "Cis72", "ND"),
    digest = "AT2_lineage_FACS",
    library_id = meta$library_raw
  )
  obj
}

load_gse180963 <- function() {
  d <- file.path(data_root, "GSE180963")
  samples <- data.frame(
    label = c("K", "KL"),
    gsm = c("GSM5481386", "GSM5481387"),
    genotype = c("K", "KL"),
    stringsAsFactors = FALSE
  )
  objs <- lapply(seq_len(nrow(samples)), function(k) {
    lab <- samples$label[k]
    counts <- Read10X(data.dir = file.path(d, lab), gene.column = 2,
                      unique.features = TRUE)
    obj <- CreateSeuratObject(counts = counts, project = "GSE180963",
                              min.cells = 0, min.features = 0)
    tag_meta(
      obj,
      dataset = "GSE180963",
      mouse = paste0("GSE180963_", lab),
      gsm = samples$gsm[k],
      genotype = samples$genotype[k],
      treatment = "untreated",
      digest = "mixed",
      library_id = paste0(samples$gsm[k], "_", lab)
    )
  })
  out <- merge(objs[[1]], y = objs[[2]], add.cell.ids = samples$label,
               project = "GSE180963")
  if (packageVersion("SeuratObject") >= "5.0.0") out <- JoinLayers(out)
  out
}

find_filtered_dir <- function(root, hint) {
  hits <- list.files(root, pattern = "filtered_feature_bc_matrix$",
                     recursive = TRUE, include.dirs = TRUE, full.names = TRUE)
  if (length(hits) == 0) stop("no filtered_feature_bc_matrix under ", root)
  pick <- hits[grepl(hint, hits, ignore.case = TRUE)]
  if (length(pick) == 0) pick <- hits
  pick[[1]]
}

load_gse165641 <- function() {
  d <- file.path(data_root, "GSE165641")
  specs <- data.frame(
    hint = c("KL1", "KL2"),
    mouse = c("GSE165641_KL1", "GSE165641_KL2"),
    gsm = c("GSM5047302", "GSM5047303"),
    stringsAsFactors = FALSE
  )
  objs <- lapply(seq_len(nrow(specs)), function(k) {
    mtx_dir <- find_filtered_dir(d, specs$hint[k])
    counts <- Read10X(mtx_dir)
    if (is.list(counts)) counts <- counts[[1]]
    obj <- CreateSeuratObject(counts = counts, project = "GSE165641",
                              min.cells = 0, min.features = 0)
    tag_meta(
      obj,
      dataset = "GSE165641",
      mouse = specs$mouse[k],
      gsm = specs$gsm[k],
      genotype = "KL",
      treatment = "untreated",
      digest = "mixed",
      library_id = specs$gsm[k]
    )
  })
  out <- merge(objs[[1]], y = objs[[2]], add.cell.ids = specs$hint,
               project = "GSE165641")
  if (packageVersion("SeuratObject") >= "5.0.0") out <- JoinLayers(out)
  out
}

spearman_safe <- function(x, y, min_n = MIN_N_SPEARMAN) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  out <- list(n = n, rho = NA_real_, p = NA_real_, usable = FALSE)
  if (n < min_n) return(out)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  out$rho <- unname(ct$estimate)
  out$p <- unname(ct$p.value)
  out$usable <- TRUE
  out
}

partial_spearman <- function(x, y, cov, min_n = MIN_N_SPEARMAN) {
  ok <- is.finite(x) & is.finite(y) & !is.na(cov)
  n <- sum(ok)
  out <- list(n = n, rho = NA_real_, p = NA_real_, usable = FALSE)
  if (n < min_n) return(out)
  cov <- droplevels(factor(cov[ok]))
  if (nlevels(cov) < 2) return(spearman_safe(x[ok], y[ok], min_n))
  rx <- tryCatch(resid(lm(x[ok] ~ cov)), error = function(e) NULL)
  ry <- tryCatch(resid(lm(y[ok] ~ cov)), error = function(e) NULL)
  if (is.null(rx) || is.null(ry)) return(out)
  spearman_safe(rx, ry, min_n)
}

q4_vs_q1 <- function(x, y, min_n = MIN_N_Q4) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  out <- list(
    n = n, n_q1 = NA_integer_, n_q4 = NA_integer_,
    median_q1 = NA_real_, median_q4 = NA_real_,
    delta_median = NA_real_, p = NA_real_, usable = FALSE
  )
  if (n < min_n) return(out)
  qs <- quantile(x[ok], probs = c(0.25, 0.75), names = FALSE, type = 7)
  q1 <- x[ok] <= qs[1]
  q4 <- x[ok] >= qs[2]
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

fmt <- function(x, d = 3) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  formatC(as.numeric(x), format = "f", digits = d)
}

rho_txt <- function(sp) {
  if (!isTRUE(sp$usable)) {
    return(sprintf("n=%d mice, Spearman locked off (need n>=%d)",
                   sp$n, MIN_N_SPEARMAN))
  }
  sprintf("n=%d mice, ρ=%.3f, p=%.3g", sp$n, sp$rho, sp$p)
}

contrast_row <- function(subset, xname, yname, x, y, cov = NULL, note = "") {
  sp <- spearman_safe(x, y)
  psp <- if (!is.null(cov)) partial_spearman(x, y, cov) else
    list(n = length(x), rho = NA_real_, p = NA_real_, usable = FALSE)
  q <- q4_vs_q1(x, y)
  data.frame(
    subset = subset,
    x = xname,
    y = yname,
    n_mice = length(x),
    spearman_n = sp$n,
    spearman_rho = sp$rho,
    spearman_p = sp$p,
    spearman_usable = sp$usable,
    dataset_partial_n = psp$n,
    dataset_partial_rho = psp$rho,
    dataset_partial_p = psp$p,
    dataset_partial_usable = psp$usable,
    q4q1_usable = q$usable,
    q4q1_delta_median = q$delta_median,
    q4q1_p = q$p,
    note = note,
    stringsAsFactors = FALSE
  )
}

# ---- load series that have a processed matrix ----
dropped <- character()
loaded <- list()
inventory <- list()

if (has_matrix("GSE154977")) {
  message("CreateSeuratObject GSE154977")
  loaded$GSE154977 <- load_gse154977()
  inventory[[length(inventory) + 1]] <- data.frame(
    dataset = "GSE154977", kept = TRUE,
    n_cells_loaded = ncol(loaded$GSE154977),
    n_genes_loaded = nrow(loaded$GSE154977),
    matrix = "GSE154977_mmLung10x_cis_dSp_rawCount.h5",
    stringsAsFactors = FALSE
  )
} else {
  dropped <- c(dropped, "GSE154977")
  inventory[[length(inventory) + 1]] <- data.frame(
    dataset = "GSE154977", kept = FALSE, n_cells_loaded = 0,
    n_genes_loaded = 0, matrix = "MISSING — dropped",
    stringsAsFactors = FALSE
  )
}

if (has_matrix("GSE180963")) {
  message("CreateSeuratObject GSE180963")
  loaded$GSE180963 <- load_gse180963()
  inventory[[length(inventory) + 1]] <- data.frame(
    dataset = "GSE180963", kept = TRUE,
    n_cells_loaded = ncol(loaded$GSE180963),
    n_genes_loaded = nrow(loaded$GSE180963),
    matrix = "GSE180963_RAW.tar K/ + KL/ MTX",
    stringsAsFactors = FALSE
  )
} else {
  dropped <- c(dropped, "GSE180963")
  inventory[[length(inventory) + 1]] <- data.frame(
    dataset = "GSE180963", kept = FALSE, n_cells_loaded = 0,
    n_genes_loaded = 0, matrix = "MISSING — dropped",
    stringsAsFactors = FALSE
  )
}

if (has_matrix("GSE165641")) {
  message("CreateSeuratObject GSE165641")
  loaded$GSE165641 <- load_gse165641()
  inventory[[length(inventory) + 1]] <- data.frame(
    dataset = "GSE165641", kept = TRUE,
    n_cells_loaded = ncol(loaded$GSE165641),
    n_genes_loaded = nrow(loaded$GSE165641),
    matrix = "GSM5047302/03 Cell Ranger filtered MTX",
    stringsAsFactors = FALSE
  )
} else {
  dropped <- c(dropped, "GSE165641")
  inventory[[length(inventory) + 1]] <- data.frame(
    dataset = "GSE165641", kept = FALSE, n_cells_loaded = 0,
    n_genes_loaded = 0, matrix = "MISSING — dropped",
    stringsAsFactors = FALSE
  )
}

inv_df <- do.call(rbind, inventory)
write.table(inv_df, file.path(tab_dir, "dataset_inventory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (length(loaded) < 2) {
  stop("Need at least two series with processed matrices to integrate. Kept: ",
       paste(names(loaded), collapse = ","),
       " Dropped: ", paste(dropped, collapse = ","))
}

# Common gene symbols only (Cldn4 must survive).
gene_lists <- lapply(loaded, rownames)
common <- Reduce(intersect, gene_lists)
message("common genes: ", length(common))
if (!(CLDN4 %in% common)) {
  stop("Cldn4 is not in the intersected gene universe — cannot score Cldn4-only.")
}
loaded <- lapply(loaded, function(o) o[common, ])

# Merge
if (length(loaded) == 2) {
  obj <- merge(loaded[[1]], y = loaded[[2]],
               add.cell.ids = names(loaded), project = "KP_KL_10x")
} else {
  obj <- merge(loaded[[1]], y = loaded[-1],
               add.cell.ids = names(loaded), project = "KP_KL_10x")
}
rm(loaded)
gc()
if (packageVersion("SeuratObject") >= "5.0.0") {
  obj <- JoinLayers(obj)
}

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
n_before <- ncol(obj)
# Liberal common QC. Author matrices are already cell-called.
obj <- subset(obj, subset = nFeature_RNA >= 200 & nCount_RNA >= 500 & percent.mt < 25)
n_after <- ncol(obj)
message("QC ", n_after, " / ", n_before)

obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, verbose = FALSE)

# Harmony: dataset is the batch. Do not Harmony-out mouse (that is the unit).
harmony_method <- "none"
if (requireNamespace("harmony", quietly = TRUE)) {
  message("RunHarmony by dataset")
  obj <- tryCatch({
    o <- harmony::RunHarmony(
      obj,
      group.by.vars = "dataset",
      reduction.use = "pca",
      dims.use = 1:30,
      reduction.save = "harmony",
      verbose = TRUE
    )
    harmony_method <<- "harmony::RunHarmony"
    o
  }, error = function(e) {
    message("RunHarmony (reduction.use) failed: ", e$message)
    tryCatch({
      o <- harmony::RunHarmony(
        obj,
        group.by.vars = "dataset",
        reduction = "pca",
        dims.use = 1:30,
        reduction.save = "harmony",
        verbose = TRUE
      )
      harmony_method <<- "harmony::RunHarmony(reduction=pca)"
      o
    }, error = function(e2) {
      message("RunHarmony (reduction=) failed: ", e2$message)
      obj
    })
  })
}
if (!("harmony" %in% Reductions(obj))) {
  if ("HarmonyIntegration" %in% getNamespaceExports("Seurat")) {
    message("IntegrateLayers HarmonyIntegration")
    obj <- IntegrateLayers(
      obj,
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "harmony",
      verbose = TRUE
    )
    harmony_method <- "Seurat::IntegrateLayers(HarmonyIntegration)"
  } else {
    stop("Harmony reduction was not created. Install r-cran-harmony.")
  }
}

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20, verbose = FALSE)

universe <- rownames(obj)
counts <- assay_mat(obj, "counts")
logn <- assay_mat(obj, "data")

inv_genes <- unique(c(CLDN4, T_NK_SCORE, IFN, MHC, EPI_CORE, HOST_LUNG, AUDIT))
gene_inv <- data.frame(
  gene = inv_genes,
  present = inv_genes %in% universe,
  stringsAsFactors = FALSE
)
write.table(gene_inv, file.path(tab_dir, "gene_inventory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Marker compartments. Cldn4 is never a caller.
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
# GSE154977 is AT2-lineage FACS / tumor-state 10x: Epcam+ (or structural+) and Ptprc-
# is accepted as epithelium when the digest is not mixed.
epi_facs <- (obj$digest == "AT2_lineage_FACS") & (ptprc == 0) &
  ((epcam > 0) | (struct > 0))
is_epi <- epi_tight | epi_facs
tnk_mark <- pos_any(counts, T_NK_CALL)
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

md <- leo.a@example.org
mice <- sort(unique(md$mouse))

mouse_rows <- lapply(mice, function(lab) {
  w <- md$mouse == lab
  we <- w & md$is_epi
  wn <- w & md$is_tnk
  n_cells <- sum(w)
  n_epi <- sum(we)
  n_tnk <- sum(wn)
  frac_tnk <- if (n_cells > 0) n_tnk / n_cells else NA_real_
  digest <- unique(as.character(md$digest[w]))[1]
  tnk_usable <- (digest == "mixed") && (n_tnk >= MIN_TNK_CELLS) &&
    is.finite(frac_tnk) && (frac_tnk >= MIN_TNK_FRAC)
  epi_usable <- n_epi >= MIN_EPI_CELLS
  data.frame(
    mouse = lab,
    dataset = unique(as.character(md$dataset[w]))[1],
    gsm = unique(as.character(md$gsm[w]))[1],
    genotype = unique(as.character(md$genotype[w]))[1],
    treatment = unique(as.character(md$treatment[w]))[1],
    digest = digest,
    n_cells = n_cells,
    n_epi = n_epi,
    n_tnk = n_tnk,
    frac_tnk = frac_tnk,
    frac_epi = if (n_cells > 0) n_epi / n_cells else NA_real_,
    frac_sftpc = mean(md$sftpc_pos[w]),
    frac_ptprc = mean(md$ptprc_pos[w]),
    Cldn4_all_mean = mean(md$Cldn4_logn[w]),
    Cldn4_all_pctpos = 100 * mean(md$Cldn4_pos[w]),
    Cldn4_epi_mean = if (any(we)) mean(md$Cldn4_logn[we]) else NA_real_,
    Cldn4_epi_pctpos = if (any(we)) 100 * mean(md$Cldn4_pos[we]) else NA_real_,
    IFN_epi_mean = if (any(we)) mean(md$ifn_score[we]) else NA_real_,
    MHC_epi_mean = if (any(we)) mean(md$mhc_score[we]) else NA_real_,
    tnk_usable = tnk_usable,
    epi_usable = epi_usable,
    tnk_exclude_reason = if (digest != "mixed") {
      "not mixed digest (GSE154977 AT2-lineage / tumor-state FACS; T/NK is design no-go)"
    } else if (!tnk_usable) {
      "mixed digest but T/NK below floor"
    } else {
      ""
    },
    stringsAsFactors = FALSE
  )
})
per_mouse <- do.call(rbind, mouse_rows)
per_mouse <- per_mouse[order(per_mouse$dataset, per_mouse$genotype, per_mouse$mouse), ]
write.table(per_mouse, file.path(tab_dir, "mouse_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

comp <- as.data.frame(table(mouse = md$mouse, compartment = md$compartment),
                      stringsAsFactors = FALSE)
write.table(comp, file.path(tab_dir, "compartment_by_mouse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- mouse-level tests ----
tnk_df <- per_mouse[per_mouse$tnk_usable, , drop = FALSE]
epi_df <- per_mouse[per_mouse$epi_usable, , drop = FALSE]
epi_nd <- epi_df[epi_df$treatment != "Cis72", , drop = FALSE]

contrasts <- list()

add_contrasts <- function(df, subset, note) {
  if (nrow(df) == 0) {
    contrasts[[length(contrasts) + 1]] <<- contrast_row(
      subset, "Cldn4_epi_mean", "frac_tnk",
      numeric(0), numeric(0), NULL, paste(note, "empty")
    )
    return(invisible(NULL))
  }
  if ("frac_tnk" %in% names(df) && any(df$tnk_usable %in% TRUE)) {
    d2 <- df[df$tnk_usable, , drop = FALSE]
    contrasts[[length(contrasts) + 1]] <<- contrast_row(
      subset, "Cldn4_epi_mean", "frac_tnk",
      d2$Cldn4_epi_mean, d2$frac_tnk, d2$dataset, note
    )
  }
  contrasts[[length(contrasts) + 1]] <<- contrast_row(
    subset, "Cldn4_epi_mean", "IFN_epi_mean",
    df$Cldn4_epi_mean, df$IFN_epi_mean, df$dataset, note
  )
  contrasts[[length(contrasts) + 1]] <<- contrast_row(
    subset, "Cldn4_epi_mean", "MHC_epi_mean",
    df$Cldn4_epi_mean, df$MHC_epi_mean, df$dataset, note
  )
}

add_contrasts(tnk_df, "TNK_all_usable",
              "mixed-digest mice only; dataset as covariate")
add_contrasts(epi_df, "EPI_all_usable",
              "mice with epithelium; dataset as covariate")
add_contrasts(epi_nd, "EPI_no_Cis72",
              "drop GSE154977 Cis72; dataset as covariate")

for (g in sort(unique(per_mouse$genotype))) {
  add_contrasts(tnk_df[tnk_df$genotype == g, , drop = FALSE],
                paste0("TNK_within_", g),
                paste0("within-genotype ", g, " only"))
  add_contrasts(epi_df[epi_df$genotype == g, , drop = FALSE],
                paste0("EPI_within_", g),
                paste0("within-genotype ", g, " only"))
}

for (ds in sort(unique(per_mouse$dataset))) {
  add_contrasts(tnk_df[tnk_df$dataset != ds, , drop = FALSE],
                paste0("TNK_LODO_drop_", ds),
                paste0("leave-one-dataset-out: drop ", ds))
  add_contrasts(epi_df[epi_df$dataset != ds, , drop = FALSE],
                paste0("EPI_LODO_drop_", ds),
                paste0("leave-one-dataset-out: drop ", ds))
}

contrast_tab <- do.call(rbind, contrasts)
write.table(contrast_tab, file.path(tab_dir, "mouse_contrasts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Primary rows
sp_tnk <- spearman_safe(tnk_df$Cldn4_epi_mean, tnk_df$frac_tnk)
psp_tnk <- partial_spearman(tnk_df$Cldn4_epi_mean, tnk_df$frac_tnk, tnk_df$dataset)
sp_ifn <- spearman_safe(epi_df$Cldn4_epi_mean, epi_df$IFN_epi_mean)
psp_ifn <- partial_spearman(epi_df$Cldn4_epi_mean, epi_df$IFN_epi_mean, epi_df$dataset)
sp_mhc <- spearman_safe(epi_df$Cldn4_epi_mean, epi_df$MHC_epi_mean)
psp_mhc <- partial_spearman(epi_df$Cldn4_epi_mean, epi_df$MHC_epi_mean, epi_df$dataset)

honest <- data.frame(
  item = c(
    "GEO series requested",
    "GEO series kept (processed matrix)",
    "GEO series dropped (no matrix)",
    "mice (unit)",
    "KP mice (GSE154977)",
    "K mice (GSE180963)",
    "KL mice (GSE180963 + GSE165641)",
    "T/NK-usable mice",
    "epithelium-usable mice",
    "cells after QC (not the unit)",
    "genes (intersected symbols)",
    "Cldn4-positive cells",
    "dual-high TACSTD2 x Cldn4",
    "private 8-KL mice used",
    "human series",
    "excluded GEO (locked)"
  ),
  n = c(
    3,
    length(unique(per_mouse$dataset)),
    length(dropped),
    nrow(per_mouse),
    sum(per_mouse$genotype == "KP"),
    sum(per_mouse$genotype == "K"),
    sum(per_mouse$genotype == "KL"),
    nrow(tnk_df),
    nrow(epi_df),
    n_after,
    nrow(obj),
    sum(md$Cldn4_pos),
    0,
    0,
    0,
    6
  ),
  note = c(
    "GSE154977 + GSE180963 + GSE165641",
    paste(sort(unique(per_mouse$dataset)), collapse = ","),
    if (length(dropped)) paste(dropped, collapse = ",") else "none",
    "honest unit; do not write cell n",
    "AT2-lineage / tumor-state FACS; T/NK design no-go",
    "mixed digest",
    "mixed digest",
    "mixed digest AND T/NK floor",
    sprintf("n_epi >= %d", MIN_EPI_CELLS),
    sprintf("nFeature>=200, nCount>=500, mt<25; before=%d", n_before),
    "symbol intersect across kept series",
    "count > 0",
    "not defined",
    "public GEO only",
    "none",
    "GSE179502, GSE154989, GSE267321, GSE127465, GSE50927, human"
  ),
  stringsAsFactors = FALSE
)
write.table(honest, file.path(tab_dir, "honest_n.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- figures ----
theme_set(theme_bw(base_size = 11))
p_ds <- DimPlot(obj, group.by = "dataset", reduction = "umap") +
  ggtitle("Harmony UMAP — dataset (batch)")
ggsave(file.path(fig_dir, "DimPlot_dataset.png"), p_ds, width = 7.2, height = 5.4, dpi = 140)
ggsave(file.path(fig_dir, "DimPlot_dataset.pdf"), p_ds, width = 7.2, height = 5.4)

p_gt <- DimPlot(obj, group.by = "genotype", reduction = "umap") +
  ggtitle("Harmony UMAP — genotype (K / KP / KL)")
ggsave(file.path(fig_dir, "DimPlot_genotype.png"), p_gt, width = 7.2, height = 5.4, dpi = 140)
ggsave(file.path(fig_dir, "DimPlot_genotype.pdf"), p_gt, width = 7.2, height = 5.4)

p_cp <- DimPlot(obj, group.by = "compartment", reduction = "umap") +
  ggtitle("Harmony UMAP — marker compartment (Cldn4 not a caller)")
ggsave(file.path(fig_dir, "DimPlot_compartment.png"), p_cp, width = 7.2, height = 5.4, dpi = 140)
ggsave(file.path(fig_dir, "DimPlot_compartment.pdf"), p_cp, width = 7.2, height = 5.4)

p_ms <- DimPlot(obj, group.by = "mouse", reduction = "umap") +
  ggtitle("Harmony UMAP — mouse (the unit)")
ggsave(file.path(fig_dir, "DimPlot_mouse.png"), p_ms, width = 8.0, height = 5.6, dpi = 140)
ggsave(file.path(fig_dir, "DimPlot_mouse.pdf"), p_ms, width = 8.0, height = 5.6)

p_c4 <- FeaturePlot(obj, features = "Cldn4", reduction = "umap", order = TRUE) +
  ggtitle("Harmony UMAP — Cldn4 (lognorm)")
ggsave(file.path(fig_dir, "FeaturePlot_Cldn4.png"), p_c4, width = 6.6, height = 5.2, dpi = 140)

if (nrow(tnk_df) > 0) {
  p_xy <- ggplot(tnk_df, aes(x = Cldn4_epi_mean, y = frac_tnk,
                             color = genotype, shape = dataset, label = mouse)) +
    geom_point(size = 3.4) +
    geom_text(vjust = -0.7, size = 2.6, show.legend = FALSE) +
    labs(title = "Mouse-level Cldn4 (epithelium) vs T/NK fraction",
         subtitle = paste0("Usable mixed-digest mice only. ", rho_txt(sp_tnk),
                           ". Dataset-adjusted: ", rho_txt(psp_tnk),
                           ". Do not sell a genotype mix as a Cldn4 effect."),
         x = "Cldn4 mean lognorm in marker epithelium",
         y = "T/NK fraction")
  ggsave(file.path(fig_dir, "fig_cldn4_vs_tnk.png"), p_xy, width = 7.4, height = 5.2, dpi = 140)
  ggsave(file.path(fig_dir, "fig_cldn4_vs_tnk.pdf"), p_xy, width = 7.4, height = 5.2)
}

if (nrow(epi_df) > 0) {
  p_ifn <- ggplot(epi_df, aes(x = Cldn4_epi_mean, y = IFN_epi_mean,
                              color = genotype, shape = dataset, label = mouse)) +
    geom_point(size = 3.4) +
    geom_text(vjust = -0.7, size = 2.6, show.legend = FALSE) +
    labs(title = "Mouse-level epithelial Cldn4 vs IFN ISG",
         subtitle = paste0(rho_txt(sp_ifn), ". Dataset-adjusted: ", rho_txt(psp_ifn)),
         x = "Cldn4 mean lognorm in marker epithelium",
         y = "IFN ISG mean lognorm in epithelium")
  ggsave(file.path(fig_dir, "fig_cldn4_vs_ifn.png"), p_ifn, width = 7.4, height = 5.2, dpi = 140)
  ggsave(file.path(fig_dir, "fig_cldn4_vs_ifn.pdf"), p_ifn, width = 7.4, height = 5.2)

  p_mhc <- ggplot(epi_df, aes(x = Cldn4_epi_mean, y = MHC_epi_mean,
                              color = genotype, shape = dataset, label = mouse)) +
    geom_point(size = 3.4) +
    geom_text(vjust = -0.7, size = 2.6, show.legend = FALSE) +
    labs(title = "Mouse-level epithelial Cldn4 vs MHC/APM",
         subtitle = paste0(rho_txt(sp_mhc), ". Dataset-adjusted: ", rho_txt(psp_mhc)),
         x = "Cldn4 mean lognorm in marker epithelium",
         y = "MHC/APM mean lognorm in epithelium")
  ggsave(file.path(fig_dir, "fig_cldn4_vs_mhc.png"), p_mhc, width = 7.4, height = 5.2, dpi = 140)
  ggsave(file.path(fig_dir, "fig_cldn4_vs_mhc.pdf"), p_mhc, width = 7.4, height = 5.2)
}

# ---- save object ----
obj_path <- file.path(obj_dir, "kpkl_10x_harmony_seurat.rds")
message("Saving ", obj_path)
saveRDS(obj, obj_path)
message("object bytes=", file.info(obj_path)$size)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(
    list(
      engine = "R_Seurat_Harmony",
      seurat_version = as.character(packageVersion("Seurat")),
      harmony_version = if (requireNamespace("harmony", quietly = TRUE))
        as.character(packageVersion("harmony")) else NA,
      harmony_method = harmony_method,
      python_primary = FALSE,
      private_8_kl = FALSE,
      cldn4_only = TRUE,
      dual_high = FALSE,
      unit = "mouse",
      datasets_kept = sort(unique(per_mouse$dataset)),
      datasets_dropped = dropped,
      n_mice = nrow(per_mouse),
      n_tnk_usable = nrow(tnk_df),
      n_epi_usable = nrow(epi_df),
      n_cells = n_after,
      n_genes = nrow(obj),
      object = obj_path
    ),
    file.path(tab_dir, "summary.json"),
    auto_unbox = TRUE, pretty = TRUE, digits = 8
  )
}

# ---- FINDING.md ----
md_mouse <- function() {
  hdr <- paste(
    "| mouse | dataset | GSM | genotype | treatment | digest | n cells | n epi | n T/NK | frac T/NK | Cldn4 epi | Cldn4 epi %pos | IFN epi | MHC epi | T/NK usable | epi usable |",
    "|---|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    sep = "\n"
  )
  rows <- apply(per_mouse, 1, function(r) {
    sprintf(
      "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
      r[["mouse"]], r[["dataset"]], r[["gsm"]], r[["genotype"]], r[["treatment"]],
      r[["digest"]], r[["n_cells"]], r[["n_epi"]], r[["n_tnk"]],
      fmt(as.numeric(r[["frac_tnk"]])),
      fmt(as.numeric(r[["Cldn4_epi_mean"]]), 4),
      fmt(as.numeric(r[["Cldn4_epi_pctpos"]]), 2),
      fmt(as.numeric(r[["IFN_epi_mean"]]), 3),
      fmt(as.numeric(r[["MHC_epi_mean"]]), 3),
      r[["tnk_usable"]], r[["epi_usable"]]
    )
  })
  paste(c(hdr, rows), collapse = "\n")
}

md_contrasts <- function(prefix) {
  sub <- contrast_tab[startsWith(contrast_tab$subset, prefix) |
                        contrast_tab$subset == prefix, , drop = FALSE]
  if (nrow(sub) == 0) return("_none_\n")
  hdr <- paste(
    "| subset | x | y | n mice | Spearman | dataset-adjusted | Q4 vs Q1 | note |",
    "|---|---|---|---:|---|---|---|---|",
    sep = "\n"
  )
  rows <- apply(sub, 1, function(r) {
    sp <- if (identical(as.character(r[["spearman_usable"]]), "TRUE") ||
              isTRUE(as.logical(r[["spearman_usable"]]))) {
      sprintf("ρ=%s p=%s", fmt(as.numeric(r[["spearman_rho"]])),
              formatC(as.numeric(r[["spearman_p"]]), format = "g", digits = 3))
    } else {
      sprintf("locked off (n=%s)", r[["spearman_n"]])
    }
    psp <- if (identical(as.character(r[["dataset_partial_usable"]]), "TRUE") ||
               isTRUE(as.logical(r[["dataset_partial_usable"]]))) {
      sprintf("ρ=%s p=%s", fmt(as.numeric(r[["dataset_partial_rho"]])),
              formatC(as.numeric(r[["dataset_partial_p"]]), format = "g", digits = 3))
    } else {
      "locked off"
    }
    q <- if (identical(as.character(r[["q4q1_usable"]]), "TRUE") ||
             isTRUE(as.logical(r[["q4q1_usable"]]))) {
      sprintf("Δmed=%s p=%s", fmt(as.numeric(r[["q4q1_delta_median"]])),
              formatC(as.numeric(r[["q4q1_p"]]), format = "g", digits = 3))
    } else {
      "locked off (need n>=8)"
    }
    sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
            r[["subset"]], r[["x"]], r[["y"]], r[["n_mice"]], sp, psp, q, r[["note"]])
  })
  paste(c(hdr, rows), collapse = "\n")
}

kept <- sort(unique(per_mouse$dataset))
finding <- paste0(
  "# FINDING — Seurat/Harmony integrate KP+K+KL 10x (Cldn4-only)\n\n",
  "**ADDITIVE. Public mouse. Cldn4-only. No dual-high.** Thesis is already correct and is not rewritten. ",
  "Primary engine is **R + Seurat ", as.character(packageVersion("Seurat")),
  " + Harmony** (`CreateSeuratObject` on public processed matrices, then `RunHarmony` by **dataset**). ",
  "No Python-only primary. Private 8-KL matrices were not opened. ",
  "Human series were not opened. Locked-out GEO (GSE179502 FACS epi / no T/NK, GSE154989 plate epi, ",
  "GSE267321 subq Cldn4 floor, GSE127465 CD45+, GSE50927) were not added.\n\n",
  "**Requested trio.** GSE154977 (KP 30w 10x) + GSE180963 (K/KL) + GSE165641 (KL). ",
  "Each series is kept **only** if a processed matrix exists. ",
  "Dropped for no matrix: ",
  if (length(dropped)) paste(dropped, collapse = ", ") else "**none**",
  ". Kept: **", paste(kept, collapse = ", "), "**.\n\n",
  "**Honest unit = mouse**, not cell n. Dataset is a covariate. ",
  "Within-genotype and leave-one-dataset-out are reported so a KP/K/KL mix is not sold as a Cldn4 effect.\n\n",
  "---\n\n",
  "## Decision\n\n",
  "| Question | Answer |\n",
  "|---|---|\n",
  "| Public processed matrices | GSE154977 **h5 COO**; GSE180963 **10x MTX**; GSE165641 **Cell Ranger MTX** |\n",
  "| CreateSeuratObject | **yes** — Seurat ", as.character(packageVersion("Seurat")), " |\n",
  "| Harmony | **yes** — ", harmony_method, ", `group.by.vars = dataset` |\n",
  "| Integrated object | `objects/kpkl_10x_harmony_seurat.rds` |\n",
  "| Mice (unit) | **", nrow(per_mouse), "** |\n",
  "| T/NK-usable mice | **", nrow(tnk_df), "** (mixed digest only; GSE154977 is AT2-lineage FACS / tumor-state) |\n",
  "| Epithelium-usable mice | **", nrow(epi_df), "** |\n",
  "| Cldn4 vs T/NK (mouse) | ", rho_txt(sp_tnk), " |\n",
  "| Cldn4 vs T/NK, dataset-adjusted | ", rho_txt(psp_tnk), " |\n",
  "| Cldn4 vs epithelial IFN | ", rho_txt(sp_ifn), " |\n",
  "| Cldn4 vs epithelial IFN, dataset-adjusted | ", rho_txt(psp_ifn), " |\n",
  "| Cldn4 vs epithelial MHC | ", rho_txt(sp_mhc), " |\n",
  "| Cldn4 vs epithelial MHC, dataset-adjusted | ", rho_txt(psp_mhc), " |\n",
  "| Dual-high TACSTD2 × Cldn4 | **not defined** |\n",
  "| Private 8 KL | **not used** |\n",
  "| Unit | **mouse** |\n\n",
  "Do not write n = ", n_after, " cells.\n\n",
  "---\n\n",
  "## Honest n\n\n",
  "| item | n | note |\n",
  "|---|---:|---|\n"
)

for (k in seq_len(nrow(honest))) {
  finding <- paste0(
    finding,
    "| ", honest$item[k], " | **", honest$n[k], "** | ", honest$note[k], " |\n"
  )
}

finding <- paste0(
  finding,
  "\n---\n\n",
  "## Per-mouse table (the actual unit)\n\n",
  md_mouse(), "\n\n",
  "T/NK fraction is a real fraction only in **mixed digest** (GSE180963, GSE165641). ",
  "GSE154977 GEO source is *Lung tissue, AT2 cells* with a FACS step (Marjanovic / Cancer Cell 2020 cisplatin 10x). ",
  "Those KP mice are epithelium for IFN/MHC and a **design no-go** for T/NK fraction. ",
  "Cis72 (m5, m6) is a 72 h cisplatin arm — kept in the object, dropped in the `EPI_no_Cis72` sensitivity.\n\n",
  "---\n\n",
  "## Cldn4 vs T/NK (mouse unit)\n\n",
  "Primary: epithelial Cldn4 mean vs T/NK fraction on **tnk_usable** mice. ",
  "Spearman requires n≥4. Q4 vs Q1 requires n≥8.\n\n",
  "**Pooled mixed-digest:** ", rho_txt(sp_tnk), ".\n\n",
  "**Dataset as covariate** (Spearman of residuals after `~ dataset`): ",
  rho_txt(psp_tnk), ".\n\n",
  "Pooled K+KL is a genotype mix. Do not sell it as a Cldn4 law. ",
  "Within-genotype and leave-one-dataset-out:\n\n",
  md_contrasts("TNK_"), "\n\n",
  "---\n\n",
  "## Cldn4 vs epithelial IFN / MHC (mouse unit)\n\n",
  "Primary: epithelial Cldn4 mean vs mean lognorm of the locked IFN or MHC sets, ",
  "mice with n_epi ≥ ", MIN_EPI_CELLS, ".\n\n",
  "**Pooled epithelium:** IFN ", rho_txt(sp_ifn), "; MHC ", rho_txt(sp_mhc), ".\n\n",
  "**Dataset as covariate:** IFN ", rho_txt(psp_ifn), "; MHC ", rho_txt(psp_mhc), ".\n\n",
  "Pooled KP+K+KL is a genotype mix and a dataset mix. Within-genotype and LODO:\n\n",
  md_contrasts("EPI_"), "\n\n",
  "---\n\n",
  "## Epithelium and T/NK (honest)\n\n",
  "- **Tight epithelium (mixed digest)** = Epcam+ **and** (Cdh1 or Krt8/18/19 or Cldn18)+ **and** Ptprc−.\n",
  "- **GSE154977 epithelium** = Ptprc− and (Epcam+ or structural+); AT2-lineage / tumor-state FACS, not a mixed digest.\n",
  "- **T/NK** = Cd3d / Cd3e / Cd3g / Cd8a / Nkg7 / Ncr1 / Klrb1c. Epithelium wins if both fire.\n",
  "- **Cldn4 is not a caller.** Tacstd2 is inventory only. No dual-high gate.\n",
  "- **T/NK usable** = mixed digest AND n_T/NK ≥ ", MIN_TNK_CELLS,
  " AND frac ≥ ", MIN_TNK_FRAC, ".\n\n",
  "Lineage genes present: ",
  paste(present_in(c(EPI_CORE, HOST_LUNG, T_NK_CALL, CLDN4, AUDIT), universe), collapse = ", "),
  ".\nMissing: ",
  {
    miss <- setdiff(c(EPI_CORE, HOST_LUNG, T_NK_CALL, CLDN4, AUDIT), universe)
    if (length(miss) == 0) "none of the core set" else paste(miss, collapse = ", ")
  },
  ".\n\n",
  "---\n\n",
  "## Methods (short)\n\n",
  "1. Public only. Processed matrices from NCBI GEO FTP. No SRA / FASTQ. No private 8-KL object.\n",
  "2. GSE154977: author COO `rawCount.h5` + smp/gene tables. GSE180963: `Read10X` on K/ and KL/ MTX. ",
  "GSE165641: `Read10X` on KL1/KL2 `filtered_feature_bc_matrix`. ",
  "A series with no matrix is dropped.\n",
  "3. Gene universe = **symbol intersect**. Cldn4 must be present.\n",
  "4. `CreateSeuratObject` → light QC (nFeature≥200, nCount≥500, percent.mt<25) → ",
  "LogNormalize 1e4 → VST 2000 → ScaleData on HVG → PCA 30.\n",
  "5. **Harmony by dataset** (not by mouse). Neighbors / clusters / UMAP on Harmony dims 1:20.\n",
  "6. Cldn4-only. Positive = raw count > 0. Module scores = mean lognorm of present locked genes.\n",
  "7. Unit = mouse. Spearman n≥4. Q4 vs Q1 n≥8. Dataset-adjusted = Spearman of `lm(~ dataset)` residuals. ",
  "Within-genotype and leave-one-dataset-out are written even when locked off.\n",
  "8. Thesis is not rewritten.\n\n",
  "```bash\n",
  "bash methods/seurat_integrate_kpkl_10x_cldn4/scripts/install_r.sh\n",
  "bash methods/seurat_integrate_kpkl_10x_cldn4/scripts/download.sh /tmp/kpkl_10x\n",
  "Rscript methods/seurat_integrate_kpkl_10x_cldn4/scripts/analyze.R --data /tmp/kpkl_10x --out methods/seurat_integrate_kpkl_10x_cldn4\n",
  "```\n\n",
  "---\n\n",
  "## How to read this\n\n",
  "- **Additive public mouse**, not a human concordant-pool join.\n",
  "- **GSE154977 is KP epithelium**, not a T/NK fraction series.\n",
  "- **GSE180963 is 1 K + 1 KL.** **GSE165641 is 2 KL.** Together they give mixed-digest T/NK mice.\n",
  "- **A pooled Spearman that mixes genotypes is not a Cldn4 effect.** Read the within-genotype and LODO rows.\n",
  "- **Honest n is mice.** Cell counts are inventory.\n",
  "- **No dual-high. No private 8 KL. No locked-out series. Thesis unchanged.**\n\n",
  "## Files\n\n",
  "- `objects/kpkl_10x_harmony_seurat.rds` — integrated Seurat (gitignored if huge; exists after the run)\n",
  "- `tables/mouse_table.tsv` — unit-level table\n",
  "- `tables/mouse_contrasts.tsv` — pooled / dataset-adjusted / within-genotype / LODO\n",
  "- `tables/honest_n.tsv`, `tables/dataset_inventory.tsv`, `tables/gene_inventory.tsv`, `tables/summary.json`\n",
  "- `figures/DimPlot_dataset.png` (and genotype / compartment / mouse)\n",
  "- `figures/fig_cldn4_vs_tnk.png`, `fig_cldn4_vs_ifn.png`, `fig_cldn4_vs_mhc.png`\n",
  "- `scripts/analyze.R` — R + Seurat / Harmony primary\n\n",
  "## 结论\n\n",
  "用 **R + Seurat / Harmony** 整合了公开小鼠肺 GEMM 10x：GSE154977（KP）、GSE180963（K/KL）、GSE165641（KL），",
  "缺矩阵的系列已丢。诚实单位是 **鼠**，不是细胞数。数据集作为协变量；同时报告基因型内和 leave-one-dataset-out，",
  "避免把基因型混合物写成 Cldn4 效应。Cldn4-only，无 dual-high，无私有 8 只 KL，不改 thesis。\n"
)

writeLines(finding, file.path(out_dir, "FINDING.md"))
message("wrote FINDING.md")
message("n_mice=", nrow(per_mouse), " n_tnk_usable=", nrow(tnk_df),
        " n_epi_usable=", nrow(epi_df), " n_cells=", n_after)
message("T/NK ", rho_txt(sp_tnk), " | IFN ", rho_txt(sp_ifn), " | MHC ", rho_txt(sp_mhc))
