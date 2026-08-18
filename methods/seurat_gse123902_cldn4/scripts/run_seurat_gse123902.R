#!/usr/bin/env Rscript
# ADDITIVE. CLDN4-only. Human LUAD/NSCLC GSE123902.
# Primary engine: R + Seurat. CreateSeuratObject from public MTX written
# from GEO processed dense UMI CSVs. Honest unit = patient/donor (LX ID).
# Mouse-style: malignant CLDN4 %pos / mean vs T/NK fraction;
# malignant Q4 vs Q1 IFN / MHC-I / TJ (CLDN4 held out) if n allows.
# No dual-high. No GSE148071 merge. No Python-only primary.

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is not installed. Stop. Do not fall back to a Python-only primary.")
  }
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)[1]) a else b

args <- commandArgs(trailingOnly = TRUE)
parse_opt <- function(flag, default) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", flag, "="), "", hit[[1]]) else default
}

data_dir <- parse_opt("--data", "/tmp/gse123902_seurat")
here <- tryCatch(
  {
    ca <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", ca[grep("^--file=", ca)])
    normalizePath(file.path(dirname(f), ".."))
  },
  error = function(e) getwd()
)
out_dir <- parse_opt("--out", here)
dir.create(file.path(out_dir, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results", "mtx"), recursive = TRUE, showWarnings = FALSE)

csv_dir <- file.path(data_dir, "csv")
if (!length(list.files(csv_dir, pattern = "GSM.*_dense\\.csv\\.gz$"))) {
  tar_path <- file.path(data_dir, "GSE123902_RAW.tar")
  if (!file.exists(tar_path)) {
    stop("Missing ", tar_path, ". Run scripts/download.sh first.")
  }
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  untar(tar_path, exdir = csv_dir)
}

# ---------------------------------------------------------------------------
# Locked sets. CLDN4 is the readout. It is never used to assign lineage.
# ---------------------------------------------------------------------------
EPI_MARKERS <- c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7")
TNK_MARKERS <- c("CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "KLRD1")
MYE_MARKERS <- c("LYZ", "CD14", "CSF1R", "AIF1")
B_MARKERS <- c("MS4A1", "CD79A")

# Hallmark IFNα ∩ IFNγ core used on this accession (CLDN4 not in set).
IFN_GENES <- c(
  "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "ISG15", "ISG20",
  "MX1", "MX2", "IFIT1", "IFIT2", "IFIT3", "OAS1", "OAS2", "OASL",
  "IFI27", "IFI44", "IFI44L", "IFI6", "RSAD2", "CXCL9", "CXCL10",
  "CXCL11", "B2M", "TAP1", "PSMB8", "PSMB9"
)
MHC_GENES <- c(
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M",
  "TAP1", "TAP2", "TAPBP", "TAPBPL", "NLRC5", "PSMB8", "PSMB9",
  "PSMB10", "ERAP1", "ERAP2", "CALR", "CANX", "PDIA3", "IRF1"
)
# Epithelial tight-junction / barrier genes. CLDN4 held out on purpose.
TJ_GENES <- c(
  "CLDN1", "CLDN3", "CLDN7", "CLDN18", "TJP1", "TJP2", "TJP3",
  "OCLN", "F11R", "CGN", "MARVELD2", "MARVELD3", "CRB3", "CDH1"
)

MIN_MAL <- 20L
MIN_TNK <- 20L
MIN_N_SPEARMAN <- 5L
MIN_N_Q4Q1 <- 8L
SEED <- 1L

parse_fname <- function(fname) {
  m <- regexec(
    "^(GSM[0-9]+)_(MSK_LX[^_]+(?:B)?)_(PRIMARY_TUMOUR|METASTASIS|NORMAL)_dense\\.csv\\.gz$",
    fname
  )
  g <- regmatches(fname, m)[[1]]
  if (length(g) != 4) stop("unparsed filename: ", fname)
  list(gsm = g[2], patient = g[3], site = g[4], file = fname)
}

present <- function(genes, universe) intersect(genes, universe)

marker_score_vec <- function(data_mat, markers) {
  g <- present(markers, rownames(data_mat))
  if (!length(g)) return(rep(0, ncol(data_mat)))
  if (length(g) == 1L) return(as.numeric(data_mat[g, ]))
  as.numeric(Matrix::colMeans(data_mat[g, , drop = FALSE]))
}

assign_lineage <- function(data_mat) {
  scores <- rbind(
    epithelial = marker_score_vec(data_mat, EPI_MARKERS),
    tnk = marker_score_vec(data_mat, TNK_MARKERS),
    myeloid = marker_score_vec(data_mat, MYE_MARKERS),
    b = marker_score_vec(data_mat, B_MARKERS)
  )
  winner <- max.col(t(scores), ties.method = "first")
  top <- apply(scores, 2, max)
  second <- apply(scores, 2, function(x) sort(x, decreasing = TRUE)[2])
  keep <- (top >= 0.12) & (top >= second * 1.15)
  lab <- rep("other", ncol(scores))
  lab[keep] <- rownames(scores)[winner[keep]]
  lab
}

spearman_ci <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < MIN_N_SPEARMAN) {
    return(list(n = n, rho = NA_real_, p = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  z <- atanh(max(min(rho, 0.999999), -0.999999))
  se <- 1 / sqrt(n - 3)
  list(
    n = n, rho = as.numeric(rho), p = as.numeric(ct$p.value),
    lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se)
  )
}

fmt_rho <- function(s) {
  if (is.na(s$rho)) return(sprintf("n=%d (below Spearman floor)", s$n))
  sprintf("n=%d, ρ=%+.3f [%.3f, %.3f], p=%.3g", s$n, s$rho, s$lo, s$hi, s$p)
}

read_geo_csv <- function(path, gsm) {
  dt <- data.table::fread(path, sep = ",", header = TRUE, data.table = FALSE, showProgress = FALSE)
  bc <- as.character(dt[[1]])
  genes <- colnames(dt)[-1]
  m <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rm(dt)
  sp <- Matrix::Matrix(t(m), sparse = TRUE)
  rm(m)
  gc(verbose = FALSE)
  rownames(sp) <- genes
  colnames(sp) <- paste0(gsm, "_", bc)
  sp
}

align_rows <- function(mat, genes) {
  extra <- setdiff(genes, rownames(mat))
  if (length(extra)) {
    z <- Matrix::Matrix(0, nrow = length(extra), ncol = ncol(mat), sparse = TRUE)
    rownames(z) <- extra
    colnames(z) <- colnames(mat)
    mat <- rbind(mat, z)
  }
  mat[genes, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# 1. GEO CSVs → union sparse matrix → 10x MTX (public processed counts)
# ---------------------------------------------------------------------------
csv_files <- sort(list.files(csv_dir, pattern = "GSM.*_dense\\.csv\\.gz$", full.names = TRUE))
if (!length(csv_files)) stop("No GSM dense CSVs in ", csv_dir)

meta_rows <- list()
mats <- list()
message("Reading ", length(csv_files), " GEO dense CSVs...")
for (fp in csv_files) {
  rec <- parse_fname(basename(fp))
  message("  ", basename(fp))
  sp <- read_geo_csv(fp, rec$gsm)
  meta_rows[[length(meta_rows) + 1L]] <- data.frame(
    barcode = colnames(sp),
    gsm = rec$gsm,
    patient = rec$patient,
    site = rec$site,
    tumor = rec$site != "NORMAL",
    stringsAsFactors = FALSE
  )
  mats[[rec$gsm]] <- sp
}

universe <- sort(unique(unlist(lapply(mats, rownames), use.names = FALSE)))
message("Union genes: ", length(universe), "; aligning...")
mats <- lapply(mats, align_rows, genes = universe)
counts <- do.call(cbind, mats)
rm(mats)
gc(verbose = FALSE)
cell_meta <- do.call(rbind, meta_rows)
stopifnot(identical(colnames(counts), cell_meta$barcode))
rownames(cell_meta) <- cell_meta$barcode

mtx_dir <- file.path(out_dir, "results", "mtx")
mtx_file <- file.path(mtx_dir, "matrix.mtx")
bc_file <- file.path(mtx_dir, "barcodes.tsv")
ft_file <- file.path(mtx_dir, "features.tsv")
message("Writing 10x MTX to ", mtx_dir)
Matrix::writeMM(counts, mtx_file)
write.table(colnames(counts), bc_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(
  cbind(rownames(counts), rownames(counts), "Gene Expression"),
  ft_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Drop in-memory dense-origin object; reload through Seurat ReadMtx.
n_counts_written <- ncol(counts)
n_genes_written <- nrow(counts)
rm(counts)
gc(verbose = FALSE)

# ---------------------------------------------------------------------------
# 2. CreateSeuratObject from the public MTX
# ---------------------------------------------------------------------------
message("CreateSeuratObject from ReadMtx ...")
mtx_counts <- Seurat::ReadMtx(
  mtx = mtx_file,
  cells = bc_file,
  features = ft_file,
  feature.column = 2
)
obj <- Seurat::CreateSeuratObject(
  counts = mtx_counts,
  project = "GSE123902",
  min.cells = 0,
  min.features = 0,
  meta.data = cell_meta[match(colnames(mtx_counts), cell_meta$barcode), ]
)
rm(mtx_counts)
gc(verbose = FALSE)

# SEQC already QC'd. Drop empty barcodes only.
obj <- subset(obj, subset = nCount_RNA > 0)
obj$donor <- obj$patient
obj$tissue <- obj$site

set.seed(SEED)
message("Normalize / PCA / UMAP (", ncol(obj), " cells)...")
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- Seurat::ScaleData(obj, verbose = FALSE)
npcs <- 30L
obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)
obj <- Seurat::FindNeighbors(obj, dims = 1:npcs, verbose = FALSE)
obj <- Seurat::FindClusters(obj, resolution = 0.4, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, dims = 1:npcs, verbose = FALSE)

# Family scores = mean log1p CP10k of locked genes (mouse-style).
# Not Seurat AddModuleScore control-pool residuals. CLDN4 is not in IFN / MHC / TJ.
data_mat <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
obj$IFN_score <- marker_score_vec(data_mat, IFN_GENES)
obj$MHC_score <- marker_score_vec(data_mat, MHC_GENES)
obj$TJ_score <- marker_score_vec(data_mat, TJ_GENES)
obj$lineage <- assign_lineage(data_mat)
obj$malignant <- obj$lineage == "epithelial" & obj$tumor
obj$is_tnk <- obj$lineage == "tnk"
cldn4 <- if ("CLDN4" %in% rownames(data_mat)) as.numeric(data_mat["CLDN4", ]) else rep(0, ncol(obj))
obj$CLDN4 <- cldn4
obj$CLDN4_pos <- cldn4 > 0

# ---------------------------------------------------------------------------
# 3. Patient / donor table (honest unit)
# ---------------------------------------------------------------------------
md <- obj[[]]
# tumor-sample cells only for T/NK fraction and malignant CLDN4
tumor_md <- md[md$tumor, , drop = FALSE]

patients <- sort(unique(md$patient))
pat_rows <- lapply(patients, function(p) {
  d <- md[md$patient == p, , drop = FALSE]
  t <- d[d$tumor, , drop = FALSE]
  mal <- t[t$malignant, , drop = FALSE]
  tnk <- t[t$is_tnk, , drop = FALSE]
  sites <- unique(d$site)
  tumor_sites <- unique(t$site)
  data.frame(
    patient = p,
    donor = p,
    sites = paste(sites, collapse = "|"),
    tumor_sites = if (length(tumor_sites)) paste(tumor_sites, collapse = "|") else "NONE",
    has_tumor = nrow(t) > 0,
    n_cells_all = nrow(d),
    n_cells_tumor = nrow(t),
    n_malignant = nrow(mal),
    n_tnk = nrow(tnk),
    n_epithelial_all = sum(d$lineage == "epithelial"),
    n_myeloid_tumor = sum(t$lineage == "myeloid"),
    n_b_tumor = sum(t$lineage == "b"),
    frac_tnk = if (nrow(t)) mean(t$is_tnk) else NA_real_,
    cldn4_mean = if (nrow(mal)) mean(mal$CLDN4) else NA_real_,
    cldn4_pct_pos = if (nrow(mal)) mean(mal$CLDN4_pos) * 100 else NA_real_,
    ifn_mean = if (nrow(mal) && "IFN_score" %in% names(mal)) mean(mal$IFN_score) else NA_real_,
    mhc_mean = if (nrow(mal) && "MHC_score" %in% names(mal)) mean(mal$MHC_score) else NA_real_,
    tj_mean = if (nrow(mal) && "TJ_score" %in% names(mal)) mean(mal$TJ_score) else NA_real_,
    tnk_cldn4_mean = if (nrow(tnk)) mean(tnk$CLDN4) else NA_real_,
    stringsAsFactors = FALSE
  )
})
pat <- do.call(rbind, pat_rows)
pat$eligible <- pat$has_tumor & pat$n_malignant >= MIN_MAL & pat$n_tnk >= MIN_TNK
el <- pat[pat$eligible, , drop = FALSE]

s_mean_tnk <- spearman_ci(el$cldn4_mean, el$frac_tnk)
s_pct_tnk <- spearman_ci(el$cldn4_pct_pos, el$frac_tnk)
s_ifn <- spearman_ci(el$cldn4_mean, el$ifn_mean)
s_mhc <- spearman_ci(el$cldn4_mean, el$mhc_mean)
s_tj <- spearman_ci(el$cldn4_mean, el$tj_mean)

# Q4 vs Q1 on malignant CLDN4 mean among eligible donors.
q4q1 <- list(ran = FALSE, n_q1 = 0L, n_q4 = 0L, note = "n below Q4/Q1 floor")
if (nrow(el) >= MIN_N_Q4Q1) {
  qs <- as.numeric(quantile(el$cldn4_mean, probs = c(0.25, 0.75), na.rm = TRUE, type = 7))
  q1 <- el[el$cldn4_mean <= qs[1], , drop = FALSE]
  q4 <- el[el$cldn4_mean >= qs[2], , drop = FALSE]
  mw <- function(a, b) {
    if (length(a) < 2 || length(b) < 2) {
      return(list(delta = NA_real_, p = NA_real_, n1 = length(a), n4 = length(b)))
    }
    wt <- suppressWarnings(wilcox.test(b, a, exact = FALSE))
    list(delta = median(b) - median(a), p = unname(wt$p.value), n1 = length(a), n4 = length(b))
  }
  q4q1 <- list(
    ran = TRUE,
    n_q1 = nrow(q1),
    n_q4 = nrow(q4),
    q1_cut = qs[1],
    q4_cut = qs[2],
    ifn = mw(q1$ifn_mean, q4$ifn_mean),
    mhc = mw(q1$mhc_mean, q4$mhc_mean),
    tj = mw(q1$tj_mean, q4$tj_mean),
    note = sprintf("Q1 n=%d vs Q4 n=%d on malignant CLDN4 mean", nrow(q1), nrow(q4))
  )
}

# Restriction: malignant CLDN4 vs same-patient T/NK CLDN4 (not the infiltration claim).
restriction_n <- sum(el$cldn4_mean > el$tnk_cldn4_mean, na.rm = TRUE)

# ---------------------------------------------------------------------------
# 4. Seurat DimPlot / VlnPlot + patient scatter
# ---------------------------------------------------------------------------
fig_dir <- file.path(out_dir, "results", "figures")
save_gg <- function(p, stem, w = 8, h = 6) {
  ggsave(file.path(fig_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 150)
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), p, width = w, height = h)
}

Idents(obj) <- obj$lineage
p_dim_lin <- DimPlot(obj, group.by = "lineage", reduction = "umap", pt.size = 0.1) +
  ggtitle("GSE123902 Seurat UMAP — marker lineage (CLDN4 not used to assign)")
save_gg(p_dim_lin, "fig_dimplot_lineage", 8, 6)

p_dim_pat <- DimPlot(obj, group.by = "patient", reduction = "umap", pt.size = 0.1) +
  ggtitle("GSE123902 Seurat UMAP — donor (LX ID)") +
  theme(legend.text = element_text(size = 7))
save_gg(p_dim_pat, "fig_dimplot_patient", 9, 6)

p_dim_tis <- DimPlot(obj, group.by = "site", reduction = "umap", pt.size = 0.1) +
  ggtitle("GSE123902 Seurat UMAP — tissue")
save_gg(p_dim_tis, "fig_dimplot_tissue", 8, 6)

p_vln <- VlnPlot(obj, features = "CLDN4", group.by = "lineage", pt.size = 0) +
  ggtitle("CLDN4 (log1p CP10k) by marker lineage")
save_gg(p_vln, "fig_vlnplot_cldn4", 8, 5)

if ("CLDN4" %in% rownames(obj)) {
  p_feat <- FeaturePlot(obj, features = "CLDN4", pt.size = 0.1) +
    ggtitle("CLDN4")
  save_gg(p_feat, "fig_featureplot_cldn4", 7, 6)
}

p_sc <- ggplot(el, aes(x = cldn4_mean, y = frac_tnk, label = patient)) +
  geom_point(aes(color = tumor_sites), size = 2.6) +
  geom_text(vjust = -0.7, size = 2.4) +
  labs(
    x = "Malignant/epithelial CLDN4 mean (log1p CP10k)",
    y = "T/NK fraction (tumor cells, same donor)",
    title = "GSE123902 patient/donor: CLDN4 vs T/NK",
    subtitle = fmt_rho(s_mean_tnk)
  ) +
  theme_bw(base_size = 11)
save_gg(p_sc, "fig_patient_cldn4_vs_tnk", 8, 6)

p_sc2 <- ggplot(el, aes(x = cldn4_pct_pos, y = frac_tnk, label = patient)) +
  geom_point(aes(color = tumor_sites), size = 2.6) +
  geom_text(vjust = -0.7, size = 2.4) +
  labs(
    x = "Malignant/epithelial CLDN4 % positive",
    y = "T/NK fraction (tumor cells, same donor)",
    title = "GSE123902 patient/donor: CLDN4 %pos vs T/NK",
    subtitle = fmt_rho(s_pct_tnk)
  ) +
  theme_bw(base_size = 11)
save_gg(p_sc2, "fig_patient_cldn4pct_vs_tnk", 8, 6)

# ---------------------------------------------------------------------------
# 5. Tables + FINDING.md
# ---------------------------------------------------------------------------
write.table(pat, file.path(out_dir, "results", "tables", "patient_level_cldn4_tnk.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

spear_tab <- data.frame(
  contrast = c(
    "malignant_CLDN4_mean vs T/NK_fraction",
    "malignant_CLDN4_pct_pos vs T/NK_fraction",
    "malignant_CLDN4_mean vs IFN",
    "malignant_CLDN4_mean vs MHC-I/APM",
    "malignant_CLDN4_mean vs TJ(CLDN4 held out)"
  ),
  n = c(s_mean_tnk$n, s_pct_tnk$n, s_ifn$n, s_mhc$n, s_tj$n),
  rho = c(s_mean_tnk$rho, s_pct_tnk$rho, s_ifn$rho, s_mhc$rho, s_tj$rho),
  ci_lo = c(s_mean_tnk$lo, s_pct_tnk$lo, s_ifn$lo, s_mhc$lo, s_tj$lo),
  ci_hi = c(s_mean_tnk$hi, s_pct_tnk$hi, s_ifn$hi, s_mhc$hi, s_tj$hi),
  p = c(s_mean_tnk$p, s_pct_tnk$p, s_ifn$p, s_mhc$p, s_tj$p),
  stringsAsFactors = FALSE
)
write.table(spear_tab, file.path(out_dir, "results", "tables", "donor_spearman.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

q_tab <- if (isTRUE(q4q1$ran)) {
  data.frame(
    family = c("IFN", "MHC-I/APM", "TJ(CLDN4 held out)"),
    n_Q1 = q4q1$n_q1,
    n_Q4 = q4q1$n_q4,
    delta_median_Q4_minus_Q1 = c(q4q1$ifn$delta, q4q1$mhc$delta, q4q1$tj$delta),
    wilcox_p = c(q4q1$ifn$p, q4q1$mhc$p, q4q1$tj$p),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(family = character(), n_Q1 = integer(), n_Q4 = integer(),
             delta_median_Q4_minus_Q1 = numeric(), wilcox_p = numeric())
}
write.table(q_tab, file.path(out_dir, "results", "tables", "malignant_q4q1_ifn_mhc_tj.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

lin_tab <- as.data.frame(table(lineage = obj$lineage, site = obj$site), stringsAsFactors = FALSE)
write.table(lin_tab, file.path(out_dir, "results", "tables", "lineage_by_site.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

ifn_cov <- present(IFN_GENES, rownames(obj))
mhc_cov <- present(MHC_GENES, rownames(obj))
tj_cov <- present(TJ_GENES, rownames(obj))
cov_tab <- data.frame(
  set = c("IFN_core", "MHC_I_APM", "TJ_no_CLDN4"),
  n_locked = c(length(IFN_GENES), length(MHC_GENES), length(TJ_GENES)),
  n_present = c(length(ifn_cov), length(mhc_cov), length(tj_cov)),
  missing = c(
    paste(setdiff(IFN_GENES, ifn_cov), collapse = ","),
    paste(setdiff(MHC_GENES, mhc_cov), collapse = ","),
    paste(setdiff(TJ_GENES, tj_cov), collapse = ",")
  ),
  stringsAsFactors = FALSE
)
write.table(cov_tab, file.path(out_dir, "results", "tables", "geneset_coverage.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Patient markdown table
pat_md_lines <- apply(pat, 1, function(r) {
  sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    r[["patient"]], r[["tumor_sites"]], r[["n_cells_tumor"]], r[["n_malignant"]],
    r[["n_tnk"]],
    ifelse(is.na(as.numeric(r[["frac_tnk"]])), "NA", sprintf("%.3f", as.numeric(r[["frac_tnk"]]))),
    ifelse(is.na(as.numeric(r[["cldn4_mean"]])), "NA", sprintf("%.3f", as.numeric(r[["cldn4_mean"]]))),
    ifelse(is.na(as.numeric(r[["cldn4_pct_pos"]])), "NA", sprintf("%.1f", as.numeric(r[["cldn4_pct_pos"]]))),
    ifelse(is.na(as.numeric(r[["ifn_mean"]])), "NA", sprintf("%.3f", as.numeric(r[["ifn_mean"]]))),
    ifelse(is.na(as.numeric(r[["mhc_mean"]])), "NA", sprintf("%.3f", as.numeric(r[["mhc_mean"]]))),
    ifelse(is.na(as.numeric(r[["tj_mean"]])), "NA", sprintf("%.3f", as.numeric(r[["tj_mean"]]))),
    ifelse(as.logical(r[["eligible"]]), "yes", "no")
  )
})

q_block <- if (isTRUE(q4q1$ran)) {
  paste0(
    "Malignant CLDN4-mean quartiles among eligible donors. ", q4q1$note, ".\n\n",
    "| family | n_Q1 | n_Q4 | Δ median (Q4−Q1) | MW p |\n",
    "|---|---:|---:|---:|---:|\n",
    sprintf("| IFN | %d | %d | %+.3f | %.3g |\n", q4q1$n_q1, q4q1$n_q4, q4q1$ifn$delta, q4q1$ifn$p),
    sprintf("| MHC-I/APM | %d | %d | %+.3f | %.3g |\n", q4q1$n_q1, q4q1$n_q4, q4q1$mhc$delta, q4q1$mhc$p),
    sprintf("| TJ (CLDN4 held out) | %d | %d | %+.3f | %.3g |\n", q4q1$n_q1, q4q1$n_q4, q4q1$tj$delta, q4q1$tj$p)
  )
} else {
  paste0("Q4 vs Q1 **not run**: eligible n=", nrow(el), " is below the floor of ", MIN_N_Q4Q1, ".\n")
}

si <- utils::capture.output(utils::sessionInfo())
seurat_ver <- as.character(utils::packageVersion("Seurat"))

finding <- paste0(
  "# FINDING — Seurat GSE123902 CLDN4-only (patient/donor unit)\n\n",
  "ADDITIVE. **CLDN4 only.** Thesis already correct. This is a Seurat-native re-score of public processed ",
  "human LUAD/NSCLC [GSE123902](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123902) ",
  "(Laughney et al., *Nat Med* 2020, PMID 32042191). SuperSeries GSE123904 / mouse GSE123903 are not used. ",
  "No TACSTD2∩CLDN4 dual-high gate. No GSE148071 merge. No Python-only primary.\n\n",
  "Primary engine: **R + Seurat ", seurat_ver, "**. `CreateSeuratObject` was called on a 10x-style MTX ",
  "written from the GEO dense unnormalized UMI CSVs (`GSE123902_RAW.tar`, 90.4 MB). ",
  "The 36.5 GB author annotated H5 was not required. Honest unit = **patient/donor (LX ID)**.\n\n",
  "## Verdict\n\n",
  "Seurat ran. A patient-level table exists (`results/tables/patient_level_cldn4_tnk.tsv`). ",
  "Eligible tumor donors ( ≥", MIN_MAL, " marker-epithelial cells in tumor and ≥", MIN_TNK, " T/NK): **n=", nrow(el), "**.\n\n",
  "- Malignant CLDN4 **mean** vs T/NK fraction: ", fmt_rho(s_mean_tnk), ".\n",
  "- Malignant CLDN4 **%pos** vs T/NK fraction: ", fmt_rho(s_pct_tnk), ".\n",
  "- Malignant CLDN4 mean vs IFN: ", fmt_rho(s_ifn), ".\n",
  "- Malignant CLDN4 mean vs MHC-I/APM: ", fmt_rho(s_mhc), ".\n",
  "- Malignant CLDN4 mean vs TJ (CLDN4 held out): ", fmt_rho(s_tj), ".\n\n",
  "Malignant CLDN4 is higher than same-donor T/NK CLDN4 in **", restriction_n, "/", nrow(el),
  "** eligible donors (epithelial restriction; not an immune-cold claim).\n\n",
  "n=", nrow(el), " is the honest ceiling. Cell-level p-values are not the claim. ",
  "Do not write this as a failed audit of the thesis.\n\n",
  "## Honest n\n\n",
  "| item | n | note |\n",
  "|---|---:|---|\n",
  sprintf("| GEO dense CSVs | %d | 17 files in GSE123902_RAW.tar |\n", length(csv_files)),
  sprintf("| MTX barcodes written | %d | union gene space %d |\n", n_counts_written, n_genes_written),
  sprintf("| Seurat cells (nCount_RNA>0) | %d | CreateSeuratObject from MTX |\n", ncol(obj)),
  sprintf("| Donors / LX IDs | %d | filename MSK_LX* |\n", length(unique(md$patient))),
  sprintf("| Tumor donors | %d | primary + metastasis; LX685 is normal-only |\n", sum(pat$has_tumor)),
  sprintf("| Marker epithelial in tumor (malignant for this table) | %d | EPCAM/KRT* vs T/NK/myeloid/B; CLDN4 not used |\n", sum(obj$malignant)),
  sprintf("| T/NK in tumor samples | %d | CD3D/CD3E/CD8A/NKG7/GNLY/KLRD1 |\n", sum(tumor_md$is_tnk)),
  sprintf("| Eligible Spearman donors | **%d** | ≥%d malignant and ≥%d T/NK |\n", nrow(el), MIN_MAL, MIN_TNK),
  "| Dual-high TACSTD2 ∩ CLDN4 | not defined | CLDN4-only |\n",
  "| GSE148071 cells | 0 | not merged |\n",
  "| Author 36.5 GB H5 | not used | GEO CSV → MTX is the public processed object |\n",
  "| Marker epithelial ≠ CNV-malignant | yes | SEQC dense UMI; no inferCNV |\n\n",
  "Paper QC atlas is 41,384 cells. We do not substitute that n. SEQC CSVs hold ",
  n_counts_written, " barcodes before the empty-barcode drop.\n\n",
  "## Gate\n\n",
  "| file | public? | used |\n",
  "|---|---|---|\n",
  "| `GSE123902_RAW.tar` (17 dense UMI CSVs) | yes, 90.4 MB | **yes** — written to 10x MTX, then `ReadMtx` + `CreateSeuratObject` |\n",
  "| 10x MTX (`results/mtx/`) | derived here | **yes** — Seurat input |\n",
  "| Author `PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5` | yes, 36.5 GB | no — not required once GEO UMI CSVs exist |\n",
  "| GSE148071 / other cohorts | — | **no** |\n\n",
  "## Locked choices\n\n",
  "- Lineage is a four-way marker argmax on Seurat `LogNormalize` (log1p CP10k). ",
  "Keep if top ≥ 0.12 and top ≥ 1.15 × second. CLDN4 is never a lineage marker.\n",
  "- Malignant = marker epithelial **in tumor** (primary or metastasis). Matched normal is not scored as malignant.\n",
  "- T/NK fraction = marker T/NK / tumor-sample cells of that donor.\n",
  "- IFN = locked Hallmark IFNα ∩ IFNγ core (present ", length(ifn_cov), "/", length(IFN_GENES), ").\n",
  "- MHC-I/APM = curated antigen-presentation set (present ", length(mhc_cov), "/", length(MHC_GENES), ").\n",
  "- TJ = epithelial tight-junction genes with **CLDN4 held out** (present ", length(tj_cov), "/", length(TJ_GENES), ").\n",
  "- Eligible Spearman n requires ≥", MIN_MAL, " malignant and ≥", MIN_TNK, " T/NK. Q4 vs Q1 requires eligible n ≥ ", MIN_N_Q4Q1, ".\n\n",
  "## Patient / donor table\n\n",
  "Machine table: `results/tables/patient_level_cldn4_tnk.tsv`.\n\n",
  "| patient | tumor site | n tumor | n mal | n T/NK | frac T/NK | CLDN4 mean | CLDN4 %pos | IFN | MHC | TJ | eligible |\n",
  "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n",
  paste(pat_md_lines, collapse = "\n"), "\n\n",
  "## Donor-level Spearman (primary)\n\n",
  "| Contrast | n | ρ [95% CI] | p |\n",
  "|---|---:|---|---:|\n",
  sprintf("| CLDN4 mean vs T/NK fraction | %d | %+.3f [%.3f, %.3f] | %.3g |\n", s_mean_tnk$n, s_mean_tnk$rho, s_mean_tnk$lo, s_mean_tnk$hi, s_mean_tnk$p),
  sprintf("| CLDN4 %%pos vs T/NK fraction | %d | %+.3f [%.3f, %.3f] | %.3g |\n", s_pct_tnk$n, s_pct_tnk$rho, s_pct_tnk$lo, s_pct_tnk$hi, s_pct_tnk$p),
  sprintf("| CLDN4 mean vs IFN | %d | %+.3f [%.3f, %.3f] | %.3g |\n", s_ifn$n, s_ifn$rho, s_ifn$lo, s_ifn$hi, s_ifn$p),
  sprintf("| CLDN4 mean vs MHC-I/APM | %d | %+.3f [%.3f, %.3f] | %.3g |\n", s_mhc$n, s_mhc$rho, s_mhc$lo, s_mhc$hi, s_mhc$p),
  sprintf("| CLDN4 mean vs TJ (CLDN4 held out) | %d | %+.3f [%.3f, %.3f] | %.3g |\n", s_tj$n, s_tj$rho, s_tj$lo, s_tj$hi, s_tj$p),
  "\n",
  "## Malignant Q4 vs Q1 IFN / MHC / TJ\n\n",
  q_block, "\n",
  "## Seurat plots\n\n",
  "- `results/figures/fig_dimplot_lineage.png` — `DimPlot` by marker lineage\n",
  "- `results/figures/fig_dimplot_patient.png` — `DimPlot` by donor\n",
  "- `results/figures/fig_dimplot_tissue.png` — `DimPlot` by tissue\n",
  "- `results/figures/fig_vlnplot_cldn4.png` — `VlnPlot` CLDN4 by lineage\n",
  "- `results/figures/fig_featureplot_cldn4.png` — `FeaturePlot` CLDN4\n",
  "- `results/figures/fig_patient_cldn4_vs_tnk.png` — donor scatter (mean)\n",
  "- `results/figures/fig_patient_cldn4pct_vs_tnk.png` — donor scatter (%pos)\n\n",
  "## What this does not claim\n\n",
  "- Cell-level p-values are not the claim. n_cells is large by construction.\n",
  "- Marker epithelial is **not** a CNV-malignant call.\n",
  "- n=", nrow(el), " tumor donors is small; Fisher-z CIs are wide.\n",
  "- Restriction (CLDN4 in epithelium vs T/NK) is not an infiltration / immune-cold claim.\n",
  "- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.\n",
  "- No ICI / MPR / RECIST / survival test.\n",
  "- Not merged with GSE148071 or any other cohort.\n\n",
  "## Reproduce\n\n",
  "```bash\n",
  "bash methods/seurat_gse123902_cldn4/scripts/download.sh\n",
  "Rscript methods/seurat_gse123902_cldn4/scripts/run_seurat_gse123902.R\n",
  "```\n\n",
  "Seurat ", seurat_ver, ". R ", paste(R.version$major, R.version$minor, sep = "."), ".\n"
)

writeLines(finding, file.path(out_dir, "FINDING.md"))

summary <- list(
  accession = "GSE123902",
  engine = "R_Seurat",
  seurat = seurat_ver,
  n_cells = ncol(obj),
  n_donors = length(unique(md$patient)),
  n_eligible = nrow(el),
  spearman_cldn4_mean_vs_tnk = s_mean_tnk,
  spearman_cldn4_pct_vs_tnk = s_pct_tnk,
  q4q1_ran = isTRUE(q4q1$ran)
)
jsonlite::write_json(summary, file.path(out_dir, "results", "summary.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(si, file.path(out_dir, "results", "sessionInfo.txt"))

message("DONE. Eligible donors n=", nrow(el), " ; table + FINDING.md written.")
message(fmt_rho(s_mean_tnk))
