#!/usr/bin/env Rscript
# ADDITIVE. Cldn4-only. GSE179502 FACS neoplastic epithelium, Lkb1 restore vs not.
# MUST use Seurat::CreateSeuratObject on the public GEO mtx.
# Honest unit = mouse (n=3 Restored vs 3 NonRestored). No dual-high. No private 8 KL.

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop(paste0(
      "Seurat is not installed. Stopped without switching to Python.\n",
      "Install R package Seurat, then re-run: Rscript methods/seurat_gse179502_cldn4/analyze.R"
    ), call. = FALSE)
  }
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(jsonlite)
})

library(patchwork)

HERE <- tryCatch(
  {
    ca <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", ca[grepl("^--file=", ca)])
    if (length(f)) dirname(normalizePath(f[1])) else getwd()
  },
  error = function(e) getwd()
)
source(file.path(HERE, "gene_sets.R"))

DATA <- file.path(HERE, "data")
TABLES <- file.path(HERE, "tables")
FIGS <- file.path(HERE, "figures")
dir.create(DATA, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGS, recursive = TRUE, showWarnings = FALSE)

GEO_BASE <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179502"
GEO_FILES <- list(
  "barcodes.tsv.gz" = paste0(GEO_BASE, "/suppl/GSE179502_XTR_sorted_scRNAseq_barcodes.tsv.gz"),
  "features.tsv.gz" = paste0(GEO_BASE, "/suppl/GSE179502_XTR_sorted_scRNAseq_features.tsv.gz"),
  "matrix.mtx.gz" = paste0(GEO_BASE, "/suppl/GSE179502_XTR_sorted_scRNAseq_matrix.mtx.gz")
)

MIN_GENES <- 200L
MIN_UMI <- 500
MAX_MITO <- 20
MIN_N_Q4 <- 8L

fmt <- function(x, nd = 3) {
  if (length(x) != 1 || is.null(x) || !is.finite(as.numeric(x))) return("—")
  sprintf(paste0("%.", nd, "f"), as.numeric(x))
}

ensure_geo <- function() {
  for (nm in names(GEO_FILES)) {
    dest <- file.path(DATA, nm)
    if (file.exists(dest) && file.info(dest)$size > 1000) next
    url <- GEO_FILES[[nm]]
    message("DOWNLOAD ", url)
    utils::download.file(url, dest, mode = "wb", quiet = FALSE)
  }
  bc <- file.path(DATA, "barcodes.tsv.gz")
  # GEO sometimes ships barcodes uncompressed under a .gz name.
  raw <- readBin(bc, what = "raw", n = 2L)
  if (!(length(raw) >= 2 && raw[1] == as.raw(0x1f) && raw[2] == as.raw(0x8b))) {
    tmp <- file.path(DATA, "barcodes.tsv")
    file.rename(bc, tmp)
    R.utils <- requireNamespace("R.utils", quietly = TRUE)
    if (R.utils) {
      R.utils::gzip(tmp, destname = bc, overwrite = TRUE, remove = TRUE)
    } else {
      system2("gzip", c("-n", "-f", tmp))
      if (!file.exists(bc) && file.exists(paste0(tmp, ".gz"))) {
        file.rename(paste0(tmp, ".gz"), bc)
      }
    }
  }
}

assay_layer <- function(obj, layer) {
  # Seurat v5 prefers layer=; v4 uses slot=.
  out <- try(GetAssayData(obj, assay = "RNA", layer = layer), silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  GetAssayData(obj, assay = "RNA", slot = layer)
}

find_gene <- function(obj, symbol) {
  feats <- rownames(obj)
  if (symbol %in% feats) return(symbol)
  hit <- feats[sub("\\.\\d+$", "", feats) == symbol]
  if (length(hit)) return(hit[1])
  NA_character_
}

present_on_object <- function(obj, genes) {
  vapply(genes, function(g) !is.na(find_gene(obj, g)), logical(1))
}

resolve_genes <- function(obj, genes) {
  out <- unique(na.omit(vapply(genes, function(g) find_gene(obj, g), character(1))))
  as.character(out)
}

gene_vector <- function(mat, gene) {
  if (is.na(gene) || !gene %in% rownames(mat)) {
    return(rep(NA_real_, ncol(mat)))
  }
  as.numeric(mat[gene, ])
}

set_mean <- function(mat, genes) {
  genes <- genes[genes %in% rownames(mat)]
  if (!length(genes)) return(rep(NA_real_, ncol(mat)))
  if (length(genes) == 1L) return(as.numeric(mat[genes, ]))
  as.numeric(Matrix::colMeans(mat[genes, , drop = FALSE]))
}

welch_mwu <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  out <- list(
    n_a = length(a), n_b = length(b),
    mean_a = if (length(a)) mean(a) else NA_real_,
    mean_b = if (length(b)) mean(b) else NA_real_,
    median_a = if (length(a)) median(a) else NA_real_,
    median_b = if (length(b)) median(b) else NA_real_,
    delta_mean = NA_real_, delta_median = NA_real_,
    welch_p = NA_real_, mwu_p = NA_real_, r_rb = NA_real_
  )
  if (!length(a) || !length(b)) return(out)
  out$delta_mean <- out$mean_b - out$mean_a
  out$delta_median <- out$median_b - out$median_a
  if (length(a) >= 2 && length(b) >= 2 && (sd(a) + sd(b)) > 0) {
    out$welch_p <- tryCatch(t.test(b, a)$p.value, error = function(e) NA_real_)
  }
  if (length(a) >= 1 && length(b) >= 1) {
    wt <- tryCatch(wilcox.test(b, a, alternative = "two.sided", exact = FALSE), error = function(e) NULL)
    if (!is.null(wt)) {
      out$mwu_p <- wt$p.value
      # rank-biserial from U: 2U/(n1 n2) - 1. wilcox.test statistic is U for the first sample.
      u <- as.numeric(wt$statistic)
      out$r_rb <- (2 * u) / (length(b) * length(a)) - 1
    }
  }
  out
}

q4_vs_q1 <- function(rank_on, y, min_n = MIN_N_Q4) {
  c <- as.numeric(rank_on); y <- as.numeric(y)
  m <- is.finite(c) & is.finite(y)
  c <- c[m]; y <- y[m]
  n <- length(c)
  out <- list(
    n = n, n_q1 = NA_real_, n_q4 = NA_real_,
    median_q1 = NA_real_, median_q4 = NA_real_,
    delta_median = NA_real_, r_rb = NA_real_, p = NA_real_, usable = FALSE
  )
  n_tail <- n %/% 4L
  if (n < min_n || n_tail < 2) return(out)
  ord <- order(c, method = "radix")
  q1 <- y[ord[seq_len(n_tail)]]
  q4 <- y[ord[(n - n_tail + 1):n]]
  d <- welch_mwu(q1, q4)
  out$n_q1 <- d$n_a
  out$n_q4 <- d$n_b
  out$median_q1 <- d$median_a
  out$median_q4 <- d$median_b
  out$delta_median <- d$delta_median
  out$r_rb <- d$r_rb
  out$p <- d$mwu_p
  out$usable <- TRUE
  out
}

pos_vs_neg <- function(counts, y) {
  d <- welch_mwu(y[counts <= 0], y[counts > 0])
  list(
    n = sum(is.finite(counts) & is.finite(y)),
    n_q1 = d$n_a, n_q4 = d$n_b,
    median_q1 = d$median_a, median_q4 = d$median_b,
    delta_median = d$delta_median, r_rb = d$r_rb, p = d$mwu_p,
    usable = d$n_a >= 2 && d$n_b >= 2
  )
}

spearman_xy <- function(x, y) {
  xa <- as.numeric(x); ya <- as.numeric(y)
  m <- is.finite(xa) & is.finite(ya)
  n <- sum(m)
  if (n < 4) return(list(rho = NA_real_, p = NA_real_, n = n))
  ct <- suppressWarnings(cor.test(xa[m], ya[m], method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = n)
}

partial_spearman <- function(x, y, z) {
  xa <- as.numeric(x); ya <- as.numeric(y); za <- as.numeric(z)
  m <- is.finite(xa) & is.finite(ya) & is.finite(za)
  n <- sum(m)
  if (n < 4) return(list(rho = NA_real_, p = NA_real_, n = n))
  rx <- rank(xa[m]); ry <- rank(ya[m]); rz <- rank(za[m])
  rx <- residuals(lm(rx ~ rz))
  ry <- residuals(lm(ry ~ rz))
  ct <- suppressWarnings(cor.test(rx, ry, method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = n)
}

list_to_row <- function(lst) {
  as.data.frame(lst, stringsAsFactors = FALSE)
}

save_plot <- function(p, stem, w = 7, h = 5) {
  png <- paste0(stem, ".png")
  pdf <- paste0(stem, ".pdf")
  ggsave(png, p, width = w, height = h, dpi = 150, limitsize = FALSE)
  ggsave(pdf, p, width = w, height = h, limitsize = FALSE)
  message("WROTE ", png)
}

# ------------------------------------------------------------------------------
message("=== GSE179502 Seurat Cldn4-only ===")
message("Seurat ", as.character(packageVersion("Seurat")),
        " / SeuratObject ", as.character(packageVersion("SeuratObject")))

ensure_geo()
samples <- read.csv(file.path(DATA, "samples.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(samples) == 6L)

message("READ10X + CreateSeuratObject")
counts <- Read10X(DATA, gene.column = 2)
obj <- CreateSeuratObject(
  counts = counts,
  project = "GSE179502",
  min.cells = 0,
  min.features = 0
)
stopifnot(inherits(obj, "Seurat"))
message("CreateSeuratObject OK: ", ncol(obj), " cells x ", nrow(obj), " genes")

bc <- colnames(obj)
mouse <- sub("_.*$", "", bc)
if (!all(mouse %in% samples$mouse)) {
  stop("Unexpected barcode prefixes: ", paste(setdiff(unique(mouse), samples$mouse), collapse = ", "))
}
map <- samples
rownames(map) <- map$mouse
obj$mouse <- mouse
obj$gsm <- map[mouse, "gsm"]
obj$genotype <- map[mouse, "genotype"]
obj$sex <- map[mouse, "sex"]
obj$treatment <- map[mouse, "treatment"]
obj$cohort <- map[mouse, "cohort"]
obj$restorable <- map[mouse, "restorable"]
obj$lkb1_restored <- map[mouse, "lkb1_restored"]
obj$orig.ident <- mouse

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

n_raw <- ncol(obj)
obj <- subset(obj, subset = nFeature_RNA >= MIN_GENES & nCount_RNA >= MIN_UMI & percent.mt <= MAX_MITO)
message("QC: ", n_raw, " -> ", ncol(obj), " cells")

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

present_symbols <- unique(sub("\\.\\d+$", "", rownames(obj)))
a8 <- load_a8_mouse(file.path(DATA, "a8_families.json"), present_symbols)

modules <- list(
  AT2 = present_subset(AT2, present_symbols),
  IFN_CORE = present_subset(IFN_CORE, present_symbols),
  MHC_CORE = present_subset(MHC_CORE, present_symbols),
  TJ_CORE = present_subset(TJ_CORE, present_symbols),
  IFN_A8 = a8$IFN_A8,
  MHC_A8 = a8$MHC_A8,
  TJ_A8 = a8$TJ_A8
)

cov <- data.frame(
  set = names(modules),
  n_requested = c(length(AT2), length(IFN_CORE), length(MHC_CORE), length(TJ_CORE),
                  NA_integer_, NA_integer_, NA_integer_),
  n_present = vapply(modules, length, integer(1)),
  genes = vapply(modules, function(g) paste(g, collapse = ","), character(1)),
  stringsAsFactors = FALSE
)
write.table(cov[, c("set", "n_requested", "n_present")], file.path(TABLES, "gene_coverage.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Seurat AddModuleScore (control-gene subtracted) — primary module scores
mod_features <- modules[vapply(modules, length, integer(1)) >= 2]
obj <- AddModuleScore(
  obj,
  features = lapply(mod_features, function(g) resolve_genes(obj, g)),
  name = "AMS_",
  ctrl = 50,
  seed = 1
)
# AddModuleScore names AMS_1, AMS_2, ... in feature-list order
ams_names <- paste0("AMS_", seq_along(mod_features))
ams_df <- FetchData(obj, vars = ams_names)
for (i in seq_along(mod_features)) {
  obj[[paste0(names(mod_features)[i], "_AMS")]] <- ams_df[, i]
}

# Mean log-normalized module (same scale as Cldn4 mean)
data_mat <- assay_layer(obj, "data")
count_mat <- assay_layer(obj, "counts")
cldn4_g <- find_gene(obj, "Cldn4")
stk11_g <- find_gene(obj, "Stk11")
tacstd2_g <- find_gene(obj, "Tacstd2")
epcam_g <- find_gene(obj, "Epcam")
stopifnot(!is.na(cldn4_g), !is.na(stk11_g))

obj$Cldn4 <- gene_vector(data_mat, cldn4_g)
obj$Cldn4_count <- gene_vector(count_mat, cldn4_g)
obj$Stk11 <- gene_vector(data_mat, stk11_g)
obj$Tacstd2 <- gene_vector(data_mat, tacstd2_g)
obj$Epcam <- gene_vector(data_mat, epcam_g)
for (nm in names(modules)) {
  obj[[paste0(nm, "_mean")]] <- set_mean(data_mat, resolve_genes(obj, modules[[nm]]))
}

md <- FetchData(obj, vars = unique(c(
  "mouse", "gsm", "cohort", "treatment", "sex", "restorable", "lkb1_restored",
  "nCount_RNA", "nFeature_RNA", "percent.mt",
  "Cldn4", "Cldn4_count", "Stk11", "Tacstd2", "Epcam",
  paste0(names(modules), "_mean"),
  paste0(names(mod_features), "_AMS")
)))

# Mouse-level table
mouse_rows <- lapply(split(md, md$mouse), function(d) {
  srow <- samples[samples$mouse == d$mouse[1], ]
  data.frame(
    mouse = d$mouse[1],
    gsm = srow$gsm,
    cohort = srow$cohort,
    treatment = srow$treatment,
    sex = srow$sex,
    restorable = srow$restorable,
    lkb1_restored = srow$lkb1_restored,
    n_cells_qc = nrow(d),
    n_cells_raw = sum(sub("_.*$", "", bc) == d$mouse[1]),
    frac_Cldn4_pos = mean(d$Cldn4_count > 0),
    median_n_umi = median(d$nCount_RNA),
    Cldn4 = mean(d$Cldn4),
    Stk11 = mean(d$Stk11),
    Tacstd2 = mean(d$Tacstd2),
    Epcam = mean(d$Epcam),
    AT2 = mean(d$AT2_mean),
    IFN_CORE = mean(d$IFN_CORE_mean),
    MHC_CORE = mean(d$MHC_CORE_mean),
    TJ_CORE = mean(d$TJ_CORE_mean),
    IFN_A8 = mean(d$IFN_A8_mean),
    MHC_A8 = mean(d$MHC_A8_mean),
    TJ_A8 = mean(d$TJ_A8_mean),
    AT2_AMS = mean(d$AT2_AMS),
    IFN_CORE_AMS = mean(d$IFN_CORE_AMS),
    MHC_CORE_AMS = mean(d$MHC_CORE_AMS),
    TJ_CORE_AMS = mean(d$TJ_CORE_AMS),
    IFN_A8_AMS = mean(d$IFN_A8_AMS),
    MHC_A8_AMS = mean(d$MHC_A8_AMS),
    TJ_A8_AMS = mean(d$TJ_A8_AMS),
    stringsAsFactors = FALSE
  )
})
mice <- do.call(rbind, mouse_rows)
mice <- mice[match(c("CM0875", "ZR1932", "ZR1966", "CM0879", "CM0884", "ZR1969"), mice$mouse), ]
rownames(mice) <- NULL
write.table(mice, file.path(TABLES, "mouse_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mice, file.path(TABLES, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("WROTE mouse-level table: ", nrow(mice), " mice")

features_mouse <- c(
  "Cldn4", "frac_Cldn4_pos", "Stk11", "Tacstd2", "Epcam",
  "AT2", "IFN_CORE", "MHC_CORE", "TJ_CORE", "IFN_A8", "MHC_A8", "TJ_A8",
  "AT2_AMS", "IFN_CORE_AMS", "MHC_CORE_AMS", "TJ_CORE_AMS",
  "IFN_A8_AMS", "MHC_A8_AMS", "TJ_A8_AMS"
)
contrast_rows <- list()
nr <- mice$lkb1_restored == FALSE
rs <- mice$lkb1_restored == TRUE
for (ft in features_mouse) {
  d <- welch_mwu(mice[[ft]][nr], mice[[ft]][rs])
  contrast_rows[[length(contrast_rows) + 1]] <- cbind(
    data.frame(contrast = "Restored_vs_NonRestored", unit = "mouse", feature = ft, stringsAsFactors = FALSE),
    list_to_row(d)
  )
}
# cell-level descriptive (not the n)
cell_feats <- c("Cldn4", "Stk11", "AT2_mean", "IFN_A8_mean", "MHC_A8_mean", "TJ_A8_mean", "AT2_AMS", "IFN_A8_AMS")
for (ft in cell_feats) {
  d <- welch_mwu(md[[ft]][!md$lkb1_restored], md[[ft]][md$lkb1_restored])
  contrast_rows[[length(contrast_rows) + 1]] <- cbind(
    data.frame(contrast = "Restored_vs_NonRestored", unit = "cell_descriptive", feature = ft, stringsAsFactors = FALSE),
    list_to_row(d)
  )
}
contrasts <- do.call(rbind, contrast_rows)
write.table(contrasts, file.path(TABLES, "restore_contrasts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

families <- c("IFN_A8_mean", "MHC_A8_mean", "TJ_A8_mean", "AT2_mean",
              "IFN_CORE_mean", "MHC_CORE_mean", "TJ_CORE_mean",
              "IFN_A8_AMS", "MHC_A8_AMS", "TJ_A8_AMS", "AT2_AMS", "Stk11")
fam_label <- function(x) sub("_mean$", "", x)

hl_rows <- list()
strata <- list(
  all_cells = md,
  NonRestored = md[!md$lkb1_restored, ],
  Restored = md[md$lkb1_restored, ]
)
for (sn in names(strata)) {
  sub <- strata[[sn]]
  for (fam in families) {
    q <- q4_vs_q1(sub$Cldn4, sub[[fam]])
    det <- pos_vs_neg(sub$Cldn4_count, sub[[fam]])
    sp <- spearman_xy(sub$Cldn4, sub[[fam]])
    psp <- partial_spearman(sub$Cldn4, sub[[fam]], log1p(sub$nCount_RNA))
    hl_rows[[length(hl_rows) + 1]] <- cbind(
      data.frame(stratum = sn, gate = "Q4_vs_Q1_rank_tails", unit = "cell_descriptive",
                 family = fam_label(fam), stringsAsFactors = FALSE),
      list_to_row(q),
      data.frame(spearman_rho = sp$rho, spearman_p = sp$p, spearman_n = sp$n,
                 partial_rho_umi = psp$rho, partial_p_umi = psp$p)
    )
    hl_rows[[length(hl_rows) + 1]] <- cbind(
      data.frame(stratum = sn, gate = "Cldn4pos_vs_neg", unit = "cell_descriptive",
                 family = fam_label(fam), stringsAsFactors = FALSE),
      list_to_row(det),
      data.frame(spearman_rho = sp$rho, spearman_p = sp$p, spearman_n = sp$n,
                 partial_rho_umi = psp$rho, partial_p_umi = psp$p)
    )
  }
}

per_mouse_hl <- list()
for (ms in unique(md$mouse)) {
  sub <- md[md$mouse == ms, ]
  for (fam in families) {
    q <- q4_vs_q1(sub$Cldn4, sub[[fam]])
    det <- pos_vs_neg(sub$Cldn4_count, sub[[fam]])
    per_mouse_hl[[length(per_mouse_hl) + 1]] <- cbind(
      data.frame(mouse = ms, cohort = sub$cohort[1], lkb1_restored = sub$lkb1_restored[1],
                 gate = "Q4_vs_Q1_rank_tails", family = fam_label(fam), stringsAsFactors = FALSE),
      list_to_row(q)
    )
    per_mouse_hl[[length(per_mouse_hl) + 1]] <- cbind(
      data.frame(mouse = ms, cohort = sub$cohort[1], lkb1_restored = sub$lkb1_restored[1],
                 gate = "Cldn4pos_vs_neg", family = fam_label(fam), stringsAsFactors = FALSE),
      list_to_row(det)
    )
  }
}
per_mouse_hl_df <- do.call(rbind, per_mouse_hl)
write.table(per_mouse_hl_df, file.path(TABLES, "cldn4_highlow_by_mouse.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

for (gate in c("Q4_vs_Q1_rank_tails", "Cldn4pos_vs_neg")) {
  for (fam in unique(vapply(families, fam_label, character(1)))) {
    dlt <- per_mouse_hl_df$delta_median[
      per_mouse_hl_df$family == fam & per_mouse_hl_df$gate == gate & per_mouse_hl_df$usable
    ]
    dlt <- dlt[is.finite(dlt)]
    if (length(dlt) >= 4) {
      wt <- suppressWarnings(wilcox.test(dlt, alternative = "two.sided", exact = TRUE))
      w <- as.numeric(wt$statistic); p <- wt$p.value
    } else {
      w <- NA_real_; p <- NA_real_
    }
    hl_rows[[length(hl_rows) + 1]] <- data.frame(
      stratum = "within_mouse_then_wilcoxon", gate = gate, unit = "mouse",
      family = fam, n = length(dlt), n_q1 = NA_real_, n_q4 = NA_real_,
      median_q1 = NA_real_, median_q4 = NA_real_,
      delta_median = if (length(dlt)) median(dlt) else NA_real_,
      r_rb = if (length(dlt)) mean(dlt > 0) - mean(dlt < 0) else NA_real_,
      p = p, usable = length(dlt) >= 4,
      spearman_rho = NA_real_, spearman_p = NA_real_, spearman_n = length(dlt),
      partial_rho_umi = NA_real_, partial_p_umi = NA_real_,
      wilcoxon_W = w,
      n_mice_delta_neg = sum(dlt < 0), n_mice_delta_pos = sum(dlt > 0),
      stringsAsFactors = FALSE
    )
  }
}

for (fam in c("AT2", "IFN_A8", "MHC_A8", "TJ_A8", "IFN_CORE", "MHC_CORE", "TJ_CORE", "Stk11")) {
  col <- if (fam == "Stk11") "Stk11" else fam
  sp <- spearman_xy(mice$Cldn4, mice[[col]])
  hl_rows[[length(hl_rows) + 1]] <- data.frame(
    stratum = "mouse_mean", gate = "mouse_mean_spearman", unit = "mouse",
    family = fam, n = sp$n, n_q1 = NA_real_, n_q4 = NA_real_,
    median_q1 = NA_real_, median_q4 = NA_real_,
    delta_median = NA_real_, r_rb = NA_real_, p = sp$p, usable = sp$n >= 4,
    spearman_rho = sp$rho, spearman_p = sp$p, spearman_n = sp$n,
    partial_rho_umi = NA_real_, partial_p_umi = NA_real_,
    stringsAsFactors = FALSE
  )
}

# bind with fill
all_cols <- unique(unlist(lapply(hl_rows, names)))
hl_rows <- lapply(hl_rows, function(d) {
  miss <- setdiff(all_cols, names(d))
  for (m in miss) d[[m]] <- NA
  d[all_cols]
})
hl <- do.call(rbind, hl_rows)
write.table(hl, file.path(TABLES, "cldn4_highlow.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
ifn_tab <- hl[
  hl$family %in% c("IFN_A8", "IFN_CORE", "MHC_A8", "MHC_CORE", "TJ_A8", "TJ_CORE", "AT2",
                   "IFN_A8_AMS", "MHC_A8_AMS", "TJ_A8_AMS", "AT2_AMS") &
    hl$gate %in% c("Q4_vs_Q1_rank_tails", "Cldn4pos_vs_neg", "mouse_mean_spearman"),
]
write.table(ifn_tab, file.path(TABLES, "cldn4_highlow_ifn.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

grab <- function(feature, unit = "mouse") {
  contrasts[contrasts$feature == feature & contrasts$unit == unit, ][1, ]
}
one <- data.frame(
  dataset = "GSE179502",
  model = "KT;Lkb1XTR KRAS lung, FACS neoplastic scRNA",
  contrast = "Restored vs NonRestored (author GEO cohort)",
  engine = paste0("Seurat ", packageVersion("Seurat"), " CreateSeuratObject"),
  n_mice_nonrestored = sum(!mice$lkb1_restored),
  n_mice_restored = sum(mice$lkb1_restored),
  n_cells_qc = nrow(md),
  Cldn4_delta_mouse = grab("Cldn4")$delta_mean,
  Cldn4_welch_p_mouse = grab("Cldn4")$welch_p,
  Cldn4_mwu_p_mouse = grab("Cldn4")$mwu_p,
  frac_Cldn4_pos_delta_mouse = grab("frac_Cldn4_pos")$delta_mean,
  AT2_delta_mouse = grab("AT2")$delta_mean,
  AT2_welch_p_mouse = grab("AT2")$welch_p,
  IFN_A8_delta_mouse = grab("IFN_A8")$delta_mean,
  IFN_A8_welch_p_mouse = grab("IFN_A8")$welch_p,
  MHC_A8_delta_mouse = grab("MHC_A8")$delta_mean,
  MHC_A8_welch_p_mouse = grab("MHC_A8")$welch_p,
  TJ_A8_delta_mouse = grab("TJ_A8")$delta_mean,
  Stk11_delta_mouse = grab("Stk11")$delta_mean,
  Stk11_welch_p_mouse = grab("Stk11")$welch_p,
  stringsAsFactors = FALSE
)
write.table(one, file.path(TABLES, "one_row.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Reductions for DimPlot
message("PCA + UMAP")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
umap_ok <- TRUE
tryCatch(
  {
    obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  },
  error = function(e) {
    message("UMAP failed: ", conditionMessage(e))
    umap_ok <<- FALSE
  }
)

# Cldn4 high/low factor for plots (global rank tails)
ord <- order(md$Cldn4)
n_tail <- nrow(md) %/% 4L
gate <- rep("mid", nrow(md))
gate[ord[seq_len(n_tail)]] <- "Cldn4_Q1"
gate[ord[(nrow(md) - n_tail + 1):nrow(md)]] <- "Cldn4_Q4"
obj$Cldn4_gate <- gate[match(colnames(obj), rownames(md))]
obj$Cldn4_pos <- ifelse(obj$Cldn4_count > 0, "Cldn4+", "Cldn4-")
obj$cohort <- factor(obj$cohort, levels = c("NonRestored", "Restored"))
obj$mouse <- factor(obj$mouse, levels = c("CM0875", "ZR1932", "ZR1966", "CM0879", "CM0884", "ZR1969"))

Idents(obj) <- "cohort"
p_vln_cldn4 <- VlnPlot(obj, features = "Cldn4", group.by = "cohort", pt.size = 0) +
  ggtitle("Cldn4 (log-normalized) by Lkb1 restore")
save_plot(p_vln_cldn4, file.path(FIGS, "vln_cldn4_by_restore"), 5.5, 4.2)

p_vln_mouse <- VlnPlot(obj, features = "Cldn4", group.by = "mouse", pt.size = 0) +
  ggtitle("Cldn4 by mouse (honest n = 6)")
save_plot(p_vln_mouse, file.path(FIGS, "vln_cldn4_by_mouse"), 7.5, 4.2)

p_vln_mods <- VlnPlot(
  obj,
  features = c("AT2_AMS", "IFN_A8_AMS", "MHC_A8_AMS", "TJ_A8_AMS"),
  group.by = "cohort",
  pt.size = 0,
  ncol = 2
) + plot_annotation(title = "Seurat AddModuleScore by Lkb1 restore")
save_plot(p_vln_mods, file.path(FIGS, "vln_modules_by_restore"), 8.5, 6.5)

p_vln_hl <- VlnPlot(
  obj,
  features = c("AT2_AMS", "IFN_A8_AMS", "TJ_A8_AMS"),
  group.by = "Cldn4_gate",
  pt.size = 0,
  ncol = 3
) + plot_annotation(title = "Modules in Cldn4 Q1 vs Q4 (rank tails; zeros in Q1)")
save_plot(p_vln_hl, file.path(FIGS, "vln_modules_cldn4_highlow"), 9.5, 4.2)

dot_genes <- resolve_genes(obj, c(
  "Cldn4", "Stk11", "Sftpc", "Sftpb", "Lamp3", "Lyz2", "Nkx2-1",
  "Stat1", "Irf7", "B2m", "H2-K1", "H2-D1", "Cldn3", "Ocln", "Tjp1"
))
p_dot <- DotPlot(obj, features = dot_genes, group.by = "cohort") +
  RotatedAxis() +
  ggtitle("Cldn4 / Stk11 / AT2 / IFN / MHC-I / TJ genes")
save_plot(p_dot, file.path(FIGS, "dotplot_restore"), 8.5, 3.8)

p_pca_coh <- DimPlot(obj, reduction = "pca", group.by = "cohort") +
  ggtitle("PCA: Lkb1 restore cohort")
p_pca_mouse <- DimPlot(obj, reduction = "pca", group.by = "mouse") +
  ggtitle("PCA: mouse (honest n)")
save_plot(p_pca_coh + p_pca_mouse, file.path(FIGS, "dimplot_pca"), 10, 4.4)

if (umap_ok && "umap" %in% Reductions(obj)) {
  p_umap_coh <- DimPlot(obj, reduction = "umap", group.by = "cohort") + ggtitle("UMAP: restore")
  p_umap_gate <- DimPlot(obj, reduction = "umap", group.by = "Cldn4_gate") + ggtitle("UMAP: Cldn4 Q1/Q4")
  p_umap_feat <- FeaturePlot(obj, features = "Cldn4", reduction = "umap") + ggtitle("UMAP: Cldn4")
  save_plot(p_umap_coh + p_umap_gate, file.path(FIGS, "dimplot_umap"), 10, 4.4)
  save_plot(p_umap_feat, file.path(FIGS, "featureplot_cldn4_umap"), 5.2, 4.4)
}

# persist a slim RDS proof that CreateSeuratObject ran
meta_out <- FetchData(obj, vars = c(
  "mouse", "cohort", "lkb1_restored", "Cldn4", "Cldn4_count", "Stk11",
  "AT2_mean", "IFN_A8_mean", "MHC_A8_mean", "TJ_A8_mean",
  "AT2_AMS", "IFN_A8_AMS", "MHC_A8_AMS", "TJ_A8_AMS",
  "Cldn4_gate", "nCount_RNA", "nFeature_RNA", "percent.mt"
))
write.table(meta_out, file.path(TABLES, "cell_meta.tsv"), sep = "\t", quote = FALSE, col.names = NA)
saveRDS(
  list(
    seurat_class = class(obj),
    seurat_version = as.character(packageVersion("Seurat")),
    seuratobject_version = as.character(packageVersion("SeuratObject")),
    n_cells_raw = n_raw,
    n_cells_qc = ncol(obj),
    n_genes = nrow(obj),
    create_seurat_object = TRUE,
    reductions = Reductions(obj),
    mice = mice$mouse
  ),
  file.path(TABLES, "seurat_object_proof.rds")
)

summary <- list(
  accession = "GSE179502",
  pmid = 35228570,
  species = "mouse",
  compartment = "FACS-sorted neoplastic epithelium",
  engine = paste0("Seurat ", packageVersion("Seurat"), " CreateSeuratObject + AddModuleScore"),
  n_mice = nrow(mice),
  n_restored = sum(mice$lkb1_restored),
  n_nonrestored = sum(!mice$lkb1_restored),
  n_cells_raw = n_raw,
  n_cells_qc = ncol(obj),
  qc = list(min_genes = MIN_GENES, min_umi = MIN_UMI, max_mito = MAX_MITO),
  mice = mice,
  restore_mouse = lapply(
    c("Cldn4", "frac_Cldn4_pos", "Stk11", "AT2", "IFN_A8", "MHC_A8", "TJ_A8", "Tacstd2"),
    function(k) grab(k)
  ),
  gene_coverage = cov[, c("set", "n_present")],
  reductions = Reductions(obj),
  create_seurat_object = TRUE,
  notes = c(
    "Cldn4 only. Tacstd2 audit. No dual-high. No human. No GSE179501. No private 8 KL.",
    "Author GEO cohort Restored = restorable + tamoxifen (3 mice).",
    "NonRestored = non-restorable +/- tamoxifen/vehicle or restorable + vehicle (3 mice).",
    "MWU at 3 vs 3 cannot go below p=0.1.",
    "Module scores: *_mean = mean log1p CP10k of member genes; *_AMS = Seurat AddModuleScore."
  )
)
# fix restore_mouse names
names(summary$restore_mouse) <- c("Cldn4", "frac_Cldn4_pos", "Stk11", "AT2", "IFN_A8", "MHC_A8", "TJ_A8", "Tacstd2")
write(
  toJSON(summary, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 8),
  file.path(TABLES, "summary.json")
)

message("DONE cells_qc=", ncol(obj), " mice=", nrow(mice))
print(mice[, c("mouse", "cohort", "n_cells_qc", "frac_Cldn4_pos", "Cldn4", "Stk11", "AT2", "IFN_A8")])
print(contrasts[contrasts$unit == "mouse" & contrasts$feature %in% c("Cldn4", "frac_Cldn4_pos", "AT2", "IFN_A8", "Stk11"), ])
