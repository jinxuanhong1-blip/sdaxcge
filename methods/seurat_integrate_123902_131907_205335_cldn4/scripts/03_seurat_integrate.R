#!/usr/bin/env Rscript
# Seurat v5 + Harmony integration of the triple that differs:
# GSE123902 + GSE131907 + GSE205335. CLDN4-only.
# Honest unit = patient. Dataset covariate. Thesis already correct; ADDITIVE.
# No dual-high. Not +GSE148071. Not +GSE189357. No Python-only primary.

.libPaths(c("/tmp/r_libs", .libPaths()))
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat is not installed. Stop. Do not fall back to a Python-only primary.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

options(warn = 1)
set.seed(1)

ROOT <- "/workspace/methods/seurat_integrate_123902_131907_205335_cldn4"
GEO <- "/tmp/triple_geo"
DATA <- file.path(ROOT, "data")
EXT <- file.path(ROOT, "extract")
RES <- file.path(ROOT, "results")
FIG <- file.path(RES, "figures")
TAB <- file.path(RES, "tables")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

MHC_II <- c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
            "HLA-DQA1", "HLA-DQB1", "CD74", "CIITA")
COHORTS <- c("GSE123902", "GSE131907", "GSE205335")

spearman_df <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 5) {
    return(data.frame(n = n, rho = NA_real_, p = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  z <- atanh(max(min(rho, 0.999999), -0.999999))
  se <- 1 / sqrt(n - 3)
  data.frame(
    n = n, rho = rho, p = ct$p.value,
    lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se)
  )
}

mw_df <- function(x, g_high) {
  ok <- is.finite(x) & !is.na(g_high)
  x <- x[ok]; g_high <- g_high[ok]
  n1 <- sum(!g_high); n2 <- sum(g_high)
  if (n1 < 3 || n2 < 3) {
    return(data.frame(n_low = n1, n_high = n2, r = NA_real_, p = NA_real_))
  }
  wt <- suppressWarnings(wilcox.test(x[g_high], x[!g_high], exact = FALSE))
  u <- as.numeric(wt$statistic)
  r_rb <- (2 * u) / (n1 * n2) - 1
  data.frame(n_low = n1, n_high = n2, r = r_rb, p = wt$p.value)
}

fisher_z_dl <- function(rho, n) {
  keep <- is.finite(rho) & is.finite(n) & n > 3 & abs(rho) < 1
  rho <- rho[keep]; n <- n[keep]
  k <- length(rho)
  if (k == 0) {
    return(data.frame(k = 0, N = 0, rho = NA_real_, p = NA_real_, I2 = NA_real_))
  }
  if (k == 1) {
    z <- atanh(rho)
    se <- 1 / sqrt(n - 3)
    p <- 2 * pnorm(-abs(z / se))
    return(data.frame(k = 1, N = n, rho = rho, p = p, I2 = 0))
  }
  z <- atanh(rho)
  w <- n - 3
  zbar <- sum(w * z) / sum(w)
  q <- sum(w * (z - zbar)^2)
  df <- k - 1
  tau2 <- max(0, (q - df) / (sum(w) - sum(w^2) / sum(w)))
  wstar <- 1 / (1 / w + tau2)
  zdl <- sum(wstar * z) / sum(wstar)
  se <- sqrt(1 / sum(wstar))
  p <- 2 * pnorm(-abs(zdl / se))
  i2 <- if (q <= 0) 0 else max(0, (q - df) / q)
  data.frame(k = k, N = sum(n), rho = tanh(zdl), p = p, I2 = 100 * i2)
}

assign_q <- function(x) {
  rnk <- rank(x, ties.method = "average", na.last = "keep")
  as.character(cut(
    rnk,
    breaks = quantile(rnk, probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE),
    include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4")
  ))
}

read_keep_mtx <- function(ds) {
  keep <- read.delim(file.path(EXT, paste0("keep_", ds, ".tsv")), check.names = FALSE)
  mat <- readMM(file.path(EXT, ds, "matrix.mtx"))
  genes <- readLines(file.path(EXT, ds, "genes.tsv"))
  bc <- readLines(file.path(EXT, ds, "barcodes.tsv"))
  if (nrow(mat) != length(genes) || ncol(mat) != length(bc)) {
    stop(ds, " MTX dim mismatch: ", paste(dim(mat), collapse = "x"),
         " vs genes=", length(genes), " bc=", length(bc))
  }
  rownames(mat) <- make.unique(genes)
  colnames(mat) <- bc
  mat <- as(mat, "CsparseMatrix")
  meta <- keep
  rownames(meta) <- meta$barcode
  miss <- setdiff(bc, rownames(meta))
  if (length(miss)) {
    stop(ds, " barcodes missing from keep table: ", length(miss))
  }
  meta <- meta[bc, , drop = FALSE]
  obj <- CreateSeuratObject(counts = mat, meta.data = meta, project = ds)
  rm(mat)
  gc(verbose = FALSE)
  obj
}

# ---- gene sets ----
js <- jsonlite::fromJSON(file.path(DATA, "a8_sets.json"))
ifn <- unique(c(js$sets$HALLMARK_INTERFERON_ALPHA_RESPONSE,
                js$sets$HALLMARK_INTERFERON_GAMMA_RESPONSE))
ifn <- setdiff(ifn, "CLDN4")
mhc_i <- unique(js$sets$CUSTOM_MHC_I_ANTIGEN_PRESENTATION)
mhc <- unique(c(mhc_i, MHC_II))
mhc <- setdiff(mhc, "CLDN4")
mhc_i <- setdiff(mhc_i, "CLDN4")

# ---- locked patient tables (PR #459; T/NK not re-audited) ----
u123 <- read.delim(file.path(DATA, "GSE123902_marker_units.tsv"), check.names = FALSE)
u123 <- u123[u123$tissue %in% c("PRIMARY", "METASTASIS") &
               tolower(as.character(u123$eligible)) == "true", ]
u123 <- u123[order(u123$patient, u123$tissue != "PRIMARY"), ]
u123 <- u123[!duplicated(u123$patient), ]
# locked %pos is a fraction on this table
u123$locked_pct <- 100 * as.numeric(u123$mal_CLDN4_pct)

s131 <- read.delim(file.path(DATA, "GSE131907_samples.tsv"), check.names = FALSE)
p205 <- read.delim(file.path(DATA, "GSE205335_patients.tsv"), check.names = FALSE)

lock <- rbind(
  data.frame(
    dataset = "GSE123902",
    unit_id = as.character(u123$patient),
    locked_frac_tnk = as.numeric(u123$frac_tnk),
    locked_mal_cldn4_pct = u123$locked_pct,
    locked_n_mal = as.integer(u123$n_malignant),
    locked_n_tnk = as.integer(u123$n_tnk),
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset = "GSE131907",
    unit_id = as.character(s131$sample),
    locked_frac_tnk = as.numeric(s131$frac_tnk),
    locked_mal_cldn4_pct = as.numeric(s131$mal_CLDN4_pct),
    locked_n_mal = as.integer(s131$n_malignant),
    locked_n_tnk = as.integer(s131$n_tnk),
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset = "GSE205335",
    unit_id = as.character(p205$patient),
    locked_frac_tnk = as.numeric(p205$frac_tnk),
    locked_mal_cldn4_pct = as.numeric(p205$mal_CLDN4_pct_pos),
    locked_n_mal = as.integer(p205$n_malignant),
    locked_n_tnk = as.integer(p205$n_tnk),
    stringsAsFactors = FALSE
  )
)

# ---- load three objects ----
cat("[load] GSE123902 MTX\n")
obj123 <- read_keep_mtx("GSE123902")
cat("[load] GSE131907 MTX\n")
obj131 <- read_keep_mtx("GSE131907")
cat("[load] GSE205335 RDS subset\n")
keep205 <- read.delim(file.path(EXT, "keep_GSE205335.tsv"), check.names = FALSE)
rds_path <- file.path(GEO, "GSE205335_Lung_IO_UMI_matrix.rds")
if (!file.exists(rds_path)) {
  stop("missing decompressed GSE205335 RDS; run scripts/00_download.sh")
}
umi205 <- readRDS(rds_path)
cat("[load] GSE205335 dim", paste(dim(umi205), collapse = "x"), "\n")
missing <- setdiff(keep205$barcode, colnames(umi205))
if (length(missing)) {
  cat("[warn] ", length(missing), " GSE205335 barcodes missing from matrix\n", sep = "")
  keep205 <- subset(keep205, barcode %in% colnames(umi205))
}
mat205 <- umi205[, keep205$barcode, drop = FALSE]
rm(umi205)
gc(verbose = FALSE)
meta205 <- keep205
rownames(meta205) <- meta205$barcode
obj205 <- CreateSeuratObject(counts = mat205, meta.data = meta205, project = "GSE205335")
rm(mat205)
gc(verbose = FALSE)

# ---- merge + Harmony ----
cat("[seurat] merge triple\n")
obj <- merge(obj123, y = c(obj131, obj205))
rm(obj123, obj131, obj205)
gc(verbose = FALSE)
obj$dataset <- as.character(obj$dataset)
obj$compartment <- as.character(obj$compartment)
obj$unit_id <- as.character(obj$unit_id)
obj$patient_id <- as.character(obj$patient_id)
DefaultAssay(obj) <- "RNA"
# merge() already splits Assay5 layers by project (= dataset). Join+resplit
# only if a single counts layer remains.
ly <- Layers(obj[["RNA"]])
cat("[seurat] RNA layers after merge:", paste(ly, collapse = ", "), "\n")
if (length(ly) <= 1) {
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$dataset)
} else {
  cat("[seurat] keeping pre-split dataset layers\n")
}

cat("[seurat] normalize / HVG / PCA\n")
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

integ_method <- "HarmonyIntegration"
integ_reduction <- "harmony"
ok_int <- FALSE
if (exists("HarmonyIntegration", mode = "function") &&
    requireNamespace("harmony", quietly = TRUE)) {
  cat("[seurat] IntegrateLayers HarmonyIntegration (dataset layers)\n")
  ok_int <- tryCatch({
    obj <<- IntegrateLayers(
      object = obj,
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "harmony",
      verbose = FALSE
    )
    TRUE
  }, error = function(e) {
    cat("[seurat] HarmonyIntegration failed:", conditionMessage(e), "\n")
    FALSE
  })
}
if (!ok_int) {
  cat("[seurat] fallback RunHarmony\n")
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Harmony is not available after Seurat load. Stop.")
  }
  obj <- JoinLayers(obj)
  obj <- harmony::RunHarmony(obj, group.by.vars = "dataset", reduction.use = "pca",
                             reduction.save = "harmony", verbose = FALSE)
  integ_method <- "RunHarmony"
  integ_reduction <- "harmony"
  ok_int <- TRUE
}

obj <- FindNeighbors(obj, reduction = integ_reduction, dims = 1:20, verbose = FALSE)
obj <- RunUMAP(obj, reduction = integ_reduction, dims = 1:20, verbose = FALSE)

cat("[seurat] JoinLayers + module scores\n")
if (length(Layers(obj[["RNA"]])) > 1) {
  obj <- JoinLayers(obj)
}
features_ifn <- intersect(ifn, rownames(obj))
features_mhc <- intersect(mhc, rownames(obj))
features_mhci <- intersect(mhc_i, rownames(obj))
if (length(features_ifn) < 10) stop("too few IFN genes in object")
if (length(features_mhc) < 5) stop("too few MHC genes in object")
obj <- AddModuleScore(
  obj,
  features = list(IFN = features_ifn, MHC = features_mhc, MHCI = features_mhci),
  name = "mod",
  ctrl = 50,
  seed = 1
)
obj$IFN_score <- obj$mod1
obj$MHC_score <- obj$mod2
obj$MHCI_score <- obj$mod3

if (!("CLDN4" %in% rownames(obj))) stop("CLDN4 absent after merge")
cldn4_counts <- as.numeric(LayerData(obj, layer = "counts")["CLDN4", ])
cldn4_data <- as.numeric(LayerData(obj, layer = "data")["CLDN4", ])
obj$CLDN4_counts <- cldn4_counts
obj$CLDN4_log1p <- cldn4_data
obj$CLDN4_pos <- as.integer(cldn4_counts > 0)

# ---- patient-level scores ----
md <- slot(obj, "meta.data")
units <- unique(md[, c("dataset", "unit_id", "patient_id")])
score_rows <- vector("list", nrow(units))
for (i in seq_len(nrow(units))) {
  ds <- units$dataset[i]
  uid <- units$unit_id[i]
  sub <- md[md$dataset == ds & md$unit_id == uid, , drop = FALSE]
  mal <- sub[sub$compartment == "malignant", , drop = FALSE]
  tnk <- sub[sub$compartment == "TNK", , drop = FALSE]
  score_rows[[i]] <- data.frame(
    dataset = ds,
    unit_id = uid,
    patient_id = units$patient_id[i],
    n_mal_in_object = nrow(mal),
    n_tnk_in_object = nrow(tnk),
    mal_CLDN4_mean = if (nrow(mal)) mean(mal$CLDN4_log1p, na.rm = TRUE) else NA_real_,
    mal_CLDN4_pct = if (nrow(mal)) 100 * mean(mal$CLDN4_pos, na.rm = TRUE) else NA_real_,
    mal_IFN = if (nrow(mal)) mean(mal$IFN_score, na.rm = TRUE) else NA_real_,
    mal_MHC = if (nrow(mal)) mean(mal$MHC_score, na.rm = TRUE) else NA_real_,
    mal_MHCI = if (nrow(mal)) mean(mal$MHCI_score, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}
scores <- do.call(rbind, score_rows)
scores <- merge(scores, lock, by = c("dataset", "unit_id"), all.x = TRUE)
scores <- scores[scores$n_mal_in_object >= 10 & is.finite(scores$locked_frac_tnk), ]
scores$dataset <- factor(scores$dataset, levels = COHORTS)
scores$q_seurat <- NA_character_
for (ds in levels(scores$dataset)) {
  idx <- scores$dataset == ds
  if (sum(idx) >= 4) scores$q_seurat[idx] <- assign_q(scores$mal_CLDN4_pct[idx])
}
scores$cldn4_z <- NA_real_
for (ds in levels(scores$dataset)) {
  idx <- scores$dataset == ds
  scores$cldn4_z[idx] <- as.numeric(scale(scores$mal_CLDN4_pct[idx]))
}
write.table(scores, file.path(TAB, "patient_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- tests: within-dataset, DL, residual (dataset covariate), OLS ----
contrasts <- list(
  cldn4_vs_tnk = list(x = "mal_CLDN4_pct", y = "locked_frac_tnk"),
  cldn4_vs_ifn = list(x = "mal_CLDN4_pct", y = "mal_IFN"),
  cldn4_vs_mhc = list(x = "mal_CLDN4_pct", y = "mal_MHC"),
  cldn4_vs_mhci = list(x = "mal_CLDN4_pct", y = "mal_MHCI"),
  locked_cldn4_vs_tnk = list(x = "locked_mal_cldn4_pct", y = "locked_frac_tnk")
)

test_rows <- list()
ols_rows <- list()
for (nm in names(contrasts)) {
  xcol <- contrasts[[nm]]$x
  ycol <- contrasts[[nm]]$y
  per <- list()
  for (ds in COHORTS) {
    sub <- scores[as.character(scores$dataset) == ds, ]
    sp <- spearman_df(sub[[xcol]], sub[[ycol]])
    qh <- sub$q_seurat == "Q4"
    ql <- sub$q_seurat == "Q1"
    mw <- mw_df(sub[[ycol]][ql | qh], qh[ql | qh])
    per[[ds]] <- data.frame(
      dataset = ds, contrast = nm, model = "within_cohort_spearman",
      n = sp$n, rho = sp$rho, p = sp$p, lo = sp$lo, hi = sp$hi,
      n_low = mw$n_low, n_high = mw$n_high, r = mw$r, p_q4q1 = mw$p,
      I2 = NA_real_, k = 1, beta = NA_real_,
      stringsAsFactors = FALSE
    )
    test_rows[[length(test_rows) + 1]] <- per[[ds]]
  }
  dl <- fisher_z_dl(vapply(per, function(z) z$rho, 1), vapply(per, function(z) z$n, 1))
  stack <- scores[scores$q_seurat %in% c("Q1", "Q4"), ]
  mw_stack <- mw_df(stack[[ycol]], stack$q_seurat == "Q4")
  test_rows[[length(test_rows) + 1]] <- data.frame(
    dataset = "combined_DL", contrast = nm, model = "dersimonian_laird",
    n = dl$N, rho = dl$rho, p = dl$p, lo = NA_real_, hi = NA_real_,
    n_low = mw_stack$n_low, n_high = mw_stack$n_high, r = mw_stack$r,
    p_q4q1 = mw_stack$p, I2 = dl$I2, k = dl$k, beta = NA_real_,
    stringsAsFactors = FALSE
  )
  # residual Spearman: both axes residualized on dataset
  ok <- is.finite(scores[[xcol]]) & is.finite(scores[[ycol]])
  rx <- residuals(lm(scores[[xcol]][ok] ~ scores$dataset[ok]))
  ry <- residuals(lm(scores[[ycol]][ok] ~ scores$dataset[ok]))
  sp_res <- spearman_df(rx, ry)
  test_rows[[length(test_rows) + 1]] <- data.frame(
    dataset = "dataset_residual", contrast = nm, model = "spearman_residual_dataset",
    n = sp_res$n, rho = sp_res$rho, p = sp_res$p, lo = sp_res$lo, hi = sp_res$hi,
    n_low = mw_stack$n_low, n_high = mw_stack$n_high, r = mw_stack$r,
    p_q4q1 = mw_stack$p, I2 = NA_real_, k = 3, beta = NA_real_,
    stringsAsFactors = FALSE
  )
  # OLS with dataset covariate; predictor z-scored within dataset
  xz <- scores$cldn4_z
  if (xcol != "mal_CLDN4_pct") {
    xz <- rep(NA_real_, nrow(scores))
    for (ds in levels(scores$dataset)) {
      idx <- scores$dataset == ds
      xz[idx] <- as.numeric(scale(scores[[xcol]][idx]))
    }
  }
  ok2 <- ok & is.finite(xz)
  fit <- lm(scores[[ycol]][ok2] ~ xz[ok2] + scores$dataset[ok2])
  sm <- summary(fit)
  beta <- unname(coef(fit)[2])
  p_beta <- unname(sm$coefficients[2, 4])
  test_rows[[length(test_rows) + 1]] <- data.frame(
    dataset = "dataset_covariate_OLS", contrast = nm, model = "ols_cldn4z_plus_dataset",
    n = sum(ok2), rho = NA_real_, p = p_beta, lo = NA_real_, hi = NA_real_,
    n_low = NA_integer_, n_high = NA_integer_, r = NA_real_, p_q4q1 = NA_real_,
    I2 = NA_real_, k = 3, beta = beta,
    stringsAsFactors = FALSE
  )
  ols_rows[[nm]] <- data.frame(
    contrast = nm,
    n = sum(ok2),
    beta_cldn4z = beta,
    se = unname(sm$coefficients[2, 2]),
    t = unname(sm$coefficients[2, 3]),
    p = p_beta,
    stringsAsFactors = FALSE
  )
}

bind_rows <- function(dfs) {
  cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    miss <- setdiff(cols, names(d))
    for (m in miss) d[[m]] <- NA
    d[, cols, drop = FALSE]
  })
  do.call(rbind, dfs)
}
tests <- bind_rows(test_rows)
write.table(tests, file.path(TAB, "patient_level_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
ols <- do.call(rbind, ols_rows)
write.table(ols, file.path(TAB, "dataset_covariate_ols.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

n_tab <- data.frame(
  item = c(
    "GSE123902 locked tumor donors",
    "GSE131907 locked tumor samples (patient-equal on this gate)",
    "GSE205335 locked patients",
    "combined patients in Seurat scores",
    "cells in Harmony object",
    "malignant cells in object",
    "T/NK cells in object",
    "IFN genes scored",
    "MHC genes scored"
  ),
  n = c(
    sum(as.character(scores$dataset) == "GSE123902"),
    sum(as.character(scores$dataset) == "GSE131907"),
    sum(as.character(scores$dataset) == "GSE205335"),
    nrow(scores),
    ncol(obj),
    sum(obj$compartment == "malignant"),
    sum(obj$compartment == "TNK"),
    length(features_ifn),
    length(features_mhc)
  ),
  note = c(
    "Laughney donor; marker-malignant; PRIMARY preferred",
    "Kim GEO Sample; unique patient numbers on the gated set",
    "Hu patient after collapsing tumor GSMs",
    "patient is the unit; cells are not n",
    "capped <=80 mal + <=80 T/NK per patient",
    "123902 marker-mal; 131907/205335 author-mal",
    "T/NK fraction uses locked full-sample frac_tnk",
    "Hallmark IFNa U IFNg, CLDN4 excluded",
    "custom MHC-I + MHC-II, CLDN4 excluded"
  )
)
write.table(n_tab, file.path(TAB, "n_honest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- figures ----
umap <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap) <- c("UMAP_1", "UMAP_2")
umap <- cbind(umap, slot(obj, "meta.data"))
set.seed(1)
plot_n <- min(nrow(umap), 8000)
umap_p <- umap[sample.int(nrow(umap), plot_n), ]

ggsave_both <- function(p, name, w = 6.2, h = 5.0) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 140)
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
}
theme_ok <- theme_bw(base_size = 11) + theme(legend.position = "bottom")

ggsave_both(
  ggplot(umap_p, aes(UMAP_1, UMAP_2, color = dataset)) +
    geom_point(size = 0.25, alpha = 0.7) + theme_ok +
    ggtitle("Seurat Harmony UMAP — dataset"),
  "umap_dataset"
)
ggsave_both(
  ggplot(umap_p, aes(UMAP_1, UMAP_2, color = compartment)) +
    geom_point(size = 0.25, alpha = 0.7) + theme_ok +
    ggtitle("Seurat Harmony UMAP — compartment"),
  "umap_compartment"
)
ggsave_both(
  ggplot(umap_p, aes(UMAP_1, UMAP_2, color = CLDN4_log1p)) +
    geom_point(size = 0.25, alpha = 0.8) +
    scale_color_viridis_c(option = "C") + theme_ok +
    ggtitle("Seurat Harmony UMAP — CLDN4 (log-norm)"),
  "umap_cldn4"
)
ggsave_both(
  ggplot(umap_p[umap_p$compartment == "malignant", ],
         aes(UMAP_1, UMAP_2, color = IFN_score)) +
    geom_point(size = 0.35, alpha = 0.85) +
    scale_color_viridis_c(option = "D") + theme_ok +
    ggtitle("Malignant cells — IFN module (CLDN4 excluded)"),
  "umap_mal_ifn"
)
ggsave_both(
  ggplot(scores, aes(mal_CLDN4_pct, locked_frac_tnk, color = dataset)) +
    geom_point(size = 2.4) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    theme_ok +
    xlab("Malignant CLDN4 %pos (Seurat object)") +
    ylab("Author / locked T/NK fraction (full sample)") +
    ggtitle("Patient unit: CLDN4 vs T/NK"),
  "scatter_cldn4_tnk"
)
ggsave_both(
  ggplot(scores, aes(mal_CLDN4_pct, mal_IFN, color = dataset)) +
    geom_point(size = 2.4) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    theme_ok +
    xlab("Malignant CLDN4 %pos") +
    ylab("Malignant IFN module (AddModuleScore)") +
    ggtitle("Patient unit: CLDN4 vs malignant IFN"),
  "scatter_cldn4_ifn"
)
ggsave_both(
  ggplot(scores, aes(mal_CLDN4_pct, mal_MHC, color = dataset)) +
    geom_point(size = 2.4) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
    theme_ok +
    xlab("Malignant CLDN4 %pos") +
    ylab("Malignant MHC module (AddModuleScore)") +
    ggtitle("Patient unit: CLDN4 vs malignant MHC"),
  "scatter_cldn4_mhc"
)
qdf <- scores[scores$q_seurat %in% c("Q1", "Q4"), ]
qdf$q_seurat <- factor(qdf$q_seurat, levels = c("Q1", "Q4"))
ggsave_both(
  ggplot(qdf, aes(q_seurat, locked_frac_tnk, fill = dataset)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               size = 1.6) +
    theme_ok + ylab("T/NK fraction") + xlab("Within-dataset CLDN4 %pos quartile") +
    ggtitle("Q4 vs Q1 T/NK (unit = patient)"),
  "box_q4q1_tnk"
)
ggsave_both(
  ggplot(qdf, aes(q_seurat, mal_IFN, fill = dataset)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               size = 1.6) +
    theme_ok + ylab("Malignant IFN") + xlab("Within-dataset CLDN4 %pos quartile") +
    ggtitle("Q4 vs Q1 malignant IFN"),
  "box_q4q1_ifn"
)
ggsave_both(
  ggplot(qdf, aes(q_seurat, mal_MHC, fill = dataset)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               size = 1.6) +
    theme_ok + ylab("Malignant MHC") + xlab("Within-dataset CLDN4 %pos quartile") +
    ggtitle("Q4 vs Q1 malignant MHC"),
  "box_q4q1_mhc"
)
ggsave_both(
  ggplot(n_tab, aes(x = item, y = n)) +
    geom_col(fill = "#4C78A8") + coord_flip() + theme_ok +
    ggtitle("Honest n (cells are not the test unit)"),
  "n_honest", w = 8.5, h = 4.2
)

summary <- list(
  seurat_version = as.character(packageVersion("Seurat")),
  seuratobject_version = as.character(packageVersion("SeuratObject")),
  harmony_version = if (requireNamespace("harmony", quietly = TRUE))
    as.character(packageVersion("harmony")) else NA_character_,
  r_version = R.version.string,
  integration_method = integ_method,
  integration_reduction = integ_reduction,
  n_cells = ncol(obj),
  n_units = nrow(scores),
  n_gse123902 = sum(as.character(scores$dataset) == "GSE123902"),
  n_gse131907 = sum(as.character(scores$dataset) == "GSE131907"),
  n_gse205335 = sum(as.character(scores$dataset) == "GSE205335"),
  n_ifn_genes = length(features_ifn),
  n_mhc_genes = length(features_mhc),
  n_mhci_genes = length(features_mhci),
  tests = tests,
  ols = ols
)
jsonlite::write_json(summary, file.path(RES, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
saveRDS(slot(obj, "meta.data"), file.path(RES, "cell_metadata.rds"))
write.table(slot(obj, "meta.data"), file.path(TAB, "cell_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)
writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo.txt"))

cat("[done] units=", nrow(scores), " cells=", ncol(obj),
    " method=", integ_method, "\n", sep = "")
print(tests[tests$contrast %in% c("cldn4_vs_tnk", "cldn4_vs_ifn", "cldn4_vs_mhc") &
              tests$model %in% c("within_cohort_spearman", "dersimonian_laird",
                                 "spearman_residual_dataset", "ols_cldn4z_plus_dataset"),
            c("dataset", "contrast", "model", "n", "rho", "p", "beta")])
