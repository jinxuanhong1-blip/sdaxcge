#!/usr/bin/env Rscript
# Seurat v5 RPCA integration of the winning pair (GSE131907 + GSE205335).
# CLDN4-only. Honest unit = sample (GSE131907) / patient (GSE205335).
# Primary: CLDN4 vs T/NK; malignant IFN/MHC. Thesis already correct; ADDITIVE.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

options(warn = 1)
set.seed(1)

ROOT <- "/workspace/methods/seurat_winpair_131907_205335_cldn4"
GEO <- "/tmp/winpair_geo"
DATA <- file.path(ROOT, "data")
EXT <- file.path(ROOT, "extract")
RES <- file.path(ROOT, "results")
FIG <- file.path(RES, "figures")
TAB <- file.path(RES, "tables")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

MHC_II <- c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
            "HLA-DQA1", "HLA-DQB1", "CD74", "CIITA")

spearman_df <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 5) {
    return(data.frame(n = n, rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(n = n, rho = unname(ct$estimate), p = ct$p.value)
}

mw_df <- function(x, g_high) {
  ok <- is.finite(x) & !is.na(g_high)
  x <- x[ok]; g_high <- g_high[ok]
  n1 <- sum(!g_high); n2 <- sum(g_high)
  if (n1 < 3 || n2 < 3) {
    return(data.frame(n_low = n1, n_high = n2, r = NA_real_, p = NA_real_))
  }
  wt <- suppressWarnings(wilcox.test(x[g_high], x[!g_high], exact = FALSE))
  # rank-biserial: 2U/(n1 n2) - 1, high minus low
  rnk <- rank(x)
  r <- (mean(rnk[g_high]) - mean(rnk[!g_high])) / (n1 + n2) * 2
  # simpler signed rank-biserial from U
  u <- as.numeric(wt$statistic)
  r_rb <- 1 - (2 * u) / (n1 * n2)
  # wilcox.test statistic is W for first sample = high group
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
  as.character(cut(rnk, breaks = quantile(rnk, probs = seq(0, 1, 0.25), na.rm = TRUE),
                   include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4")))
}

# ---- gene sets ----
js <- jsonlite::fromJSON(file.path(DATA, "a8_sets.json"))
ifn <- unique(c(js$sets$HALLMARK_INTERFERON_ALPHA_RESPONSE,
                js$sets$HALLMARK_INTERFERON_GAMMA_RESPONSE))
ifn <- setdiff(ifn, "CLDN4")
mhc <- unique(c(js$sets$CUSTOM_MHC_I_ANTIGEN_PRESENTATION, MHC_II))
mhc <- setdiff(mhc, "CLDN4")

# ---- locked unit tables ----
s131 <- read.delim(file.path(DATA, "GSE131907_samples.tsv"), check.names = FALSE)
p205 <- read.delim(file.path(DATA, "GSE205335_patients.tsv"), check.names = FALSE)
lock <- rbind(
  data.frame(
    dataset = "GSE131907",
    unit_id = s131$sample,
    locked_frac_tnk = s131$frac_tnk,
    locked_mal_cldn4_pct = s131$mal_CLDN4_pct,
    locked_n_mal = s131$n_malignant,
    locked_n_tnk = s131$n_tnk,
    stringsAsFactors = FALSE
  ),
  data.frame(
    dataset = "GSE205335",
    unit_id = p205$patient,
    locked_frac_tnk = p205$frac_tnk,
    locked_mal_cldn4_pct = p205$mal_CLDN4_pct_pos,
    locked_n_mal = p205$n_malignant,
    locked_n_tnk = p205$n_tnk,
    stringsAsFactors = FALSE
  )
)

# ---- GSE131907 MTX ----
cat("[load] GSE131907 MTX\n")
keep131 <- read.delim(file.path(EXT, "keep_GSE131907.tsv"), check.names = FALSE)
mat131 <- readMM(file.path(EXT, "GSE131907", "matrix.mtx"))
genes131 <- readLines(file.path(EXT, "GSE131907", "genes.tsv"))
bc131 <- readLines(file.path(EXT, "GSE131907", "barcodes.tsv"))
if (nrow(mat131) != length(genes131) || ncol(mat131) != length(bc131)) {
  stop("GSE131907 MTX dim mismatch")
}
rownames(mat131) <- make.unique(genes131)
colnames(mat131) <- bc131
mat131 <- as(mat131, "CsparseMatrix")
meta131 <- keep131
rownames(meta131) <- meta131$barcode
meta131 <- meta131[bc131, , drop = FALSE]
obj131 <- CreateSeuratObject(counts = mat131, meta.data = meta131, project = "GSE131907")
rm(mat131); gc()

# ---- GSE205335 RDS subset ----
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
rm(umi205); gc()
meta205 <- keep205
rownames(meta205) <- meta205$barcode
obj205 <- CreateSeuratObject(counts = mat205, meta.data = meta205, project = "GSE205335")
rm(mat205); gc()

# ---- merge + integrate ----
cat("[seurat] merge\n")
obj <- merge(obj131, obj205)
rm(obj131, obj205); gc()
obj$dataset <- as.character(obj$dataset)
obj$compartment <- as.character(obj$compartment)
obj$unit_id <- as.character(obj$unit_id)
DefaultAssay(obj) <- "RNA"
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$dataset)

cat("[seurat] normalize / HVG / PCA\n")
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

cat("[seurat] IntegrateLayers RPCA\n")
integ_method <- "RPCAIntegration"
integ_reduction <- "integrated.rpca"
ok_int <- tryCatch({
  obj <<- IntegrateLayers(
    object = obj,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    verbose = FALSE
  )
  TRUE
}, error = function(e) {
  cat("[seurat] RPCA failed:", conditionMessage(e), "\n")
  FALSE
})
if (!ok_int) {
  cat("[seurat] fallback CCAIntegration\n")
  integ_method <- "CCAIntegration"
  integ_reduction <- "integrated.cca"
  obj <- IntegrateLayers(
    object = obj,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE
  )
}

obj <- FindNeighbors(obj, reduction = integ_reduction, dims = 1:20, verbose = FALSE)
obj <- RunUMAP(obj, reduction = integ_reduction, dims = 1:20, verbose = FALSE)

cat("[seurat] JoinLayers + module scores\n")
obj <- JoinLayers(obj)
features_ifn <- list(IFN = intersect(ifn, rownames(obj)))
features_mhc <- list(MHC = intersect(mhc, rownames(obj)))
if (length(features_ifn$IFN) < 10) stop("too few IFN genes in object")
if (length(features_mhc$MHC) < 5) stop("too few MHC genes in object")
obj <- AddModuleScore(obj, features = features_ifn, name = "IFN", ctrl = 50, seed = 1)
obj <- AddModuleScore(obj, features = features_mhc, name = "MHC", ctrl = 50, seed = 1)
# AddModuleScore appends 1: IFN1, MHC2 if run separately; run together:
# re-score as one call for stable names
obj <- AddModuleScore(
  obj,
  features = list(IFN = features_ifn$IFN, MHC = features_mhc$MHC),
  name = "mod",
  ctrl = 50,
  seed = 1
)
# mod1 = IFN, mod2 = MHC
obj$IFN_score <- obj$mod1
obj$MHC_score <- obj$mod2

if (!("CLDN4" %in% rownames(obj))) stop("CLDN4 absent after merge")
cldn4_counts <- as.numeric(LayerData(obj, layer = "counts")["CLDN4", ])
cldn4_data <- as.numeric(LayerData(obj, layer = "data")["CLDN4", ])
obj$CLDN4_counts <- cldn4_counts
obj$CLDN4_log1p <- cldn4_data
obj$CLDN4_pos <- as.integer(cldn4_counts > 0)

# ---- unit-level scores (malignant cells) ----
md <- slot(obj, "meta.data")
units <- unique(md[, c("dataset", "unit_id")])
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
    n_mal_in_object = nrow(mal),
    n_tnk_in_object = nrow(tnk),
    mal_CLDN4_mean = if (nrow(mal)) mean(mal$CLDN4_log1p, na.rm = TRUE) else NA_real_,
    mal_CLDN4_pct = if (nrow(mal)) 100 * mean(mal$CLDN4_pos, na.rm = TRUE) else NA_real_,
    mal_IFN = if (nrow(mal)) mean(mal$IFN_score, na.rm = TRUE) else NA_real_,
    mal_MHC = if (nrow(mal)) mean(mal$MHC_score, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}
scores <- do.call(rbind, score_rows)
scores <- merge(scores, lock, by = c("dataset", "unit_id"), all.x = TRUE)
# primary unit gate: locked n_mal>=20 (already), plus malignant cells in object
scores <- scores[scores$n_mal_in_object >= 10 & is.finite(scores$locked_frac_tnk), ]
scores$q_seurat <- NA_character_
for (ds in unique(scores$dataset)) {
  idx <- scores$dataset == ds
  scores$q_seurat[idx] <- assign_q(scores$mal_CLDN4_pct[idx])
}
scores$q_locked <- NA_character_
for (ds in unique(scores$dataset)) {
  idx <- scores$dataset == ds
  scores$q_locked[idx] <- assign_q(scores$locked_mal_cldn4_pct[idx])
}
write.table(scores, file.path(TAB, "patient_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- tests ----
contrasts <- list(
  cldn4_vs_tnk = list(x = "mal_CLDN4_pct", y = "locked_frac_tnk"),
  cldn4_vs_ifn = list(x = "mal_CLDN4_pct", y = "mal_IFN"),
  cldn4_vs_mhc = list(x = "mal_CLDN4_pct", y = "mal_MHC"),
  locked_cldn4_vs_tnk = list(x = "locked_mal_cldn4_pct", y = "locked_frac_tnk")
)

test_rows <- list()
for (nm in names(contrasts)) {
  xcol <- contrasts[[nm]]$x
  ycol <- contrasts[[nm]]$y
  per <- list()
  for (ds in c("GSE131907", "GSE205335")) {
    sub <- scores[scores$dataset == ds, ]
    sp <- spearman_df(sub[[xcol]], sub[[ycol]])
    qh <- sub$q_seurat == "Q4"
    ql <- sub$q_seurat == "Q1"
    mw <- mw_df(sub[[ycol]][ql | qh], qh[ql | qh])
    names(mw) <- c("n_low", "n_high", "r", "p_q4q1")
    per[[ds]] <- data.frame(
      dataset = ds, contrast = nm,
      n = sp$n, rho = sp$rho, p = sp$p,
      n_low = mw$n_low, n_high = mw$n_high, r = mw$r, p_q4q1 = mw$p_q4q1,
      I2 = NA_real_, k = 1,
      stringsAsFactors = FALSE
    )
    test_rows[[length(test_rows) + 1]] <- per[[ds]]
  }
  dl <- fisher_z_dl(vapply(per, function(z) z$rho, 1), vapply(per, function(z) z$n, 1))
  # stacked Q4 vs Q1 (within-cohort quartiles)
  stack <- scores[scores$q_seurat %in% c("Q1", "Q4"), ]
  mw_stack <- mw_df(stack[[ycol]], stack$q_seurat == "Q4")
  test_rows[[length(test_rows) + 1]] <- data.frame(
    dataset = "combined_DL",
    contrast = nm,
    n = dl$N,
    rho = dl$rho,
    p = dl$p,
    n_low = mw_stack$n_low,
    n_high = mw_stack$n_high,
    r = mw_stack$r,
    p_q4q1 = mw_stack$p,
    I2 = dl$I2,
    k = dl$k,
    stringsAsFactors = FALSE
  )
}
# bind with NAs for extra cols
tests <- do.call(function(...) {
  dfs <- list(...)
  cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(d) {
    miss <- setdiff(cols, names(d))
    for (m in miss) d[[m]] <- NA
    d[, cols, drop = FALSE]
  })
  do.call(rbind, dfs)
}, test_rows)
# fix combined p column name
if (!("p_q4q1" %in% names(tests))) tests$p_q4q1 <- NA
write.table(tests, file.path(TAB, "patient_level_tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

n_tab <- data.frame(
  item = c(
    "GSE131907 locked samples n_mal>=20",
    "GSE131907 unique patient numbers on those samples",
    "GSE205335 locked patients n_mal>=20",
    "combined units in Seurat scores",
    "cells in integrated object",
    "malignant cells in object",
    "T/NK cells in object",
    "IFN genes scored",
    "MHC genes scored"
  ),
  n = c(
    length(unique(scores$unit_id[scores$dataset == "GSE131907"])),
    length(unique(keep131$patient_id[keep131$unit_id %in% scores$unit_id[scores$dataset == "GSE131907"]])),
    length(unique(scores$unit_id[scores$dataset == "GSE205335"])),
    nrow(scores),
    ncol(obj),
    sum(obj$compartment == "malignant"),
    sum(obj$compartment == "TNK"),
    length(features_ifn$IFN),
    length(features_mhc$MHC)
  ),
  note = c(
    "GEO Sample is the locked unit (PR #320)",
    "same n as samples on the gated tumor set",
    "patient after collapsing tumor GSMs",
    "sample + patient; cells are not n",
    "capped <=120 mal + <=120 T/NK per unit",
    "author malignant labels",
    "author T/NK labels; fraction uses locked full-sample frac_tnk",
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
# subsample for plotting
set.seed(1)
plot_n <- min(nrow(umap), 8000)
umap_p <- umap[sample.int(nrow(umap), plot_n), ]

ggsave_both <- function(p, name, w = 6.2, h = 5.0) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 140)
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
}

theme_ok <- theme_bw(base_size = 11) + theme(legend.position = "bottom")

p1 <- ggplot(umap_p, aes(UMAP_1, UMAP_2, color = dataset)) +
  geom_point(size = 0.25, alpha = 0.7) + theme_ok +
  ggtitle("Seurat RPCA UMAP — dataset")
ggsave_both(p1, "umap_dataset")

p2 <- ggplot(umap_p, aes(UMAP_1, UMAP_2, color = compartment)) +
  geom_point(size = 0.25, alpha = 0.7) + theme_ok +
  ggtitle("Seurat RPCA UMAP — author compartment")
ggsave_both(p2, "umap_compartment")

p3 <- ggplot(umap_p, aes(UMAP_1, UMAP_2, color = CLDN4_log1p)) +
  geom_point(size = 0.25, alpha = 0.8) +
  scale_color_viridis_c(option = "C") + theme_ok +
  ggtitle("Seurat RPCA UMAP — CLDN4 (log-norm)")
ggsave_both(p3, "umap_cldn4")

p4 <- ggplot(umap_p[umap_p$compartment == "malignant", ],
             aes(UMAP_1, UMAP_2, color = IFN_score)) +
  geom_point(size = 0.35, alpha = 0.85) +
  scale_color_viridis_c(option = "D") + theme_ok +
  ggtitle("Malignant cells — IFN module (CLDN4 excluded)")
ggsave_both(p4, "umap_mal_ifn")

p5 <- ggplot(scores, aes(mal_CLDN4_pct, locked_frac_tnk, color = dataset)) +
  geom_point(size = 2.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  theme_ok +
  xlab("Malignant CLDN4 %pos (Seurat object)") +
  ylab("Author T/NK fraction (locked full sample)") +
  ggtitle("Patient/sample unit: CLDN4 vs T/NK")
ggsave_both(p5, "scatter_cldn4_tnk")

p6 <- ggplot(scores, aes(mal_CLDN4_pct, mal_IFN, color = dataset)) +
  geom_point(size = 2.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  theme_ok +
  xlab("Malignant CLDN4 %pos") +
  ylab("Malignant IFN module (Seurat AddModuleScore)") +
  ggtitle("Patient/sample unit: CLDN4 vs malignant IFN")
ggsave_both(p6, "scatter_cldn4_ifn")

p7 <- ggplot(scores, aes(mal_CLDN4_pct, mal_MHC, color = dataset)) +
  geom_point(size = 2.4) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  theme_ok +
  xlab("Malignant CLDN4 %pos") +
  ylab("Malignant MHC module (Seurat AddModuleScore)") +
  ggtitle("Patient/sample unit: CLDN4 vs malignant MHC")
ggsave_both(p7, "scatter_cldn4_mhc")

qdf <- scores[scores$q_seurat %in% c("Q1", "Q4"), ]
qdf$q_seurat <- factor(qdf$q_seurat, levels = c("Q1", "Q4"))
p8 <- ggplot(qdf, aes(q_seurat, locked_frac_tnk, fill = dataset)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1.6) +
  theme_ok + ylab("T/NK fraction") + xlab("Within-cohort CLDN4 %pos quartile") +
  ggtitle("Q4 vs Q1 T/NK (unit = patient/sample)")
ggsave_both(p8, "box_q4q1_tnk")

p9 <- ggplot(qdf, aes(q_seurat, mal_IFN, fill = dataset)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 1.6) +
  theme_ok + ylab("Malignant IFN") + xlab("Within-cohort CLDN4 %pos quartile") +
  ggtitle("Q4 vs Q1 malignant IFN")
ggsave_both(p9, "box_q4q1_ifn")

p10 <- ggplot(n_tab, aes(x = item, y = n)) +
  geom_col(fill = "#4C78A8") + coord_flip() + theme_ok +
  ggtitle("Honest n (cells are not the test unit)")
ggsave_both(p10, "n_honest", w = 8.5, h = 4.2)

# ---- summary json-ish ----
pull <- function(ds, contrast) {
  row <- tests[tests$dataset == ds & tests$contrast == contrast, ]
  if (!nrow(row)) return(NULL)
  as.list(row[1, ])
}

summary <- list(
  seurat_version = as.character(packageVersion("Seurat")),
  seuratobject_version = as.character(packageVersion("SeuratObject")),
  r_version = R.version.string,
  integration_method = integ_method,
  integration_reduction = integ_reduction,
  n_cells = ncol(obj),
  n_units = nrow(scores),
  n_gse131907 = sum(scores$dataset == "GSE131907"),
  n_gse205335 = sum(scores$dataset == "GSE205335"),
  n_ifn_genes = length(features_ifn$IFN),
  n_mhc_genes = length(features_mhc$MHC),
  tests = tests
)
# write a compact json via write
jsonlite::write_json(summary, file.path(RES, "summary.json"), pretty = TRUE, auto_unbox = TRUE)
saveRDS(slot(obj, "meta.data"), file.path(RES, "cell_metadata.rds"))
write.table(slot(obj, "meta.data"), file.path(TAB, "cell_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE)
writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo.txt"))

cat("[done] units=", nrow(scores), " cells=", ncol(obj),
    " method=", integ_method, "\n", sep = "")
print(tests[tests$contrast %in% c("cldn4_vs_tnk", "cldn4_vs_ifn", "cldn4_vs_mhc"),
            c("dataset", "contrast", "n", "rho", "p", "n_low", "n_high", "r")])
