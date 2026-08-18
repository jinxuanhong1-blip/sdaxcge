#!/usr/bin/env Rscript
# GSE207422 — Seurat CreateSeuratObject, Cldn4-only malignant, patient unit.
# Additive. Not dual-high. Not merged with GSE148071.

suppressPackageStartupMessages({
  .libPaths(c("/usr/lib/R/site-library", "/usr/lib/R/library"))
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
tenx_dir <- if (length(args) >= 1) args[[1]] else "data/GSE207422/tenx"
out_dir  <- if (length(args) >= 2) args[[2]] else "methods/seurat_gse207422_cldn4"
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

MIN_MAL <- 20L
MIN_TAIL <- 8L
MIN_MAL_SENS <- 10L

LINEAGES <- list(
  epithelial  = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  T           = c("CD3D", "CD3E", "CD2"),
  NK          = c("NKG7", "GNLY", "FGFBP2"),
  B           = c("CD79A", "MS4A1"),
  plasma      = c("IGHG1", "MZB1"),
  myeloid     = c("LYZ", "CD68", "CD14"),
  neutrophil  = c("CSF3R"),
  fibroblast  = c("COL1A1", "DCN"),
  endothelial = c("VWF", "PECAM1"),
  mast        = c("KIT")
)

A3_NORMAL <- c("SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3")

ISG_CORE <- c(
  "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
  "IFIT5", "IFITM1", "IFITM2", "IFITM3", "MX1", "MX2", "OAS1", "OAS2",
  "OAS3", "OASL", "RSAD2", "USP18", "BST2", "XAF1", "STAT1", "STAT2",
  "IRF7", "IRF9", "DDX58", "IFIH1", "SAMD9", "SAMD9L", "HERC5", "HERC6",
  "EPSTI1", "CMPK2", "PARP9", "DTX3L", "LY6E", "SP100", "SP110", "PLSCR1"
)

MHC1_APM <- c(
  "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1", "TAP2",
  "TAPBP", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "NLRC5",
  "ERAP1", "ERAP2", "CALR", "PDIA3", "CANX", "IRF1"
)

IFN_CORE6 <- c("IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A")

CTRL_OXPHOS <- c(
  "NDUFA1", "NDUFB3", "COX5A", "COX7A2", "UQCRC1", "SDHA",
  "ATP5F1A", "ATP5F1B", "CYCS", "VDAC1"
)

present <- function(genes, universe) intersect(genes, universe)

module_mean <- function(data, genes) {
  genes <- present(genes, rownames(data))
  if (!length(genes)) return(rep(NA_real_, ncol(data)))
  if (length(genes) == 1L) return(as.numeric(data[genes, ]))
  Matrix::colMeans(data[genes, , drop = FALSE])
}

spearman_row <- function(x, y, contrast) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 4L) {
    return(data.frame(contrast = contrast, n = n, rho = NA_real_, p = NA_real_,
                      note = "too_few_samples", stringsAsFactors = FALSE))
  }
  s <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(contrast = contrast, n = n, rho = unname(s$estimate), p = s$p.value,
             note = "", stringsAsFactors = FALSE)
}

wilcoxon_paired <- function(high, low, contrast) {
  ok <- is.finite(high) & is.finite(low)
  high <- high[ok]; low <- low[ok]
  n <- length(high)
  if (n < 3L) {
    return(data.frame(contrast = contrast, n = n,
                      mean_high = if (n) mean(high) else NA_real_,
                      mean_low  = if (n) mean(low) else NA_real_,
                      median_delta = NA_real_, W = NA_real_, p = NA_real_,
                      n_high_gt_low = if (n) sum(high > low) else 0L,
                      note = "too_few_samples", stringsAsFactors = FALSE))
  }
  d <- high - low
  if (all(abs(d) < 1e-12)) {
    w <- 0; p <- 1
  } else {
    wt <- wilcox.test(high, low, paired = TRUE, exact = TRUE)
    w <- unname(wt$statistic); p <- wt$p.value
  }
  data.frame(contrast = contrast, n = n, mean_high = mean(high), mean_low = mean(low),
             median_delta = median(d), W = w, p = p,
             n_high_gt_low = sum(high > low), note = "", stringsAsFactors = FALSE)
}

mwu_row <- function(a, b, contrast) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  na <- length(a); nb <- length(b)
  if (na < 2L || nb < 2L) {
    return(data.frame(contrast = contrast, n_a = na, n_b = nb,
                      mean_a = NA_real_, mean_b = NA_real_, U = NA_real_,
                      p = NA_real_, note = "too_few_samples", stringsAsFactors = FALSE))
  }
  wt <- wilcox.test(a, b, paired = FALSE, exact = TRUE)
  data.frame(contrast = contrast, n_a = na, n_b = nb, mean_a = mean(a),
             mean_b = mean(b), U = unname(wt$statistic), p = wt$p.value,
             note = "", stringsAsFactors = FALSE)
}

fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  if (p < 0.001) return(sprintf("%.2e", p))
  sprintf("%.3f", p)
}

# ------------------------------------------------------------------------------
# 1. Ingest public UMI via Seurat
# ------------------------------------------------------------------------------
feat <- file.path(tenx_dir, "features.tsv.gz")
bc   <- file.path(tenx_dir, "barcodes.tsv.gz")
mtx  <- file.path(tenx_dir, "matrix.mtx.gz")
if (!file.exists(mtx) || !file.exists(feat) || !file.exists(bc)) {
  stop("10x MTX missing under ", tenx_dir, " — stop. Need features/barcodes/matrix.")
}

message("ReadMtx ...")
counts <- ReadMtx(
  mtx = mtx, cells = bc, features = feat,
  feature.column = 1, cell.column = 1, mtx.transpose = FALSE
)
message("CreateSeuratObject ...")
obj <- CreateSeuratObject(
  counts = counts, project = "GSE207422",
  min.cells = 0, min.features = 0, assay = "RNA"
)
rm(counts); gc()

obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)

meta_path <- file.path(out_dir, "sample_metadata.tsv")
smeta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
rownames(smeta) <- smeta$Sample

bcodes <- colnames(obj)
sample_id <- sub("_[^_]+$", "", bcodes)
if (!all(sample_id %in% smeta$Sample)) {
  missing <- sort(unique(sample_id[!sample_id %in% smeta$Sample]))
  stop("barcode samples not in metadata: ", paste(missing, collapse = ","))
}
obj$Sample <- sample_id
obj$Patient <- smeta[sample_id, "Patient"]
obj$timing <- smeta[sample_id, "timing"]
obj$ici_group <- smeta[sample_id, "ici_group"]
obj$Pathologic_Response <- smeta[sample_id, "Pathologic_Response"]
obj$Pathology <- smeta[sample_id, "Pathology"]
obj$RECIST <- smeta[sample_id, "RECIST"]
obj$PD1_Antibody <- smeta[sample_id, "PD1_Antibody"]

genes <- rownames(obj)
data <- GetAssayData(obj, assay = "RNA", layer = "data")
raw  <- GetAssayData(obj, assay = "RNA", layer = "counts")

# ------------------------------------------------------------------------------
# 2. Marker lineage + Cldn4-only A3-malignant (not dual-high)
# ------------------------------------------------------------------------------
lin_names <- names(LINEAGES)
lin_mat <- do.call(rbind, lapply(LINEAGES, function(g) module_mean(data, g)))
rownames(lin_mat) <- lin_names
best_i <- max.col(t(lin_mat), ties.method = "first")
best_v <- lin_mat[cbind(best_i, seq_len(ncol(lin_mat)))]
second <- apply(lin_mat, 2, function(v) sort(v, decreasing = TRUE)[2])
lineage <- lin_names[best_i]
lineage[best_v < 0.15 | (best_v - second) < 0.05] <- "other"
obj$lineage <- lineage
obj$is_epi <- lineage == "epithelial"
obj$is_tnk <- lineage %in% c("T", "NK")

a3_umi <- rep(0, ncol(obj))
for (g in present(A3_NORMAL, genes)) a3_umi <- a3_umi + as.numeric(raw[g, ])
obj$a3_normal_umi <- a3_umi
obj$is_malig_a3 <- obj$is_epi & (a3_umi == 0)
# Cldn4-only: TACSTD2 is never a gate. Dual-high is not computed as a filter.

obj$CLDN4 <- if ("CLDN4" %in% genes) as.numeric(data["CLDN4", ]) else stop("CLDN4 missing")
obj$TACSTD2 <- if ("TACSTD2" %in% genes) as.numeric(data["TACSTD2", ]) else NA_real_

ifn_genes <- present(ISG_CORE, genes)
mhc_genes <- present(MHC1_APM, genes)
core6_genes <- present(IFN_CORE6, genes)
ox_genes <- present(CTRL_OXPHOS, genes)
obj$mod_ifn_isg <- module_mean(data, ifn_genes)
obj$mod_mhc1 <- module_mean(data, mhc_genes)
obj$mod_ifn_core6 <- module_mean(data, core6_genes)
obj$mod_oxphos <- module_mean(data, ox_genes)

# Seurat-native AddModuleScore (sensitivity; not the primary mean-log1p score)
mod_list <- list(
  IFN_ISG = ifn_genes,
  MHC1 = mhc_genes,
  IFN_CORE6 = core6_genes,
  OXPHOS = ox_genes
)
mod_list <- mod_list[vapply(mod_list, length, 1L) >= 3L]
set.seed(42)
obj <- AddModuleScore(obj, features = mod_list, name = "AMS_", ctrl = 100, seed = 42)
ams_cols <- grep("^AMS_", colnames(obj[[]]), value = TRUE)
# Rename positional AMS_1.. to module names
if (length(ams_cols) == length(mod_list)) {
  newn <- paste0("ams_", names(mod_list))
  for (i in seq_along(ams_cols)) obj[[newn[i]]] <- obj[[ams_cols[i]]]
}

# ------------------------------------------------------------------------------
# 3. Patient-level table (honest unit)
# ------------------------------------------------------------------------------
md <- obj@meta.data
md$cell <- colnames(obj)

patient_rows <- list()
paired_rows <- list()
samples <- sort(unique(md$Sample))
for (sid in samples) {
  g <- md[md$Sample == sid, ]
  rec <- data.frame(
    Sample = sid,
    Patient = g$Patient[1],
    timing = g$timing[1],
    ici_group = g$ici_group[1],
    Pathologic_Response = g$Pathologic_Response[1],
    Pathology = g$Pathology[1],
    RECIST = g$RECIST[1],
    PD1_Antibody = g$PD1_Antibody[1],
    n_cells = nrow(g),
    n_epithelial = sum(g$is_epi),
    n_malig_a3 = sum(g$is_malig_a3),
    n_T = sum(g$lineage == "T"),
    n_NK = sum(g$lineage == "NK"),
    n_tnk = sum(g$is_tnk),
    frac_tnk = mean(g$is_tnk),
    stringsAsFactors = FALSE
  )
  for (tag in c("a3", "epi")) {
    sub <- if (tag == "a3") g[g$is_malig_a3, ] else g[g$is_epi, ]
    rec[[paste0("mean_CLDN4_", tag)]] <- if (nrow(sub)) mean(sub$CLDN4) else NA_real_
    rec[[paste0("mean_TACSTD2_", tag)]] <- if (nrow(sub)) mean(sub$TACSTD2) else NA_real_
    rec[[paste0("pct_CLDN4pos_", tag)]] <- if (nrow(sub)) mean(sub$CLDN4 > 0) * 100 else NA_real_
    rec[[paste0("mean_ifn_isg_", tag)]] <- if (nrow(sub)) mean(sub$mod_ifn_isg) else NA_real_
    rec[[paste0("mean_mhc1_", tag)]] <- if (nrow(sub)) mean(sub$mod_mhc1) else NA_real_
    rec[[paste0("mean_ifn_core6_", tag)]] <- if (nrow(sub)) mean(sub$mod_ifn_core6) else NA_real_
    rec[[paste0("mean_oxphos_", tag)]] <- if (nrow(sub)) mean(sub$mod_oxphos) else NA_real_
    if ("ams_IFN_ISG" %in% colnames(sub)) {
      rec[[paste0("ams_ifn_isg_", tag)]] <- if (nrow(sub)) mean(sub$ams_IFN_ISG) else NA_real_
      rec[[paste0("ams_mhc1_", tag)]] <- if (nrow(sub)) mean(sub$ams_MHC1) else NA_real_
    }
  }
  patient_rows[[sid]] <- rec

  for (spec in list(
    list(flag = "is_malig_a3", tag = "a3", min_n = MIN_MAL, min_tail = MIN_TAIL),
    list(flag = "is_malig_a3", tag = "a3_sens", min_n = MIN_MAL_SENS, min_tail = 5L),
    list(flag = "is_epi", tag = "epi", min_n = MIN_MAL, min_tail = MIN_TAIL)
  )) {
    sub <- g[g[[spec$flag]], ]
    n_sub <- nrow(sub)
    if (n_sub < spec$min_n) {
      paired_rows[[length(paired_rows) + 1L]] <- data.frame(
        Sample = sid, Patient = g$Patient[1], timing = g$timing[1],
        ici_group = g$ici_group[1], compartment = spec$tag,
        n_cells = n_sub, n_high = 0L, n_low = 0L, eligible = FALSE,
        reason = sprintf("n_cells<%d", spec$min_n), stringsAsFactors = FALSE
      )
      next
    }
    q25 <- as.numeric(quantile(sub$CLDN4, 0.25, names = FALSE))
    q75 <- as.numeric(quantile(sub$CLDN4, 0.75, names = FALSE))
    hi <- sub[sub$CLDN4 >= q75, ]
    lo <- sub[sub$CLDN4 <= q25, ]
    if (nrow(hi) < spec$min_tail || nrow(lo) < spec$min_tail) {
      paired_rows[[length(paired_rows) + 1L]] <- data.frame(
        Sample = sid, Patient = g$Patient[1], timing = g$timing[1],
        ici_group = g$ici_group[1], compartment = spec$tag,
        n_cells = n_sub, n_high = nrow(hi), n_low = nrow(lo),
        eligible = FALSE, reason = sprintf("tail<%d", spec$min_tail),
        stringsAsFactors = FALSE
      )
      next
    }
    prow <- data.frame(
      Sample = sid, Patient = g$Patient[1], timing = g$timing[1],
      ici_group = g$ici_group[1], compartment = spec$tag,
      n_cells = n_sub, n_high = nrow(hi), n_low = nrow(lo),
      eligible = TRUE, reason = "",
      cldn4_high = mean(hi$CLDN4), cldn4_low = mean(lo$CLDN4),
      ifn_isg_high = mean(hi$mod_ifn_isg), ifn_isg_low = mean(lo$mod_ifn_isg),
      mhc1_high = mean(hi$mod_mhc1), mhc1_low = mean(lo$mod_mhc1),
      ifn_core6_high = mean(hi$mod_ifn_core6), ifn_core6_low = mean(lo$mod_ifn_core6),
      oxphos_high = mean(hi$mod_oxphos), oxphos_low = mean(lo$mod_oxphos),
      umi_high = mean(hi$nCount_RNA), umi_low = mean(lo$nCount_RNA),
      stringsAsFactors = FALSE
    )
    prow$ifn_isg_delta <- prow$ifn_isg_high - prow$ifn_isg_low
    prow$mhc1_delta <- prow$mhc1_high - prow$mhc1_low
    paired_rows[[length(paired_rows) + 1L]] <- prow
  }
}

patients <- do.call(rbind, patient_rows)
# Harmonize paired columns (eligible vs dropped rows differ)
all_cols <- unique(unlist(lapply(paired_rows, names)))
paired <- do.call(rbind, lapply(paired_rows, function(d) {
  miss <- setdiff(all_cols, names(d))
  for (m in miss) d[[m]] <- if (identical(m, "eligible")) FALSE else NA
  d[, all_cols, drop = FALSE]
}))

write.table(patients, file.path(out_dir, "tables", "per_patient.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(paired, file.path(out_dir, "tables", "paired_high_low.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Lineage counts
lin_tab <- as.data.frame(table(md$Sample, md$lineage), stringsAsFactors = FALSE)
names(lin_tab) <- c("Sample", "lineage", "n")
write.table(lin_tab, file.path(out_dir, "tables", "lineage_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Tests — patient is the unit
# ------------------------------------------------------------------------------
tests <- list()
add <- function(kind, df) {
  df$kind <- kind
  tests[[length(tests) + 1L]] <<- df
}

post <- patients[patients$timing == "post", ]
post_a3 <- post[post$n_malig_a3 >= MIN_MAL, ]
post_epi <- post[post$n_epithelial >= MIN_MAL, ]

add("sample_spearman", spearman_row(post_a3$mean_CLDN4_a3, post_a3$frac_tnk,
                                    sprintf("post_a3_CLDN4_vs_tnk_frac_min%d", MIN_MAL)))
add("sample_spearman", spearman_row(post_a3$mean_CLDN4_a3, post_a3$mean_ifn_isg_a3,
                                    sprintf("post_a3_CLDN4_vs_ifn_isg_min%d", MIN_MAL)))
add("sample_spearman", spearman_row(post_a3$mean_CLDN4_a3, post_a3$mean_mhc1_a3,
                                    sprintf("post_a3_CLDN4_vs_mhc1_min%d", MIN_MAL)))
add("sample_spearman", spearman_row(post_a3$mean_CLDN4_a3, post_a3$mean_ifn_core6_a3,
                                    sprintf("post_a3_CLDN4_vs_ifn_core6_min%d", MIN_MAL)))
add("sample_spearman", spearman_row(post_a3$mean_CLDN4_a3, post_a3$mean_oxphos_a3,
                                    sprintf("post_a3_CLDN4_vs_oxphos_min%d", MIN_MAL)))
add("sample_spearman", spearman_row(post_a3$pct_CLDN4pos_a3, post_a3$frac_tnk,
                                    sprintf("post_a3_CLDN4pct_vs_tnk_frac_min%d", MIN_MAL)))

add("sample_spearman", spearman_row(post$mean_CLDN4_a3, post$frac_tnk,
                                    "post_a3_CLDN4_vs_tnk_frac_all_post_incl_tiny"))
add("sample_spearman", spearman_row(post_epi$mean_CLDN4_epi, post_epi$frac_tnk,
                                    "post_epi_CLDN4_vs_tnk_frac"))
add("sample_spearman", spearman_row(post_epi$mean_CLDN4_epi, post_epi$mean_ifn_isg_epi,
                                    "post_epi_CLDN4_vs_ifn_isg"))
add("sample_spearman", spearman_row(post_epi$mean_CLDN4_epi, post_epi$mean_mhc1_epi,
                                    "post_epi_CLDN4_vs_mhc1"))

if ("ams_ifn_isg_a3" %in% names(post_a3)) {
  add("sample_spearman", spearman_row(post_a3$mean_CLDN4_a3, post_a3$ams_ifn_isg_a3,
                                      "post_a3_CLDN4_vs_ams_ifn_isg"))
  add("sample_spearman", spearman_row(post_a3$mean_CLDN4_a3, post_a3$ams_mhc1_a3,
                                      "post_a3_CLDN4_vs_ams_mhc1"))
}

for (comp in c("a3", "a3_sens", "epi")) {
  sub <- paired[paired$compartment == comp & paired$eligible %in% TRUE & paired$timing == "post", ]
  for (mod in c("ifn_isg", "mhc1", "ifn_core6", "oxphos")) {
    add("paired_high_vs_low", {
      r <- wilcoxon_paired(sub[[paste0(mod, "_high")]], sub[[paste0(mod, "_low")]],
                           sprintf("post_%s_%s_CLDN4high_vs_low", comp, mod))
      r$compartment <- comp
      r$module <- mod
      r$n_patients_eligible <- nrow(sub)
      r$patients <- paste(sub$Patient, collapse = ",")
      r
    })
  }
  add("paired_high_vs_low", {
    r <- wilcoxon_paired(sub$umi_high, sub$umi_low, sprintf("post_%s_UMI_CLDN4high_vs_low", comp))
    r$compartment <- comp
    r$module <- "total_UMI"
    r$n_patients_eligible <- nrow(sub)
    r$patients <- paste(sub$Patient, collapse = ",")
    r
  })
}

# ICI labels present: NMPR vs MPR on post patients (pCR folded into MPR)
nmpr <- post_a3$mean_CLDN4_a3[post_a3$ici_group == "NMPR"]
mpr  <- post_a3$mean_CLDN4_a3[post_a3$ici_group == "MPR"]
add("ici_nmpr_vs_mpr", mwu_row(nmpr, mpr, sprintf("post_a3_CLDN4_NMPR_vs_MPR_min%d", MIN_MAL)))
add("ici_nmpr_vs_mpr", mwu_row(post$mean_CLDN4_a3[post$ici_group == "NMPR"],
                               post$mean_CLDN4_a3[post$ici_group == "MPR"],
                               "post_a3_CLDN4_NMPR_vs_MPR_all_post"))
add("ici_nmpr_vs_mpr", mwu_row(post_epi$mean_CLDN4_epi[post_epi$ici_group == "NMPR"],
                               post_epi$mean_CLDN4_epi[post_epi$ici_group == "MPR"],
                               "post_epi_CLDN4_NMPR_vs_MPR"))
add("ici_nmpr_vs_mpr", mwu_row(post$frac_tnk[post$ici_group == "NMPR"],
                               post$frac_tnk[post$ici_group == "MPR"],
                               "post_tnk_frac_NMPR_vs_MPR"))

# Bind tests (union of columns)
test_cols <- unique(unlist(lapply(tests, names)))
tests_df <- do.call(rbind, lapply(tests, function(d) {
  miss <- setdiff(test_cols, names(d))
  for (m in miss) d[[m]] <- NA
  d[, test_cols, drop = FALSE]
}))
write.table(tests_df, file.path(out_dir, "tables", "tests.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

elig <- paired[paired$compartment == "a3" & paired$timing == "post",
               c("Sample", "Patient", "ici_group", "n_cells", "n_high", "n_low",
                 "eligible", "reason")]
write.table(elig, file.path(out_dir, "tables", "eligibility.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Figures
# ------------------------------------------------------------------------------
post$eligible_a3 <- post$n_malig_a3 >= MIN_MAL
post$ici_group <- factor(post$ici_group, levels = c("MPR", "NMPR"))
col_ici <- c(MPR = "#d1495b", NMPR = "#2c6eaf")

p1 <- ggplot(post, aes(mean_CLDN4_a3, frac_tnk, color = ici_group)) +
  geom_point(aes(size = pmax(n_malig_a3, 1)), alpha = 0.9) +
  geom_text(aes(label = Patient), vjust = -0.8, size = 3) +
  scale_color_manual(values = col_ici, drop = FALSE) +
  labs(x = "Malignant CLDN4 mean log1p(CP10k)", y = "T/NK fraction",
       title = "GSE207422 post: malignant CLDN4 vs T/NK (patient)",
       color = "ICI", size = "n malig") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "figures", "fig_cldn4_vs_tnk.png"), p1, width = 6.2, height = 4.6, dpi = 160)
ggsave(file.path(out_dir, "figures", "fig_cldn4_vs_tnk.pdf"), p1, width = 6.2, height = 4.6)

long_mod <- rbind(
  data.frame(Patient = post$Patient, ici_group = post$ici_group,
             CLDN4 = post$mean_CLDN4_a3, module = "IFN ISG", score = post$mean_ifn_isg_a3,
             eligible = post$eligible_a3, stringsAsFactors = FALSE),
  data.frame(Patient = post$Patient, ici_group = post$ici_group,
             CLDN4 = post$mean_CLDN4_a3, module = "MHC-I APM", score = post$mean_mhc1_a3,
             eligible = post$eligible_a3, stringsAsFactors = FALSE)
)
p2 <- ggplot(long_mod[long_mod$eligible, ], aes(CLDN4, score, color = ici_group)) +
  geom_point(size = 3) +
  geom_text(aes(label = Patient), vjust = -0.7, size = 3) +
  facet_wrap(~module, scales = "free_y") +
  scale_color_manual(values = col_ici) +
  labs(x = "Malignant CLDN4 mean log1p(CP10k)", y = "Module mean log1p(CP10k)",
       title = "GSE207422 post: malignant CLDN4 vs IFN / MHC (patient)",
       color = "ICI") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "figures", "fig_cldn4_vs_ifn_mhc.png"), p2, width = 8.2, height = 4.4, dpi = 160)
ggsave(file.path(out_dir, "figures", "fig_cldn4_vs_ifn_mhc.pdf"), p2, width = 8.2, height = 4.4)

paired_a3 <- paired[paired$compartment == "a3" & paired$timing == "post" & paired$eligible %in% TRUE, ]
if (nrow(paired_a3)) {
  pair_long <- rbind(
    data.frame(Patient = paired_a3$Patient, ici_group = paired_a3$ici_group,
               tail = "Q1 low", module = "IFN ISG", score = paired_a3$ifn_isg_low),
    data.frame(Patient = paired_a3$Patient, ici_group = paired_a3$ici_group,
               tail = "Q4 high", module = "IFN ISG", score = paired_a3$ifn_isg_high),
    data.frame(Patient = paired_a3$Patient, ici_group = paired_a3$ici_group,
               tail = "Q1 low", module = "MHC-I APM", score = paired_a3$mhc1_low),
    data.frame(Patient = paired_a3$Patient, ici_group = paired_a3$ici_group,
               tail = "Q4 high", module = "MHC-I APM", score = paired_a3$mhc1_high)
  )
  p3 <- ggplot(pair_long, aes(tail, score, group = Patient, color = ici_group)) +
    geom_line(alpha = 0.7) + geom_point(size = 2.4) +
    facet_wrap(~module, scales = "free_y") +
    scale_color_manual(values = col_ici) +
    labs(x = "Within-patient malignant CLDN4 tail", y = "Module mean log1p(CP10k)",
         title = "Same cells: CLDN4 Q4 vs Q1 IFN / MHC (eligible post patients)",
         color = "ICI") +
    theme_bw(base_size = 11)
  ggsave(file.path(out_dir, "figures", "fig_paired_ifn_mhc.png"), p3, width = 8.2, height = 4.4, dpi = 160)
  ggsave(file.path(out_dir, "figures", "fig_paired_ifn_mhc.pdf"), p3, width = 8.2, height = 4.4)
}

p4 <- ggplot(post, aes(Patient, n_malig_a3, fill = ici_group)) +
  geom_col() +
  geom_hline(yintercept = MIN_MAL, linetype = 2) +
  scale_fill_manual(values = col_ici) +
  labs(y = "A3-malignant cells", title = "Honest n: A3-malignant occupancy (floor 20)",
       fill = "ICI") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "figures", "fig_honest_n.png"), p4, width = 6.4, height = 4.2, dpi = 160)
ggsave(file.path(out_dir, "figures", "fig_honest_n.pdf"), p4, width = 6.4, height = 4.2)

# ------------------------------------------------------------------------------
# 6. Summary JSON + session
# ------------------------------------------------------------------------------
pick <- function(contrast) {
  row <- tests_df[tests_df$contrast == contrast, ]
  if (!nrow(row)) return(NULL)
  as.list(row[1, ])
}

summary <- list(
  dataset = "GSE207422",
  paper = "Hu et al. Genome Med 2023 PMID 36869384",
  seurat_version = as.character(packageVersion("Seurat")),
  seuratobject_version = as.character(packageVersion("SeuratObject")),
  create_seurat_object = TRUE,
  n_cells = ncol(obj),
  n_genes = nrow(obj),
  n_epithelial = as.integer(sum(obj$is_epi)),
  n_malig_a3 = as.integer(sum(obj$is_malig_a3)),
  n_tnk = as.integer(sum(obj$is_tnk)),
  n_T = as.integer(sum(obj$lineage == "T")),
  n_NK = as.integer(sum(obj$lineage == "NK")),
  dual_high = FALSE,
  tacstd2_used_as_gate = FALSE,
  merged_gse148071 = FALSE,
  unit = "patient",
  ici_labels_present = TRUE,
  ici_label_source = "GEO GSE207422_NSCLC_scRNAseq_metadata.xlsx Pathologic Response (pCR folded into MPR)",
  n_post_patients = nrow(post),
  n_post_mpr = sum(post$ici_group == "MPR"),
  n_post_nmpr = sum(post$ici_group == "NMPR"),
  n_post_a3_min20 = nrow(post_a3),
  eligible_paired_a3 = as.character(paired_a3$Patient),
  ifn_genes_used = ifn_genes,
  mhc_genes_used = mhc_genes,
  missing_ifn = setdiff(ISG_CORE, genes),
  missing_mhc = setdiff(MHC1_APM, genes),
  missing_a3_normal = setdiff(A3_NORMAL, genes),
  tests = tests_df
)

jsonlite_ok <- requireNamespace("jsonlite", quietly = TRUE)
if (jsonlite_ok) {
  jsonlite::write_json(summary, file.path(out_dir, "summary.json"),
                       auto_unbox = TRUE, pretty = TRUE, na = "null")
} else {
  # minimal fallback
  con <- file(file.path(out_dir, "summary.json"), "w")
  on.exit(close(con), add = TRUE)
  writeLines(sprintf('{"dataset":"GSE207422","n_cells":%d,"n_genes":%d,"seurat":"%s"}',
                     ncol(obj), nrow(obj), as.character(packageVersion("Seurat"))), con)
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

# compact primary table
prim <- tests_df[tests_df$contrast %in% c(
  sprintf("post_a3_CLDN4_vs_tnk_frac_min%d", MIN_MAL),
  sprintf("post_a3_CLDN4_vs_ifn_isg_min%d", MIN_MAL),
  sprintf("post_a3_CLDN4_vs_mhc1_min%d", MIN_MAL),
  "post_a3_ifn_isg_CLDN4high_vs_low",
  "post_a3_mhc1_CLDN4high_vs_low",
  sprintf("post_a3_CLDN4_NMPR_vs_MPR_min%d", MIN_MAL)
) | grepl("post_a3_(ifn_isg|mhc1)_CLDN4high_vs_low", tests_df$contrast), ]
write.table(prim, file.path(out_dir, "tables", "primary_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("DONE cells=", ncol(obj), " genes=", nrow(obj),
        " malig=", sum(obj$is_malig_a3), " tnk=", sum(obj$is_tnk))
message("post patients=", nrow(post), " a3>=20=", nrow(post_a3),
        " paired eligible=", nrow(paired_a3))
