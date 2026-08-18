#!/usr/bin/env Rscript
# Seurat Cldn4-only scoring of leftover public mouse lung tumor scRNA.
# Processed GEO matrices only. No FASTQ. No private 8 KL.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

out_root <- "/workspace/methods/seurat_mouse_kl_hunt"
data_root <- "/tmp/geo_dl"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

tnk_genes <- c(
  "Cd3d", "Cd3e", "Cd3g", "Cd2", "Cd8a", "Cd8b1", "Cd4",
  "Nkg7", "Gzma", "Gzmb", "Prf1", "Klrb1c", "Ncr1", "Klrd1", "Klrc1", "Ifng"
)
ifn_genes <- c(
  "Stat1", "Stat2", "Irf1", "Irf7", "Irf9", "Isg15",
  "Ifit1", "Ifit2", "Ifit3", "Mx1", "Oasl2", "Rsad2", "Ifih1", "Ddx58", "Ifnb1",
  "B2m", "H2-K1", "H2-D1", "H2-Q4", "H2-Q6", "H2-Q7",
  "H2-Aa", "H2-Ab1", "H2-Eb1", "Tap1", "Tap2", "Psmb8", "Psmb9", "Nlrc5", "Ciita"
)
epi_genes <- c("Epcam", "Krt8", "Krt18", "Krt19", "Sftpc", "Nkx2-1", "Cdh1")

present <- function(obj, genes) intersect(genes, rownames(obj))

expr_vec <- function(obj, gene) {
  if (!gene %in% rownames(obj)) return(rep(0, ncol(obj)))
  as.numeric(GetAssayData(obj, layer = "data")[gene, ])
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) {
    return(list(r = NA_real_, p = NA_real_, n = sum(ok)))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

safe_mwu <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(a, b, exact = FALSE)$p.value)
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

qc_filter <- function(obj) {
  obj$nFeature_RNA <- obj$nFeature_RNA
  keep <- obj$nFeature_RNA >= 200 & obj$nCount_RNA >= 500 & obj$nFeature_RNA <= 10000
  list(obj = obj[, keep], n_before = ncol(obj), n_after = sum(keep))
}

annotate_and_score <- function(obj, mouse, group, accession) {
  obj$mouse <- mouse
  obj$group <- group
  obj$accession <- accession
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

  tnk_use <- present(obj, tnk_genes)
  ifn_use <- present(obj, ifn_genes)
  if (length(tnk_use) >= 2) {
    obj <- AddModuleScore(obj, features = list(tnk_use), name = "TNK", ctrl = max(8, length(tnk_use)), verbose = FALSE)
    obj$tnk_score <- obj$TNK1
  } else {
    obj$tnk_score <- 0
  }
  if (length(ifn_use) >= 2) {
    obj <- AddModuleScore(obj, features = list(ifn_use), name = "IFN", ctrl = max(8, length(ifn_use)), verbose = FALSE)
    obj$ifn_score <- obj$IFN1
  } else {
    obj$ifn_score <- 0
  }

  obj$Cldn4 <- expr_vec(obj, "Cldn4")
  obj$Epcam <- expr_vec(obj, "Epcam")
  obj$Ptprc <- expr_vec(obj, "Ptprc")
  obj$Cd3e <- expr_vec(obj, "Cd3e")
  obj$Cd3d <- expr_vec(obj, "Cd3d")
  obj$Nkg7 <- expr_vec(obj, "Nkg7")
  obj$Ncr1 <- expr_vec(obj, "Ncr1")
  obj$Krt8 <- expr_vec(obj, "Krt8")

  obj$is_epi <- obj$Epcam > 0 | (obj$Krt8 > 0 & obj$Ptprc == 0)
  obj$is_tnk <- obj$Cd3e > 0 | obj$Cd3d > 0 | obj$Nkg7 > 0 | obj$Ncr1 > 0
  obj$is_immune <- obj$Ptprc > 0
  obj$cell_class <- ifelse(obj$is_tnk, "TNK",
                    ifelse(obj$is_epi & !obj$is_immune, "epithelial",
                    ifelse(obj$is_immune, "other_immune", "other")))

  list(
    obj = obj,
    genes = data.frame(
      accession = accession,
      set = c(rep("Cldn4", as.integer("Cldn4" %in% rownames(obj))),
              rep("TNK", length(tnk_use)),
              rep("IFN_MHC", length(ifn_use)),
              rep("epithelial_marker", length(present(obj, epi_genes)))),
      gene = c(if ("Cldn4" %in% rownames(obj)) "Cldn4" else character(0),
               tnk_use, ifn_use, present(obj, epi_genes)),
      stringsAsFactors = FALSE
    )
  )
}

mouse_table <- function(obj) {
  md <- obj[[]]
  pieces <- split(md, md$mouse, drop = TRUE)
  rows <- lapply(pieces, function(x) {
    epi <- x[x$is_epi, , drop = FALSE]
    data.frame(
      accession = unique(x$accession),
      mouse = unique(x$mouse),
      group = unique(x$group),
      n_cells = nrow(x),
      n_epithelial = sum(x$is_epi),
      n_tnk = sum(x$is_tnk),
      n_immune = sum(x$is_immune),
      frac_epithelial = mean(x$is_epi),
      frac_tnk = mean(x$is_tnk),
      frac_immune = mean(x$is_immune),
      mean_Cldn4_all = mean(x$Cldn4),
      mean_Cldn4_epithelial = if (nrow(epi) >= 10) mean(epi$Cldn4) else NA_real_,
      pct_Cldn4_pos_epithelial = if (nrow(epi) >= 10) mean(epi$Cldn4 > 0) else NA_real_,
      mean_tnk_score = mean(x$tnk_score),
      mean_ifn_epithelial = if (nrow(epi) >= 10) mean(epi$ifn_score) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

read_counts <- function(path) {
  mat <- if (grepl("\\.h5$", path)) Read10X_h5(path) else Read10X(data.dir = path)
  if (is.list(mat) && !is.data.frame(mat) && !inherits(mat, "dgCMatrix")) {
    if ("Gene Expression" %in% names(mat)) mat <- mat[["Gene Expression"]] else mat <- mat[[1]]
  }
  mat
}

contrast_block <- function(mt, accession) {
  cldn <- mt$mean_Cldn4_epithelial
  use <- is.finite(cldn)
  mt2 <- mt[use, , drop = FALSE]
  n <- nrow(mt2)
  if (n < 3) {
    return(data.frame(
      accession = accession,
      endpoint = c("frac_tnk", "mean_ifn_epithelial"),
      n_mice = n,
      spearman_r = NA_real_,
      spearman_p = NA_real_,
      high_n = NA_integer_,
      low_n = NA_integer_,
      high_vs_low_mwu_p = NA_real_,
      high_mean = NA_real_,
      low_mean = NA_real_,
      note = "n<3 mice with >=10 epithelial cells",
      stringsAsFactors = FALSE
    ))
  }
  med <- median(mt2$mean_Cldn4_epithelial)
  hi <- mt2$mean_Cldn4_epithelial > med
  # if median split is empty on one side (ties), use >= median vs < median
  if (sum(hi) == 0 || sum(!hi) == 0) {
    hi <- mt2$mean_Cldn4_epithelial >= med
    if (sum(hi) == n) hi[which.min(mt2$mean_Cldn4_epithelial)] <- FALSE
  }
  rows <- list()
  for (ep in c("frac_tnk", "mean_ifn_epithelial")) {
    sc <- safe_cor(mt2$mean_Cldn4_epithelial, mt2[[ep]])
    rows[[ep]] <- data.frame(
      accession = accession,
      endpoint = ep,
      n_mice = n,
      spearman_r = sc$r,
      spearman_p = sc$p,
      high_n = sum(hi),
      low_n = sum(!hi),
      high_vs_low_mwu_p = safe_mwu(mt2[[ep]][hi], mt2[[ep]][!hi]),
      high_mean = mean(mt2[[ep]][hi], na.rm = TRUE),
      low_mean = mean(mt2[[ep]][!hi], na.rm = TRUE),
      note = sprintf("median split of epithelial Cldn4 (median=%.4f)", med),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

write_accession <- function(accession, mt, genes, honest, contrasts, extra_note) {
  d <- file.path(out_root, accession)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  write_tsv(mt, file.path(d, "mouse_level_scores.tsv"))
  write_tsv(genes, file.path(d, "genes_used.tsv"))
  write_tsv(honest, file.path(d, "honest_n.tsv"))
  write_tsv(contrasts, file.path(d, "contrasts.tsv"))
  readme <- c(
    paste0("# ", accession, " Seurat Cldn4-only scores"),
    "",
    extra_note,
    "",
    "Unit: one row per biological mouse.",
    "Cldn4-only. Processed GEO matrix. No FASTQ.",
    "Seurat NormalizeData + AddModuleScore. Epithelial = Epcam>0 or (Krt8>0 and Ptprc==0).",
    "T/NK cells = Cd3e/Cd3d/Nkg7/Ncr1 > 0. Mouse-level Cldn4 is mean log-normalized Cldn4 in epithelial cells (require >=10 epi cells)."
  )
  writeLines(readme, file.path(d, "README.md"))
}

message("=== GSE179501 ===")
m179 <- read_counts(file.path(data_root, "GSE179501"))
bc <- colnames(m179)
mouse179 <- sub("_.*", "", bc)
map179 <- c(
  CM2260 = "Restored",
  CM2319 = "Restored",
  CM2324 = "Non-Restored",
  CM2328 = "Non-Restored"
)
obj179 <- CreateSeuratObject(m179, project = "GSE179501", min.cells = 0, min.features = 0)
obj179$mouse_raw <- mouse179
qc179 <- qc_filter(obj179)
obj179 <- qc179$obj
ann179 <- annotate_and_score(
  obj179,
  mouse = obj179$mouse_raw,
  group = unname(map179[obj179$mouse_raw]),
  accession = "GSE179501"
)
obj179 <- ann179$obj
mt179 <- mouse_table(obj179)
honest179 <- data.frame(
  accession = "GSE179501",
  n_mice_deposited = 4L,
  n_mice_scored = nrow(mt179),
  n_mice_with_ge10_epi = sum(is.finite(mt179$mean_Cldn4_epithelial)),
  n_cells_before_qc = qc179$n_before,
  n_cells_after_qc = qc179$n_after,
  cldn4_present = "Cldn4" %in% rownames(obj179),
  unit = "biological mouse (barcode prefix CM2260/CM2319/CM2324/CM2328)",
  stringsAsFactors = FALSE
)
ct179 <- contrast_block(mt179, "GSE179501")
# genotype contrast
rest <- mt179$group == "Restored"
geno179 <- data.frame(
  accession = "GSE179501",
  endpoint = "mean_Cldn4_epithelial",
  n_mice = sum(is.finite(mt179$mean_Cldn4_epithelial)),
  spearman_r = NA_real_,
  spearman_p = NA_real_,
  high_n = sum(rest),
  low_n = sum(!rest),
  high_vs_low_mwu_p = safe_mwu(mt179$mean_Cldn4_epithelial[rest], mt179$mean_Cldn4_epithelial[!rest]),
  high_mean = mean(mt179$mean_Cldn4_epithelial[rest], na.rm = TRUE),
  low_mean = mean(mt179$mean_Cldn4_epithelial[!rest], na.rm = TRUE),
  note = "Restored (high_n) vs Non-Restored (low_n) epithelial Cldn4; not a Cldn4 median split",
  stringsAsFactors = FALSE
)
ct179 <- rbind(ct179, geno179)
write_accession(
  "GSE179501", mt179, ann179$genes, honest179, ct179,
  "Lkb1-XTR total-viable leftover (sister of assigned GSE179502). Restored n=2 vs Non-Restored n=2."
)
message("GSE179501 mice:\n")
print(mt179)

message("=== GSE201247 ===")
h5dir <- file.path(data_root, "GSE201247", "h5")
meta201 <- data.frame(
  file = c(
    "GSM6056008_1_ATTAC_filtered_feature_bc_matrix.h5",
    "GSM6056009_2_Kras_filtered_feature_bc_matrix.h5",
    "GSM6056010_3_ATTAC_Kras_filtered_feature_bc_matrix.h5",
    "GSM6056011_4_ATTAC_filtered_feature_bc_matrix.h5",
    "GSM6056012_5_Kras_filtered_feature_bc_matrix.h5",
    "GSM6056013_6_ATTAC_Kras_filtered_feature_bc_matrix.h5"
  ),
  mouse = c("ATTAC_1", "Kras_2", "ATTAC_Kras_3", "ATTAC_4", "Kras_5", "ATTAC_Kras_6"),
  group = c("WT_ATTAC", "Kras", "ATTAC_Kras", "WT_ATTAC", "Kras", "ATTAC_Kras"),
  stringsAsFactors = FALSE
)
objs201 <- list()
n_before201 <- 0L
n_after201 <- 0L
genes201 <- NULL
for (i in seq_len(nrow(meta201))) {
  mat <- read_counts(file.path(h5dir, meta201$file[i]))
  o <- CreateSeuratObject(mat, project = meta201$mouse[i], min.cells = 0, min.features = 0)
  n_before201 <- n_before201 + ncol(o)
  qc <- qc_filter(o)
  n_after201 <- n_after201 + qc$n_after
  ann <- annotate_and_score(qc$obj, mouse = meta201$mouse[i], group = meta201$group[i], accession = "GSE201247")
  objs201[[i]] <- ann$obj
  genes201 <- rbind(genes201, ann$genes)
}
obj201 <- merge(objs201[[1]], y = objs201[-1], add.cell.ids = meta201$mouse, project = "GSE201247")
# merge can drop module columns inconsistently; rebuild mouse table from list
mt201 <- do.call(rbind, lapply(objs201, mouse_table))
genes201 <- unique(genes201)
honest201 <- data.frame(
  accession = "GSE201247",
  n_mice_deposited = 6L,
  n_mice_scored = nrow(mt201),
  n_mice_with_ge10_epi = sum(is.finite(mt201$mean_Cldn4_epithelial)),
  n_cells_before_qc = n_before201,
  n_cells_after_qc = n_after201,
  cldn4_present = any(genes201$gene == "Cldn4"),
  unit = "biological mouse (one filtered h5 library per mouse)",
  stringsAsFactors = FALSE
)
ct201 <- contrast_block(mt201, "GSE201247")
# Kras-bearing (n=4) only
kb <- mt201$group %in% c("Kras", "ATTAC_Kras")
ct201 <- rbind(ct201, contrast_block(mt201[kb, , drop = FALSE], "GSE201247_Kras_bearing_only"))
write_accession(
  "GSE201247", mt201, genes201, honest201, ct201,
  "Kras vs ATTAC whole-lung leftover. Honest n=6 mice (2 WT ATTAC, 2 Kras, 2 ATTAC;Kras)."
)
message("GSE201247 mice:\n")
print(mt201)

message("=== GSE266323 ===")
meta266 <- data.frame(
  dir = c("d10_1", "d10_2", "dTom_1", "dTom_2"),
  mouse = c("d10_1", "d10_2", "dTom_1", "dTom_2"),
  group = c("Malat1_CRISPRa", "Malat1_CRISPRa", "Tomato_control", "Tomato_control"),
  stringsAsFactors = FALSE
)
objs266 <- list()
n_before266 <- 0L
n_after266 <- 0L
genes266 <- NULL
for (i in seq_len(nrow(meta266))) {
  mat <- read_counts(file.path(data_root, "GSE266323", meta266$dir[i]))
  o <- CreateSeuratObject(mat, project = meta266$mouse[i], min.cells = 0, min.features = 0)
  n_before266 <- n_before266 + ncol(o)
  qc <- qc_filter(o)
  n_after266 <- n_after266 + qc$n_after
  ann <- annotate_and_score(qc$obj, mouse = meta266$mouse[i], group = meta266$group[i], accession = "GSE266323")
  objs266[[i]] <- ann$obj
  genes266 <- rbind(genes266, ann$genes)
}
mt266 <- do.call(rbind, lapply(objs266, mouse_table))
genes266 <- unique(genes266)
honest266 <- data.frame(
  accession = "GSE266323",
  n_mice_deposited = 4L,
  n_mice_scored = nrow(mt266),
  n_mice_with_ge10_epi = sum(is.finite(mt266$mean_Cldn4_epithelial)),
  n_cells_before_qc = n_before266,
  n_cells_after_qc = n_after266,
  cldn4_present = any(genes266$gene == "Cldn4"),
  unit = "biological mouse (one GEX library per mouse)",
  stringsAsFactors = FALSE
)
ct266 <- contrast_block(mt266, "GSE266323")
mal <- mt266$group == "Malat1_CRISPRa"
geno266 <- data.frame(
  accession = "GSE266323",
  endpoint = "mean_Cldn4_epithelial",
  n_mice = sum(is.finite(mt266$mean_Cldn4_epithelial)),
  spearman_r = NA_real_,
  spearman_p = NA_real_,
  high_n = sum(mal),
  low_n = sum(!mal),
  high_vs_low_mwu_p = safe_mwu(mt266$mean_Cldn4_epithelial[mal], mt266$mean_Cldn4_epithelial[!mal]),
  high_mean = mean(mt266$mean_Cldn4_epithelial[mal], na.rm = TRUE),
  low_mean = mean(mt266$mean_Cldn4_epithelial[!mal], na.rm = TRUE),
  note = "Malat1 CRISPRa (high_n) vs Tomato control (low_n) epithelial Cldn4; not a Cldn4 median split",
  stringsAsFactors = FALSE
)
ct266 <- rbind(ct266, geno266)
write_accession(
  "GSE266323", mt266, genes266, honest266, ct266,
  "KP LUAD TME leftover. Malat1 CRISPRa n=2 vs Tomato n=2."
)
message("GSE266323 mice:\n")
print(mt266)

all_mt <- rbind(mt179, mt201, mt266)
all_ct <- rbind(ct179, ct201, ct266)
all_h <- rbind(honest179, honest201, honest266)
write_tsv(all_mt, file.path(out_root, "all_mouse_level_scores.tsv"))
write_tsv(all_ct, file.path(out_root, "all_contrasts.tsv"))
write_tsv(all_h, file.path(out_root, "all_honest_n.tsv"))
message("DONE")
message(paste(capture.output(print(all_h)), collapse = "\n"))
message(paste(capture.output(print(all_ct)), collapse = "\n"))
