## 03_de_compartment.R
## Genome-wide DE with a patient block. Primary: limma-voom + duplicateCorrelation.
## Sensitivity: patient x segment pseudobulk, and (optional) dream / lmer.
##
## 以患者为 block 的全基因组差异表达。主分析：limma-voom + duplicateCorrelation。
## 敏感性：患者 x 分区 pseudobulk，以及可选的 dream / lmer。
##
## Inputs:  results/02_spe.rds
## Outputs: results/03_de_<contrast>.tsv
##          results/03_consensus_correlation.txt
##          results/figures/03_volcano_<contrast>.pdf

source("templates/R/00_setup.R")
need(c("SpatialExperiment", "SummarizedExperiment", "edgeR", "limma"), hard = TRUE)
need(c("variancePartition", "lme4", "lmerTest"), hard = FALSE)

log_step("03  compartment / group DE")

spe <- readRDS(out_path("02_spe.rds"))
cd <- as.data.frame(SummarizedExperiment::colData(spe))
cd$segment <- factor_segment(cd$segment)
cd$patient <- as.factor(cd$patient)
cd$slide <- as.factor(cd$slide)
if (!"group" %in% names(cd)) {
  ## Optional between-patient factor. If absent, within-tissue contrasts still run.
  ## 可选的患者间分组。缺失时仍可做组织内分区对比。
  cd$group <- factor("all")
}

counts <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
## Cross-compartment contrasts use the INTERSECTION gene list (playbook §3.3).
## 跨分区对比使用基因表交集（手册 §3.3）。
genes_inter <- S4Vectors::metadata(spe)$genes_intersection
if (is.null(genes_inter)) genes_inter <- rownames(counts)
counts_x <- counts[intersect(rownames(counts), genes_inter), , drop = FALSE]

## ---------------------------------------------------------------------------
## Design: group-means of segment (and segment x group if group has >1 level)
## 设计矩阵：分区的 group-means（若 group 多于 1 个水平则用 segment x group）
## ---------------------------------------------------------------------------
multi_group <- nlevels(droplevels(as.factor(cd$group))) > 1
if (multi_group) {
  cd$sg <- make_segment_group(cd$segment, cd$group)
  design <- stats::model.matrix(~ 0 + sg, data = cd)
  colnames(design) <- levels(cd$sg)
} else {
  design <- stats::model.matrix(~ 0 + segment, data = cd)
  colnames(design) <- levels(cd$segment)
}

## ---------------------------------------------------------------------------
## Primary: voom + two-round duplicateCorrelation blocked on patient
## 主分析：voom + 以患者为 block 的两轮 duplicateCorrelation
## ---------------------------------------------------------------------------
dge <- edgeR::DGEList(counts = counts_x)
dge <- edgeR::calcNormFactors(dge, method = "TMM")

## filterByExpr on the segment (not patient) to avoid dropping compartment markers.
## 按分区而非患者做 filterByExpr，以免丢掉分区标记基因。
keep <- edgeR::filterByExpr(dge, group = cd$segment)
dge <- dge[keep, , keep.lib.sizes = FALSE]
log_step("voom input: ", nrow(dge), " genes x ", ncol(dge), " AOIs")

v <- limma::voom(dge, design, plot = FALSE)
corfit1 <- limma::duplicateCorrelation(v, design, block = cd$patient)
log_step("consensus.correlation round 1: ", signif(corfit1$consensus.correlation, 3))
if (is.finite(corfit1$consensus.correlation) && corfit1$consensus.correlation > 0.5) {
  corfit2 <- limma::duplicateCorrelation(v, design, block = cd$patient,
                                         correlation = corfit1$consensus.correlation)
  consensus <- corfit2$consensus.correlation
  log_step("consensus.correlation round 2: ", signif(consensus, 3))
} else {
  consensus <- corfit1$consensus.correlation
}
writeLines(c(
  sprintf("round1\t%s", corfit1$consensus.correlation),
  sprintf("used\t%s", consensus),
  "Values near 0: patient clustering is negligible.",
  "Values > 0.3: a naive t-test would have been badly anti-conservative."
), out_path("03_consensus_correlation.txt"))

fit <- limma::lmFit(v, design, block = cd$patient, correlation = consensus)

## Contrasts. Names are taken from the config; formulas are built from segment levels.
## 对比名称来自配置；公式由分区水平自动构建。
cm_list <- list()
if (multi_group) {
  ## Pooled Tumor vs CD45: (Tumor.g1 + Tumor.g2 + ...) - (CD45.g1 + ...)
  segs <- SEGMENT_LEVELS
  grps <- levels(droplevels(as.factor(cd$group)))
  for (i in seq_len(length(segs) - 1)) {
    for (j in seq(i + 1, length(segs))) {
      a <- segs[i]; b <- segs[j]
      plus  <- paste(sanitize_level(a), sanitize_level(grps), sep = ".")
      minus <- paste(sanitize_level(b), sanitize_level(grps), sep = ".")
      plus  <- intersect(plus, colnames(design))
      minus <- intersect(minus, colnames(design))
      if (length(plus) && length(minus)) {
        expr <- paste(paste(plus, collapse = " + "), "-",
                      paste(minus, collapse = " - "))
        cm_list[[paste0(a, "_vs_", b)]] <- expr
      }
    }
  }
  ## Within-segment group contrast for the first two group levels, if present.
  if (length(grps) >= 2) {
    g1 <- sanitize_level(grps[1]); g2 <- sanitize_level(grps[2])
    for (s in segs) {
      a <- paste(sanitize_level(s), g1, sep = ".")
      b <- paste(sanitize_level(s), g2, sep = ".")
      if (all(c(a, b) %in% colnames(design))) {
        cm_list[[paste0(s, "_", grps[1], "_vs_", grps[2])]] <- paste(a, "-", b)
      }
    }
  }
} else {
  segs <- SEGMENT_LEVELS
  for (i in seq_len(length(segs) - 1)) {
    for (j in seq(i + 1, length(segs))) {
      a <- segs[i]; b <- segs[j]
      if (all(c(a, b) %in% colnames(design))) {
        cm_list[[paste0(a, "_vs_", b)]] <- paste(a, "-", b)
      }
    }
  }
}

if (!length(cm_list)) stop("No estimable contrasts. Check segment / group levels.", call. = FALSE)
cm <- limma::makeContrasts(contrasts = cm_list, levels = design)
fit2 <- limma::contrasts.fit(fit, cm)
lfc_thr <- cfg$differential_expression$lfc_threshold %||% 0
if (lfc_thr > 0) {
  fit2 <- limma::treat(fit2, lfc = lfc_thr)
  tt_fun <- function(fit, coef) limma::topTreat(fit, coef = coef, number = Inf, sort.by = "none")
} else {
  fit2 <- limma::eBayes(fit2)
  tt_fun <- function(fit, coef) limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
}

## BH is applied WITHIN each contrast, never across the stacked set.
## BH 在每个对比内部进行，不把所有对比堆在一起校正。
for (cn in colnames(cm)) {
  tt <- tt_fun(fit2, cn)
  tt$gene <- rownames(tt)
  tt$contrast <- cn
  tt$method <- "voom_dupcor"
  write_tsv(tt, out_path(sprintf("03_de_%s.tsv", cn)))
  pdf(fig_path(sprintf("03_volcano_%s.pdf", cn)), width = 6, height = 5)
  adj <- tt$adj.P.Val
  plot(tt$logFC, -log10(pmax(tt$P.Value, 1e-300)),
       pch = 16, cex = 0.4, col = ifelse(adj < cfg$differential_expression$fdr, "red", "grey50"),
       xlab = "log2 FC", ylab = "-log10 p", main = cn)
  targets <- intersect(unlist(cfg$targets$genes), tt$gene)
  if (length(targets)) {
    points(tt[targets, "logFC"], -log10(pmax(tt[targets, "P.Value"], 1e-300)),
           pch = 1, cex = 1.4, col = "blue", lwd = 2)
    text(tt[targets, "logFC"], -log10(pmax(tt[targets, "P.Value"], 1e-300)),
         labels = targets, pos = 3, cex = 0.8, col = "blue")
  }
  dev.off()
}

## ---------------------------------------------------------------------------
## Sensitivity: patient x segment pseudobulk + limma
## 敏感性：患者 x 分区 pseudobulk + limma
## ---------------------------------------------------------------------------
if ("pseudobulk" %in% unlist(cfg$differential_expression$sensitivity_methods)) {
  log_step("sensitivity: patient x segment pseudobulk")
  key <- paste(cd$patient, cd$segment, sep = "||")
  pb <- t(rowsum(t(counts_x), group = key))
  pb_meta <- do.call(rbind, strsplit(colnames(pb), "||", fixed = TRUE))
  pb_cd <- data.frame(patient = pb_meta[, 1],
                      segment = factor_segment(pb_meta[, 2]),
                      stringsAsFactors = FALSE)
  dge_pb <- edgeR::DGEList(counts = pb)
  dge_pb <- edgeR::calcNormFactors(dge_pb)
  des_pb <- stats::model.matrix(~ 0 + segment, data = pb_cd)
  colnames(des_pb) <- levels(pb_cd$segment)
  v_pb <- limma::voom(dge_pb, des_pb)
  fit_pb <- limma::eBayes(limma::lmFit(v_pb, des_pb))
  for (i in seq_len(length(SEGMENT_LEVELS) - 1)) {
    for (j in seq(i + 1, length(SEGMENT_LEVELS))) {
      a <- SEGMENT_LEVELS[i]; b <- SEGMENT_LEVELS[j]
      if (!all(c(a, b) %in% colnames(des_pb))) next
      cm_pb <- limma::makeContrasts(contrasts = paste(a, "-", b), levels = des_pb)
      tt <- limma::topTable(limma::contrasts.fit(fit_pb, cm_pb), number = Inf, sort.by = "none")
      tt$gene <- rownames(tt)
      tt$contrast <- paste0(a, "_vs_", b)
      tt$method <- "pseudobulk"
      write_tsv(tt, out_path(sprintf("03_de_pseudobulk_%s_vs_%s.tsv", a, b)))
    }
  }
}

## ---------------------------------------------------------------------------
## Optional: dream (crossed random effects) / 可选 dream（交叉随机效应）
## ---------------------------------------------------------------------------
if (identical(cfg$differential_expression$primary_method, "dream") ||
    "dream" %in% unlist(cfg$differential_expression$sensitivity_methods)) {
  if (requireNamespace("variancePartition", quietly = TRUE)) {
    log_step("dream mixed model (patient + slide)")
    form <- ~ segment + (1 | patient) + (1 | slide)
    vobj <- variancePartition::voomWithDreamWeights(dge, form, cd)
    fit_d <- variancePartition::dream(vobj, form, cd)
    fit_d <- variancePartition::eBayes(fit_d)
    ## Extract Tumor - CD45 etc. from the coefficient table if estimable.
    ## 若可估计则从系数表提取 Tumor - CD45 等。
    saveRDS(fit_d, out_path("03_dream_fit.rds"))
  } else {
    warning("variancePartition not installed; dream skipped.")
  }
}

save_session_info("03")
log_step("03 done. Contrasts: ", paste(names(cm_list), collapse = ", "))
