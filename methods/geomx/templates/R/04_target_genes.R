## 04_target_genes.R
## Prespecified targets (default TACSTD2, CLDN4): level, specificity, heterogeneity,
## and the spillover controls that decide whether immune-AOI signal is biology.
##
## 事先指定的靶点（默认 TACSTD2、CLDN4）：表达水平、分区特异性、异质性，
## 以及判定免疫 AOI 信号是否为生物学的渗漏对照。
##
## Public-data note: TACSTD2 / CLDN4 are on the Human WTA panel. This script
## does not use any private probe annotation.
## 公开数据说明：TACSTD2 / CLDN4 均在 Human WTA panel 上，不使用任何非公开探针注释。
##
## Inputs:  results/02_spe.rds
##          results/03_de_Tumor_vs_CD45.tsv (optional, for the genome-wide row)
## Outputs: results/04_targets_aoi.tsv
##          results/04_targets_patient.tsv
##          results/04_targets_lmm.tsv
##          results/04_spillover.tsv
##          results/figures/04_targets.pdf

source("templates/R/00_setup.R")
need(c("SpatialExperiment", "SummarizedExperiment"), hard = TRUE)
need(c("lme4", "lmerTest", "emmeans"), hard = FALSE)

log_step("04  TACSTD2 / CLDN4 compartment analysis")

spe <- readRDS(out_path("02_spe.rds"))
cd <- as.data.frame(SummarizedExperiment::colData(spe))
cd$segment <- factor_segment(cd$segment)
cd$patient <- as.factor(cd$patient)
logmat <- as.matrix(SummarizedExperiment::assay(spe, "log2_norm"))
loq <- cd$LOQ
if (is.null(loq) || all(is.na(loq))) loq <- rep(NA_real_, ncol(spe))

targets <- intersect(unlist(cfg$targets$genes), rownames(logmat))
if (!length(targets)) {
  stop("None of the configured target genes are in the matrix: ",
       paste(cfg$targets$genes, collapse = ", "),
       "\nCheck WTA target names (TACSTD2, CLDN4 on Human WTA).",
       call. = FALSE)
}
log_step("targets present: ", paste(targets, collapse = ", "))

spill_genes <- unlist(cfg$targets$spillover_genes)
imm_genes   <- unlist(cfg$targets$immune_control_genes)
cd$spillover <- spillover_index(logmat, spill_genes, transform = "mean_z")
cd$immune_ctrl <- spillover_index(logmat, imm_genes, transform = "mean_z")

## ---------------------------------------------------------------------------
## Q1 Level: detection above LOQ per segment
## Q1 水平：各分区高于 LOQ 的检出
## ---------------------------------------------------------------------------
counts <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
aoi_rows <- lapply(targets, function(g) {
  data.frame(
    gene = g,
    sample_id = colnames(spe),
    patient = as.character(cd$patient),
    segment = as.character(cd$segment),
    slide = as.character(cd$slide),
    expr = as.numeric(logmat[g, ]),
    count = as.numeric(counts[g, ]),
    loq = loq,
    above_loq = if (all(is.na(loq))) NA else as.numeric(counts[g, ]) > loq,
    spillover = cd$spillover,
    immune_ctrl = cd$immune_ctrl,
    area_um2 = cd$area_um2,
    nuclei = cd$nuclei,
    stringsAsFactors = FALSE
  )
})
aoi_df <- do.call(rbind, aoi_rows)
write_tsv(aoi_df, out_path("04_targets_aoi.tsv"))

level_tab <- do.call(rbind, lapply(targets, function(g) {
  d <- aoi_df[aoi_df$gene == g, ]
  do.call(rbind, lapply(SEGMENT_LEVELS, function(s) {
    ds <- d[d$segment == s, ]
    data.frame(
      gene = g, segment = s, n_aoi = nrow(ds),
      n_above_loq = sum(ds$above_loq, na.rm = TRUE),
      pct_above_loq = mean(ds$above_loq, na.rm = TRUE),
      median_log2 = stats::median(ds$expr, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}))
write_tsv(level_tab, out_path("04_targets_level.tsv"))

## ---------------------------------------------------------------------------
## Q2 Specificity: paired LMM with random slope, plus spillover controls
## Q2 特异性：带随机斜率的配对 LMM，外加渗漏对照
## ---------------------------------------------------------------------------
lmm_rows <- list()
if (requireNamespace("lmerTest", quietly = TRUE)) {
  for (g in targets) {
    d <- aoi_df[aoi_df$gene == g, ]
    d$segment <- factor_segment(d$segment)
    ## Co-existing segments -> random intercept AND random slope (playbook §5.1).
    ## 共存分区 -> 随机截距 + 随机斜率（手册 §5.1）。
    form <- expr ~ segment + (1 + segment | patient)
    fit <- try(lmerTest::lmer(form, data = d, REML = TRUE), silent = TRUE)
    if (inherits(fit, "try-error") || lme4::isSingular(fit, tol = 1e-4)) {
      warning(g, ": random-slope model singular or failed; falling back to random intercept.")
      fit <- try(lmerTest::lmer(expr ~ segment + (1 | patient), data = d, REML = TRUE), silent = TRUE)
    }
    if (inherits(fit, "try-error")) {
      warning(g, ": lmer failed entirely.")
      next
    }
    vc <- as.data.frame(lme4::VarCorr(fit))
    icc <- vc$vcov[vc$grp == "patient"][1] / sum(vc$vcov, na.rm = TRUE)
    if (requireNamespace("emmeans", quietly = TRUE)) {
      emm <- emmeans::emmeans(fit, ~ segment)
      pairs <- as.data.frame(emmeans::contrast(emm, "revpairwise"))
      pairs$gene <- g
      pairs$icc_patient <- icc
      pairs$model <- deparse(stats::formula(fit))
      lmm_rows[[g]] <- pairs
    } else {
      cf <- as.data.frame(summary(fit)$coefficients)
      cf$gene <- g
      cf$term <- rownames(summary(fit)$coefficients)
      cf$icc_patient <- icc
      lmm_rows[[g]] <- cf
    }

    ## Spillover: among immune AOIs, does the target track the epithelial index?
    ## 渗漏：在免疫 AOI 中，靶点是否跟随上皮指数？
    imm <- d[is_immune(d$segment), ]
    if (nrow(imm) >= 8 && any(is.finite(imm$spillover))) {
      r <- suppressWarnings(stats::cor(imm$expr, imm$spillover, use = "pairwise.complete.obs"))
      ## Tertile restriction: lowest spillover immune AOIs only.
      ## 三分位限制：仅保留渗漏指数最低的免疫 AOI。
      q <- stats::quantile(imm$spillover, probs = c(1/3), na.rm = TRUE)
      imm_low <- imm[is.finite(imm$spillover) & imm$spillover <= q, ]
      lmm_rows[[paste0(g, "_spill")]] <- data.frame(
        gene = g,
        contrast = "immune_expr_vs_spillover_index",
        estimate = r,
        r_squared = r^2,
        n_immune = nrow(imm),
        n_low_spill = nrow(imm_low),
        median_expr_all_immune = stats::median(imm$expr, na.rm = TRUE),
        median_expr_low_spill = stats::median(imm_low$expr, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
} else {
  warning("lmerTest not installed; writing rank-based paired summaries only.")
  for (g in targets) {
    d <- aoi_df[aoi_df$gene == g, ]
    ## Patient-level paired median: tumor - each immune segment.
    ## 患者层面配对中位数：肿瘤 - 各免疫分区。
    for (imm_s in IMMUNE_SEGMENTS) {
      wide <- stats::reshape(
        d[d$segment %in% c(TUMOR_SEGMENT, imm_s),
          c("patient", "segment", "expr")],
        idvar = "patient", timevar = "segment", direction = "wide"
      )
      tn <- paste0("expr.", TUMOR_SEGMENT)
      im <- paste0("expr.", imm_s)
      if (!all(c(tn, im) %in% names(wide))) next
      delta <- wide[[tn]] - wide[[im]]
      lmm_rows[[paste(g, imm_s, sep = "_")]] <- data.frame(
        gene = g, contrast = paste(TUMOR_SEGMENT, "-", imm_s),
        estimate = stats::median(delta, na.rm = TRUE),
        n_patients = sum(is.finite(delta)),
        pct_positive = mean(delta > 0, na.rm = TRUE),
        wilcox_p = if (sum(is.finite(delta)) >= 6)
          stats::wilcox.test(delta)$p.value else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
}
lmm_tab <- do.call(rbind, lapply(lmm_rows, function(x) {
  x[] <- lapply(x, function(col) if (is.list(col)) sapply(col, as.character) else col)
  x
}))
rownames(lmm_tab) <- NULL
write_tsv(lmm_tab, out_path("04_targets_lmm.tsv"))

## ---------------------------------------------------------------------------
## Q3 Heterogeneity: within-patient ROI variance, ICC, ROIs needed for R >= 0.8
## Q3 异质性：患者内 ROI 方差、ICC、达到可靠性 0.8 所需 ROI 数
## ---------------------------------------------------------------------------
het_rows <- lapply(targets, function(g) {
  do.call(rbind, lapply(SEGMENT_LEVELS, function(s) {
    ds <- aoi_df[aoi_df$gene == g & aoi_df$segment == s, ]
    rho <- icc_oneway(ds$expr, ds$patient)
    m <- mean(table(ds$patient))
    data.frame(
      gene = g, segment = s,
      n_patient = length(unique(ds$patient)),
      mean_m = m,
      icc = rho,
      reliability = spearman_brown(rho, m),
      n_roi_for_R08 = n_roi_for_reliability(rho, 0.8),
      n_eff = n_eff(length(unique(ds$patient)), m, rho),
      stringsAsFactors = FALSE
    )
  }))
})
d_ok <- function(df, g) TRUE
het_tab <- do.call(rbind, het_rows)
write_tsv(het_tab, out_path("04_targets_heterogeneity.tsv"))

## Patient-level paired deltas (used by 05 / 06 as specificity_score).
## 患者层面配对差值（供 05 / 06 作为 specificity_score）。
pat_rows <- lapply(targets, function(g) {
  d <- aoi_df[aoi_df$gene == g, ]
  agg <- aggregate_aoi(d$expr, d$patient, d$segment, d$nuclei,
                       method = cfg$aggregation$method %||% "median",
                       min_n = cfg$aggregation$min_aoi_per_patient %||% 1)
  agg$gene <- g
  agg
})
pat_tab <- do.call(rbind, pat_rows)
write_tsv(pat_tab, out_path("04_targets_patient.tsv"))

## ---------------------------------------------------------------------------
## Figures / 图
## ---------------------------------------------------------------------------
pdf(fig_path("04_targets.pdf"), width = 9, height = 6)
for (g in targets) {
  d <- aoi_df[aoi_df$gene == g, ]
  d$segment <- factor_segment(d$segment)
  graphics::boxplot(expr ~ segment, data = d, ylab = "log2 normalized",
                    main = sprintf("%s  by segment (AOI-level)", g),
                    col = c("#4C78A8", "#F58518", "#54A24B")[seq_along(SEGMENT_LEVELS)])
  ## Paired patient means
  pm <- stats::aggregate(expr ~ patient + segment, d, mean)
  wide <- stats::reshape(pm, idvar = "patient", timevar = "segment", direction = "wide")
  if (all(c(paste0("expr.", TUMOR_SEGMENT),
            paste0("expr.", IMMUNE_SEGMENTS[1])) %in% names(wide))) {
    graphics::stripchart(
      list(tumor = wide[[paste0("expr.", TUMOR_SEGMENT)]],
           immune = wide[[paste0("expr.", IMMUNE_SEGMENTS[1])]]),
      vertical = TRUE, method = "jitter", pch = 16, col = c("#4C78A8", "#F58518"),
      ylab = "patient-mean log2",
      main = sprintf("%s  paired patient means  %s vs %s",
                     g, TUMOR_SEGMENT, IMMUNE_SEGMENTS[1])
    )
  }
  imm <- d[is_immune(d$segment), ]
  if (nrow(imm) && any(is.finite(imm$spillover))) {
    plot(imm$spillover, imm$expr, pch = 16, col = as.integer(factor(imm$segment)),
         xlab = "epithelial spillover index (z)", ylab = paste(g, "log2 in immune AOIs"),
         main = sprintf("%s  immune AOI vs spillover  r = %.2f",
                        g, suppressWarnings(stats::cor(imm$expr, imm$spillover,
                                                       use = "pairwise"))))
    legend("topleft", legend = unique(imm$segment),
           col = seq_along(unique(imm$segment)), pch = 16, bty = "n")
  }
}
dev.off()

## Reporting block (fill-in numbers for the methods / results paragraph).
## 报告块：给方法学/结果段落填空的数字。
sink(out_path("04_reporting_block.txt"))
for (g in targets) {
  cat("==== ", g, " ====\n", sep = "")
  print(level_tab[level_tab$gene == g, ])
  cat("\nheterogeneity\n")
  print(het_tab[het_tab$gene == g, ])
  cat("\n")
}
sink()

save_session_info("04")
log_step("04 done")
