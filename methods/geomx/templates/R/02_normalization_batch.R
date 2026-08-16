## 02_normalization_batch.R
## Compare Q3 / TMM / quantile (optional GeoDiff), pick the PRIMARY method declared
## in the config BEFORE looking at DE, then optional RUV4. Writes a comparison
## figure and the chosen assay onto the SpatialExperiment.
##
## 比较 Q3 / TMM / quantile（可选 GeoDiff），使用配置中事先指定的主方法，
## 再可选做 RUV4。写出比较图，并把选定 assay 写回 SpatialExperiment。
##
## Inputs:  results/01_spe.rds
## Outputs: results/02_spe.rds
##          results/figures/02_rle_pca.pdf
##          results/02_norm_choice.txt

source("templates/R/00_setup.R")
need(c("SpatialExperiment", "SummarizedExperiment", "S4Vectors"), hard = TRUE)
need(c("edgeR", "limma"), hard = TRUE)
need(c("standR", "ggplot2"), hard = FALSE)

log_step("02  normalization + batch")

spe <- readRDS(out_path("01_spe.rds"))
counts <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
cd <- as.data.frame(SummarizedExperiment::colData(spe))

## ---------------------------------------------------------------------------
## Size-factor / rank normalizations / 缩放与秩归一化
## ---------------------------------------------------------------------------
log2p <- function(x) log2(x + 1)

norm_q3 <- function(m) {
  q3 <- apply(m, 2, stats::quantile, probs = 0.75, na.rm = TRUE)
  q3[q3 <= 0] <- stats::median(q3[q3 > 0], na.rm = TRUE)
  sweep(m, 2, q3, "/") * stats::median(q3)
}

norm_tmm <- function(m) {
  dge <- edgeR::DGEList(counts = m)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  edgeR::cpm(dge, log = FALSE)
}

norm_quantile <- function(m) {
  ## Rank-based quantile normalization on log2(count+1), then back-transform
  ## so the stored assay stays on a count-like scale for voom.
  ## 在 log2(count+1) 上做分位数归一化，再反变换，使 assay 保持类计数尺度供 voom 使用。
  lg <- log2p(m)
  qn <- limma::normalizeQuantiles(lg)
  pmax(2^qn - 1, 0)
}

methods_run <- unique(c(cfg$normalization$primary, unlist(cfg$normalization$sensitivity)))
assays <- list()
for (meth in methods_run) {
  log_step("normalizing: ", meth)
  assays[[meth]] <- switch(meth,
    q3 = norm_q3(counts),
    tmm = norm_tmm(counts),
    quantile = norm_quantile(counts),
    geodiff = {
      warning("GeoDiff path is a stub in this template: falling back to TMM. ",
              "See the GeoDiff WTA vignette for Poisson-background + NB-threshold size factors.")
      norm_tmm(counts)
    },
    stop("Unknown normalization: ", meth)
  )
}

primary <- cfg$normalization$primary
SummarizedExperiment::assay(spe, "counts_norm") <- assays[[primary]]
SummarizedExperiment::assay(spe, "log2_norm") <- log2p(assays[[primary]])
for (nm in setdiff(names(assays), primary)) {
  SummarizedExperiment::assay(spe, paste0("log2_", nm)) <- log2p(assays[[nm]])
}

## ---------------------------------------------------------------------------
## Diagnostics: RLE + PCA coloured by segment / slide / area / detection
## 诊断：RLE 与按分区/玻片/面积/检出率着色的 PCA
## ---------------------------------------------------------------------------
rle_matrix <- function(logmat) {
  med <- apply(logmat, 1, stats::median, na.rm = TRUE)
  sweep(logmat, 1, med, "-")
}

pca_df <- function(logmat, cd, ntop = 2000) {
  rv <- apply(logmat, 1, stats::var, na.rm = TRUE)
  keep <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pc <- stats::prcomp(t(logmat[keep, , drop = FALSE]), scale. = TRUE)
  data.frame(cd,
             PC1 = pc$x[, 1],
             PC2 = pc$x[, 2],
             var1 = summary(pc)$importance[2, 1],
             var2 = summary(pc)$importance[2, 2],
             stringsAsFactors = FALSE)
}

pdf(fig_path("02_rle_pca.pdf"), width = 10, height = 8)
on.exit(dev.off(), add = TRUE)
for (meth in names(assays)) {
  lg <- log2p(assays[[meth]])
  rle <- rle_matrix(lg)
  graphics::boxplot(rle, outline = FALSE, las = 2, cex.axis = 0.4,
                    main = sprintf("RLE  %s", meth), ylab = "RLE (log2)")
  pcs <- pca_df(lg, cd)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  col_seg <- as.integer(factor(pcs$segment, levels = SEGMENT_LEVELS))
  plot(pcs$PC1, pcs$PC2, col = col_seg, pch = 16,
       xlab = sprintf("PC1 (%.0f%%)", 100 * pcs$var1[1]),
       ylab = sprintf("PC2 (%.0f%%)", 100 * pcs$var2[1]),
       main = sprintf("PCA %s  by segment", meth))
  legend("topright", legend = SEGMENT_LEVELS, col = seq_along(SEGMENT_LEVELS),
         pch = 16, cex = 0.7, bty = "n")
  plot(pcs$PC1, pcs$PC2, col = as.integer(factor(pcs$slide)), pch = 16,
       main = sprintf("PCA %s  by slide", meth), xlab = "PC1", ylab = "PC2")
  if ("area_um2" %in% names(pcs) && any(is.finite(pcs$area_um2))) {
    plot(pcs$PC1, pcs$PC2,
         col = grDevices::hcl.colors(10, "viridis")[
           as.integer(cut(log10(pmax(pcs$area_um2, 1)), 10))],
         pch = 16, main = sprintf("PCA %s  by log10 area", meth),
         xlab = "PC1", ylab = "PC2")
  } else {
    plot.new(); title("no area_um2")
  }
  if ("GeneDetectionRate" %in% names(pcs)) {
    plot(pcs$PC1, pcs$PC2,
         col = grDevices::hcl.colors(10, "plasma")[
           as.integer(cut(pcs$GeneDetectionRate, 10))],
         pch = 16, main = sprintf("PCA %s  by detection rate", meth),
         xlab = "PC1", ylab = "PC2")
  } else {
    plot.new(); title("no GeneDetectionRate")
  }
  par(op)
}

## ---------------------------------------------------------------------------
## Batch correction (visualization only) / 批次校正（仅用于可视化）
## Inference keeps the batch term in the model (see 03). The corrected matrix
## is stored as log2_batchcorr and must not be reused for DE.
## 推断时把批次项留在模型里（见 03）。校正矩阵存为 log2_batchcorr，不得用于 DE。
## ---------------------------------------------------------------------------
bc <- cfg$normalization$batch_correction
if (!is.null(bc) && !identical(bc$method, "none")) {
  log_step("batch correction: ", bc$method, " (visualization / clustering only)")
  logmat <- SummarizedExperiment::assay(spe, "log2_norm")
  if (identical(bc$method, "limma")) {
    ## removeBatchEffect is appropriate when slide is not confounded with biology.
    ## 仅当玻片与生物学因素不混杂时才适合 removeBatchEffect。
    chk <- confound_check(cd, "slide", "cohort", "slide vs cohort (batch correction)")
    report_confound(chk)
    design <- stats::model.matrix(~ segment, data = cd)
    corr <- limma::removeBatchEffect(logmat, batch = cd$slide, design = design)
    SummarizedExperiment::assay(spe, "log2_batchcorr") <- corr
  } else if (identical(bc$method, "ruv4")) {
    if (requireNamespace("standR", quietly = TRUE)) {
      ## standR::findNCGs + geomxBatchCorrection on a SpatialExperiment.
      ## Requires biology factors that should be *kept*.
      tryCatch({
        spe2 <- spe
        SummarizedExperiment::assay(spe2, "logcounts") <- logmat
        ncg <- standR::findNCGs(spe2, n_assay = "logcounts",
                                batch_name = "slide", n = bc$n_ncg %||% 300)
        spe2 <- standR::geomxBatchCorrection(spe2, factors = "segment",
                                             NCGs = ncg, k = bc$k %||% 3)
        if ("logcounts" %in% SummarizedExperiment::assayNames(spe2)) {
          SummarizedExperiment::assay(spe, "log2_batchcorr") <-
            SummarizedExperiment::assay(spe2, "logcounts")
        }
      }, error = function(e) {
        warning("standR RUV4 failed: ", conditionMessage(e),
                "  Falling back to limma::removeBatchEffect.")
        design <- stats::model.matrix(~ segment, data = cd)
        SummarizedExperiment::assay(spe, "log2_batchcorr") <-
          limma::removeBatchEffect(logmat, batch = cd$slide, design = design)
      })
    } else {
      warning("standR not installed; RUV4 skipped, using limma::removeBatchEffect.")
      design <- stats::model.matrix(~ segment, data = cd)
      SummarizedExperiment::assay(spe, "log2_batchcorr") <-
        limma::removeBatchEffect(logmat, batch = cd$slide, design = design)
    }
  }
}

S4Vectors::metadata(spe)$normalization <- list(
  primary = primary,
  sensitivity = cfg$normalization$sensitivity,
  batch = bc,
  note = paste("Primary method was declared in the config before DE.",
               "Batch-corrected assay is for visualization only.")
)

saveRDS(spe, out_path("02_spe.rds"))
writeLines(c(
  sprintf("primary_normalization\t%s", primary),
  sprintf("sensitivity\t%s", paste(cfg$normalization$sensitivity, collapse = ",")),
  sprintf("batch_method\t%s", bc$method %||% "none"),
  "DE must use counts + voom (or log2_norm), NOT log2_batchcorr."
), out_path("02_norm_choice.txt"))

save_session_info("02")
log_step("02 done -> ", out_path("02_spe.rds"))
