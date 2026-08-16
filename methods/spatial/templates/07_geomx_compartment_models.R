#!/usr/bin/env Rscript
# 07 · GeoMx compartment models. GeoMx 区室（compartment）模型。
# ---------------------------------------------------------------------------
# GeoMx AOIs are collected per compartment (e.g. tumour / immune / stroma, or
# PanCK+ vs PanCK-) within ROIs on slides from patients. Differential expression
# between compartments must account for this nested structure: multiple AOIs per
# ROI, multiple ROIs per slide/patient. A linear mixed model (LMM) with a random
# effect for slide/patient is the GeoMx-recommended approach.
# GeoMx 的 AOI 按区室采集（如 肿瘤/免疫/基质，或 PanCK+ vs PanCK-），嵌套于
# ROI 与患者切片。区室间差异表达需考虑该嵌套结构（每 ROI 多个 AOI、每切片多个
# ROI），采用带切片/患者随机效应的线性混合模型（LMM），此为 GeoMx 推荐做法。
#
# Two interchangeable engines (pick with --engine):
#   mixed   : GeomxTools::mixedModelDE  (LMM on log2 Q3 counts)  [default]
#   standr  : standR + limma-voom with duplicateCorrelation (block = patient)
# 两种可选引擎：GeomxTools 的 LMM；或 standR + limma-voom（以患者为 block）。
#
# Tools: GeomxTools, standR, limma, edgeR (Bioconductor 3.18+, 2024+).
# No results are fabricated; contrasts/statistics are produced by the model.
# 工具：GeomxTools、standR、limma、edgeR。不伪造结果；统计量由模型产生。
#
# Usage:
#   Rscript 07_geomx_compartment_models.R --rds <geomx_target_qnorm.rds> \
#       --compartment segment --contrast Tumor,Immune --subject patient \
#       --engine mixed --out methods/spatial/demo/out/geomx
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(GeomxTools)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--rds", type = "character", help = "output of template 06"),
  make_option("--compartment", type = "character", default = "segment",
              help = "pData column naming the compartment (e.g. segment)"),
  make_option("--contrast", type = "character",
              help = "two compartment levels 'A,B' to compare (A vs B)"),
  make_option("--subject", type = "character", default = "patient",
              help = "random-effect grouping (slide or patient id)"),
  make_option("--engine", type = "character", default = "mixed",
              help = "mixed (GeomxTools LMM) or standr (limma-voom)"),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--out", type = "character", default = "methods/spatial/demo/out/geomx")
)))

dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
target_geomx <- readRDS(opt$rds)
levels_ab <- strsplit(opt$contrast, ",")[[1]]
stopifnot(length(levels_ab) == 2)

# log2 of Q3-normalized counts for linear modelling. 对 Q3 归一化计数取 log2。
assayDataElement(target_geomx, "log_q") <-
  assayDataApply(target_geomx, 2, FUN = log, base = 2, elt = "q_norm")

pd <- pData(target_geomx)
comp <- opt$compartment
subj <- opt$subject
stopifnot(comp %in% colnames(pd), subj %in% colnames(pd))

# Restrict to the two compartments being contrasted. 仅保留待比较的两个区室。
keep <- pd[[comp]] %in% levels_ab
target_geomx <- target_geomx[, keep]
pData(target_geomx)[[comp]] <- factor(pData(target_geomx)[[comp]], levels = levels_ab)

if (opt$engine == "mixed") {
  # LMM: expression ~ compartment + (1 | subject). Random intercept per subject
  # absorbs patient/slide-level correlation between AOIs.
  # LMM：表达 ~ 区室 + (1 | 患者)。患者随机截距吸收 AOI 间的相关性。
  fml <- as.formula(sprintf("~ %s + (1 | %s)", comp, subj))
  res <- mixedModelDE(target_geomx, elt = "log_q", modelFormula = fml,
                      groupVar = comp, nCores = 1, multiCore = FALSE)

  tab <- do.call(rbind, lapply(seq_len(ncol(res)), function(i) {
    x <- res["lsmeans", i][[1]]
    data.frame(gene = colnames(res)[i],
               contrast = rownames(x)[1],
               log2FC = as.numeric(x[1, "Estimate"]),
               pval = as.numeric(x[1, "Pr(>|t|)"]))
  }))
  tab$FDR <- p.adjust(tab$pval, method = "BH")

} else if (opt$engine == "standr") {
  suppressPackageStartupMessages({ library(standR); library(limma); library(edgeR) })
  spe <- as.SpatialExperiment(target_geomx, normData = "q_norm", forceRaw = FALSE)
  design <- model.matrix(as.formula(sprintf("~ 0 + %s", comp)), data = colData(spe))
  colnames(design) <- make.names(levels(colData(spe)[[comp]]))
  logmat <- log2(assay(spe, "q_norm") + 1)
  block <- colData(spe)[[subj]]
  dupcor <- duplicateCorrelation(logmat, design, block = block)
  fit <- lmFit(logmat, design, block = block, correlation = dupcor$consensus)
  contr <- makeContrasts(contrasts = paste(make.names(levels_ab), collapse = "-"),
                         levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contr))
  tab <- topTable(fit2, number = Inf)
  tab$gene <- rownames(tab)
  tab$log2FC <- tab$logFC; tab$pval <- tab$P.Value; tab$FDR <- tab$adj.P.Val

} else {
  stop("unknown --engine; use 'mixed' or 'standr'")
}

tab <- tab[order(tab$FDR), ]
out_csv <- file.path(opt$out, sprintf("compartment_DE_%s_vs_%s.csv",
                                      levels_ab[1], levels_ab[2]))
write.csv(tab, out_csv, row.names = FALSE)
message(sprintf("[07_geomx] %s vs %s: %d genes at FDR < %.2f -> %s",
                levels_ab[1], levels_ab[2], sum(tab$FDR < opt$fdr, na.rm = TRUE),
                opt$fdr, out_csv))
