#!/usr/bin/env Rscript
# =============================================================================
# edger_template.R  -  edgeR quasi-likelihood for KD/KO RNA-seq
# -----------------------------------------------------------------------------
# INPUT : RAW integer counts (genes x samples) + sample sheet (see DESeq2 tmpl).
# WHY QL: edgeR's quasi-likelihood F-test (glmQLFit/glmQLFTest) controls the
#         type-I error better than the LRT for the small n (2-3) typical of
#         perturbation screens.
# USAGE : Rscript edger_template.R counts.tsv coldata.tsv out_dir
# =============================================================================

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
counts_file  <- ifelse(length(args) >= 1, args[1], "counts.tsv")
coldata_file <- ifelse(length(args) >= 2, args[2], "coldata.tsv")
out_dir      <- ifelse(length(args) >= 3, args[3], "edger_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cts <- as.matrix(read.delim(counts_file, row.names = 1, check.names = FALSE))
coldata <- read.delim(coldata_file, row.names = 1, check.names = FALSE)
coldata <- coldata[colnames(cts), , drop = FALSE]
stopifnot(all(rownames(coldata) == colnames(cts)))

condition <- relevel(factor(coldata$condition), ref = "control")
has_batch <- "batch" %in% colnames(coldata) && length(unique(coldata$batch)) > 1

# Design: block on batch if present (edgeR removes batch as a nuisance term).
if (has_batch) {
  batch <- factor(coldata$batch)
  design <- model.matrix(~ batch + condition)
} else {
  design <- model.matrix(~ condition)
}
colnames(design) <- make.names(colnames(design))
coef <- tail(grep("condition", colnames(design), value = TRUE), 1)
message("Testing coefficient: ", coef)

y <- DGEList(counts = cts, group = condition)
keep <- filterByExpr(y, design = design)          # edgeR's recommended filter
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")           # TMM normalization
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)                         # quasi-likelihood
qlf <- glmQLFTest(fit, coef = coef)

tt <- topTags(qlf, n = Inf)$table
tt$gene <- rownames(tt)
tt <- tt[, c("gene", setdiff(colnames(tt), "gene"))]
write.table(tt, file.path(out_dir, "edger_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Ranked list for GSEA pre-rank: signed -log10(PValue) is a robust metric.
tt$rank_metric <- sign(tt$logFC) * -log10(pmax(tt$PValue, .Machine$double.xmin))
rnk <- tt[order(-tt$rank_metric), c("gene", "rank_metric")]
write.table(rnk, file.path(out_dir, "ranked_for_gsea.rnk"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

focus <- c("CLDN4", "TACSTD2", "EPCAM", "CLDN1", "CLDN3", "CLDN7")
print(tt[tt$gene %in% focus, c("gene", "logFC", "FDR")])
message("Done. Outputs in: ", normalizePath(out_dir))
