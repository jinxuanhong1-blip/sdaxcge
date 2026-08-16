#!/usr/bin/env Rscript
# =============================================================================
# batch_correction.R  -  handling batch effects in KD/KO RNA-seq
# -----------------------------------------------------------------------------
# GOLDEN RULE: model the batch, do not "clean" the counts before DE testing.
#   * For DE: put `batch` in the DESeq2/edgeR design (~ batch + condition).
#     This is statistically correct and keeps the raw counts intact.
#   * ComBat-seq: only when you must MERGE datasets that were quantified
#     separately and there is no other way to model the batch in the design.
#   * removeBatchEffect(): ONLY for visualization (PCA/heatmap) on log-CPM/VST,
#     NEVER feed its output into DESeq2/edgeR.
#   * If batch is confounded with condition (all KD in batch1, all control in
#     batch2) NO method can rescue it -> the experiment must be re-run/balanced.
#
# This script shows all three, plus a confounding check.
# =============================================================================

suppressPackageStartupMessages({
  library(sva)          # ComBat_seq
  library(limma)        # removeBatchEffect
  library(edgeR)        # cpm
})

args <- commandArgs(trailingOnly = TRUE)
counts_file  <- ifelse(length(args) >= 1, args[1], "counts.tsv")
coldata_file <- ifelse(length(args) >= 2, args[2], "coldata.tsv")

cts <- as.matrix(read.delim(counts_file, row.names = 1, check.names = FALSE))
meta <- read.delim(coldata_file, row.names = 1, check.names = FALSE)
meta <- meta[colnames(cts), , drop = FALSE]

# ---- 0. confounding check (do this FIRST) ---------------------------------
tab <- table(meta$batch, meta$condition)
print(tab)
if (any(rowSums(tab > 0) == 1) && ncol(tab) > 1 && all(colSums(tab > 0) == 1)) {
  stop("batch is CONFOUNDED with condition; correction is impossible. Re-balance.")
}

# ---- 1. PREFERRED: model batch in the DE design ---------------------------
# DESeq2 : DESeqDataSetFromMatrix(..., design = ~ batch + condition)
# edgeR  : model.matrix(~ batch + condition)
# (see deseq2_template.R / edger_template.R -- nothing else to do here.)

# ---- 2. ComBat-seq: adjust RAW counts when merging separate datasets ------
adj_counts <- ComBat_seq(cts, batch = meta$batch,
                         group = meta$condition)   # preserves condition signal
write.table(adj_counts, "counts_combatseq.tsv", sep = "\t", quote = FALSE, col.names = NA)

# ---- 3. removeBatchEffect: FOR PLOTTING ONLY ------------------------------
logcpm <- cpm(cts, log = TRUE, prior.count = 2)
design0 <- model.matrix(~ condition, data = meta)   # protect condition
vis <- removeBatchEffect(logcpm, batch = meta$batch, design = design0)
write.table(vis, "logcpm_batchcorrected_for_plots.tsv",
            sep = "\t", quote = FALSE, col.names = NA)
message("Wrote ComBat-seq counts (for merging) and batch-corrected logCPM (for plots only).")
