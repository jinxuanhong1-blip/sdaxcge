#!/usr/bin/env Rscript
# =============================================================================
# deseq2_template.R  -  DESeq2 for KD/KO RNA-seq (si / sh / CRISPRi / CRISPR-KO)
# -----------------------------------------------------------------------------
# INPUT  : a RAW integer count matrix (genes x samples) + a sample sheet.
#          DESeq2 must receive un-normalized counts (NOT FPKM/TPM/RPKM).
# DESIGN : perturbation (KD/KO) vs control (NTC / scramble / safe-harbor sgRNA).
# USAGE  : Rscript deseq2_template.R counts.tsv coldata.tsv out_dir
#
# coldata.tsv (tab-sep, header) minimally:
#   sample     condition   batch
#   NTC_1      control     b1
#   NTC_2      control     b2
#   siCLDN4_1  knockdown   b1
#   siCLDN4_2  knockdown   b2
#
# Key choices you MUST make (see playbook.md):
#   * reference level = the CONTROL guide/siRNA, not "alphabetical first".
#   * include `batch` in the design if samples were processed in batches.
#   * shrink LFCs (apeglm) before ranking / GSEA.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)      # LFC shrinkage
})

args <- commandArgs(trailingOnly = TRUE)
counts_file <- ifelse(length(args) >= 1, args[1], "counts.tsv")
coldata_file <- ifelse(length(args) >= 2, args[2], "coldata.tsv")
out_dir <- ifelse(length(args) >= 3, args[3], "deseq2_out")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- load ------------------------------------------------------------------
cts <- as.matrix(read.delim(counts_file, row.names = 1, check.names = FALSE))
storage.mode(cts) <- "integer"           # DESeq2 needs integers
coldata <- read.delim(coldata_file, row.names = 1, check.names = FALSE)
coldata <- coldata[colnames(cts), , drop = FALSE]   # align order (critical!)
stopifnot(all(rownames(coldata) == colnames(cts)))

coldata$condition <- relevel(factor(coldata$condition), ref = "control")
has_batch <- "batch" %in% colnames(coldata) && length(unique(coldata$batch)) > 1
if (has_batch) coldata$batch <- factor(coldata$batch)

design <- if (has_batch) ~ batch + condition else ~ condition
message("Design: ", deparse(design))

# ---- build + pre-filter ----------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = design)
# Keep genes with >=10 counts in at least "smallest group size" samples.
smallest_grp <- min(table(coldata$condition))
keep <- rowSums(counts(dds) >= 10) >= smallest_grp
dds <- dds[keep, ]
message(sprintf("Genes kept after filtering: %d", nrow(dds)))

# ---- fit -------------------------------------------------------------------
dds <- DESeq(dds)                # size factors + dispersions + Wald test
res_name <- "condition_knockdown_vs_control"
if (!res_name %in% resultsNames(dds)) {
  # fall back to whatever the non-intercept, non-batch coefficient is
  res_name <- grep("condition", resultsNames(dds), value = TRUE)[1]
}
message("Coefficient: ", res_name)

# Shrunken LFC (recommended for ranking / visualization / GSEA input).
res <- lfcShrink(dds, coef = res_name, type = "apeglm")
res <- res[order(res$padj), ]

# ---- QC artifacts ----------------------------------------------------------
vsd <- vst(dds, blind = FALSE)
write.table(as.data.frame(assay(vsd)),
            file.path(out_dir, "vst_matrix.tsv"), sep = "\t", quote = FALSE, col.names = NA)

pca <- plotPCA(vsd, intgroup = intersect(c("condition", "batch"), colnames(coldata)),
               returnData = TRUE)
write.table(pca, file.path(out_dir, "pca.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- results table ---------------------------------------------------------
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]
write.table(res_df, file.path(out_dir, "deseq2_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Ranked file for GSEA pre-rank (rank by shrunken LFC; ties broken by stat).
rnk <- res_df[is.finite(res_df$log2FoldChange), c("gene", "log2FoldChange")]
rnk <- rnk[order(-rnk$log2FoldChange), ]
write.table(rnk, file.path(out_dir, "ranked_for_gsea.rnk"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- on-target + opposite-gene sanity check --------------------------------
focus <- c("CLDN4", "TACSTD2", "EPCAM", "CLDN1", "CLDN3", "CLDN7")
print(res_df[res_df$gene %in% focus, c("gene", "log2FoldChange", "padj")])

message("Done. Outputs in: ", normalizePath(out_dir))
