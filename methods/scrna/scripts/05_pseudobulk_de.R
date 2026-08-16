## PRIMARY differential expression = PSEUDOBULK with edgeR-QLF (and DESeq2 cross-check).
## Input: a samples x genes pseudobulk COUNT matrix + a sample design table.
## (Produce these from 05_pseudobulk_de.py, muscat::aggregateData, or manual rowsum.)
## 主要差异分析 = 伪bulk（edgeR-QLF，DESeq2交叉验证）。样本=重复单位。
suppressPackageStartupMessages({library(edgeR); library(DESeq2)})

args <- commandArgs(trailingOnly = TRUE)
counts_csv <- args[[1]]   # rows = samples, cols = genes (raw summed counts)
design_csv <- args[[2]]   # rows = samples; must contain the group column
group_col  <- ifelse(length(args) >= 3, args[[3]], "path_response")
ref_level  <- ifelse(length(args) >= 4, args[[4]], "NMPR")
out_prefix <- ifelse(length(args) >= 5, args[[5]], "pseudobulk_de")

cnt <- as.matrix(read.csv(counts_csv, row.names = 1, check.names = FALSE))
cnt <- t(cnt)                                   # -> genes x samples for edgeR
meta <- read.csv(design_csv, row.names = 1, check.names = FALSE)
meta <- meta[colnames(cnt), , drop = FALSE]
grp  <- factor(meta[[group_col]])
grp  <- relevel(grp, ref = ref_level)
stopifnot(nlevels(grp) == 2)

## Optional covariates present in meta, e.g. histology / batch:
covars <- intersect(c("histology", "Pathology", "batch"), colnames(meta))
form   <- as.formula(paste("~", paste(c(covars, "grp"), collapse = " + ")))
design <- model.matrix(form, data = cbind(meta, grp = grp))

## ---- edgeR QLF ----
y <- DGEList(counts = cnt, group = grp)
keep <- filterByExpr(y, design = design)
message(sprintf("Genes kept by filterByExpr: %d / %d", sum(keep), nrow(y)))
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)                          # TMM
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = ncol(design))      # last coef = grp effect
res_edger <- topTags(qlf, n = Inf)$table
write.csv(res_edger, paste0(out_prefix, "_edgeR.csv"))
message("edgeR top hits:"); print(head(res_edger, 15))

## ---- DESeq2 cross-check ----
dds <- DESeqDataSetFromMatrix(countData = cnt,
                              colData   = cbind(meta, grp = grp),
                              design    = form)
dds <- DESeq(dds)
res_deseq <- as.data.frame(results(dds, name = tail(resultsNames(dds), 1)))
res_deseq <- res_deseq[order(res_deseq$padj), ]
write.csv(res_deseq, paste0(out_prefix, "_DESeq2.csv"))
message("DESeq2 top hits:"); print(head(res_deseq, 15))

## Report log2FC + FDR, and ALWAYS plot per-sample points for top genes.
## 报告 log2FC 与 FDR，并对top基因绘制每样本点图。
message("Done -> ", out_prefix, "_{edgeR,DESeq2}.csv")
