## Demo pseudobulk DE in R (edgeR + DESeq2) on the exported epithelial pseudobulk.
## Equivalent to pseudobulk_de.py. EXPLORATORY only (small n). Run after run_demo.py.
## 与 pseudobulk_de.py 等价的 R 版本；小样本，仅作探索。
suppressPackageStartupMessages({library(edgeR); library(DESeq2)})

counts <- as.matrix(read.csv("pseudobulk_epithelial_counts.csv",
                             row.names = 1, check.names = FALSE))
design <- read.csv("pseudobulk_design.csv", row.names = 1, check.names = FALSE)

resp_map <- c(pCR = "MPR_like", MPR = "MPR_like", NMPR = "NMPR")
design$grp <- resp_map[as.character(design[["Pathologic.Response"]])]
if (all(is.na(design$grp))) {
  design$grp <- resp_map[as.character(design[["Pathologic Response"]])]
}
keep_s <- rownames(design)[!is.na(design$grp) &
                             design$n_cells_epithelial >= 20 &
                             rownames(design) %in% rownames(counts)]
counts <- t(counts[keep_s, , drop = FALSE])
design <- design[keep_s, , drop = FALSE]
grp <- factor(design$grp, levels = c("NMPR", "MPR_like"))
cat("Group sizes:\n"); print(table(grp))

mm <- model.matrix(~ grp)
y <- DGEList(counts = counts, group = grp)
keep <- filterByExpr(y, design = mm)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, mm)
fit <- glmQLFit(y, mm)
qlf <- glmQLFTest(fit, coef = 2)
res <- topTags(qlf, n = Inf)$table
write.csv(res, "de_epithelial_MPRlike_vs_NMPR_edgeR.csv")

goi <- c("TACSTD2","CLDN4","CLDN3","CLDN7","CLDN18","CDH1","EPCAM",
         "TJP1","OCLN","F11R","ELF3","CXCL13")
cat("\nGenes of interest:\n"); print(res[rownames(res) %in% goi, ])
cat("\nEXPLORATORY: small n, no external replication here.\n")
