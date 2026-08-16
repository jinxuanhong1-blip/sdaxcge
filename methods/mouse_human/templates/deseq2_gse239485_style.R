# Template only: GSE239485-style nested combination RNA-seq.
# Does not download data. Fill `counts` and `coldata` from GEO/SRA.
# Contrast class is treatment-arm DE. Do not emit R/NR labels.

suppressPackageStartupMessages({
  library(DESeq2)
})

# coldata columns expected (see gse239485_style_metadata.tsv):
#   sample_id, arm, factor_polyIC, factor_anti_PD1, factor_anti_C5aR1
# arm levels: vehicle | polyIC_aPD1 | polyIC_aPD1_aC5aR1

stopifnot(exists("counts"), exists("coldata"))
stopifnot(all(c("arm") %in% colnames(coldata)))

coldata$arm <- factor(
  coldata$arm,
  levels = c("vehicle", "polyIC_aPD1", "polyIC_aPD1_aC5aR1")
)

# Nested encoding (equivalent, useful when adding covariates later):
# combo = 0 vehicle, 1 either treated arm
# add_on = 1 only triple
coldata$combo <- factor(ifelse(coldata$arm == "vehicle", "none", "polyIC_aPD1"))
coldata$add_on <- factor(ifelse(coldata$arm == "polyIC_aPD1_aC5aR1", "aC5aR1", "none"))

# Primary 3-level arm model. Do NOT fit a 2x2x2 factorial: empty cells.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ arm
)
dds <- DESeq(dds)

contrast_D_vs_C <- results(dds, contrast = c("arm", "polyIC_aPD1", "vehicle"))
contrast_T_vs_D <- results(dds, contrast = c("arm", "polyIC_aPD1_aC5aR1", "polyIC_aPD1"))
contrast_T_vs_C <- results(dds, contrast = c("arm", "polyIC_aPD1_aC5aR1", "vehicle"))

# Focal genes: look up by NCBI / Ensembl, not by case-folded symbol alone.
focal_mouse <- c(Tacstd2 = "56753", Cldn4 = "12740", Cldn7 = "53624", Epcam = "17075")

# After DE:
# 1) Rank genes by stat or signed -log10(pvalue) * sign(log2FoldChange).
# 2) Map 1:1 orthologs (NCBI / Ensembl Compara / MGI). Drop many-to-one.
# 3) FGSEA / RRHO only against a human contrast of the SAME class.
# 4) Never rbind these counts with a human matrix.

invisible(list(
  D_vs_C = contrast_D_vs_C,
  T_vs_D = contrast_T_vs_D,
  T_vs_C = contrast_T_vs_C,
  focal_mouse = focal_mouse
))
