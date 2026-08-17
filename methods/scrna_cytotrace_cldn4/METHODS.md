# Methods — CLDN4 vs CytoTRACE-like potency (GSE131907)

Additive module. Lives under `methods/scrna_cytotrace_cldn4/`. **CLDN4 only.**

## What was (and was not) run

**CytoTRACE2 was not run.** The CytoTRACE2 R package (Kang et al. *Nat Commun* 2024) is not installed. The prompt allowed CytoTRACE2 **or** a potency score. This module uses the latter.

Stemness proxy (Gulati et al. *Science* 2020 CytoTRACE idea, not the R package):

1. Per cell, count genes with UMI > 0 (`n_genes`) and total UMI while streaming the public UMI TSV.
2. Among epithelial cells of that contrast, OLS: `n_genes ~ log1p(total_UMI)`.
3. Residual = observed − fitted. Rank-scale residuals to `[0, 1]`. Higher = more stem-like.

CLDN4: `log1p(CP10k)` from the same streamed totals.

## Data

**GSE131907** (Kim et al. *Nat Commun* 2020, PMID 32385277). Processed GEO UMI matrix (0.39 GB gzip) + author cell annotation. Skip the 2.86 GB log2TPM text and EGA FASTQ.

Epithelial = author `Cell_type == Epithelial cells`. Primary = tLung. Sensitivity = tLung + tL/B + mLN + mBrain + PE. Samples with <20 epithelial cells are excluded from sample-level tests.

No ICI / pathologic-response labels in this atlas.

## Tests

Unit = sample. Cell-level Spearman is stored and labeled exploratory.

1. Spearman of sample-mean epithelial CLDN4 vs CytoTRACE-like.
2. Extra: paired Wilcoxon of CytoTRACE-like in within-sample CLDN4 Q4 vs Q1 (min 8 cells/stratum).

Two-sided p. Honest n. Every test is in `results/stats.tsv`.
