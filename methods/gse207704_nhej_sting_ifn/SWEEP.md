# FINDING — GSE207704 replicate-level DE sweep

Additive counts analysis. The earlier FPKM scores stand. This file uses the eight deposited runs, not the collapsed FPKM table.

## Design that was actually fit

Kallisto 0.51.1, Ensembl 110 cDNA, single-end, fragment length 200 (sd 30). Transcript estimates were summed to gene symbols and rounded. n = 2 per genotype per line. Positive log2FC = higher in CLDN4−/− than in parental WT.

Methods: DESeq2 Wald (inmoose 0.9.1), edgeR quasi-likelihood F on a common dispersion with TMM factors, limma-trend (`eBayes(trend=True)` on log2 CPM). Robust empirical-Bayes squeezing is not implemented in this inmoose build, so both edgeR QL and limma use the standard squeeze. Genes with fewer than 2 samples at count ≥ 10 were dropped before each fit. FDR within a method is Benjamini-Hochberg on the genes tested in that contrast.

CLDN4 across fits: DESeq2 T47D -1.048 (p 2.41e-22), edgeR_QL T47D -1.054 (p 1.62e-10), limma_trend T47D -1.049 (p 5.17e-10), DESeq2 MCF7 -0.763 (p 2.07e-32), edgeR_QL MCF7 -0.762 (p 3.8e-11), limma_trend MCF7 -0.761 (p 1.48e-10), DESeq2 pooled -0.888 (p 4.63e-28), edgeR_QL pooled -0.900 (p 5.62e-07), limma_trend pooled -0.902 (p 4.99e-07).

Library sizes are in `tables/kallisto_sample_qc.tsv`. Gene counts: `tables/kallisto_gene_counts.tsv.gz` (30123 genes × 8 samples).

## Pre-declared thesis and the rule for “best”

Thesis, fixed before looking at these p-values: NHEJ sets go **down** after CLDN4 loss; STING, IFN, and APM sets go **up**. A row matches when the set’s mean log2FC has that sign.

Primary family: GSEA on the signed DE statistic, effect-size floor 0, two-sided gene-set permutation, three methods × three contrasts (T47D, MCF7, pooled `line + genotype`). The best primary match is the smallest nominal p inside that family among sign-matching rows. Sensitivity rows (floors 0.5 and 1, one-sided Wilcoxon, mean-statistic permutation) are in the same table and are not substituted for the primary result.

Global BH-FDR is computed across every test in the sweep, because the sweep itself is many looks.

## What is most significant

The strongest primary result opposes the thesis. IFN sets go down after CLDN4 loss. No specification in this sweep — method, contrast, test, or effect-size floor — puts a thesis-matching set at within-family FDR < 0.05. Of 1188 tests, 0 thesis-matching rows have within-family FDR < 0.05 and 0 have global FDR < 0.05.

Gene-set permutation p-values cannot be smaller than 1/1001 with 1000 permutations.

- Strongest primary GSEA row, any direction: DESeq2 / pooled / IFN / gsea / floor 0.0: mean log2FC -0.338, effect -1.787, nominal p 0.000999, within-family FDR 0.01099, global FDR 0.04929, n=50
- Strongest primary GSEA row whose sign matches the thesis: limma_trend / MCF7 / REACTOME_CLASS_I_MHC_PEPTIDE_LOADING / gsea / floor 0.0: mean log2FC +0.151, effect +1.534, nominal p 0.01698, within-family FDR 0.1758, global FDR 0.2211, n=26
- Strongest primary GSEA row whose sign opposes the thesis: DESeq2 / pooled / IFN / gsea / floor 0.0: mean log2FC -0.338, effect -1.787, nominal p 0.000999, within-family FDR 0.01099, global FDR 0.04929, n=50
- Strongest thesis-matching row anywhere in the sensitivity sweep: limma_trend / MCF7 / APM / gsea / floor 0.5: mean log2FC +0.705, effect +1.577, nominal p 0.01299, within-family FDR 0.07692, global FDR 0.2028, n=5

The sensitivity winner is an effect-size floor that keeps only genes with |log2FC| ≥ 0.5. It is not the primary result.

## Primary GSEA direction (floor 0)

Nine fits = DESeq2, edgeR QL, and limma-trend, each on T47D, MCF7, and the line-adjusted pool. The smallest nominal p and the smallest within-family FDR are taken over those nine whatever the sign. A small FDR on a down IFN set is evidence against the thesis.

| Set | thesis | fits with that sign | median mean log2FC | smallest nominal p (any sign) | best thesis-match p | smallest within-family FDR (any sign) |
|---|---|---:|---:|---:|---:|---:|
| NHEJ_CORE | down | 0/9 | +0.046 | 0.04096 | no row | 0.1126 |
| REACTOME_NONHOMOLOGOUS_END_JOINING | down | 9/9 | -0.013 | 0.3247 | 0.3247 | 0.3457 |
| STING_AXIS | up | 9/9 | +0.195 | 0.1039 | 0.1039 | 0.1525 |
| REACTOME_STING_MEDIATED_INDUCTION | up | 5/9 | +0.002 | 0.1389 | 0.2128 | 0.2546 |
| IFN | up | 0/9 | -0.285 | 0.000999 | no row | 0.01099 |
| APM | up | 6/9 | +0.029 | 0.06094 | 0.06094 | 0.1363 |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | up | 3/9 | -0.191 | 0.006993 | 0.2867 | 0.03846 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | up | 0/9 | -0.139 | 0.006993 | no row | 0.03846 |
| REACTOME_INTERFERON_ALPHA_BETA_SIGNALING | up | 0/9 | -0.229 | 0.02797 | no row | 0.07692 |
| REACTOME_INTERFERON_GAMMA_SIGNALING | up | 3/9 | -0.047 | 0.2018 | 0.2018 | 0.2747 |
| REACTOME_CLASS_I_MHC_PEPTIDE_LOADING | up | 9/9 | +0.089 | 0.01698 | 0.01698 | 0.1429 |

NHEJ_CORE mean log2FC is positive in every fit, so the 6-gene panel does not go down. Reactome NHEJ is slightly negative in every fit and is not significant. Reactome STING-mediated induction contains PRKDC, XRCC5, XRCC6, and MRE11, so it is not a pure CGAS–STING1–TBK1–IRF3 test. STING1 itself is present in the cDNA counts but the counts are small (see the QC table); the positive log2FC is not a significant gene-level result.

IFN, Reactome IFN-α/β, and Hallmark IFN-γ are negative in every primary fit. Hallmark IFN-α is negative in the fits that reach nominal p < 0.05.

## Named genes, DESeq2

### T47D

| Gene | log2FC | nominal p | FDR |
|---|---:|---:|---:|
| CLDN4 | -1.048 | 2.41e-22 | 2.38e-20 |
| TACSTD2 | -0.642 | 2.89e-15 | 1.62e-13 |
| PRKDC | +0.258 | 2.92e-06 | 4.47e-05 |
| LIG4 | +0.125 | 0.584 | 0.806 |
| XRCC4 | -0.098 | 0.554 | 0.786 |
| XRCC5 | +0.168 | 0.000628 | 0.00507 |
| XRCC6 | -0.195 | 0.000124 | 0.00123 |
| NHEJ1 | +0.046 | 0.768 | 0.911 |
| CGAS | +0.365 | 0.00551 | 0.0311 |
| STING1 | +0.429 | 0.606 | 0.82 |
| TBK1 | +0.028 | 0.718 | 0.886 |
| IRF3 | -0.075 | 0.517 | 0.761 |

### MCF7

| Gene | log2FC | nominal p | FDR |
|---|---:|---:|---:|
| CLDN4 | -0.763 | 2.07e-32 | 1.68e-30 |
| TACSTD2 | -1.042 | 5.85e-48 | 9e-46 |
| PRKDC | -0.157 | 3.06e-05 | 0.000237 |
| LIG4 | -0.026 | 0.887 | 0.954 |
| XRCC4 | +0.294 | 0.0339 | 0.106 |
| XRCC5 | +0.162 | 8.95e-08 | 1.08e-06 |
| XRCC6 | +0.095 | 0.0035 | 0.0162 |
| NHEJ1 | -0.126 | 0.417 | 0.651 |
| CGAS | +0.191 | 0.165 | 0.354 |
| STING1 | +0.403 | 0.527 | 0.745 |
| TBK1 | +0.050 | 0.355 | 0.592 |
| IRF3 | +0.186 | 0.0498 | 0.145 |

edgeR and limma log2FC values are in `tables/gene_de_replicate.tsv` and in `figures/fig_replicate_de_named_genes.png`. Primary GSEA NES values are in `figures/fig_de_sweep_gsea_nes.png`. Set mean log2FC is in `figures/fig_de_sweep_mean_log2fc.png`. The full grid is `tables/de_sweep.tsv`.
