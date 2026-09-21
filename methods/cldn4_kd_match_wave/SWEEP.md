# Cancer KD sweep: IFN, APM, and NHEJ after CLDN4 loss

Question: does any honest re-scoring of GSE207704 or GSE22493 show IFN, APM, or NHEJ going up after CLDN4 loss?

A concordant hit requires two things at once. CLDN4 log2FC is negative in that same scoring, and the panel is called up. For gene-median scores the up call is the Claim C4 rule (median above 0 and one-sided Wilcoxon P ≤ 0.05, or the sign test when fewer than 5 genes are observed). Pathway and array-unit rows use the test named in the method.

GSE207704 in GEO is a four-column Cufflinks FPKM table (condition means). Replicate FASTQs exist (two per genotype). Kallisto 0.52.0 on Ensembl 111 cDNA, single-end fragment length 200 (sd 20). Estimated counts summed to gene_symbol. Median-of-ratios size factors. log2((mean KO + 1) / (mean WT + 1)), n=2 vs 2. Pseudoalignment of the scored runs: 37.0–42.1%. A row is a concordant hit only when CLDN4 log2FC is negative in this same scoring.

## Primary estimates still match the locked table

GSE207704 sum of loci, pseudocount 0.5, mean FPKM at least 0.5, and GSE22493 median-probe then mean across the three arrays, reproduce the locked IFN/ISG and MHC-I/APM medians. Those calls stay discordant or unsigned, and every primary median is negative.

## What still stands

The combined cancer contrasts do not open IFN or APM. Concordant hits on those combined contrasts, excluding single arrays and max-probe rules: 0. Primary IFN/ISG and MHC-I/APM medians stay negative in T47D, MCF-7, and SKOV-3, and they still match the locked table. Prerank GSEA, the three-array t-test, and Stouffer's combination of the gene-wise t-tests do not call those panels up. Raw ScanArray log2(KD/control), with CLDN4 down, calls IFN/ISG, MHC-I/APM, and NHEJ_CORE down.

NHEJ_CORE is not a second version of the IFN split. T47D: median +0.092, 5/7 up, Wilcoxon P(greater) = 0.109, call weak_up. MCF-7: median -0.046, 3/7 up, Wilcoxon P(greater) = 0.812, call weak_down. SKOV-3: median -0.062, 2/5 up, Wilcoxon P(greater) = 0.688, call weak_down. Lung baseline: median +0.049, 4/7 up, Wilcoxon P(greater) = 0.344, call weak_up.

## The one array that leans the other way

SKOV-3 is three two-colour arrays. GSM558701 IFN/ISG: median +0.683, 14/24 up, Wilcoxon P(greater) = 0.0604, call weak_up. GSM558701 MHC-I/APM: median +0.964, 10/14 up, Wilcoxon P(greater) = 0.0765, call weak_up. GSM558701 Hallmark IFN-α: median +0.345, 39/68 up, Wilcoxon P(greater) = 0.0178, call up. GSM558701 Hallmark IFN-γ: median +0.508, 92/145 up, Wilcoxon P(greater) = 0.000104, call up. GSM558700 has no CLDN4 measurement and keeps IFN and APM negative. GSM558702 has CLDN4 down and keeps IFN and APM negative. One array out of three is not the combined SKOV-3 contrast, and it does not move T47D or MCF-7.

## Replicate quantification

Kallisto 0.52.0 on Ensembl 111 cDNA, single-end fragment length 200 (sd 20). Estimated counts summed to gene_symbol. Median-of-ratios size factors. log2((mean KO + 1) / (mean WT + 1)), n=2 vs 2. Pseudoalignment of the scored runs: 37.0–42.1%. A row is a concordant hit only when CLDN4 log2FC is negative in this same scoring.

T47D KO vs WT: CLDN4 log2FC -1.047. IFN/ISG: median -0.251, 13/35 up, Wilcoxon P(greater) = 0.987, call down. MHC-I/APM: median -0.232, 2/17 up, Wilcoxon P(greater) = 1, call down. NHEJ_CORE: median +0.076, 5/8 up, Wilcoxon P(greater) = 0.422, call weak_up. Hallmark IFN-α: median -0.175, 29/89 up, Wilcoxon P(greater) = 0.999, call down. Hallmark IFN-γ: median -0.095, 69/171 up, Wilcoxon P(greater) = 1, call down.

MCF-7 KO vs WT: CLDN4 log2FC -0.760. IFN/ISG: median -0.205, 11/34 up, Wilcoxon P(greater) = 0.993, call down. MHC-I/APM: median -0.060, 7/17 up, Wilcoxon P(greater) = 0.445, call weak_down. NHEJ_CORE: median +0.037, 4/8 up, Wilcoxon P(greater) = 0.23, call weak_up. Hallmark IFN-α: median +0.026, 44/86 up, Wilcoxon P(greater) = 0.496, call weak_up. Hallmark IFN-γ: median -0.044, 82/177 up, Wilcoxon P(greater) = 0.811, call weak_down.

Prerank GSEA on those kallisto rankings does not call IFN/ISG, MHC-I/APM, or NHEJ_CORE up.

## Concordant hits

Single-array hits: 4. Other non-optimistic hits: 0. Optimistic max-probe or max-locus hits: 5.

| Family | Method | Contrast | Panel | Estimate | Call |
|---|---|---|---|---:|---|
| optimistic | fpkm_max_logfc_pc0.5_min0.5 | T47D KO vs WT | REACTOME_NHEJ | +0.087 | up |
| optimistic | deposited_max_probe | SKOV-3 KD | NHEJ_EXT | +0.601 | up |
| optimistic | deposited_max_probe | SKOV-3 KD | HALLMARK_IFNA | +0.090 | up |
| optimistic | deposited_max_probe | SKOV-3 KD | HALLMARK_IFNG | +0.072 | up |
| optimistic | deposited_max_probe+stouffer_t | SKOV-3 KD | MHC_APM | +0.453 | up |
| single_array | deposited_median_probe_GSM558701 | SKOV-3 KD | NHEJ_EXT | +0.343 | up |
| single_array | deposited_median_probe_GSM558701 | SKOV-3 KD | HALLMARK_IFNA | +0.345 | up |
| single_array | deposited_median_probe_GSM558701 | SKOV-3 KD | HALLMARK_IFNG | +0.508 | up |
| single_array | deposited_median_probe_GSM558701 | SKOV-3 KD | REACTOME_NHEJ | +0.275 | up |

Optimistic rows pick the highest probe or locus per gene before taking the panel median. They are a hunt, not the primary estimate.

## How to read the figure

Filled points are the primary medians. The bar is the min-to-max of the other non-optimistic sensitivities in which CLDN4 fell (pseudocount, filter, mean or median collapse, leave-one-array-out, single arrays, raw ScanArray log2(KD/control), outlier trim, and kallisto when that line was quantified). The open circle is the highest median from a max-probe or max-locus rule that also has CLDN4 down. Lung IFN/APM points are the locked baseline and were not refit. Lung NHEJ, when drawn, is scored on the author edgeR table.

Raw ScanArray uses GEO channel 2 as knockdown and channel 1 as the overexpression control, log2((Ch2 background-subtracted median + 1) / (Ch1 + 1)). A scoring is left out of the hit list and out of the sensitivity bar when CLDN4 does not fall.

Reactome NHEJ contains histone genes, so it is a sensitivity panel. NHEJ_CORE is XRCC5, XRCC6, PRKDC, LIG4, XRCC4, NHEJ1, DCLRE1C, and PAXX.

Tables: `sweep/sweep_results.tsv`, `sweep/concordant_hits.tsv`.

Reproduce: `python3 scripts/cldn4_kd_match_wave/sweep_cancer_kd.py`
