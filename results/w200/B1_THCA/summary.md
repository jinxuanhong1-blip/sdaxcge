# B1 analog — TCGA-THCA: TACSTD2 × CLDN4 surface co-expression rank

**Question:** Is `CLDN4` the top co-expression partner of `TACSTD2` among surface genes in TCGA-THCA?

**Answer:** NO — CLDN4 is NOT the single top surface-gene co-expression partner of TACSTD2.

- Cohort: TCGA-THCA primary tumors, n=505 (sample-type 01, one vial per patient).
- Expression: UCSC Xena GDC hub `TCGA-THCA.star_tpm` = log2(TPM+1), GENCODE v36.
- Surface-gene universe: 2674 ranked genes (2799 unique SURFY symbols in the workbook; `TACSTD2` and zero-variance genes removed from the co-expression ranking).
- Primary metric: Spearman correlation (honest full-universe rank; no pre-selected neighbourhood).

## Co-expression rank (claim B1 analog)

| metric | value |
| --- | --- |
| CLDN4 Spearman rank | **#30 of 2674** (98.915th percentile) |
| CLDN4 Spearman rho | 0.676 (FDR q=7.91e-67) |
| CLDN4 Pearson rank | #49 of 2674 |
| CLDN4 Pearson rho | 0.680 |
| Actual #1 (Spearman) | DCSTAMP (rho=0.826) |
| TACSTD2–CLDN4 Spearman (bootstrap 95% CI) | 0.676 [0.62, 0.72] |

## Surface-expression rank (median log2(TPM+1))

This is the original “surface rank” readout: how highly each gene is expressed among the surfaceome, independent of co-expression.

| gene | expression rank |
| --- | --- |
| TACSTD2 | #14 of 2686 (median log2(TPM+1) = 9.39) |
| CLDN4 | #23 of 2686 (median log2(TPM+1) = 8.53) |

Median log2(TPM+1) in the analysis cohort: TACSTD2 = 9.39, CLDN4 = 8.53. Both are highly expressed in THCA.

## Sensitivity

| analysis | n | Spearman rho | p |
| --- | --- | --- | --- |
| all_primary_tumors | 505 | 0.676 [0.62, 0.72] | 8.9e-69 |
| classic_papillary_8260_3 | 355 | 0.594 [0.52, 0.66] | 3.4e-35 |
| follicular_variant_8340_3 | 105 | 0.586 [0.43, 0.71] | 5.0e-11 |

Purity (ABSOLUTE, DNA-based):

| quantity | n | rho | p |
| --- | --- | --- | --- |
| spearman_unadjusted_purity_subset | 463 | 0.677 | 2.9e-63 |
| spearman_partial_given_ABSOLUTE_purity | 463 | 0.688 | 3.3e-66 |
| spearman_TACSTD2_vs_purity | 463 | -0.085 | 6.9e-02 |
| spearman_CLDN4_vs_purity | 463 | 0.081 | 8.3e-02 |

## Top 15 surface-gene partners of TACSTD2 (Spearman)

| rank | gene | spearman_rho | pearson_rho | FDR q |
| --- | --- | --- | --- | --- |
| 1 | DCSTAMP | 0.826 | 0.873 | 1.08e-123 |
| 2 | ERBB3 | 0.779 | 0.827 | 6.88e-101 |
| 3 | TGFBR1 | 0.763 | 0.754 | 2.60e-94 |
| 4 | MET | 0.761 | 0.803 | 1.35e-93 |
| 5 | GJB3 | 0.759 | 0.842 | 4.49e-93 |
| 6 | ITGA3 | 0.759 | 0.782 | 4.82e-93 |
| 7 | TMPRSS6 | 0.754 | 0.836 | 4.06e-91 |
| 8 | QSOX1 | 0.753 | 0.769 | 4.49e-91 |
| 9 | GABRB2 | 0.742 | 0.839 | 8.06e-87 |
| 10 | CRLF2 | 0.739 | 0.682 | 4.71e-86 |
| 11 | MXRA8 | 0.738 | 0.774 | 1.39e-85 |
| 12 | ITGB8 | 0.737 | 0.795 | 2.12e-85 |
| 13 | SLC34A2 | 0.727 | 0.864 | 5.19e-82 |
| 14 | GJB4 | 0.715 | 0.628 | 8.77e-78 |
| 15 | LY6G6C | 0.712 | 0.714 | 3.74e-77 |

## Honest caveats

1. **This is not a pan-cancer test of claim B1.** Claim B1 is stated for a TCGA pan-cancer surfaceome ranking. THCA is one analog cohort. A top (or near-top) rank here does not prove the pan-cancer claim; a mid-pack rank here does not disprove it outside thyroid.
2. **Bulk co-expression is not proof of a protein complex or of co-occurrence on the same cell.** THCA is epithelium-rich, so two epithelial surface genes can correlate through shared epithelial content. ABSOLUTE purity adjustment is reported; DNA purity is an imperfect proxy for the epithelial mRNA fraction.
3. **mRNA ≠ protein.** ADC-target relevance needs IHC / protein. This analysis cannot replace that.
4. **Xena GDC `star_tpm` is treated as log2(TPM+1).** The matrix floor is 0.0 (not −9.97), which is the log2(TPM+1) signature. Ranks are invariant to any strictly increasing transform, so a mistaken invert-to-linear-TPM step would not change ranks, only the reported TPM numbers.
5. **TCGA-THCA is almost all papillary thyroid carcinoma** (classic 8260/3 + follicular variant 8340/3). Follicular, poorly differentiated, and anaplastic carcinomas are essentially absent. Nothing here speaks to those histologies, or to ICI / ADC response. Pooled rho is higher than either histology subset (see Sensitivity): between-subtype mean shift inflates the headline number. Quote the within-histology figures if the claim is about a within-tumor-type relationship.
6. Surfaceome membership is the 2018 SURFY in-silico set (2,886 proteins). Gene-symbol matching misses retired / renamed symbols; unmatched genes are dropped, not imputed.

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/w200/B1_THCA/download_data.py
python3 scripts/w200/B1_THCA/run_analysis.py
```

Outputs: `coexpression_TACSTD2_surfaceome.csv`, `top200.csv`, `surface_expression_rank.csv`, `coexpression_main.csv`, `purity_adjusted.csv`, `summary.json`, `summary.md`, `scatter_TACSTD2_vs_CLDN4.png`, `top_partners_TACSTD2.png`.
