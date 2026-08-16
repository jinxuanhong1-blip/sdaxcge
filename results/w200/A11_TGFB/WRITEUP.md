# A11 TGF-β genes vs TACSTD2-high public lung

> Honest report. Numbers are computed from public matrices.
> TCGA patients were not ICI-treated. Correlations are not causation.
> CLDN4 is a junction positive control, not a TGF-β gene.
> CD8A is immune context only.

## TL;DR

**Honest answer: not a TGF-β signaling class effect. SMAD3 yes; TGFB1 modest and LUAD-stronger; HALLMARK / F-TBRS are LUAD-only weak scores; receptors and the other SMADs no.**

The pre-specified ligand rule is `LIGAND_CLASS_EFFECT_SUPPORTED` (2/3 TGF-β ligands are purity+histology-adjusted positive with TACSTD2 and none are significantly negative.) That rule is too easy to trip: TGFB2 is only WEAK_POSITIVE (pooled partial ρ = 0.11) and **fails the TACSTD2 Q4 vs Q1 test** (Δmedian = +0.34, p = 0.056). In LUSC, TGFB2 is negative (partial ρ = −0.11; Q4 vs Q1 Δ = −0.36, p = 0.004). The pathway rule is `PATHWAY_CLAIM_SUPPORTED` because HALLMARK is WEAK_POSITIVE (ρ = 0.15) and TGFB1 is not negative. That score is **LUAD-only** (LUAD 0.32, LUSC 0.01, LUSC Q4 vs Q1 p = 0.89). The named A11 question is TACSTD2-**high** public lung, not a weak LUAD-driven pooled ρ.

Primary cohort: TCGA LUAD n=516 + LUSC n=501 primary tumors (one sample/patient); 996 / 1017 have ABSOLUTE purity and enter the primary partial correlation. HALLMARK score uses 54/54 genes; F-TBRS score uses 25 genes.

## Pre-specified design

- Ligands (class rule): TGFB1, TGFB2, TGFB3.
- Receptors and SMADs reported separately: TGFBR1/2/3, SMAD2/3/4/6/7.
- Pathway score: z-mean of MSigDB HALLMARK_TGF_BETA_SIGNALING genes present in the extract.
- Sensitivity score: z-mean of Mariathasan 2018 F-TBRS (fibroblast TGF-β response) genes present.
- Expression: Xena GDC STAR TPM, log2(TPM+1), GENCODE v36.
- Primary statistic: rank-partial Spearman of TACSTD2 vs each gene/score; covariates = ABSOLUTE purity + histology in the pooled analysis.
- TACSTD2-high: Q4 vs Q1 Mann-Whitney U (median split is sensitivity).
- FDR: BH within each cohort across the 3 ligands (family) and across all partners (descriptive).
- Pathway claim requires HALLMARK score positive and TGFB1 not negative.
- CLDN4 = junction control. CD8A = immune context.

## Primary result — pooled partial Spearman (purity + histology)

| gene | role | partial ρ | partial p | FDR | label | LUAD partial ρ | LUSC partial ρ |
|---|---|---:|---:|---:|---|---:|---:|
| TGFB1 | primary_ligand | 0.215 | 6.85e-12 | 2.05e-11 | ASSOCIATED_POSITIVE | 0.291 | 0.168 |
| TGFB2 | primary_ligand | 0.111 | 4.55e-04 | 6.83e-04 | WEAK_POSITIVE | 0.354 | -0.114 |
| TGFB3 | primary_ligand | -0.080 | 1.20e-02 | 1.20e-02 | NULL | 0.035 | -0.179 |
| TGFBR1 | receptor | -0.061 | 5.40e-02 | 6.75e-02 | NULL | 0.023 | -0.133 |
| TGFBR2 | receptor | 0.056 | 7.83e-02 | 9.04e-02 | NULL | 0.205 | -0.089 |
| TGFBR3 | receptor | 0.037 | 2.39e-01 | 2.56e-01 | NULL | 0.119 | -0.002 |
| SMAD2 | smad | -0.067 | 3.39e-02 | 4.63e-02 | NULL | -0.011 | -0.119 |
| SMAD3 | smad | 0.340 | 2.56e-28 | 1.92e-27 | ASSOCIATED_POSITIVE | 0.342 | 0.341 |
| SMAD4 | smad | -0.129 | 4.33e-05 | 1.08e-04 | WEAK_NEGATIVE | -0.031 | -0.216 |
| SMAD6 | smad | -0.015 | 6.27e-01 | 6.27e-01 | NULL | 0.065 | -0.084 |
| SMAD7 | smad | -0.081 | 1.02e-02 | 1.70e-02 | NULL | 0.080 | -0.208 |
| HALLMARK_TGFB_ZMEAN | pathway_score | 0.153 | 1.21e-06 | 3.63e-06 | WEAK_POSITIVE | 0.316 | 0.012 |
| FTBRS_ZMEAN | pathway_score | 0.107 | 7.03e-04 | 1.32e-03 | WEAK_POSITIVE | 0.261 | -0.027 |
| CLDN4 | positive_control | 0.436 | 2.07e-47 | 3.11e-46 | ASSOCIATED_POSITIVE | 0.525 | 0.388 |
| CD8A | immune_context | -0.206 | 5.87e-11 | 2.20e-10 | ASSOCIATED_NEGATIVE | -0.085 | -0.295 |

## TACSTD2-high vs TACSTD2-low (pooled Q4 vs Q1)

| gene | median Q4 | median Q1 | Δmedian | MWU p |
|---|---:|---:|---:|---:|
| TGFB1 | 6.583 | 6.050 | 0.533 | 1.35e-11 |
| TGFB2 | 2.300 | 1.957 | 0.344 | 5.61e-02 |
| TGFB3 | 3.709 | 3.853 | -0.144 | 1.88e-01 |
| TGFBR1 | 4.806 | 4.885 | -0.079 | 2.14e-02 |
| TGFBR2 | 5.464 | 5.556 | -0.092 | 3.36e-02 |
| TGFBR3 | 2.825 | 2.462 | 0.363 | 4.67e-03 |
| SMAD2 | 3.202 | 3.210 | -0.008 | 8.29e-01 |
| SMAD3 | 5.562 | 4.489 | 1.073 | 1.47e-30 |
| SMAD4 | 3.853 | 3.977 | -0.124 | 4.69e-03 |
| SMAD6 | 1.993 | 2.069 | -0.076 | 2.49e-02 |
| SMAD7 | 3.882 | 4.167 | -0.285 | 6.89e-06 |
| HALLMARK_TGFB_ZMEAN | 0.091 | -0.080 | 0.171 | 6.93e-06 |
| FTBRS_ZMEAN | 0.036 | -0.152 | 0.188 | 6.30e-03 |
| CLDN4 | 8.029 | 7.204 | 0.825 | 3.61e-18 |
| CD8A | 3.362 | 4.118 | -0.756 | 8.13e-11 |

## Honest reading of each object

- **SMAD3**: the only robust, histology-replicated TGF-β-cassette partner. Pooled partial ρ = 0.34 (LUAD 0.34, LUSC 0.34). Q4 vs Q1 Δmedian = +1.07 (p = 1e-30). DepMap all-lung ρ = 0.48; NSCLC-only drops to 0.17 (FDR-null). Strength is in the same range as the CLDN4 control (partial ρ = 0.44). This is a SMAD3 association, not proof of TGF-β pathway activity.
- **TGFB1**: real modest co-expression, stronger in LUAD (0.29) than LUSC (0.17). Pooled Q4 vs Q1 Δmedian = +0.53 (p = 1e-11). Weaker than CLDN4 and SMAD3. DepMap NSCLC is **negative** (ρ = −0.23). Do not treat TGFB1 as a CLDN4-like TACSTD2-high partner.
- **TGFB2**: LUAD-only (0.35); LUSC negative (−0.11). Pooled Q4 vs Q1 is null (p = 0.056). Do not call this a TACSTD2-high TGF-β ligand.
- **TGFB3**: pooled |ρ| < 0.10. LUSC Q4 vs Q1 is lower, not higher (Δ = −0.54, p = 8e-4). Not a TACSTD2-high partner.
- **TGFBR1 / TGFBR2 / TGFBR3**: pooled null. Q4 vs Q1 is slightly *lower* for TGFBR1 and TGFBR2. No receptor class effect.
- **SMAD2 / SMAD6 / SMAD7**: pooled null. SMAD7 is *lower* in TACSTD2 Q4 (Δ = −0.29, p = 7e-6), the opposite of a TGF-β-feedback-on signature.
- **SMAD4**: weakly negative (pooled ρ = −0.13), LUSC-driven (−0.22). Opposite of the claim.
- **HALLMARK TGF-β score**: pooled WEAK_POSITIVE (0.15) is LUAD-only (0.32 vs LUSC 0.01). Dropping SMAD3 from the score barely changes the pooled ρ (0.15 → 0.15). LUSC Q4 vs Q1 p = 0.89. Not a lung-wide TGF-β program.
- **F-TBRS**: same pattern (LUAD 0.26, LUSC −0.03). A fibroblast response signature in bulk RNA can track LUAD stroma, not malignant-cell TGF-β output.
- **CLDN4 control**: recovered (partial ρ = 0.44). The pipeline can see a junction association; that does not make TGFB2/3 or the receptors positive.
- **CD8A**: pooled partial ρ = −0.21 (LUSC −0.30, LUAD −0.09). TACSTD2-high tumors are immune-colder on this marker. That is context, not TGF-β evidence.

## What this does **not** show

- It does not show that TGF-β causes TACSTD2-high tumors, or the reverse.
- It does not show ICI response, ADC response, or protein-level TGF-β activity.
- Bulk TGF-β ligand mRNA is often stromal. A purity-adjusted null or negative is the more honest tumor-cell test than a raw Spearman.
- HALLMARK_TGF_BETA_SIGNALING is a mixed transcriptional set, not phospho-SMAD activity.
- F-TBRS is a fibroblast response signature; in bulk tumor it can track stroma, not malignant-cell TGF-β output.
- Histology can disagree. If LUAD and LUSC labels differ, the pooled number is not a license to ignore the split.
- DepMap lung lines are models, not tumors; they are sensitivity only.

## Sensitivity — DepMap 24Q4 lung cell lines (no purity adjustment)

| cohort | gene | Spearman ρ | p | n | label |
|---|---|---:|---:|---:|---|
| depmap_lung_cell_lines | TGFB1 | 0.168 | 1.41e-02 | 214 | WEAK_POSITIVE |
| depmap_lung_cell_lines | TGFB2 | 0.231 | 6.74e-04 | 214 | ASSOCIATED_POSITIVE |
| depmap_lung_cell_lines | TGFB3 | -0.159 | 2.01e-02 | 214 | WEAK_NEGATIVE |
| depmap_lung_cell_lines | TGFBR1 | -0.362 | 5.07e-08 | 214 | ASSOCIATED_NEGATIVE |
| depmap_lung_cell_lines | TGFBR2 | 0.468 | 4.64e-13 | 214 | ASSOCIATED_POSITIVE |
| depmap_lung_cell_lines | TGFBR3 | -0.100 | 1.44e-01 | 214 | NULL |
| depmap_lung_cell_lines | SMAD2 | -0.267 | 7.47e-05 | 214 | ASSOCIATED_NEGATIVE |
| depmap_lung_cell_lines | SMAD3 | 0.481 | 9.22e-14 | 214 | ASSOCIATED_POSITIVE |
| depmap_lung_cell_lines | SMAD4 | -0.417 | 2.08e-10 | 214 | ASSOCIATED_NEGATIVE |
| depmap_lung_cell_lines | SMAD6 | 0.001 | 9.93e-01 | 214 | NULL |
| depmap_lung_cell_lines | SMAD7 | -0.091 | 1.83e-01 | 214 | NULL |
| depmap_lung_cell_lines | CLDN4 | 0.607 | 5.64e-23 | 214 | ASSOCIATED_POSITIVE |
| depmap_lung_cell_lines | CD8A | -0.229 | 7.56e-04 | 214 | ASSOCIATED_NEGATIVE |
| depmap_lung_cell_lines | HALLMARK_TGFB_ZMEAN | 0.233 | 5.92e-04 | 214 | ASSOCIATED_POSITIVE |
| depmap_lung_cell_lines | FTBRS_ZMEAN | 0.352 | 1.19e-07 | 214 | ASSOCIATED_POSITIVE |
| depmap_NSCLC_cell_lines | TGFB1 | -0.235 | 4.81e-03 | 143 | ASSOCIATED_NEGATIVE |
| depmap_NSCLC_cell_lines | TGFB2 | 0.006 | 9.43e-01 | 143 | NULL |
| depmap_NSCLC_cell_lines | TGFB3 | -0.149 | 7.56e-02 | 143 | NULL |
| depmap_NSCLC_cell_lines | TGFBR1 | -0.261 | 1.67e-03 | 143 | ASSOCIATED_NEGATIVE |
| depmap_NSCLC_cell_lines | TGFBR2 | 0.122 | 1.47e-01 | 143 | NULL |
| depmap_NSCLC_cell_lines | TGFBR3 | -0.171 | 4.14e-02 | 143 | NULL |
| depmap_NSCLC_cell_lines | SMAD2 | -0.255 | 2.15e-03 | 143 | ASSOCIATED_NEGATIVE |
| depmap_NSCLC_cell_lines | SMAD3 | 0.169 | 4.32e-02 | 143 | NULL |
| depmap_NSCLC_cell_lines | SMAD4 | -0.207 | 1.29e-02 | 143 | ASSOCIATED_NEGATIVE |
| depmap_NSCLC_cell_lines | SMAD6 | -0.156 | 6.24e-02 | 143 | NULL |
| depmap_NSCLC_cell_lines | SMAD7 | -0.138 | 1.00e-01 | 143 | NULL |
| depmap_NSCLC_cell_lines | CLDN4 | 0.667 | 8.81e-20 | 143 | ASSOCIATED_POSITIVE |
| depmap_NSCLC_cell_lines | CD8A | -0.064 | 4.46e-01 | 143 | NULL |
| depmap_NSCLC_cell_lines | HALLMARK_TGFB_ZMEAN | 0.069 | 4.11e-01 | 143 | NULL |
| depmap_NSCLC_cell_lines | FTBRS_ZMEAN | 0.043 | 6.09e-01 | 143 | NULL |

## Tumor vs adjacent normal (context, not the A11 claim)

| cohort | gene | median tumor | median normal | Δ | unpaired p |
|---|---|---:|---:|---:|---:|
| LUAD | TACSTD2 | 9.013 | 8.637 | 0.376 | 7.97e-03 |
| LUAD | TGFB1 | 6.232 | 6.750 | -0.518 | 3.47e-09 |
| LUAD | TGFB2 | 2.368 | 4.115 | -1.748 | 2.51e-15 |
| LUAD | TGFB3 | 3.581 | 3.780 | -0.199 | 1.31e-01 |
| LUAD | TGFBR1 | 4.880 | 5.129 | -0.249 | 4.04e-04 |
| LUAD | TGFBR2 | 6.331 | 8.176 | -1.845 | 2.21e-33 |
| LUAD | TGFBR3 | 2.378 | 5.388 | -3.010 | 2.21e-33 |
| LUAD | SMAD2 | 3.075 | 3.421 | -0.347 | 9.54e-08 |
| LUAD | SMAD3 | 4.680 | 4.825 | -0.145 | 2.34e-01 |
| LUAD | SMAD4 | 3.911 | 4.403 | -0.492 | 1.68e-13 |
| LUAD | SMAD6 | 2.317 | 4.638 | -2.320 | 2.11e-29 |
| LUAD | SMAD7 | 4.421 | 5.436 | -1.016 | 5.70e-23 |
| LUAD | CLDN4 | 8.181 | 7.016 | 1.164 | 6.59e-18 |
| LUAD | CD8A | 3.847 | 4.196 | -0.350 | 1.70e-02 |
| LUAD | HALLMARK_TGFB_ZMEAN | 0.044 | 0.013 | 0.031 | 8.26e-01 |
| LUAD | FTBRS_ZMEAN | -0.031 | 0.105 | -0.136 | 8.93e-01 |
| LUSC | TACSTD2 | 9.575 | 8.699 | 0.876 | 4.83e-09 |
| LUSC | TGFB1 | 6.548 | 6.846 | -0.298 | 2.17e-03 |
| LUSC | TGFB2 | 2.101 | 3.932 | -1.831 | 3.23e-17 |
| LUSC | TGFB3 | 4.012 | 4.079 | -0.066 | 7.62e-01 |
| LUSC | TGFBR1 | 4.900 | 5.067 | -0.166 | 4.55e-02 |
| LUSC | TGFBR2 | 5.102 | 8.126 | -3.024 | 7.55e-32 |
| LUSC | TGFBR3 | 2.967 | 5.099 | -2.132 | 1.29e-26 |
| LUSC | SMAD2 | 3.365 | 3.297 | 0.068 | 4.74e-01 |
| LUSC | SMAD3 | 5.442 | 4.893 | 0.549 | 9.34e-09 |
| LUSC | SMAD4 | 3.963 | 4.244 | -0.281 | 7.16e-05 |
| LUSC | SMAD6 | 1.845 | 4.248 | -2.403 | 2.79e-28 |
| LUSC | SMAD7 | 3.857 | 5.543 | -1.686 | 5.64e-27 |
| LUSC | CLDN4 | 7.206 | 6.808 | 0.398 | 1.36e-02 |
| LUSC | CD8A | 3.670 | 4.523 | -0.853 | 2.99e-07 |
| LUSC | HALLMARK_TGFB_ZMEAN | 0.037 | 0.071 | -0.034 | 9.72e-01 |
| LUSC | FTBRS_ZMEAN | -0.003 | -0.143 | 0.140 | 8.45e-01 |

## Files

- `correlations.csv` — marginal and partial Spearman
- `tacstd2_high_vs_low.csv` — Q4 vs Q1 and median split
- `tumor_vs_normal.csv` — adjacent-normal context
- `summary.json` — machine-readable verdict
- `fig1_partial_rho_heatmap.png` / `fig2_scatter_*.png` / `fig3_*boxplots.png`

