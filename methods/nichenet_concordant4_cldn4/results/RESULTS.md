# RESULTS — NicheNet ligand activity, concordant-4

CLDN4-high vs CLDN4-low **malignant senders** to **T/NK receivers**.
Cohorts are only the locked concordant four: GSE123902, GSE131907, GSE205335, GSE189357.
GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526, and GSE207422 are not in this run.
TACSTD2 is not a gate. There is no dual-high call.

Activity is `predict_ligand_activities` from saeyslab/nichenetr
`66f90d5eeafef280b2b2f339b3fd70ffec1781dd` (2026-06-01), executed in R:
`get_single_ligand_importances` → `evaluate_target_prediction` →
`classification_evaluation_continuous_pred` (ROCR curves, `caTools::trapz` AUPR).
The prior is NicheNet v2 `ligand_target_matrix_nsga2r_final.rds` (Zenodo 7074291).
Official rank is **AUPR corrected** (AUPR minus the positive-class prevalence). Pearson and AUROC are reported beside it.

## What the two arms say

On 59 eligible units, CLDN4-high malignant cells express more of the barrier ligands than CLDN4-low malignant cells from the same unit. The core set F11R, NECTIN2, CDH1, LGALS9 has mean Δ **+0.259** (Wilcoxon p = 2.65×10⁻¹¹). All four cohorts are positive, and dropping any one cohort leaves the sign unchanged. LGALS9 is the weak member of that set (mean Δ +0.052).

The core IFN/recruit ligands do not rise in CLDN4-high cells. CXCL9, CXCL10, and IFNG are too sparse to be potential ligands. CCL5 is lower in CLDN4-high cells in three cohorts and higher in GSE131907, so the family random-effects p is 0.073. The extended IFN/recruit list (CXCL16, ICAM1, IL15, and the rest) is higher, not lower, in CLDN4-high cells. HLA-A/B/C are higher in CLDN4-high cells, in the same direction as the earlier CellChat HLA–CD8 result, except in GSE189357.

NicheNet activity does not rank those barrier ligands as the regulators of T/NK programs. Against 286 potential ligands, the four core barrier ligands are no higher than a random draw on the a priori IFN, cytotoxicity, or exhaustion genesets (one-sided p = 0.90, 0.82, 0.46). CD274, not F11R, is the ligand whose prior targets best match the exhaustion genes (rank 1 of 286), and its sender Δ is only +0.010. Absolute AUPR-corrected values are small.

## Honest n

| cohort | units loaded | eligible quartile (n_mal≥40, n_tnk≥20, CLDN4 high>low, ≥25% CLDN4+) |
|---|---:|---:|
| GSE123902 | 13 | 10 |
| GSE131907 | 21 | 19 |
| GSE205335 | 22 | 21 |
| GSE189357 | 9 | 9 |
| **all** | **65** | **59** |

The unit is the locked patient or sample (GSE131907 stays at the sample id used in the concordant-4 inventory, not a further patient collapse).
Cells are not the independent unit. Primary sender split is **within-unit**: malignant cells ranked by CLDN4 `log1p(CP10k)`, top quartile vs bottom quartile (`rank` ties broken by first occurrence, same rule as the CellChat concordant-4 run).
A unit enters the quartile test only when CLDN4 itself is higher in the top quartile than the bottom quartile, and at least 25% of malignant cells express CLDN4. Below 25% positive, zero-ties would fill the high arm. Those units are kept for the percent-positive split.

## Sender expression — within-patient Q4 minus Q1

For each eligible unit, the family score is the mean ligand Δ across the pre-specified ligands measured in that unit.
The test is a Wilcoxon signed-rank test of that patient-level score against 0.
Cohort means are also combined with DerSimonian–Laird random effects. I² is that meta-analysis, not a cell-pooled test.

| family | expect | n patients | mean Δ | Wilcoxon p | DL mean | DL p | I² | agrees with expect |
|---|---|---:|---:|---:|---:|---:|---:|---|
| core_barrier | high>low | 59 | +0.259 | 2.65e-11 | +0.257 | 4.06e-31 | 44% | yes |
| ext_barrier | high>low | 59 | +0.131 | 2.39e-11 | +0.131 | 6.03e-09 | 83% | yes |
| core_ifn_recruit | low>high | 59 | -0.040 | 9.07e-06 | -0.043 | 0.0725 | 95% | yes |
| ext_ifn_recruit | low>high | 59 | +0.054 | 1.79e-09 | +0.056 | 7.34e-19 | 0% | no |
| mhci | low>high | 59 | +0.289 | 1.75e-05 | +0.261 | 0.0213 | 75% | no |

Core barrier ligands are F11R, NECTIN2, CDH1, LGALS9.
Core IFN/recruit ligands are CXCL9, CXCL10, CCL5, IFNG.
MHC-I (HLA-A/B/C) is its own row because the CellChat concordant-4 run found HLA–CD8 communication higher, not lower, from CLDN4-high senders. That opposite direction is what this quartile split shows as well, except GSE189357, where the MHC-I family mean is negative.

Core barrier stays positive in every leave-one-cohort-out (see `family_sender.tsv`). The core IFN/recruit family mean is negative in GSE123902, GSE205335, and GSE189357 and positive in GSE131907, so I² is high and the random-effects p does not match the pooled Wilcoxon.
CXCL9, CXCL10, and IFNG are detected in ≥10% of malignant cells in fewer than 10% of eligible units, so they are not potential ligands. Their within-patient deltas are near zero. CCL5 is the only core IFN/recruit ligand that passes the expression filter.

### Called ligands

| ligand | family | n | mean Δ Q4−Q1 | p | potential ligand | AUPR-corrected rank on T/NK genes up with CLDN4 | AUPR-corrected rank on T/NK genes down with CLDN4 |
|---|---|---:|---:|---:|---|---:|---:|
| F11R | core_barrier | 59 | +0.350 | 2.52e-11 | yes | 63.0 | 134.0 |
| NECTIN2 | core_barrier | 59 | +0.253 | 2.65e-11 | yes | 248.0 | 237.0 |
| CDH1 | core_barrier | 59 | +0.382 | 3.25e-11 | yes | 34.0 | 263.0 |
| LGALS9 | core_barrier | 59 | +0.052 | 0.00192 | yes | 218.0 | 56.0 |
| CXCL9 | core_ifn_recruit | 59 | -0.000 | 0.171 | no | not scored | not scored |
| CXCL10 | core_ifn_recruit | 59 | -0.005 | 0.98 | no | not scored | not scored |
| CCL5 | core_ifn_recruit | 59 | -0.121 | 1.16e-05 | yes | 241.0 | 3.0 |
| IFNG | core_ifn_recruit | 59 | -0.032 | 4.52e-06 | no | not scored | not scored |
| HLA-A | mhci | 59 | +0.329 | 4.14e-06 | yes | 87.0 | 69.0 |
| HLA-B | mhci | 59 | +0.280 | 5.93e-05 | yes | 128.0 | 229.0 |
| HLA-C | mhci | 59 | +0.258 | 1.01e-04 | no | not scored | not scored |
| CXCL11 | ext_ifn_recruit | 59 | -0.000 | 0.745 | no | not scored | not scored |
| PVR | ext_barrier | 59 | +0.047 | 1.51e-08 | yes | 50.0 | 223.0 |
| CD274 | ext_barrier | 59 | +0.010 | 0.0282 | yes | 187.0 | 20.0 |
| TGFB1 | ext_barrier | 59 | -0.067 | 0.149 | yes | 7.0 | 115.0 |
| HLA-E | ext_barrier | 59 | +0.186 | 1.11e-04 | yes | 198.0 | 77.0 |

Full per-ligand sender table: `results/tables/ligand_sender_delta.tsv`.

## Receiver geneset (patient-aware)

The empirical geneset is fit on author-annotated T/NK only (GSE131907 and GSE205335). Marker-gated T/NK calls in GSE123902 and GSE189357 are not used to choose the geneset.
Within each author cohort, Spearman correlation of the T/NK pseudobulk (`mean log1p(CP10k)`) with malignant CLDN4 % positive. The two correlations are Fisher-z combined (DerSimonian–Laird). A gene must keep the same sign in both cohorts.
Genes up with CLDN4: **12** genes. Rule: FALLBACK unadjusted p<=0.01, same sign in every contributing cohort, epithelial genes removed (FDR set had <15 genes).
Genes down with CLDN4: **62** genes. Rule: FALLBACK unadjusted p<=0.01, same sign in every contributing cohort, epithelial genes removed (FDR set had <15 genes).
Lung epithelial, secretory, and ciliated markers are removed before the cutoff (EPCAM, TACSTD2, CDH1, KRTs, CLDNs, surfactant and secretoglobin genes, ELF3, CEACAM5/6, WFDC2, FOXJ1, TMC5, and the rest of the list in METHODS).
The all-four correlation table is `tnk_gene_vs_cldn4_all4.tsv`. It is not the activity geneset.

Highest-confidence T/NK genes **up** with malignant CLDN4 %pos (meta ρ, p):

ATP1B1 (ρ=+0.70, p=3.47e-07), APP (ρ=+0.52, p=7.17e-04), MDK (ρ=+0.50, p=0.00151), TNFRSF12A (ρ=+0.48, p=0.00202), LAPTM4B (ρ=+0.47, p=0.00282), ICA1 (ρ=+0.46, p=0.00399), SNRNP25 (ρ=+0.45, p=0.00443), KNOP1 (ρ=+0.43, p=0.00704), CD59 (ρ=+0.43, p=0.00774), TSPAN13 (ρ=+0.42, p=0.00929), CYSTM1 (ρ=+0.47, p=0.00954), PHLDA2 (ρ=+0.42, p=0.00979)

Highest-confidence T/NK genes **down** with malignant CLDN4 %pos:

AKAP7 (ρ=-0.63, p=1.25e-05), UBE2V1 (ρ=-0.62, p=2.77e-05), SH2D1A (ρ=-0.60, p=6.25e-05), CPNE1 (ρ=-0.57, p=1.83e-04), HCST (ρ=-0.56, p=2.26e-04), TMA7 (ρ=-0.52, p=7.79e-04), S100A8 (ρ=-0.64, p=9.39e-04), SAMD3 (ρ=-0.51, p=9.66e-04), GPR174 (ρ=-0.51, p=0.00105), P2RX4 (ρ=-0.51, p=0.00107), NOP10 (ρ=-0.51, p=0.00118), NDUFA3 (ρ=-0.50, p=0.00137)

Background expressed in T/NK (≥10% of eligible units with detection in ≥10% of T/NK cells, gene measured in that unit's cohort): see `background_genes.tsv`.
Potential ligands (sender detected in ≥10% of malignant cells in ≥10% of eligible units, and ≥1 NicheNet-v2 receptor detected in T/NK at the same threshold): **286**.

## NicheNet activity

Activity asks which malignant-expressed ligands' NicheNet-v2 target profiles match the T/NK geneset.
It is not a second copy of the sender Δ. A ligand can be higher in CLDN4-high cells and still be a weak predictor of the T/NK program, or the reverse.

Core IFN/recruit contributes only CCL5 as a potential ligand, so a barrier-versus-IFN activity contrast inside that pair of families is not estimated.
The comparison that is estimated is whether the four core barrier ligands, as a set, have higher AUPR-corrected than a random draw of four potential ligands (one-sided, 4000 draws).

| geneset | n genes | core barrier mean AUPR-corrected | one-sided p vs random ligands |
|---|---:|---:|---:|
| a_priori_cytotoxicity | 14 | 0.005 | 0.817 |
| a_priori_exhaustion | 13 | -0.000 | 0.459 |
| a_priori_ifn | 12 | 0.015 | 0.905 |
| empirical_high | 12 | 0.003 | 0.467 |
| empirical_low | 62 | -0.001 | 0.728 |

Ranks below are among the potential ligands (lower rank = higher AUPR-corrected). Absolute AUPR-corrected values on these genesets are small; the rank is the comparison.

Prioritization on the empirical-high geneset, in the differential-NicheNet style, is the sum of z-scored AUPR-corrected, z-scored within-patient sender Δ (high − low), and z-scored best-receptor detection. Stored in `ligand_priority.tsv`.
AUPR-corrected remains the activity rank. The sum is a joint score, not a substitute for either arm.

Top potential ligands by AUPR-corrected on T/NK genes **up** with CLDN4:

| rank | ligand | AUPR-corrected | AUROC | Pearson |
|---:|---|---:|---:|---:|
| 1 | CD59 | 0.024 | 0.674 | +0.069 |
| 2 | CD58 | 0.019 | 0.667 | +0.068 |
| 3 | S100A9 | 0.012 | 0.622 | +0.079 |
| 4 | BMP2 | 0.011 | 0.675 | +0.073 |
| 5 | IL18 | 0.010 | 0.703 | +0.038 |
| 6 | ADAM17 | 0.009 | 0.690 | +0.087 |
| 7 | TGFB1 | 0.008 | 0.784 | +0.067 |
| 8 | EXOC3-AS1 | 0.006 | 0.684 | +0.044 |

Top potential ligands by AUPR-corrected on T/NK genes **down** with CLDN4:

| rank | ligand | AUPR-corrected | AUROC | Pearson |
|---:|---|---:|---:|---:|
| 1 | C5 | 0.017 | 0.458 | +0.049 |
| 2 | APP | 0.005 | 0.453 | +0.016 |
| 3 | CCL5 | 0.004 | 0.484 | +0.048 |
| 4 | LTB | 0.003 | 0.538 | +0.032 |
| 5 | IL15 | 0.002 | 0.508 | +0.020 |
| 6 | TNF | 0.002 | 0.505 | +0.021 |
| 7 | HBEGF | 0.002 | 0.441 | +0.019 |
| 8 | TNFSF12 | 0.002 | 0.497 | +0.015 |

A priori receiver programs are scored with the same function. They are not the CLDN4 contrast. Ranks among 286 potential ligands:

| ligand | IFN program | exhaustion program | cytotoxicity program |
|---|---:|---:|---:|
| F11R | 167 | 104 | 127 |
| NECTIN2 | 185 | 69 | 265 |
| CDH1 | 192 | 250 | 236 |
| LGALS9 | 256 | 55 | 85 |
| CCL5 | 38 | 42 | 50 |
| CD274 | 16 | 1 | 2 |
| HLA-A | 22 | 146 | 8 |
| HLA-B | 20 | 230 | 68 |

CD274 and the classical HLA ligands sit nearer the IFN and cytotoxicity programs than F11R, NECTIN2, or CDH1 do. That is prior structure plus receiver expression, not a claim that PD-L1 rose with CLDN4.

## Other patient-aware splits

Same eligible units. Median split: malignant CLDN4 above vs at-or-below the unit median. Percent-positive split: CLDN4 UMI > 0 vs 0. Family Wilcoxon:

| family | median-split mean Δ | p | %pos-split mean Δ | p |
|---|---:|---:|---:|---:|
| core_barrier | +0.199 | 2.94e-11 | +0.280 | 3.53e-12 |
| core_ifn_recruit | -0.029 | 1.33e-05 | -0.033 | 3.80e-05 |
| mhci | +0.216 | 3.89e-05 | +0.343 | 6.07e-07 |

Leave-one-cohort-out means for the primary Q4−Q1 family score are columns `loco_mean_drop_*` in `family_sender.tsv`.
A cohort that flips the family sign when it is the only one left out is visible there. Cohorts are not dropped to manufacture agreement, and no discordant accession is added.

## What this does not claim

- It does not re-estimate the locked concordant-4 CLDN4 %pos vs T/NK fraction correlation.
- It does not turn a NicheNet target score into a statement that a barrier ligand excludes T cells in space. That evidence is the CosMx contact result, not this table.
- It does not claim CXCL9/10 are abundant. If they fail the potential-ligand filter, their activity is not scored.
- Cell-pooled p-values are not reported.
- The ligand–target matrix is a published prior, not a model fit on these four cohorts.

## Reproduce

```bash
bash methods/nichenet_concordant4_cldn4/scripts/download.sh /tmp/concordant4_raw /tmp/nichenet_prior /tmp/nichenet_src
Rscript methods/nichenet_concordant4_cldn4/scripts/extract_gse205335.R --raw=/tmp/concordant4_raw --here=methods/nichenet_concordant4_cldn4 --out=/tmp/nichenet_work/extract
python3 methods/nichenet_concordant4_cldn4/scripts/extract_rest.py
python3 methods/nichenet_concordant4_cldn4/scripts/analyze.py
```

