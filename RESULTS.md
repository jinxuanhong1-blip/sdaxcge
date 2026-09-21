# Malignant CLDN4 Q4 vs Q1, patient-level pseudo-bulk

Concordant-4 only: **GSE123902, GSE131907, GSE205335, GSE189357**. Malignant cells, CLDN4 percent-positive Q4 versus Q1, within each cohort. Positive log fold-change means higher in CLDN4-high.

The count matrices are already malignant UMI sums, one column per donor, sample, or patient. That is the muscat `pbDS` collapse. Gene tests use the two engines `pbDS` calls: **limma-voom** with robust empirical Bayes, and **edgeR** quasi-likelihood (`glmQLFit`, robust). Family scores are the mean edgeR log2 CPM of a shared gene panel, then a cohort linear mixed model (random intercept) and a DerSimonian–Laird random-effects meta with a Knapp–Hartung test. Quartile labels are the locked percent-positive split. They were not re-cut on the DE subset.

Script: `methods/concordant4_pb_voom_edger_cldn4/run.sh`.

## Honest n

The quartile vector has **65** units (13 donors + 21 samples + 22 patients + 9 patients). The DE uses only the Q1 and Q4 units that are in the UMI-sum: **18 vs 16**.

| Cohort | Unit | Vector n | Q1 / Q4 in the vector | Q1 / Q4 in DE | DE ids |
|---|---|---:|---|---|---|
| GSE123902 | donor | 13 | 4 / 3 | 4 / 3 | Q1 LX675, LX682, LX699, LX701; Q4 LX653, LX680, LX684 |
| GSE131907 | sample | 21 | 6 / 5 | 6 / 5 | Q1 EBUS_13, EBUS_15, EBUS_49, NS_02, NS_06, NS_16; Q4 EBUS_19, EBUS_28, NS_03, NS_04, NS_07 |
| GSE205335 | patient | 22 | 6 / 6 | 5 / 6 | Q1 P1015, P1062, P1063, P1090, P1119; Q4 P1016, P1025, P1037, P1084, P1089, P1115 |
| GSE189357 | patient | 9 | 3 / 2 | 3 / 2 | Q1 TD2, TD4, TD7; Q4 TD6, TD9 |

GSE131907 stays a **sample**. The locked table has no patient id, so these 21 tumor-bearing samples (n_malignant ≥ 20) are not collapsed to patients.

**P4001** is Q1 on the 22-patient GSE205335 vector (27 malignant cells) and is absent from the UMI-sum. It is not a DE unit.

**GSE189357 Q4 is 2 patients.** That cohort is inside the pooled fits and is flagged thin. A pre-specified sensitivity drops it.

Smallest malignant-cell count inside a DE unit: GSE123902 **46** (LX699), GSE131907 79, GSE205335 131, GSE189357 491. Dropping LX699 leaves the stacked IFN result in place (fry p = 4.2×10⁻⁴, 17 vs 16).

## Family scores

Shared panel: genes with a count of at least 10 in at least two DE samples of every cohort. IFN 203, MHC-I/APM 21, chemokine 14, tight junction 168 (CLDN4 held out), keratin 9. The chemokine score is the 14 genes present in all four cohorts, out of the 26-gene panel.

### Within cohort

Family-score logFC (Q4 − Q1) and ordinary t p-value.

| Family | GSE123902 (4 vs 3) | GSE131907 (6 vs 5) | GSE205335 (5 vs 6) | GSE189357 (3 vs 2) |
|---|---|---|---|---|
| IFN | −0.626 (0.16) | +0.146 (0.71) | **−1.620 (0.0029)** | −0.240 (0.40) |
| MHC-I/APM | −0.748 (0.13) | −0.071 (0.91) | **−1.770 (0.020)** | −0.339 (0.43) |
| chemokine | −1.164 (0.14) | −0.618 (0.24) | **−2.613 (0.0032)** | −0.624 (0.53) |
| TJ | +0.181 (0.57) | **+0.456 (0.0043)** | −0.173 (0.27) | +0.014 (0.88) |
| keratin | −0.545 (0.62) | +0.899 (0.62) | −1.402 (0.097) | +0.724 (0.042) |

GSE189357’s keratin p-value sits on two Q4 patients. Its within-cohort FDR across five families is 0.21.

### Pooled

| Family | Mixed model (n = 18 vs 16) | Cohort-adjusted OLS | RE meta, 4 cohorts (I²) |
|---|---|---|---|
| IFN | **−0.665, p = 0.0083, FDR = 0.021** | −0.643, p = 0.010 | −0.560, p = 0.23 (75%) |
| MHC-I/APM | **−0.782, p = 0.020, FDR = 0.033** | −0.802, p = 0.020 | −0.671, p = 0.13 (37%) |
| chemokine | **−1.398, p = 4.0×10⁻⁴, FDR = 0.0020** | −1.381, p = 4.5×10⁻⁴ | −1.256, p = 0.076 (54%) |
| TJ | +0.129, p = 0.19 | +0.131, p = 0.18 | +0.115, p = 0.47 (77%) |
| keratin | −0.152, p = 0.81 | −0.172, p = 0.80 | −0.106, p = 0.86 (65%) |

The mixed model and the cohort-adjusted regression use the 34 units. The random-effects meta uses four cohort estimates, so a Knapp–Hartung p-value with 3 degrees of freedom is a heterogeneity check, not a second copy of the unit-level test. IFN and tight junction are the heterogeneous scores (I² 75% and 77%). GSE205335 carries the IFN, MHC, and chemokine magnitude. GSE131907’s IFN score is flat.

The same signs were already in the earlier cohort-adjusted logCPM regression (IFN −0.58, MHC −0.78, chemokine −0.96, TJ +0.04). This pass replaces that OLS gene test with voom and edgeR.

## Gene-level voom and edgeR

Stacked limma-voom, design `~ cohort + Q4`, 12,193 genes after `filterByExpr`. edgeR on the same design agrees: correlation of logFC **0.84**, same sign for **86%** of genes.

| Family | Voom genes down / up | Median logFC | FDR < 0.05 | fry direction (p) | edgeR down / up |
|---|---|---:|---:|---|---|
| IFN | 176 / 30 | −0.813 | 37 | Down (3.3×10⁻⁴) | 167 / 39 |
| MHC-I/APM | 19 / 2 | −0.967 | 3 | Down (0.0044) | 18 / 3 |
| chemokine | 14 / 2 | −2.02 | 7 | Down (1.7×10⁻⁵) | 13 / 3 |
| TJ | 91 / 71 | −0.081 | 6 | Down (0.36) | 83 / 79 |
| keratin | 4 / 6 | +0.111 | 0 | Down (0.88) | 7 / 3 |

Competitive camera on the stacked voom fit gives the same directions: IFN p = 1.5×10⁻²⁵, chemokine p = 1.2×10⁻¹¹, MHC-I/APM p = 2.9×10⁻¹⁰, TJ p = 0.20, keratin p = 0.69.

Random-effects meta of the cohort voom logFCs (genes in at least 3 cohorts): IFN 169 down / 34 up, median −0.67; MHC 18 / 3 down-majority, median −0.81; chemokine 9 / 3, median −1.63; TJ 87 down / 76 up, median −0.03; keratin 5 up / 3 down, median +0.30. edgeR meta signs match (IFN 160 down / 43 up).

CLDN4 itself, held out of the tight-junction set, is higher in Q4: stacked voom logFC **+1.68** (p = 0.041); four-cohort voom meta logFC **+1.06** (p_t = 0.043).

Genes with FDR < 0.05 in the stacked voom fit include CCL5 (−3.25), TAP2 (−1.91), TAP1 (−1.73), STAT1 (−1.10), CXCL10 (−4.22), PSMB9 (−1.69), and **CLDN3 (+1.80)**. CDH1 (+0.95), OCLN (+0.59), and F11R (+0.47) point up with p between 0.07 and 0.09. TJP1 and CLDN7 are near zero. KRT8, KRT18, and KRT19 are near zero (logFC +0.16, +0.40, +0.06).

## Where the IFN signal sits

Leave-one-cohort-out of the stacked voom fit:

| Dropped | DE n (Q1 vs Q4) | IFN median logFC | IFN genes down / up | IFN fry p |
|---|---|---:|---|---:|
| GSE123902 | 14 vs 13 | −0.72 | 170 / 32 | 0.0030 |
| GSE131907 | 12 vs 11 | −1.12 | 186 / 20 | 4.3×10⁻⁵ |
| GSE205335 | 13 vs 10 | −0.36 | 140 / 65 | 0.097 |
| GSE189357 | 15 vs 14 | −0.87 | 177 / 27 | 4.2×10⁻⁴ |

The IFN, MHC, and chemokine down-shift is reproducible across engines and is still there after dropping the thin cohort. **GSE205335 supplies the magnitude.** Without it, most IFN genes still have a negative logFC (140 vs 65) and the fry p-value moves to 0.097. Within GSE131907 the IFN score is +0.15 and the gene median is −0.06 (108 genes down, 90 up).

## Tight junction and keratin

The locked tight-junction list is KEGG tight junction plus GO tight-junction organization, with CLDN4 removed. On that broad set the pooled score is a small positive number (+0.13, mixed-model p = 0.19) and the stacked fry test points **Down** with p = 0.36.

The cohort split is the result:

- **GSE131907** tight-junction score **+0.46** (p = 0.0043, within-cohort FDR 0.021). Gene median logFC +0.04, 81 genes up and 76 down, fry direction Up (p = 0.27).
- **GSE205335** score **−0.17** (p = 0.27). Gene counts 65 up and 109 down.

GSE205335 Q4 versus Q1 is histology-imbalanced. In the DE tails, Q1 is 4 adenocarcinoma and 1 squamous; Q4 is 2 adenocarcinoma, **3 small-cell**, and 1 squamous. Small-cell enrichment in the CLDN4-high arm pulls epithelial programs down and works against a tight-junction or keratin increase. Keratin follows that split: GSE205335 score −1.40 (p = 0.097), GSE189357 +0.72 on two Q4 patients, pooled score −0.15 (p = 0.81).

CLDN3 is the tight-junction gene that survives the stacked model (logFC +1.80, FDR = 0.014). The broad set does not.

## Robustness

Pre-specified in `analyze_robust.R` and run after the numbers above were locked. The Q4-versus-Q1 interferon score in this script matches the cohort-adjusted regression already reported (−0.643, n = 34, Q1 = 18, Q4 = 16). Family-score FDR values in this section are Benjamini–Hochberg across IFN, MHC-I/APM, and chemokine inside that model. They are a smaller family set than the five-family FDR in the pooled table above.

Continuous coefficients are the change in the family score per within-cohort standard deviation of malignant CLDN4 percent-positive. Q4-versus-Q1 coefficients are the group difference. cGAS is MB21D1 in GSE123902 and GSE131907 and CGAS in GSE205335 and GSE189357. STING_core is those eight genes with IRF7 and ZBP1 left out, because those two symbols are already in the Hallmark interferon set. NHEJ is XRCC4, XRCC5, XRCC6, LIG4, PRKDC, NHEJ1, DCLRE1C, APLF, POLL, and POLM. PAXX is missing from two cohorts and is not in the score. NHEJ had no pre-specified direction.

The keratin covariate is the mean log2 CPM of KRT8, KRT18, and KRT19, z-scored inside each cohort. The immune-leak covariate is the mean log2 CPM of PTPRC, CD3D, CD3E, CD68, NKG7, and MS4A1. Malignant fraction is n_malignant / n_cells on the biopsy, which is composition of the biopsy and not purity of an already-malignant UMI sum. GSE123902, GSE131907, and GSE189357 are lung-adenocarcinoma spectrum, so a LUAD cut that also drops GSE205335 is the same 13-versus-10 set as leaving GSE205335 out. Those two labels are one contrast.

### Locked Q1 = 18 versus Q4 = 16, with STING and NHEJ

| Score | Mixed model | Cohort-adjusted OLS |
|---|---|---|
| STING_core (8 genes) | **−0.803, p = 0.0094** | −0.797, p = 0.013 |
| STING plus IRF7 and ZBP1 | −0.762, p = 0.018 | −0.752, p = 0.025 |
| NHEJ (10 genes) | +0.040, p = 0.66 | +0.037, p = 0.69 |

STING_core moves with the interferon score on the locked contrast. NHEJ sits on zero. Within cohort, STING_core is −2.08 in GSE205335 (p = 0.0025) and −1.25 in GSE123902 (p = 0.046). GSE131907 is +0.51 (p = 0.21).

Adjusting the locked contrast for KRT8/KRT18/KRT19 leaves the three claim scores down and slightly steeper: IFN −0.697 (p = 0.0015, FDR = 0.0022), MHC-I/APM −0.828 (p = 0.0071), chemokine −1.452 (p = 2.9×10⁻⁵). The interferon shift is not a low-keratin artifact of the small-cell cases in Q4.

### Leaving GSE205335 out

Q1 = 13, Q4 = 10, three cohorts. Mixed model:

| Score | logFC | p | Three-family FDR |
|---|---:|---:|---:|
| IFN | −0.177 | 0.45 | 0.45 |
| MHC-I/APM | −0.340 | 0.32 | 0.45 |
| chemokine | **−0.791** | **0.034** | 0.10 |
| STING_core | −0.167 | 0.56 |  |
| NHEJ | +0.063 | 0.46 |  |
| TJ | **+0.277** | **0.022** |  |

Chemokine is the claim score that stays negative with a unit-level p-value below 0.05. Its FDR across the three claim scores is 0.10. The three cohort slopes are all negative (GSE123902 −1.16, GSE131907 −0.62, GSE189357 −0.62) and none is significant on its own. A voom refit on these 23 units gives chemokine median logFC −1.74, 10 genes down and 4 up, fry p = 0.0070. IFN on that refit is median −0.31, 136 down and 67 up, fry p = 0.14. MHC-I/APM is 17 down and 4 up, fry p = 0.17. The earlier stacked-matrix leave-one-out (IFN fry p = 0.097, 140 versus 65) is a different fit and is unchanged.

The three-cohort Knapp–Hartung meta of the chemokine slopes is −0.78 (p = 0.047, I² = 0). With I² at zero the Knapp–Hartung scale factor drops below 1 and the interval gets narrower than the usual random-effects interval. Flooring that factor at 1 (`test = "adhoc"`) gives p = 0.16. The unit-level mixed-model p-value, 0.034, is the number for this contrast.

The same Q4-versus-Q1 model with the keratin covariate, still without GSE205335: chemokine −0.925 (p = 0.010, FDR = 0.031), MHC-I/APM −0.462 (p = 0.17), IFN −0.254 (p = 0.26).

Continuous percent-positive on the 43 units outside GSE205335 is flat for IFN (+0.005, p = 0.96). MHC-I/APM is −0.083 (p = 0.52) and chemokine is −0.168 (p = 0.33). Voom fry for IFN is p = 0.50 (118 down, 85 up). Chemokine fry is p = 0.065 (9 down, 5 up). edgeR signs on the same continuous design are IFN 109 down / 94 up and chemokine 8 / 6. Spearman meta-analysis of the ranks is IFN rho −0.09 (p = 0.65). GSE123902 alone has a continuous chemokine slope of −0.61 (p = 0.050, 13 donors). GSE131907 and GSE189357 are near zero.

Tight junction is the score that rises once GSE205335 is out: +0.277 on the Q4-versus-Q1 mixed model (p = 0.022) and +0.135 per SD on the continuous model (p = 8.0×10⁻⁴). GSE131907 carries the continuous tight-junction slope (+0.205, p = 2.5×10⁻⁵).

### Histology

Dropping the three small-cell Q4 patients (P1016, P1025, P1115) leaves Q1 = 18 and Q4 = 13. Squamous cases stay, and GSE205335 is still in the fit. Mixed model: chemokine −1.059 (p = 0.0035, FDR = 0.011), IFN −0.423 (p = 0.058, FDR = 0.072), MHC-I/APM −0.596 (p = 0.072). Voom fry: IFN p = 0.0098 (163 down, 40 up), MHC-I/APM p = 0.039, chemokine p = 4.3×10⁻⁴, STING_core p = 0.047. edgeR signs agree (IFN 158 / 45, chemokine 11 / 3).

Restricting GSE205335 to adenocarcinoma and keeping the three adenocarcinoma-spectrum cohorts gives Q1 = 17 and Q4 = 12. Chemokine −0.898 (p = 0.013, FDR = 0.038). IFN −0.350 (p = 0.13) with fry p = 0.033 (160 down, 43 up). MHC-I/APM −0.560 (p = 0.11). Tight junction +0.233 (p = 0.022).

GSE205335 adenocarcinoma only, continuous, n = 13: IFN −0.297 (p = 0.21), MHC-I/APM −0.497 (p = 0.13), chemokine −0.346 (p = 0.32). The sign matches the full GSE205335 slope. Thirteen patients do not separate that slope from the small-cell cases.

### Continuous percent-positive on all 64 units

This includes GSE205335. Mixed model, per within-cohort SD: IFN −0.175 (p = 0.067, FDR = 0.067), MHC-I/APM −0.273 (p = 0.024, FDR = 0.036), chemokine −0.400 (p = 0.014, FDR = 0.036), STING_core −0.238 (p = 0.036), NHEJ +0.002 (p = 0.97). Voom fry: IFN p = 0.0074 (166 down, 37 up), MHC-I/APM p = 0.010, chemokine p = 9.8×10⁻⁴, STING_core p = 0.049.

Keratin adjustment on that full continuous model (cohort-adjusted OLS): IFN −0.269 (p = 0.0013), MHC-I/APM −0.381 (p = 6.2×10⁻⁴), chemokine −0.542 (p = 2.9×10⁻⁴). GSE205335 is inside those three p-values. The same keratin-adjusted slopes fit inside each cohort and then meta-analyzed without GSE205335 are all negative (IFN −0.17, MHC-I/APM −0.25, chemokine −0.42). The Knapp–Hartung p-value for MHC-I/APM is 0.0059 because the scale factor collapses (I² = 0). The floored test is p = 0.17. The unit-level keratin-adjusted mixed model on the same 43 units is MHC-I/APM −0.220 (p = 0.096) and chemokine −0.305 (p = 0.096). Voom fry for that keratin-adjusted continuous model without GSE205335: chemokine p = 0.010 (10 down, 4 up), MHC-I/APM p = 0.049 (17 down, 4 up), IFN p = 0.084 (138 down, 65 up).

### Immune leak and malignant fraction

Adding the six immune-lineage genes to the locked Q4-versus-Q1 model moves the claim scores to IFN −0.261 (p = 0.46), MHC-I/APM −0.497 (p = 0.33), and chemokine −0.157 (p = 0.72). Keratin plus immune leak on the same 34 units: chemokine −0.359 (p = 0.32), IFN −0.426 (p = 0.16). Without GSE205335, keratin plus immune leak gives chemokine +0.006 (p = 0.99). The chemokine association is sensitive to residual immune transcripts in the malignant UMI sum.

Malignant fraction on the locked contrast leaves chemokine at −0.978 (p = 0.055) and IFN at −0.489 (p = 0.15). Without GSE205335, IFN is −0.248 (p = 0.50).

### Reading the grid

Chemokine stays lower in CLDN4-high malignant cells after GSE205335 is removed, on the Q4-versus-Q1 mixed model and on the keratin-adjusted version of that model, and the voom gene set agrees. IFN and MHC-I/APM stay lower on the full four-cohort continuous model and on the histology cuts that still contain GSE205335 adenocarcinoma. Their family scores move to about zero on the continuous model once GSE205335 is removed, with gene signs still mostly negative and fry p-values of 0.14 and 0.17 on the Q4-versus-Q1 refit. STING_core follows interferon on the locked 18-versus-16 contrast and in GSE123902, and it does not survive leaving GSE205335 out. NHEJ does not move. Tight junction is higher in CLDN4-high cells in the fits that are free of the GSE205335 small-cell Q4 arm. Immune-leak adjustment removes the chemokine difference.

## What this does not say

- The quartile-vector n of 65 is not the DE n.
- GSE131907 contributes 11 samples, not 11 patients.
- This is not a T/NK re-analysis and does not reopen the concordant-4 Spearman.
- GSE148071, GSE127465, GSE207422, GSE154826, and GSE200563 stay out.
- No TACSTD2×CLDN4 dual-high score.
- The association is not a claim that CLDN4 causes the interferon change.
- Genome-wide FDR outside these families is thin; CLDN4’s own stacked FDR is 0.29.

## Files

- `methods/concordant4_pb_voom_edger_cldn4/tables/n_honest.tsv`
- `methods/concordant4_pb_voom_edger_cldn4/tables/family_effects.tsv`
- `methods/concordant4_pb_voom_edger_cldn4/tables/family_direction.tsv`
- `methods/concordant4_pb_voom_edger_cldn4/tables/geneset_fry_camera.tsv`
- `methods/concordant4_pb_voom_edger_cldn4/tables/sensitivity_loo.tsv`
- `methods/concordant4_pb_voom_edger_cldn4/tables/gse205335_histology_qtails.tsv`
- `methods/concordant4_pb_voom_edger_cldn4/tables/robust_models.tsv`
- `methods/concordant4_pb_voom_edger_cldn4/tables/robust_fry.tsv`
- `methods/concordant4_pb_voom_edger_cldn4/tables/robust_edger_signs.tsv`
- `methods/concordant4_pb_voom_edger_cldn4/tables/robust_n.tsv`
- `methods/concordant4_pb_voom_edger_cldn4/figures/forest_family_score.png`
- `methods/concordant4_pb_voom_edger_cldn4/figures/forest_ifn_robust.png`
- `methods/concordant4_pb_voom_edger_cldn4/figures/heatmap_robust_grid.png`
- `methods/concordant4_pb_voom_edger_cldn4/figures/heatmap_family_median_logfc.png`
- `methods/concordant4_pb_voom_edger_cldn4/figures/volcano_stacked_voom.png`
- `methods/concordant4_pb_voom_edger_cldn4/figures/n_honest.png`
