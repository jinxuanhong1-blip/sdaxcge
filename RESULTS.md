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
- `methods/concordant4_pb_voom_edger_cldn4/figures/forest_family_score.png`
- `methods/concordant4_pb_voom_edger_cldn4/figures/heatmap_family_median_logfc.png`
- `methods/concordant4_pb_voom_edger_cldn4/figures/volcano_stacked_voom.png`
- `methods/concordant4_pb_voom_edger_cldn4/figures/n_honest.png`
