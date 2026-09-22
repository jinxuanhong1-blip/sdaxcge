# FINDING — GSE131907 max effect, TJ Δ and TACSTD2 vs T

Kim et al., *Nat Commun* 2020, GSE131907. Public processed UMI matrix only (208,506 cells). The log2TPM text and EGA FASTQ were not used.

Inferential unit is the **patient**. Each TJ Δ is a within-patient contrast (TACSTD2-high minus TACSTD2-low cells in the gate). Each T association is one TACSTD2 summary and one T fraction per patient. Sample-level tests that repeat patients are not eligible to win. p-values from the sweep are descriptive.

The pipeline first reproduces the locked author-malignant CLDN4 % positive versus T/NK row (n=21): sample ρ=-0.522 (locked -0.522), patient ρ=-0.478 (locked -0.478).

## TJ Δ

Score scale is locked: mean of log1p(UMI / full library × 10,000). TACSTD2 is not in any TJ score. A panel needs ≥3 genes, n≥8 patients, and both the mean and the median paired Δ > 0. Epithelial controls cannot win.

**Maximum mean paired Δ = 0.542** on `claudins` (4 genes). Gate `epi` (all tumor-site epithelial cells), site `met_lung` (tL/B, mLN, and pleural effusion), pool `patient`, split `q4q1`, score `cell_mean`, min cells/arm 10. n=10 patients. Median Δ=0.556. Fraction of patients with Δ>0 = 1.000. Two-sided Wilcoxon p=0.001953.

Best author-malignant panel (metastatic cells only, not the broader epithelial gate): `claudins`, site `no_mBrain` (tumor origins except brain metastasis), split `pos`, mean Δ=0.513, n=10, Wilcoxon p=0.001953. Author-malignant cells are absent from primary lung and from pleural effusion, so the author `no_mBrain` and author `met_lung` rows are the same cells.

At the same contrast, the epithelial control EPCAM/KRT8/18/19 mean Δ is 0.514 (n=10) and KRT8/18/19 alone is 0.489. The winning TJ module is larger than that epithelial control on this contrast. After a cell-level linear residual of the winning score on EPCAM/KRT8/18/19, the mean paired Δ is 0.398 (median 0.379, Wilcoxon p=0.001953, n=10).

Pre-specified core (CLDN1/3/4/7, OCLN, TJP1/2/3, F11R, CGN, MARVELD2, CRB3), median split, cell-mean score, patient pool: author-malignant mean Δ=0.131 (n=20, p=9.537e-06); primary tS mean Δ=0.141 (n=10, p=0.001953).

Sweep size: 6624 module rows, 4080 eligible. Largest single-gene mean Δ is CLDN4 (epi, no_mBrain, q4q1): mean Δ=0.851, median Δ=0.835, n=21, Wilcoxon p=9.537e-07. That row is not a multi-gene TJ module.

### Pre-specified TJ rows

| gate | site | pool | split | score | module | genes | n | mean Δ | median Δ | frac>0 | Wilcoxon p |
|---|---|---|---|---|---|---|---:|---:|---:|---:|---:|
| author | all | patient | median | cell_mean | pre_core | 12 | 20 | 0.131 | 0.124 | 0.900 | 9.537e-06 |
| author | all | patient | q4q1 | cell_mean | pre_core | 12 | 20 | 0.209 | 0.223 | 0.900 | 9.537e-06 |
| author | all | patient | median | cell_mean | pre_core_no_cldn4 | 11 | 20 | 0.103 | 0.104 | 0.900 | 9.537e-06 |
| ts | tLung | patient | median | cell_mean | pre_core | 12 | 10 | 0.141 | 0.137 | 1.000 | 0.001953 |
| ts | tLung | patient | q4q1 | cell_mean | pre_core | 12 | 10 | 0.240 | 0.243 | 1.000 | 0.001953 |
| broad | all | patient | median | cell_mean | pre_core | 12 | 30 | 0.135 | 0.126 | 0.933 | 9.313e-09 |

### Eight largest eligible mean Δ

| gate | site | pool | split | score | module | genes | n | mean Δ | median Δ | frac>0 | Wilcoxon p |
|---|---|---|---|---|---|---|---:|---:|---:|---:|---:|
| epi | met_lung | patient | q4q1 | cell_mean | claudins | 4 | 10 | 0.542 | 0.556 | 1.000 | 0.001953 |
| epi | met_lung | sample_then_patient | pos | cell_mean | claudins | 4 | 10 | 0.538 | 0.493 | 1.000 | 0.001953 |
| epi | met_lung | sample_then_patient | q4q1 | cell_mean | claudins | 4 | 10 | 0.538 | 0.556 | 1.000 | 0.001953 |
| epi | met_lung | patient | pos | cell_mean | claudins | 4 | 10 | 0.537 | 0.493 | 1.000 | 0.001953 |
| author | no_mBrain | patient | pos | cell_mean | claudins | 4 | 10 | 0.513 | 0.483 | 1.000 | 0.001953 |
| author | no_mBrain | sample_then_patient | pos | cell_mean | claudins | 4 | 10 | 0.513 | 0.483 | 1.000 | 0.001953 |
| author | met_lung | patient | pos | cell_mean | claudins | 4 | 10 | 0.513 | 0.483 | 1.000 | 0.001953 |
| author | met_lung | sample_then_patient | pos | cell_mean | claudins | 4 | 10 | 0.513 | 0.483 | 1.000 | 0.001953 |

## TACSTD2-high versus T depletion

Pre-specified primary: author-malignant cells, samples that contain those cells, TACSTD2 % (UMI>0) versus T-cell fraction, ≥20 malignant cells. ρ=-0.157, Spearman p=0.4963, n=21. Primary tS, same predictor versus T fraction: ρ=-0.188, p=0.6032, n=10.

Smallest nominal depletion p is 0.03681 (spearman, author, mBrain, gate_samples, tac_pct2 vs cyto_frac, n=10 patients, ρ=-0.663). Spearman tests in the sweep: 1932. Depletion-direction and n≥8: 605, of which 12 have nominal p<0.05. Those p<0.05 rows collapse to: cyto_frac at mBrain, tac_pct2, ρ=-0.663, n=10 (12 duplicate rows). Total T fraction (T_frac) never reaches p<0.05 in this sweep (rows below 0.05: 0). Q4 vs Q1 tests: 420, of which 0 are eligible and p<0.05. Bonferroni 0.05 line for the Spearman family: 2.588e-05. Permutation p for this single panel is 0.0406 (4999 shuffles). Partial Spearman on the gate-cell fraction is ρ=-0.537, p=0.1097 (n=10). After the gate-fraction partial, that nominal p is no longer below 0.05. Across these patients the outcome totals about 29 cells (2 patients are zero). Leave-one-out Spearman p ranges from 0.009106 to 0.1373. The searched minimum is not a confirmatory p-value.

T outcomes that can win: T fraction of the sample, T fraction of immune cells, CD8 fraction, CD8 fraction of immune cells, cytotoxic CD8 fraction, exhausted CD8 fraction. T/NK is reported in the sweep and cannot win. Denominators are either the samples that contain the gate cells, or every tumor sample in the site (including a capture with no gate cells). The second denominator is the patient pool used for the locked CLDN4 row.

### Pre-specified T rows

| gate | site | denominator | predictor | outcome | n | ρ | Spearman p | Q4−Q1 | Q4 vs Q1 p |
|---|---|---|---|---|---:|---:|---:|---:|---:|
| author | all | gate_samples | tac_pct | T_frac | 21 | -0.157 | 0.4963 | -0.113 | 0.6623 |
| author | all | all_tumor_in_site | tac_pct | T_frac | 21 | -0.025 | 0.9154 | -0.039 | 0.7922 |
| author | all | gate_samples | tac_mean | T_frac | 21 | -0.022 | 0.9243 | -0.031 | 0.9307 |
| author | all | gate_samples | tac_pct | CD8_frac | 21 | -0.078 | 0.7371 | -0.035 | 0.5368 |
| ts | tLung | gate_samples | tac_pct | T_frac | 10 | -0.188 | 0.6032 | -0.025 | 0.7 |
| ts | tLung | gate_samples | tac_mean | T_frac | 10 | 0.103 | 0.777 | 0.210 | 1 |
| ts | tLung | gate_samples | tac_pct | CD8_frac | 10 | 0.588 | 0.07388 | 0.071 | 0.7 |
| broad | all | gate_samples | tac_pct | T_frac | 31 | 0.221 | 0.2331 | 0.016 | 0.5054 |
| epi | all | gate_samples | tac_pct | T_frac | 32 | 0.184 | 0.3134 | 0.021 | 0.3823 |

## What this is not

- Not a cell-level Wilcoxon. The MPE and primary-tS cell-level TJ p-values in the paper-funnel pass are a different unit.
- Not a claim that TACSTD2 causes T-cell loss, and not a spatial exclusion test. This atlas is dissociated.
- Not a confirmatory p-value for the searched minimum.
- Not GSE148071, not the four-cohort pool, and not a mouse Tacstd2 result.
- Author `Malignant cells` are metastatic (mBrain, tL/B, mLN), not primary tLung. Primary tumor epithelium in this matrix is tS1/tS2/tS3.

## Files

- `results/tables/tj_sweep.tsv` — module sweep
- `results/tables/tj_patient_deltas.tsv` — one Δ per patient
- `results/tables/t_sweep.tsv` — T/CD8 sweep
- `results/tables/calibration_cldn4_tnk.tsv` — lock check
- `results/figures/fig_tj_patient_deltas.png`
- `results/figures/fig_tacstd2_vs_t.png`
- `results/figures/fig_honest_n.png`
