# Concordant-4 malignant cells: patient-level GSVA / ssGSEA and fgsea by CLDN4

ADDITIVE. **CLDN4-only.** Cohorts are the locked concordant four:
**GSE123902 + GSE131907 + GSE205335 + GSE189357**. Not GSE148071, GSE127465,
GSE207422, GSE154826, or CD45+/T-only extracts. The T/NK result
(malignant CLDN4 %pos vs T/NK, DL rho about -0.53, N=65) is not re-derived.

The unit is the patient, donor, or sample malignant pseudobulk. Cell-level
p-values are not computed. Quartiles are the locked within-cohort CLDN4 %pos
labels. P4001 is out of the malignant UMI sum, so expression n is not 65.

## Prespecified sets

Scored separately. IFN-alpha and IFN-gamma are not merged. Antigen presentation
is the frozen 21-gene classical MHC-I / APM panel (MHC-II excluded; not a
Hallmark set). Junction sets are Hallmark apical junction and KEGG tight
junction. Hallmark EMT is included as its own set. CLDN4 is removed from every
set before scoring. Two GO tight-junction sets are secondary.

Shared genes for GSVA/ssGSEA (count >= 10 in >= 3 units inside every cohort): **12285**. Within-cohort gene counts: GSE123902 13549; GSE131907 17372; GSE205335 20393; GSE189357 15607.

## Honest n

| cohort | unit in the locked table | in the count matrix | Q1 | Q4 |
|---|---:|---:|---:|---:|
| GSE123902 | 13 | 13 | 4 | 3 |
| GSE131907 | 21 | 21 | 6 | 5 |
| GSE205335 | 22 | 21 | 5 | 6 |
| GSE189357 | 9 | 9 | 3 | 2 |
| pooled expression | 65 | **64** | **18** | **16** |

GSE189357 Q4 has 2 units, so its within-cohort binary test is skipped. Those
two units stay in the pooled Q4 vs Q1 model. Dropped from the matrix: P4001.

Direction check, same pooled OLS as the ranking (positive = higher in CLDN4-high):
CLDN4 Q4 vs Q1 logFC **1.673** (p 0.012, n=34, genes 21604). Continuous logFC per SD of %pos **1.116** (p 8.38e-10).

## 1. fgsea NES (headline)

Ranking for the pool is the patient-level OLS t on log2(TMM-CPM+1),
`~ cohort + exposure`, built the same way as the concordant-4 malignant DE
(zero-fill genes absent from a cohort, then one TMM). Positive NES = the set
sits toward genes that are higher in CLDN4-high malignant pseudobulks.

p-values are fgseaMultilevel. The interval is a patient bootstrap (B=200, fgseaSimple nperm=1000) on that fixed logCPM matrix; normalization
factors are not re-estimated inside the bootstrap. FDR is Benjamini-Hochberg
across the six primary sets, within the contrast.

### Pooled Q4 vs Q1 (n=34, 18 vs 16)

| set | NES | bootstrap 95% | p | FDR | size |
|---|---:|---|---:|---:|---:|
| Hallmark IFN-alpha | -3.42 | -3.74 to -1.99 | 8.16e-24 | 1.63e-23 | 96 |
| Hallmark IFN-gamma | -3.85 | -4.01 to -2.70 | 1.05e-49 | 6.33e-49 | 198 |
| MHC-I antigen presentation | -2.75 | -3.05 to -1.76 | 4.54e-08 | 6.82e-08 | 21 |
| Hallmark EMT | -3.42 | -3.66 to -2.19 | 1.38e-34 | 4.15e-34 | 199 |
| Hallmark apical junction | -1.62 | -2.28 to 0.95 | 3.94e-05 | 4.73e-05 | 189 |
| KEGG tight junction | 1.13 | -1.12 to 1.42 | 0.244 | 0.244 | 149 |

### Pooled continuous CLDN4 %pos (n=64, t per SD)

| set | NES | bootstrap 95% | p | FDR | size |
|---|---:|---|---:|---:|---:|
| Hallmark IFN-alpha | -2.61 | -3.39 to -0.80 | 8.74e-12 | 1.75e-11 | 96 |
| Hallmark IFN-gamma | -3.08 | -3.74 to -1.69 | 3.59e-28 | 2.15e-27 | 198 |
| MHC-I antigen presentation | -2.44 | -2.99 to -1.15 | 2.21e-06 | 3.32e-06 | 21 |
| Hallmark EMT | -3.09 | -3.41 to -1.77 | 1.70e-27 | 5.11e-27 | 199 |
| Hallmark apical junction | -1.59 | -2.10 to 1.04 | 2.45e-04 | 2.95e-04 | 193 |
| KEGG tight junction | 1.32 | -1.18 to 1.67 | 0.045 | 0.045 | 153 |

Pooled Q4 intervals that stay negative: IFN-alpha (−3.74 to −1.99), IFN-gamma
(−4.01 to −2.70), MHC-I antigen presentation (−3.05 to −1.76), and Hallmark
EMT (−3.66 to −2.19). Hallmark apical junction (−2.28 to 0.95) and KEGG tight
junction (−1.12 to 1.42) cross zero. The fgsea p-value conditions on the
observed ranking. The bootstrap resamples patients. Where they disagree, the
bootstrap is the uncertainty that matches the unit of analysis. A cohort point
can sit just outside its interval when the multilevel NES and the
simple-permutation bootstrap differ slightly.

Leave-one-cohort-out keeps IFN-alpha, IFN-gamma, MHC-I, and Hallmark EMT
negative, including the drop of GSE205335 (NES −1.99, −2.37, −1.99, −3.35).
GSE131907 alone is flat for IFN. The pool is not one cohort, and it is not
uniform.

Hallmark apical junction is a broad MSigDB set, not a pure tight-junction
list. KEGG tight junction and the two GO tight-junction sets are the junction
lists. CLDN4 is held out of every set. Putting it back does not flip the sign.
This run is not tight-junction up.

### Cohort NES (same sets; binary skipped when an arm has n<3)

| set | cohort | n | NES | bootstrap 95% | p | FDR |
|---|---|---:|---:|---|---:|---:|
| Hallmark IFN-alpha | GSE123902 | 7 | -2.50 | -2.62 to 0.83 | 7.06e-10 | 1.41e-09 |
| Hallmark IFN-alpha | GSE131907 | 11 | -0.98 | -2.70 to 1.68 | 0.499 | 0.598 |
| Hallmark IFN-alpha | GSE205335 | 11 | -3.54 | -3.51 to -2.18 | 4.25e-32 | 1.28e-31 |
| Hallmark IFN-alpha | pooled | 34 | -3.42 | -3.74 to -1.99 | 8.16e-24 | 1.63e-23 |
| Hallmark IFN-gamma | GSE123902 | 7 | -2.53 | -2.72 to -1.02 | 1.89e-14 | 5.66e-14 |
| Hallmark IFN-gamma | GSE131907 | 11 | -1.32 | -2.51 to 1.37 | 0.019 | 0.056 |
| Hallmark IFN-gamma | GSE205335 | 11 | -3.88 | -3.80 to -2.56 | 2.25e-60 | 1.35e-59 |
| Hallmark IFN-gamma | pooled | 34 | -3.85 | -4.01 to -2.70 | 1.05e-49 | 6.33e-49 |
| MHC-I antigen presentation | GSE123902 | 7 | -2.31 | -2.30 to -1.06 | 6.28e-05 | 9.42e-05 |
| MHC-I antigen presentation | GSE131907 | 11 | -1.10 | -2.31 to 1.62 | 0.303 | 0.455 |
| MHC-I antigen presentation | GSE205335 | 11 | -2.83 | -2.82 to -1.64 | 1.91e-12 | 3.82e-12 |
| MHC-I antigen presentation | pooled | 34 | -2.75 | -3.05 to -1.76 | 4.54e-08 | 6.82e-08 |
| Hallmark EMT | GSE123902 | 7 | -3.12 | -3.14 to -1.14 | 3.09e-25 | 1.86e-24 |
| Hallmark EMT | GSE131907 | 11 | -1.68 | -2.16 to 0.91 | 2.81e-05 | 1.69e-04 |
| Hallmark EMT | GSE205335 | 11 | -1.69 | -2.29 to 0.83 | 8.87e-05 | 1.33e-04 |
| Hallmark EMT | pooled | 34 | -3.42 | -3.66 to -2.19 | 1.38e-34 | 4.15e-34 |
| Hallmark apical junction | GSE123902 | 7 | -1.77 | -2.21 to 2.14 | 1.86e-04 | 2.23e-04 |
| Hallmark apical junction | GSE131907 | 11 | -0.81 | -1.36 to 1.09 | 0.952 | 0.952 |
| Hallmark apical junction | GSE205335 | 11 | -1.51 | -2.07 to 1.24 | 0.002 | 0.003 |
| Hallmark apical junction | pooled | 34 | -1.62 | -2.28 to 0.95 | 3.94e-05 | 4.73e-05 |
| KEGG tight junction | GSE123902 | 7 | 1.21 | -1.25 to 2.32 | 0.134 | 0.134 |
| KEGG tight junction | GSE131907 | 11 | 1.21 | -0.85 to 1.47 | 0.137 | 0.274 |
| KEGG tight junction | GSE205335 | 11 | -1.34 | -1.61 to 0.98 | 0.032 | 0.032 |
| KEGG tight junction | pooled | 34 | 1.13 | -1.12 to 1.42 | 0.244 | 0.244 |

### Continuous cohort NES

| set | cohort | n | NES | bootstrap 95% | p | FDR |
|---|---|---:|---:|---|---:|---:|
| Hallmark IFN-alpha | GSE123902 | 13 | -2.38 | -2.67 to 0.85 | 2.20e-09 | 4.40e-09 |
| Hallmark IFN-alpha | GSE131907 | 21 | 0.90 | -2.17 to 1.91 | 0.665 | 0.767 |
| Hallmark IFN-alpha | GSE205335 | 21 | -3.52 | -3.61 to -0.97 | 5.46e-29 | 1.64e-28 |
| Hallmark IFN-alpha | GSE189357 | 9 | -1.18 | -1.90 to 1.29 | 0.168 | 0.202 |
| Hallmark IFN-alpha | pooled | 64 | -2.61 | -3.39 to -0.80 | 8.74e-12 | 1.75e-11 |
| Hallmark IFN-gamma | GSE123902 | 13 | -2.76 | -2.86 to -1.15 | 5.16e-20 | 1.55e-19 |
| Hallmark IFN-gamma | GSE131907 | 21 | 0.87 | -2.05 to 1.81 | 0.767 | 0.767 |
| Hallmark IFN-gamma | GSE205335 | 21 | -3.92 | -3.85 to -1.81 | 7.61e-57 | 4.56e-56 |
| Hallmark IFN-gamma | GSE189357 | 9 | -1.37 | -2.46 to 1.08 | 0.015 | 0.030 |
| Hallmark IFN-gamma | pooled | 64 | -3.08 | -3.74 to -1.69 | 3.59e-28 | 2.15e-27 |
| MHC-I antigen presentation | GSE123902 | 13 | -2.30 | -2.37 to 0.70 | 2.10e-05 | 3.15e-05 |
| MHC-I antigen presentation | GSE131907 | 21 | -0.90 | -2.25 to 1.58 | 0.609 | 0.767 |
| MHC-I antigen presentation | GSE205335 | 21 | -3.01 | -3.01 to -0.81 | 4.42e-14 | 8.83e-14 |
| MHC-I antigen presentation | GSE189357 | 9 | -0.61 | -2.12 to 1.59 | 0.958 | 0.958 |
| MHC-I antigen presentation | pooled | 64 | -2.44 | -2.99 to -1.15 | 2.21e-06 | 3.32e-06 |
| Hallmark EMT | GSE123902 | 13 | -3.03 | -3.30 to 0.66 | 7.42e-26 | 4.45e-25 |
| Hallmark EMT | GSE131907 | 21 | -1.84 | -2.42 to 0.75 | 4.83e-06 | 2.90e-05 |
| Hallmark EMT | GSE205335 | 21 | -1.43 | -2.16 to 0.99 | 0.005 | 0.006 |
| Hallmark EMT | GSE189357 | 9 | -2.70 | -2.85 to -1.54 | 4.96e-20 | 2.97e-19 |
| Hallmark EMT | pooled | 64 | -3.09 | -3.41 to -1.77 | 1.70e-27 | 5.11e-27 |
| Hallmark apical junction | GSE123902 | 13 | -1.69 | -2.42 to 1.14 | 1.01e-04 | 1.21e-04 |
| Hallmark apical junction | GSE131907 | 21 | 0.94 | -1.41 to 1.23 | 0.613 | 0.767 |
| Hallmark apical junction | GSE205335 | 21 | -1.45 | -2.00 to 1.36 | 0.002 | 0.003 |
| Hallmark apical junction | GSE189357 | 9 | -1.51 | -1.60 to 1.06 | 0.003 | 0.010 |
| Hallmark apical junction | pooled | 64 | -1.59 | -2.10 to 1.04 | 2.45e-04 | 2.95e-04 |
| KEGG tight junction | GSE123902 | 13 | 1.01 | -1.38 to 1.50 | 0.446 | 0.446 |
| KEGG tight junction | GSE131907 | 21 | 1.68 | 0.72 to 1.87 | 2.69e-04 | 8.06e-04 |
| KEGG tight junction | GSE205335 | 21 | -1.28 | -1.69 to 1.07 | 0.041 | 0.041 |
| KEGG tight junction | GSE189357 | 9 | -1.24 | -1.35 to 1.09 | 0.080 | 0.120 |
| KEGG tight junction | pooled | 64 | 1.32 | -1.18 to 1.67 | 0.045 | 0.045 |

### Leave-one-cohort-out, pooled Q4 ranking (re-TMM)

| dropped | set | n | NES | p |
|---|---|---:|---:|---:|
| GSE123902 | Hallmark IFN-alpha | 27 | -3.24 | 6.37e-19 |
| GSE123902 | Hallmark IFN-gamma | 27 | -3.77 | 1.48e-44 |
| GSE123902 | MHC-I antigen presentation | 27 | -2.81 | 1.50e-08 |
| GSE123902 | Hallmark EMT | 27 | -2.81 | 2.88e-20 |
| GSE123902 | Hallmark apical junction | 27 | -1.45 | 8.63e-04 |
| GSE123902 | KEGG tight junction | 27 | 0.95 | 0.604 |
| GSE131907 | Hallmark IFN-alpha | 23 | -3.48 | 1.05e-27 |
| GSE131907 | Hallmark IFN-gamma | 23 | -4.01 | 7.20e-64 |
| GSE131907 | MHC-I antigen presentation | 23 | -2.92 | 2.40e-12 |
| GSE131907 | Hallmark EMT | 23 | -3.29 | 6.34e-34 |
| GSE131907 | Hallmark apical junction | 23 | -1.87 | 7.97e-07 |
| GSE131907 | KEGG tight junction | 23 | -1.33 | 0.025 |
| GSE205335 | Hallmark IFN-alpha | 23 | -1.99 | 7.48e-06 |
| GSE205335 | Hallmark IFN-gamma | 23 | -2.37 | 2.25e-14 |
| GSE205335 | MHC-I antigen presentation | 23 | -1.99 | 9.66e-04 |
| GSE205335 | Hallmark EMT | 23 | -3.35 | 1.22e-34 |
| GSE205335 | Hallmark apical junction | 23 | -1.59 | 7.13e-05 |
| GSE205335 | KEGG tight junction | 23 | 1.35 | 0.037 |
| GSE189357 | Hallmark IFN-alpha | 29 | -3.37 | 8.36e-25 |
| GSE189357 | Hallmark IFN-gamma | 29 | -3.79 | 1.33e-49 |
| GSE189357 | MHC-I antigen presentation | 29 | -2.81 | 2.51e-10 |
| GSE189357 | Hallmark EMT | 29 | -3.11 | 5.32e-28 |
| GSE189357 | Hallmark apical junction | 29 | -1.70 | 3.61e-05 |
| GSE189357 | KEGG tight junction | 29 | 1.08 | 0.314 |

### CLDN4 put back into junction sets (same ranking)

| set | contrast | NES held out | NES with CLDN4 | delta |
|---|---|---:|---:|---:|
| Hallmark apical junction | q4_vs_q1 | -1.62 | -1.59 | 0.02 |
| Hallmark apical junction | continuous_z | -1.59 | -1.53 | 0.06 |
| KEGG tight junction | q4_vs_q1 | 1.13 | 1.16 | 0.03 |
| KEGG tight junction | continuous_z | 1.32 | 1.40 | 0.08 |
| GO TJ organization | q4_vs_q1 | 1.44 | 1.48 | 0.04 |
| GO TJ organization | continuous_z | 1.60 | 1.71 | 0.11 |
| GO bicellular TJ assembly | q4_vs_q1 | 1.51 | 1.60 | 0.09 |
| GO bicellular TJ assembly | continuous_z | 1.64 | 1.80 | 0.16 |

## 2. GSVA and ssGSEA on the shared-gene malignant pseudobulks

One GSVA (`kcdf=Gaussian`, `maxDiff=TRUE`) and one ssGSEA (`alpha=0.25`,
`normalize=TRUE`) on the within-cohort log2(TMM-CPM+1) matrix restricted to
the 12285 shared genes. Scores are per unit. Tests:

- Within-cohort Spearman of the score vs CLDN4 %pos, then DerSimonian-Laird
  on Fisher z (all four cohorts).
- Within-cohort Welch difference Q4 minus Q1 where both arms have n>=3, and
  the same contrast as Hedges g.
- Pooled OLS `score ~ cohort + Q4` and `score ~ cohort + z(CLDN4 %pos)`.

FDR is BH across the six primary sets inside method x scope x contrast.

### Pooled OLS, Q4 vs Q1 (score units; positive = higher in Q4)

| method | set | scope | effect | 95% CI | p | FDR | I2 |
|---|---|---|---:|---|---:|---:|---|
| GSVA | Hallmark IFN-alpha | pooled_OLS | -0.401 | -0.618 to -0.183 | 7.43e-04 | 0.001 | — |
| GSVA | Hallmark IFN-gamma | pooled_OLS | -0.390 | -0.575 to -0.204 | 1.78e-04 | 5.35e-04 | — |
| GSVA | MHC-I antigen presentation | pooled_OLS | -0.456 | -0.737 to -0.176 | 0.002 | 0.004 | — |
| GSVA | Hallmark EMT | pooled_OLS | -0.301 | -0.443 to -0.159 | 1.58e-04 | 5.35e-04 | — |
| GSVA | Hallmark apical junction | pooled_OLS | -0.148 | -0.289 to -0.008 | 0.039 | 0.046 | — |
| GSVA | KEGG tight junction | pooled_OLS | -0.045 | -0.140 to 0.051 | 0.348 | 0.348 | — |
| SSGSEA | Hallmark IFN-alpha | pooled_OLS | -0.150 | -0.244 to -0.057 | 0.003 | 0.006 | — |
| SSGSEA | Hallmark IFN-gamma | pooled_OLS | -0.142 | -0.220 to -0.065 | 7.54e-04 | 0.005 | — |
| SSGSEA | MHC-I antigen presentation | pooled_OLS | -0.131 | -0.215 to -0.047 | 0.003 | 0.006 | — |
| SSGSEA | Hallmark EMT | pooled_OLS | -0.144 | -0.239 to -0.049 | 0.004 | 0.006 | — |
| SSGSEA | Hallmark apical junction | pooled_OLS | -0.044 | -0.104 to 0.015 | 0.140 | 0.168 | — |
| SSGSEA | KEGG tight junction | pooled_OLS | 0.009 | -0.022 to 0.039 | 0.565 | 0.565 | — |

### Pooled OLS, per SD of CLDN4 %pos

| method | set | scope | effect | 95% CI | p | FDR | I2 |
|---|---|---|---:|---|---:|---:|---|
| GSVA | Hallmark IFN-alpha | pooled_OLS | -0.111 | -0.214 to -0.007 | 0.036 | 0.054 | — |
| GSVA | Hallmark IFN-gamma | pooled_OLS | -0.113 | -0.202 to -0.024 | 0.014 | 0.034 | — |
| GSVA | MHC-I antigen presentation | pooled_OLS | -0.148 | -0.269 to -0.027 | 0.017 | 0.034 | — |
| GSVA | Hallmark EMT | pooled_OLS | -0.105 | -0.164 to -0.047 | 6.11e-04 | 0.004 | — |
| GSVA | Hallmark apical junction | pooled_OLS | -0.047 | -0.098 to 0.004 | 0.071 | 0.086 | — |
| GSVA | KEGG tight junction | pooled_OLS | -0.006 | -0.044 to 0.033 | 0.779 | 0.779 | — |
| SSGSEA | Hallmark IFN-alpha | pooled_OLS | -0.045 | -0.086 to -0.004 | 0.030 | 0.046 | — |
| SSGSEA | Hallmark IFN-gamma | pooled_OLS | -0.043 | -0.077 to -0.008 | 0.016 | 0.033 | — |
| SSGSEA | MHC-I antigen presentation | pooled_OLS | -0.050 | -0.084 to -0.016 | 0.005 | 0.014 | — |
| SSGSEA | Hallmark EMT | pooled_OLS | -0.051 | -0.086 to -0.016 | 0.005 | 0.014 | — |
| SSGSEA | Hallmark apical junction | pooled_OLS | -0.014 | -0.033 to 0.006 | 0.167 | 0.201 | — |
| SSGSEA | KEGG tight junction | pooled_OLS | 0.006 | -0.006 to 0.018 | 0.299 | 0.299 | — |

### DL meta of within-cohort Spearman rho

| method | set | scope | effect | 95% CI | p | FDR | I2 |
|---|---|---|---:|---|---:|---:|---|
| GSVA | Hallmark IFN-alpha | DL_meta | -0.319 | -0.539 to -0.059 | 0.017 | 0.034 | 0% |
| GSVA | Hallmark IFN-gamma | DL_meta | -0.399 | -0.643 to -0.082 | 0.015 | 0.034 | 32% |
| GSVA | MHC-I antigen presentation | DL_meta | -0.249 | -0.483 to 0.017 | 0.066 | 0.080 | 0% |
| GSVA | Hallmark EMT | DL_meta | -0.563 | -0.814 to -0.136 | 0.013 | 0.034 | 68% |
| GSVA | Hallmark apical junction | DL_meta | -0.269 | -0.499 to -0.004 | 0.047 | 0.070 | 0% |
| GSVA | KEGG tight junction | DL_meta | -0.084 | -0.341 to 0.186 | 0.546 | 0.546 | 0% |
| SSGSEA | Hallmark IFN-alpha | DL_meta | -0.370 | -0.600 to -0.083 | 0.013 | 0.073 | 18% |
| SSGSEA | Hallmark IFN-gamma | DL_meta | -0.461 | -0.776 to 0.036 | 0.068 | 0.101 | 72% |
| SSGSEA | MHC-I antigen presentation | DL_meta | -0.348 | -0.607 to -0.022 | 0.037 | 0.074 | 33% |
| SSGSEA | Hallmark EMT | DL_meta | -0.359 | -0.607 to -0.048 | 0.024 | 0.073 | 28% |
| SSGSEA | Hallmark apical junction | DL_meta | -0.119 | -0.373 to 0.151 | 0.387 | 0.464 | 0% |
| SSGSEA | KEGG tight junction | DL_meta | 0.079 | -0.191 to 0.337 | 0.570 | 0.570 | 0% |

### DL meta of Hedges g (Q4 vs Q1; cohorts with both arms n>=3)

| method | set | scope | effect | 95% CI | p | FDR | I2 |
|---|---|---|---:|---|---:|---:|---|
| GSVA | Hallmark IFN-alpha | DL_meta | -1.396 | -2.586 to -0.207 | 0.021 | 0.032 | 57% |
| GSVA | Hallmark IFN-gamma | DL_meta | -1.631 | -3.016 to -0.245 | 0.021 | 0.032 | 65% |
| GSVA | MHC-I antigen presentation | DL_meta | -1.249 | -2.193 to -0.306 | 0.009 | 0.028 | 35% |
| GSVA | Hallmark EMT | DL_meta | -1.135 | -1.857 to -0.414 | 0.002 | 0.012 | 0% |
| GSVA | Hallmark apical junction | DL_meta | -0.652 | -1.330 to 0.026 | 0.059 | 0.071 | 0% |
| GSVA | KEGG tight junction | DL_meta | -0.171 | -0.973 to 0.630 | 0.675 | 0.675 | 29% |
| SSGSEA | Hallmark IFN-alpha | DL_meta | -1.245 | -2.567 to 0.077 | 0.065 | 0.117 | 66% |
| SSGSEA | Hallmark IFN-gamma | DL_meta | -1.439 | -3.038 to 0.161 | 0.078 | 0.117 | 75% |
| SSGSEA | MHC-I antigen presentation | DL_meta | -1.235 | -2.272 to -0.198 | 0.020 | 0.059 | 46% |
| SSGSEA | Hallmark EMT | DL_meta | -0.839 | -1.532 to -0.147 | 0.018 | 0.059 | 0% |
| SSGSEA | Hallmark apical junction | DL_meta | -0.415 | -1.084 to 0.254 | 0.224 | 0.269 | 0% |
| SSGSEA | KEGG tight junction | DL_meta | 0.398 | -0.616 to 1.413 | 0.442 | 0.442 | 53% |

## How to read the two layers

fgsea asks whether set members sit at one end of the patient-level t ranking.
GSVA and ssGSEA ask whether that patient's enrichment score is lower when
malignant CLDN4 %pos is high.

On the pooled Q4 contrast, GSVA and ssGSEA both put IFN-alpha, IFN-gamma,
MHC-I antigen presentation, and Hallmark EMT lower in CLDN4-high (FDR < 0.05
inside each method). That matches the negative fgsea NES. The Hedges g
meta-analysis is wider: ssGSEA IFN intervals cross zero and I2 is high, so
the score shift is not the same size in every cohort.

KEGG tight junction does not clear the score tests. Its Q4 fgsea NES is
positive and not FDR < 0.05. This run is not tight-junction up. Hallmark
apical junction is negative on fgsea and on the GSVA Q4 model, and not
FDR < 0.05 on ssGSEA. It is not a junction-barrier score.

## What this is

A second look at the same malignant pseudobulks, with rank enrichment
(fgsea) and sample-wise enrichment scores (GSVA, ssGSEA) instead of a mean
logFC across a merged IFN list. It does not replace the T/NK Spearman, and
it is not a cell-level test.

## What this is not

- Not a mega-merge and not TACSTD2.
- Not proof that CLDN4 causes the pathway shift.
- Not a DoRothEA / PROGENy activity estimate.
- Not cell-level GSEA p-values.
- Genome-wide significance is not claimed. The tested family is the six
  primary sets.

## Figures

- `results/figures/forest_fgsea_nes_pooled.png` — pooled NES, both contrasts
- `results/figures/forest_fgsea_nes_q4q1.png` — NES by cohort
- `results/figures/forest_fgsea_nes_continuous.png`
- `results/figures/forest_gsva_q4q1.png` / `forest_ssgsea_q4q1.png`
- `results/figures/forest_gsva_spearman.png` / `forest_ssgsea_spearman.png`
- `results/figures/forest_score_hedges_g.png`
- `results/figures/box_ssgsea_q4q1_patient.png`

## Reproduce

```bash
Rscript methods/concordant4_cldn4_gsva_fgsea/analyze.R
```

Requires R packages GSVA (>= 2.0), fgsea, ggplot2, data.table, jsonlite.
`N_BOOT` and `NPERM_BOOT` override the patient-bootstrap size.
Seed 20260921.

