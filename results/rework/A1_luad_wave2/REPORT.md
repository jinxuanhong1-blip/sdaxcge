# REWORK A1 wave 2 — TCGA-LUAD double residual (ABSOLUTE + epithelium) and ImmuneScore residual

**Self-contained. Public data only. No fabricated numbers.**

**Cohort:** TCGA-LUAD primary tumors (`-01`) only. Not LUSC. Not OncoSG. Not the pooled n ≈ 1031.

**Question:** Prior LUAD-only A1 was weak after ABSOLUTE (|ρ| ≤ 0.22 markers, PR #81; |ρ| ≤ 0.13 xCell/MCP/GEP18/ESTIMATE, PR #169). OncoSG and LUSC *did* support a negative TACSTD2–immune association after purity. Can a second residual — epithelium or ESTIMATE ImmuneScore — recover the user's implied **strong inverse** on LUAD?

## Verdict: STILL_WEAK

- One-sided double residual (requested): most negative ρ = -0.072; |ρ| stays inside the prior LUAD marker bound of 0.22 and far from OncoSG −0.30.
- Two-sided partial | ABS+KRT/EPCAM: most negative ρ = -0.082; |ρ| stays inside the prior LUAD marker bound of 0.22 and far from OncoSG −0.30.
- ImmuneScore residual, one-sided: most negative ρ = -0.140; |ρ| stays inside the prior LUAD marker bound of 0.22 and far from OncoSG −0.30.
- ImmuneScore residual, two-sided: most negative ρ = -0.211; |ρ| stays inside the prior LUAD marker bound of 0.22 and far from OncoSG −0.30.

ABSOLUTE-only two-sided partials **reproduce** the histology / deconv reworks: CD8 −0.104, CYT −0.116, GEP18 z-mean −0.013, xCell CD8 −0.126 (n=502). TACSTD2 ⟂ ABSOLUTE (ρ = 0.007, p = 0.869).

The epithelial double residual **weakens** those already-small inverses (TACSTD2 vs KRT/EPCAM ρ = 0.364). GEP18 even flips to a small positive. Epithelium was not masking a strong LUAD inverse.

Scan-wide strongest primary ρ is `ABS_IMMUNE` two-sided CYT = −0.217 (n=502). That is still `WEAK_NEGATIVE` under the locked −0.22 / −0.30 rules. It is not hidden; it is not OncoSG.

**Implied strong inverse, locked before results.** OncoSG A1 after published purity (PR #139, n=169) is the user's working example of a supported inverse: CD8A partial ρ = −0.309, GEP18 = −0.349. Thresholds:

| Call | Rule |
|---|---|
| `STRONG_INVERSE` | ρ ≤ **−0.30** (OncoSG CD8/GEP) |
| `REACHES_LUSC_BOUND` | −0.30 < ρ ≤ **−0.22** (LUSC CD8 / prior LUAD marker bound) |
| `WEAK_NEGATIVE` | −0.22 < ρ < 0 |
| `NOT_INVERSE` | ρ ≥ 0 |

n ≈ 1031 in the original claim is the **pooled TCGA LUAD+LUSC** ESTIMATE-complete set (PR #89), not LUAD-only. This file stays LUAD-only. Complete-case n is in every table cell (~502 with ABSOLUTE).

## Why one-sided and two-sided are both shown

The request is: residualize **TACSTD2** on (1) ABSOLUTE and (2) epithelial fraction / KRT/EPCAM, **then** correlate with CD8 / GEP18 / CYT / xCell CD8.

That is a **one-sided** residual. Immune scores still carry epithelial / purity composition. Because KRT/EPCAM is expected to anti-correlate with CD8 in bulk RNA, one-sided residualization can *inflate* a negative ρ. The **two-sided** rank residual (standard multi-covariate partial Spearman) is the fair test of “TACSTD2 vs immune beyond tumor content and epithelium.” Both are reported. The overall tag uses the two-sided result as the claim-level answer and states the one-sided number separately.

Residualizing TACSTD2 on ESTIMATE ImmuneScore and then correlating with CD8/GEP18/CYT/xCell CD8 is valid (different features). Correlating that residual with ImmuneScore itself would be circular and is not used as a primary readout. ImmuneScore is strongly anti-correlated with ABSOLUTE (ρ = -0.616) by construction of leukocyte content.

TIMER2 xCell has **no Epithelial-cells column** in this freeze (keratinocyte column: also absent). Epithelial fraction is therefore the KRT/EPCAM gene score and ESTIMATE TumorPurity, not an xCell epithelial fraction. That absence is reported, not filled in.

## Analysis set

| Filter | n |
|---|---:|
| Xena HiSeqV2 primary tumors (`-01`) | 515 |
| + ABSOLUTE purity | 502 |
| + KRT/EPCAM score (same matrix) | 502 |
| + official ESTIMATE ImmuneScore | 515 |
| + TIMER2 xCell CD8 | 515 |

Primary tumors only. One row per 15-character barcode. Replicate aliquots averaged.

## Definitions (locked)

### TACSTD2
UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` **log2(RSEM normalized_count + 1)**. Gene symbol `TACSTD2`. Not protein. Not ADC. Not ICI.

### Samples
TCGA-LUAD only. Sample type `01`. No LUSC, no OncoSG, no normals.

### Covariates
1. **ABSOLUTE purity** — PanCanAtlas GDC `4f277128-f793-4354-a13d-30cc7fe9f6b5`, column `purity`.
2. **KRT/EPCAM score (primary epithelium)** — mean log2(RSEM+1) of `EPCAM, KRT7, KRT8, KRT18, KRT19`. Present: EPCAM, KRT7, KRT8, KRT18, KRT19. Missing (not imputed): none. **TACSTD2 is not in the score.**
3. **ESTIMATE TumorPurity (sensitivity epithelium)** — Yoshihara 2013 transform `cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)` from the official MD Anderson LUAD RNAseqV2 table. This is RNA tumor fraction, not ABSOLUTE.
4. **ESTIMATE ImmuneScore** — official `Immune_score` from the same table. Used only as the requested TACSTD2 residual covariate, not as a primary immune *outcome*.
5. **xCell epithelial** — used only if a TIMER2 `*_XCELL` epithelial column exists; otherwise skipped (see above).

Sensitivity epithelium: `EPCAM, KRT7, KRT8, KRT18, KRT19, CDH1` mean; and mean of EPCAM + every `KRT[0-9]+` gene in HiSeqV2 (56 genes).

### Immune outcomes (primary)
| Name | Definition |
|---|---|
| CD8 | `CD8A` log2(RSEM+1) on the same Xena matrix |
| CYT | Rooney *Cell* 2015 cytolytic score = mean(`GZMA`, `PRF1`) on log2 data |
| GEP18 | Ayers *JCI* 2017 18 genes, **unweighted within-LUAD z-mean** (histology-rework definition, PR #107) |
| xCell CD8 | TIMER2.0 immunedeconv column `T cell CD8+_XCELL` |

Sensitivity outcomes (not in the verdict): `CD8B`, GEP18 ssGSEA (Barbie/GSVA tau=0.25, ranks 1..10000), GEP18 log-mean, official ESTIMATE ImmuneScore as a *Y* (only for ABS/KRT models).

### Estimators
- **onesided_rank (user-requested primary):** OLS of rank(TACSTD2) on rank(covariates); Pearson of that residual vs rank(immune). df = n − 2 − k.
- **onesided_value:** OLS of TACSTD2 on raw covariates; Spearman of residual vs raw immune. Same df.
- **twosided_rank:** both sides residualized on ranked covariates (partial Spearman). Same df.
- Unadjusted Spearman is the no-covariate baseline.
- 95% Fisher-z CI uses variance `1/(n − k − 3)`.
- BH-FDR across the 4 primary immune features, within each model × estimator.

## Covariate context (unadjusted Spearman)

| Pair | ρ | p | n |
|---|---:|---:|---:|
| TACSTD2 vs ABSOLUTE | 0.007 | 0.869 | 502 |
| TACSTD2 vs KRT/EPCAM | 0.364 | 1.47e-17 | 515 |
| TACSTD2 vs ESTIMATE ImmuneScore | 0.036 | 0.41 | 515 |
| ABSOLUTE vs KRT/EPCAM | 0.023 | 0.61 | 502 |
| ABSOLUTE vs ImmuneScore | -0.616 | 7.71e-54 | 502 |

If TACSTD2 is nearly orthogonal to ABSOLUTE (as in PR #81 / #169), purity-only adjustment cannot create a strong inverse. A second residual can only help if TACSTD2 shares variance with epithelium (or ImmuneScore) that was *masking* an inverse with CD8/GEP/CYT.

## Primary results

### Unadjusted

| Feature | ρ | p | FDR (4) | n | 95% CI | Call |
|---|---:|---:|---:|---:|---|---|
| `CD8` | -0.096 | 0.0294 | 0.0392 | 515 | -0.181 to -0.010 | WEAK_NEGATIVE |
| `CYT` | -0.112 | 0.0108 | 0.0217 | 515 | -0.197 to -0.026 | WEAK_NEGATIVE |
| `GEP18_zmean` | -0.022 | 0.62 | 0.62 | 515 | -0.108 to 0.065 | WEAK_NEGATIVE |
| `xCell_T_cell_CD8` | -0.117 | 0.00809 | 0.0217 | 515 | -0.201 to -0.030 | WEAK_NEGATIVE |

### ABSOLUTE only (prior A1 replication, one-sided rank residual)

| Feature | ρ | p | FDR (4) | n | 95% CI | Call |
|---|---:|---:|---:|---:|---|---|
| `CD8` | -0.092 | 0.0397 | 0.0529 | 502 | -0.178 to -0.004 | WEAK_NEGATIVE |
| `CYT` | -0.103 | 0.0216 | 0.0432 | 502 | -0.189 to -0.015 | WEAK_NEGATIVE |
| `GEP18_zmean` | -0.010 | 0.817 | 0.817 | 502 | -0.098 to 0.077 | WEAK_NEGATIVE |
| `xCell_T_cell_CD8` | -0.119 | 0.00744 | 0.0298 | 502 | -0.205 to -0.032 | WEAK_NEGATIVE |

### Double residual: ABSOLUTE + KRT/EPCAM — one-sided (requested)

| Feature | ρ | p | FDR (4) | n | 95% CI | Call |
|---|---:|---:|---:|---:|---|---|
| `CD8` | -0.059 | 0.19 | 0.255 | 502 | -0.146 to 0.029 | WEAK_NEGATIVE |
| `CYT` | -0.072 | 0.106 | 0.255 | 502 | -0.159 to 0.015 | WEAK_NEGATIVE |
| `GEP18_zmean` | 0.024 | 0.591 | 0.591 | 502 | -0.064 to 0.112 | NOT_INVERSE |
| `xCell_T_cell_CD8` | -0.059 | 0.192 | 0.255 | 502 | -0.145 to 0.029 | WEAK_NEGATIVE |

### Double residual: ABSOLUTE + KRT/EPCAM — two-sided partial (fair test)

| Feature | ρ | p | FDR (4) | n | 95% CI | Call |
|---|---:|---:|---:|---:|---|---|
| `CD8` | -0.067 | 0.137 | 0.215 | 502 | -0.153 to 0.021 | WEAK_NEGATIVE |
| `CYT` | -0.082 | 0.066 | 0.215 | 502 | -0.169 to 0.005 | WEAK_NEGATIVE |
| `GEP18_zmean` | 0.030 | 0.503 | 0.503 | 502 | -0.058 to 0.117 | NOT_INVERSE |
| `xCell_T_cell_CD8` | -0.063 | 0.161 | 0.215 | 502 | -0.150 to 0.025 | WEAK_NEGATIVE |

### ESTIMATE ImmuneScore residual of TACSTD2 — one-sided

| Feature | ρ | p | FDR (4) | n | 95% CI | Call |
|---|---:|---:|---:|---:|---|---|
| `CD8` | -0.122 | 0.00549 | 0.00732 | 515 | -0.207 to -0.036 | WEAK_NEGATIVE |
| `CYT` | -0.140 | 0.00152 | 0.00319 | 515 | -0.223 to -0.054 | WEAK_NEGATIVE |
| `GEP18_zmean` | -0.054 | 0.222 | 0.222 | 515 | -0.140 to 0.033 | WEAK_NEGATIVE |
| `xCell_T_cell_CD8` | -0.139 | 0.0016 | 0.00319 | 515 | -0.223 to -0.053 | WEAK_NEGATIVE |

### ESTIMATE ImmuneScore residual — two-sided partial

| Feature | ρ | p | FDR (4) | n | 95% CI | Call |
|---|---:|---:|---:|---:|---|---|
| `CD8` | -0.177 | 5.60e-05 | 8.45e-05 | 515 | -0.259 to -0.092 | WEAK_NEGATIVE |
| `CYT` | -0.211 | 1.42e-06 | 5.67e-06 | 515 | -0.292 to -0.127 | WEAK_NEGATIVE |
| `GEP18_zmean` | -0.114 | 0.00995 | 0.00995 | 515 | -0.198 to -0.027 | WEAK_NEGATIVE |
| `xCell_T_cell_CD8` | -0.175 | 6.34e-05 | 8.45e-05 | 515 | -0.258 to -0.090 | WEAK_NEGATIVE |

## All models, one-sided rank residual (primary four)

| Model | CD8 | CYT | GEP18 z-mean | xCell CD8 | most negative |
|---|---:|---:|---:|---:|---:|
| `unadjusted` | -0.096 | -0.112 | -0.022 | -0.117 | -0.117 |
| `ABS_only` | -0.092 | -0.103 | -0.010 | -0.119 | -0.119 |
| `KRT_only` | -0.055 | -0.074 | 0.023 | -0.049 | -0.074 |
| `ABS_KRT` | -0.059 | -0.072 | 0.024 | -0.059 | -0.072 |
| `ABS_KRT_CDH1` | -0.046 | -0.060 | 0.031 | -0.045 | -0.060 |
| `ABS_allKRT` | -0.046 | -0.068 | -0.008 | -0.026 | -0.068 |
| `ABS_ESTPUR` | -0.108 | -0.120 | -0.032 | -0.134 | -0.134 |
| `IMMUNE` | -0.122 | -0.140 | -0.054 | -0.139 | -0.140 |
| `ABS_IMMUNE` | -0.132 | -0.145 | -0.058 | -0.157 | -0.157 |

## All models, two-sided rank partial (primary four)

| Model | CD8 | CYT | GEP18 z-mean | xCell CD8 | most negative |
|---|---:|---:|---:|---:|---:|
| `unadjusted` | -0.096 | -0.112 | -0.022 | -0.117 | -0.117 |
| `ABS_only` | -0.104 | -0.116 | -0.013 | -0.126 | -0.126 |
| `KRT_only` | -0.055 | -0.075 | 0.023 | -0.050 | -0.075 |
| `ABS_KRT` | -0.067 | -0.082 | 0.030 | -0.063 | -0.082 |
| `ABS_KRT_CDH1` | -0.052 | -0.068 | 0.039 | -0.049 | -0.068 |
| `ABS_allKRT` | -0.052 | -0.077 | -0.010 | -0.028 | -0.077 |
| `ABS_ESTPUR` | -0.136 | -0.154 | -0.052 | -0.152 | -0.154 |
| `IMMUNE` | -0.177 | -0.211 | -0.114 | -0.175 | -0.211 |
| `ABS_IMMUNE` | -0.190 | -0.217 | -0.121 | -0.199 | -0.217 |

## Does |ρ| reach the implied strong inverse?

| Test | Most negative primary ρ | Strong (≤ −0.30)? | LUSC bound (≤ −0.22)? |
|---|---:|---|---|
| Unadjusted | -0.117 | False | False |
| ABS only, one-sided | -0.119 | False | False |
| ABS+KRT, one-sided | -0.072 | False | False |
| ABS+KRT, two-sided | -0.082 | False | False |
| ImmuneScore residual, one-sided | -0.140 | False | False |
| ImmuneScore residual, two-sided | -0.211 | False | False |

Largest |ρ| among primary four, ABS+KRT one-sided: 0.072. Two-sided: 0.082. Scan-wide strongest (any model, primary four): -0.217 (`ABS_IMMUNE` two-sided CYT). Still above −0.22.

## Honest limits

1. **LUAD-only.** OncoSG (n=169, East-Asian surgical LUAD) and TCGA-LUSC are other cohorts. Their support does not transfer by residualization.
2. **n is not 1031.** Pooled LUAD+LUSC ESTIMATE n ≈ 1014–1034 (PR #89). This table is ~502 LUAD primaries with ABSOLUTE.
3. **KRT/EPCAM is a transcript score, not a cell fraction.** ESTIMATE TumorPurity is the published RNA tumor-fraction transform; ABSOLUTE is DNA. They are not interchangeable.
4. **xCell CD8 is TIMER2, not a from-scratch Xena xCell run.**
5. **GEP18 z-mean is not the Merck NanoString TIS.** ssGSEA is a sensitivity row.
6. **One-sided residual vs raw immune is not a partial correlation.** If it looks stronger than two-sided, composition in Y is the first explanation.
7. **No causality, no protein, no ICI outcome.**

## Reproduce

```
pip install -r requirements.txt
python scripts/rework_A1_luad_wave2.py
```

Downloads public tables into `data/` (gitignored) on first run. Signature lists are bundled under `data/signatures/`.

## Files

- `REPORT.md` — this writeup
- `DEFINITIONS.md` — locked definitions
- `correlations.tsv` — every model × estimator × feature
- `primary_results.tsv` — unadjusted / ABS / ABS+KRT / ImmuneScore × primary four
- `context_correlations.tsv` — TACSTD2 / purity / epithelium / ImmuneScore pairwise
- `sample_table.tsv` — per-sample values
- `gene_coverage.tsv`
- `summary.json` / `provenance.json`
- `figures/bar_residual_models.png`
- `figures/forest_double_residual.png`
- `figures/scatter_double_residual.png`
- `figures/scatter_tacstd2_vs_covariates.png`

## Data

- Expression: https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz
- Purity: https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5
- ESTIMATE: https://ibl.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt
- TIMER2: https://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz
