# GSE205335 — TACSTD2 (TROP2) in malignant cells: ICI responders vs non-responders

## TL;DR (honest summary)

- **No significant difference in malignant-cell TACSTD2 between responders and non-responders** (per-patient Mann-Whitney: p = 0.96 for mean normalized expression, p = 0.71 for % positive cells, p = 0.64 for pseudobulk CPM; PR-vs-PD-only sensitivity p = 0.95). With n = 6 R vs 10 NR patients, only a large effect could have been detected; this is an underpowered null, not proof of no association.
- **TACSTD2 is strongly malignant-restricted**: malignant cells vastly exceed T/NK cells in every one of 22 evaluable patients (paired Wilcoxon p = 4.8e-7; median mean-log1p(CP10K) 0.93 vs 0.02; median 58.4% vs 1.2% positive cells). TROP2 is a clean tumor-vs-lymphocyte marker in this cohort regardless of ICI response.
- **Important dropout asymmetry**: 3 of 9 responder patients (P2001, P2009, P3032) and 1 non-responder (P2016) had **zero** cells annotated as malignant, so they could not contribute to the tumor-cell comparison. Absence of captured tumor cells in responders may itself reflect treatment effect or sampling and shrinks the R arm.

## Dataset

- **GSE205335**: "Single-cell transcriptome profiles of tumor tissues from lung cancer patients receiving immune checkpoint inhibitors" (Ahn / Lee, Catholic Univ. of Korea; public Nov 2024). 33 samples (metastatic LN, liver, lung/bronchus, effusion) from 26 lung cancer patients (ADC/SQ, mostly stage IV) on ICI.
- Inputs used (author-processed, from GEO supplementary files):
  - `GSE205335_Lung_IO_UMI_matrix.rds.gz` — dgCMatrix, 33,714 genes x 96,505 cells (raw UMI).
  - `GSE205335_Lung_IO_CellIdentity.txt.gz` — author cell annotations; malignant cells are the **authors' own call** (`lineage.sub == "Malignant cells"`, n = 28,512), T/NK from `lineage.total == "T/NK cells"` (n = 39,875).
- Response from GEO sample characteristics (RECIST): **R = PR** (no CR in cohort), **NR = SD/PD**, **NE excluded** from the response comparison (6 patients).

## Methods

1. Per-cell TACSTD2 UMI and total UMI extracted from the RDS in R; normalized as log1p(counts/total x 10^4) ("mean_log1p_cp10k").
2. Samples pooled per patient (patient = unit of analysis; several patients have 2 samples). Patients required >= 20 cells in a compartment to contribute.
3. Metrics per patient x compartment: mean log1p(CP10K), % TACSTD2+ cells, pseudobulk CPM.
4. Tests: two-sided Mann-Whitney U (R vs NR, malignant cells); two-sided Wilcoxon signed-rank (malignant vs T/NK, paired within patient, all 22 evaluable patients incl. NE).
5. Sanity checks passed: EPCAM detected in 69.3% of malignant vs 2.5% of T/NK cells; PTPRC (CD45) in 4.2% vs 80.9%.

## Files

| File | Contents |
|---|---|
| `stats_results.txt` | All test statistics and exclusion notes |
| `per_patient_tacstd2.csv` | Per-patient x compartment metrics (n cells, % positive, mean log1p CP10K, pseudobulk CPM) |
| `per_sample_tacstd2.csv` | Same at sample level (GSM, tissue, platform) |
| `tacstd2_summary.png` | Panels: malignant TACSTD2 by response; % positive by response; paired malignant vs T/NK |
| `tacstd2_per_patient_bars.png` | Per-patient malignant TACSTD2 bar chart, colored by response |

Reproduction: `scripts/extract_tacstd2.R` (RDS -> per-cell CSV; note the GEO file is double-gzipped — `gunzip` once before `readRDS`) then `scripts/analyze_tacstd2.py`.

## Caveats (read before citing)

- **Small n**: 6 R vs 10 NR patients after exclusions. The R-vs-NR comparison is descriptive; do not over-interpret the null.
- **Responder tumor-cell dropout**: 3/9 R patients had no malignant cells captured, a potential selection bias in the R arm.
- **Mixed tissues and platforms**: LN, liver, lung, effusion samples; 10x 3' and 5' chemistries are pooled. Platform/tissue are confounded with patient and not adjusted for (n too small to model them).
- **Treatment timing not exposed in GEO metadata**: samples may vary in timing relative to ICI; this cannot be controlled here (raw data are at EGA EGAD00001008703 under controlled access).
- **Malignant calls are the authors'**; no independent CNV inference was performed here. Marker sanity checks (EPCAM/PTPRC) are consistent with correct labels.
- RECIST for P1084/P1018/P1063/P1079/P1072/P0031 is "NE"; they are shown in plots for context but excluded from R-vs-NR tests. P0031's normal-lung sample (LUNG_N31) contributes no malignant cells by definition; its cells enter only via the authors' annotations.
