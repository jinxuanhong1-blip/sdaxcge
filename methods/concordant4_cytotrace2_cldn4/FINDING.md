# Finding — concordant-4 malignant CLDN4 vs CytoTRACE2 and stemness

ADDITIVE. **CLDN4 only.** Same four public sets as the locked concordant-4 patient table: **GSE123902 + GSE131907 + GSE205335 + GSE189357**. Not GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526. No TACSTD2∩CLDN4 dual-high. This does **not** re-audit the T/NK exclusion ρ.

Primary potency = **CytoTRACE2** (`cytotrace2-py` 1.1.0.4, Kang et al. *Nat Methods* 2025), fit **inside each dataset**. Score 0 = differentiated, 1 = totipotent. Stemness is paired, not a substitute: **Ben-Porath ES1** and **Wong ESC core** (MSigDB c2.cgp v2023.2.Hs), mean z of log1p(CP10k), with CLDN4, EPCAM, TACSTD2, and the barrier/keratin genes removed. Barrier/keratin = KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**). Inferential unit = patient / donor / sample. Primary number = DerSimonian–Laird meta of the four within-cohort Spearman ρ. Cell-level ρ is not the claim.

**Question.** Are CLDN4-high malignant cells more differentiated and more barrier-like? **Answer:** CLDN4-high is barrier-like (DL ρ=0.350, p=0.0078, q=0.0273, I²=0.0%). The CytoTRACE2 Differentiated fraction rises with CLDN4 (ρ=0.395, p=0.0024, q=0.0166, I²=0.0%). Continuous CytoTRACE2 is lower where CLDN4 is higher, but the four-cohort meta is not significant (ρ=-0.303, p=0.0703, q=0.1072, I²=33.0%, 95% CI includes 0). Embryonic stemness does not fall with CLDN4 (Ben-Porath ES1 ρ=0.239, p=0.0766; Wong ESC ρ=0.069, p=0.6158).

## Verdict

DL meta CLDN4 vs CytoTRACE2: k=4, N=65, ρ=-0.303, 95% CI [-0.573, 0.026], p=0.0703, q=0.1072, I²=33.0%.
DL meta CLDN4 vs Ben-Porath ES1 stemness: ρ=0.239, p=0.0766, q=0.1072, I²=0.0%.
DL meta CLDN4 vs Wong ESC stemness: ρ=0.069, p=0.6158, q=0.6158, I²=0.0%.
DL meta CLDN4 vs barrier/keratin (CLDN4 excluded): ρ=0.350, p=0.0078, q=0.0273, I²=0.0%.
CytoTRACE2 vs ES1 (do the two potency readouts agree): ρ=0.367, p=0.0728, q=0.1072.
CytoTRACE2 vs barrier: ρ=-0.243, p=0.3439, q=0.4012.
CLDN4 vs fraction Differentiated: ρ=0.395, p=0.0024, q=0.0166.

Predicted signs if CLDN4-high is differentiated / barrier-like: CytoTRACE2 and stemness **negative**, barrier and Differentiated fraction **positive**, CytoTRACE2 vs stemness **positive**, CytoTRACE2 vs barrier **negative**.

**What holds.** Between patients, CLDN4 tracks barrier/keratin with CLDN4 held out of the score (CLDN4 vs barrier/keratin (no CLDN4): ρ=0.350, p=0.0078, q=0.0273, I²=0.0%). Between patients, CLDN4 tracks a higher CytoTRACE2 Differentiated fraction (CLDN4 vs frac Differentiated: ρ=0.395, p=0.0024, q=0.0166, I²=0.0%). GSE131907 carries the cohort-level signal; the other three cohorts are the same sign and not significant alone. Within a tumor, CLDN4-high cells are barrier-higher (n=65, Δmed=+0.324, p=3.53e-12).

**What does not hold.** Continuous CytoTRACE2 vs CLDN4 is negative but the four-cohort meta CI includes 0 (CLDN4 vs CytoTRACE2: ρ=-0.303, p=0.0703, q=0.1072, I²=33.0%). Cohorts: GSE123902 ρ=-0.467 (p=0.1076, n=13); GSE131907 ρ=-0.512 (p=0.0177, n=21); GSE205335 ρ=-0.248 (p=0.2660, n=22); GSE189357 ρ=0.367 (p=0.3317, n=9). Leave-one-out dropping GSE189357 is ρ=-0.402, p=0.0035. That sensitivity is not the primary. The within-patient CytoTRACE2 shift is null (Δmed=-0.001, p=0.5384). Do not quote it as a within-tumor potency drop. Stemness does not mark CLDN4-high cells as less stem-like. Ben-Porath ES1 vs CLDN4 is positive (CLDN4 vs Ben-Porath ES1 stemness: ρ=0.239, p=0.0766, q=0.1072, I²=0.0%). Wong ESC vs CLDN4 is null (CLDN4 vs Wong ESC stemness: ρ=0.069, p=0.6158, q=0.6158, I²=0.0%). Within a tumor, CLDN4-high cells score higher on both (ES1 Δmed=+0.072, p=8.17e-08; Wong Δmed=+0.069, p=7.21e-06). Wong stemness does track CytoTRACE2 (sensitivity DL ρ=0.495, p=7.65e-05), so the stemness score is coupled to potency and CLDN4 is not on that axis. Gulati 2020 gene-count CytoTRACE goes the other way on the paired test (Δmed=+0.179, p=6.14e-11). Same disagreement as the winning-pair potency folder. This is why Gulati is not the primary. CLDN4 %pos vs CytoTRACE2 is null. Do not substitute the T/NK %pos score for the mean used here.

## Honest n

- Analysis cells after QC and cap ≤200/unit: **n_cells = 11902** GSE123902 2103, GSE131907 3857, GSE205335 4142, GSE189357 1800.
- Units: **n_units = 65** GSE123902 13, GSE131907 21, GSE205335 22, GSE189357 9. Eligible (≥20 malignant cells): **n = 65**.
- Paired within-unit CLDN4-high vs low (≥8 cells/arm): **n = 65**.
- GSE131907 / GSE205335 malignant = **author label**, not a new CNV call. GSE131907 author-malignant catalog 24784; tLung author-malignant = 0. GSE205335 author-malignant non-normal catalog 28512.
- GSE123902 and GSE189357 malignant = marker gate `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`. **Not CNV.** NORMAL lung from GSE123902 is dropped. Gate counts are in `extract_audit.tsv`.
- CytoTRACE2 potency categories (cells): {'Differentiated': 8897, 'Unipotent': 2160, 'Oligopotent': 625, 'Multipotent': 220}.
- Cells with <500 genes (CytoTRACE2 prefers 500–1000): 0.068.
- Stemness gene coverage: {'GSE123902': {'es1_used': 358, 'wong_used': 325, 'core9_detected_frac': 0.761, 'frac_genes_lt_500': 0.114}, 'GSE131907': {'es1_used': 361, 'wong_used': 325, 'core9_detected_frac': 0.883, 'frac_genes_lt_500': 0.019}, 'GSE205335': {'es1_used': 373, 'wong_used': 331, 'core9_detected_frac': 0.874, 'frac_genes_lt_500': 0.08}, 'GSE189357': {'es1_used': 359, 'wong_used': 328, 'core9_detected_frac': 0.59, 'frac_genes_lt_500': 0.089}}.
- GSE148071 not used. Dual-high not used. T/NK ρ not re-fit.

## Locked choices

- CytoTRACE2 per dataset, human model, `max_cores=1`, seed 14. Not a joint KNN.
- Gulati 2020 gene-count CytoTRACE is sensitivity only.
- CLDN4 = patient-mean log1p(CP10k). %pos is sensitivity.
- Stemness z-scores are within dataset.
- Primary family = seven DL metas, BH inside that list.
- Paired test = within-unit CLDN4 tertile (median split if the unit is small).
- Cap ≤200 cells/unit after QC (author sets: cap then QC, matching the winning-pair potency run).

## Primary (DL meta of within-cohort Spearman, BH inside this list)

| Contrast | k | N | ρ | 95% CI | p | q | I² |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: |
| CLDN4 vs CytoTRACE2 | 4 | 65 | -0.303 | [-0.573, 0.026] | 0.0703 | 0.1072 | 33.0% |
| CLDN4 vs Ben-Porath ES1 stemness | 4 | 65 | 0.239 | [-0.026, 0.472] | 0.0766 | 0.1072 | 0.0% |
| CLDN4 vs Wong ESC stemness | 4 | 65 | 0.069 | [-0.198, 0.326] | 0.6158 | 0.6158 | 0.0% |
| CLDN4 vs barrier/keratin (no CLDN4) | 4 | 65 | 0.350 | [0.096, 0.561] | 0.0078 | 0.0273 | 0.0% |
| CytoTRACE2 vs Ben-Porath ES1 | 4 | 65 | 0.367 | [-0.036, 0.668] | 0.0728 | 0.1072 | 55.7% |
| CytoTRACE2 vs barrier/keratin (no CLDN4) | 4 | 65 | -0.243 | [-0.642, 0.259] | 0.3439 | 0.4012 | 70.1% |
| CLDN4 vs frac Differentiated | 4 | 65 | 0.395 | [0.147, 0.596] | 0.0024 | 0.0166 | 0.0% |

## Within-cohort Spearman (not BH)

| Contrast | dataset | n | ρ | p |
| --- | --- | ---: | ---: | ---: |
| CLDN4 vs CytoTRACE2 | GSE123902 | 13 | -0.467 | 0.1076 |
| CLDN4 vs CytoTRACE2 | GSE131907 | 21 | -0.512 | 0.0177 |
| CLDN4 vs CytoTRACE2 | GSE205335 | 22 | -0.248 | 0.2660 |
| CLDN4 vs CytoTRACE2 | GSE189357 | 9 | 0.367 | 0.3317 |
| CLDN4 vs Ben-Porath ES1 stemness | GSE123902 | 13 | 0.445 | 0.1275 |
| CLDN4 vs Ben-Porath ES1 stemness | GSE131907 | 21 | 0.106 | 0.6459 |
| CLDN4 vs Ben-Porath ES1 stemness | GSE205335 | 22 | 0.295 | 0.1821 |
| CLDN4 vs Ben-Porath ES1 stemness | GSE189357 | 9 | 0.067 | 0.8647 |
| CLDN4 vs Wong ESC stemness | GSE123902 | 13 | 0.099 | 0.7479 |
| CLDN4 vs Wong ESC stemness | GSE131907 | 21 | 0.058 | 0.8013 |
| CLDN4 vs Wong ESC stemness | GSE205335 | 22 | 0.074 | 0.7436 |
| CLDN4 vs Wong ESC stemness | GSE189357 | 9 | 0.033 | 0.9322 |
| CLDN4 vs barrier/keratin (no CLDN4) | GSE123902 | 13 | 0.571 | 0.0413 |
| CLDN4 vs barrier/keratin (no CLDN4) | GSE131907 | 21 | 0.116 | 0.6178 |
| CLDN4 vs barrier/keratin (no CLDN4) | GSE205335 | 22 | 0.392 | 0.0708 |
| CLDN4 vs barrier/keratin (no CLDN4) | GSE189357 | 9 | 0.450 | 0.2242 |
| CytoTRACE2 vs Ben-Porath ES1 | GSE123902 | 13 | -0.027 | 0.9290 |
| CytoTRACE2 vs Ben-Porath ES1 | GSE131907 | 21 | 0.365 | 0.1038 |
| CytoTRACE2 vs Ben-Porath ES1 | GSE205335 | 22 | 0.202 | 0.3683 |
| CytoTRACE2 vs Ben-Porath ES1 | GSE189357 | 9 | 0.850 | 0.0037 |
| CytoTRACE2 vs barrier/keratin (no CLDN4) | GSE123902 | 13 | -0.390 | 0.1876 |
| CytoTRACE2 vs barrier/keratin (no CLDN4) | GSE131907 | 21 | -0.073 | 0.7541 |
| CytoTRACE2 vs barrier/keratin (no CLDN4) | GSE205335 | 22 | -0.669 | 6.61e-04 |
| CytoTRACE2 vs barrier/keratin (no CLDN4) | GSE189357 | 9 | 0.483 | 0.1875 |
| CLDN4 vs frac Differentiated | GSE123902 | 13 | 0.511 | 0.0743 |
| CLDN4 vs frac Differentiated | GSE131907 | 21 | 0.525 | 0.0144 |
| CLDN4 vs frac Differentiated | GSE205335 | 22 | 0.255 | 0.2512 |
| CLDN4 vs frac Differentiated | GSE189357 | 9 | 0.167 | 0.6682 |

## Sensitivity (not in the BH family)

| Contrast | k or n | ρ or Δ | p | note |
| --- | ---: | ---: | ---: | --- |
| pooled Spearman CLDN4 vs CytoTRACE2 (not meta) | 65 | -0.255 | 0.0408 | all eligible units stacked; cohort is not the weight |
| DL CLDN4 %pos vs CytoTRACE2 | 65 | -0.041 | 0.8705 | sensitivity; %pos is the T/NK score |
| DL CLDN4 vs Gulati2020 | 65 | 0.232 | 0.0855 | gene-count CytoTRACE, per dataset |
| DL CLDN4 vs preKNN CytoTRACE2 | 65 | -0.296 | 0.0780 | unsmoothed CytoTRACE2 |
| DL CLDN4 vs Ben-Porath core nine | 65 | 0.167 | 0.2209 | n=9 TF set |
| DL CytoTRACE2 vs Wong ESC | 65 | 0.495 | 7.65e-05 | second stemness agreement |
| LOO drop GSE123902: CLDN4 vs CytoTRACE2 | 52 | -0.233 | 0.3013 | leave-one-cohort-out |
| LOO drop GSE131907: CLDN4 vs CytoTRACE2 | 44 | -0.191 | 0.3742 | leave-one-cohort-out |
| LOO drop GSE205335: CLDN4 vs CytoTRACE2 | 43 | -0.298 | 0.2455 | leave-one-cohort-out |
| LOO drop GSE189357: CLDN4 vs CytoTRACE2 | 56 | -0.402 | 0.0035 | leave-one-cohort-out |
| GSE205335 ADC-only CLDN4 vs CytoTRACE2 | 14 | -0.323 | 0.2599 | not meta |
| GSE123902 primary-only CLDN4 vs CytoTRACE2 | 8 | -0.310 | 0.4556 | donor tissue; n may be small |
| GSE123902 met-only CLDN4 vs CytoTRACE2 | 5 | -0.100 | 0.8729 | donor tissue; n may be small |

## Paired within-unit CLDN4-high minus CLDN4-low

Negative potency or stemness Δ = high arm more differentiated. Positive barrier Δ = high arm more barrier-like.

| Contrast | n_units | Δ median | p |
| --- | ---: | ---: | ---: |
| CytoTRACE2 high vs low | 65 | -0.001 | 0.5384 |
| Ben-Porath ES1 high vs low | 65 | 0.072 | 8.17e-08 |
| Wong ESC high vs low | 65 | 0.069 | 7.21e-06 |
| barrier/keratin high vs low | 65 | 0.324 | 3.53e-12 |
| Gulati2020 high vs low | 65 | 0.179 | 6.14e-11 |
| preKNN CytoTRACE2 high vs low | 65 | -0.003 | 0.3526 |

## What this does not claim

- Not a rewrite of concordant-4 T/NK exclusion (malignant CLDN4 %pos vs T/NK, n=65).
- Not an ICI / RECIST / stage test. GSE189357 stage n=3 per bin is metadata only.
- Marker-malignant is not copy-number malignant.
- GSE205335 mixes ADC / SQ / SCLC / NUT. ADC-only is in the sensitivity table.
- A tiny paired Δ is not a large within-tumor stemness collapse.
- Do not call CLDN4 a stemness marker from this table.
- Visium same-spot correlation is not in this folder.

## Outputs

- `results/tables/patient_cldn4_vs_potency.tsv` — done criterion
- `results/tables/extract_audit.tsv` — marker-gate counts vs the locked pseudobulk n
- `results/tables/patient_means.tsv`
- `results/tables/cohort_spearman.tsv`
- `results/tables/stats.tsv`
- `results/figures/fig_patient_cldn4_vs_potency.png`
- `results/figures/fig_patient_cldn4_vs_stemness.png`
- `results/figures/fig_patient_cldn4_vs_barrier.png`
- `results/figures/fig_forest_cohort_rho.png`
- `results/figures/fig_paired_high_vs_low.png`
- `results/figures/fig_honest_n.png`

## Reproduce

```bash
python3 -m venv /tmp/ct2venv
/tmp/ct2venv/bin/pip install -r methods/concordant4_cytotrace2_cldn4/requirements.txt
/tmp/ct2venv/bin/python methods/concordant4_cytotrace2_cldn4/scripts/download.py --out /tmp/concordant4_ct2
/tmp/ct2venv/bin/python methods/concordant4_cytotrace2_cldn4/scripts/extract_malignant.py --data /tmp/concordant4_ct2
/tmp/ct2venv/bin/python methods/concordant4_cytotrace2_cldn4/scripts/analyze.py \
  --h5ad /tmp/concordant4_ct2/h5ad --outdir methods/concordant4_cytotrace2_cldn4/results
```

