# FINDING — merged GSE131907+GSE205335 NicheNet-style CLDN4-only ligand activity

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. **GSE207422 is not re-run**
(PR #334: n=12 signed tests were NS). Sender = CLDN4-high vs low **malignant**.
Receiver = **same-patient T/NK**. Patient is the unit. GSE131907 GEO extract is
sample-level; cells from eligible samples are collapsed to unique `patient_id`.

Prior = published **NicheNet-v2** ligand–target matrix (Browaeys et al.;
Zenodo 10.5281/zenodo.7074291, `ligand_target_matrix_nsga2r_final.rds`).
Converted and scored in **Python**. R / `nichenetr` was **not** installed and was **not** run.

## Honest n

| Item | n | Note |
| --- | ---: | --- |
| GSE131907 eligible samples | 31 | tumor-origin, ≥20 author-malig, ≥20 T/NK |
| GSE131907 unique patients | **31** | **unit** after collapse |
| GSE205335 patients | **22** | locked ≥20 author-malig + T/NK (PR #320) |
| Merged patients | **53** | two cohorts; not one mixed bag of cells |
| Malignant cells (eligible) | 59,643 | author labels |
| CLDN4-high malignant | 29,823 | log1p(CP10k) ≥ cohort malignant median |
| Dual-high companion (not used) | 23,367 | TACSTD2 and CLDN4 ≥ median |
| T/NK cells (eligible) | 92,147 | same-patient receivers |
| Potential ligands (merged) | **136** | ≥10% CLDN4-high + T/NK receptor in v2 LR |
| Background genes | 264 | T/NK-expressed ∩ prior targets |
| IFN / cytotoxicity genes in prior | 16 / 14 | a priori lists |

Cells are counts, not replicates. p-values are descriptive.

## Signed patient-level T/NK programs (not GSE207422)

CLDN4 score = malignant **%pos** (same primary cut as PR #320). Programs = mean
log1p(CP10k) of a priori IFN / cytotoxicity lists in same-patient T/NK.

| test | n | result |
| --- | ---: | --- |
| Spearman CLDN4 %pos vs T/NK IFN (DL merge) | 53 | ρ=-0.243 p=0.0897 I²=0% |
| Spearman CLDN4 %pos vs T/NK cytotoxicity (DL merge) | 53 | ρ=-0.161 p=0.449 I²=52% |
| Spearman CLDN4 %pos vs T/NK fraction (DL merge) | 53 | ρ=-0.362 p=0.00928 I²=0% |
| Q4 vs Q1 T/NK IFN (pooled r, tails only) | 28 | ρ=-0.482 p=0.0138 I²=0% |
| Q4 vs Q1 T/NK cytotoxicity (pooled r, tails only) | 28 | ρ=-0.228 p=0.276 I²=0% |

Within-cohort rows are in `results/patient_level_tests.tsv`.
GSE205335 Q4 mixes SCLC with ADC (PR #362); that mix is reported, not hidden.

## Ligand activity (unsigned NicheNet-v2 prior recovery)

Primary sender = pooled CLDN4-high malignant from eligible patients in both
cohorts. Receiver background = T/NK-expressed genes present in the prior.
Rank = Pearson of the ligand’s prior target scores vs gene-set membership.
These ranks recover **prior structure** (ISG / MHC ligands sit next to IFN
genes). They are **not** by themselves a signed CLDN4→IFN inductive axis.

**A priori IFN (16 genes in prior). Top 8 of 136:**

| Rank | Ligand | Pearson | AUROC | pearson p |
| --- | --- | ---: | ---: | ---: |
| 1 | **IFITM1** | 0.515 | 0.934 | 2.57e-19 |
| 2 | **LIF** | 0.457 | 0.895 | 5.10e-15 |
| 3 | **CLCF1** | 0.416 | 0.892 | 1.77e-12 |
| 4 | **HLA-F** | 0.416 | 0.884 | 1.81e-12 |
| 5 | **CRLF1** | 0.410 | 0.884 | 3.87e-12 |
| 6 | **VSIG10** | 0.409 | 0.885 | 4.47e-12 |
| 7 | **HLA-DRB5** | 0.386 | 0.874 | 8.38e-11 |
| 8 | **HLA-B** | 0.375 | 0.867 | 3.03e-10 |

**A priori cytotoxicity (14 genes). Top 8:**

| Rank | Ligand | Pearson | AUROC | pearson p |
| --- | --- | ---: | ---: | ---: |
| 1 | **LRPAP1** | 0.321 | 0.353 | 9.94e-08 |
| 2 | **HLA-DRA** | 0.286 | 0.400 | 2.26e-06 |
| 3 | **HLA-DQB1** | 0.251 | 0.406 | 3.63e-05 |
| 4 | **SPP1** | 0.189 | 0.344 | 0.00209 |
| 5 | **TNC** | 0.182 | 0.347 | 0.00306 |
| 6 | **HLA-A** | 0.181 | 0.471 | 0.00316 |
| 7 | **HLA-DRB1** | 0.162 | 0.402 | 0.00836 |
| 8 | **ITGB2** | 0.152 | 0.359 | 0.0135 |

## Signed ligand vs same-patient T/NK program (patient-level)

Per-patient mean ligand log1p(CP10k) in **CLDN4-high malignant** cells vs
same-patient T/NK IFN or cytotoxicity. Within-cohort Spearman, then
DerSimonian–Laird on Fisher-z. This is the signed, patient-unit test.

**vs T/NK IFN (lowest p, unadjusted):**

| Ligand | k | N | pooled ρ | p | I² |
| --- | ---: | ---: | ---: | ---: | ---: |
| LAMA5 | 2 | 51 | +0.435 | 0.00179 | 0% |
| SDC2 | 2 | 51 | -0.427 | 0.00221 | 0% |
| CD99 | 2 | 51 | +0.409 | 0.00361 | 0% |
| S100A8 | 2 | 51 | +0.405 | 0.00398 | 0% |
| NECTIN4 | 1 | 22 | +0.523 | 0.0113 | 0% |
| DUSP18 | 2 | 51 | +0.341 | 0.017 | 0% |
| MIF | 2 | 51 | +0.430 | 0.0182 | 40% |
| GPI | 2 | 51 | +0.329 | 0.0218 | 0% |

**vs T/NK cytotoxicity (lowest p, unadjusted):**

| Ligand | k | N | pooled ρ | p | I² |
| --- | ---: | ---: | ---: | ---: | ---: |
| NUCB2 | 2 | 51 | -0.425 | 0.0023 | 0% |
| MMP14 | 2 | 51 | +0.425 | 0.00234 | 0% |
| TGFA | 2 | 51 | +0.417 | 0.0029 | 0% |
| CLCF1 | 2 | 51 | +0.409 | 0.00355 | 0% |
| HLA-E | 2 | 51 | +0.400 | 0.00449 | 0% |
| HLA-F | 2 | 51 | +0.393 | 0.00531 | 0% |
| JAG1 | 2 | 51 | +0.393 | 0.00535 | 0% |
| ANGPTL4 | 2 | 51 | +0.371 | 0.00889 | 0% |

## What this is not

- Not GSE207422 NicheNet (PR #286 dual-high; PR #334 CLDN4-only n=12 NS).
- Not a full-transcriptome `nichenetr` R run. Background is T/NK-expressed genes ∩ prior.
- Not TACSTD2∩CLDN4 dual-high senders. Companion dual-high count is reported only.
- Not CellChat / LIANA (those are other PRs on the single cohorts).
- Cells are not n. GSE131907 sample n and unique-patient n are both stated.

## Files

- `results/ligand_activity_table.tsv` — primary table (prior Pearson/AUROC + setting)
- `results/ligand_activity_all.tsv` — all settings × gene sets
- `results/ligand_vs_tnk_program_merged.tsv` — signed patient-level ligand tests
- `results/patient_scores.tsv`, `n_table.tsv`, `patient_level_tests.tsv`
- Extra figures: `figures/`

## Reproduce

```bash
python3 methods/merge_131907_205335_nichenet_cldn4/scripts/00_download.py
python3 methods/merge_131907_205335_nichenet_cldn4/scripts/01_convert_prior.py
python3 methods/merge_131907_205335_nichenet_cldn4/scripts/02_extract.py
python3 methods/merge_131907_205335_nichenet_cldn4/scripts/03_analyze.py
```

