# FINDING — QUAD CLDN4-only AUCell (IFN / MHC / TJ regulons)

**Verdict: `PARTIAL`.** Patient-level CLDN4-high vs low on malignant cells from **GSE207422 + GSE131907 + GSE148071 + GSE205335** (207422 included). A10 ELF3–CLDN4 is taken as given and is not a success criterion. No TACSTD2 dual-high split.

Pooled n = **92 patients** (within-cohort z-score of patient means; min 20 malignant cells). **TJ regulons are higher in CLDN4-high** (TJ_prior_union Δ=+0.43, p=0.033; TJ_STRUCT Δ=+0.47, p=0.013; GRHL2_prior Δ=+0.57, p=0.001). **IFN and MHC are not supported** (ISG_CORE p=0.36; MHC1_APM p=0.35; MHC Spearman ρ≈−0.19, p≈0.07). Controls (OXPHOS/ribo) did not move.

Primary unit is the **patient**. Cell counts below are inventory only.

Patients with malignant cells but **&lt;20 cells were dropped** (honest): GSE207422 2/15; GSE131907 1/26; GSE148071 5/42; GSE205335 0/17. NLRC5 has **no public prior edges** in TRRUST/DoRothEA/CollecTRI (n_targets=0).

## Pre-specified rules

- Split: within-cohort median of patient-mean malignant CLDN4 (high ≥ median).
- Pool: z-score patient means within cohort, then Wilcoxon high vs low and Spearman CLDN4 vs regulon.
- Expected: TJ up in CLDN4-high; IFN/MHC down; OXPHOS/ribo controls null.
- SUPPORT on an axis: expected direction and p<0.05 on a primary regulon (axis prior-union or program set).
- Overall SUPPORT: (IFN or MHC) **and** TJ. PARTIAL: one axis. NOT_SUPPORTED: none. UNDERPOWERED: n<12.
- ELF3_prior is A10-given and excluded from the verdict.
- AUCell = Aibar recovery curve on public priors / hardcoded programs. Not cisTarget. Not ChIP.

## Honest n

| cohort | patients (≥20 mal. cells) | malignant cells | CLDN4-high | CLDN4-low | malignant definition |
| --- | ---: | ---: | ---: | ---: | --- |
| GSE207422 | 13 | 10846 | 7 | 6 | marker epithelial AND normal-lung score ≤ p75 (CopyKAT not on GEO) |
| GSE131907 | 25 | 32749 | 13 | 12 | author Epithelial cells in tumor/met tissue (tLung/tL/B/mLN/mBrain/PE); no CopyKAT |
| GSE148071 | 37 | 48082 | 19 | 18 | TISCH2 major-lineage = Malignant (Wu et al. 2021 via TISCH2) |
| GSE205335 | 17 | 15635 | 9 | 8 | author lineage.sub = Malignant cells; ADC+SQ; drop normal tissues |
| **QUAD pool** | **92** | **107312** | 48 | 44 | within-cohort z |

## Primary regulon table (pooled, within-cohort z)

| regulon | axis | expected | n_high | n_low | Δ (high−low) | Wilcoxon p | Spearman ρ (p) | call |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| IFN_prior_union | IFN | down | 48 | 44 | 0.002 | 0.935 | -0.071 (0.504) | NS |
| ISG_CORE | IFN | down | 48 | 44 | -0.197 | 0.363 | -0.103 (0.329) | NS |
| IRF1_prior | IFN | down | 48 | 44 | 0.217 | 0.396 | 0.042 (0.688) | NS |
| STAT1_prior | IFN | down | 48 | 44 | 0.058 | 0.699 | -0.035 (0.739) | NS |
| MHC_prior_union | MHC | down | 48 | 44 | -0.124 | 0.342 | -0.190 (0.070) | NS |
| MHC1_APM | MHC | down | 48 | 44 | -0.182 | 0.350 | -0.192 (0.067) | NS |
| NLRC5_prior | MHC | down | 0 | 0 | NA | NA | NA (NA) | UNDERPOWERED |
| CIITA_prior | MHC | down | 48 | 44 | -0.094 | 0.401 | -0.189 (0.071) | NS |
| TJ_prior_union | TJ | up | 48 | 44 | 0.431 | 0.033 | 0.316 (0.002) | SUPPORT |
| TJ_STRUCT | TJ | up | 48 | 44 | 0.468 | 0.013 | 0.299 (0.004) | SUPPORT |
| GRHL2_prior | TJ | up | 48 | 44 | 0.574 | 0.001 | 0.449 (7.11e-06) | SUPPORT |
| KLF4_prior | TJ | up | 48 | 44 | 0.370 | 0.054 | 0.261 (0.012) | NS |
| OVOL1_prior | TJ | up | 48 | 44 | -0.375 | 0.013 | -0.215 (0.039) | WRONG_SIGN |
| ELF3_prior | A10 | given | 48 | 44 | 0.325 | 0.043 | 0.289 (0.005) | A10_GIVEN |
| CTRL_OXPHOS | CTRL | null | 48 | 44 | -0.319 | 0.113 | -0.181 (0.085) | NULL_OK |
| CTRL_RIBO | CTRL | null | 48 | 44 | -0.191 | 0.193 | -0.095 (0.366) | NULL_OK |

Full tests: `tables/tests.tsv`. Regulon gene lists: `tables/regulons.tsv`. Patient scores: `tables/patient_scores.tsv`.

## What this is / is not

- **Is** additive CLDN4-only AUCell on a four-cohort malignant merge, patient-level, 207422 included.
- **Is not** a Harmony/scVI joint embedding (patients are z-scored within cohort; no cell mixing).
- **Is not** pySCENIC cisTarget or a binding map. Priors are TRRUST/DoRothEA/CollecTRI.
- **Is not** a re-test of A10 ELF3–CLDN4. ELF3_prior is reported and tagged `A10_GIVEN`.
- **Is not** a TACSTD2 / dual-high analysis.
- GSE131907 malignant is tumor-tissue epithelium (author label), not CopyKAT.
- GSE148071 expression is TISCH2 `log2(TPM/10+1)`; others are `log1p(CP10k)` from UMI. Ranking AUCell is within-matrix; pooling uses within-cohort z.
- GSE205335 SCLC/NUT and normal tissues were excluded from the primary NSCLC malignant set.

## Regulon inventory

| regulon | axis | kind | n_targets | expected |
| --- | --- | --- | ---: | --- |
| IRF1_prior | IFN | public_prior | 146 | down |
| IRF7_prior | IFN | public_prior | 20 | down |
| IRF9_prior | IFN | public_prior | 30 | down |
| STAT1_prior | IFN | public_prior | 806 | down |
| STAT2_prior | IFN | public_prior | 72 | down |
| NLRC5_prior | MHC | public_prior | 0 | down |
| CIITA_prior | MHC | public_prior | 28 | down |
| RFX5_prior | MHC | public_prior | 31 | down |
| ELF3_prior | A10 | a10_given_prior | 26 | given |
| GRHL1_prior | TJ | public_prior | 1 | up |
| GRHL2_prior | TJ | public_prior | 10 | up |
| KLF4_prior | TJ | public_prior | 110 | up |
| OVOL1_prior | TJ | public_prior | 9 | up |
| OVOL2_prior | TJ | public_prior | 1 | up |
| IFN_prior_union | IFN | combinatorial_public_prior | 921 | down |
| MHC_prior_union | MHC | combinatorial_public_prior | 45 | down |
| TJ_prior_union | TJ | combinatorial_public_prior | 127 | up |
| ISG_CORE | IFN | program_ISG | 40 | down |
| MHC1_APM | MHC | program_APM | 21 | down |
| TJ_STRUCT | TJ | program_TJ | 27 | up |
| CTRL_OXPHOS | CTRL | program_control | 15 | null |
| CTRL_RIBO | CTRL | program_control | 16 | null |

