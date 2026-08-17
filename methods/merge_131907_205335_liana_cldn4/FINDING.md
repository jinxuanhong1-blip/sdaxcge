# FINDING — merged GSE131907 + GSE205335 LIANA/LR from CLDN4-high malignant to T/NK

**Verdict:** On 40 paired patients (19 GSE131907 + 21 GSE205335), CLDN4-high vs CLDN4-low malignant → same-patient T/NK focus pairs do **not** support a coordinated T-recruit / MHC-I drop (5/22 median Δ < 0; **0 FDR < 0.05** weaker from high). The LR table is a ranked co-expression list, not a causal claim.

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. No GSE207422. Patient is the unit. p-values are descriptive.

The between-patient abundance cut on this same slice is **given** and is not re-audited: PR #320 Q4 vs Q1 author %pos vs T/NK **n=23 (12/11) r=-0.705**. This folder tests a different contrast: outgoing CellPhoneDB-style scores from CLDN4-high vs CLDN4-low malignant cells to **same-patient** T/NK.

---

## Honest n

| Item | n | Note |
| --- | ---: | --- |
| GSE131907 cells in UMI matrix | 208,506 | author barcodes |
| GSE205335 cells in UMI matrix | 96,505 | author barcodes |
| GSE131907 locked samples (PR #320) | 21 | tumor-origin, n_malignant ≥20 |
| GSE131907 unique patients in locked samples | 21 | GEO patient_id; extract is sample-level |
| GSE205335 locked patients (PR #320) | 22 | ≥20 malignant and ≥20 T/NK |
| Author malignant cells (locked slice) | 53,296 | 24,784 + 28,512 |
| CLDN4-high / low malignant | 26,648 / 26,648 | cohort-specific global medians |
| T/NK cells (locked slice) | 48,150 | 15,150 + 33,000 |
| **Paired patients (Wilcoxon unit)** | **40** | ≥10 high, ≥10 low, ≥20 T/NK |
|  … GSE131907 / GSE205335 | 19 / 21 | GSE131907 cells pooled by patient_id |
| Pairs scored / ranked | 2873 / 109 | CellPhoneDB v5 + overlay |
| LIANA | ran_cellphonedb_method | secondary; not the patient test |
| CellChat | not_run_R_unavailable | not run |

GSE131907 CLDN4 median log1p(CP10k) among locked malignant cells = 1.819. GSE205335 median = 1.371. Thresholds are not shared across platforms.

Per-patient counts: `results/n_cells_patients.tsv`.

### Locked patients that failed the paired LR gate

| cohort | patient | n_high | n_low | n_T/NK | reason |
|---|---|---:|---:|---:|---|
| GSE131907 | P3016 | 0 | 79 | 385 | high<10 |
| GSE131907 | P1013 | 0 | 376 | 1150 | high<10 |
| GSE205335 | P4001 | 6 | 21 | 1801 | high<10 |

## Method

Documented CellPhoneDB-style score on log1p(CP10k): partner expression = **min of subunit means**; pair score = **mean of the two partner means** (Efremova et al. 2020; Garcia-Alonso et al. 2022). `pass_expr_prop` requires both partners in ≥10% of cells in their group. The test is a **paired Wilcoxon** of high vs low **per patient**. FDR is Benjamini–Hochberg within the outgoing contrast. This is **not** a CellChat communication probability and does **not** observe secretion or spatial contact.

LIANA `mt.cellphonedb` was run. If present, LIANA p-values are within-object specificity, not patient-level tests.

## Primary LR table — outgoing CLDN4-high malignant → T/NK

Median Δ = high − low. Negative = weaker from CLDN4-high. Paired n = **40** (19 GSE131907 + 21 GSE205335).

### MHC-I and T-recruit focus pairs

| pathway | ligand | receptor | n_patients | n_GSE131907 | n_GSE205335 | median_delta | pval | padj |
|---|---|---|---|---|---|---|---|---|
| MHC_I | HLA-F | KIR3DL2 | 3 | 1 | 2 | +0.026 | 0.25 | 0.296 |
| MHC_I | HLA-C | CD8A | 35 | 17 | 18 | +0.112 | 9.99e-05 | 0.000419 |
| MHC_I | HLA-B | CD8A | 35 | 17 | 18 | +0.118 | 0.000262 | 0.000894 |
| MHC_I | HLA-A | CD8A | 35 | 17 | 18 | +0.122 | 8.43e-05 | 0.000383 |
| MHC_I | HLA-A | CD8B | 31 | 14 | 17 | +0.132 | 0.000453 | 0.00145 |
| MHC_I | HLA-E | KLRC1 | 13 | 8 | 5 | +0.137 | 0.000244 | 0.000858 |
| MHC_I | HLA-C | CD8B | 31 | 14 | 17 | +0.137 | 0.000495 | 0.0015 |
| MHC_I | HLA-E | KLRC1+KLRD1 | 12 | 7 | 5 | +0.138 | 0.000488 | 0.0015 |
| MHC_I | HLA-E | KLRD1 | 34 | 16 | 18 | +0.139 | 5.79e-07 | 7.01e-06 |
| MHC_I | HLA-B | CD8B | 31 | 14 | 17 | +0.151 | 0.00107 | 0.00298 |
| MHC_I | HLA-E | KLRC3+KLRD1 | 6 | 0 | 6 | +0.159 | 0.0625 | 0.0885 |
| MHC_I | HLA-E | KLRC2 | 15 | 4 | 11 | +0.174 | 0.000122 | 0.000475 |
| MHC_I | HLA-E | KLRC2+KLRD1 | 15 | 4 | 11 | +0.174 | 0.000122 | 0.000475 |
| MHC_I | HLA-E | KLRK1 | 19 | 0 | 19 | +0.174 | 0.0141 | 0.0264 |
| MHC_I | HLA-B | KIR3DL2 | 3 | 1 | 2 | +0.184 | 0.25 | 0.296 |
| T_recruit | CCL5 | CCR5 | 5 | 2 | 3 | -0.084 | 0.125 | 0.166 |
| T_recruit | CCL5 | CCR1 | 4 | 2 | 2 | -0.050 | 0.25 | 0.296 |
| T_recruit | CCL4 | CCR5 | 4 | 3 | 1 | -0.030 | 0.125 | 0.166 |
| T_recruit | CXCL10 | CXCR3 | 6 | 4 | 2 | -0.012 | 0.562 | 0.607 |
| T_recruit | CCL3 | CCR5 | 4 | 3 | 1 | -0.005 | 0.625 | 0.661 |
| T_recruit | CXCL16 | CXCR6 | 21 | 10 | 11 | +0.031 | 0.00118 | 0.0032 |
| T_recruit | CX3CL1 | CX3CR1 | 7 | 4 | 3 | +0.036 | 0.0156 | 0.0279 |

Focus pairs with median Δ < 0: **5/22**. FDR < 0.05 and Δ < 0: **0/22**. FDR < 0.05 in either direction: **14/22**.

### Top paired hits (all pathways, by Wilcoxon p)

| pathway | ligand | receptor | n_patients | n_GSE131907 | n_GSE205335 | median_delta | pval | padj |
|---|---|---|---|---|---|---|---|---|
| other | CDH1 | ITGAE+ITGB7 | 32 | 15 | 17 | +0.152 | 3.26e-09 | 3.55e-07 |
| other | CD55 | ADGRE5 | 37 | 18 | 19 | +0.125 | 1.11e-08 | 5.41e-07 |
| other | F11R | ITGAL+ITGB2 | 28 | 12 | 16 | +0.147 | 1.49e-08 | 5.41e-07 |
| other | ICAM1 | SPN | 32 | 17 | 15 | +0.105 | 7.87e-08 | 1.62e-06 |
| other | APP | CD74 | 38 | 18 | 20 | +0.063 | 8.27e-08 | 1.62e-06 |
| other | CDH1 | KLRG1 | 26 | 10 | 16 | +0.155 | 8.94e-08 | 1.62e-06 |
| other | APP | SORL1 | 27 | 11 | 16 | +0.082 | 1.49e-07 | 2.32e-06 |
| other | ALCAM | CD6 | 35 | 17 | 18 | +0.044 | 5.78e-07 | 7.01e-06 |
| MHC_I | HLA-E | KLRD1 | 34 | 16 | 18 | +0.139 | 5.79e-07 | 7.01e-06 |
| other | PTGES | PTGER4 | 25 | 12 | 13 | +0.023 | 1.13e-06 | 1.23e-05 |
| other | ICAM1 | ITGAL | 25 | 12 | 13 | +0.121 | 3.28e-06 | 2.86e-05 |
| other | ICAM1 | ITGAL+ITGB2 | 25 | 12 | 13 | +0.121 | 3.28e-06 | 2.86e-05 |
| other | ITGAV+ITGB1 | ADGRE5 | 38 | 18 | 20 | +0.030 | 3.41e-06 | 2.86e-05 |
| other | IGFBP3 | TMEM219 | 27 | 14 | 13 | +0.097 | 4.57e-06 | 3.56e-05 |
| other | PLAUR | ITGA4+ITGB1 | 26 | 9 | 17 | +0.089 | 5.04e-06 | 3.66e-05 |

Full ranked table: `results/lr_table.tsv`. Per-patient scores: `results/patient_outgoing.tsv.gz`.

## LIANA CellPhoneDB method (secondary, pooled / downsampled)

LIANA `mt.cellphonedb` ran separately per cohort (`expr_prop=0.10`, 50 permutations, ≤2,000 cells/group). These p-values are within-object specificity, not the patient Wilcoxon.

| Cohort | groups | outgoing CLDN4-high → T/NK edges |
|---|---|---:|
| GSE131907 | high/low/T/NK = 2000 each | 50 |
| GSE205335 | high/low/T = 2000 each (author T/NK is one label) | 34 |

Focus-ish edges that cleared LIANA `expr_prop` from CLDN4-high malignant to T or NK:

| Cohort | Pair | target | `lr_means` | cellphone_p |
|---|---|---|---:|---:|
| GSE131907 | CXCL16–CXCR6 | T | 0.270 | 0 |
| GSE131907 | CX3CL1–CX3CR1 | NK | 0.356 | 0 |
| GSE131907 | CCL20–CXCR3 | T | 0.400 | 0 |
| GSE205335 | CXCL16–CXCR6 | T | 0.205 | 0 |
| GSE205335 | CCL20–CXCR3 | T | 0.293 | 0 |
| GSE205335 | HLA-E–KLRK1 | T | 0.623 | 1 |

CXCL9/10–CXCR3 typically fail the 10% filter. LIANA does not replace the paired n=40 table.

## What this is not

- Not a re-audit of PR #320 Q4 vs Q1 n=23 r=−0.705.
- Not dual-high TACSTD2×CLDN4.
- Not GSE207422.
- Not inferCNV/CopyKAT recomputed malignant IDs.
- Not CellChat. LIANA permutation p-values (if present) are not patient tests.
- Not spatial proximity or protein secretion.
- Cells are not the sample size. GSE131907 PR #320 extract is sample-level; this Wilcoxon uses unique patients.

## Extra figures

- `figures/n_cells_by_patient.png`
- `figures/cldn4_outgoing_focus.png`
- `figures/top_hits_delta.png`
- `figures/focus_by_cohort.png`

## Reproduce

```bash
python3 methods/merge_131907_205335_liana_cldn4/scripts/00_download.py
python3 methods/merge_131907_205335_liana_cldn4/scripts/01_analyze.py
```
