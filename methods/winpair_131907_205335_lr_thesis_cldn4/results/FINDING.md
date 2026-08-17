# FINDING — winning-pair thesis LR (CLDN4-only, both families)

**Verdict:** On 40 paired patients (19 GSE131907 + 21 GSE205335), the barrier/inhibitory family is **8/8 pairs in the thesis direction** (high > low; composite Δ=+0.082, p=1.82e-11; 8 pairs p<0.05 and Δ>0). The IFN/T-recruit/MHC-I family is **2/10 pairs in the thesis direction** (low > high; composite Δ=+0.069, p=4.73e-05; 1 pair p<0.05 and Δ<0). Barrier-up-in-high is a primary arm, not a recruit leftover. Scores are CellPhoneDB-style co-expression, not secretion or contact.

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. No GSE148071. Patient is the unit. p-values are descriptive.

The between-patient abundance cut on this same slice is **given** and is not re-audited: PR #320 Q4 vs Q1 author %pos vs T/NK **n=23 (12/11) r=-0.705**. This folder tests a different contrast: outgoing CellPhoneDB-style scores from CLDN4-high vs CLDN4-low malignant cells to **same-patient** T/NK, on two **pre-specified** pair families.

---

## Thesis (pre-specified; both arms)

1. **CLDN4-high** malignant cells send **more barrier / inhibitory** pairs to T/NK: F11R–ITGAL/ITGB2, NECTIN2–TIGIT, CDH1–ITGAE, LGALS9–CD45/CD44.
2. **CLDN4-low (KD-like)** malignant cells send **more IFN / T-recruit / MHC-I** pairs: CXCL9/10/11–CXCR3, CCL5–CCR5, HLA-A/B/C–CD8.

Both families are scored and tabulated. Barrier-up-in-high hits are not filed as "recruit-up, skip."

## Honest n

| Item | n | Note |
| --- | ---: | --- |
| GSE131907 cells in UMI matrix | 208,506 | author barcodes |
| GSE205335 cells in UMI matrix | 96,505 | author barcodes |
| GSE131907 locked samples (PR #320) | 21 | tumor-origin, n_malignant ≥20 |
| GSE131907 unique patients in locked samples | 21 | GEO patient_id; extract is sample-level |
| GSE205335 locked patients (PR #320) | 22 | tumor extract |
| Author malignant cells (locked slice) | 53,296 | 24,784 + 28,512 |
| CLDN4-high / low malignant | 26,648 / 26,648 | cohort-specific global medians |
| T/NK cells (locked slice) | 48,150 | 15,150 + 33,000 |
| **Paired patients (Wilcoxon unit)** | **40** | ≥10 high, ≥10 low, ≥20 T/NK |
|  … GSE131907 / GSE205335 | 19 / 21 | GSE131907 cells pooled by patient_id |
| Pre-specified pairs | 18 | 8 barrier/inhibitory + 10 IFN/recruit/MHC-I |

GSE131907 CLDN4 median log1p(CP10k) among locked malignant cells = 1.819. GSE205335 median = 1.371. Thresholds are not shared across platforms.

Per-patient counts: `results/n_cells_patients.tsv`.

### Locked patients that failed the paired LR gate

| cohort | patient | n_high | n_low | n_T/NK | reason |
|---|---|---:|---:|---:|---|
| GSE131907 | P3016 | 0 | 79 | 385 | high<10 |
| GSE131907 | P1013 | 0 | 376 | 1150 | high<10 |
| GSE205335 | P4001 | 6 | 21 | 1801 | high<10 |

## Method

Documented CellPhoneDB-style score on log1p(CP10k): partner expression = **min of subunit means**; pair score = **mean of the two partner means** (Efremova et al. 2020; Garcia-Alonso et al. 2022). Every pre-specified pair is scored on **every eligible patient** (no expr_prop drop from the primary n). `n_pass_expr_prop` is the number of patients in whom both partners are in ≥10% of cells in high or low; it is a sparsity flag, not a license to skip the pair. The test is a **paired Wilcoxon** of high vs low **per patient**. Δ = median(high − low). Positive Δ = stronger from CLDN4-high. Family FDR is Benjamini–Hochberg **within family**. Family composite = mean of that family's pair scores per patient, then the same Wilcoxon. This is **not** a CellChat communication probability and does **not** observe secretion or spatial contact. GSE131907 `PVRL2` is aliased to `NECTIN2`. Because the T/NK receiver is the same patient pool for high and low senders, Δ = ½(ligand_high − ligand_low). Pairs that share a ligand (F11R–ITGAL vs F11R–ITGB2; HLA-A–CD8A vs HLA-A–CD8B) therefore share Δ. That is the tumor-cell program, not a copy error.

## Primary table 1 — barrier / inhibitory (thesis: Δ > 0)

Family composite: n=40, Δ=+0.082, p=1.82e-11 (one-sided thesis p=9.09e-12).

Pairs in the thesis direction (Δ > 0) are **hits for this arm**, not a reason to skip to recruit.

| pair | n | n_GSE131907 | n_GSE205335 | delta | p | padj_family | agrees_thesis | n_pass_expr_prop | sparse_sender |
|---|---|---|---|---|---|---|---|---|---|
| F11R–ITGAL+ITGB2 | 40 | 19 | 21 | +0.100 | 9.09e-12 | 2.43e-11 | True | 28 | False |
| F11R–ITGAL | 40 | 19 | 21 | +0.100 | 9.09e-12 | 2.43e-11 | True | 28 | False |
| F11R–ITGB2 | 40 | 19 | 21 | +0.100 | 9.09e-12 | 2.43e-11 | True | 40 | False |
| NECTIN2–TIGIT | 40 | 19 | 21 | +0.083 | 7.82e-11 | 1.56e-10 | True | 35 | False |
| CDH1–ITGAE | 40 | 19 | 21 | +0.130 | 3.07e-10 | 4.1e-10 | True | 36 | False |
| CDH1–ITGAE+ITGB7 | 40 | 19 | 21 | +0.130 | 3.07e-10 | 4.1e-10 | True | 32 | False |
| LGALS9–PTPRC | 40 | 19 | 21 | +0.010 | 0.0106 | 0.0106 | True | 31 | False |
| LGALS9–CD44 | 40 | 19 | 21 | +0.010 | 0.0106 | 0.0106 | True | 32 | False |

Full table: `results/barrier_inhibitory_high.tsv`.

## Primary table 2 — IFN / T-recruit / MHC-I (thesis: Δ < 0)

Family composite: n=40, Δ=+0.069, p=4.73e-05 (one-sided thesis p=1).

A positive Δ here is a **miss for this arm** (higher from CLDN4-high, opposite the KD-like prediction). Sparse CXCR3 ligands keep their honest n; they are not dropped.

| pair | n | n_GSE131907 | n_GSE205335 | delta | p | padj_family | agrees_thesis | n_pass_expr_prop | sparse_sender |
|---|---|---|---|---|---|---|---|---|---|
| CXCL9–CXCR3 | 40 | 19 | 21 | +0.000 | 0.663 | 0.663 | False | 0 | True |
| CXCL10–CXCR3 | 40 | 19 | 21 | -0.000 | 0.446 | 0.558 | True | 6 | True |
| CXCL11–CXCR3 | 40 | 19 | 21 | +0.000 | 0.586 | 0.651 | False | 1 | True |
| CCL5–CCR5 | 40 | 19 | 21 | -0.003 | 0.0224 | 0.0319 | True | 5 | False |
| HLA-A–CD8A | 40 | 19 | 21 | +0.114 | 1.13e-05 | 3.3e-05 | False | 35 | False |
| HLA-A–CD8B | 40 | 19 | 21 | +0.114 | 1.13e-05 | 3.3e-05 | False | 31 | False |
| HLA-B–CD8A | 40 | 19 | 21 | +0.101 | 5.08e-05 | 8.46e-05 | False | 35 | False |
| HLA-B–CD8B | 40 | 19 | 21 | +0.101 | 5.08e-05 | 8.46e-05 | False | 31 | False |
| HLA-C–CD8A | 40 | 19 | 21 | +0.105 | 1.32e-05 | 3.3e-05 | False | 35 | False |
| HLA-C–CD8B | 40 | 19 | 21 | +0.105 | 1.32e-05 | 3.3e-05 | False | 31 | False |

Full table: `results/ifn_recruit_low.tsv`.

## Family composites

| family | thesis_direction | n | n_GSE131907 | n_GSE205335 | delta | p | p_onesided_thesis | agrees_thesis |
|---|---|---|---|---|---|---|---|---|
| barrier_inhibitory | high>low | 40 | 19 | 21 | +0.082 | 1.82e-11 | 9.09e-12 | True |
| ifn_recruit | low>high | 40 | 19 | 21 | +0.069 | 4.73e-05 | 1 | False |

## What this is not

- Not a re-audit of PR #320 Q4 vs Q1 n=23 r=−0.705.
- Not dual-high TACSTD2×CLDN4.
- Not GSE148071.
- Not inferCNV/CopyKAT recomputed malignant IDs.
- Not CellChat / LIANA permutation p-values. Those are not patient tests.
- Not spatial proximity or protein secretion.
- Cells are not the sample size. GSE131907 PR #320 extract is sample-level; this Wilcoxon uses unique patients.

## Extra figures

- `figures/n_cells_by_patient.png`
- `figures/barrier_inhibitory_high.png`
- `figures/ifn_recruit_low.png`
- `figures/thesis_two_arms.png`
- `figures/barrier_paired_boxes.png`
- `figures/recruit_paired_boxes.png`
- `figures/family_by_cohort.png`

## Reproduce

```bash
python3 methods/winpair_131907_205335_lr_thesis_cldn4/scripts/00_download.py
python3 methods/winpair_131907_205335_lr_thesis_cldn4/scripts/01_analyze.py
```
