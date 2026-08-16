# scRNA compositional analysis playbook (TACSTD2 / T/NK / MPR)

**Focus.** Patient-level cell-type fractions after neoadjuvant ICI, tested
against (i) malignant TACSTD2 and (ii) MPR vs NMPR. This folder is an
**analysis**, not a generic template.

中文：本目录是公开 scRNA 组成分析（患者为重复单位），不是把细胞当独立样本。

## Why not cells-as-replicates

Composition is a **sample-level** multinomial. Treating 90k cells as
replicates of an MPR label is pseudoreplication. scCODA / sccomp / milo
and CLR-ILR all operate on **samples**.

组成是样本层 multinomial。细胞不能当 MPR 的独立重复。

## Default 2024–2026 stack

| Step | Default | This run |
|---|---|---|
| Differential abundance | **scCODA** or **sccomp** | scCODA **not installable** (rpy2/R) |
| Fallback | CLR + ILR + Dirichlet-multinomial | **used** |
| Neighborhood DA | miloR | out of scope |
| Unit | patient / sample | patient |
| Parts | 5–8 lineages; collapse T+NK | Epi, T/NK, B, Myeloid, Stromal, Other |

## Hypotheses (pre-specified)

1. **H1.** T/NK fraction (and CLR/ILR T/NK) is lower when patient malignant
   TACSTD2 is high. Continuous Spearman + median-split Wilcoxon.
2. **H2.** NMPR has higher epithelial and lower T/NK than MPR (pCR = MPR).

Do not re-define “match” after seeing p-values. Report direction **and** p.

## Residual-tumor confound

Post-neoadjuvant MPR/pCR samples have few leftover epithelial cells **by
definition**. An NMPR > MPR epithelial-fraction test on resection scRNA is
partly that fact. State it every time H2 is shown.

术后 MPR 残存上皮少是定义，不是独立的“排斥”证据。

## Annotation

- Prefer **author** `major.cell.type` when GEO has it (GSE241934).
- Else a documented marker hierarchy (GSE207422, GSE291670). Hu et al.
  CopyKAT barcodes are **not** on GEO — do not invent them.
- Malignant TACSTD2 = mean `log1p(UMI)` in epithelial cells of that patient.
  Minimum 20 epithelial cells for H1.

## Combinations

Test each series, then each union. For unions, also fit

`CLR(T/NK) ~ TACSTD2 + C(cohort)`

and a cohort-centered Spearman. Do not pool IIT (EGFR-mutant, uniform
sintilimab) with REAL (EGFR-WT, mixed PD-1) without saying so.

## What this does not do

- No FASTQ, no private EGA/dbGaP objects.
- No cell-level Wilcoxon on MPR.
- No scCODA inclusion probabilities (package missing).
- No claim that a null ρ “supports” H1 because the sign is negative.
