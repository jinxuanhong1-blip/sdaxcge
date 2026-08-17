# Methods — GSE123902 + GSE205335 CLDN4-only high-end CellChat

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. Patient (GSE123902
donor) is the unit. Do not add GSE148071. Do not pile the seven-cohort pool.

## Locked combo (not re-audited)

PR #459 pair that **differs**:

| combo | score | N | ρ (p, I²) | Q4 vs Q1 r (p) |
|---|---|---:|---|---|
| GSE123902+GSE205335 | %pos | 35 | −0.522 (0.002, 0%) | −0.802 (0.005; 9 vs 9) |

| unit | malignant def | unit of analysis | n |
|---|---|---|---:|
| GSE123902 | marker-malignant (EPCAM\|KRT8\|KRT18\|KRT19 > 0 and PTPRC == 0) | donor | 13 |
| GSE205335 | author `Malignant cells` | patient | 22 |

T/NK: GSE123902 = (CD3D\|CD3E\|CD8A\|NKG7\|GNLY\|KLRD1) > 0 and not
malignant; GSE205335 = author `T/NK cells`. This folder does not re-rank
or re-pool that Spearman.

## What is new

CellChat-style outgoing communication **CLDN4-high malignant → same-patient
T/NK**, then a meta-analysis of **patient ΔP** across the two units.

High-end split (primary): within each patient, malignant cells are cut on
`log1p(CP10k)` **CLDN4** into Q4 vs Q1 (middle half dropped). Sensitivity:
median and tertile. TACSTD2 is not used.

ΔP = P_high − P_low. One delta per patient per ligand–receptor pair.
Unit means are DerSimonian–Laird random-effects pooled. A patient-level
Wilcoxon on the stacked deltas is also reported. p-values are descriptive.

Extra between-patient cut: quartiles of the **given** locked CLDN4 %pos
scores (same vectors as the n=35 row). Per-patient P uses all malignant
cells → that patient’s T/NK. This is not a Spearman re-audit.

## CellChat-style probability

Jin et al. 2021 Hill / mass-action on CellChatDB v2 protein pairs
(same rule as sibling CLDN4 CellChat folders):

1. Keep a pair only if every ligand and receptor subunit is in the matrix.
2. Per cell group: 10% truncated mean of `log1p(CP10k)`; complexes =
   geometric mean (0 if any subunit mean is 0).
3. Detected if each complex has expressing-cell fraction ≥ 0.10.
4. \(P = (L \cdot R) / (0.5 + L \cdot R)\).
5. Arms with fewer than 10 cells are not scored. Within-patient Q4 vs Q1
   also requires ≥40 malignant cells.

Outgoing only is the ligand table (Mal → T/NK). CellChat R was not run.

## Honest n and skipped files

Given combo n=35 is not the CellChat n. Patients that fail the malignant
or T/NK cell floor are dropped from the LR table and are reported in
`results/patient_coverage.tsv`.

GSE148071 is not downloaded. The other five PR #459 cohorts are not
piled. No FASTQ / EGA / SRA. No file ≥ 2 GB.

## Software

Python (numpy / pandas / scipy / matplotlib / rdata). Public GEO
processed files only: `GSE123902_RAW.tar` (90 MB dense CSV) and
`GSE205335_Lung_IO_UMI_matrix.rds.gz` (524 MB).
