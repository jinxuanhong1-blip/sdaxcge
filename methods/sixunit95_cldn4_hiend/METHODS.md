# Methods — 6-unit CLDN4-only high-end CellChat

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. Patient is the unit.

## Locked combo (not re-audited)

The strict-malignant 6-unit mean-vs-T/NK Spearman is taken as given from
PR #312 (`methods/strict_malig_tnk`): **n=95, ρ=−0.260, p=0.0195, I²=0%**.

| unit | malignant def | n |
|---|---|---:|
| GSE207422 | author DRMref | 12 |
| GSE205335 | author malignant | 22 |
| GSE291670 | marker malignant | 6 |
| GSE253013 | marker malignant-like (Tumor) | 9 |
| GSE131907 | author malignant, n_mal≥20 | 21 |
| GSE325414 | author malignant | 25 |

All-epithelial rows were already dropped in that recut. This folder does
not re-rank or re-pool that Spearman.

## What is new

CellChat-style outgoing communication **CLDN4-high malignant → same-patient T/NK**,
then a meta-analysis of **patient ΔP** across the six units.

High-end split (primary): within each patient, malignant cells are sorted
on `log1p(CP10k)` **CLDN4** and the top vs bottom quarter are kept
(middle half dropped). Quantile cuts collapse when CLDN4 is zero-inflated,
so the split is a stable rank slice, not `qcut` on raw values. Sensitivity:
median (halves) and tertile (top vs bottom third). TACSTD2 is not used.

ΔP = P_high − P_low. One delta per patient per ligand–receptor pair.
Unit means are DerSimonian–Laird random-effects pooled. A patient-level
Wilcoxon on the stacked deltas is also reported. p-values are descriptive.

Extra between-patient cut: quartiles of the **given** locked CLDN4 scores
(same vectors as the n=95 row). Per-patient P uses all malignant cells →
that patient’s T/NK. This is not a Spearman re-audit.

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

GSE253013’s only public processed matrix is `GSE253013_all_luad_garnett_temp.rds.gz`
(9.3 GB). It is **not downloaded**. That unit stays in the given n=95 row
and contributes 0 patients to the LR meta.

GSE131907 2.86 GB log2TPM and all FASTQ/EGA/SRA files are not used.
GSE207422 DRMref barcodes are not public; epithelial lineage is the
malignant proxy on the same 12 locked samples.

## Software

Python (numpy / pandas / scipy / matplotlib / rdata). Public GEO
processed files only.
