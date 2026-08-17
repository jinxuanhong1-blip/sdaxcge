# Methods — thesis-aligned CLDN4-only ligand families (GSE189357 + GSE205335)

ADDITIVE. CLDN4 only. No dual-high. Patient is the unit.
The PR #459 %pos Spearman is not re-audited.

## Given cut

PR #459 pair that differs: GSE189357+GSE205335 %pos n=31, ρ=−0.478, Q4 vs Q1 r=−0.750.
Locked units: `data/GSE189357_marker_units.tsv` (patient, n=9) and
`data/GSE205335_patients.tsv` (patient, n=22 with author malignant + T/NK).

## Matrices

Public processed GEO only (<2 GB):

- GSE189357 `GSE189357_RAW.tar` (Zhu/Wang AIS–IAC 10x MTX, ~624 MB).
- GSE205335 UMI RDS + author cell identity + GEO SOFT. No FASTQ / EGA.

## Malignant / T/NK

- GSE189357: marker-malignant = EPCAM|KRT8|KRT18|KRT19 count>0 and PTPRC==0.
  T/NK = CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1 count>0 and not malignant.
  Patient-level (TD1–TD9).
- GSE205335: author `lineage.sub` == `Malignant cells`;
  T/NK = `lineage.total` == `T/NK cells`. Locked patients only.

Malignant definitions are **not the same**. TACSTD2 is never a gate.

## Pre-specified families

Not a discovery screen. CellChatDB v2 protein pairs only:

1. Barrier/inhibitory (expect CLDN4-high > low outgoing → T/NK):
   F11R–ITGAL/ITGB2, NECTIN2–TIGIT, CDH1–ITGAE/ITGB7, CDH1–KLRG1,
   LGALS9–HAVCR2, LGALS9–CD44, LGALS9–PTPRC.
2. IFN / T-recruit / MHC-I (expect CLDN4-low / KD-like > high):
   CXCL9–CXCR3, CXCL10–CXCR3, CCL5–CCR5, CCL5–CCR1,
   HLA-A–CD8A, HLA-B–CD8A, HLA-C–CD8A.

## Score

Within-unit CLDN4-high vs low malignant cells (primary: median split).
Outgoing Hill probability (Jin et al. 2021): 10% truncated mean, Kh=0.5,
detected if ligand and receptor expr_prop ≥ 0.10. ΔP = P_high − P_low.
Each unit contributes one delta. Family aggregate = mean of detected pair
ΔPs in that unit. Honest n = number of units that pass the cell floor
and have a detected pair.

Cell floors: ≥10 cells/arm and ≥10 T/NK. Q4 vs Q1 extra also requires
n_mal ≥ 40. Tertile extra uses exclusive outer thirds.

## Not done

CellChat R, LIANA, dual-high, T/NK rho re-audit, cell-pooled tests.
