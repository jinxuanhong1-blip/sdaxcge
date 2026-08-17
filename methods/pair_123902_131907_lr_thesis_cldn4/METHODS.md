# Methods — thesis-aligned CLDN4-only ligand families (GSE123902 + GSE131907)

ADDITIVE. CLDN4 only. No dual-high. No GSE148071.
Patient/donor is the unit. The PR #459 %pos Spearman is not re-audited.

## Given cut

PR #459 pair that differs: GSE123902+GSE131907 %pos n=34, ρ=−0.575, Q4 vs Q1 r=−0.700.
Locked units: `data/GSE123902_marker_units.tsv` (donor, n=13) and
`data/GSE131907_samples.tsv` (sample, n=21 with author malignant + T/NK).

## Matrices

Public processed GEO only (<2 GB):

- GSE123902 `GSE123902_RAW.tar` (Laughney 2020 dense UMI CSVs).
- GSE131907 raw UMI + author cell annotation. The 2.86 GB log2TPM text and
  EGA FASTQ are skipped. The 36.5 GB Laughney H5 is skipped.

## Malignant / T/NK

- GSE123902: marker-malignant = EPCAM|KRT8|KRT18|KRT19 count>0 and PTPRC==0.
  T/NK = CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1 count>0 and not malignant.
  Donor-level; PRIMARY preferred when a donor also has METASTASIS. Normals dropped.
- GSE131907: author `Cell_subtype` in {Malignant cells, tS1, tS2, tS3};
  T/NK = `Cell_type` in {T lymphocytes, NK cells}. Locked samples only.

TACSTD2 is never a gate.

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

CellChat R, LIANA, dual-high, GSE148071, T/NK rho re-audit, cell-pooled tests.
