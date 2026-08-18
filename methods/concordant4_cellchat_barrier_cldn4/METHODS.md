# Methods — concordant-four CLDN4-only CellChat barrier call

ADDITIVE. Thesis already correct. CLDN4 only. No dual-high.
Concordant four: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Do not add GSE148071 / GSE127465 / CD45-only. Not a full-pool.
Patient/donor/locked-sample is the unit.
PR #424 / #474 numbers are comparison rows only.

## Matrices

Public processed GEO only (<2 GB):

- GSE123902 `GSE123902_RAW.tar` (Laughney 2020 dense UMI CSVs).
- GSE131907 raw UMI + author cell annotation. The 2.86 GB log2TPM text
  and EGA FASTQ are skipped. `PVRL2` is aliased to `NECTIN2`.
- GSE205335 UMI RDS + CellIdentity + family SOFT (Hu 2023).
- GSE189357 `GSE189357_RAW.tar` (Zhu/Wang AIS–IAC 10x MTX; TD1–TD9).

## Malignant / T/NK

- GSE123902 and GSE189357: author malignant labels thin → epithelium
  marker-malignant = EPCAM|KRT8|KRT18|KRT19 count>0 and PTPRC==0.
  T/NK = CD3D|CD3E|CD8A|NKG7|GNLY|KLRD1 count>0 and not malignant.
  GSE123902 is donor-level (PRIMARY preferred over METASTASIS).
- GSE131907: author `Cell_subtype` in {Malignant cells, tS1, tS2, tS3};
  T/NK = `Cell_type` in {T lymphocytes, NK cells}. Locked samples only.
- GSE205335: author `lineage.sub` == Malignant cells;
  T/NK = `lineage.total` == T/NK cells. Locked patients only.

TACSTD2 is never a gate.

## Gates

CLDN4-only. Primary: within-unit malignant **Q4 vs Q1**.
Extra: %pos (CLDN4>0 vs =0), median, tertile.
Honest n requires both compartments: ≥10 cells/arm and T/NK ≥20.
Q4 vs Q1 also requires n_mal ≥ 40.

## Pre-specified families

Not a discovery screen. CellChatDB v2 protein pairs:

1. Barrier/inhibitory (expect CLDN4-high > low outgoing → T/NK):
   F11R–ITGAL/ITGB2, NECTIN2–TIGIT, CDH1–ITGAE/ITGB7, CDH1–KLRG1,
   LGALS9–HAVCR2, LGALS9–CD44, LGALS9–PTPRC.
2. IFN / T-recruit / MHC-I (expect CLDN4-low / KD-like > high):
   CXCL9–CXCR3, CXCL10–CXCR3, CCL5–CCR5, CCL5–CCR1,
   HLA-A–CD8A, HLA-B–CD8A, HLA-C–CD8A.

Barrier-up-in-high is ON-thesis. CXCL9/10 are often undetected at
expr_prop≥0.10; say so in one line and still report the high-side
barrier ligands.

## Score

Outgoing Hill probability (Jin et al. 2021): 10% truncated mean,
Kh=0.5, detected if ligand and receptor expr_prop ≥ 0.10.
ΔP = P_high − P_low. Family aggregate = mean of detected pair ΔPs
in that unit. Honest n = units that pass the cell floor and have a
detected pair. CellChat R / LIANA R were not run.

## Not done

Dual-high, GSE148071, GSE127465, CD45-only, 7-pool, T/NK rho re-audit,
PR #424/#474 re-audit, cell-pooled tests, discovery screen.
