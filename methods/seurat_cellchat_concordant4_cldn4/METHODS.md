# Methods — Seurat + CellChat concordant-four CLDN4-high vs low senders

ADDITIVE. Thesis already correct. CLDN4 only. No dual-high.
Concordant four: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Do not add GSE148071 / GSE127465 / CD45-only. Not a full-pool.
Patient / locked sample is the unit.

## Engine

R + Seurat `CreateSeuratObject` + CellChat `createCellChat` /
`computeCommunProb`. This is not a Python reimplementation of the Hill
score. If Seurat or CellChat cannot install, the analysis stops.

CellChat settings: `type = "truncatedMean"`, `trim = 0.1`,
`population.size = TRUE`. CellChatDB is subset to the pre-specified
protein pairs below.

## Matrices

Public processed GEO only. Prefer MTX / dense CSV / UMI TXT over RDS.

- GSE123902 `GSE123902_RAW.tar` (Laughney 2020 dense UMI CSVs).
- GSE131907 raw UMI TXT + author cell annotation. The 2.9 GB log2TPM
  text and both GEO RDS files are skipped. `PVRL2` is aliased to
  `NECTIN2`.
- GSE205335 UMI RDS + CellIdentity + family SOFT (Hu 2023). GEO does
  not ship MTX/H5/CSV for this series. The RDS is read with a
  multi-layer gzip peeler (files may be double-gzipped). Genes are
  subset immediately after `readRDS`.
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
  Normal-tissue biopsies are excluded.

TACSTD2 is never a gate.

## Gates

CLDN4-only. Primary: within-unit malignant **Q4 vs Q1**.
Extra: median split and %pos (CLDN4>0 vs =0).
Honest n requires both compartments: ≥10 cells/arm and T/NK ≥20.
Q4 vs Q1 also requires n_mal ≥ 40.

Per unit, CellChat is run on three identities: `CLDN4_high`,
`CLDN4_low`, `TNK`. The readout is outgoing probability
high→TNK minus low→TNK.

## Pre-specified families

Not a discovery screen.

1. Barrier/inhibitory (expect CLDN4-high > low outgoing → T/NK):
   F11R–ITGAL/ITGB2 (`JAM1_ITGAL_ITGB2`), NECTIN2–TIGIT, CDH1–ITGAE/ITGB7,
   CDH1–KLRG1, LGALS9–HAVCR2, LGALS9–CD44, LGALS9–PTPRC (`LGALS9_CD45`).
2. IFN / T-recruit / MHC-I (expect CLDN4-low / KD-like > high):
   CXCL9–CXCR3, CXCL10–CXCR3, CCL5–CCR5, CCL5–CCR1,
   HLA-A–CD8A, HLA-B–CD8A, HLA-C–CD8A.

Barrier-up-in-high is ON-thesis. CXCL9/10 are often undetected.

## Score

ΔP = P_high→TNK − P_low→TNK from CellChat. Family aggregate = mean of
detected pair ΔPs in that unit. Honest n = units that pass the cell
floor and have a detected pair. Wilcoxon signed-rank vs 0 across units.

## Not done

Dual-high, GSE148071, GSE127465, CD45-only, 7-pool, cell-pooled tests,
discovery screen, Python Hill-score primary.
