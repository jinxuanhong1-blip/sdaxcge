# Methods — GSE207422 CLDN4-only ligand–receptor scoring

## Data

- Accession: **GSE207422** (Hu et al., *Genome Medicine* 2023, PMID 36869384).
- Input: author-processed UMI matrix
  `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` and
  `GSE207422_NSCLC_scRNAseq_metadata.xlsx` from GEO Series supplementary files.
- Platform: BD Rhapsody WTA (not 10x). 15 samples, neoadjuvant PD-1 + chemo.
- Raw GSA-Human HRA001033 is not used. Author per-cell labels are not on GEO.

## A3 rules (taken as given)

Lineage = argmax of mean log1p marker scores, with the A3 floor
(`best < 0.15` or `best − second < 0.05` → `other`):

| Lineage | Markers |
|---|---|
| epithelial | EPCAM, KRT8, KRT18, KRT19 |
| T | CD3D, CD3E, CD2 |
| NK | NKG7, GNLY, FGFBP2 |
| B / plasma / myeloid / neutrophil / fibroblast / endothelial / mast | as in A3 |

**Malignant-like** = epithelial AND zero UMI of SFTPA2, AGER, SCGB1A1, SCGB3A1, TPPP3.
T/NK = T or NK. This is a marker proxy, not public CopyKAT calls.

## State split (CLDN4 only)

Among malignant-like cells, CLDN4 is scored as log1p(CP10k) using the cell’s
full UMI total. **High** = at or above the global median; **low** = below.
TACSTD2 is not used to define dual-high or to filter patients.

A patient enters the paired test if it has ≥10 malignant-like cells in both
bins and ≥20 T/NK cells.

## Pair list

CellPhoneDB v5 (`ventolab/cellphonedb-data`, `interaction_input.csv`).
The `interactors` string is split into ligand and receptor only when both
sides match the CellPhoneDB gene list, so hyphenated symbols such as HLA-A
are preserved. Complexes use `+` subunits. A small curated overlay adds
canonical T-recruit / IFN / MHC-I pairs if a parse miss would drop them.

Pathway labels (T_recruit, IFN, MHC-I) organize the ranks. They do not
subset the score.

## Documented score (primary)

For sender group S and receiver group R, each partner’s expression is the
**minimum of subunit means** on log1p(CP10k) (CellPhoneDB complex rule).
The pair score is the **arithmetic mean of the two partner means**
(Efremova et al., *Nat Protoc* 2020; Garcia-Alonso et al., *Nat Protoc* 2022).
A pair **passes `expr_prop`** when both partners are detected in ≥10% of cells
in their group.

Pooled (all patients together) scores are descriptive. The test is a
**paired Wilcoxon** of the high vs low score across patients. FDR is
Benjamini–Hochberg within each contrast (outgoing or incoming).

This is **not** a CellChat communication probability.

## LIANA (secondary, if installed)

`liana.mt.cellphonedb` with `resource_name="cellphonedb"`, `expr_prop=0.10`,
`n_perms=50`. Groups: `Malig_CLDN4high` / `Malig_CLDN4low` / `T` / `NK`.
Each group is downsampled to ≤2,000 cells (`random_seed=0`).
LIANA p-values are within-object specificity, not patient-level tests.

## CellChat

Not run. This environment has no R. No CellChat RDS, probability matrix, or
pathway barplot is produced.

## Software

Python 3.12; numpy, pandas, scipy, statsmodels, matplotlib; optional
scanpy + liana. See `requirements.txt`.
