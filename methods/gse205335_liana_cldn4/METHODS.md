# Methods — GSE205335 CLDN4-only ligand–receptor scoring

ADDITIVE. **CLDN4 only.** TACSTD2 is not a gate and does not define dual-high.
Patient is the unit. This is not GSE207422 (that LIANA is PR #344) and is not
the between-patient Q4 vs Q1 CellChat-style table (PR #362).

## Data

- Accession: **GSE205335** (Hu / Ahn / Lee; GEO public 2024). Palliative ICI
  biopsy / effusion / LN / liver scRNA.
- Input: author-processed UMI `GSE205335_Lung_IO_UMI_matrix.rds.gz`
  (`dgCMatrix`, 33,714 genes × 96,505 cells),
  `GSE205335_Lung_IO_CellIdentity.txt.gz`, and `GSE205335_family.soft.gz`.
- Controlled EGA raw (`EGAD00001008703`) was not accessed.
- Platform: 10x 3' / 5' mix (Illumina HiSeq 2500 / NovaSeq 6000).

## Labels (author, not recomputed)

| Set | Rule |
|---|---|
| Tumor samples | GEO `tissue` does **not** start with `Normal` |
| Malignant | tumor sample **and** `lineage.sub == "Malignant cells"` |
| T | `lineage.sub` ∈ {CD4+ T cells, CD8+ T cells} |
| NK | `lineage.sub == "NK cells"` |
| T/NK | `lineage.total == "T/NK cells"` (exactly T + NK here) |

Normal-only samples (P2001 / P2009 / P2016 Normal LN; P3032 Normal Brain;
P0031 Normal Lung) are excluded from the malignant split and from that
patient’s T/NK pool. P0031 tumor lung is kept. This is the authors’ malignant
label, not an independent CNV re-call.

MPR/NMPR is not labeled. RECIST is recorded and is **not** substituted for MPR.

## State split (CLDN4 only)

Among author-malignant cells on tumor samples, CLDN4 is `log1p(CP10k)` using
the cell’s full UMI total. **High** = at or above the **global median**;
**low** = below. TACSTD2 is not used.

A **patient** enters the paired test if, after pooling that patient’s tumor
samples, it has ≥10 malignant cells in both bins and ≥20 T/NK cells. Cells
are not n. Patients with multiple sites (LN + liver, lung + effusion) are
one unit.

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

Python 3; numpy, pandas, scipy, statsmodels, matplotlib, seaborn, pyyaml,
rdata. Optional scanpy + liana. See `requirements.txt`.
The ~500 MB GEO matrix is not stored in git.
