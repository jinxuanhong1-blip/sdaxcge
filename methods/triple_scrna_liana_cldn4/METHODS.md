# Methods — triple-merge CLDN4-only LIANA / CellPhoneDB

## Scope

Additive merge of three public lung scRNA-seq series:

| Accession | Citation | Labels used |
|---|---|---|
| GSE131907 | Kim et al., *Nat Commun* 2020 | Author `Cell_type` / `Cell_subtype` |
| GSE148071 | Wu et al., *Nat Commun* 2021 | Marker-argmax (no GEO labels) |
| GSE205335 | Hui et al. / Lung IO | Author `lineage.sub` / `lineage.total` |

This is **not** GSE207422 and **not** the 131907+205335-only pair.
**CLDN4 only.** TACSTD2 is not a gate and dual-high is not defined.

## Malignant and T/NK

- **GSE131907:** tumor-site cells only (`tLung`, `tL/B`, `mLN`, `mBrain`, `PE`).
  Malignant = `Cell_type == Epithelial cells` and
  `Cell_subtype ∈ {Malignant cells, tS1, tS2, tS3}`.
  T/NK = `T lymphocytes` or `NK cells`. Normal-lung epithelium is not malignant.
- **GSE148071:** A3 marker-argmax lineage. Malignant-like = epithelial AND
  zero UMI of `SFTPA2, AGER, SCGB1A1, SCGB3A1, TPPP3`. T/NK = T or NK.
  This is a marker proxy, not CopyKAT.
- **GSE205335:** `lineage.sub == Malignant cells`; `lineage.total == T/NK cells`.

Patient IDs are prefixed by accession so they cannot collide.

## State split

Among malignant / malignant-like cells in the merge, CLDN4 is
`log1p(CP10k)` using each cell’s full UMI total.
**High** = at or above the **global** median; **low** = below.
A patient enters the paired test if it has ≥10 malignant cells in **both**
bins and ≥20 T/NK cells. The header *n* is that paired count, not the
number of deposited samples.

## Pair list and score

CellPhoneDB v5 (`ventolab/cellphonedb-data`). Hyphenated symbols (HLA-A)
are split only when both sides match the gene list. Complexes use `+`.
A curated T-recruit / IFN / MHC-I overlay is added if a parse miss
would drop a canonical pair.

Partner expression = **minimum subunit mean** on log1p(CP10k).
Pair score = **mean of the two partner means**
(Efremova et al., *Nat Protoc* 2020; Garcia-Alonso et al., *Nat Protoc* 2022).
`pass_expr_prop` = both partners detected in ≥10% of cells in their group.

Primary test = **paired Wilcoxon** of high vs low scores across patients.
FDR = Benjamini–Hochberg within each contrast (outgoing or incoming).
This is **not** a CellChat communication probability.

## Secondary

`liana.mt.cellphonedb` (`resource_name=cellphonedb`, `expr_prop=0.10`,
50 permutations) on groups downsampled to ≤2,000 cells. LIANA p-values
are within-object specificity, not patient-level tests.

CellChat is not run (no R in this environment).

## Software

Python 3.12; numpy, pandas, scipy, statsmodels, matplotlib, rdata.
Optional: anndata + liana.
