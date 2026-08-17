# Methods — GSE131907 LIANA/LR from CLDN4-high malignant to T/NK

## Data

- Accession: **GSE131907** (Kim et al., *Nat Commun* 2020, PMID 32385277).
- Input: author `GSE131907_Lung_Cancer_cell_annotation.txt.gz` and processed
  `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz` (~0.38 GB gzip) from GEO
  Series supplementary files. Series matrix supplies `patient_id`.
- Skipped: 2.86 GB log2TPM text; EGA raw FASTQ (`EGAD00001005054`).
- Platform: 10x Chromium. Treatment-naive LUAD atlas. **No ICI / MPR labels.**

## Labels (author, not recomputed)

| Set | Rule |
|---|---|
| Tumor-origin samples | `Sample_Origin` ∈ {tLung, tL/B, mLN, mBrain, PE} |
| Malignant | tumor-origin **and** `Cell_subtype` ∈ {Malignant cells, tS1, tS2, tS3} |
| T/NK | `Cell_type` ∈ {T lymphocytes, NK cells} |

Honest limits on that malignant call:

- Primary tLung uses Kim **tS1/tS2/tS3** (tumor-specific epithelial clusters).
  The author string `Malignant cells` is used in metastases / tL-B / mLN / mBrain.
- PE epithelial cells are present but **unlabeled**; they are **not** counted as malignant.
- nLung AT1/AT2/Club/Ciliated are **not** malignant.
- This is not a re-run of inferCNV/CopyKAT.

## CLDN4 split

Among malignant cells, CLDN4 is `log1p(CP10k)` using the cell’s full UMI total
(one pass over the processed matrix). **High** = at or above the **global
median** of malignant CLDN4; **low** = below. TACSTD2 is stored but does not
define the split.

A tumor sample enters the paired test if it has ≥10 malignant cells in both
bins and ≥20 T/NK cells. The inferential unit is the **sample**. Unique
`patient_id` in that set is reported separately (GEO has 58 samples / 44
patients; some patients have more than one site). Cells are not n.

## Pair list

CellPhoneDB v5 interactors parsed as in the GSE207422 folder
(`resources/cellphonedb_v5_lr_pairs.tsv`), plus a small curated T-recruit /
IFN / MHC-I overlay. Pathway labels organize ranks; they do not subset the
score.

## Documented score (primary)

For sender group S and receiver group R, each partner’s expression is the
**minimum of subunit means** on log1p(CP10k). The pair score is the
**arithmetic mean of the two partner means** (Efremova et al., *Nat Protoc*
2020; Garcia-Alonso et al., *Nat Protoc* 2022). A pair **passes `expr_prop`**
when both partners are detected in ≥10% of cells in their group.

- `lr_table.tsv` — pooled CLDN4-high malignant → tumor-origin T/NK (descriptive).
- `lr_table_high_vs_low.tsv` — paired Wilcoxon of high vs low per sample; FDR-BH
  within the outgoing contrast.

This is **not** a CellChat communication probability.

## LIANA (secondary, if installed)

`liana.mt.cellphonedb` with `resource_name="cellphonedb"`, `expr_prop=0.10`,
`n_perms=50`. Groups: `Malig_CLDN4high` / `Malig_CLDN4low` / `T` / `NK`.
Each group is downsampled to ≤2,000 cells (`random_seed=0`). LIANA p-values
are within-object specificity, not patient-level tests.

## CellChat

Not run. This environment has no R.

## Software

Python 3; numpy, pandas, scipy, statsmodels, matplotlib, pyyaml. Optional
scanpy + liana. See `requirements.txt`.
