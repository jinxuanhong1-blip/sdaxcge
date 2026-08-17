# Methods — GSE148071 CLDN4-high malignant → T/NK ligand–receptor

Additive package. Does not edit other `methods/` folders.

## Data

- Accession: **GSE148071** (Wu et al., *Nature Communications* 2021, PMID 33953163).
- Object used: public TISCH2 `NSCLC_GSE148071` (`expression.h5` + `CellMetainfo_table.tsv`) from
  `https://tisch.compbio.cn/static/data/NSCLC_GSE148071/`.
- Values in the h5 are the TISCH2/MAESTRO **log2(TPM/10 + 1)** matrix. This is not a
  CellRanger reprocess of GEO `GSE148071_RAW.tar`.
- All 82,267 TISCH2 cells are Tissue=`Tumor`. 42 patients, one Sample each.

## Labels (taken as given)

| Role | Rule |
|---|---|
| Malignant | TISCH2 `Celltype (major-lineage)` == `Malignant` |
| T/NK | TISCH2 major-lineage in {CD8T, Tprolif, CD4Tconv, Treg, Tcell, NK, NKT, ILC, TMKI67, CD8Tex} |
| CLDN4-high | malignant cells with CLDN4 ≥ **global median** among malignant cells |

In this object the T/NK set collapses to **CD8T + Tprolif**. NK / CD4 / Treg bins are empty.
That is reported, not filled in.

## Pair list

CellPhoneDB v5 `interactors` parsed as in `methods/scrna_cellchat_gse207422`
(`resources/cellphonedb_v5_lr_pairs.tsv`). Complexes use `+` subunits.
Pathway labels (T_recruit, IFN, MHC_I, checkpoint) organize ranks; they do not subset the score.

## Documented score (primary)

For sender group S and receiver group R, each partner’s expression is the **minimum of
subunit means** on the TISCH2 log matrix (CellPhoneDB complex rule). The pair score is the
**arithmetic mean of the two partner means** (Efremova et al., *Nat Protoc* 2020;
Garcia-Alonso et al., *Nat Protoc* 2022). A pair **passes `expr_prop`** when both partners
are detected (>0) in ≥10% of cells in their group.

Pooled scores (all patients together) are descriptive. The test, when n allows, is a
**paired Wilcoxon** of the high vs low score across patients that have
≥10 CLDN4-high malignant, ≥10 CLDN4-low malignant, and ≥20 T/NK cells.
FDR is Benjamini–Hochberg within that contrast. Cells are not replicates.

This is **not** a CellChat communication probability.

## LIANA (secondary, if installed)

`liana.mt.cellphonedb` with `resource_name="cellphonedb"`, `expr_prop=0.10`, `n_perms=50`.
Groups: `Malig_CLDN4high` / `Malig_CLDN4low` / `TNK`. Each group is downsampled to
≤2,000 cells (`random_seed=0`). LIANA p-values are within-object specificity, not
patient-level tests.

## CellChat

Not run. This environment has no R. No CellChat tables are written.

## Software

Python 3; numpy, pandas, scipy, statsmodels, matplotlib, h5py, pyyaml.
Optional: scanpy + liana. See `requirements.txt`.
