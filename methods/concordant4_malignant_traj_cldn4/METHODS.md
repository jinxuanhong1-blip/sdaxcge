# Methods

Public processed matrices only. No FASTQ. The Laughney 36.5 GB H5 and the GSE131907 log2TPM text matrix are not used.

## Cohorts and malignant gates

| Dataset | Unit | Malignant gate |
| --- | --- | --- |
| GSE123902 | donor | Tumor/metastasis files. `(EPCAM\|KRT8\|KRT18\|KRT19)>0` and `PTPRC==0`. NORMAL libraries dropped. |
| GSE131907 | sample | `Cell_subtype == Malignant cells`. tS1/tS2/tS3 and nLung AT2 are outside this locked gate and are not added. |
| GSE205335 | patient | `lineage.sub == Malignant cells` and tissue does not start with `Normal`. |
| GSE189357 | patient | Same marker gate as GSE123902. |

QC: ≥200 genes, ≥500 counts, mitochondrial percent < 20. Then ≤180 cells per unit (seed 1). The inferential unit stays the patient or sample. Cell counts are not n.

## Embedding

Log-normalize (10,000), 2,000 HVGs (Seurat, `batch_key=dataset`), 30 PCs, Harmony on `dataset`, 30-neighbor graph, UMAP, Leiden resolution 0.6.

## Poles (fixed before the T/NK test)

- AT2 score: SFTPC, SFTPB, SFTPA1, NAPSA, LAMP3, ABCA3. SFTPC and SFTPA1 are absent from the GSE123902 public dense matrix, so they drop out of the shared gene set. The score uses the members that remain (SFTPB, NAPSA, LAMP3, ABCA3). That absence is reported, not filled in.
- Barrier score: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR. **CLDN4 is excluded.**
- IFN score: compact ISG core. **CLDN4 is excluded.**
- Root cluster: highest mean AT2 among Leiden clusters with ≥40 cells, after dropping the cluster with the highest mean CLDN4.
- Early cell: within the root cluster, among cells at or below the cluster median CLDN4, the cell nearest the median AT2 score.
- Barrier cluster: highest mean barrier score among the remaining clusters with ≥40 cells.
- Barrier terminal cell: nearest the 80th percentile of the barrier score in that cluster. CLDN4 is not used.
- Alternate terminal: highest mean IFN among clusters that are neither root nor barrier. Two terminals are required so a fate probability is not identically 1.

## Three algorithms

1. **PAGA** on Leiden (`scanpy.tl.paga`). The reported path is the connectivity shortest path from the root cluster to the barrier cluster.
2. **Palantir** 1.4 on Harmony PCs: diffusion maps (10 components), multiscale space, waypoints, early cell as above, terminals `{barrier cell, alternate cell}`. The primary quantity is the absorption probability of the barrier terminal.
3. **Slingshot** (Street et al. 2018): `pyslingshot.core.getLineages` builds the minimum spanning tree on Leiden centroids in the first 10 Harmony PCs and enumerates root-to-leaf lineages. Principal curves (LOWESS / Hastie–Stuetzle) are fit on the cells of each lineage. The barrier lineage is the curve that ends on the barrier cluster, or the shortest curve that passes through it. This is not diffusion pseudotime. The package `__init__` imports `ggplot2_py`, which is not installed, so the core module is loaded directly. R `slingshot` was not on PATH.

## Association

For each unit with ≥10 malignant cells in the object, take the mean Palantir barrier-fate probability and join `frac_tnk` from `data/locked_patient_units.tsv` (PR #539). That T/NK column is not recomputed.

Report a pooled Spearman and a Fisher-z DerSimonian–Laird meta-analysis across the four cohorts. The same unit-level tests are reported for Slingshot barrier-lineage pseudotime and weight, and for CLDN4 / AT2 / barrier score versus pseudotime.

No RECIST, MPR, or treatment-time model is fit.
