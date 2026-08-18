# FINDING — Seurat GSE123902 CLDN4-only (patient/donor unit)

ADDITIVE. **CLDN4 only.** Thesis already correct. This is a Seurat-native re-score of public processed human LUAD/NSCLC [GSE123902](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123902) (Laughney et al., *Nat Med* 2020, PMID 32042191). SuperSeries GSE123904 / mouse GSE123903 are not used. No TACSTD2∩CLDN4 dual-high gate. No GSE148071 merge. No Python-only primary.

Primary engine: **R + Seurat 5.5.1**. `CreateSeuratObject` was called on a 10x-style MTX written from the GEO dense unnormalized UMI CSVs (`GSE123902_RAW.tar`, 90.4 MB). The 36.5 GB author annotated H5 was not required. Honest unit = **patient/donor (LX ID)**.

## Verdict

Seurat ran. A patient-level table exists (`results/tables/patient_level_cldn4_tnk.tsv`). Eligible tumor donors ( ≥20 marker-epithelial cells in tumor and ≥20 T/NK): **n=13**.

- Malignant CLDN4 **mean** vs T/NK fraction: n=13, ρ=-0.132 [-0.637, 0.452], p=0.668.
- Malignant CLDN4 **%pos** vs T/NK fraction: n=13, ρ=-0.434 [-0.795, 0.154], p=0.138.
- Malignant CLDN4 mean vs IFN: n=13, ρ=-0.440 [-0.797, 0.147], p=0.133.
- Malignant CLDN4 mean vs MHC-I/APM: n=13, ρ=-0.071 [-0.599, 0.499], p=0.817.
- Malignant CLDN4 mean vs TJ (CLDN4 held out): n=13, ρ=+0.654 [0.161, 0.886], p=0.0153.

Malignant CLDN4 is higher than same-donor T/NK CLDN4 in **13/13** eligible donors (epithelial restriction; not an immune-cold claim).

n=13 is the honest ceiling. Cell-level p-values are not the claim. Do not write this as a failed audit of the thesis.

## Honest n

| item | n | note |
|---|---:|---|
| GEO dense CSVs | 17 | 17 files in GSE123902_RAW.tar |
| MTX barcodes written | 42847 | union gene space 23266 |
| Seurat cells (nCount_RNA>0) | 42847 | CreateSeuratObject from MTX |
| Donors / LX IDs | 14 | filename MSK_LX* |
| Tumor donors | 13 | primary + metastasis; LX685 is normal-only |
| Marker epithelial in tumor (malignant for this table) | 4433 | EPCAM/KRT* vs T/NK/myeloid/B; CLDN4 not used |
| T/NK in tumor samples | 13305 | CD3D/CD3E/CD8A/NKG7/GNLY/KLRD1 |
| Eligible Spearman donors | **13** | ≥20 malignant and ≥20 T/NK |
| Dual-high TACSTD2 ∩ CLDN4 | not defined | CLDN4-only |
| GSE148071 cells | 0 | not merged |
| Author 36.5 GB H5 | not used | GEO CSV → MTX is the public processed object |
| Marker epithelial ≠ CNV-malignant | yes | SEQC dense UMI; no inferCNV |

Paper QC atlas is 41,384 cells. We do not substitute that n. SEQC CSVs hold 42847 barcodes before the empty-barcode drop.

## Gate

| file | public? | used |
|---|---|---|
| `GSE123902_RAW.tar` (17 dense UMI CSVs) | yes, 90.4 MB | **yes** — written to 10x MTX, then `ReadMtx` + `CreateSeuratObject` |
| 10x MTX (`results/mtx/`) | derived here | **yes** — Seurat input |
| Author `PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5` | yes, 36.5 GB | no — not required once GEO UMI CSVs exist |
| GSE148071 / other cohorts | — | **no** |

## Locked choices

- Lineage is a four-way marker argmax on Seurat `LogNormalize` (log1p CP10k). Keep if top ≥ 0.12 and top ≥ 1.15 × second. CLDN4 is never a lineage marker.
- Malignant = marker epithelial **in tumor** (primary or metastasis). Matched normal is not scored as malignant.
- T/NK fraction = marker T/NK / tumor-sample cells of that donor.
- IFN = locked Hallmark IFNα ∩ IFNγ core (present 27/27).
- MHC-I/APM = curated antigen-presentation set (present 21/21).
- TJ = epithelial tight-junction genes with **CLDN4 held out** (present 14/14).
- Eligible Spearman n requires ≥20 malignant and ≥20 T/NK. Q4 vs Q1 requires eligible n ≥ 8.

## Patient / donor table

Machine table: `results/tables/patient_level_cldn4_tnk.tsv`.

| patient | tumor site | n tumor | n mal | n T/NK | frac T/NK | CLDN4 mean | CLDN4 %pos | IFN | MHC | TJ | eligible |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| MSK_LX255B | METASTASIS | 3886 |  205 | 1472 | 0.379 | 0.718 | 40.0 | 0.293 | 0.619 | 0.163 | yes |
| MSK_LX653 | PRIMARY_TUMOUR |  392 |  154 |   34 | 0.087 | 1.122 | 69.5 | 0.328 | 0.709 | 0.316 | yes |
| MSK_LX661 | PRIMARY_TUMOUR | 4566 |  117 | 3023 | 0.662 | 0.942 | 65.0 | 0.396 | 0.895 | 0.322 | yes |
| MSK_LX666 | METASTASIS | 1316 |  790 |  140 | 0.106 | 0.538 | 67.7 | 0.431 | 1.040 | 0.097 | yes |
| MSK_LX675 | PRIMARY_TUMOUR | 3435 |  514 | 1033 | 0.301 | 0.594 | 39.9 | 0.400 | 0.888 | 0.194 | yes |
| MSK_LX676 | PRIMARY_TUMOUR | 4014 |  140 | 2388 | 0.595 | 0.901 | 47.9 | 0.376 | 0.896 | 0.192 | yes |
| MSK_LX679 | PRIMARY_TUMOUR | 2017 |  464 |  568 | 0.282 | 0.748 | 48.3 | 0.426 | 0.907 | 0.372 | yes |
| MSK_LX680 | PRIMARY_TUMOUR |  941 |  265 |  298 | 0.317 | 1.475 | 76.6 | 0.268 | 0.834 | 0.305 | yes |
| MSK_LX681 | METASTASIS | 1630 | 1249 |   70 | 0.043 | 0.860 | 61.9 | 0.363 | 0.473 | 0.269 | yes |
| MSK_LX682 | PRIMARY_TUMOUR | 4092 |  259 | 1970 | 0.481 | 0.661 | 29.0 | 0.619 | 1.243 | 0.138 | yes |
| MSK_LX684 | PRIMARY_TUMOUR |  378 |  182 |   43 | 0.114 | 1.095 | 81.3 | 0.314 | 0.864 | 0.362 | yes |
| MSK_LX685 | NONE |    0 |    0 |    0 | NA | NA | NA | NA | NA | NA | no |
| MSK_LX699 | METASTASIS | 2671 |   37 | 2106 | 0.788 | 0.472 | 24.3 | 0.285 | 0.678 | 0.103 | yes |
| MSK_LX701 | METASTASIS |  788 |   57 |  160 | 0.203 | 0.458 | 22.8 | 0.399 | 0.715 | 0.246 | yes |

## Donor-level Spearman (primary)

| Contrast | n | ρ [95% CI] | p |
|---|---:|---|---:|
| CLDN4 mean vs T/NK fraction | 13 | -0.132 [-0.637, 0.452] | 0.668 |
| CLDN4 %pos vs T/NK fraction | 13 | -0.434 [-0.795, 0.154] | 0.138 |
| CLDN4 mean vs IFN | 13 | -0.440 [-0.797, 0.147] | 0.133 |
| CLDN4 mean vs MHC-I/APM | 13 | -0.071 [-0.599, 0.499] | 0.817 |
| CLDN4 mean vs TJ (CLDN4 held out) | 13 | +0.654 [0.161, 0.886] | 0.0153 |

## Malignant Q4 vs Q1 IFN / MHC / TJ

Malignant CLDN4-mean quartiles among eligible donors. Q1 n=4 vs Q4 n=4 on malignant CLDN4 mean.

| family | n_Q1 | n_Q4 | Δ median (Q4−Q1) | MW p |
|---|---:|---:|---:|---:|
| IFN | 4 | 4 | -0.078 | 0.194 |
| MHC-I/APM | 4 | 4 | +0.048 | 1 |
| TJ (CLDN4 held out) | 4 | 4 | +0.170 | 0.0304 |

## Seurat plots

- `results/figures/fig_dimplot_lineage.png` — `DimPlot` by marker lineage
- `results/figures/fig_dimplot_patient.png` — `DimPlot` by donor
- `results/figures/fig_dimplot_tissue.png` — `DimPlot` by tissue
- `results/figures/fig_vlnplot_cldn4.png` — `VlnPlot` CLDN4 by lineage
- `results/figures/fig_featureplot_cldn4.png` — `FeaturePlot` CLDN4
- `results/figures/fig_patient_cldn4_vs_tnk.png` — donor scatter (mean)
- `results/figures/fig_patient_cldn4pct_vs_tnk.png` — donor scatter (%pos)

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- Marker epithelial is **not** a CNV-malignant call.
- n=13 tumor donors is small; Fisher-z CIs are wide.
- Restriction (CLDN4 in epithelium vs T/NK) is not an infiltration / immune-cold claim.
- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.
- No ICI / MPR / RECIST / survival test.
- Not merged with GSE148071 or any other cohort.

## Reproduce

```bash
bash methods/seurat_gse123902_cldn4/scripts/download.sh
Rscript methods/seurat_gse123902_cldn4/scripts/run_seurat_gse123902.R
```

Seurat 5.5.1. R 4.6.1.

