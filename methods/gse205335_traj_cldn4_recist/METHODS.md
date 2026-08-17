# Methods — GSE205335 PAGA / DPT, CLDN4 and RECIST

ADDITIVE. **CLDN4 only.** GSE205335 advanced NSCLC ICI (Ahn / Lee, *eLife* 2024).
Public processed UMI only. **No GSE148071. No dual-high.** Does not re-audit
PR #320 T/NK ρ.

**Question.** Among malignant / epithelial cells, does CLDN4 sit on a
trajectory that also separates RECIST (**PR vs PD/SD**)? Patient is the unit.

## Dataset

| Item | Choice |
| --- | --- |
| Accession | GSE205335 |
| Matrix | `GSE205335_Lung_IO_UMI_matrix.rds.gz` (dgCMatrix UMI) |
| Labels | `GSE205335_Lung_IO_CellIdentity.txt.gz` + GEO SOFT |
| Epithelium | `lineage.total == Epithelial cells` on non-normal tissues |
| Dropped | Normal Lung / LN / Brain; EGA FASTQ; GSE148071; GSE131907 |
| RECIST | GEO `recist`; test = PR vs PD+SD; NE excluded from that table |
| Histology | GEO `cancer subtype` (ADC / SQ / SCLC / NUT), reported not hidden |

## Trajectory clock

Slingshot (Street et al. 2018) if `Rscript` + `slingshot` are present.
Otherwise the **documented fallback** is scanpy diffusion pseudotime
(Haghverdi et al. 2016). Palantir is not required.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → k-NN 30. No Harmony (single cohort). |
| Clusters | Leiden 0.6; PAGA on Leiden (geometry only) |
| Root | Biological: author Non-malignant if that pool is not CLDN4-high; else the Leiden cluster with highest AT2 score among clusters with mean CLDN4 ≤ object median. **Never CLDN4-high.** |
| Score | CLDN4 continuous + tertiles |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| Inferential n | Patient. Cell-level ρ is descriptive. |
| Cap | ≤500 cells / patient (catalog vs analysis n both reported) |
| Extra figures | Within-patient CLDN4-high vs low; patient CLDN4 vs DPT colored by RECIST + histology |

## Tests

1. Patient-level Spearman of mean CLDN4 vs mean DPT / AT2 / barrier (CLDN4 excluded). BH inside that list only.
2. Patient-level Mann–Whitney of mean DPT (and CLDN4) in **PR vs PD/SD**. Honest n = patients with evaluable RECIST.
3. Sensitivity (not BH): malignant-only means, ADC+SQ, drop-SCLC, RECIST-evaluable-only Spearman.

If CLDN4 vs DPT is null, that is **not** lineage proof. The RECIST table is still delivered.
