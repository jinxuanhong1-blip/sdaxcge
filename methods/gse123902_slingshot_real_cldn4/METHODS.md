# Methods — GSE123902 real Slingshot/PAGA, CLDN4 only

ADDITIVE. Does **not** rewrite pair folders that use GSE123902 as one arm of a combo.
Mouse GSE123903 and SuperSeries GSE123904 are not used. TACSTD2 does not define groups.
No dual-high gate.

**Question.** On Laughney 2020 **tumor epithelium**, where do **CLDN4**, a CLDN4-excluded
barrier/keratin program, and an IFN program sit on a Slingshot/PAGA trajectory whose
root is **not** CLDN4-high?

## Dataset (public processed UMI only)

| File | Public? | Size | Used? |
| --- | --- | ---: | --- |
| GSE123902_RAW.tar (17 dense UMI CSVs) | yes | **90.4 MB** | **yes** |
| Author `PATIENT_LUNG_ADENOCARCINOMA_ANNOTATED.h5` | yes | **36.5 GB** | **no** |

Epithelium gate (locked, PR #459): `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`.
This is **not** a CNV-malignant call. PRIMARY + METASTASIS are the analysis cells.
Matched NORMAL epithelium is kept in the graph only so the root can be AT2-like.

## Trajectory

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → k-NN 30 |
| Clusters | Leiden 0.6; PAGA on Leiden (geometry only) |
| Clock | Real Slingshot (Street 2018). R/Bioconductor if installed; otherwise the same MST + polyline projection. DPT is a companion, not the claim if Slingshot lineages exist. |
| Root | Matched-normal AT2-like, **never CLDN4-high**. Fallback: max AT2 among CLDN4-not-high cells. |
| Scores | CLDN4 continuous + tertiles; barrier/keratin **without CLDN4**; Hallmark IFNα∩IFNγ core |
| Inferential n | **Donor** (LX ID). Cell-level ρ is descriptive. |
| Cap | ≤350 cells / GEO sample after protecting matched-normal AT2-like cells |
| Extra figures | CLDN4 + barrier + IFN along PT; within-donor CLDN4-high vs low |
| Unused | 36.5 GB H5; dual-high TACSTD2∩CLDN4; ICI / MPR / RECIST as a trajectory label |

## Tests

Primary: donor-level Spearman of mean CLDN4 / barrier / IFN vs mean Slingshot PT (tumor epithelium).
BH inside that list only.

Sensitivity (not BH): primary-only, met-only, include-normal-cells.

The lineage table (`results/tables/lineage.tsv`) is the done criterion.
