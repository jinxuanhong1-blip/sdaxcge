# Methods — GSE127465 REAL Slingshot/PAGA, CLDN4 only

ADDITIVE. Does **not** rewrite the GSE127465+GSE148071 T/NK pair folder.
**CLDN4 only.** No TACSTD2∩CLDN4 dual-high gate.

**Question.** On GSE127465 (Zilionis et al., *Immunity* 2019, PMID 30979687) **tumor epithelium**, if author malignant cells with CLDN4>0 exist, where does **CLDN4** sit on a Type II–rooted **real Slingshot** lineage relative to AT2 and a CLDN4-excluded barrier/keratin program?

If malignant+CLDN4 is missing: **stop**, write honest n=0, do not invent a lineage.

## Dataset

| Item | Choice |
| --- | --- |
| Accession | GSE127465 human inDrops (mouse MTX unused) |
| Matrix | `GSE127465_human_counts_normalized_54773x41861.mtx.gz` (author total-count normalized) |
| Tumor epithelium | Type I, Type II, club, ciliated, `PatientN-specific` |
| Malignant | `Major cell type` contains `specific`, Tissue=tumor |
| Dropped | blood; immune; fibroblasts; endothelium; smooth muscle |
| Patients | p1–p7. **n=7 is thin.** |

## Trajectory

| Item | Choice |
| --- | --- |
| Clock | **REAL** R `slingshot` (Street et al. 2018) on PCA + Leiden. Not a DPT stand-in. |
| Graph | Seurat HVG 3000 → PCA → k-NN 30. No Harmony (malignant is PatientN-specific). |
| Clusters | Leiden 0.6; PAGA on Leiden (geometry only) |
| Start / root | Leiden with max Type II fraction / AT2 score. **Never the CLDN4-high cluster.** |
| Companion | scanpy DPT rooted on a Type II cell in that start cluster |
| Score | CLDN4 continuous + tertiles. Barrier/keratin **excludes CLDN4**. |
| Inferential n | patient (n=7, thin). Cell-level ρ is descriptive. |
| Extra figures | paired tertile; lineage CLDN4; PAGA; patient scatter; honest n |

Slingshot MST connects cluster centers even when patient tumors are discrete. PAGA connectivity is the honesty check.
