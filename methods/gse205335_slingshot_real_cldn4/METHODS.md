# Methods — GSE205335 malignant-only REAL Slingshot/PAGA, CLDN4 only

ADDITIVE. Does **not** rewrite the winning-pair epithelium trajectory
(`methods/winpair_131907_205335_slingshot_cldn4/`) or any TACSTD2 hunt.
GSE131907 and GSE207422 are not added. TACSTD2 does not define groups.
No dual-high gate.

**Question.** On GSE205335 **author-malignant cells only**, where do
**CLDN4**, a CLDN4-excluded barrier/keratin program, and a compact IFN
ISG score sit on a real Slingshot/PAGA lineage, and does patient-mean
diffusion / Slingshot pseudotime differ by RECIST?

## Dataset

GEO [GSE205335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205335)
(Ahn / Lee, *eLife* 2024, palliative ICI biopsy/effusion scRNA). Open
processed UMI `dgCMatrix` + author identity + GEO SOFT. Controlled EGA
raw (`EGAD00001008703`) was not accessed.

| Item | Choice |
| --- | --- |
| Compartment | `lineage.sub == Malignant cells` on non-normal tissues |
| Dropped | Normal Lung / LN / Brain; four normal-only patients (0 malignant) |
| Unit | patient (biopsies from the same patient pooled) |
| Cap | ≤400 malignant cells / patient after catalog (honest n reports both) |
| Score | CLDN4 continuous + tertiles on the malignant object |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| IFN | compact Hallmark-like ISG panel (locked in `scripts/gene_sets.py`) |
| Unused | TACSTD2∩CLDN4 dual-high; MPR (unlabeled); GSE131907; GSE207422 |

MPR/NMPR is not labelled. RECIST is not substituted for MPR.
RECIST R = PR; NR = SD+PD; NE dropped from R vs NR.

## Trajectory clock

Slingshot (Street et al. 2018) is installed (Bioconductor, including
`DelayedMatrixStats`) and run on Harmony-PCA (Leiden clusters,
`start.clus` = lowest-CLDN4 eligible cluster).
PAGA (Wolf et al. 2019) is the cluster graph. Scanpy DPT (Haghverdi
et al. 2016) is the companion clock, same root. Direction is an
external arrow (lowest-CLDN4 Leiden, never the CLDN4-high cluster),
not a UMAP.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → Harmony (`batch = patient_id`) → k-NN 30 |
| Clusters | Leiden 0.6; PAGA on Leiden |
| Root | Lowest-CLDN4 Leiden with ≥40 cells and ≥3 patients. Never CLDN4-high. |
| Slingshot | PCA (Harmony) for pseudotime; UMAP principal curves for the lineage plot only |
| Inferential n | patient. Cell-level ρ is descriptive. |
| Extra figures | Within-patient CLDN4-high vs low program scores (min 8 cells/arm); patient CLDN4/barrier/IFN vs DPT |

## Tests

Primary: patient-level Spearman of mean CLDN4 vs mean DPT / Slingshot PT /
barrier (CLDN4 excluded) / IFN, plus barrier vs DPT and IFN vs DPT.
BH inside that list only.

RECIST (not BH): Mann–Whitney of patient-mean DPT, Slingshot PT, CLDN4,
barrier, and IFN in R vs NR.

Sensitivity (not BH): ADC+SQ only; drop patients with <50 analysis cells;
CLDN4 %pos vs DPT.

Harmony on patient removes a main-effect batch. It does not prove a
within-tumor differentiation axis. Author malignant is not CNV.
