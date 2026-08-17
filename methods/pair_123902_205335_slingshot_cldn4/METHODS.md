# Methods — pair GSE123902+GSE205335 Palantir/PAGA, CLDN4 only

ADDITIVE. Does **not** re-audit PR #459 T/NK infiltrate. Does **not** re-run
PR #473 malignant IFN DE. GSE148071 is not added. TACSTD2 does not define
groups. No dual-high gate.

**Question.** On the PR #459 pair that differs (GSE123902+GSE205335), where do
**CLDN4**, an IFN program, and a CLDN4-excluded barrier/keratin program sit on
a real Palantir trajectory whose root is **not** CLDN4-high?

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix | Epithelium kept | Dropped |
| --- | --- | --- | --- | --- |
| GSE123902 | Laughney et al., *Nat Med* 2020, PMID 32066974 | dense CSV in `GSE123902_RAW.tar` (~90 MB) | marker-epi: (EPCAM\|KRT8\|KRT18\|KRT19)>0 and PTPRC==0; locked tumor/met donors (PRIMARY preferred) plus NORMAL libraries for the root | unused donors; FASTQ |
| GSE205335 | Ahn / Lee, *eLife* 2024 | UMI dgCMatrix + author identity + SOFT | `lineage.total == Epithelial cells` on non-normal tissues | Normal Lung / LN / Brain; EGA FASTQ |

Pair membership is taken as given from PR #459. This folder asks a **different**
question (trajectory on epithelium), not a T/NK re-rank.

## Trajectory clock

Slingshot R is missing (`Rscript` not on PATH). **Palantir** (Setty et al. 2019)
is installed and is the primary clock. PAGA on Leiden is geometry only. scanpy
DPT is a companion, not the claim.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → Harmony (`batch = dataset`) → k-NN 30 |
| Clusters | Leiden 0.6; PAGA on Leiden |
| Root | GSE123902 NORMAL AT2-like (median AT2 in the top AT2 tercile). **Never CLDN4-high.** |
| Score | CLDN4 continuous + tertiles; IFN = Hallmark IFNα ∪ IFNγ; barrier/keratin **no CLDN4** |
| Inferential n | GSE123902 `donor` (tumor library) + GSE205335 `patient`. NORMAL libraries are root material, not primary units. |
| Cap | ≤350 cells / unit after protecting NORMAL SFTPC+ |
| Extra figure | CLDN4 + IFN + barrier along Palantir PT; within-unit CLDN4-high vs low |
| Unused | GSE148071; dual-high TACSTD2∩CLDN4; ICI / MPR / RECIST as a trajectory label |

## Tests

Primary: sample-level Spearman of mean CLDN4 / IFN / barrier vs mean Palantir PT,
plus CLDN4 vs IFN and CLDN4 vs barrier. BH inside that list only.

Sensitivity (not BH): per-cohort, include-NORMAL units, GSE205335 ADC+SQ, DPT companion.

The pooled PT correlation mixes cohorts. It is **not** a within-tumor progression test.
