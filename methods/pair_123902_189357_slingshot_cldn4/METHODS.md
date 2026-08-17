# Methods — REAL Slingshot/PAGA on GSE123902+GSE189357, CLDN4 only

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2∩CLDN4. Does **not** redo the winning-pair
Slingshot/DPT (PR #449) or GSE131907-only PAGA (PR #325).

**Question.** On the PR #459 pair that **differs** (GSE123902+GSE189357 %pos, n=22,
ρ=−0.638 vs T/NK), where do **CLDN4**, a CLDN4-excluded **barrier/keratin** score, and
an **IFN** score sit on a **real Slingshot** lineage? PAGA is geometry, not the clock.

The given T/NK Spearman is **not re-audited**. Q4 vs Q1 tails **7/5 are thin** — said so.

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix | Cells kept | Dropped |
| --- | --- | --- | --- | --- |
| GSE123902 | Laughney et al., *Nat Med* 2020 | dense CSV in `GSE123902_RAW.tar` (~90 MB) | marker-epithelial on tumor/met (given units) + NORMAL (root pool only) | FASTQ; ineligible libraries |
| GSE189357 | Zhu/Wang AIS–IAC atlas | 10x MTX in `GSE189357_RAW.tar` (~624 MB) | marker-epithelial, one tumor sample per patient (TD1–TD9) | FASTQ |

Marker-epithelial gate (same as PR #459): `(EPCAM\|KRT8\|KRT18\|KRT19)>0` and `PTPRC==0`.
Patient/donor is the unit. Cap ≤350 cells / unit after protecting NORMAL SFTPC+ cells.

## Trajectory clock

**Slingshot** (Street et al. 2018; Bioconductor) is required. Tools are installed
(`scripts/install_r_slingshot.sh`). This is **not** a DPT-only fallback.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → Harmony (`batch = dataset`) → k-NN 30 |
| Clusters | Leiden 0.6; **PAGA** on Leiden (geometry only) |
| Root | GSE123902 NORMAL, **not CLDN4-high**, median AT2 score. Never a CLDN4-high cell. |
| Clock | Slingshot on Harmony PCs, `start.clus` = root Leiden |
| Scores | CLDN4 continuous + tertiles; barrier/keratin **without CLDN4**; Hallmark IFNα mean |
| Inferential n | given 22 tumor/met units. Cell-level ρ is descriptive. |
| Extra figures | CLDN4+barrier+IFN along PT; within-unit high vs low; honest n |
| Unused | dual-high; CellChat; ICI/MPR as a trajectory label; 7-cohort pool |

## Tests

Primary: given-unit Spearman of mean CLDN4 / barrier / IFN vs mean Slingshot PT. BH inside that list only.

Sensitivity (not BH): per-cohort.

Tails may be thin. Q4 vs Q1 7/5 is **not** the Slingshot n.

SFTPC/SFTPA1 and the club panel were absent from the GSE123902∩GSE189357 gene intersection
on this extract. The AT2 score is therefore incomplete. Report that, do not hide it.
