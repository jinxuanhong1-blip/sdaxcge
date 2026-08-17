# Methods — winning-pair Slingshot/Palantir trajectory, CLDN4 only

ADDITIVE. Does **not** rewrite `methods/scrna_paga_cldn4/` (PR #325, GSE131907-only PAGA).
GSE207422 is not added. TACSTD2 does not define groups. No dual-high gate.

**Question.** On the winning public pair GSE131907+GSE205335 **epithelium only**, where does **CLDN4** sit on an AT2-rooted trajectory relative to AT2 and a CLDN4-excluded barrier/keratin program?

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix | Epithelium kept | Dropped |
| --- | --- | --- | --- | --- |
| GSE131907 | Kim et al., *Nat Commun* 2020, PMID 32385277 | raw UMI text (~0.39 GB) + author annotation | `Cell_type == Epithelial cells` and `Sample_Origin` ∈ {tLung, nLung, tL/B, mLN, mBrain} | PE unlabeled epithelium; 2.86 GB log2TPM; EGA FASTQ |
| GSE205335 | Ahn / Lee, *eLife* 2024 | UMI dgCMatrix + author identity + SOFT | `lineage.total == Epithelial cells` on non-normal tissues | Normal Lung / LN / Brain; EGA FASTQ |

Winning-pair membership is taken as given from PR #290 (author %pos GSE131907+GSE205335 vs T/NK). This folder asks a **different** question (trajectory on epithelium), not a T/NK re-rank.

## Trajectory clock

Slingshot (Street et al. 2018) is the requested method. If `Rscript` / `slingshot` is missing, the **documented fallback** is scanpy diffusion pseudotime (Haghverdi et al. 2016) rooted on GSE131907 nLung author AT2 (median AT2 score). Palantir is optional and is not required for the primary table. Direction is an external arrow, not a UMAP.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → Harmony (`batch = dataset`) → k-NN 30 |
| Clusters | Leiden 0.6; PAGA on Leiden (geometry only) |
| Root | GSE131907 nLung author AT2. Never CLDN4-high. |
| Score | CLDN4 continuous + tertiles on the epithelial object |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| Inferential n | GSE131907 `Sample` + GSE205335 `patient`. Cell-level ρ is descriptive. |
| Cap | ≤350 cells / unit after protecting nLung AT2 (honest n reports catalog vs analysis) |
| Extra figure | Within-unit CLDN4-high vs low program scores (min 8 cells/arm) |
| Unused | GSE207422; dual-high TACSTD2∩CLDN4; ICI / MPR / RECIST as a trajectory label |

## Tests

Primary: sample-level Spearman of mean CLDN4 vs mean DPT / AT2 / barrier (CLDN4 excluded). BH inside that list only.

Sensitivity (not BH): per-cohort, tLung-only, nLung-only, tumor-only (drop nLung), GSE205335 ADC+SQ.

The pooled DPT correlation mixes cohorts and nLung vs tumor. It is **not** a within-tumor progression test.
