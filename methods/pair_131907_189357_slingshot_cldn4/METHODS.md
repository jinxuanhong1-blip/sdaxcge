# Methods — pair GSE131907+GSE189357 REAL Slingshot/PAGA, CLDN4 only

ADDITIVE. Does **not** rewrite `methods/scrna_paga_cldn4/` (PR #325) or
`methods/winpair_131907_205335_slingshot_cldn4/` (PR #449).
GSE148071 is not added. TACSTD2 does not define groups. No dual-high gate.

**Question.** On the PR #459 pair that already differs (GSE131907+GSE189357
malignant CLDN4 %pos vs T/NK, n=30, ρ=−0.542; **not re-audited here**),
where do **CLDN4**, a CLDN4-excluded barrier/keratin program, and an IFN
program sit on a real Slingshot lineage rooted off CLDN4-high?

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix | Epithelium kept | Dropped |
| --- | --- | --- | --- | --- |
| GSE131907 | Kim et al., *Nat Commun* 2020, PMID 32385277 | raw UMI text (~0.39 GB) + author annotation | `Cell_type == Epithelial cells` and `Sample_Origin` ∈ {tLung, nLung, tL/B, mLN, mBrain} | PE unlabeled epithelium; 2.86 GB log2TPM; EGA FASTQ |
| GSE189357 | Zhu / Fan / Jiang, GEO GSE189357, PMID 36434043 | 10x MTX in `GSE189357_RAW.tar` (~624 MB) | marker epithelium: `(EPCAM\|KRT8\|KRT18\|KRT19)>0` and `PTPRC==0` | FASTQ / SRA |

GSE189357 histology from the series matrix (field typo `histolgical type`):
TD1 IAC, TD2 IAC, TD3 MIA, TD4 MIA, TD5 AIS, TD6 MIA, TD7 AIS, TD8 AIS, TD9 IAC.

## Trajectory clock

**Real Slingshot** (Street et al. 2018, Bioconductor `slingshot`) is required.
PAGA (Wolf et al. 2019) is the Leiden graph. Diffusion pseudotime is a companion
only. Direction is an external arrow, not a UMAP.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → Harmony (`batch = dataset`) → k-NN 30 |
| Clusters | Leiden 0.6; PAGA on Leiden (geometry only) |
| Root / start | GSE131907 nLung author AT2. Never CLDN4-high. If that cell's Leiden is the highest-CLDN4 cluster, start is reassigned to the next AT2-bearing cluster. |
| Score | CLDN4 continuous + tertiles; barrier/keratin without CLDN4; Hallmark IFNα ∪ IFNγ mean |
| Inferential n | GSE131907 `Sample` + GSE189357 `patient`. Cell-level ρ is descriptive. |
| Cap | ≤350 cells / unit after protecting nLung AT2 |
| Extra figures | CLDN4 + barrier + IFN along Slingshot PT; within-unit CLDN4-high vs low |
| Unused | GSE148071; dual-high TACSTD2∩CLDN4; T/NK re-rank of PR #459 |

## Tests

Primary: sample-level Spearman of mean CLDN4 / barrier / IFN vs mean Slingshot
pseudotime (primary lineage). BH inside that list only.

Sensitivity (not BH): per-cohort, tLung-only, nLung-only, tumor-only (drop nLung),
GSE189357 AIS/MIA/IAC, companion DPT.

The pooled Slingshot correlation mixes cohorts and nLung vs tumor. It is **not**
a within-tumor progression test.

**Done when** `results/tables/lineage_table.tsv` (copy of `slingshot_lineages.tsv`) exists.
