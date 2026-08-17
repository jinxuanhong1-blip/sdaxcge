# Methods — pair GSE123902+GSE131907 REAL Slingshot/PAGA, CLDN4 only

ADDITIVE. Does **not** rewrite `methods/scrna_paga_cldn4/` (PR #325) or
`methods/winpair_131907_205335_slingshot_cldn4/` (PR #449).
GSE148071 is not added. TACSTD2 does not define groups. No dual-high gate.

**Question.** On the PR #459 pair that **differs** (GSE123902+GSE131907),
where do **CLDN4**, a CLDN4-excluded **barrier/keratin** program, and a compact
**IFN** program sit on a real Slingshot lineage?

PR #459 T/NK ρ is taken as given and is **not re-audited**.

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix | Epithelium kept | Unit | Dropped |
| --- | --- | --- | --- | --- | --- |
| GSE123902 | Laughney et al., *Nat Med* 2020 (human arm of GSE123904) | GEO RAW.tar dense CSVs (~90 MB) | marker epithelium on PRIMARY_TUMOUR / METASTASIS: (EPCAM\|KRT8\|KRT18\|KRT19)>0 and PTPRC==0 | **donor** | NORMAL samples; 36.5 GB author H5 |
| GSE131907 | Kim et al., *Nat Commun* 2020, PMID 32385277 | raw UMI text (~0.39 GB) + author annotation | `Cell_type == Epithelial cells` and `Sample_Origin` ∈ {tLung, nLung, tL/B, mLN, mBrain} | **sample** | PE unlabeled epithelium; 2.86 GB log2TPM; EGA FASTQ |

GSE148071 is not added. Cap ≤350 cells / unit after protecting GSE131907 nLung AT2.

## Trajectory clock

**Slingshot** (Street et al. 2018) is required. This folder does **not** stop
if R is missing: install R + Bioconductor `slingshot` first
(`scripts/install_r_slingshot.sh`). DPT is a companion ordering only.
PAGA (Wolf et al. 2019) is geometry on Leiden vertices.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → Harmony (`batch = dataset`) → k-NN 30 |
| Clusters | Leiden 0.6; PAGA on Leiden |
| Root | GSE131907 nLung author AT2. **Never CLDN4-high** (start cluster refused if it ranks in the top 15% of cluster-mean CLDN4). |
| Score | CLDN4 continuous + tertiles; barrier/keratin **without CLDN4**; compact IFN/ISG core **without CLDN4** |
| Inferential n | GSE123902 `donor` + GSE131907 `Sample`. Cell-level ρ is descriptive. |
| Extra figures | Programs along Slingshot PT; within-unit CLDN4-high vs low |
| Unused | GSE148071; dual-high TACSTD2∩CLDN4; ICI / MPR / RECIST as a trajectory label |

## Tests

Primary: sample-level Spearman of mean CLDN4 / barrier / IFN vs mean Slingshot
pseudotime, plus CLDN4 vs barrier and CLDN4 vs IFN. BH inside that list only.

Sensitivity (not BH): per-cohort, tLung-only, nLung-only, tumor-only (drop nLung),
DPT companion.

The pooled Slingshot correlation mixes cohorts and nLung vs tumor. It is **not**
a within-tumor progression test.
