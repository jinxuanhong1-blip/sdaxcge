# Methods — triple-that-differs Slingshot/PAGA, CLDN4 only

ADDITIVE. GSE123902 + GSE131907 + GSE205335 epithelium.
Does **not** add GSE148071. This is **not** the 7-pool.
TACSTD2 does not define groups. No dual-high gate.
The PR #459 T/NK Spearman is not re-audited.

**Question.** On the triple that differs, where does **CLDN4** sit on an
AT2-rooted trajectory relative to AT2 and a CLDN4-excluded barrier/keratin
program? Inferential unit = sample / patient / donor-sample.

## Datasets (public processed UMI only)

| Cohort | Citation | Matrix | Epithelium kept | Dropped |
| --- | --- | --- | --- | --- |
| GSE123902 | Laughney et al., *Nat Med* 2020, PMID 31959988 | RAW.tar dense UMI CSVs (90.4 MB) | marker epithelium `(EPCAM\|KRT8\|KRT18\|KRT19)>0`; NORMAL kept as units for an AT2-like root | 36.5 GB annotated H5; SRA FASTQ |
| GSE131907 | Kim et al., *Nat Commun* 2020, PMID 32385277 | raw UMI text (~0.39 GB) + author annotation | `Cell_type == Epithelial cells` and `Sample_Origin` ∈ {tLung, nLung, tL/B, mLN, mBrain} | PE unlabeled epithelium; 2.86 GB log2TPM; EGA FASTQ |
| GSE205335 | Ahn / Lee, *eLife* 2024 | UMI dgCMatrix + author identity + SOFT | `lineage.total == Epithelial cells` on non-normal tissues | Normal Lung / LN / Brain; EGA FASTQ |

## Trajectory clock

Slingshot (Street et al. 2018) is the requested method. Tools are installed
in this run (R 4.3 + slingshot 2.10.0). On this ~24k-cell object, lineages
are fit on a stratified subsample (≤80 cells / Leiden cluster,
`approx_points=150`) and projected to all cells by 5-NN in Harmony-PCA.
If `Rscript` / `slingshot` is still missing after install, the
**documented fallback** is scanpy diffusion pseudotime (Haghverdi et al. 2016).
Root is **never CLDN4-high**. Preferred root: GSE131907 nLung author AT2
(median AT2 score among not-CLDN4-high). Palantir is optional.

| Item | Choice |
| --- | --- |
| Graph | Seurat-v3 HVG 3000 → PCA → Harmony (`batch = dataset`) → k-NN 30 |
| Harmony death | Per-dataset graphs (no Harmony), then **stack** patient CLDN4 vs that dataset's pseudotime |
| Clusters | Leiden 0.6; PAGA on Leiden (geometry only) |
| Root | Never CLDN4-high. Prefer GSE131907 nLung author AT2. Else GSE123902 NORMAL. Else max-AT2 Leiden. |
| Score | CLDN4 continuous + tertiles on the epithelial object |
| Barrier/keratin | KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**no CLDN4**) |
| Inferential n | GSE123902 `donor:tissue` + GSE131907 `Sample` + GSE205335 `patient` |
| Cap | ≤350 cells / unit after protecting nLung AT2 / NORMAL SFTPC+ |
| Extra figure | Within-unit CLDN4-high vs low program scores (min 8 cells/arm) |
| Unused | GSE148071; 7-pool; dual-high TACSTD2∩CLDN4; ICI / MPR / RECIST as a trajectory label |

## Tests

Primary: stacked sample-level Spearman of mean CLDN4 vs mean pseudotime / AT2 /
barrier (CLDN4 excluded). BH inside that list only.

Sensitivity (not BH): per-cohort, tLung/PRIMARY-only, nLung/NORMAL-only,
tumor-only (drop nLung/NORMAL), GSE205335 ADC+SQ.

The stacked pseudotime correlation mixes cohorts and normal vs tumor. It is
**not** a within-tumor progression test.

Done when `results/tables/stacked_sample_cldn4_pseudotime.tsv` exists.
