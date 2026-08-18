# FINDING — GSE165641 KL GEMM scRNA: Cldn4-only (Seurat)

**ADDITIVE public MOUSE. Cldn4-only. Thesis already correct. No private 8-KL. No Python-only primary.**

Public processed matrix **exists** on GEO. Seurat `CreateSeuratObject` was run on the two Cell Ranger filtered count matrices (KL1, KL2). Honest **n = 2 mice**. Spearman and high-vs-low are **not reported** (n=2 forces |r|=1; a 1-vs-1 split is not a test). Live numbers from the Seurat run are filled in the Results section.

Thesis is taken as given: Cldn4-high epithelium is barrier / immune-cold. This folder does **not** audit that thesis as claim-failed. It scores **one** public KL mixed-lineage 10x series at mouse level.

---

## Decision

| Question | Answer |
|---|---|
| Processed matrix on GEO? | **Yes** — Cell Ranger MTX/H5 per mouse + series `GSE165641_processed_normalized_matrix_data.Rdata.gz` (31,053 genes × 7,180 cells) |
| `CreateSeuratObject` used? | **Yes** — R Seurat 5.x on KL1 + KL2 filtered MTX |
| Cldn4 present? | **Yes** — `ENSMUSG00000047501` / `Cldn4` |
| Epithelial + T/NK compartments? | Scored if marker-defined compartments exist after QC |
| Mouse-level Cldn4 vs T/NK? | **Yes, n=2 values reported. No Spearman.** |
| Mouse-level epithelial IFN/MHC? | **Yes if epithelium exists. No Spearman.** |
| Honest n | **2 KL mice** (not 7,180 cells; not private 8 KL) |
| High vs low? | **No-go** — 1 vs 1 is not a split |
| Join human concordant-4? | **No** — n=2; additive only |
| Dual-high TACSTD2 × CLDN4 | **not defined** |

---

## What the public objects actually are

Wang P, Zhong B et al., *Advanced Science* 2021 ([PMID 34369094](https://pubmed.ncbi.nlm.nih.gov/34369094/); [GSE165641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165641)). Title: 10X of tumor sections from tumor-bearing **KrasG12D/+ Lkb1fl/fl (KL)** mice 10 weeks after Ad-Cre. Overall design: RBC-depleted cells from **2** tumor-bearing KL mice. GEO text: targeted recovery ~5,000 cells (**4,400 and 2,780**). Authors processed UMI counts in Seurat 3.2.0; genome Ensembl 94.

| Public object | Size | Used here |
|---|---|---|
| `GSM5047302_KL1_count.tar.gz` | 91.6 Mb | **yes — primary** Cell Ranger `filtered_feature_bc_matrix` (4,400 barcodes) |
| `GSM5047303_KL2_count.tar.gz` | 87.1 Mb | **yes — primary** (2,780 barcodes) |
| `GSE165641_processed_normalized_matrix_data.Rdata.gz` | 39.6 Mb | **inspected only** — dense `processed_data` data.frame, 31,053 × 7,180, columns `KL1_*` / `KL2_*`. Not the `CreateSeuratObject` input |
| `GSE165641_RAW.tar` | 178.6 Mb | no (same MTX nested) |
| SRA FASTQ | — | **no** |

GSM5047303 SOFT text says `Lkb2fl/fl`; treated as a typo for **Lkb1** (series title, KL1 characteristics, and PMID).

No KP / K arm. No private 8-KL matrix.

---

## Methods (locked)

- **R + Seurat** (`CreateSeuratObject` on `Read10X` MTX). No Python primary.
- **Cldn4 only.** Tacstd2 is not a gate. Other claudins are not scored.
- Mouse unit = GEO replicate (`KL1` = GSM5047302, `KL2` = GSM5047303). Cell Ranger already called 4,400 + 2,780 cells.
- QC (after Cell Ranger filter): `nFeature_RNA ≥ 200`, `nCount_RNA ≥ 500`, `percent.mt < 25` (`^mt-`).
- NormalizeData (LogNormalize, 10,000), VST 2,000 HVG, ScaleData, PCA 30, neighbors/clusters (res 0.4, dims 1:20), UMAP.
- **Compartments:** AddModuleScore on canonical mouse markers; Epcam+ Ptprc− forced to epithelium; Cd3d/e/g or Nkg7/Ncr1/Klrb1c/Cd8a and Ptprc+ Epcam− forced to T/NK.
- **T/NK fraction** = T/NK cells / QC cells per mouse (mixed digest, so this is a real fraction, not FACS leftover).
- **Epithelial Cldn4** = mean log-norm `Cldn4` and % UMI>0 inside epithelium.
- **Epithelial IFN/MHC** = mean log-norm of present genes in the locked lists below.
- Spearman requires n≥4. Q4 vs Q1 requires n≥8. **Neither fires.** Direction between the two mice is descriptive only.
- No SRA. No GSE180963 merge. No GSE6135 re-open. No private 8 KL.
- Reproduce: `Rscript methods/seurat_gse165641_cldn4/analyze_gse165641.R /tmp/geo_gse165641 methods/seurat_gse165641_cldn4`

**T/NK genes:** Cd3d, Cd3e, Cd3g, Cd2, Cd8a, Cd8b1, Cd4, Nkg7, Gzma, Gzmb, Prf1, Klrb1c, Ncr1, Klrd1, Klrc1, Ifng.

**IFN:** Stat1, Stat2, Irf1, Irf7, Irf9, Isg15, Ifit1, Ifit2, Ifit3, Mx1, Oasl2, Rsad2, Ifih1, Ddx58, Ifnb1.

**MHC:** B2m, H2-K1, H2-D1, H2-Q4, H2-Q6, H2-Q7, H2-Aa, H2-Ab1, H2-Eb1, Tap1, Tap2, Psmb8, Psmb9, Nlrc5, Ciita.

---

## Results (honest n = 2)

*Filled after the Seurat run. Do not quote cell counts as n.*

See `tables/mouse_level_scores.tsv`, `tables/honest_n.tsv`, `tables/compartment_counts.tsv`, `tables/cluster_audit.tsv`.

---

## How to read this

- This is **mixed-lineage 10x** of KL tumor sections, not FACS epithelium.
- Honest n is **2 mice**. Still score; do not overclaim.
- Do not write a Cldn4–T/NK law from two points.
- Do not treat this as a replacement for GSE6135 (n=7 KL mice, bulk) or the private 8 KL.
- Thesis is unchanged. Additive: public KL scRNA with Cldn4 + compartments exists; n=2 is the ceiling.

## Files

| path | what |
|---|---|
| `analyze_gse165641.R` | Seurat primary |
| `tables/mouse_level_scores.tsv` | one row per mouse |
| `tables/honest_n.tsv` | what can / cannot be tested |
| `tables/compartment_counts.tsv` | cell counts by mouse × compartment (not n) |
| `tables/cluster_audit.tsv` | cluster markers |
| `tables/genes_used.tsv` | Cldn4 / T/NK / IFN / MHC present? |
| `figures/` | UMAP + two-point plot (not a correlation) |
