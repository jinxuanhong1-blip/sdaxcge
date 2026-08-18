# FINDING — GSE165641 KL GEMM scRNA: Cldn4-only (Seurat)

**ADDITIVE public MOUSE. Cldn4-only. Thesis already correct. No private 8-KL. No Python-only primary.**

Public processed matrix **exists**. Seurat `CreateSeuratObject` was run on the two Cell Ranger filtered count matrices. Epithelium and T/NK both exist. Honest **n = 2 KL mice**. Spearman and high-vs-low are **not reported** (n=2 forces |r|=1; a 1-vs-1 split is not a test).

Thesis is taken as given: Cldn4-high epithelium is barrier / immune-cold. This folder does **not** audit that thesis as claim-failed. It scores **one** public KL mixed-lineage 10x series at mouse level.

---

## Decision

| Question | Answer |
|---|---|
| Processed matrix on GEO? | **Yes** — Cell Ranger MTX/H5 per mouse + series `GSE165641_processed_normalized_matrix_data.Rdata.gz` (31,053 × 7,180) |
| `CreateSeuratObject` used? | **Yes** — R 4.3.3, Seurat 5.5.1, SeuratObject 5.4.0 |
| Cldn4 present? | **Yes** — `ENSMUSG00000047501` / `Cldn4`. All locked T/NK, IFN, and MHC genes present |
| Epithelial + T/NK compartments? | **Yes** — epithelium clusters 3/10/13; T/NK cluster 7 (98% T/NK-marker) |
| Mouse-level Cldn4 vs T/NK? | **Yes, two values. No Spearman.** Direction: higher epithelial Cldn4, lower T/NK fraction |
| Mouse-level epithelial IFN/MHC? | **Yes.** Combined IFN/MHC and MHC are lower in the Cldn4-higher mouse; IFN alone is not |
| Honest n | **2 KL mice** (not 7,180 cells; not 6,693 QC cells; not private 8 KL) |
| High vs low? | **No-go** — 1 vs 1 is not a split |
| Join human concordant-4? | **No** — n=2; additive only |
| Dual-high TACSTD2 × CLDN4 | **not defined** |

---

## What the public objects actually are

Wang P, Zhong B et al., *Advanced Science* 2021 ([PMID 34369094](https://pubmed.ncbi.nlm.nih.gov/34369094/); [GSE165641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165641)). Title: 10X of tumor sections from tumor-bearing **KrasG12D/+ Lkb1fl/fl (KL)** mice 10 weeks after Ad-Cre. Overall design: RBC-depleted cells from **2** tumor-bearing KL mice. GEO: targeted recovery ~5,000 cells (**4,400 and 2,780**). Authors processed UMI counts in Seurat 3.2.0; genome Ensembl 94.

| Public object | Size | Used here |
|---|---|---|
| `GSM5047302_KL1_count.tar.gz` | 91.6 Mb | **yes — primary** Cell Ranger `filtered_feature_bc_matrix` (4,400 barcodes) |
| `GSM5047303_KL2_count.tar.gz` | 87.1 Mb | **yes — primary** (2,780 barcodes) |
| `GSE165641_processed_normalized_matrix_data.Rdata.gz` | 39.6 Mb | **inspected only** — dense `processed_data` data.frame, 31,053 × 7,180, columns `KL1_*` / `KL2_*`. Not the `CreateSeuratObject` input |
| `GSE165641_RAW.tar` | 178.6 Mb | no (same MTX nested) |
| SRA FASTQ | — | **no** |

GSM5047303 SOFT text says `Lkb2fl/fl`; treated as a typo for **Lkb1** (series title, KL1 characteristics, and PMID).

No KP / K arm. No private 8-KL matrix.

Cell Ranger metrics (depositor):

| mouse | GSM | estimated cells | median genes | median UMI |
|---|---|---:|---:|---:|
| KL1 | GSM5047302 | 4,400 | 1,206 | 3,115 |
| KL2 | GSM5047303 | 2,780 | 2,425 | 7,193 |

KL2 is deeper per cell. That is a technical difference, not a genotype contrast.

---

## Methods (locked)

- **R + Seurat** (`Read10X` → `CreateSeuratObject` on each mouse → `merge` → `JoinLayers`). No Python primary.
- **Cldn4 only.** Tacstd2 is not a gate. Other claudins are not scored.
- Mouse unit = GEO replicate (`KL1` = GSM5047302, `KL2` = GSM5047303).
- QC after Cell Ranger filter: `nFeature_RNA ≥ 200`, `nCount_RNA ≥ 500`, `percent.mt < 25` (`^mt-`). 7,180 → **6,693** cells.
- `NormalizeData` (LogNormalize, 10,000), VST 2,000 HVG, `ScaleData`, PCA 30, neighbors/clusters (res 0.4, dims 1:20), UMAP (uwot).
- **Compartments:** `AddModuleScore` on canonical mouse markers; Epcam+ Ptprc− forced to epithelium; Cd3d/e/g or Nkg7/Ncr1/Klrb1c/Cd8a and Ptprc+ Epcam− forced to T/NK.
- **T/NK fraction** = T/NK cells / QC cells per mouse (mixed digest — a real fraction, not FACS leftover).
- **Epithelial Cldn4** = mean log-norm `Cldn4` and % UMI>0 inside epithelium.
- **Epithelial IFN / MHC** = mean log-norm of present genes in the locked lists.
- Spearman requires n≥4. Q4 vs Q1 requires n≥8. **Neither fires.** The two-mouse direction is descriptive only.
- No SRA. No GSE180963 merge. No GSE6135 re-open. No private 8 KL.
- Reproduce: `Rscript methods/seurat_gse165641_cldn4/analyze_gse165641.R /tmp/geo_gse165641 methods/seurat_gse165641_cldn4`

**T/NK genes:** Cd3d, Cd3e, Cd3g, Cd2, Cd8a, Cd8b1, Cd4, Nkg7, Gzma, Gzmb, Prf1, Klrb1c, Ncr1, Klrd1, Klrc1, Ifng.

**IFN:** Stat1, Stat2, Irf1, Irf7, Irf9, Isg15, Ifit1, Ifit2, Ifit3, Mx1, Oasl2, Rsad2, Ifih1, Ddx58, Ifnb1.

**MHC:** B2m, H2-K1, H2-D1, H2-Q4, H2-Q6, H2-Q7, H2-Aa, H2-Ab1, H2-Eb1, Tap1, Tap2, Psmb8, Psmb9, Nlrc5, Ciita.

All 1 + 16 + 15 + 15 genes are present (`tables/genes_used.tsv`).

---

## Results (honest n = 2)

Do **not** quote 6,693 or 7,180 as n.

### Compartments exist

Cldn4 is epithelial-restricted (sanity, not a claim):

| compartment | mean log-norm Cldn4 | % Cldn4+ |
|---|---:|---:|
| epithelial | 0.614 | 41.4 |
| T/NK | 0.010 | 0.6 |
| myeloid | 0.007 | 1.1 |
| neutrophil | 0.008 | 0.6 |
| B / endo / fibro / other | ≤0.007 | ≤0.4 |

T/NK is a discrete cluster (cluster 7: 269 cells, 98% T/NK-marker, 96% Ptprc). Epithelium is clusters 3 + 10 + 13 (Cldn4 %pos 37 / 56 / 36). UMAP: `figures/umap_compartment.png`, `figures/umap_cldn4.png`.

### Mouse-level scores

From `tables/mouse_level_scores.tsv`.

| mouse | GSM | n_QC | n_epi | n_T/NK | T/NK fraction | epi Cldn4 mean | epi Cldn4 %pos | epi IFN | epi MHC | epi IFN/MHC |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| KL1 | GSM5047302 | 4,199 | 358 | 94 | **0.022** | **0.781** | 0.455 | 0.114 | **0.350** | **0.232** |
| KL2 | GSM5047303 | 2,494 | 260 | 242 | **0.097** | **0.384** | 0.358 | 0.065 | **0.618** | **0.342** |

Descriptive two-point direction (not a test):

- The Cldn4-higher mouse (KL1) has the **lower T/NK fraction** (2.2% vs 9.7%) and the **lower epithelial MHC / IFN/MHC**.
- Epithelial IFN alone is the **opposite** sign (higher in KL1).
- Whole-digest T/NK gene score tracks the fraction (0.011 vs 0.075).

That is the thesis-consistent direction for Cldn4 vs T/NK and vs MHC. It is **two points**. Do not write ρ, p, or a KL Cldn4–T/NK law.

### Composition confound (why n=2 cannot carry this)

| mouse | neutrophil fraction | T/NK among Ptprc+ |
|---|---:|---:|
| KL1 | **0.602** (2,528 / 4,199) | 0.028 |
| KL2 | **0.057** (142 / 2,494) | 0.152 |

KL1 is a neutrophil-heavy digest with lower per-cell complexity (Cell Ranger median genes 1,206 vs 2,425). T/NK fraction is a valid mixed-tumor readout, but the two mice are not exchangeable tissue samples. The T/NK difference survives restriction to Ptprc+ cells, so it is not only neutrophil dilution — still n=2.

---

## How to read this

- This is **mixed-lineage 10x** of KL tumor sections, not FACS epithelium. Both compartments exist; Cldn4 sits in epithelium.
- Honest n is **2 mice**. Still score; do not overclaim.
- Do not write a Cldn4–T/NK or Cldn4–IFN/MHC law from two points.
- Do not treat this as a replacement for GSE6135 (public KL bulk, larger mouse n) or the private 8 KL.
- IFN and MHC should not be collapsed into one sentence: MHC/IFN/MHC follow the cold direction; IFN alone does not.
- Thesis is unchanged. Additive: public KL scRNA with Cldn4 + T/NK + epithelium exists; n=2 is the ceiling.

## Files

| path | what |
|---|---|
| `analyze_gse165641.R` | Seurat primary (`CreateSeuratObject`) |
| `tables/mouse_level_scores.tsv` | one row per mouse |
| `tables/honest_n.tsv` | what can / cannot be tested |
| `tables/compartment_counts.tsv` | cell counts by mouse × compartment (not n) |
| `tables/cluster_audit.tsv` | cluster markers |
| `tables/genes_used.tsv` | Cldn4 / T/NK / IFN / MHC present? |
| `tables/cell_metadata.tsv` | per-cell scores (do not quote as n) |
| `figures/umap_compartment.png` | UMAP split by mouse |
| `figures/umap_cldn4.png` | Cldn4 feature plot |
| `figures/mouse_two_point.png` | two-point plot, labeled **not a Spearman** |
