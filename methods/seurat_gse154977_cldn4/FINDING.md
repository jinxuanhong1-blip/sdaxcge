# FINDING — GSE154977 KP 30w 10x: Cldn4-only (Seurat)

**ADDITIVE public MOUSE. Cldn4-only. Thesis already correct. No dual-high.**

Primary analysis is **R + Seurat 5.5.1** (not Python). Public processed raw-count matrix **exists** on GEO. Epithelial Cldn4 and IFN / MHC / TJ were scored at **biological mouse** level. **T/NK fraction cannot be scored.** These are tumor-only 10x libraries: live `tdTomato+ / CD45− / CD11b− / TER119− / CD31−` (Marjanovic et al., STAR Methods). Residual `Cd3d` / `Nkg7` / `Ptprc` is leak, not a compartment. T/NK Spearman is **empty**. This accession does **not** join the human concordant-4 pool.

Primary KP (n=4 mice): IFN and MHC-I/APM are **not** down in Cldn4-high. Spearman is at the lock floor (n=4). Exact p for ρ=1 is 0.083; asymptotic p=0 is invalid here and is **not** quoted. Q4−Q1 is locked off (needs n_units≥8). Two mice are cisplatin 72h. Do not over-read sign on n=4.

---

## Decision

| Question | Answer |
|---|---|
| 30w 10x h5 on GEO? | **Yes** — `GSE154977_mmLung10x_cis_dSp_rawCount.h5` (custom COO, not Cell Ranger feature-barcode) |
| Seurat object? | **Yes** — `objects/gse154977_kp30w_seurat.rds` (CreateSeuratObject + QC + UMAP; 254 Mb, gitignored) |
| Usable for epithelial Cldn4? | **Yes** — 5,610 / 11,017 cells Cldn4+; 10,361 Epcam+ |
| Usable for T/NK fraction? | **No** — CD45− FACS tumor-only libraries. Design no-go, not a missing-file no-go |
| High vs low IFN / MHC / TJ? | **Mouse-level table written.** Spearman n=4 at floor; Q4−Q1 locked |
| IFN/MHC down in Cldn4-high? | **No** — IFN ρ=+0.80 exact p=0.333; MHC ρ=+1.00 exact p=0.083 |
| Honest primary n | **4 KP mice** (not 11,017 cells; not 4 libraries as if they were extra mice) |
| Join human concordant-4? | **No** — T/NK inverse not testable; n=4 too thin |
| Dual-high TACSTD2 × CLDN4 | **not defined** |

T/NK stop rule: *score epithelial/malignant Cldn4 vs T/NK fraction at mouse level*. There is no T/NK fraction. Do not invent one from leftover `Cd3d`.

---

## What the public objects actually are

Marjanovic et al., *Cancer Cell* 2020 ([PMID 32707077](https://pubmed.ncbi.nlm.nih.gov/32707077/); SuperSeries [GSE152607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152607)). Series [GSE154977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154977): “Single-cell RNAseq by 10x of advanced (30week) mouse lung tumors, comparing response to cisplatin at 72h of treatment.”

GEO lists four GSM (HTML brief view sometimes truncates the fourth):

| GSM | title | treatment | cells in h5 |
|---|---|---|---:|
| GSM4685281 | KP_30w_ND_m3_PT | no-drug | 3945 |
| GSM4685282 | KP_30w_ND_m4_PT | no-drug | 4008 |
| GSM4685283 | KP_30w_Cis72_m5_PT | cisplatin 72h | 2065 |
| GSM4685284 | KP_30w_Cis72_m6_PT | cisplatin 72h | 999 |

Source name on each GSM: “Lung tissue, AT2 cells.” 10x V2. Depositor already ran Cell Ranger + DropletUtils empty-drop filter.

| Public object | Size | Used here |
|---|---|---|
| `GSE154977_mmLung10x_cis_dSp_rawCount.h5` | 110.9 Mb | **yes** — sparse COO (`i`,`j`,`v`); Seurat counts |
| `GSE154977_mmLung10x_cis_dSp_normTPM.h5` | 118.8 Mb | no (raw counts used) |
| `GSE154977_mmLung10x_cis_smpTable.csv.gz` | 72.6 Kb | **yes** — library / barcode |
| `GSE154977_mmLung10x_cis_geneTable.csv.gz` | 212.0 Kb | **yes** — mouse symbols (`Cldn4` present) |
| `GSE154977_mmLung10x_cis_dZ_annot_annot_smpTable.csv.gz` | 271.0 Kb | **yes** — `timecourse_pred_cluster` 1–12 audit |
| `GSE154977_mmLung10x_cis_dZ_QCstat_QCstat_smpTable.csv.gz` | 161.1 Kb | **yes** — author `mitoPct` |
| SRA FASTQ | — | **no** — public processed only |

H5 is **not** a 10x `filtered_feature_bc_matrix.h5`. It is the same custom COO dump as sister Smart-seq2 GSE154989. Sister GSE154989 was **not** merged. No private 8 KL.

---

## T/NK is a design no-go

STAR Methods (same paper as GSE154989): live `tdTomato+ / CD45− / CD11b− / TER119− / CD31−`. These 30w 10x libraries are tumor epithelium, not a whole-tumor immune atlas.

| marker | cells with count>0 / 11017 | note |
|---|---:|---|
| Epcam | 10361 | epithelium present |
| Krt8 | 10714 | epithelium present |
| Sftpc | 11017 | AT2/tumor lineage present |
| Cldn4 | 5610 | epithelial Cldn4 present |
| Ptprc (CD45) | 67 | ~0.6% leak |
| Cd3d | 6 | not a T compartment |
| Nkg7 | 7 | not an NK compartment |
| Ptprc ∩ (Cd3d ∪ Cd3e ∪ Nkg7) | **1** | leak; not a fraction |

Author `timecourse_pred_cluster` is the paper’s 12 epithelial / tumor states (projected onto the 10x cisplatin libraries). Cluster 5 is Cldn4-high (91% Cldn4+; Highly mixed / HPCS marker in the paper). Cluster 11 is Epcam-lower (EMT-like in the paper); Ptprc+/Cd3d+ rate is 0 there. No T/NK cluster exists.

Do not write n=11,017 as a T/NK n. Do not write a T/NK Spearman.

---

## Honest n

| item | n | note |
|---|---:|---|
| GEO QC cells | 11017 | do not quote as analysis n |
| cells after Seurat QC (`nFeature≥200`, `percent.mt<20`) | 11017 | all depositor cells already passed |
| libraries / GSM | 4 | one 10x library per mouse |
| biological mice | **4** | honest unit |
| ND mice | 2 | m3, m4 |
| Cis72 mice | 2 | m5, m6 |
| PRIMARY mice (n_cells≥20) | **4** | all pass |
| T/NK-fraction mice | **0** | FACS CD45− |
| Q4 vs Q1 tails | **locked** | needs n_units≥8 |

Primary unit = biological KP mouse. Do not write n=11,017 cells. Do not treat four libraries as extra mice.

---

## Mouse-level Cldn4 vs IFN / MHC / TJ

Seurat `NormalizeData` (log1p of 10k-normalized counts). Family score = mean of mapped mouse genes in the data slot. IFN = Hallmark IFNα∪IFNγ (215 genes mapped). MHC-I/APM = curated mouse set (24 genes). TJ = KEGG∪GOBP minus `Cldn4` (196 genes). AddModuleScore was computed as a Seurat-native check; the claim table uses the simple family mean (same lock as GSE154989).

### Primary: 4 KP mice (2 ND + 2 Cis72)

| family | n | Spearman ρ (exact p) | Q4−Q1 | note |
|---|---:|---|---|---|
| IFN | 4 | +0.800 (0.333) | locked | nearly flat (0.197–0.218); not down |
| MHC-I/APM | 4 | +1.000 (0.083) | locked | perfect rank; exact p is **not** 0; not down |
| TJ (Cldn4 held out) | 4 | +1.000 (0.083) | locked | expected junction sanity; same n=4 limit |

IFN/MHC **down** would be ρ<0. Observed sign is the opposite and **not** significant on the exact test. Do not quote asymptotic p=0. Do not claim MHC-up. ND-only n=2 is below the Spearman lock (empty). Wilcoxon Cldn4 ~ treatment: exact p=0.667 (2 vs 2).

### Per-mouse primary units

| mouse | GSM | treatment | n_cells | Cldn4 mean | Cldn4 %pos | IFN | MHC-I/APM | TJ | T/NK leak |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| KP_30w_ND_m3 | GSM4685281 | ND | 3945 | 0.311 | 0.360 | 0.197 | 0.328 | 0.324 | 0 |
| KP_30w_ND_m4 | GSM4685282 | ND | 4008 | 0.624 | 0.598 | 0.201 | 0.409 | 0.350 | 0 |
| KP_30w_Cis72_m5 | GSM4685283 | Cis72 | 2065 | 0.535 | 0.562 | 0.202 | 0.385 | 0.337 | 1 cell |
| KP_30w_Cis72_m6 | GSM4685284 | Cis72 | 999 | 0.823 | 0.636 | 0.218 | 0.435 | 0.362 | 0 |

Honest Spearman n = **4** (not 11,017 cells). Fisher-z / exact permutation on n=4 is the whole test. Cisplatin is a real confound; the table mixes it because that is the deposited series.

---

## Seurat methods (locked)

- Cldn4 only. Tacstd2 is not a gate.
- Matrix: GEO `GSE154977_mmLung10x_cis_dSp_rawCount.h5` only. Undetected = 0.
- `CreateSeuratObject` → `nFeature_RNA≥200` and `percent.mt<20` → `NormalizeData` / `FindVariableFeatures` / `ScaleData` / `RunPCA` (30) / `FindNeighbors` (dims 1–20) / `FindClusters` (res 0.4) / `RunUMAP`.
- Mouse unit = library with trailing `_PT` stripped (`KP_30w_ND_m3`, …). One mouse per library.
- Primary cohort: all 4 KP mice (all ≥20 cells).
- Family scores: mean of mapped mouse genes on the Seurat data slot. Spearman uses the **exact** test when n≤10.
- Spearman requires n≥4 finite pairs. Q4 vs Q1 requires n_units≥8 — **not computed**.
- T/NK fraction is **not** computed.
- No SRA. No GSE154989 merge. No human concordant-4 merge. No private 8 KL.
- Reproduce: `/home/ubuntu/micromamba/envs/seurat/bin/Rscript methods/seurat_gse154977_cldn4/analyze.R`
- Plots: Seurat `DimPlot` (library, clusters, author cluster, treatment) and `VlnPlot` (Cldn4 / Epcam / Ptprc / Cd3d / Nkg7).

---

## How to read this

- This is **epithelial-only 10x**, not a whole-tumor atlas.
- T/NK are truly absent (tumor-only libraries). Only epithelial Cldn4 vs IFN / MHC / TJ is scored.
- Do not quote n=11,017 cells as the mouse n.
- Do not treat leftover `Cd3d` as T/NK fraction.
- Do not quote MHC/TJ ρ=1 as a result with p=0. Exact p=0.083 on n=4.
- Thesis is unchanged. Additive mouse result: Cldn4 is present and the libraries are tumor epithelium; T/NK is unscorable; IFN/MHC is not down; n=4 is honest and thin.
