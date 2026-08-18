# FINDING — Seurat GSE154989 FACS epithelium: Cldn4-only

**ADDITIVE public MOUSE. Cldn4-only. Thesis already correct. No dual-high.**

Primary analysis is **R + Seurat 5.5.1** (`CreateSeuratObject`) on the public GEO processed matrix. Epithelial Cldn4 vs IFN / MHC / TJ was scored at **biological mouse** level. **T/NK fraction cannot be scored.** Cells were FACS `tdTomato+ / CD45− / CD11b− / TER119− / CD31−`. Residual `Cd3d` / `Nkg7` / `Ptprc` is leak, not a compartment. T/NK Spearman is **empty**. This accession does **not** join the human concordant-4 pool. No private 8 KL.

Primary KP (n=15 mice, ≥20 cells): Seurat `AddModuleScore` IFN and MHC-I/APM are **not** significantly down in Cldn4-high (IFN ρ=−0.100 p=0.723; MHC ρ=−0.214 p=0.443). Mean depositor TPM of the same families is the opposite sign and also ns (IFN ρ=+0.339 p=0.216; MHC ρ=+0.139 p=0.621). TJ (Cldn4 held out) tracks Cldn4 in both scores (module ρ=+0.818 p=0.0002; TPM ρ=+0.757 p=0.001) — sanity that Cldn4 is a real tight-junction gene here, not a join key.

---

## Decision

| Question | Answer |
|---|---|
| Seurat installed and ran? | **Yes** — Seurat 5.5.1 / SeuratObject 5.4.0; `CreateSeuratObject` used |
| Processed matrix on GEO? | **Yes** — `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5` |
| Usable for epithelial Cldn4? | **Yes** — 1,960 / 3,891 cells Cldn4+; 3,569 Epcam+ |
| Usable for T/NK fraction? | **No** — CD45− FACS. Design no-go, not a missing-file no-go |
| High vs low IFN / MHC / TJ? | **Yes, mouse-level** — primary KP n=**15** |
| IFN/MHC down in Cldn4-high (KP)? | **No** — Seurat module scores ns (p>0.4); mean TPM not down |
| Honest primary n | **15 KP mice** (not 3,891 cells; not 39 deposited `mouseID`) |
| Join human concordant-4? | **No** — T/NK inverse not testable; KP IFN/MHC not a significant down |
| Dual-high TACSTD2 × CLDN4 | **not defined** |
| Private 8 KL | **not used** |

T/NK stop rule: *score epithelial/malignant Cldn4 vs T/NK fraction at mouse level*. There is no T/NK fraction. Do not invent one from leftover `Cd3d`.

---

## What the public objects actually are

Marjanovic et al., *Cancer Cell* 2020 ([PMID 32707077](https://pubmed.ncbi.nlm.nih.gov/32707077/); SuperSeries [GSE152607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152607)). Series [GSE154989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154989): “Emergence of a high-plasticity cell state during lung cancer evolution, single-cell RNAseq timecourse.”

STAR Methods: live `tdTomato+ / CD45− / CD11b− / TER119− / CD31−` cells, modified SMART-Seq2. 3,891 QC full-length transcriptomes from K / KP / T lungs at defined time points. Paper text says “39 mice”; GEO `mouseID` has **39** values because **KP 30w tumors are split** (`m1_T1`…`m3_T9`). Collapsing `_T#` gives **30 biological animals**.

| Public object | Size | Used here |
|---|---|---|
| `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5` | 190.9 Mb | **yes** — MATLAB v7.3 COO (`i`=gene 1…52638, `j`=cell 1…3891, `v`=depositor normTPM) |
| `GSE154989_mmLungPlate_fQC_dSp_rawCount.h5` | 75.6 Mb | no (TPM used) |
| `GSE154989_mmLungPlate_fQC_smpTable.csv.gz` | 26.1 Kb | **yes** — `mouseID`, stage |
| `GSE154989_mmLungPlate_fQC_dZ_annot_smpTable.csv.gz` | 153.9 Kb | **yes** — `clusterK12` audit + author tSNE |
| `GSE154989_mmLungPlate_fQC_geneTable.csv.gz` | 762.1 Kb | **yes** — mouse symbols (`Cldn4` present) |
| SRA FASTQ | — | **no** — public processed only |

Sister 10x series GSE154977 (3 cisplatin libraries) was **not** opened.

H5 is not 10x. It is a custom COO dump. Undetected genes are absent from the triplet and are scored as 0.

Seurat was given that sparse TPM matrix via `CreateSeuratObject`. The `data` layer is `log1p(TPM)`. `NormalizeData` was **not** run as if these were UMI counts.

---

## T/NK is a design no-go

| marker | cells with v>0 / 3891 | note |
|---|---:|---|
| Epcam | 3569 | epithelium present |
| Cldn4 | 1960 | epithelial Cldn4 present |
| Ptprc (CD45) | 122 | ~3% leak |
| Cd3d | 31 | not a T compartment |
| Nkg7 | 36 | not an NK compartment |
| any T/NK audit gene | 202 | still leak; two KP 12w plates are dirtier |

`clusterK12` is 12 epithelial / tumor states. Ptprc+ rate is ≤8.3% in every cluster. Cluster 11 (89 cells) is Epcam-low; it is not a T/NK cluster and was not used to mint a fraction.

Do not write n=3,891 as a T/NK n. Do not write a T/NK Spearman.

---

## Honest n

| item | n | note |
|---|---:|---|
| GEO QC cells | 3891 | do not quote as analysis n |
| deposited `mouseID` | 39 | paper’s “39 mice”; KP 30w = 12 tumors from 3 mice |
| biological animals (strip `_T#`) | 30 | true mouse unit |
| T / normal AT2 | 5 | not KP; Cldn4-low, MHC-high |
| K-only | 9 | Kras; Trp53 WT |
| KP animals | 16 | named GEMM |
| KP animals with ≥20 cells (**PRIMARY**) | **15** | dropped `KP_2w_ND_m1` (4 cells) |
| K+KP animals with ≥20 cells | 24 | sensitivity only |
| T/NK-fraction mice | **0** | FACS CD45− |

Primary unit = biological KP mouse with ≥20 epithelial cells. Q4 vs Q1 tails are **4 vs 4**. That is honest and thin.

---

## Mouse-level Cldn4 vs IFN / MHC / TJ (Seurat)

Family genes: IFN = Hallmark IFNα∪IFNγ mapped to mouse (217 / 224 present). MHC-I/APM = curated mouse set (24 genes: `H2-K1`, `H2-D1`, `H2-Q*`, `H2-T23`, `B2m`, TAP/PSMB/ERAP machinery). TJ = KEGG∪GOBP minus `Cldn4` (200 genes).

Two scores per family, both aggregated to the mouse mean:

1. **Seurat `AddModuleScore`** (control-gene normalized; primary Seurat arm).
2. **Mean depositor normTPM** of the same genes (companion; matches the public matrix as deposited).

The `ifn_mhc_down` column in `family_q4q1.tsv` is a mechanical sign check (ρ<0 or Q4−Q1 Δ<0). It is **not** a significance call.

### Primary: KP, n=15 — Seurat AddModuleScore

| family | n | n_Q1 / n_Q4 | Spearman ρ (p) | Q4−Q1 Δ median | rank-biserial r | MW p |
|---|---:|---|---|---:|---:|---:|
| IFN | 15 | 4 / 4 | −0.100 (0.723) | −0.013 | −0.125 | 0.885 |
| MHC-I/APM | 15 | 4 / 4 | −0.214 (0.443) | −0.066 | −0.375 | 0.470 |
| TJ (Cldn4 held out) | 15 | 4 / 4 | **+0.818 (0.0002)** | +0.071 | +1.000 | 0.030 |

IFN/MHC **down** would require a real negative effect, not a 15-mouse ns wobble. Primary KP module scores do not support that claim.

### Primary: KP, n=15 — mean depositor TPM (companion)

| family | n | n_Q1 / n_Q4 | Spearman ρ (p) | Q4−Q1 Δ median | rank-biserial r | MW p |
|---|---:|---|---|---:|---:|---:|
| IFN | 15 | 4 / 4 | +0.339 (0.216) | +0.228 | +0.625 | 0.194 |
| MHC-I/APM | 15 | 4 / 4 | +0.139 (0.621) | +0.193 | +0.250 | 0.665 |
| TJ (Cldn4 held out) | 15 | 4 / 4 | **+0.757 (0.001)** | +0.550 | +1.000 | 0.030 |

Mean TPM is **not** down. Sign disagreement vs `AddModuleScore` is expected: module scores subtract random control genes; raw family means do not. Neither arm is a significant IFN/MHC-down result.

### Per-mouse primary units

| mouse | week | n_cells | n_tumors | Cldn4 mean | Cldn4 %pos | IFN module | MHC module | TJ module |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| KP_2w_ND_m2 | 2w | 27 | 1 | 1.622 | 0.296 | 0.116 | 0.360 | 0.231 |
| KP_12w_ND_m1 | 12w | 131 | 1 | 4.012 | 0.687 | 0.015 | −0.102 | 0.223 |
| KP_12w_ND_m2 | 12w | 124 | 1 | 4.942 | 0.758 | 0.102 | 0.033 | 0.261 |
| KP_12w_ND_m3 | 12w | 74 | 1 | 4.575 | 0.770 | 0.029 | −0.094 | 0.237 |
| KP_12w_ND_m4 | 12w | 62 | 1 | 3.713 | 0.435 | 0.118 | 0.212 | 0.243 |
| KP_12w_ND_m5 | 12w | 48 | 1 | 4.506 | 0.667 | 0.033 | −0.017 | 0.239 |
| KP_12w_ND_m6 | 12w | 52 | 1 | 3.452 | 0.519 | 0.097 | 0.170 | 0.227 |
| KP_20w_ND_m2 | 20w | 43 | 1 | 4.556 | 0.744 | 0.092 | 0.133 | 0.243 |
| KP_20w_ND_m3 | 20w | 60 | 1 | 5.958 | 0.817 | 0.122 | 0.160 | 0.265 |
| KP_20w_ND_m4 | 20w | 173 | 1 | 3.014 | 0.561 | 0.047 | −0.061 | 0.208 |
| KP_20w_ND_m5 | 20w | 126 | 1 | 1.495 | 0.317 | 0.062 | −0.001 | 0.143 |
| KP_20w_ND_m6 | 20w | 136 | 1 | 1.598 | 0.331 | 0.096 | 0.083 | 0.176 |
| KP_30w_ND_m1 | 30w | 520 | 4 | 2.945 | 0.448 | 0.034 | −0.026 | 0.183 |
| KP_30w_ND_m2 | 30w | 619 | 5 | 5.226 | 0.775 | 0.026 | −0.084 | 0.240 |
| KP_30w_ND_m3 | 30w | 415 | 3 | 4.392 | 0.667 | 0.041 | −0.047 | 0.224 |

Honest Spearman n = **15** (not 2,610 KP cells, not 12 KP-30w tumors, not 39 `mouseID`).

### Sensitivity (not the claim)

| contrast | score | family | n | ρ (p) | Q4−Q1 Δ | note |
|---|---|---|---:|---|---:|---|
| K+KP ≥20 | module | IFN | 24 | −0.302 (0.152) | −0.006 | ns |
| K+KP ≥20 | module | MHC-I/APM | 24 | −0.526 (0.008) | −0.218 | mixes K; MW p=0.128 |
| K+KP ≥20 | TPM | MHC-I/APM | 24 | −0.288 (0.173) | −0.171 | ns |
| K+KP ≥20 | either | TJ | 24 | +0.676 / +0.737 | + | same TJ sanity |
| all animals ≥20 (incl. T AT2) | TPM | MHC-I/APM | 28 | −0.500 (0.007) | −0.665 | **confounded** — normal AT2 are Cldn4-low / MHC-high |
| KP including 4-cell mouse | TPM | IFN | 16 | +0.103 (0.704) | +0.023 | still not down |
| KP 30w tumors as units | TPM | IFN | 12 | +0.566 (0.055) | +0.418 | not mice; same-mouse tumors |

Do **not** quote the all-animal MHC-down as support. That sign is T AT2 vs tumor, not Cldn4 vs IFN inside KP.

---

## Seurat plots

| file | what |
|---|---|
| `figures/seurat_umap_genotype.png` | `DimPlot` UMAP by genotype |
| `figures/seurat_umap_week.png` | `DimPlot` UMAP by week |
| `figures/seurat_featureplot_cldn4_modules.png` | `FeaturePlot` Cldn4 + IFN/MHC/TJ modules |
| `figures/seurat_author_tsne_cldn4.png` | author tSNE, `FeaturePlot` Cldn4 |
| `figures/seurat_vln_cldn4_genotype.png` | `VlnPlot` Cldn4 by genotype |
| `figures/kp_mouse_cldn4_vs_modules.png` | primary n=15 module scatter |
| `figures/kp_mouse_cldn4_vs_tpm_families.png` | primary n=15 mean-TPM scatter |
| `figures/kp_q4q1_modules.png` | Cldn4 Q4 vs Q1 (4 vs 4) |
| `figures/honest_n.png` | do not quote 3891 cells |

UMAP is for display of this already-QC FACS epithelium. It is not used to invent a T/NK cluster.

---

## Methods (locked)

- Cldn4 only. Tacstd2 is not a gate.
- Matrix: GEO `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5` only. Undetected = 0.
- R + Seurat: `CreateSeuratObject(counts = sparse TPM)`. `data = log1p(TPM)`. `AddModuleScore` for IFN / MHC / TJ. PCA/UMAP for Seurat plots only.
- Mouse unit = `mouseID` with trailing `_T#` stripped. KP 30w tumors from the same `m#` are one mouse.
- Primary cohort: genotype `KP`, n_cells ≥ 20. Dropped `KP_2w_ND_m1` (4 cells).
- Family scores: Hallmark IFNα∪IFNγ; curated mouse MHC-I/APM; KEGG∪GOBP TJ minus Cldn4.
- Spearman requires n≥4 finite pairs. Q4 vs Q1 requires n_units≥8 and two tails.
- T/NK fraction is **not** computed.
- No SRA. No GSE154977. No human concordant-4 merge. No private 8 KL.
- Reproduce: `Rscript methods/seurat_gse154989_cldn4/analyze.R`

---

## How to read this

- This is **epithelial-only Smart-seq2**, not a whole-tumor 10x atlas.
- Do not quote n=3,891 cells or n=39 `mouseID` as the mouse n.
- Do not treat leftover `Cd3d` as T/NK fraction.
- Do not treat T-AT2–driven MHC-down as a KP finding.
- Do not treat a non-significant negative `AddModuleScore` ρ as IFN/MHC-down.
- Thesis is unchanged. Additive mouse result: Seurat ran; Cldn4 is present and tracks TJ; T/NK is unscorable; KP IFN/MHC is not a significant down.
