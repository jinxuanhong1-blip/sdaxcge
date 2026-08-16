# REWORK A1 deconvolution — TCGA-LUAD TACSTD2 vs xCell / MCP-counter / ESTIMATE / GEP18 ssGSEA / Danaher

**Self-contained. Public data only. Written to be read without the rest of the repo.**

**Cohort:** TCGA-LUAD primary tumors (`-01`). 515 with RNA-seq; **n = 502** with RNA-seq + ABSOLUTE purity (partial-correlation set).

**Question:** The prior LUAD-only A1 scan was weak (|ρ| ≤ 0.22 on CIBERSORT / Wolf / marker genes). Does a different immune *definition* — official deconvolution and signature methods — change that?

## Verdict: STILL_WEAK

STILL WEAK for the named immune-infiltration / T-cell-inflamed scores. xCell ImmuneScore, ESTIMATE ImmuneScore, GEP18 ssGSEA, MCP T/CD8/cytotoxicity, and Danaher T/CD8/cytotoxic/Total TILs all stay inside the prior LUAD |ρ|≤0.22 bound (largest primary |ρ| = 0.126). Two non-primary xCell subsets exceed 0.22 and are reported, not hidden: xCell_T_cell_NK ρ=0.250; xCell_B_cell_plasma ρ=-0.226. They do not form a coherent immune-cold pattern (one is positive).

TACSTD2 vs ABSOLUTE purity: Spearman ρ = 0.007, p = 0.869, n = 502. Purity adjustment is nearly a no-op for TACSTD2 itself (same observation as the prior LUAD-only A1). It is **not** optional for ESTIMATE / xCell ImmuneScore / MCP totals, which track leukocyte content and are strongly (negatively) correlated with purity.

## Why this rework exists

Prior LUAD-only A1 (`results/w200/A1_LUAD/`): 32 features, all |ρ| ≤ 0.22 after ABSOLUTE; no association with total leukocyte fraction (ρ = 0.03). Histology rework (`results/rework/A1_histology/`): LUAD GEP18 ρ ≈ 0, ESTIMATE ImmuneScore weakly positive / NS; CD8/CYT small negatives. Those used CD8A, CYT, unweighted GEP18 z-mean, and official ESTIMATE. They did **not** run xCell, MCP-counter, GEP18 *ssGSEA*, or Danaher.

This file runs those five named methods on the same TCGA-LUAD Xena freeze + ABSOLUTE, and reports **every population** each method emits.

## Analysis set

| Filter | n |
|---|---:|
| Xena HiSeqV2 primary tumors (`-01`) | 515 |
| + ABSOLUTE purity | 502 |
| + official ESTIMATE ImmuneScore | 515 |
| + TIMER2 xCell / MCP-counter | 515 |

Primary tumors only. One row per 15-character barcode. Replicate aliquots averaged. Partial n is TACSTD2 + feature + ABSOLUTE (exact n is in every table cell).

## ALL definitions

### TACSTD2
- Gene symbol `TACSTD2` (TROP2, GA733-1, EGP-1). Ensembl `ENSG00000184292`.
- Measure: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` **log2(RSEM normalized_count + 1)**.
- Not protein. Not ADC response. Not ICI outcome.

### Samples
- TCGA-LUAD only. No LUSC. No OncoSG.
- Keep barcode sample-type `01` (primary solid tumor). Drop `11` normals and other types.
- Identifier: first 15 characters (`TCGA-XX-XXXX-01`).

### Covariate
- PanCanAtlas **ABSOLUTE** `purity` from GDC file `4f277128-f793-4354-a13d-30cc7fe9f6b5` (`TCGA_mastercalls.abs_tables_JSedit.fixed.txt`).
- Not ESTIMATE TumorPurity. Not methylation leukocyte fraction.

### Partial Spearman (primary statistic)
Two equivalent one-covariate estimators; both are reported.

1. **Algebraic (primary, matches original pooled A1 / histology rework):**
   \( r_{xy\cdot z} = (r_{xy} - r_{xz} r_{yz}) / \sqrt{(1-r_{xz}^2)(1-r_{yz}^2)} \)
   with pairwise Spearman rhos. t-test, **df = n − 3**.
2. **Residual-rank (sensitivity, matches LUAD-only A1):** Pearson of rank residuals after OLS of ranked x and y on ranked purity.

Also reported: unadjusted Spearman; 95% Fisher-z CI with variance `1/(n−4)`; feature vs purity Spearman (so ImmuneScore-like composites are not silently called “immune”).

BH-FDR is shown two ways: across **all** reported features, and across the **15** pre-specified primary features. Raw p-values are always shown. No feature is hidden.

### 1. xCell (Aran *Genome Biol* 2017)
- **Source used:** TIMER2.0 precomputed scores (`infiltration_estimation_for_tcga.csv.gz`, Li *NAR* 2020), produced by **immunedeconv** running the official xCell algorithm on TCGA TPM.
- **Not recomputed** on Xena log2(RSEM+1). xCell’s power/calibration/spillover steps expect TPM-like input; the published TIMER2 table is the standard TCGA xCell freeze.
- **Primary features:** `immune score`, `T cell CD8+`, `microenvironment score`.
- **All other xCell columns** in TIMER2 are reported (B, CD4 subsets, Tcm/Tem/naive CD8, NK, macrophages M1/M2, Tregs, stroma score, progenitors, etc.).
- xCell ImmuneScore is a composite of immune cell scores after spillover compensation. It is expected to anti-correlate with ABSOLUTE purity.

### 2. MCP-counter (Becht *Genome Biol* 2016)
Two implementations, both reported (they are not identical).

**A. TIMER2 / immunedeconv (official algorithm on TCGA TPM).** Columns: T cell, T cell CD8+, cytotoxicity score, NK cell, B cell, Monocyte, Macrophage/Monocyte, myeloid DC, Neutrophil, Endothelial cell, CAF.

**B. Xena recompute (definition-transparent, same matrix as TACSTD2).** Score = **mean of log2(RSEM+1)** of the official Becht marker genes for each of the 10 populations (T cells, CD8 T cells, Cytotoxic lymphocytes, B lineage, NK cells, Monocytic lineage, Myeloid dendritic cells, Neutrophils, Endothelial cells, Fibroblasts). Missing genes are **not imputed**; coverage is in `gene_coverage.tsv`.

Xena symbol aliases (do not treat as new genes): `ADGRL4` → `ELTD1`; `DIPK2B` → `CXorf36`. Absent from HiSeqV2 and dropped: `CHRM3-AS2`, `IGKC`, `KIR3DS1`, `WFDC21P`.

CD8 T cells is a **single gene** (`CD8B`) in the official MCP signature. Cytotoxic lymphocytes is a 7-gene set that includes `CD8A`.

### 3. ESTIMATE ImmuneScore (Yoshihara *Nat Commun* 2013)
- Official MD Anderson **RNAseqV2** table `lung_adenocarcinoma_RNAseqV2.txt`.
- Primary: `Immune_score`. Secondary: `Stromal_score`, `ESTIMATE_score` (= Immune + Stromal).
- **Not recomputed** by ssGSEA in this slice (official table returned HTTP 200).
- ImmuneScore is built to track leukocyte content. ABSOLUTE partialling is required before calling it an “immune” association.

### 4. GEP18 ssGSEA (Ayers *JCI* 2017)
18 genes (all present in Xena):
`CCL5, CD27, CD274, CD276, CD8A, CMKLR1, CXCL9, CXCR6, HLA-DQA1, HLA-DRB1, HLA-E, IDO1, LAG3, NKG7, PDCD1LG2, PSMB10, STAT1, TIGIT`.

- **Primary score:** Barbie/GSVA **ssGSEA**, tau = 0.25, running-sum ES, ranks scaled to 1..10000 over the full HiSeqV2 gene universe (20,530 genes). No Merck NanoString TIS weights (not public). No inter-gene-set `ssgsea.norm` (one set).
- **Sensitivity:** unweighted mean of within-LUAD z-scores of log2(RSEM+1) — the histology-rework definition — and simple mean of the 18 log2 values.
- This is **not** the clinical NanoString TIS.

### 5. Danaher (Danaher *JITC* 2017 Table 1)
14 populations, 60 marker genes. Score = **average log2 expression** of that population’s markers (paper Methods). Xena is already log2(RSEM+1), so this is the paper’s estimator.

| Population | Genes |
|---|---|
| B-cells | BLK, CD19, FCRL2, MS4A1, KIAA0125, TNFRSF17, TCL1A, SPIB, PNOC |
| CD45 | PTPRC (paper typo: PTRPC) |
| Cytotoxic cells | PRF1, GZMA, GZMB, NKG7, GZMH, KLRK1, KLRB1, KLRD1, CTSW, GNLY |
| DC | CCL13, CD209, HSD11B1 |
| Exhausted CD8 | LAG3, CD244, EOMES, PTGER4 |
| Macrophages | CD68, CD84, CD163, MS4A4A |
| Mast cells | TPSB2, TPSAB1, CPA3, MS4A2, HDC |
| Neutrophils | FPR1, SIGLEC5, CSF3R, FCAR, FCGR3B, CEACAM3, S100A12 |
| NK CD56dim | KIR2DL3, KIR3DL1, KIR3DL2, IL21R |
| NK cells | XCL1, XCL2, NCR1 |
| T-cells | CD6, CD3D, CD3E, SH2D1A, TRAT1, CD3G |
| Th1 cells | TBX21 |
| Treg | FOXP3 |
| CD8 T cells | CD8A, CD8B |

**Total TILs (primary):** mean of the 11 populations the paper retained after requiring correlation with PTPRC > 0.6, i.e. exclude DC, Treg, and mast cells. Sensitivity: `Danaher_Total_TILs_empirical` recomputes that cutoff on this LUAD set (populations used: Danaher_B_cells, Danaher_Cytotoxic_cells, Danaher_DC, Danaher_Exhausted_CD8, Danaher_Macrophages, Danaher_NK_CD56dim, Danaher_T_cells, Danaher_Th1_cells, Danaher_Treg, Danaher_CD8_T_cells).

Th1 and Treg are single-gene scores; the paper says so. They are not multi-gene deconvolution.

## Primary result — ABSOLUTE partial Spearman

| Feature | Unadj ρ | Unadj p | Partial ρ | Partial p | FDR (15) | n | vs purity ρ |
|---|---:|---:|---:|---:|---:|---:|---:|
| `ESTIMATE_ImmuneScore` | 0.036 | 0.41 | 0.072 | 0.109 | 0.16 | 502 | -0.616 |
| `GEP18_ssGSEA` | -0.039 | 0.381 | -0.035 | 0.429 | 0.495 | 502 | -0.550 |
| `xCell_immune_score` | 0.002 | 0.965 | 0.022 | 0.625 | 0.669 | 502 | -0.562 |
| `xCell_T_cell_CD8` | -0.117 | 0.00809 | -0.126 | 0.00479 | 0.036 | 502 | -0.315 |
| `xCell_microenvironment_score` | -0.008 | 0.852 | 0.017 | 0.704 | 0.704 | 502 | -0.606 |
| `TIMER2_MCP_T_cell` | -0.070 | 0.111 | -0.079 | 0.077 | 0.128 | 502 | -0.508 |
| `TIMER2_MCP_T_cell_CD8` | -0.114 | 0.00937 | -0.118 | 0.00808 | 0.036 | 502 | -0.388 |
| `TIMER2_MCP_cytotoxicity_score` | -0.111 | 0.0121 | -0.117 | 0.00851 | 0.036 | 502 | -0.424 |
| `Xena_MCP_T_cells` | -0.059 | 0.18 | -0.066 | 0.14 | 0.175 | 502 | -0.552 |
| `Xena_MCP_CD8_T_cells` | -0.105 | 0.017 | -0.106 | 0.0181 | 0.0387 | 502 | -0.375 |
| `Xena_MCP_Cytotoxic_lymphocytes` | -0.106 | 0.0162 | -0.112 | 0.012 | 0.036 | 502 | -0.462 |
| `Danaher_T_cells` | -0.065 | 0.143 | -0.070 | 0.117 | 0.16 | 502 | -0.533 |
| `Danaher_CD8_T_cells` | -0.110 | 0.0127 | -0.115 | 0.01 | 0.036 | 502 | -0.430 |
| `Danaher_Cytotoxic_cells` | -0.103 | 0.0189 | -0.109 | 0.0146 | 0.0364 | 502 | -0.513 |
| `Danaher_Total_TILs` | -0.093 | 0.0344 | -0.100 | 0.0254 | 0.0477 | 502 | -0.591 |

Dashed lines on the forest plot mark the prior LUAD-only |ρ| = 0.22 bound. 8 / 15 primary tests have BH-FDR < 0.05.

## All populations, by method

### ESTIMATE
| Feature | Unadj ρ | Unadj p | Partial ρ | Partial p | FDR (all) | n | vs purity ρ |
|---|---:|---:|---:|---:|---:|---:|---:|
| `ESTIMATE_StromalScore` | -0.020 | 0.65 | 0.007 | 0.869 | 0.902 | 502 | -0.596 |
| `ESTIMATE_Score` | 0.006 | 0.899 | 0.041 | 0.364 | 0.523 | 502 | -0.659 |
| `ESTIMATE_ImmuneScore` | 0.036 | 0.41 | 0.072 | 0.109 | 0.218 | 502 | -0.616 |

### GEP18
| Feature | Unadj ρ | Unadj p | Partial ρ | Partial p | FDR (all) | n | vs purity ρ |
|---|---:|---:|---:|---:|---:|---:|---:|
| `GEP18_ssGSEA` | -0.039 | 0.381 | -0.035 | 0.429 | 0.593 | 502 | -0.550 |
| `GEP18_logmean` | -0.033 | 0.454 | -0.028 | 0.538 | 0.671 | 502 | -0.582 |
| `GEP18_zmean` | -0.022 | 0.62 | -0.013 | 0.775 | 0.836 | 502 | -0.588 |

### xCell (TIMER2, all columns)
| Feature | Unadj ρ | Unadj p | Partial ρ | Partial p | FDR (all) | n | vs purity ρ |
|---|---:|---:|---:|---:|---:|---:|---:|
| `xCell_B_cell_plasma` | -0.210 | 1.50e-06 | -0.226 | 3.20e-07 | 1.31e-05 | 502 | -0.163 |
| `xCell_T_cell_CD8` | -0.117 | 0.00809 | -0.126 | 0.00479 | 0.0476 | 502 | -0.315 |
| `xCell_Endothelial_cell` | -0.135 | 0.00217 | -0.117 | 0.00874 | 0.051 | 502 | -0.303 |
| `xCell_Plasmacytoid_dendritic_cell` | -0.099 | 0.0242 | -0.101 | 0.0232 | 0.0732 | 502 | -0.478 |
| `xCell_T_cell_CD4_Th2` | -0.096 | 0.0301 | -0.095 | 0.0331 | 0.0936 | 502 | -0.110 |
| `xCell_Common_lymphoid_progenitor` | -0.071 | 0.109 | -0.093 | 0.038 | 0.1 | 502 | 0.264 |
| `xCell_T_cell_CD4_Th1` | -0.077 | 0.0827 | -0.089 | 0.0459 | 0.114 | 502 | 0.063 |
| `xCell_T_cell_CD8_naive` | -0.084 | 0.0575 | -0.087 | 0.0521 | 0.125 | 502 | -0.001 |
| `xCell_T_cell_CD4_memory` | -0.077 | 0.0808 | -0.084 | 0.06 | 0.137 | 502 | -0.161 |
| `xCell_T_cell_CD4_non_regulatory` | -0.077 | 0.0789 | -0.073 | 0.103 | 0.216 | 502 | -0.076 |
| `xCell_B_cell` | -0.062 | 0.157 | -0.065 | 0.144 | 0.262 | 502 | -0.446 |
| `xCell_T_cell_regulatory_Tregs` | -0.057 | 0.195 | -0.060 | 0.177 | 0.31 | 502 | -0.049 |
| `xCell_Eosinophil` | -0.052 | 0.237 | -0.055 | 0.218 | 0.357 | 502 | -0.079 |
| `xCell_T_cell_CD8_central_memory` | -0.044 | 0.32 | -0.050 | 0.267 | 0.421 | 502 | -0.436 |
| `xCell_T_cell_CD8_effector_memory` | -0.045 | 0.308 | -0.043 | 0.332 | 0.505 | 502 | -0.271 |
| `xCell_B_cell_naive` | -0.041 | 0.351 | -0.035 | 0.434 | 0.593 | 502 | -0.176 |
| `xCell_Neutrophil` | -0.022 | 0.621 | -0.025 | 0.574 | 0.703 | 502 | 0.027 |
| `xCell_B_cell_memory` | -0.029 | 0.511 | -0.024 | 0.596 | 0.718 | 502 | -0.346 |
| `xCell_NK_cell` | -0.020 | 0.653 | -0.021 | 0.633 | 0.721 | 502 | -0.003 |
| `xCell_T_cell_gamma_delta` | 0.005 | 0.904 | 0.004 | 0.926 | 0.937 | 502 | 0.100 |
| `xCell_Hematopoietic_stem_cell` | -0.007 | 0.874 | 0.009 | 0.84 | 0.884 | 502 | -0.136 |
| `xCell_Macrophage_M2` | -0.002 | 0.964 | 0.012 | 0.794 | 0.845 | 502 | -0.212 |
| `xCell_T_cell_CD4_central_memory` | 0.020 | 0.652 | 0.014 | 0.746 | 0.816 | 502 | 0.102 |
| `xCell_microenvironment_score` | -0.008 | 0.852 | 0.017 | 0.704 | 0.79 | 502 | -0.606 |
| `xCell_T_cell_CD4_naive` | 0.015 | 0.734 | 0.021 | 0.632 | 0.721 | 502 | -0.347 |
| `xCell_immune_score` | 0.002 | 0.965 | 0.022 | 0.625 | 0.721 | 502 | -0.562 |
| `xCell_Macrophage_M1` | 0.008 | 0.85 | 0.030 | 0.497 | 0.657 | 502 | -0.458 |
| `xCell_stroma_score` | 0.001 | 0.975 | 0.031 | 0.494 | 0.657 | 502 | -0.406 |
| `xCell_Class_switched_memory_B_cell` | 0.040 | 0.37 | 0.038 | 0.392 | 0.554 | 502 | -0.342 |
| `xCell_Mast_cell` | 0.033 | 0.451 | 0.043 | 0.34 | 0.507 | 502 | -0.104 |
| `xCell_Macrophage` | 0.027 | 0.547 | 0.048 | 0.286 | 0.442 | 502 | -0.420 |
| `xCell_Monocyte` | 0.061 | 0.164 | 0.093 | 0.0379 | 0.1 | 502 | -0.459 |
| `xCell_Myeloid_dendritic_cell_activated` | 0.074 | 0.0941 | 0.097 | 0.0304 | 0.0889 | 502 | -0.554 |
| `xCell_Granulocyte_monocyte_progenitor` | 0.117 | 0.00792 | 0.114 | 0.0106 | 0.051 | 502 | -0.164 |
| `xCell_Cancer_associated_fibroblast` | 0.080 | 0.0713 | 0.115 | 0.01 | 0.051 | 502 | -0.380 |
| `xCell_Myeloid_dendritic_cell` | 0.094 | 0.0329 | 0.115 | 0.00999 | 0.051 | 502 | -0.457 |
| `xCell_Common_myeloid_progenitor` | 0.110 | 0.0128 | 0.115 | 0.0097 | 0.051 | 502 | 0.031 |
| `xCell_T_cell_CD4_effector_memory` | 0.169 | 0.000113 | 0.158 | 0.000374 | 0.00613 | 502 | -0.057 |
| `xCell_T_cell_NK` | 0.256 | 3.94e-09 | 0.250 | 1.32e-08 | 1.08e-06 | 502 | -0.177 |

### MCP-counter (TIMER2 + Xena recompute)
| Feature | Unadj ρ | Unadj p | Partial ρ | Partial p | FDR (all) | n | vs purity ρ |
|---|---:|---:|---:|---:|---:|---:|---:|
| `TIMER2_MCP_NK_cell` | -0.166 | 0.000155 | -0.165 | 0.000212 | 0.00454 | 502 | -0.369 |
| `Xena_MCP_NK_cells` | -0.165 | 0.000169 | -0.164 | 0.000221 | 0.00454 | 502 | -0.410 |
| `TIMER2_MCP_Endothelial_cell` | -0.142 | 0.00124 | -0.125 | 0.00522 | 0.0476 | 502 | -0.192 |
| `TIMER2_MCP_T_cell_CD8` | -0.114 | 0.00937 | -0.118 | 0.00808 | 0.051 | 502 | -0.388 |
| `TIMER2_MCP_cytotoxicity_score` | -0.111 | 0.0121 | -0.117 | 0.00851 | 0.051 | 502 | -0.424 |
| `Xena_MCP_Endothelial_cells` | -0.129 | 0.00335 | -0.113 | 0.0115 | 0.0518 | 502 | -0.310 |
| `Xena_MCP_Cytotoxic_lymphocytes` | -0.106 | 0.0162 | -0.112 | 0.012 | 0.0518 | 502 | -0.462 |
| `Xena_MCP_CD8_T_cells` | -0.105 | 0.017 | -0.106 | 0.0181 | 0.0667 | 502 | -0.375 |
| `TIMER2_MCP_B_cell` | -0.082 | 0.0631 | -0.092 | 0.0391 | 0.1 | 502 | -0.492 |
| `Xena_MCP_B_lineage` | -0.079 | 0.0719 | -0.086 | 0.0534 | 0.125 | 502 | -0.479 |
| `TIMER2_MCP_T_cell` | -0.070 | 0.111 | -0.079 | 0.077 | 0.166 | 502 | -0.508 |
| `Xena_MCP_T_cells` | -0.059 | 0.18 | -0.066 | 0.14 | 0.262 | 502 | -0.552 |
| `Xena_MCP_Neutrophils` | -0.033 | 0.45 | -0.030 | 0.508 | 0.661 | 502 | -0.120 |
| `Xena_MCP_Monocytic_lineage` | -0.027 | 0.54 | -0.004 | 0.92 | 0.937 | 502 | -0.566 |
| `TIMER2_MCP_Monocyte` | 0.001 | 0.986 | 0.027 | 0.54 | 0.671 | 502 | -0.468 |
| `TIMER2_MCP_Macrophage_Monocyte` | 0.001 | 0.986 | 0.027 | 0.54 | 0.671 | 502 | -0.468 |
| `TIMER2_MCP_Cancer_associated_fibroblast` | 0.014 | 0.743 | 0.041 | 0.355 | 0.519 | 502 | -0.468 |
| `TIMER2_MCP_Neutrophil` | 0.057 | 0.195 | 0.053 | 0.236 | 0.38 | 502 | 0.136 |
| `Xena_MCP_Fibroblasts` | 0.031 | 0.479 | 0.059 | 0.186 | 0.318 | 502 | -0.505 |
| `TIMER2_MCP_Myeloid_dendritic_cell` | 0.064 | 0.148 | 0.070 | 0.115 | 0.224 | 502 | -0.311 |
| `Xena_MCP_Myeloid_dendritic_cells` | 0.141 | 0.00129 | 0.153 | 0.00058 | 0.00793 | 502 | -0.323 |

### Danaher
| Feature | Unadj ρ | Unadj p | Partial ρ | Partial p | FDR (all) | n | vs purity ρ |
|---|---:|---:|---:|---:|---:|---:|---:|
| `Danaher_CD8_T_cells` | -0.110 | 0.0127 | -0.115 | 0.01 | 0.051 | 502 | -0.430 |
| `Danaher_Cytotoxic_cells` | -0.103 | 0.0189 | -0.109 | 0.0146 | 0.0597 | 502 | -0.513 |
| `Danaher_NK_CD56dim` | -0.105 | 0.0166 | -0.106 | 0.0171 | 0.0667 | 502 | -0.470 |
| `Danaher_Total_TILs_empirical` | -0.093 | 0.0355 | -0.105 | 0.0193 | 0.0667 | 502 | -0.609 |
| `Danaher_B_cells` | -0.094 | 0.0329 | -0.104 | 0.0195 | 0.0667 | 502 | -0.500 |
| `Danaher_Exhausted_CD8` | -0.099 | 0.0249 | -0.103 | 0.0214 | 0.0702 | 502 | -0.556 |
| `Danaher_Total_TILs` | -0.093 | 0.0344 | -0.100 | 0.0254 | 0.0772 | 502 | -0.591 |
| `Danaher_Th1_cells` | -0.076 | 0.0844 | -0.081 | 0.0693 | 0.154 | 502 | -0.501 |
| `Danaher_NK_cells` | -0.083 | 0.0588 | -0.072 | 0.106 | 0.216 | 502 | -0.340 |
| `Danaher_T_cells` | -0.065 | 0.143 | -0.070 | 0.117 | 0.224 | 502 | -0.533 |
| `Danaher_DC` | -0.064 | 0.146 | -0.061 | 0.173 | 0.308 | 502 | -0.516 |
| `Danaher_CD45` | -0.062 | 0.159 | -0.057 | 0.2 | 0.335 | 502 | -0.566 |
| `Danaher_Macrophages` | -0.038 | 0.387 | -0.022 | 0.625 | 0.721 | 502 | -0.486 |
| `Danaher_Neutrophils` | -0.038 | 0.388 | -0.016 | 0.713 | 0.79 | 502 | -0.411 |
| `Danaher_Treg` | 0.000 | 0.991 | -0.000 | 0.998 | 0.998 | 502 | -0.534 |
| `Danaher_Mast_cells` | 0.127 | 0.00388 | 0.148 | 0.000886 | 0.0104 | 502 | -0.331 |

## Gene coverage (Xena recomputes)

Missing genes are listed, not imputed.

### MCP-counter (Xena)
```
                      method              population  n_official  n_present                                                                                                                                                                             present_genes missing_genes                     score_column
MCP-counter (Xena recompute)                 T cells          15         14                                                                                                             CD28,CD3D,CD3G,CD5,CD6,CTLA4,FLT3LG,ICOS,MAL,PBX4,SIRPG,THEMIS,TNFRSF25,TRAT1     CHRM3-AS2                 Xena_MCP_T_cells
MCP-counter (Xena recompute)             CD8 T cells           1          1                                                                                                                                                                                      CD8B                           Xena_MCP_CD8_T_cells
MCP-counter (Xena recompute)   Cytotoxic lymphocytes           7          7                                                                                                                                                  CD8A,EOMES,FGFBP2,GNLY,KLRC3,KLRC4,KLRD1                 Xena_MCP_Cytotoxic_lymphocytes
MCP-counter (Xena recompute)               B lineage           9          8                                                                                                                                                BANK1,CD19,CD22,CD79A,CR2,FCRL2,MS4A1,PAX5          IGKC               Xena_MCP_B_lineage
MCP-counter (Xena recompute)                NK cells           9          8                                                                                                                                   CD160,KIR2DL1,KIR2DL3,KIR2DL4,KIR3DL1,NCR1,PTGDR,SH2D1B       KIR3DS1                Xena_MCP_NK_cells
MCP-counter (Xena recompute)       Monocytic lineage           7          7                                                                                                                                                  ADAP2,CSF1R,FPR3,KYNU,PLA2G7,RASSF4,TFEC                     Xena_MCP_Monocytic_lineage
MCP-counter (Xena recompute) Myeloid dendritic cells           6          5                                                                                                                                                              CD1A,CD1B,CD1E,CLEC10A,CLIC2       WFDC21P Xena_MCP_Myeloid_dendritic_cells
MCP-counter (Xena recompute)             Neutrophils          15         15                                                                                         CA4,CEACAM3,CXCR1,CXCR2,CYP4F3,FCGR3B,HAL,KCNJ15,MEGF9,SLC25A37,STEAP4,TECPR2,TLE3,TNFRSF10C,VNN3                           Xena_MCP_Neutrophils
MCP-counter (Xena recompute)       Endothelial cells          33         33 ACVRL1,APLN,BCL6B,BMP6,BMX,CDH5,CLEC14A,CXorf36,EDN1,ELTD1,EMCN,ESAM,ESM1,FAM124B,HECW2,HHIP,KDR,MMRN1,MMRN2,MYCT1,PALMD,PEAR1,PGF,PLXNA2,PTPRB,ROBO4,SOX18,SHANK3,SHE,TEK,TIE1,VEPH1,VWF                     Xena_MCP_Endothelial_cells
MCP-counter (Xena recompute)             Fibroblasts           8          8                                                                                                                                         COL1A1,COL3A1,COL6A1,COL6A2,DCN,GREM1,PAMR1,TAGLN                           Xena_MCP_Fibroblasts
```

### Danaher
```
      method      population  n_official  n_present                                          present_genes missing_genes            score_column
Danaher 2017         B-cells           9          9 BLK,CD19,FCRL2,MS4A1,KIAA0125,TNFRSF17,TCL1A,SPIB,PNOC                       Danaher_B_cells
Danaher 2017            CD45           1          1                                                  PTPRC                          Danaher_CD45
Danaher 2017 Cytotoxic cells          10         10   PRF1,GZMA,GZMB,NKG7,GZMH,KLRK1,KLRB1,KLRD1,CTSW,GNLY               Danaher_Cytotoxic_cells
Danaher 2017              DC           3          3                                    CCL13,CD209,HSD11B1                            Danaher_DC
Danaher 2017   Exhausted CD8           4          4                                LAG3,CD244,EOMES,PTGER4                 Danaher_Exhausted_CD8
Danaher 2017     Macrophages           4          4                                 CD68,CD84,CD163,MS4A4A                   Danaher_Macrophages
Danaher 2017      Mast cells           5          5                            TPSB2,TPSAB1,CPA3,MS4A2,HDC                    Danaher_Mast_cells
Danaher 2017     Neutrophils           7          7         FPR1,SIGLEC5,CSF3R,FCAR,FCGR3B,CEACAM3,S100A12                   Danaher_Neutrophils
Danaher 2017      NK CD56dim           4          4                          KIR2DL3,KIR3DL1,KIR3DL2,IL21R                    Danaher_NK_CD56dim
Danaher 2017        NK cells           3          3                                         XCL1,XCL2,NCR1                      Danaher_NK_cells
Danaher 2017         T-cells           6          6                        CD6,CD3D,CD3E,SH2D1A,TRAT1,CD3G                       Danaher_T_cells
Danaher 2017       Th1 cells           1          1                                                  TBX21                     Danaher_Th1_cells
Danaher 2017            Treg           1          1                                                  FOXP3                          Danaher_Treg
Danaher 2017     CD8 T cells           2          2                                              CD8A,CD8B                   Danaher_CD8_T_cells
```

GEP18: 18/18 present.

## Honest interpretation

1. **This is still LUAD-only.** A strong LUSC TACSTD2–immune signal (histology rework) is irrelevant here. Do not pool.
2. **Effect sizes.** Even where p is small, |ρ| near 0.1–0.2 in bulk RNA is a weak rank association in mixed tissue, not “immune desert because of TROP2”.
3. **Purity.** ESTIMATE ImmuneScore, xCell ImmuneScore, and MCP/Danaher totals are compositionally entangled with purity. The ABSOLUTE-partial number is the one that is allowed to be called “immune” rather than “not tumor”. TACSTD2 itself is orthogonal to ABSOLUTE (ρ = 0.007).
4. **Two MCP implementations.** TIMER2 ran official MCP-counter on TPM. Xena recompute is mean log2 of the same markers on RSEM. They should agree in sign; they will not match coefficient-for-coefficient. Both are shown.
5. **xCell is TIMER2, not a from-scratch Xena port.** Re-running xCell on log2(RSEM+1) would be a different assay.
6. **GEP18 ssGSEA is not the Merck assay.** Unweighted ssGSEA / z-mean of the 18 genes is the public approximation.
7. **Not protein, not ICI, not OncoSG.** TROP2 protein (the ADC target) was not measured.
8. **No causality.** Bulk correlation, even purity-adjusted, does not say TACSTD2 excludes T cells.

## Reproduce

```
pip install -r requirements.txt
python scripts/rework_A1_deconv.py
```

Downloads ~40 MB of public tables into `data/` (gitignored) on first run. Signature gene lists are bundled under `data/signatures/`.

## Files

- `REPORT.md` — this writeup
- `DEFINITIONS.md` — the definitions section, standalone
- `correlations.tsv` — unadjusted + both partial methods + CIs + FDR + feature-vs-purity
- `primary_partial.tsv` — the 15 pre-specified features
- `sample_table.tsv` — per-sample TACSTD2, purity, and all scores
- `gene_coverage.tsv` — present/missing genes per signature
- `feature_definitions.tsv` — method / role / source per column
- `summary.json` / `provenance.json`
- `figures/forest_primary_partial_rho.png`
- `figures/bar_all_partial_rho.png`
- `figures/scatter_primary.png`

## Data

- Expression: https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz
- Purity: https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5
- ESTIMATE: https://ibl.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt
- TIMER2: https://timer.cistrome.org/infiltration_estimation_for_tcga.csv.gz
