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


## Pointer

Numbers live in `REPORT.md` and `correlations.tsv`. This file is the locked definition list.
