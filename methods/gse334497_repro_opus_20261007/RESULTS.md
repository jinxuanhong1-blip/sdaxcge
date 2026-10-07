# GSE334497 — reproduction and assessment of "TROP2/claudin program mediates immune exclusion to impede checkpoint blockade in breast cancer"

Prepared for the SKB264 (TROP2-ADC) + ICI lung-adenocarcinoma paper (jin xuanhong). All numbers below are computed by the scripts in `scripts/` from public files. The output tables are in `tables/` and the figures in `figures/`. Nothing was re-signed or tuned. Null and opposite results are reported as they came out. Contrast convention: **log2FC = KO vs WT (positive = higher in Trop2-KO)**. The exceptions are GSEA NES, which follow the paper's orientation (positive = enriched in Trop2-WT), and GSE241876, where the difference is responders minus non-responders.

## TL;DR

* **What GSE334497 is.** Bulk RNA-seq of 10 mouse tumours: 4T1 Trop2-KO (n=5) vs Trop2-WT (n=5) in BALB/c mice, harvested at 3 weeks.
  * **There is no anti-PD-1, hRS7 or any other treatment arm.** All 10 samples are "in vivo growth".
  * The only possible contrast is KO vs WT, so "response to checkpoint blockade" cannot be tested in this dataset.
* **Reproduced:**
  * Tacstd2 is lost in KO (log2FC −3.87, DESeq2 p=1e-12).
  * The paper's Fig 4A / S4B GSEA has the same sign and nearly the same NES for all 12 gene sets. For example, KEGG_TIGHT_JUNCTION is +1.53 here vs +1.49 in the paper, and positive regulation of cell killing is −1.82 vs −1.72.
  * Cytotoxic-lymphocyte transcripts are higher in KO: Immune core (CD8A/B, GZMB, PRF1, NKG7, IFNG) +0.91 log2, exact permutation p=0.048. Gene-level changes are +0.7 to +1.1 log2, nominal p 0.006–0.04. All 19 gene sets with within-collection GSEA FDR<0.25 are KO-enriched; 18 of them are immune sets.
* **Not reproduced, or weaker than presented:**
  1. **Nothing in the paper's GSEA passes FDR.** The paper reports nominal p only. The WT-enriched tight-junction and cell-contact sets have FDR = 1.00, and the immune sets have FDR 0.07–0.91. edgeR-QL and limma-voom find 1 gene at FDR<0.05 (Cobl).
  2. **IFN and MHC-I/APM are null.**
     * APM core: +0.22, p=0.29. IFN core: +0.27, p=0.25.
     * H2-K1 and H2-D1: +0.12 and +0.05.
     * Type-I ISGs: about 0.
     * Only Cxcl9 stands out (+1.10, padj 0.023).
  3. **Cldn4 is not significantly changed** (−0.64, p=0.20).
  4. **The Cldn1 drop (−2.48) is explained by skin carry-over.**
     * 3/5 WT and 0/5 KO tumours contain epidermis (Krt1, Lor, Flg2, Krt77, Dsg1a, Dsc3).
     * Dropping those 3 WT tumours moves Cldn1 from −2.48 (p=0.0002) to −0.53 (p=0.46).
     * Cldn7 survives this and gets stronger: −1.62, padj 0.015.
  5. **The claudin decrease is not specific.** It goes with a broader loss of epithelial differentiation (Epcam, Cdh1, Esrp1, Grhl2, Rab25, Krt7/19).
     * Epithelial-differentiation score: −1.22, p=0.024.
     * The claudin score correlates with it at r=0.94.
     * After adjusting for it, the genotype effect on Cldn1/3/4/7 is −0.07 (p=0.85).
  6. **The WT-enriched "tight junction" GSEA is partly skin and skeletal-muscle signal.**
     * Myosin heavy chains, actinins and Mylpf sit in the KEGG leading edge.
     * Removing contamination-tracking genes drops KEGG_TIGHT_JUNCTION from NES 1.53 (p=0.008) to 1.35 (p=0.11).
     * The immune sets are unaffected by this pruning.
* **Usability for the lung thesis (TROP2/claudin → immune exclusion → ICI resistance): low.**
  * At most, GSEA cites as preclinical mechanistic rationale from the original authors, with caveats.
  * It is **not** independent support:
    * same authors and same data
    * mouse TNBC, not LUAD
    * genetic KO, not an ADC
    * no ICB arm
    * n=5 vs 5 with nominal-only statistics
  * Of the three steps in the thesis, it supports only "TROP2 loss goes with lower junctional/epithelial transcripts and more cytotoxic-lymphocyte transcripts". It says nothing about checkpoint response.
  * The closest lung-specific human evidence is Bessede et al. 2024 (OAK/POPLAR, 891 NSCLC; TROP2-high predicts worse survival on atezolizumab but not on chemotherapy). That is a separate, controlled-access dataset not reanalysed here.

## 1. The paper

**Citation.** Wu B, Thant W, Bitman E, … Bardia A, Ellisen LW. *TROP2/claudin program mediates immune exclusion to impede checkpoint blockade in breast cancer.* J Immunother Cancer 2026;14(4):e012265.
* doi:10.1136/jitc-2025-012265, PMID 41932810, PMC13052784.
* Preprint: bioRxiv 2024-12-05, PMID 39677819.
* GEO: GSE334497, submitted by Bogang Wu and Leif Ellisen; public 2026-06-09; BioProject PRJNA1475339.

**Claims, mapped to figures.** Items marked † use GSE334497.

| # | Claim | Evidence in paper | Data type |
|---|---|---|---|
| 1 | Immune-cold TNBC tumour cells express apical-junction genes incl. CLDN1, CLDN7, TACSTD2 | Fig 1B–C: GeoMx, 88 ROIs, 8 treatment-naïve TNBC; cold = CD8 <100/mm² | Human spatial |
| 2 | TROP2 is the only gene shared by cold-vs-hot GeoMx and non-responders in two ICB cohorts | Fig 1D–E | Human |
| 3 | TROP2 expression correlates with tight-junction pathways (positively) and with immunity pathways and T-cell markers (negatively); it predicts poor outcome in basal BC | Fig 2A–C (METABRIC GSVA, TIMER/TCGA, KM-plotter), S2 | Human public |
| 4 | Trop2 loss slows 4T1/AT-3 growth only in immunocompetent hosts | Fig 2D–G, S3A–C | Mouse growth |
| 5 | Trop2-KO tumours have more CD3/CD4/CD8, effector-memory, PD-1⁺, GZMB⁺ and 4-1BB⁺/CD69⁺ T cells | Fig 3A–H, S3D–S (IHC/flow) | Mouse IHC/flow |
| 5† | Bulk RNA-seq deconvolution (CIBERSORT, TIMER) shows more CD8 T cells in KO | S4A | GSE334497 |
| 6 | CD8 depletion abolishes the KO growth disadvantage | Fig 3I–J | Mouse growth |
| 7† | GSEA, WT vs KO: cell–cell contact zone (NES 1.61) and tight junction (1.49) up in WT. Inflammation, immune response to tumour, cytotoxicity, PD-1 signalling, NK, Th1 cytotoxic module and cell killing (NES −1.51 to −1.72) up in KO. **Nominal p only (0.008–0.027); no FDR reported** | Fig 4A, S4B | GSE334497 |
| 8 | TROP2 co-IPs with claudin-7 (4T1, MDA-MB-468, HCC1806); KO loses membrane claudin-7 and occludin protein | Fig 4B–D, S4C, S4E–F | Mouse/human cells, IF |
| 8† | KO tumours have lower Cldn7, Cldn1 and Ocln mRNA ("normalized gene counts", n=5/5, each marked \*; test not stated) | S4D | GSE334497 |
| 9 | Cldn7 shRNA phenocopies Trop2 loss: slower growth in BALB/c but not nude mice, and more CD3 | Fig 4E–F, S5 | Mouse |
| 10 | The TROP2 ECD-TM domain (no ICD) is sufficient to restore claudin-7, tumour growth and T-cell exclusion | Fig 5 | Mouse |
| 11 | hRS7 (the naked SG antibody) + anti-PD-1 slows hTROP2-4T1 growth; neither alone does. hRS7 lowers claudin-7 and raises CD3; the combination raises TIM3⁻PD-1⁺CD8 | Fig 6, S6 | Mouse growth/IF/flow; **no RNA-seq** |
| 12 | Dox-inducible Cldn7 KD sensitises 4T1 to anti-PD-1 | S7 | Mouse growth/flow; **no RNA-seq** |
| 13 | High TACSTD2 predicts non-response to pembrolizumab. Cohort 1 (GSE241876, carboplatin/nab-paclitaxel/pembrolizumab): OR 0.06 (95% CI 0.0–0.61). Cohort 2 (Bassez 2021, EGAD00001006608, single-agent pembrolizumab, TCR expansion): OR 0.06 (0.01–0.36). Epithelial scRNA: OR 0.15 (0.02–0.84). EPCAM is not associated. Cut-points are in-sample Youden-optimal | Fig 7A, S8 | Human |

The RNA-seq methods describe 4T1 WT/KO tumour frozen sections with total RNA, rRNA depletion, NEBNext Ultra II and UMIs, analysed with GSEA desktop v4.3.2. The paper gives no n and no FDR. The GEO matrix holds 5 WT and 5 KO.

## 2. Data, QC and design problems

Input files:
* `raw/GSE334497_normalized_counts.csv.gz`: 18,047 Ensembl genes × 10 samples, the only processed file deposited. GSE334497 has no raw or count file.
* `raw/GSE334497_series_matrix.txt.gz`

MD5 checksums are in `raw/MD5SUMS`.

* **Raw counts recovered exactly.** Each column is DESeq2 counts divided by a size factor, and every value is an integer multiple of the column's minimum. So size factor = 1/min, and multiplying back gives integer counts (maximum deviation from an integer is 0). DESeq2 size factors recomputed on those counts match the deposited ones. Recovered library sizes are 8.6–19.9 M. See `tables/qc_size_factors.csv` and `raw/GSE334497_raw_counts_recovered.csv.gz`.
* **Annotation.** Ensembl 102 / GRCm38 names 98.9% of genes, and 76.2% have an MSigDB human ortholog.
* **Library wave.** Five libraries carry a `RESUB-` prefix, although GEO lists batch = 1 for all. The split is KO 4 original / 1 RESUB vs WT 1 original / 4 RESUB (Fisher p=0.21). This is a partial confound, so wave-adjusted estimates are reported alongside.
* **PCA** (top 500 variable genes, log2 CPM): no PC separates genotype (PC1 p=0.33; PC2 p=0.051) or wave. Pairwise sample correlations are 0.88–0.98.
* **RESUB-170R** is an outlier: PC1 −46, highest Hbb-bs, very high skeletal muscle (Acta1 27k, Ckm 22k normalized counts). It shares tumour number "170" with control170, but their r of 0.93 is within the range of unrelated pairs. Both are WT.
* **Tissue carry-over (post hoc; `scripts/07_tissue_contamination.py`, `figures/fig6_tissue_contamination.png`).**
  * The top DESeq2 genes include Krt1 (log2FC −25), Krtdap (−8.2) and Lypd5 (−8.0). These are suprabasal-epidermis genes that 4T1 cells do not make.
  * Rule: a tumour is skin-positive when Krt1+Flg2+Krt77+Dsg1a+Dsc3 normalized counts sum to ≥50.
  * Skin-positive tumours: control170, RESUB-170R and RESUB-169R, i.e. **3/5 WT and 0/5 KO** (Fisher p=0.17).
  * The gap between groups is clean: skin-positive tumours sum to 443–5,165, all others to 0–2.
  * Skeletal muscle (panniculus carnosus or chest wall) is present in both genotypes. It is highest in RESUB-170R.

| tumour | genotype | wave | skin flag sum | skin+ | Krt1 | Lor | Acta1 | Ckm | Tacstd2 | Cldn1 | Cldn4 | Cldn7 | Ocln |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| control170 | WT | original | 2884 | True | 343 | 72 | 1419 | 1480 | 95 | 128 | 1591 | 158 | 202 |
| RESUB-171R | WT | RESUB | 1 | False | 0 | 0 | 3885 | 1933 | 124 | 22 | 2386 | 215 | 228 |
| RESUB-170R | WT | RESUB | 443 | True | 45 | 10 | 27352 | 22248 | 66 | 53 | 1148 | 249 | 96 |
| RESUB-169R | WT | RESUB | 5165 | True | 1152 | 80 | 2 | 0 | 361 | 256 | 2843 | 159 | 289 |
| RESUB-168R | WT | RESUB | 2 | False | 0 | 0 | 1452 | 686 | 150 | 29 | 5056 | 647 | 136 |
| KO162 | KO | original | 0 | False | 0 | 0 | 84 | 49 | 13 | 15 | 2050 | 126 | 133 |
| KO164 | KO | original | 0 | False | 0 | 0 | 1620 | 782 | 1 | 3 | 464 | 106 | 86 |
| KO165 | KO | original | 0 | False | 0 | 0 | 0 | 6 | 16 | 18 | 3408 | 149 | 121 |
| KO172 | KO | original | 0 | False | 0 | 2 | 1073 | 812 | 17 | 39 | 1825 | 157 | 110 |
| RESUB-KO163R | KO | RESUB | 1 | False | 0 | 1 | 1395 | 577 | 9 | 12 | 622 | 159 | 113 |

(Values are DESeq2-normalized counts.)

**Models.**

| Model | Samples | Design | Genes at padj<0.05 |
|---|---|---|---|
| M1 (primary) | 5 KO vs 5 WT | ~genotype | 122 (26 up, 96 down in KO) |
| M2 | all 10 | ~wave + genotype | 79 |
| M3 | drop RESUB-170R | ~genotype | 159 |
| M4 | drop control170 | ~genotype | 83 |
| M5 (post hoc) | drop the 3 skin-positive WT: 2 WT (both RESUB) vs 5 KO | ~genotype | 105 |

M5 is a sensitivity check only. It is small, and its WT group is entirely RESUB wave.

edgeR quasi-likelihood (robust) and limma-voom, both after filterByExpr, each find 1 gene at FDR<0.05 (Cobl).

## 3. Gene-level results (the only contrast: Trop2-KO vs Trop2-WT)

Test: DESeq2 Wald test on recovered counts, ~genotype, n = 5 KO vs 5 WT, BH adjustment over 18,047 genes. edgeR-QL and limma-voom p-values come from unfiltered fits so that low-count genes such as Tacstd2 are included. M5 compares 2 skin-free WT vs 5 KO. Human panel genes were translated to mouse as follows (`tables/human_to_mouse_panel_notes.csv`):
* HLA-A/B/C → H2-K1, H2-D1
* CD8B → Cd8b1
* GBP1 → Gbp2b
* OAS1 → Oas1a

All per-gene, per-model numbers are in `tables/focal_genes_all_models.csv`. Plots: `figures/fig2_focal_genes_per_sample.png` and `figures/fig5_focal_gene_forest_sensitivity.png`.

| group | gene | baseMean | log2FC KO/WT | p | padj | edgeR-QL p | limma-voom p | M5 log2FC (2 skin-free WT vs 5 KO) | M5 p |
|---|---|---|---|---|---|---|---|---|---|
| TROP2 / epithelial | Tacstd2 | 85 | -3.87 | 1.2e-12 | 7.2e-09 | 1.8e-05 | 6.6e-05 | -3.66 | 8.5e-11 |
| TROP2 / epithelial | Epcam | 2136 | -1.17 | 0.0009 | 0.090 | 0.005 | 0.004 | -1.50 | 0.0005 |
| Claudins / TJ | Cldn1 | 58 | -2.48 | 0.0002 | 0.037 | 0.006 | 0.019 | -0.53 | 0.464 |
| Claudins / TJ | Cldn3 | 200 | -1.35 | 0.012 | 0.354 | 0.041 | 0.026 | -1.80 | 0.053 |
| Claudins / TJ | Cldn4 | 2139 | -0.64 | 0.197 | 0.893 | 0.261 | 0.187 | -1.15 | 0.163 |
| Claudins / TJ | Cldn7 | 213 | -1.04 | 0.006 | 0.255 | 0.013 | 0.036 | -1.62 | 3.9e-05 |
| Claudins / TJ | Ocln | 151 | -0.76 | 0.011 | 0.341 | 0.018 | 0.044 | -0.69 | 0.030 |
| Claudins / TJ | Cgn | 219 | -0.94 | 0.006 | 0.260 | 0.013 | 0.022 | -1.04 | 0.004 |
| Claudins / TJ | Tjp2 | 780 | -0.71 | 0.026 | 0.489 | 0.035 | 0.062 | -1.07 | 0.002 |
| Claudins / TJ | Marveld2 | 132 | -0.70 | 0.034 | 0.541 | 0.045 | 0.032 | -1.05 | 0.009 |
| Immune (requested core) | Cd8a | 141 | +0.70 | 0.038 | 0.569 | 0.046 | 0.067 | +0.50 | 0.231 |
| Immune (requested core) | Cd8b1 | 56 | +1.03 | 0.024 | 0.467 | 0.038 | 0.035 | +0.74 | 0.149 |
| Immune (requested core) | Gzmb | 91 | +1.14 | 0.026 | 0.487 | 0.046 | 0.129 | +1.73 | 0.008 |
| Immune (requested core) | Prf1 | 58 | +1.07 | 0.009 | 0.303 | 0.017 | 0.017 | +1.23 | 0.018 |
| Immune (requested core) | Nkg7 | 72 | +1.14 | 0.006 | 0.251 | 0.013 | 0.021 | +0.84 | 0.097 |
| Immune (requested core) | Ifng | 19 | +1.08 | 0.041 | 0.582 | 0.043 | 0.055 | +1.70 | 0.014 |
| Other lymphocyte | Ptprc | 3040 | +0.38 | 0.027 | 0.491 | 0.055 | 0.044 | +0.66 | 0.003 |
| Other lymphocyte | Cd3e | 74 | +0.43 | 0.277 | 0.943 | 0.270 | 0.261 | +0.12 | 0.811 |
| Other lymphocyte | Cd4 | 39 | +0.10 | 0.821 | 0.991 | 0.803 | 0.786 | -0.36 | 0.456 |
| Other lymphocyte | Klrd1 | 58 | +1.04 | 0.002 | 0.128 | 0.006 | 0.004 | +1.63 | 6.4e-05 |
| Other lymphocyte | Ncr1 | 14 | -0.01 | 0.988 | 0.999 | 0.994 | 0.696 | +0.38 | 0.619 |
| Other lymphocyte | Pdcd1 | 55 | +0.77 | 0.062 | 0.673 | 0.072 | 0.124 | +0.53 | 0.310 |
| Other lymphocyte | Cd274 | 763 | +0.61 | 0.089 | 0.749 | 0.096 | 0.038 | +1.13 | 0.0002 |
| APM (requested core) | B2m | 31288 | +0.43 | 0.012 | 0.348 | 0.053 | 0.065 | +0.49 | 0.061 |
| APM (requested core) | H2-K1 | 6965 | +0.12 | 0.539 | 0.989 | 0.550 | 0.509 | +0.00 | 0.996 |
| APM (requested core) | H2-D1 | 11694 | +0.05 | 0.833 | 0.991 | 0.820 | 0.777 | -0.17 | 0.590 |
| APM (requested core) | Tap1 | 1749 | +0.19 | 0.492 | 0.988 | 0.465 | 0.279 | +0.11 | 0.725 |
| APM (requested core) | Tap2 | 967 | +0.08 | 0.811 | 0.991 | 0.786 | 0.478 | -0.26 | 0.491 |
| APM (requested core) | Psmb8 | 1552 | +0.21 | 0.467 | 0.986 | 0.442 | 0.585 | +0.18 | 0.638 |
| APM (requested core) | Psmb9 | 1049 | +0.26 | 0.300 | 0.954 | 0.276 | 0.332 | +0.19 | 0.568 |
| APM (requested core) | Nlrc5 | 2694 | +0.50 | 0.199 | 0.893 | 0.214 | 0.199 | +0.66 | 0.119 |
| APM (requested core) | Tapbp | 4155 | -0.11 | 0.726 | 0.991 | 0.739 | 0.975 | -0.47 | 0.226 |
| IFN (requested core) | Stat1 | 2549 | +0.30 | 0.367 | 0.967 | 0.360 | 0.180 | +0.46 | 0.128 |
| IFN (requested core) | Irf1 | 2072 | +0.23 | 0.236 | 0.921 | 0.237 | 0.170 | +0.15 | 0.543 |
| IFN (requested core) | Cxcl9 | 1709 | +1.10 | 8.5e-05 | 0.023 | 0.001 | 0.001 | +1.55 | 5.3e-06 |
| IFN (requested core) | Cxcl10 | 246 | +0.55 | 0.041 | 0.582 | 0.042 | 0.036 | +0.78 | 0.030 |
| IFN (requested core) | Cxcl11 | 214 | +0.41 | 0.225 | 0.913 | 0.212 | 0.147 | +1.24 | 0.0008 |
| IFN (requested core) | Gbp2b | 2564 | +0.55 | 0.038 | 0.568 | 0.048 | 0.045 | +0.89 | 0.003 |
| IFN (requested core) | Gbp2 | 4100 | +0.40 | 0.105 | 0.779 | 0.118 | 0.118 | +0.62 | 0.032 |
| IFN (requested core) | Gbp4 | 3529 | +0.53 | 0.160 | 0.858 | 0.175 | 0.061 | +0.98 | 0.0010 |
| IFN (requested core) | Irf7 | 466 | -0.29 | 0.362 | 0.965 | 0.355 | 0.568 | -0.27 | 0.450 |
| IFN (requested core) | Isg15 | 280 | -0.07 | 0.841 | 0.991 | 0.855 | 0.867 | -0.01 | 0.972 |
| IFN (requested core) | Mx1 | 127 | +0.05 | 0.894 | 0.994 | 0.883 | 0.585 | +0.47 | 0.160 |
| IFN (requested core) | Oas1a | 320 | -0.33 | 0.213 | 0.901 | 0.206 | 0.232 | -0.17 | 0.596 |
| IFN (requested core) | Ifit1 | 424 | -0.01 | 0.967 | 0.999 | 0.984 | 0.987 | +0.34 | 0.317 |
| IFN (requested core) | Ifit3 | 192 | +0.12 | 0.734 | 0.991 | 0.709 | 0.548 | +0.48 | 0.194 |
| IFN (requested core) | Ido1 | 56 | -0.02 | 0.965 | 0.999 | 0.977 | 0.795 | +0.31 | 0.526 |

**Full claudin family.** Same test (DESeq2 Wald, n=5/5); M5 as above. Source: `tables/claudin_family_deseq2.csv`.

| gene | baseMean | mean norm WT | mean norm KO | log2FC KO/WT | p | padj | M5 log2FC | M5 p |
|---|---|---|---|---|---|---|---|---|
| Cldn4 | 2139 | 2605 | 1674 | -0.64 | 0.197 | 0.893 | -1.15 | 0.163 |
| Cldn12 | 451 | 490 | 412 | -0.25 | 0.104 | 0.779 | -0.28 | 0.205 |
| Cldn7 | 213 | 286 | 139 | -1.04 | 0.006 | 0.255 | -1.62 | 3.9e-05 |
| Cldn3 | 200 | 288 | 113 | -1.35 | 0.012 | 0.354 | -1.80 | 0.053 |
| Cldn5 | 93 | 81 | 106 | +0.40 | 0.324 | 0.958 | +0.14 | 0.778 |
| Cldn1 | 58 | 98 | 18 | -2.48 | 0.0002 | 0.037 | -0.53 | 0.464 |
| Cldn2 | 52 | 60 | 45 | -0.43 | 0.413 | 0.981 | -0.75 | 0.239 |
| Cldn15 | 18 | 22 | 14 | -0.66 | 0.217 | 0.903 | -0.57 | 0.409 |
| Cldn23 | 17 | 22 | 12 | -0.83 | 0.062 | 0.673 | -0.99 | 0.070 |
| Cldn20 | 13 | 13 | 13 | +0.04 | 0.941 | 0.998 | +0.31 | 0.651 |
| Cldn10 | 10 | 8 | 13 | +0.78 | 0.136 | 0.827 | +0.92 | 0.213 |

The deposited matrix (pre-filtered by the submitters) does not contain Cldn6, 8, 9, 11, 13, 14, 16, 17, 18, 19, 22, 24 or the Cldn34 cluster.

**Paper S4D (Cldn7, Cldn1, Ocln, each marked \*).** The WT and KO means plotted in S4D match the deposited normalized counts. The paper does not name the test, so several plausible ones are shown. Source: `tables/paper_S4D_style_tests.csv`.

| gene | WT group | n WT/KO | mean WT | mean KO | log2 ratio | Welch p | Student p | exact MWU p | Welch p (log2) |
|---|---|---|---|---|---|---|---|---|---|
| Cldn1 | all WT | 5/5 | 98 | 18 | -2.48 | 0.143 | 0.108 | 0.032 | 0.035 |
| Cldn1 | skin-free WT | 2/5 | 25 | 18 | -0.53 | 0.317 | 0.477 | 0.381 | 0.211 |
| Cldn7 | all WT | 5/5 | 286 | 139 | -1.03 | 0.187 | 0.152 | 0.016 | 0.091 |
| Cldn7 | skin-free WT | 2/5 | 431 | 139 | -1.63 | 0.405 | 0.053 | 0.095 | 0.317 |
| Ocln | all WT | 5/5 | 190 | 113 | -0.76 | 0.085 | 0.057 | 0.095 | 0.078 |
| Ocln | skin-free WT | 2/5 | 182 | 113 | -0.69 | 0.366 | 0.053 | 0.095 | 0.311 |

What the S4D tests show:
* On the plotted normalized counts, a t-test (Welch or Student) gives p ≥ 0.057 for all three genes.
* The exact Mann–Whitney test gives Cldn7 p=0.016, Cldn1 p=0.032 and Ocln p=0.095.
* DESeq2 gives smaller p-values: 0.006, 0.0002 and 0.011.
* So the three asterisks hold for Cldn7 and Cldn1 under some tests. Ocln at p<0.05 is reproduced only by DESeq2/edgeR/limma, not by tests on the plotted values.
* Within WT, Cldn1 follows the epidermis score (Spearman ρ=0.9, n=5).

## 4. Panel scores

* **Score:** per-sample mean of log2(normalized+1) over panel genes. The KO−WT difference of this score equals the average log2FC across the panel.
* **Tests:** Welch t; exact Mann–Whitney; exact permutation over all C(10,5)=252 labelings; OLS score ~ genotype + wave; leave-one-out range.
* **Gene-set tests on voom:** limma fry (self-contained) and camera (competitive, with inter-gene correlation estimated). Camera's default fixed correlation of 0.01 gave p-values down to 3e-9; those are anti-conservative and are listed only in the CSV.
* **Skin-free comparison:** 2 WT vs 5 KO, exact permutation over 21 labelings. The smallest attainable p is 0.048.
* Panels were fixed before testing (`scripts/gene_sets.py`). The "requested" panels are the three given in the brief: APM, IFN and immune core.

Sources: `tables/panel_score_tests.csv`, `tables/panel_tests_limma_fry_camera.csv`, `tables/panel_score_tests_skin_free_WT.csv`, `figures/fig3_panel_scores.png`.

| panel | genes | diff KO-WT | Welch p | exact MWU p | exact perm p (252) | wave-adj diff | wave-adj p | LOO range | fry p (FDR) | camera p, est. cor (FDR) | skin-free diff (2v5) | skin-free perm p (21) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| APM_core | 9 | +0.22 | 0.282 | 0.421 | 0.294 | +0.31 | 0.243 | +0.15 to +0.35 | 0.280 (0.280) | 0.289 (0.289) | +0.10 | 0.810 |
| IFN_core | 15 | +0.27 | 0.241 | 0.222 | 0.254 | +0.43 | 0.149 | +0.17 to +0.44 | 0.241 (0.263) | 0.234 (0.255) | +0.49 | 0.095 |
| Immune_core | 6 | +0.91 | 0.041 | 0.056 | 0.048 | +1.47 | 0.003 | +0.65 to +1.09 | 0.037 (0.064) | 0.041 (0.083) | +0.92 | 0.095 |
| CD8_T | 5 | +0.79 | 0.075 | 0.095 | 0.071 | +1.40 | 0.005 | +0.58 to +0.96 | 0.078 (0.104) | 0.086 (0.115) | +0.59 | 0.333 |
| Effector | 6 | +0.70 | 0.026 | 0.032 | 0.024 | +1.08 | 0.003 | +0.55 to +0.87 | 0.015 (0.060) | 0.026 (0.079) | +0.87 | 0.048 |
| NK | 6 | +0.60 | 0.019 | 0.056 | 0.032 | +0.69 | 0.035 | +0.50 to +0.74 | 0.013 (0.060) | 0.028 (0.079) | +0.87 | 0.048 |
| Leukocyte | 8 | +0.36 | 0.048 | 0.095 | 0.056 | +0.38 | 0.095 | +0.29 to +0.45 | 0.025 (0.064) | 0.033 (0.079) | +0.40 | 0.143 |
| Checkpoint | 7 | +0.54 | 0.040 | 0.032 | 0.040 | +0.75 | 0.025 | +0.46 to +0.68 | 0.036 (0.064) | 0.060 (0.094) | +0.45 | 0.143 |
| MHC_II | 7 | +0.37 | 0.139 | 0.151 | 0.119 | +0.52 | 0.100 | +0.25 to +0.51 | 0.100 (0.120) | 0.128 (0.154) | +0.35 | 0.333 |
| Claudin_1_3_4_7 | 4 | -1.31 | 0.031 | 0.016 | 0.016 | -1.29 | 0.086 | -1.53 to -0.96 | 0.027 (0.064) | 0.024 (0.079) | -1.40 | 0.190 |
| Paper_TJ_trio (Cldn1/Cldn7/Ocln) | 3 | -1.22 | 0.005 | 0.008 | 0.008 | -1.28 | 0.019 | -1.36 to -1.00 | 0.002 (0.024) | 0.007 (0.079) | -0.98 | 0.095 |
| TJ_claudin_program | 18 | -0.76 | 0.039 | 0.095 | 0.071 | -0.85 | 0.076 | -0.91 to -0.61 | 0.060 (0.091) | 0.063 (0.094) | -0.93 | 0.095 |

Panel definitions (mouse symbols):
* **CD8_T:** Cd8a, Cd8b1, Cd3e, Cd3d, Cd3g
* **Effector:** Gzma, Gzmb, Gzmk, Prf1, Nkg7, Ifng, Fasl
* **NK:** Ncr1, Klrb1c, Klrd1, Klrk1, Klrc1, Il2rb
* **Leukocyte:** Ptprc, Coro1a, Laptm5, Lcp1, Cd53, Itgb2, Cd48, Lcp2
* **Checkpoint:** Pdcd1, Cd274, Pdcd1lg2, Ctla4, Lag3, Havcr2, Tigit
* **TJ program:** Cldn1/3/4/7, Ocln, Tjp1–3, F11r, Marveld2/3, Cgn, Cgnl1, Crb3, Pard3, Pard6b, Patj, Llgl2

Reading the panel results:
* After BH correction over the 12 panels, only the paper's own trio (Cldn1/Cldn7/Ocln) has fry FDR <0.05 (0.024). Every other panel has FDR 0.06–0.29.
* All leave-one-out signs are stable.
* The wave adjustment strengthens the immune panels. The confound runs against KO (KO is mostly original wave, WT mostly RESUB).
* The ratio CD8_T − Leukocyte is +0.43 (p=0.31; wave-adjusted p=0.035). There is no clean evidence that the CD8 share of leukocytes rises, as opposed to leukocyte content as a whole.

**Post hoc epithelial-state diagnostics** (`tables/tj_scores_adjusted_for_epithelial_state.csv`, `figures/fig4_claudin_vs_epithelial_state.png`):

| diagnostic (post hoc) | diff KO-WT | Welch p | exact perm p | wave-adj diff | wave-adj p |
|---|---|---|---|---|---|
| Epithelial differentiation (Cdh1, Esrp1/2, Grhl2, Ovol1, Rab25, St14, Krt7, Krt19) | -1.22 | 0.019 | 0.024 | -1.27 | 0.056 |
| Epithelial content (Krt8, Krt18) | -0.36 | 0.110 | 0.095 | -0.15 | 0.545 |
| Basal keratins (Krt5, Krt14) | -3.24 | 0.014 | 0.008 | -3.79 | 0.025 |
| Mesenchymal (Vim, Zeb1, Snai2, Twist1, Cdh2) | +0.04 | 0.817 | 0.833 | +0.21 | 0.280 |
| Proliferation (Mki67, Top2a, Ccnb1, Cdk1, Birc5, Mcm2) | -0.19 | 0.407 | 0.421 | -0.25 | 0.432 |

In OLS score ~ genotype + covariate:

| Claudin score | Covariate | r with covariate | Genotype coefficient | p |
|---|---|---|---|---|
| Claudin_1_3_4_7 | epithelial differentiation | 0.94 | −0.07 | 0.85 |
| TJ program | epithelial differentiation | 0.99 | +0.13 | 0.15 |
| Paper trio | epithelial differentiation | 0.87 | −0.57 | 0.14 |
| Claudin_1_3_4_7 | Krt8/18 | — | −1.35 | 0.063 |
| Paper trio | Krt8/18 | — | −1.36 | 0.0098 |

So the claudin drop is not explained by less tumour epithelium (Krt8/18), but it cannot be separated from a coordinated loss of the epithelial-differentiation program. Basal keratins Krt5/14 are partly skin-derived: within WT, Krt5 follows the skin score with ρ=0.9.

## 5. GSEA reproduction (paper Fig 4A / S4B)

Method (`scripts/04_gsea.py`):
* **Emulation of GSEA desktop:**
  * input: the deposited normalized counts
  * gene IDs: mouse Ensembl → human via the MSigDB 2024.1 chip, collapsed with Max_probe
  * metric: Signal2Noise; weighted ES (p=1)
  * set sizes 15–500
  * null: the **exact phenotype-permutation null of all 252 labelings**. GSEA's 1,000 random permutations can only resample these 252.
  * FDR: GSEA's NES-based formula, computed within each MSigDB sub-collection (as if each GMT were run separately) and pooled.
* **Contamination pruning (`scripts/07_tissue_contamination.py`):** set members whose expression tracks the skin or muscle score across all 10 tumours (Spearman ρ ≥ 0.8; 295 genes) are removed, the ranking is kept, and the analysis is rerun.
* **Skin-free run:** 2 WT vs 5 KO with 21 labelings. NES is reported; p-values at this resolution carry little information.

Sources: `tables/gsea_paper_fig4A_reproduction.csv`, `tables/gsea_paper_sets_skin_free_WT.csv`, `tables/gsea_paper_sets_leading_edge.csv`, `figures/fig1_gsea_paper_vs_reproduction.png`.

| paper label | MSigDB set | n | paper NES (p) | repro NES (p) | FDR within collection | FDR pooled | NES (p), contamination-tracking genes pruned | NES skin-free WT (2v5) |
|---|---|---|---|---|---|---|---|---|
| Cell_cell contact zone | GOCC_CELL_CELL_CONTACT_ZONE | 66 | +1.61 (0.008) | +1.57 (0.032) | 1.00 | 1.00 | +1.33 (0.152); −14 genes | +1.55 |
| Tight junction | KEGG_TIGHT_JUNCTION | 107 | +1.49 (0.009) | +1.53 (0.008) | 1.00 | 1.00 | +1.35 (0.112); −8 genes | +1.51 |
| Tight junction (alt.) | GOCC_TIGHT_JUNCTION | 111 | +1.49 (0.009) | +1.52 (0.056) | 1.00 | 1.00 | +1.35 (0.152); −5 genes | +1.39 |
| Leukocyte degranulation | GOBP_LEUKOCYTE_DEGRANULATION | 67 | -1.51 (0.015) | -1.54 (0.048) | 0.56 | 0.61 | -1.54 (0.048); −0 | -1.24 |
| Inflammation pathway | BIOCARTA_INFLAM_PATHWAY | 19 | -1.55 (0.025) | -1.56 (0.056) | 0.30 | 0.56 | -1.56 (0.056); −0 | -1.61 |
| Leukocyte mediated cytotoxicity | GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY | 110 | -1.62 (0.025) | -1.70 (0.040) | 0.41 | 0.46 | -1.73 (0.032); −2 | -1.56 |
| Immune response to tumor cell | GOBP_IMMUNE_RESPONSE_TO_TUMOR_CELL | 24 | -1.63 (0.014) | -1.69 (0.040) | 0.39 | 0.45 | -1.69 (0.040); −0 | -1.67 |
| Lymphocyte costimulation | GOBP_LYMPHOCYTE_COSTIMULATION | 39 | -1.64 (0.017) | -1.60 (0.056) | 0.49 | 0.52 | -1.60 (0.056); −0 | -1.37 |
| PD1 signaling | REACTOME_PD_1_SIGNALING | 16 | -1.64 (0.027) | -1.62 (0.040) | 0.91 | 0.50 | -1.62 (0.040); −0 | -1.25 |
| NK cell mediated cytotoxicity | KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY | 97 | -1.65 (0.016) | -1.69 (0.032) | 0.07 | 0.46 | -1.69 (0.032); −0 | -1.34 |
| T cell mediated cytotoxicity | GOBP_T_CELL_MEDIATED_CYTOTOXICITY | 44 | -1.66 (0.025) | -1.71 (0.024) | 0.52 | 0.49 | -1.71 (0.024); −0 | -1.77 |
| Th1 cytotoxic module | BOSCO_TH1_CYTOTOXIC_MODULE | 89 | -1.68 (0.018) | -1.69 (0.024) | 0.55 | 0.45 | -1.69 (0.024); −0 | -1.63 |
| Positive regulation of cell killing | GOBP_POSITIVE_REGULATION_OF_CELL_KILLING | 56 | -1.72 (0.015) | -1.82 (0.024) | 0.72 | 0.87 | -1.82 (0.024); −0 | -1.78 |

* **Sign and size are reproduced for all 12 paper sets.** KEGG_TIGHT_JUNCTION is the best match for the paper's "Tight junction" (NES 1.53, p=0.008, vs the paper's 1.49 and 0.009). Feeding GSEA log2-transformed input instead gives the same signs (CSV column `repro_NES_log2input`).
* **Significance is not established.**
  * Across 9,387 sets in 9 collections, 19 have within-collection FDR <0.25. **All 19 are KO-enriched, and 18 are immune sets** (e.g. KEGG allograft rejection, NK cytotoxicity, antigen processing; REACTOME immunoregulatory interactions, FDR 0.085).
  * The 19th, and the strongest (FDR 0.011), is KORKOLA_YOLK_SAC_TUMOR_UP. This is a 20-gene chromosome-12p germ-cell-tumour signature, not an immune set. It is driven mostly by Ccnd2 (+0.85, padj 7e-05), with small shifts in Slc2a3, Tpi1 and Gapdh.
  * **No WT-enriched set has FDR <0.25.** The paper's three WT-enriched sets have FDR 1.00: the positive tail is no heavier than the permutation null.
  * No GO BP/CC/MF set reaches FDR <0.25.
  * Pooled across all collections, one set passes FDR <0.25 (KORKOLA_YOLK_SAC_TUMOR_UP, KO-enriched, 0.060).
* **Leading edges of the WT-enriched sets** (`tables/gsea_paper_sets_leading_edge.csv`):
  * KEGG_TIGHT_JUNCTION (53 genes) contains skeletal-muscle sarcomere genes: Myh1, Myh2, Myh4, Myh7, Actn2, Actn3, Myl2 and Mylpf (log2FC −2.7 to −3.7, falling to −0.5 to −2.2 without the skin-positive WT). It also contains Cldn1, plus genuine epithelial junction genes (Cldn3/4/7, Ocln, Cgn, Pard6b, Llgl2, Tjp2, Myh14).
  * GOCC_CELL_CELL_CONTACT_ZONE (31 genes) is dominated by desmosome genes (Dsc2, Dsg2, Dsp, Jup, Pkp2/4) and muscle/sarcolemma genes (Myh1, Nrap, Des, Cav3, Scn4b, Scn1b, Kcnj11, Ankrd23, Fxyd1).
  * GOCC_TIGHT_JUNCTION includes Pof1b. Pof1b is −4.21 in M1 but +0.03 without the skin-positive WT, i.e. it is entirely skin-derived.
* **IFN/APM context sets** (desktop emulation, all 10 samples):

| Set | NES | p |
|---|---|---|
| HALLMARK_INTERFERON_GAMMA_RESPONSE | −1.45 | 0.17 |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | −1.14 | 0.36 |
| REACTOME_INTERFERON_GAMMA_SIGNALING | −1.22 | 0.26 |
| GO MHC-I antigen presentation | −1.12 | 0.31 |
| HALLMARK_ALLOGRAFT_REJECTION | −1.75 | 0.040 |
| HALLMARK_APICAL_JUNCTION | +1.16 | 0.21 (0.65 after pruning) |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | −0.93 | 0.55 |

* **Mouse-panel GSEA** (phenotype permutation, `tables/gsea_mouse_panels.csv`):

| Panel | NES | p |
|---|---|---|
| APM | −1.29 | 0.24 |
| IFN | −1.36 | 0.22 |
| Immune core | −1.46 | 0.008 |
| Effector | −1.64 | <0.008 (0 of 125 same-sign null labelings) |
| NK | −1.75 | 0.008 |
| TJ program | +1.59 | 0.024 |

* **Preranked GSEA on the DESeq2 Wald statistic** (gene-set permutation; columns `*_preranked_*`) gives p between 1e-4 and 0.02 for nearly every set examined, including E2F and G2M (both "WT-enriched", p≈1e-4). That null ignores inter-gene correlation, so those p-values are not used anywhere in this report.

## 6. What the paper's other RNA-linked claims look like here

* **S4A (CIBERSORT/TIMER: more CD8 T cells in KO).** CIBERSORT LM22 is a human signature and was not rerun. The marker-level results agree in direction but are only borderline:
  * CD8_T panel +0.79 (exact p=0.071); Cd8a +0.70 (p=0.038); Cd8b1 +1.03 (p=0.024)
  * Cd3e +0.43 (p=0.28); Cd4 +0.10 (p=0.82)
* **"Inflammatory / anti-tumour immune programs increased in KO" (GEO summary).** Lymphocyte-cytotoxicity transcripts are up. Interferon-response and MHC-I/APM transcripts are not (see sections 3–5). Cxcl9 is the one IFN-γ-inducible gene that rises clearly (+1.10, padj 0.023).
* **Cldn7 (the paper's mechanistic claudin).** Lower in KO: −1.04 in M1 (p=0.006, padj 0.25) and −1.62 without the skin-positive WT (padj 0.015). This is the most robust claudin result.
* **Cldn1.** Explained by epidermis in 3 WT tumours (section 2). Without them: −0.53, p=0.46.
* **Ocln.** −0.76 (p=0.011); −0.69 without skin (p=0.030).
* **Cldn4.** Not significant in any model: −0.64 (p=0.20); M2 −0.96 (p=0.23); M3 −0.82 (p=0.18); M4 −0.78 (p=0.15); M5 −1.15 (p=0.16).
* **Not testable in GSE334497** (no RNA-seq exists for these): CD8-depletion dependence, the Cldn7-KD phenocopy, ECD-TM rescue, hRS7 ± anti-PD-1 and inducible Cldn7-KD + anti-PD-1.

## 7. Related datasets and literature (same authors / same topic)

Accessions and PMIDs were all checked against NCBI on 2026-10-07.

**Checked numerically here:**

1. **GSE289287** (Vacek, Kurucová, Tichý, Souček; "Trop-2 governs anti-metastatic desmosomal integrity"; no PMID yet).
   * Design: TACSTD2-KO vs WT **T-47D (ER⁺ luminal) xenografts in NRG mice** (no T, B or NK cells). Submitters' DESeq2 table: 3 WT vs 4 KO.
   * Use: tests whether TROP2 loss lowers claudins cell-intrinsically, without adaptive immunity.
   * Result:
     * TACSTD2 −3.26 (padj 2e-61)
     * **CLDN7 +0.51 (p=0.003, padj 0.057), i.e. opposite to 4T1**
     * CLDN1 +0.19 (p=0.43), CLDN3 +0.34 (p=0.12), CLDN4 +0.28 (p=0.13), EPCAM +0.09
     * OCLN +0.47, but barely expressed (baseMean 2)
     * CGN −0.74 (p=0.003), F11R −0.31 (p=0.036)
   * Interpretation: in this line and host, TROP2 loss does not lower claudins. The 4T1 claudin decrease is therefore either specific to context or the immune environment, or downstream of the epithelial-state change; the data cannot tell these apart. Source: `tables/related_GSE289287_trop2ko_xenograft_tj_genes.csv`.
2. **GSE241876** (Wilkerson … Montero, Clin Cancer Res 2024, PMID 37882661). This is the paper's ICB "Cohort 1": metastatic TNBC treated with carboplatin + nab-paclitaxel + pembrolizumab.
   * Samples: 15 pre-treatment biopsies, 8 responders (CR/PR) vs 7 non-responders (SD/PD), the paper's definition. log2 CPM.
   * Tests: two-sided Mann–Whitney, plus the paper's approach of an in-sample Youden-optimal cut-point and an odds ratio (Haldane-corrected). New here: an exact permutation p that re-runs the cut-point search on all 6,435 relabelings. Source: `tables/related_GSE241876_pembro_response_tj_genes.csv`.

| gene | median log2CPM R (n=8) | median log2CPM NR (n=7) | mean diff R−NR | MWU p | Youden split | OR (Haldane) | naive Fisher p | cut-point-adjusted perm p |
|---|---|---|---|---|---|---|---|---|
| TACSTD2 | 7.43 | 9.38 | -1.65 | 0.072 | high:0R/5NR low:8R/2NR | 0.027 | 0.007 | 0.033 |
| EPCAM | 7.76 | 8.57 | -0.65 | 0.613 | high:2R/4NR low:6R/3NR | 0.299 | 0.315 | 0.736 |
| CLDN1 | 2.37 | 4.03 | -0.96 | 0.536 | high:2R/4NR low:6R/3NR | 0.299 | 0.315 | 0.736 |
| CLDN3 | 7.14 | 7.32 | +0.04 | 0.779 | high:7R/4NR low:1R/3NR | 3.889 | 0.282 | 0.799 |
| CLDN4 | 8.44 | 9.34 | -0.74 | 0.232 | high:1R/4NR low:7R/3NR | 0.156 | 0.119 | 0.373 |
| CLDN7 | 7.76 | 7.84 | -0.58 | 0.694 | high:1R/3NR low:7R/4NR | 0.257 | 0.282 | 0.799 |
| OCLN | 4.5 | 5.29 | +0.34 | 0.955 | high:3R/5NR low:5R/2NR | 0.289 | 0.315 | 0.660 |
| CD8B | 2.66 | 1.7 | +1.21 | 0.072 | high:5R/1NR low:3R/6NR | 6.81 | 0.119 | 0.252 |
| CXCL9 | 4.3 | 2.57 | +1.84 | 0.152 | high:5R/1NR low:3R/6NR | 6.81 | 0.119 | 0.252 |
| PDCD1 | 2.32 | 1.47 | +1.14 | 0.148 | high:5R/1NR low:3R/6NR | 6.81 | 0.119 | 0.252 |

   * TACSTD2 is lower in responders and the in-sample split is extreme (0/5 responders among TACSTD2-high). The OR of 0.027 is in line with the paper's 0.06.
   * Accounting for the cut-point search, p=0.033 (naive 0.007). The rank test gives p=0.072. This is real but fragile with n=15.
   * **No claudin is associated with response** (CLDN4 p=0.23, CLDN7 p=0.69).
   * The paper's Cohort 2 (Bassez 2021, EGAD00001006608) is controlled access and was not reanalysed.

**Same group: Ellisen lab, MGH**
* GSE118389 / GSE118390 (Karaayvaz … Michor, Ellisen; Nat Commun 2018; PMID 30181541): TNBC scRNA-seq with no ICB or TROP2-immune contrast.
* bioRxiv preprint of this paper: PMID 39677819.
* The paper's ref. 5 is the group's collagen-alignment immune-exclusion work. No GEO accession was located for it.

**TROP2 and ICB / immunity, human**

| Reference | Setting | Finding | Direction for the thesis | Data |
|---|---|---|---|---|
| Bessede … Italiano, Clin Cancer Res 2024 (PMID 38048058) | NSCLC, OAK/POPLAR trials, n=891 | TROP2-high goes with worse PFS/OS on atezolizumab but not chemotherapy, and with lower T-cell infiltration; confirmed by mIF and plasma proteomics | **Supportive, lung-specific, human. The most relevant external evidence** | Controlled; not reanalysed |
| Yang … Chuang, Am J Cancer Res 2025 (PMID 40520853) | LUAD | TACSTD2 vs CD8 infiltration ρ = −0.11 (P=0.014); TACSTD2 associated with worse outcome in immunotherapy-treated patients | Weakly supportive; very small effect | — |
| Khan … Raphael, Transl Cancer Res 2025 (PMID 40104697) | Commentary on Bessede ("TROP2 paradox") | — | Context only | — |
| Brust … Linxweiler, Front Oncol 2026 (PMID 42625594) | HNSCC, n=47 treated with anti-PD-1 (IHC) | TROP2 **not** associated with PD-1 response; in TCGA, associated with several immune populations | **Null / contrary** in another squamous epithelial cancer | — |

**Claudins / tight junctions and immunity**

| Reference | Finding | Direction for the thesis | Data |
|---|---|---|---|
| De Sanctis … Bronte, Immunity 2024 (PMID 38749447; GSE230868, mouse, n=6) | Claudin-18 on PDAC cells **promotes** T-cell infiltration and anti-tumour immunity | **Opposite** for a claudin | GEO |
| Liu … Zhang, Sci Immunol 2026 (PMID 41931598; GSE316655, n=4) | Claudins bind LILRB inhibitory receptors and drive myeloid immunosuppression | Claudin → immunosuppression via a different (non-barrier) mechanism | GEO |
| Yu … Yin, Front Oncol 2025 (PMID 40951359) | CLDN18.2 defines an intrahepatic cholangiocarcinoma subgroup with CD8 exclusion | Supportive for claudin–exclusion; IHC only | — |
| Liang … Ding, Cell Death Dis 2025 (PMID 41102168) | Claudin-7 **deficiency** reprograms neutrophils in colorectal cancer | Opposite direction for Cldn7 loss | — |
| GSE273512 (mouse, n=12, no PMID) | "Claudin 7 suppresses invasion and metastasis through repression of a smooth muscle actin program" | No immune read-out | GEO |

**Epithelial state / EMT and immunity (relevant to the epithelial-program confound)**

| Reference | Finding | Direction for the thesis | Data |
|---|---|---|---|
| Dongre … Weinberg, Cancer Discov 2021 (PMID 33328216; GSE161746/GSE161748) | **Mesenchymal**, not epithelial, carcinoma cells are the more immunosuppressive and ICB-resistant state | Contrary to "epithelial tight junctions exclude T cells" | GEO |
| GSE155577 (mouse, n=26, no PMID) | Epithelial vs mesenchymal tumour cells from E-A14 tumours grown in immunodeficient vs immunocompetent hosts | Could test immune-dependent epithelial-state shifts; not analysed | GEO |
| Vick … Serody, J Immunol 2021 (PMID 34607937; GSE182778) | Claudin-low (T11) TNBC model; anti-PD-1 enhances Treg function | "Claudin-low" is a subtype label, not a claudin mechanism | GEO |

## 8. Usability verdict for the lung thesis (TROP2/claudin → immune exclusion → resists ICI)

**What GSE334497 shows** (n=5 vs 5; 4T1 TNBC in BALB/c; nominal statistics):

1. Genetic Trop2 loss (Tacstd2 −3.9 log2) goes with lower epithelial-junction transcripts:
   * Cldn7 −1.0 to −1.6; Cldn3 −1.35; Ocln −0.76; Cgn −0.94; Epcam −1.17
   * Paper trio score −1.22 (exact p=0.008; fry FDR 0.024)
2. It goes with higher cytotoxic-lymphocyte transcripts:
   * Immune core +0.91 (p=0.048); Effector +0.70 (p=0.024); NK +0.60 (p=0.032); Cxcl9 +1.10 (padj 0.023)
   * These are the most consistent signals genome-wide: all 19 sets at FDR <0.25 are enriched in KO, and 18 of them are immune sets.
3. The paper's Fig 4A GSEA directions and effect sizes are reproducible from the deposited matrix.

**What it does not show:**

1. **Nothing about checkpoint blockade.** There is no anti-PD-1 arm, so "resists ICI" has no transcriptomic evidence here. The paper's ICB evidence is tumour growth in hTROP2-4T1 with hRS7 + anti-PD-1 and S7 (no RNA-seq), plus human cohorts (GSE241876: TACSTD2 effect real but fragile, n=15; claudins null).
2. **Nothing about lung.** The model is mouse 4T1 mammary carcinoma. The paper itself only speculates about lung.
3. **Nothing about ADC therapy or payloads.** The perturbation is a germline CRISPR KO in the tumour line. hRS7 in the paper is the naked antibody, not an ADC. SKB264's (sacituzumab tirumotecan's) TOP1-inhibitor payload effects, such as STING/IFN induction or immunogenic cell death, are not addressed.
4. **Causality from tight junctions to exclusion.** A single genetic contrast cannot separate "claudin barrier" from the broader loss of epithelial differentiation that comes with Trop2 KO in this line (claudin score r=0.94 with that program; adjusted genotype effect −0.07, p=0.85). Two further points weaken the cell-intrinsic link: in GSE289287 (T-47D, no immune system) TROP2 loss does not lower claudins and CLDN7 goes up, and in GSE241876 claudins are not associated with ICB response.
5. **An inflamed (IFN/MHC-I) conversion.** APM and IFN cores are null (+0.22, +0.27; p≈0.25–0.29), as are H2-K1/H2-D1 and type-I ISGs. The immune change looks like modestly more lymphocytes, not an interferon-driven "hot" conversion.
6. **The specific claudins used in the narrative.**
   * Cldn1: skin artefact.
   * Cldn4: never significant.
   * Cldn7 is the only claudin result that survives every model.
7. **Robust statistics.**
   * No gene set passes FDR in the direction claimed for the WT side.
   * edgeR and limma find 1 gene at FDR <0.05.
   * The tight-junction GSEA drops to p=0.11–0.15 once skin- and muscle-tracking genes are pruned.
8. **Independence.** It is the same authors' data supporting the same authors' claims. Reanalysing it does not create independent replication.

**How to use it in the lung paper.** Cite Wu et al. 2026 (with GSE334497) only as *preclinical mechanistic rationale from TNBC*, with the caveats. For example:

> "In a murine TNBC model, Trop2 deletion was accompanied by reduced tight-junction (notably claudin-7) and increased cytotoxic-lymphocyte transcripts (Wu et al. 2026); in our reanalysis of the public data (GSE334497, n=5/group) these differences were nominally significant without changes in interferon or MHC-I programs, and no checkpoint-treated arm was profiled."

Do **not** present it as independent validation, as lung evidence, or as evidence about ICB response. For lung-specific support, Bessede 2024 (OAK/POPLAR) is the relevant human dataset, but it is controlled access and was not reanalysed here.

## 9. Reproduce

The scripts were run with Python 3.12.3 and R 4.3.3 (DESeq2 1.42.0, edgeR 4.0.16, limma 3.58.1). Package versions are in `requirements.txt` and `tables/sessionInfo_R.txt`.

```bash
cd methods/gse334497_repro_opus_20261007/scripts && bash run_all.sh
```

| Script | What it does |
|---|---|
| `00_fetch.sh` | GEO, paper, MSigDB 2024.1, Ensembl 102 GTF. References go to `cache/` (git-ignored) |
| `01_prepare_qc.py` | Raw-count recovery, annotation, QC |
| `02_de.R` | DESeq2 M1–M5, edgeR-QL, limma-voom, fry/camera |
| `03_scores.py` | Panel scores and exact tests |
| `04_gsea.py` | GSEA emulation and preranked GSEA |
| `05_related_datasets.py` | GSE289287, GSE241876 |
| `07_tissue_contamination.py` | Skin/muscle diagnostics, skin-free and pruned GSEA |
| `06_figures.py` | Figures |

Everything is deterministic. The only random step is the preranked gene-set null, which uses a fixed seed.

Analyses defined after seeing the data are labelled post hoc: the epithelial-state panels, the skin/muscle panels and flag, M5, and the pruned GSEA. The primary panels and the paper-set mapping were fixed before any test was run.
