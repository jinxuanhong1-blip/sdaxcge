# FINDING — GSE316655 CLDN4-related NHEJ / STING / IFN / APM

**Species / model:** *Homo sapiens* CD45+ FACS scRNA from **SK-MEL-5** tumors and peripheral blood in **NSG-SGM3** mice humanized with human cord-blood CD34+ cells. Cell Ranger 5.0.0, `GRCh38_and_mm10-2020-A`. Liu et al. *Sci Immunol* 2026 (PMID 41931598, DOI 10.1126/sciimmunol.adt7832). GEO **GSE316655**.

The series-design sentence also names MIA PaCa-2. Every deposited library’s extract protocol says SK-MEL-5, and the four MTX files are the SK-MEL-5 CD45+ libraries only.

This file is an additive extract on top of the myeloid re-score (same QC: 15,566 cells). It does not replace that re-score.

---

## Verdict

The matrices contain every gene in the a priori panels (NHEJ 8/8, STING 10/10, IFN 40/40, APM 21/21). **STING1 is stored as TMEM173.** A CLDN4-high versus CLDN4-low program is not estimable: **CLDN4 is 6 QC cells, each with 1 UMI**, split across T (2), B (2), and residual melanoma (2). Myeloid cells with a CLDN4 UMI: **0**. CLDN18: **0**. Mouse `Cldn4` inside QC cells: **0**. LILRB2 UMI in the six CLDN4 cells: **0**.

Those six cells do not sit together on the four scores. The two residual melanoma cells, the only tumor-cell pair, land on opposite sides of their own library × lineage distribution (NHEJ / STING / IFN / APM percentiles **17 / 45 / 8 / 53** versus **94 / 68 / 98 / 95**).

The experiment that was actually deposited is anti-LILRB2 versus isotype, **1 library versus 1**. Inside residual melanoma the four program means differ by at most 0.02. Inside tumor myeloid (50 versus 22 cells) IFN and APM means are lower in the anti-LILRB2 library; that is one library against one library.

---

## Honest n

| Item | n |
| --- | ---: |
| 10x libraries | **4** (GSM9457798–801) |
| Barcodes deposited | 30,296 |
| QC cells | **15,566** |
| Tumor anti-LILRB2 / tumor isotype QC | 6,940 / 4,261 |
| PB anti-LILRB2 / PB isotype QC | 1,449 / 2,916 |
| CLDN4 UMI>0 | **6**, each UMI = 1 |
| CLDN4+ myeloid | **0** |
| Residual tumor melanoma | 310 (153 anti / 157 isotype) |
| Tumor myeloid | 72 (50 anti / 22 isotype) |
| Biological replicates per arm | **1** |

QC matches the earlier re-score: human UMI ≥ 200, human genes ≥ 100, human UMI fraction ≥ 0.80, mito fraction ≤ 0.20. Lineage is the same marker argmax (`LYZ/CD14/…`, `CD3D/E/G`, `NKG7/GNLY`, `MLANA/PMEL/TYR`, …).

---

## Gene sets

Scores are the unweighted mean of log1p(CP10k) over genes present in the matrix. Missing genes would have been dropped; none were.

| Program | Genes | In matrix |
| --- | --- | --- |
| NHEJ | XRCC6, XRCC5, PRKDC, XRCC4, LIG4, NHEJ1, PAXX, DCLRE1C | 8/8 |
| STING | CGAS, STING1, TBK1, IKBKE, IRF3, TREX1, ENPP1, IFI16, DDX41, ZBP1 | 10/10; **STING1 → TMEM173** |
| IFN | 40-gene ISG panel used in the repo CLDN4 IFN analyses (ISG15, IFIT1, MX1, OAS2, STAT1, …) | 40/40 |
| APM | HLA-A/B/C/E/F, B2M, TAP1/2, TAPBP, PSMB8/9/10, PSME1/2, NLRC5, ERAP1/2, CALR, PDIA3, CANX, IRF1 | 21/21 |

Detection in QC cells (mean of the four libraries): B2M 99.6%, HLA-B 96.3%, HLA-A 92.9%, HLA-C 88.7%. Ku genes are common (XRCC5 51%, XRCC6 36%). Ligation genes are rare (LIG4 3.9%, XRCC4 4.8%, NHEJ1 2.9%). cGAS 13%, TMEM173 17%, IFI16 49%. **TREX1 has zero UMIs in every QC cell** even though the gene is on the reference. Canonical ISGs are mostly sparse (IFIT1 5.9%, IFI27 5.9%, RSAD2 11%).

---

## The six CLDN4 cells

Percentile is the cell’s module score inside the **same library and the same lineage**.

| Library | Lineage | Human UMIs | NHEJ | STING | IFN | APM |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| PB isotype | T | 5,495 | 95 | 69 | 7 | 28 |
| Tumor anti-LILRB2 | B | 14,617 | 72 | 57 | 78 | 97 |
| Tumor anti-LILRB2 | T | 4,923 | 99 | 40 | 49 | 73 |
| Tumor anti-LILRB2 | melanoma | 10,010 | 17 | 45 | 8 | 53 |
| Tumor isotype | B | 10,881 | 70 | 15 | 90 | 80 |
| Tumor isotype | melanoma | 18,602 | 94 | 68 | 98 | 95 |

Median human UMIs: **10,446** in the six CLDN4 cells versus **2,298** in the other 15,560 QC cells. Among CLDN4-negative cells, Spearman of log1p(human UMI) versus the NHEJ score is **0.32** overall and **0.42** inside T cells. XRCC6 is detected in 6/6 CLDN4 cells and in 31.6% of other QC cells; XRCC5 in 6/6 versus 46.0%. LIG4 is detected in 0/6. The NHEJ shift is the Ku genes in deeper libraries, which is where a single extra UMI is also more likely to land.

Descriptive median score, CLDN4 UMI>0 minus UMI=0, all QC cells (n = 6 versus 15,560): NHEJ **+0.264**, STING **−0.040**, IFN **+0.004**, APM **+0.075**. The IFN median delta is essentially zero. The NHEJ median delta is the depth-sensitive Ku detection above.

---

## Residual melanoma (the only tumor-cell split)

n = 2 CLDN4 UMI>0 versus 308 CLDN4-negative tumor melanoma cells. Median score delta: NHEJ **+0.056**, STING **+0.019**, IFN **+0.100**, APM **+0.215**. Those medians are one high cell and one ordinary cell (percentiles in the table above). Gene rows are in `tables/cldn4_gene_delta.tsv` (`compartment = tumor_melanoma`). Both positive cells have UMIs for CALR, TAPBP, HLA-A, HLA-B, HLA-E, B2M, XRCC6, and IFI16. Both are **negative** for CGAS, TREX1, LIG4, IRF1, and TAP1. TMEM173 is present in 1 of the 2 cells (background detection in melanoma is 3.9%).

Anti-LILRB2 versus isotype, residual melanoma means (153 versus 157 cells, 1 library each):

| Program | Anti mean | Isotype mean | Anti − isotype |
| --- | ---: | ---: | ---: |
| NHEJ | 0.440 | 0.456 | −0.016 |
| STING | 0.158 | 0.142 | +0.016 |
| IFN | 0.199 | 0.206 | −0.007 |
| APM | 0.745 | 0.764 | −0.019 |

---

## Tumor myeloid (no CLDN4 cell to split)

Zero myeloid QC cells have a CLDN4 UMI, so there is no myeloid CLDN4 row with a delta. Program means by the deposited antibody contrast (tumor myeloid, 50 anti versus 22 isotype, 1 library each):

| Program | Anti mean | Isotype mean | Anti − isotype |
| --- | ---: | ---: | ---: |
| NHEJ | 0.231 | 0.222 | +0.009 |
| STING | 0.145 | 0.107 | +0.038 |
| IFN | 0.257 | 0.391 | −0.133 |
| APM | 1.020 | 1.133 | −0.113 |

Whole-tumor QC means move even less, except APM (+0.069 anti minus isotype). Within tumor T cells that APM delta is −0.006, within tumor B cells −0.015, within melanoma −0.019. The whole-tumor APM difference follows the T-heavy versus B-heavy library mix already described in the myeloid re-score.

---

## What this extract is

- Public MTX/TSV only. Same four libraries as the myeloid re-score.
- Gene-level CLDN4 UMI>0 versus UMI=0 deltas for all 79 panel genes, by tissue and lineage: `tables/cldn4_gene_delta.tsv`.
- The six cells, their scores, and their within-group percentiles: `tables/cldn4_positive_cells.tsv`.
- Anti-LILRB2 versus isotype program means inside tumor compartments: `tables/treatment_program_tumor.tsv`.

## What this extract is outside of

- A CLDN4 knockdown, knockout, or CLDN4-high malignant program.
- A replicated antibody effect. Cell-level tests are not reported.
- CellChat, LIANA, or a spatial CLDN–LILRB map.
- The paper’s CLDN18.2 biochemistry or the human spatial cohorts. Those measurements are not in these four matrices.

---

## Methods

`download.py` fetches the GEO MTX/TSV. Human and mouse features are split on `GRCh38_` and `mm10___`. Scores and lineage markers use log1p(CP10k) with the human UMI total as the library size. `analyze.py` writes the tables and `figures/fig1`–`fig4`.

---

## 中文

**物种/模型：** 人源化 NSG-SGM3 小鼠皮下 SK-MEL-5，FACS 人 CD45+ 单细胞。四个 10x 文库，QC **15,566** 细胞，与已有髓系重分析同一套阈值。

**CLDN4：** 6 个细胞、各 1 条 UMI（T 2、B 2、残留黑色素瘤 2）。髓系里 CLDN4 UMI>0 为 **0**。CLDN18 为 0。这 6 个细胞的 NHEJ / STING / IFN / APM 百分位并不齐：两只残留黑色素瘤细胞一个在本库本类型的 17/45/8/53 分位，另一个在 94/68/98/95 分位。IFN 全库中位差 **+0.004**。NHEJ 中位差偏正，来自更深的细胞更容易测到 Ku（XRCC5/XRCC6）；连接酶 LIG4 在这 6 个细胞里是 0。

**抗体对照（每臂 1 个文库）：** 残留黑色素瘤四个程序均值差都在 0.02 以内。肿瘤髓系 IFN / APM 均值在 anti-LILRB2 文库较低（50 对 22 个细胞），不能当成有重复的处理效应。

面板 79 个基因都在矩阵里。STING1 在参考里的名字是 **TMEM173**。TREX1 在全部 QC 细胞中 UMI 为 0。
