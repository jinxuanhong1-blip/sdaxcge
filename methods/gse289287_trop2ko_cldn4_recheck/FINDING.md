# FINDING — GSE289287 Trop-2 KO T-47D: CLDN4 padj recheck, IFN and APM

Additive public check on the deposited author DESeq2 table. Human gene symbol is **CLDN4** (the earlier “Cldn4” note is this gene). Not lung. Not SKB264. Does not replace the prerank GSEA in `methods/gse289287_trop2ko_gsea`.

Reproduce: `python3 methods/gse289287_trop2ko_cldn4_recheck/analyze.py`

---

## TL;DR

Trop-2 KO RNA-seq in this series exists **only as xenografts** (4 KO vs 3 WT, NRG mammary fat pad). There is no Trop-2 KO cell-line library and no `Trop2KO_cells` DESeq2 file on the GEO FTP.

**CLDN4 padj is non-significant.** That is not an FDR-only result: the author nominal p is already 0.129.

| Test on CLDN4 | Result |
|---|---|
| Author DESeq2 | log2FC **+0.284**, p **0.129**, padj **0.466** |
| Welch t on log2(norm+1), then BH across 16,258 protein-coding symbols | unadjusted p 0.041, padj **0.407** |
| Exact label permutation, all 35 assignments of 3 vs 4 | p **0.086** (3/35) |

KO control: TACSTD2 log2FC **−3.263**, padj **1.95×10⁻⁶¹**.

IFN collateral is a short author-significant ISG list, not a sample-separating score. APM has **no** gene at padj < 0.05.

---

## What “± xenografts” is in this accession

GEO series [GSE289287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289287), public 10 Feb 2026. Vacek et al., preprint doi:10.21203/rs.3.rs-6123457/v1.

| Arm | Samples on GEO | Trop-2 KO DESeq2 deposited? |
|---|---|---|
| In vitro T-47D | WT ×3 (GSM8788410–GSM8788412); DSG2 KO ×6 (GSM8788413–GSM8788418) | **No.** FTP has no Trop-2 KO cell table. The preprint RNA-seq cell top-hit list (supplementary file 5) is labeled DSG2 KO vs WT. |
| Xenograft | WT ×3 (animals 2808, 2810, 2812; GSM8788420, GSM8788421, GSM8788419). Trop-2 KO ×4 (animals 2807, 2815, 2817, 2818; GSM8788425, GSM8788422–424) | **Yes.** `GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz` |

The series summary says an in vitro Trop-2 KO comparison was done. The public sample list and the supplementary DESeq2 files do not contain it. This recheck does not invent that contrast from the DSG2 tables.

Hosts are NRG, so IFN / APM movement in the xenograft table is tumour-cell RNA, not an adaptive infiltrate.

---

## CLDN4 recheck

Primary numbers are the author columns, read from the table in this folder (not copied from the earlier note).

CLDN4 is ENSG00000189143, one protein-coding row, baseMean 746.

- log2FC **+0.284** (KO − WT), lfcSE 0.182, Wald stat +1.517
- p **0.129174**, padj **0.466374**, `significant_DE` = FALSE
- The 95% Wald interval 0.284 ± 1.96×0.182 crosses 0

Per-sample log2(author normalized count + 1):

- WT: 9.688, 9.935, 9.854
- KO: 10.078, 9.930, 10.157, 10.287

The groups overlap (one KO sample sits inside the WT range). Mean difference on this scale is +0.287.

Independent checks, same seven count columns:

- Welch t-test, two-sided, unequal variance: p = **0.041**. Benjamini–Hochberg across the 16,258 protein-coding symbols in the table: padj = **0.407**.
- Exact two-sided permutation of the 7 labels (C(7,3) = 35): p = **0.086**.

The unadjusted Welch p is the only number under 0.05. It does not survive the same multiple-testing standard the note was asking about, and the exact permutation does not clear 0.05 either. **padj stays non-significant.**

Nearby claudins in the same author table, all padj ≥ 0.05: CLDN1 +0.188 (padj 0.757), CLDN3 +0.343 (padj 0.445), CLDN7 +0.510 (p 0.0028, padj 0.057).

---

## IFN

Hallmark sets are the symbol lists in `gene_sets.json`. Competitive p is a one-sided Mann–Whitney of author Wald statistics (genes in the set vs other protein-coding genes). Sample p is an exact permutation of the mean log2(norm+1) score. With 35 labelings the smallest two-sided sample p is 1/35 ≈ 0.029.

| Set | Present | log2FC > 0 | Author padj < 0.05, up | Rank MW p | Sample-score perm p |
|---|---:|---:|---:|---:|---:|
| Hallmark IFN-γ | 186 / 200 | 106 | 5 | 8.2×10⁻⁷ | 0.66 |
| Hallmark IFN-α | 95 / 97 | 66 | 5 | 2.3×10⁻¹¹ | 0.40 |

The five author-significant IFN genes are the same in both sets, all up in Trop-2 KO:

| Gene | log2FC | p | padj |
|---|---:|---:|---:|
| IFI44L | +1.648 | 3.9×10⁻⁷ | 1.1×10⁻⁴ |
| ISG15 | +1.440 | 1.0×10⁻⁵ | 0.0012 |
| IFI44 | +1.372 | 1.1×10⁻⁴ | 0.0068 |
| IFITM3 | +1.052 | 8.3×10⁻⁴ | 0.026 |
| STAT2 | +0.711 | 1.1×10⁻³ | 0.031 |

No Hallmark IFN gene is significantly down. IFIT1 (+0.91, padj 0.16), MX1 (+0.93, padj 0.11), OAS2 (+0.81, padj 0.099), and STAT1 (+0.74, padj 0.15) are up in sign and not FDR-significant. CXCL10 is flat (−0.14, padj 0.93).

The rank test and the earlier prerank GSEA agree that the IFN sets sit toward the KO-up end of the gene list. The sample-mean score does not separate KO from WT (permutation p 0.40 and 0.66). n = 4 vs 3, and several of the called ISGs are noisy across animals (ISG15 exact permutation p = 0.11). Report the five author calls. Do not report a sample-level IFN program.

---

## APM (MHC-I antigen presentation, 21 genes)

All 21 symbols are in the table. 15 have log2FC > 0, 5 have log2FC < 0, PSMB10 is ~0 and was filtered (baseMean 0.22, no padj). **Zero genes have author padj < 0.05.**

| Gene | log2FC | p | padj |
|---|---:|---:|---:|
| HLA-F | +0.893 | 0.0043 | 0.074 |
| HLA-C | +0.699 | 0.023 | 0.196 |
| HLA-A | +0.628 | 0.071 | 0.355 |
| HLA-B | +0.562 | 0.054 | 0.312 |
| B2M | +0.539 | 0.049 | 0.300 |
| CANX | −0.382 | 0.0091 | 0.114 |
| ERAP1 | −0.640 | 0.020 | 0.180 |
| ERAP2 | −0.728 | 0.047 | 0.295 |

Competitive rank p = 0.031. Sample-score permutation p = 0.46. Direction of the classical MHC-I genes is up, and none clear FDR. This matches the earlier prerank result that the APM set was not FDR-significant.

---

## Figure

`figures/fig_cldn4_ifn_apm.png`

- A: TACSTD2 and CLDN4 per xenograft, log2(normalized count + 1)
- B: count of author padj < 0.05 genes in IFN-γ, IFN-α, and APM
- C: author log2FC ± 1.96×lfcSE. Red marks padj < 0.05

Tables: `tables/cldn4_recheck.tsv`, `key_genes.tsv`, `ifn_genes.tsv`, `apm_genes.tsv`, `geneset_summary.tsv`, `sample_inventory.tsv`.

---

## 中文

GSE289287 里 Trop-2 KO 的 RNA-seq **只有移植瘤**（4 vs 3，NRG）。GEO 没有 Trop-2 KO 细胞系文库，也没有对应的 DESeq2 表。体外文库是 WT 和 DSG2 KO。

CLDN4（人源，不是小鼠 Cldn4）作者 DESeq2：log2FC +0.284，p = 0.129，padj = 0.466。名义 p 已经不显著，不是“只被 FDR 压掉”。Welch 的未校正 p = 0.041，但在 16,258 个蛋白编码基因上 BH 后 padj = 0.407；7 个样本标签的精确置换 p = 0.086。**padj 仍然不显著。** TACSTD2 log2FC −3.263，padj 1.95×10⁻⁶¹，敲除成立。

IFN：作者表里显著上调的只有 IFI44L、ISG15、IFI44、IFITM3、STAT2（padj 1.1×10⁻⁴ 到 0.031）。基因排序检验偏 KO 上调，但样本均值的精确置换不分开（IFN-γ p = 0.66，IFN-α p = 0.40）。APM 21 个基因没有一个 padj < 0.05（B2M p = 0.049，padj = 0.30）。
