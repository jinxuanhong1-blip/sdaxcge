# NHEJ-STING logic wave

## 结论

公共库里没有新的 CLDN4 敲低或过表达 RNA-seq / 芯片 / RPPA。Bitler 实验室的 HGSOC / PARPi 论文做了 OVCAR3 shCLDN4 RPPA，也做了 OVCAR8 过表达和 OVCAR3/OVCA429 CRISPRi，原始矩阵都是向通讯作者索取，GEO 上没有。

能当 STING 阳性对照的是 **GSE252340**（不是 CLDN4 实验）：H1944 用 MPS1 抑制剂 BAY1217389 强开 cGAS-STING 后，预设 STING/IFN 基因集 45/46 上调，中位 log2FC **+2.45**。同一对比里 NHEJ 基因集是平的（中位 **−0.11**）。已经编目的两个人源 CLDN4 缺失矩阵（GSE207704、GSE22493）上，NHEJ 基因集和 STING 基因集都没有按「NHEJ 下降、STING 打开」一起移动。

## Hunt

Live NCBI queries on 2026-09-21 are in `search_log.md`. The GEO CLDN4 series list matches the screen already filed in PR #116 (2026-08-16). The only public CLDN4 overexpression transcriptome is still GSE22493, and that overexpression is the control channel of the SKOV-3 siRNA array. A GEO query for CLDN4 plus RPPA returned zero series. SRA text hits are 4C-seq at the CLDN4 locus after MSH2 knockout (SRP263109).

Accepted public CLDN4-loss transcriptomes remain the three already cataloged:

| Accession | What it is |
|---|---|
| GSE50927 | Mouse lung germline Cldn4 KO. IFN opens in the uninjured lung (PR #43 / #88). Not a cancer-cell or PARPi experiment, and not re-cut here. |
| GSE207704 | MCF7 and T47D CRISPR CLDN4 knockout, collapsed FPKM. |
| GSE22493 | SKOV-3 CLDN4 siRNA versus CLDN4-overexpressing control, three two-color arrays. |

## Bitler / HGSOC / PARPi

| Paper | Perturbation | Public matrix |
|---|---|---|
| Yamamoto et al., Mol Cancer Ther 2022, PMID 35373300 | OVCAR3 shCLDN4 vs shCtrl, RPPA, n=3 vs 3. The paper reports 53BP1 and XRCC1 protein down, NHEJ activity down, homology-directed repair unchanged, and PARPi sensitization. | On request. No GEO. The supplement is TCGA differential genes, BRCA status, a drug screen, and the ex vivo tumor list. |
| Villagomez et al., Sci Rep 2025, PMID 41214101 | Stable CLDN4 overexpression in OVCAR8 and CRISPRi in OVCAR3 and OVCA429. ISRE luciferase and pSTING blots. In the text, OVCAR8 overexpression raises basal ISRE and OVCAR3 knockdown lowers it. | On request. No RNA-seq and no RPPA. |
| Villagomez et al., Cancer Res Commun 2025, PMID 39625235 | Same overexpression and knockdown models. Nuclear shape, cell cycle, olaparib colonies. | TCGA only. Other data on request. |
| Villagomez et al., Cancer Res Commun 2024, PMID 38867360 | CLDN4 CRISPRi, autophagy and micronuclei. | No RNA-seq or RPPA. |
| Breed et al., Mol Cancer Res 2019, PMID 30606772 | shCLDN4, paclitaxel apoptosis. | No transcriptome. Expression context is observational GSE18521. |
| Webb, Neville, Bitler, AACR 2019 abstract, doi 10.1158/1557-3265.ovca19-a30 | Describes OVCAR3 CLDN4-knockdown RNA-seq (EMT genes). | No accession was ever assigned. |

Two non-Bitler RNA-seq claims remain without an accession, both already listed in PR #116: Kashiwagi 2025 H1688 CRISPR CLDN4 knockout (PMID 41016339) and Zheng 2025 pancreatitis CLDN4 shRNA (PMID 40892111).

## GSE252340 is a positive STING control

Series: Kitajima, Tani, Barbie and colleagues, PMID 38227896. NCI-H1944 (KRAS/LKB1 lung line) treated with DMSO or the MPS1 inhibitor BAY1217389 for 48 h, then a 24 h washout, then RNA-seq. Deposited files are read-count workbooks (four samples). Raw library sizes are 39.3, 38.9, 43.7 and 43.1 million. Counts were summed by gene symbol and normalized by median-of-ratios. log2FC is log2((mean BAY + 1) / (mean DMSO + 1)). n=2 versus 2.

Prespecified STING/IFN panel (TREX1, cGAS, STING1, TBK1, IRF3, STAT1/2, type-I ISGs, CXCL10, CCL5, and the antigen-presentation genes listed in the script):

| Set | Present | Up | Median log2FC | Wilcoxon greater p | vs expressed background |
|---|---|---|---|---|---|
| STING/IFN panel | 46 | 45 | +2.45 | 2.6×10⁻⁹ | 2.0×10⁻²⁷ (background median +0.020, n=16,409) |
| Hallmark interferon alpha | 97 | 95 | +1.77 | 9.5×10⁻¹⁸ | 4.6×10⁻⁵³ |
| Hallmark interferon gamma | 200 | 172 | +0.86 | 4.8×10⁻³⁰ | 1.1×10⁻⁵⁷ |
| NHEJ panel | 12 | 3 | −0.11 | 0.99 | 0.98 |

Focal log2FC on the acute contrast: TREX1 +2.56, STING1 +1.03, CGAS +0.30, IFNB1 +4.59, IFIT1 +4.01, ISG15 +2.89, CCL5 +2.16, CXCL10 +6.07, TAP1 +2.23. TP53BP1 +0.16 and XRCC1 −0.21. Genome-wide Benjamini-Hochberg on the Welch tests calls **zero** genes at q<0.05, which is the expected behavior at n=2, so individual-gene FDR is not claimed. The inferential result is the prespecified set shift.

An exact label permutation has only 6 two-versus-two partitions. The acute STING median is rank 1 of those 6, so the permutation p is 0.167 and cannot go lower. That floor is reported next to the Wilcoxon p, not instead of it.

Each workbook also contains a second pair of count columns labeled BAY.R (the Japanese sheet title is the resistant-line contrast). Those columns are copied into all four GEO files and are absent from the series summary. On that extra contrast the STING panel median is −0.08 (16/46 up, permutation rank 6/6). It is not a second positive control.

CLDN4 itself on the acute contrast is log2FC +0.66. This series does not perturb CLDN4.

Figure: `figures/gse252340_acute_focal_log2fc.png`.

## NHEJ / STING slice on the two cataloged human CLDN4-loss matrices

These accessions are not new. The slice asks whether CLDN4 loss moves NHEJ mRNA down and the STING set up. GSE207704 is one collapsed FPKM value per genotype, so the log2FC is the mean of the two lines with pseudocount 0.1. GSE22493 values are the deposited sample-to-control ratios (knockdown over the CLDN4-overexpressing control); gene-level values are the mean of probes, tested with a one-sample t across three arrays.

| Matrix | CLDN4 log2FC | STING set | NHEJ set | 53BP1 | XRCC1 |
|---|---|---|---|---|---|
| GSE207704 | −0.90 (MCF7 −0.73, T47D −1.07) | 9/20 up, median −0.17 | 5/11 up, median −0.027 | +0.19 in both lines | −0.034 |
| GSE22493 | −1.23 (p=0.25, q=0.999) | 14/34 up, median −0.18 | 4/9 up, median −0.062 | −0.33 (p=0.52) | +0.69 (p=0.0009, q=0.65) |

GSE207704's cufflinks table does not contain TREX1, STING1, or TMEM173, so 20 of 48 STING-panel symbols are scored. GSE22493's oligonucleotide platform also lacks STING1/TMEM173. On GSE22493, IFIT1 is −2.03 (p=0.055), the same gene the earlier IFN audits flagged.

The Yamamoto RPPA claim (53BP1 and XRCC1 protein down after shCLDN4) has no public matrix to recompute. The two public human CLDN4-loss transcriptomes do not show that pair going down together, and they do not show the STING/IFN set going up.

## Reproduce

```bash
python3 scripts/nhej_sting_wave/analyze.py
```

The script downloads the public GEO count workbooks, GSE207704 FPKM table, GSE22493 series matrix, and GPL10555 into `/tmp/nhej_sting_cache` (override with `NHEJ_CACHE`). Outputs land in `results/nhej_sting_wave/tables/`.
