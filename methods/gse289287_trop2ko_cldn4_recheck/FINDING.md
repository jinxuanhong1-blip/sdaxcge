# FINDING — GSE289287 Trop-2 KO T-47D: IFN holds; STING and NHEJ do not

Additive public check on the deposited author DESeq2 table for T-47D Trop-2 KO xenografts (4 vs 3, NRG). Not lung. Not SKB264. Human symbol **CLDN4**. The single-gene CLDN4 padj is recorded below and is not the collateral result.

Reproduce:

```bash
python3 methods/gse289287_trop2ko_cldn4_recheck/analyze.py
python3 methods/gse289287_trop2ko_cldn4_recheck/sweep.py
```

Panel: `figures/fig_sweep_collateral.png`.

---

## TL;DR

The collateral that survives the sweep is **type I / type II IFN on the author ranking**. It does not depend on one ISG. It does not separate the seven samples. STING is not opened. Core NHEJ is not opened. CLDN4 stays non-significant on the author genome-wide test.

Primary test, locked before looking: prerank GSEA on the author Wald statistic, same engine as the earlier xenograft GSEA (weighted KS, 1,000 gene-set permutations, seed 42). Positive NES = up in Trop-2 KO. Family = 19 public IFN, APM, STING/cytosolic-DNA, and NHEJ sets (`gene_sets_sweep.json`). BH is inside that family.

| Set | NES | nominal p | family FDR |
|---|---:|---:|---:|
| Hallmark IFN-α | **+2.211** | <0.001 | **0.0063** |
| Reactome IFN-α/β | **+1.990** | <0.001 | **0.0063** |
| Hallmark IFN-γ | **+1.844** | <0.001 | **0.0063** |

Nominal p is the floor of 1,000 permutations, (0+1)/1001. No APM, STING, or NHEJ set passes this FDR. TACSTD2 log2FC **−3.263**, padj **1.95×10⁻⁶¹**.

---

## Sweep (methods, ranks, thresholds, leave-one-out)

Ranks: author Wald, author log2FC, and a Welch t on log2(x+1) after three normalizations (author size-factor counts, CPM, DESeq2 median-of-ratios recomputed from the raw counts in the same file). Sample scores use the mean of log2(x+1) and all 35 label permutations. Thresholds are author padj 0.05 / 0.10 / 0.20 and nominal p 0.05 / 0.10. Symbol aliases, applied only when the current symbol is absent and the previous symbol is present: STING1→TMEM173, CGAS→MB21D1, H2AX→H2AFX. TREX1 is absent from the table.

**IFN, author ranks.** Wald and log2FC agree. The same three sets are the only FDR < 0.05 hits on both ranks (log2FC NES +2.061, +2.014, +1.619). Across the 95 rank×set tests, those three Wald rows have sweep FDR **0.016**.

**IFN, Welch ranks.** The sign stays positive (Hallmark IFN-α NES about **+1.61** on all three normalizations) but the within-rank FDR rises to **0.11**. Re-ranking by a per-gene Welch t does not carry the set at FDR < 0.05.

**Leave-one-out ISG.** Dropping each of the 95 Hallmark IFN-α genes in the Wald rank leaves NES between **2.176 and 2.270**. All 95 drops stay at the permutation floor. Dropping each of 186 Hallmark IFN-γ genes leaves NES between **1.836 and 1.902**. The enrichment is not IFI44L, ISG15, or any other single gene. The sample-mean score is a different question: Hallmark IFN-α KO−WT is +0.29 on author-norm counts, exact two-sided p **0.40** (one-sided 0.17). CPM moves it to +0.38, two-sided p **0.26**. n = 4 vs 3 does not separate.

**Thresholds.** At author padj < 0.05 the IFN calls are still the five genes from the first pass (IFI44L, ISG15, IFI44, IFITM3, STAT2). Relaxing to padj < 0.10 adds IFIT3, OAS1, UBE2L6, MOV10 on the Hallmark IFN-α list (9 up, 1 down). That is a threshold choice, not a new test.

**STING.** Not up. Reactome STING NES **−1.113**, nominal p 0.144, family FDR 0.195. The IRF3 / type I IFN Reactome set is the strongest STING-family row and it is negative on every rank (Wald NES **−1.454**, p 0.029, family FDR **0.092**). STING1 (TMEM173) log2FC **+0.076**, p 0.816, baseMean 5.2, no padj. cGAS (MB21D1) baseMean **0.18**. The only author-significant gene in the Reactome STING set is **PRKDC**, which is down (log2FC **−0.490**, p 0.0015, padj **0.039**) and is also an NHEJ gene. Leave-one-out of the STING set never produces a positive NES.

**NHEJ.** Not a core-machinery result. KEGG NHEJ Wald NES **+1.295**, p 0.076, family FDR 0.154. The genes that clear padj < 0.05 point in opposite directions: POLM **+1.102**, padj **0.018**; PRKDC **−0.490**, padj **0.039**. XRCC4/5/6 and LIG4 do not. The classical GO set (5 genes) has Wald NES **−0.61**, p 0.47. The CPM mean score for KEGG NHEJ hits the permutation floor (two-sided p **0.029**, 1/35). BH across the 57 set×normalization sample tests gives FDR **0.27**. A floor p-value in one normalization is not a called NHEJ program.

**APM.** Custom 21-gene MHC-I set: Wald NES **+1.107**, p 0.156, family FDR 0.195. Zero genes at padj < 0.05. KEGG antigen processing (MHC I and II mixed) leans the other way (NES **−1.198**, p 0.078). HLA-F is the closest single gene (log2FC +0.893, padj 0.074).

Tables: `sweep_gsea.tsv`, `sweep_sample_scores.tsv`, `sweep_thresholds.tsv`, `sweep_loo.tsv`, `sweep_focus_genes.tsv`, `sweep_cldn4_norms.tsv`, `sweep_set_coverage.tsv`.

---

## CLDN4, recorded

Author DESeq2 is unchanged by the sweep: log2FC **+0.284**, p **0.129**, padj **0.466**. The groups overlap.

| Normalization | Welch p on log2(x+1) | Exact two-sided p |
|---|---:|---:|
| Author normalized counts | 0.041 | 0.086 |
| CPM | 0.0087 | 0.057 |
| Median-of-ratios | 0.038 | 0.057 |

CPM is the strongest unadjusted row. The exact two-sided permutation does not go below 0.057, and there is no new genome-wide padj. CLDN4 is not called.

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

## Figures

`figures/fig_sweep_collateral.png` is the sweep panel: Wald NES, the five-rank heatmap, IFN-α leave-one-out, and the sensor genes next to the genes that actually move.

`figures/fig_cldn4_ifn_apm.png` is the first-pass gene plot (TACSTD2, CLDN4, author IFN/APM calls).

Tables from the first pass: `cldn4_recheck.tsv`, `key_genes.tsv`, `ifn_genes.tsv`, `apm_genes.tsv`, `geneset_summary.tsv`, `sample_inventory.tsv`.

---

## 中文

扫过的是锁定的 19 个 IFN / APM / STING / NHEJ 基因集，不是全库。主检验是作者 Wald 秩的 prerank GSEA（1000 次置换，seed 42）。过家族 FDR 的只有三个 IFN 集：Hallmark IFN-α NES +2.211，Reactome IFN-α/β +1.990，Hallmark IFN-γ +1.844，FDR 都是 0.0063。丢掉任意一个 IFN-α 基因，NES 仍在 2.176–2.270。样本均值分不开（IFN-α 精确双侧 p = 0.40）。

STING 不是打开的：STING1 log2FC +0.076（p = 0.82），cGAS 几乎不表达。Reactome STING NES −1.113（p = 0.14）。集合里唯一 padj < 0.05 的是 PRKDC 下调（−0.490，padj 0.039）。NHEJ 核心（XRCC4/5/6、LIG4）不动；POLM 上调（padj 0.018）和 PRKDC 下调方向相反。KEGG NHEJ 的 NES +1.295，家族 FDR 0.15。CPM 样本分的置换 p 落到 1/35，但 57 次样本检验的 FDR 是 0.27。APM 自定义集 NES +1.107，p = 0.16，没有基因过 padj < 0.05。

CLDN4 作者 padj 仍是 0.466（p = 0.129）。CPM 上 Welch p = 0.0087，精确双侧 p = 0.057，不构成基因组范围的显著调用。
