# CLDN4 / TACSTD2 (TROP2) true KD–KO public transcriptomes — IFN, MHC-I, tight junction, immune exclusion

English first, 中文在后. Every number below is taken from `results/kdko/*.tsv` produced by `scripts/kdko/04`–`10`. No accession was invented: each row in `notes/kdko_catalog.tsv` was retrieved from NCBI GEO E-utilities or EBI BioStudies and then opened on the public record.

Bessede 2024 (high TACSTD2 → immune exclusion via a claudin / tight-junction barrier) is treated **only as a prior**. It is not assumed. The one public series whose title matches that prior (`GSE334497`) is re-analysed from the deposited matrix.

---

## 1. What was searched, what exists

Two GEO rounds (`scripts/kdko/01_search_geo.py`, `02_search_geo_round2.py`) plus an ArrayExpress/BioStudies search.

| Query class | Unique GEO entries | Gene mentioned in title/summary | True KD/KO of CLDN4 or TACSTD2 with an open processed matrix |
|---|---:|---:|---|
| Round 1 (gene + knockdown/KO/sh/si/CRISPR) | 570 | 39 | 6 series (below) |
| Round 2 (title/description + lung-focused) | 964 | (filtered) | same 6; plus GSE50927 recovered as the only lung hit |
| ArrayExpress | 202 unique accessions | 13 | 3, all GEO mirrors (E-GEOD-334497 / 50927 / 22493) |

**Analysed true perturbations (human and mouse):**

| Accession | Gene | Perturbation | System | n | Immune host |
|---|---|---|---|---|---|
| [GSE334497](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334497) | TACSTD2 | KO | mouse 4T1 tumour | 5 vs 5 | immunocompetent |
| [GSE289287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289287) | TACSTD2 | KO | human T-47D xenograft | 4 vs 3 | NRG (no T/B/NK) |
| [GSE15212](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15212) | TACSTD2 | 2× siRNA, 72 h | human SW480 | 6 vs 9 | in vitro |
| [GSE207704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207704) | CLDN4 | CRISPR KO | human T47D + MCF7 | group means of 2+2 vs 2+2 | in vitro |
| [GSE50927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50927) | CLDN4 | germline KO | mouse whole lung | **1 vs 1** | immunocompetent |
| [GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493) | CLDN4 | siRNA vs CLDN4-OE | human SKOV-3 | 3 two-colour arrays | in vitro |

**Lung.** The only verified lung CLDN4 perturbation is `GSE50927` (mouse, unreplicated). No human lung CLDN4 or TACSTD2 KD/KO transcriptome, and no mouse lung Tacstd2 KO transcriptome, was found (see `results/kdko/missing_controlled.tsv`).

LINCS L1000 (`GSE106127` CGS metadata, `GSE92742` pert_info) does **not** contain CLDN4 or TACSTD2 perturbagens.

---

## 2. Knockdown / knockout validity

| Series | Target | Residual vs control | Test | n |
|---|---|---|---|---|
| GSE334497 | Tacstd2 | 7.1% of WT (log2FC −3.82) | Welch t = −5.30, p = 1.1×10⁻³ | 5 vs 5 |
| GSE289287 | TACSTD2 | 9.1% of WT (DESeq2 log2FC −3.26) | submitter padj = 1.9×10⁻⁶¹ | 4 vs 3 |
| GSE15212 | TACSTD2 | 11% of neg-ctrl (log2FC −3.16) | Welch p = 2.2×10⁻¹² | 6 vs 9 |
| GSE15212 siRNA#1 | TACSTD2 | 10% (log2FC −3.26) | p = 1.0×10⁻⁶ | 3 vs 9 |
| GSE15212 siRNA#4 | TACSTD2 | 12% (log2FC −3.07) | p = 3.7×10⁻⁵ | 3 vs 9 |
| GSE207704 T47D | CLDN4 | 48% of WT FPKM (42.5 → 20.4; log2FC −1.04) | no replicate p (collapsed) | 2 vs 2, collapsed |
| GSE207704 MCF7 | CLDN4 | 59% of WT FPKM (85.9 → 50.8; log2FC −0.75) | no replicate p (collapsed) | 2 vs 2, collapsed |
| GSE50927 naive lung | Cldn4 | 1.5% of WT (logFC −6.06) | submitter edgeR FDR = 4.1×10⁻²⁶ | 1 vs 1 |
| GSE22493 | CLDN4 | mean log2(KD/OE) = −1.23 (arrays: NA, −1.74, −0.71) | one-sample t not trusted (arrays disagree) | 3 arrays |

All six series therefore contain a real (or, for GSE207704, partial) loss of the targeted transcript. GSE207704 is a hypomorph at the RNA level, not a clean null.

---

## 3. Does losing one gene change the other?

`results/kdko/cross_reciprocal_CLDN4_TACSTD2.tsv`

| Perturbation | Other gene | log2FC | p | Call |
|---|---|---|---|---|
| TACSTD2 KO, 4T1 tumours (GSE334497) | Cldn4 | −0.82 | 0.25 | no significant change; Cldn1 −2.18 p=0.035, Epcam −1.26 p=0.022 |
| TACSTD2 KO, T-47D xenograft (GSE289287) | CLDN4 | **+0.28** | 0.13 (q=0.47) | not down; CLDN7 +0.51 p=0.0028 q=0.057 |
| TACSTD2 siRNA, SW480 (GSE15212) | CLDN4 | −0.008 | 0.97 | no change; CLDN1 **+0.70** q=0.004 |
| CLDN4 KO, T47D+MCF7 (GSE207704) | TACSTD2 | −0.82 (T47D −0.62, MCF7 −1.01) | no replicate p | same direction in both lines; **descriptive only** |
| CLDN4 KO, naive lung (GSE50927) | Tacstd2 | −0.14 | 0.49 (FDR=1) | no change |
| CLDN4 siRNA, SKOV-3 (GSE22493) | TACSTD2 | — | — | **not on the platform** |

**Conclusion.** Reciprocal transcriptional control of CLDN4 ↔ TACSTD2 is **not** a general result. The only hint is the GSE207704 group-mean FPKM drop of TACSTD2 after CLDN4 KO, which cannot be given a p-value from public data. TACSTD2 loss does not consistently lower CLDN4 mRNA.

---

## 4. Tight-junction / claudin / epithelial-identity programmes

| Series | Claudin-family set | Tight-junction core | Epithelial identity |
|---|---|---|---|
| GSE334497 Trop2 KO (n=5 vs 5) | Δz = −0.57, sample p=0.037, q=0.071, d=−1.58 | Δz = −0.65, p=0.035, q=0.071, d=−1.63 | Δz = −1.17, p=0.0069, q=0.050, d=−2.29 |
| GSE289287 TACSTD2 KO xg (n=4 vs 3) | Δz = +0.48, p=0.18 | Δz = +0.28, p=0.31 | Δz = +0.22, p=0.43 |
| GSE15212 TACSTD2 siRNA (n=6 vs 9) | Δz = +0.16, p=0.22 | Δz = −0.05, p=0.46 | Δz = −0.13, p=0.29 |
| GSE207704 CLDN4 KO (gene-level) | median log2FC −0.55, competitive p=0.0016, q=0.018 | median −0.14, p=0.0039, q=0.021 | median −0.09, p=0.089 |
| GSE50927 Cldn4 KO lung (n=1) | median logFC 0.00, p=0.35 | median −0.006, p=0.50 | median −0.11, p=0.037, q=0.066 |
| GSE22493 CLDN4 siRNA | median −0.22, p=0.34 | median −0.06, p=0.80 | median −0.02, p=0.70 |

The **claudin / tight-junction down-shift after TACSTD2 loss is seen only in immunocompetent 4T1 tumours** (GSE334497). It is not recovered as a cell-intrinsic effect in human T-47D xenografts or SW480 siRNA. CLDN4 KO in human breast lines (GSE207704) does lower the claudin-family gene set at the competitive-gene level, as expected when CLDN4 itself is a member, but this is not a TACSTD2-driven barrier programme.

---

## 5. IFN and MHC-I — the focus of this update

Figure: `results/kdko/fig_cross_ifn_mhci.png`. Numbers: `cross_ifn_mhci_set_summary.tsv`, `cross_ifn_mhci_gene_log2FC.tsv`.

### 5.1 TACSTD2 / TROP2 loss

**Type-I IFN is the most consistent cell-intrinsic signal.**

- **GSE15212 SW480 siRNA (best-powered in-vitro design, two independent oligos).** IFN-I set Δz = +0.47, Welch p = 6.6×10⁻⁵, BH q = 0.0011, Cohen’s d = 2.90, n=6 vs 9. Replicated per oligo: siRNA#1 p=0.012; siRNA#4 p=2.0×10⁻⁵, q=3.2×10⁻⁴. Leading genes: STAT1 log2FC +0.43, q=0.017; OAS1 +0.35, q=0.022; IRF7 +0.76, p=0.032 (q=0.22). MHC-I set is unchanged (Δz = +0.012, p=0.97). B2M −0.18 p=0.27; HLA-A −0.11 p=0.12; TAP1 +0.14 p=0.050.
- **GSE289287 T-47D TACSTD2-KO xenografts (NRG).** IFN-I genes move together (competitive Mann-Whitney on DESeq2 stat p = 6.6×10⁻¹⁸; median log2FC of the set +0.69) but the **sample-level** set test is under-powered (Δz = +0.87, p=0.21, n=4 vs 3). Per-gene DESeq2: ISG15 +1.44 q=0.0012; STAT2 +0.71 q=0.031; IFIT3 +0.96 q=0.080; OAS1 +0.86 q=0.082; MX1 +0.93 q=0.11; STAT1 +0.74 p=0.014 q=0.15. MHC-I set sample p=0.70. B2M +0.54 p=0.049 q=0.30; HLA-A +0.63 p=0.071 q=0.36; HLA-B +0.56 p=0.054.
- **GSE334497 4T1 Trop2-KO tumours (immunocompetent).** Type-I IFN is **not** up (Δz = +0.17, p=0.69; Isg15 −0.09 p=0.84; Ifit1 −0.02 p=0.94). Type-II / IFN-γ-adjacent programmes **are** up: IFN-II set Δz = +0.71, p=0.040, q=0.071, d=1.57; competitive p=4.6×10⁻⁹. Cxcl9 +1.08 p=0.0062; Cxcl10 +0.56 p=0.067; Ifng +0.89 p=0.13. MHC-I set sample p=0.48, but B2m +0.41 p=0.040. This is an **immune-infiltrate / IFN-γ** signature, not a tumour-intrinsic type-I ISG burst.

### 5.2 CLDN4 loss

- **GSE50927 naive mouse lung (n=1 vs 1 — descriptive ranking only).** IFN-I, IFN-II and MHC-I gene sets are all shifted up versus the genome-wide background (competitive Mann-Whitney p = 8.2×10⁻⁶, 1.2×10⁻⁷, 6.9×10⁻⁴; q = 4.4×10⁻⁵, 9.4×10⁻⁷, 0.0018). Individual submitter edgeR FDRs that survive 0.05: Isg15 +1.08 FDR=0.0031; Ifit1 +0.58 FDR=0.028; Oasl2 +0.54 FDR=0.018; B2m +0.66 FDR=0.034; Psmb9 +0.90 FDR=0.037; Ccl5 +2.71 FDR=0.0011; Prf1 +1.35 FDR=0.0016; Gzmb +2.08 FDR=0.046. **These p-values are not replicate-based.** The same direction is **not** stable across the two Cldn4-KO VILI mice (VILI-low actually *lowers* the T-cell-cytotoxicity set, competitive p=6×10⁻⁵).
- **GSE207704 T47D/MCF7 CLDN4-/- (group-mean FPKM).** IFN-I set moves **down** (median log2FC −0.21, competitive p=0.0071, q=0.026, 22 genes). MHC-I is unmeasured for most members (B2M, HLA-A, TAP1, PSMB9 all absent after the FPKM≥1 filter). ISG15 is −1.91 in T47D and +0.11 in MCF7 — not reproducible across lines. Cross-line transcriptome concordance is null (Pearson r = −0.019, n=12591).
- **GSE22493 SKOV-3.** No IFN or MHC-I set shift (all competitive q > 0.8). Inter-array Pearson r of log-ratios is −0.18 / +0.20 / −0.16, so this series is not used for any positive claim.

### 5.3 IFN / MHC-I take-home

1. **TACSTD2 loss induces a type-I IFN programme in human cancer cells** (SW480 siRNA, sample-level q=0.001; T-47D KO xenografts, gene-level p=6.6×10⁻¹⁸, sample-level n.s.). MHC-I / APM is at most a weak co-movement (B2M p=0.049 in GSE289287, q=0.30; null in GSE15212).
2. **In an immunocompetent tumour, TACSTD2 loss induces IFN-γ / T-cell / chemokine programmes, not type-I ISGs** (GSE334497). That is the immune-inclusion side of the Bessede prior, and it is supported by this matrix.
3. **CLDN4 loss does not reproduce the human cell-intrinsic type-I IFN induction.** The only CLDN4 series with an IFN-up ranking is unreplicated whole lung (GSE50927). The only replicated human CLDN4 KO (GSE207704) trends IFN-I *down* and has no usable MHC-I read-out.

---

## 6. Immune exclusion, ICI genes, and the Bessede 2024 prior

Prior (not assumed): high TACSTD2 builds a claudin/TJ barrier that excludes T cells and blunts checkpoint blockade.

**Supported, in the one immunocompetent tumour series (GSE334497, n=5 vs 5):**

| Set / gene | Effect (KO − WT) | p | q | n |
|---|---|---|---|---|
| ICI checkpoint gene set | Δz +0.69, d=2.50 | 0.0047 | 0.050 | 5 vs 5 |
| Ayers 6-gene IFN-γ GEP | Δz +1.04, d=2.27 | 0.0093 | 0.050 | 5 vs 5 |
| T-cell infiltration (broad) | Δz +0.61, d=1.60 | 0.037 | 0.071 | 5 vs 5 |
| T-cell cytotoxicity | Δz +0.69, d=1.19 | 0.097 | 0.15 | 5 vs 5 |
| Ayers 18-gene TIS | Δz +0.64, d=1.22 | 0.11 | 0.15 | 5 vs 5 |
| Ctla4 | log2FC +0.62 | 0.0043 | 0.71 | 5 vs 5 |
| Cxcl9 | log2FC +1.08 | 0.0062 | 0.71 | 5 vs 5 |
| Klrd1 | log2FC +1.05 | 0.0097 | 0.74 | 5 vs 5 |
| Prf1 | log2FC +1.03 | 0.043 | 0.85 | 5 vs 5 |
| Cd8a | log2FC +0.67 | 0.11 | 0.89 | 5 vs 5 |
| Cd274 (PD-L1) | log2FC +0.77 | 0.10 | 0.89 | 5 vs 5 |

Direction matches the prior (Trop2 loss → more inflamed / less excluded). Magnitude at the single-gene level does not survive genome-wide BH FDR (smallest q among these is 0.71); the claim that is statistically supported is the **set-level** ICI and IFN-γ-GEP shift.

**Not supported as a general cell-intrinsic rule:**

- Human TACSTD2 KD/KO does **not** down-regulate CLDN4 or the claudin/TJ set (GSE15212, GSE289287).
- Human TACSTD2 KD/KO does **not** induce CXCL9/10 or a TIS (GSE15212 CXCL9 −0.02 p=0.69; GSE289287 CXCL10 −0.14 p=0.79; CD274 not measured in the xenograft protein-coding table).
- GSE289287 is on an NRG host, so it cannot test exclusion.
- No public series pairs a CLDN4 or TACSTD2 perturbation with ICI treatment and a response label.

**Comparator.** DSG2 KO in T-47D cells (GSE289287, not a TACSTD2 perturbation) also raises IFN-I genes (competitive p < 10⁻⁵). Genome-wide concordance with TACSTD2-KO xenografts is r = −0.085 (p=1.0×10⁻²¹, n=12644) — the two adhesion knockouts are not the same programme.

---

## 7. Design caveats that change the interpretation

- **GSE334497** has a RESUB vs original library split that is unbalanced (4/5 WT are RESUB; 4/5 KO are not). PCA on the top 2,000 variable genes: PC1 (38% variance) is not associated with genotype (p=0.48) or prep group (p=0.28); PC2 (20%) *is* associated with genotype (p=0.025) and not with prep (p=0.63). The genotype signal is therefore not an obvious batch artefact, but n=5 vs 5 still yields no genome-wide q<0.1.
- **GSE207704** residual CLDN4 RNA is ~50% of WT. Claims are fold-changes of group means.
- **GSE50927** is n=1. Used only as a ranked list.
- **GSE22493** arrays do not agree. Used only to document that TACSTD2 cannot be tested and that no trustworthy IFN claim can be made.
- **GSE15212** two-siRNA genome-wide r=0.218. On-target TACSTD2 KD is clear; off-target load is real. The IFN-I set rising with *both* oligos is the reason that result is kept.

---

## 8. Files

| Path | Content |
|---|---|
| `notes/kdko_catalog.tsv` | Verified accessions, including excluded / missing / controlled |
| `results/kdko/missing_controlled.tsv` | Why each gap or comparator matters |
| `results/kdko/fig_cross_ifn_mhci.png` | Cross-series IFN / MHC-I figure |
| `results/kdko/fig_GSE334497_trop2ko_4T1.png` | Immunocompetent Trop2 KO |
| `results/kdko/fig_GSE289287_tacstd2ko_T47D.png` | Human TACSTD2 KO xenograft |
| `results/kdko/fig_GSE15212_tacstd2_kd_SW480.png` | Human TACSTD2 siRNA |
| `results/kdko/fig_GSE207704_cldn4ko_T47D_MCF7.png` | Human CLDN4 CRISPR KO |
| `results/kdko/fig_GSE50927_cldn4ko_lung.png` | Mouse lung Cldn4 KO (n=1) |
| `results/kdko/fig_GSE22493_cldn4_kd_SKOV3.png` | SKOV-3 CLDN4 siRNA (untrusted) |
| `scripts/kdko/` | Search, download, analysis, synthesis |

Raw matrices live under `data/kdko/` (gitignored). Checksums: `results/kdko/downloads.tsv`.

---

## 中文摘要

**检索。** GEO 两轮 + ArrayExpress。真正的 CLDN4 或 TACSTD2/TROP2 敲低/敲除且有开放处理后矩阵的系列只有 6 个：人 TACSTD2 siRNA（GSE15212，SW480，n=6 vs 9）、人 TACSTD2 KO 异种移植（GSE289287，T-47D，n=4 vs 3，NRG 免疫缺陷）、小鼠 TACSTD2 KO 肿瘤（GSE334497，4T1，n=5 vs 5，免疫健全）、人 CLDN4 CRISPR KO（GSE207704，T47D+MCF7，公开文件只有组平均 FPKM）、小鼠肺 Cldn4 KO（GSE50927，**每组 n=1**）、人 CLDN4 siRNA（GSE22493，对照是 CLDN4 过表达，三张芯片互相关为负，不可信）。**没有**人肺 CLDN4 或 TACSTD2 的 KD/KO 转录组，也没有小鼠肺 Tacstd2 KO。LINCS 里没有这两个基因的 shRNA。

**敲除是否成立。** TACSTD2 三个系列残留约 7–12%（GSE334497 p=1.1×10⁻³；GSE289287 padj=1.9×10⁻⁶¹；GSE15212 p=2.2×10⁻¹²）。GSE207704 的 CLDN4 残留 48%/59%，是弱等位而不是干净 null。GSE50927 的 Cldn4 logFC=−6.06。

**丢掉一个会不会改变另一个。** 不会作为普遍规律。TACSTD2 丢失后 CLDN4 在 GSE15212 不变（log2FC −0.008，p=0.97），在 GSE289287 略升（+0.28，p=0.13），在 GSE334497 不显著下降（−0.82，p=0.25）。CLDN4 丢失后 TACSTD2 在 GSE50927 不变（−0.14，p=0.49）；GSE207704 两组均值都下降（T47D −0.62，MCF7 −1.01）但没有可计算的重复 p；GSE22493 平台上没有 TACSTD2。

**紧密连接。** 只有免疫健全的 4T1 Trop2 KO（GSE334497）出现 claudin / TJ / 上皮身份下调（claudin 家族 Δz=−0.57，p=0.037；TJ Δz=−0.65，p=0.035；上皮身份 Δz=−1.17，p=0.0069，q=0.050）。人细胞/异种移植里 TACSTD2 丢失不下调 CLDN4 或 TJ 基因集。

**IFN / MHC-I（本次重点）。**

- 人 TACSTD2 丢失 → **I 型 IFN 上调**：GSE15212 IFN-I 基因集 Δz=+0.47，p=6.6×10⁻⁵，q=0.001，d=2.90；两条独立 siRNA 都在（p=0.012 与 2.0×10⁻⁵）。GSE289287 基因水平竞争检验 p=6.6×10⁻¹⁸（ISG15 +1.44，q=0.0012；STAT2 +0.71，q=0.031），样本水平因 n=4 vs 3 不显著（p=0.21）。MHC-I 基因集在这两套里都不显著（GSE15212 p=0.97；GSE289287 p=0.70）；GSE289287 的 B2M +0.54，p=0.049，q=0.30。
- 免疫健全肿瘤里 TACSTD2 丢失 → **II 型 / T 细胞程序上调，I 型 IFN 不上**：GSE334497 IFN-I p=0.69；IFN-II Δz=+0.71，p=0.040，q=0.071；Cxcl9 +1.08，p=0.0062。
- CLDN4 丢失 **不重复** 人细胞的 I 型 IFN 诱导。GSE207704 IFN-I 反而下降（中位 log2FC −0.21，竞争 p=0.007，q=0.026）。GSE50927 肺里 IFN-I/II 和 MHC-I 排序上调（竞争 q 分别为 4.4×10⁻⁵、9.4×10⁻⁷、0.0018），但是 n=1，不能当重复证据。

**Bessede 2024 先验（高 TACSTD2 → 免疫排斥）。** 只在 GSE334497 得到支持：ICI 检查点基因集 Δz=+0.69，p=0.0047，q=0.050；Ayers 6 基因 IFN-γ GEP Δz=+1.04，p=0.0093，q=0.050。单基因（Ctla4、Cxcl9、Cd8a、Cd274）方向一致但全基因组 BH q 都不<0.1。人细胞/免疫缺陷移植中，TACSTD2 丢失不下调 claudin、不诱导 CXCL9/10，因此“屏障 → 排斥”不是细胞自主的普遍规则。没有带 ICI 疗效标签的 CLDN4/TACSTD2 扰动公开集。

**受控对照。** GSE289287 的 DSG2 KO 也会抬高 IFN-I（竞争 p<10⁻⁵），与 TACSTD2 KO 全基因组相关 r=−0.085（n=12644），说明“粘附丢失 → IFN-I”在 T-47D 里不是 TACSTD2 特有，但两套程序并不相同。
