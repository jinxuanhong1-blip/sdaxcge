# C4 analog — GSE334497 4T1 Trop2 KO: Cldn4, Cxcl9, IFN/APM

**Verdict: does not support Claim C4.** Trop2 knockout is real. *Cldn4* does not fall with confidence. The private CORE6 IFN/APM genes are flat. *Cxcl9* is up, together with a T-cell/cytotoxicity module, in **bulk immunocompetent tumors** — the source paper’s own infiltration story, not a tumor-intrinsic IFN program after CLDN4 loss.

Slice path: `results/w200/C4_GSE334497/`. Scripts: `scripts/w200/C4_GSE334497/`.

---

# English

## Claim being tested

Private **Claim C4**: CLDN4 knockdown opens an IFN / MHC-I / APM program (core: IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A).

**This dataset is an analog, not the claim.** GSE334497 is CRISPR **Trop2 (Tacstd2) KO vs WT** in 4T1 tumors grown 3 weeks in BALB/c mice (Wu *et al.*, *JITC* 2026; GEO public 2026-06-09). It is not a CLDN4 KD/KO. It is breast, not lung. RNA is from **frozen whole-tumor sections**, so IFN/chemokine signal can be stroma and infiltrating leukocytes.

The analog logic is: if Trop2 loss lowers *Cldn4*, the contrast can stand in for CLDN4 loss. That hinge fails here.

## Design (pre-specified)

| Item | Choice |
|---|---|
| Contrast | KO − WT, n = 5 vs 5 |
| Matrix | Author-normalized counts (`GSE334497_normalized_counts.csv.gz`, GRCm38 Ensembl) |
| Transform | log2(norm + 1); no re-estimation of size factors |
| Per-gene | Welch t + Mann–Whitney; Cohen’s *d*; BH-FDR over 15,073 expressed genes |
| Set tests | (1) competitive MWU of set log2FC vs background; (2) per-sample mean z-score, exact permutation of the 252 label splits |
| Focal genes | *Cldn4*, *Cxcl9* |
| Sets | USER_CORE6, ISG_CORE, MHC1_APM, CXCL_IFN, T_CYT; negative controls OXPHOS and MYC |

KO libraries: KO162, KO164, KO165, KO172, RESUB-KO163R.  
WT libraries: control170, RESUB-171R, RESUB-170R, RESUB-169R, RESUB-168R.

## Perturbation QC

| Gene | mean norm KO | mean norm WT | log2FC | Welch *p* | MWU *p* | *d* |
|---|---:|---:|---:|---:|---:|---:|
| **Tacstd2** | 11.0 | 159.2 | **−3.82** | **0.0011** | 0.0079 | −3.35 |

Knockout worked. Residual *Tacstd2* in KO tumors is expected (stroma + incomplete CRISPR). Genome-wide BH-FDR is 0.60 for *Tacstd2* because **zero genes** pass FDR < 0.05 at n = 5+5. Pre-specified tests, not discovery FDR, are the right frame.

## Cldn4 (the C4 hinge) — not a CLDN4-loss analog

| Gene | log2FC | Welch *p* | MWU *p* | *d* | Note |
|---|---:|---:|---:|---:|---|
| **Cldn4** | −0.82 | **0.25** | 0.42 | −0.79 | Highly expressed either way (KO 1.7k, WT 2.6k). Groups overlap. |
| Cldn7 | −0.82 | 0.091 | 0.016 | −1.34 | Paper’s stated mediator (protein localization). Rank test driven in part by one high WT. |
| Cldn1 | −2.18 | 0.035 | 0.032 | −1.61 | Strongest claudin RNA drop. |
| Ocln | −0.66 | 0.078 | 0.095 | −1.40 | Directional. |
| Epcam | −1.26 | 0.021 | 0.016 | −1.89 | Directional. |
| Tjp1, F11r | ~0 | >0.5 | — | ~0 | Flat. |

*Cldn4* moves in the hypothesized direction and is **not significant**. This series cannot be sold as “CLDN4 KD.” The Wu paper itself assigns the barrier phenotype to **claudin-7 protein**, not claudin-4, and does not report a *Cldn4* RNA claim for this matrix.

## Cxcl9 — up (nominal)

| Gene | mean norm KO | mean norm WT | log2FC | Welch *p* | MWU *p* | *d* |
|---|---:|---:|---:|---:|---:|---:|
| **Cxcl9** | 2331 | 1087 | **+1.08** | **0.0062** | 0.016 | +2.34 |
| Cxcl10 | 292 | 200 | +0.56 | 0.067 | 0.032 | +1.36 |
| Cxcl11 | — | — | — | — | — | absent from matrix |

*Cxcl9* is the cleanest pre-specified single-gene hit after *Tacstd2*. Genome-wide FDR = 0.66 (rank 136 by Welch *p*). CXCL_IFN set (n = 2 present): sample-score Δ = +1.33, exact two-sided perm *p* = **0.016**.

## IFN / APM — CORE6 does not open

Private CORE6 mouse orthologs (KO − WT):

| Gene | log2FC | Welch *p* | *d* |
|---|---:|---:|---:|
| Ifi27 | +0.22 | 0.42 | +0.55 |
| Ifi27l2a | +0.21 | 0.69 | +0.26 |
| Oas2 | −0.10 | 0.81 | −0.16 |
| Ifit1 | −0.02 | 0.94 | −0.05 |
| Mx1 | +0.16 | 0.67 | +0.29 |
| Isg15 | −0.09 | 0.84 | −0.13 |
| H2-K1 | +0.13 | 0.46 | +0.49 |
| H2-D1 | +0.06 | 0.78 | +0.19 |
| Ifi27l2b | absent | | |

These are null. *Isg15* / *Ifit1* / *Mx1* / *Oas2* do not move.

Set-level (exact permutation is the honest n = 5+5 test):

| Set | n tested | median log2FC | competitive MWU *p* | sample Δ | Welch *p* | **perm *p* (two)** |
|---|---:|---:|---:|---:|---:|---:|
| USER_CORE6 | 8 | +0.10 | 0.19 | +0.18 | 0.65 | **0.63** |
| ISG_CORE | 38 | +0.08 | 0.028 | +0.17 | 0.68 | **0.71** |
| MHC1_APM | 21 | +0.14 | 0.0012 | +0.38 | 0.38 | **0.36** |
| CXCL_IFN | 2 | +0.82 | 0.0021 | +1.33 | 0.018 | **0.016** |
| T_CYT | 11 | +0.76 | 2.9e-5 | +0.77 | 0.044 | **0.040** |
| CTRL_OXPHOS | 16 | −0.18 | 0.40 | −0.23 | 0.67 | 0.67 |
| CTRL_MYC | 19 | −0.15 | 0.25 | −0.20 | 0.57 | 0.58 |

Competitive *p* for ISG/APM is **not** a win. Background median log2FC is −0.03; a +0.08 to +0.14 set median beats background without any sample-level IFN score shift. Negative-control sets are null, as they should be.

T_CYT (*Cd3d/e, Cd8a, Gzmb, Prf1, Ifng, Nkg7, Klrd1*, …) **does** rise at the sample level (perm *p* = 0.040). *Prf1* +1.03 (*p* = 0.043); *Klrd1* +1.05 (*p* = 0.0097); *Cd8a* +0.67 (*p* = 0.11). That is the parsimonious explanation of *Cxcl9*: more lymphoid/cytotoxic RNA in the bulk slice, not an epithelial APM program.

## What this is and is not

**Is**

- A real Trop2-loss contrast (5 vs 5) with processed public counts.
- Evidence that *Cxcl9* (and a T-cell/cytotoxicity module) is higher in Trop2-KO 4T1 tumors.
- Consistent with the depositing paper’s GSEA (tight-junction programs in WT; inflammatory / T-cell cytotoxicity in KO).

**Is not**

- A CLDN4 knockdown.
- A lung model.
- Tumor-cell-intrinsic IFN/APM (no sorted epithelium, no in-vitro 4T1 arm on GEO).
- An independent test of Wu *et al.* — this **is** their matrix.
- Support for private C4 (CORE6 / MHC-I opening after CLDN4 loss).

## Bottom line

| Question | Answer |
|---|---|
| Did Trop2 KO work? | Yes. *Tacstd2* log2FC −3.82, *p* = 0.001. |
| Does *Cldn4* fall? | No confident yes. log2FC −0.82, *p* = 0.25. |
| Does *Cxcl9* rise? | Yes, nominally. log2FC +1.08, *p* = 0.006; set perm *p* = 0.016. |
| Does IFN/APM CORE6 open? | **No.** |
| C4 analog? | **Fail.** Wrong perturbation, hinge gene not significant, CORE6 null. *Cxcl9* up is infiltration-compatible bulk RNA, same paper. |

## Reproduce

```bash
pip install pandas numpy scipy statsmodels matplotlib
python3 scripts/w200/C4_GSE334497/analyze.py
```

Uses the shipped counts + `raw/ensembl_to_symbol.tsv` if present; otherwise re-downloads GEO and NCBI gene maps.

---

# 中文

## 要检验的主张

私有 **Claim C4**：CLDN4 敲低打开 IFN / MHC-I / APM（核心：IFI27、OAS2、IFIT1、MX1、ISG15、HLA-A）。

**本数据集是类比，不是该主张本身。** GSE334497 是 4T1 **Trop2 (Tacstd2) CRISPR KO 对 WT**，BALB/c 体内生长 3 周后的冰冻全瘤 RNA-seq（Wu 等，*JITC* 2026）。不是 CLDN4 KD/KO，不是肺，也不是纯化上皮。IFN/趋化因子可以来自间质和浸润免疫细胞。

类比成立的前提是 Trop2 丢失会拉低 *Cldn4*。这个铰链在这里断了。

## 设计（预先指定）

KO − WT，n = 5 vs 5。作者归一化 counts，log2(norm+1)。单基因 Welch + Mann–Whitney，全基因组 BH-FDR。基因集：竞争性 MWU（set log2FC vs 背景）+ 样本均值 z 分数的精确置换（252 种标签）。焦点基因：*Cldn4*、*Cxcl9*。

## 结果（诚实）

1. **敲除成立。** *Tacstd2* log2FC = −3.82，*p* = 0.001。n = 5+5 下全基因组 FDR < 0.05 的基因数为 **0**，所以只报告预先指定检验。
2. ***Cldn4* 没有可靠下降。** log2FC = −0.82，Welch *p* = 0.25。两组重叠，表达量仍然很高。不能把这个系列写成 CLDN4 丢失。论文自己把屏障功能归到 **claudin-7 蛋白**，不是 claudin-4。RNA 上掉得更明显的是 *Cldn1*（log2FC −2.18，*p* = 0.035）。
3. ***Cxcl9* 上升（名义显著）。** log2FC = +1.08，*p* = 0.006，*d* = +2.34。CXCL_IFN 精确置换 *p* = 0.016。全基因组 FDR = 0.66。
4. **IFN/APM CORE6 没有打开。** *Isg15 / Ifit1 / Mx1 / Oas2 / H2-K1* 基本不动。USER_CORE6、ISG_CORE、MHC1_APM 的样本水平置换 *p* 分别为 0.63 / 0.71 / 0.36。竞争性检验的“显著”只是相对略负的背景多了一点正 log2FC，不是 IFN 程序被打开。
5. **T 细胞/细胞毒模块上升**（T_CYT 置换 *p* = 0.040）。这是 *Cxcl9* 更省事的解释：全瘤切片里淋巴/细胞毒 RNA 更多，不是上皮 APM。

## 一句话

GSE334497 验证了 Trop2 丢失后全瘤 *Cxcl9* 和 T 细胞毒签名升高（与原文一致），**不能**当作 CLDN4 KD 后肿瘤内在 IFN/APM 开放的公开类比。Claim C4：**不支持**。
