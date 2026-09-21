# NHEJ-STING wave — Bitler MCT CLDN4-knockdown RPPA

No PRIDE or MassIVE accession. The public numbers are the supplement tables, not an RPPA matrix.

Paper: Yamamoto, Webb, Bitler et al., *Mol Cancer Ther* 2022;21:647–657. DOI [10.1158/1535-7163.MCT-21-0827](https://doi.org/10.1158/1535-7163.MCT-21-0827). PMC8988515.

---

## 中文

### 检索

OVCAR3 shCtrl（n=3）对 shCLDN4#1（n=3）的反相蛋白阵列在 MD Anderson RPPA Core（Yiling Lu）完成。正文数据声明是 “Data are available upon request from the corresponding author.” RPPA 不是质谱，本来就不会以 PXD/MSV 的形式入库，除非作者另交了一份 MS。

| 查询 | 结果 |
| --- | --- |
| PRIDE v2 keyword：`CLDN4`、`claudin-4`、`claudin 4`、`shCLDN4`、`Bitler` | 0 |
| PRIDE v2 keyword：`claudin` | 30 个项目，标题/摘要里没有 CLDN4 knockdown，也没有 RPPA |
| PRIDE v2 keyword：`OVCAR3` | 22，均为质谱项目，不是这篇 RPPA |
| EBI Search pride 短语 `"CLDN4"` / `"claudin-4"` / `"claudin 4"` / `shCLDN4` | hitCount = 0 |
| EBI Search `Bitler` | 9，是 bitter/bitterness 词干误伤（奶酪苦肽、番茄苦味、denatonium），不是 Benjamin Bitler |
| OmicsDI `(source:pride OR source:massive) AND (CLDN4 OR "claudin-4" OR shCLDN4)` | 0 |
| OmicsDI `(source:massive) AND (CLDN4 OR claudin OR shCLDN4)` | 0 |
| OmicsDI `RPPA AND (CLDN4 OR claudin-4 OR OVCAR3)` | 2 篇文献记录，不是蛋白组 accession |

**Accession：无。** 近邻但不是本实验：PXD031094（claudin 家族 CoIP，不是敲低）、PXD066158（CLDN6 缺失）、PXD003651（CLDN3 外泌体）、PXD005292（Claudin-low 乳腺癌分型）。

### 补充表

Figshare collection [10.1158/1535-7163.c.6543310](https://doi.org/10.1158/1535-7163.c.6543310)。四张表都下了。没有 RPPA 矩阵，没有抗体清单，没有 DNA-PKcs / 53BP1 / XRCC1 的 fold。

- Table S1：TCGA 卵巢肿瘤转录组，CLDN4 高 vs 低。19,834 基因，q<0.05 正好 1,582，与正文一致。这不是敲低蛋白。
- Table S2：细胞系 BRCA 状态。
- Table S3：OVCAR3 shCtrl vs shCLDN4 药物屏 + DepMap 相关。化合物名单里没有 NU7441、peposertib、nedisertib、M3814、AZD7648。
- Table S4：15 例原代瘤临床注释。

### 三个蛋白怎么打分

`kd_logic_score` 只认敲低蛋白证据：正文写明显著下降 = −1；RPPA/补充表都没点名 = 0。不编 fold。Figure 3 的星号没有逐条写进图注，所以不把 \*/\*\*/\*\*\* 安到某一根柱子上。

| 蛋白 | 基因 | 敲低蛋白 | kd_logic_score | Table S1 log2(CLDN4高/低) | q | TCGA FDR<0.05 |
| --- | --- | --- | --- | --- | --- | --- |
| DNA-PKcs | PRKDC | 未报告 | **0** | **−0.27**（低表达组更高） | 0.00445 | 是 |
| 53BP1 | TP53BP1 | 下降（Fig 3A–B；两个 shRNA） | **−1** | −0.15 | 0.13 | 否 |
| XRCC1 | XRCC1 | 下降（Fig 3A、3C；两个 shRNA） | **−1** | +0.04 | 0.76 | 否 |

读法：

- NHEJ 因子里，这篇 RPPA **只撑得住 53BP1 蛋白下降**。XRCC1 也下降，但它是单链断裂/BER 支架，不是 NHEJ 核心酶。DNA-PKcs 没有蛋白读数。
- 同一篇的功能实验（不是这张表）：HDR 不变，NHEJ 下降，olaparib 诱导的 53BP1 foci 被削弱（Fig 4）。
- Table S1 的方向不要和 RPPA 焊在一起。三个基因里只有 **PRKDC** 进了 FDR<0.05，而且是 CLDN4 低的肿瘤转录更高（log2 −0.27）。53BP1、XRCC1 的转录不显著。这张 RNA 表不支持 “CLDN4 高 = DNA-PKcs/53BP1/XRCC1 转录更高”。

### STING 侧（同一张 Table S1，不是 RPPA）

RPPA 正文没有点名 STING、cGAS、TBK1。Table S1：STING1 log2 +0.37，p=0.0045，q=0.053，卡在 1,582 个 FDR<0.05 基因门外，方向是 CLDN4 高组略高。CGAS、TBK1、IRF3、CCL5、CXCL10 的 q 都 >0.1。NHEJ 结构基因 XRCC4/5/6、LIG4、NHEJ1 也不显著。这张公共表没有 “CLDN4 低 → STING/ISG 升高” 的转录支持。

---

## English

### Search

The OVCAR3 shCtrl (n=3) vs shCLDN4#1 (n=3) array was run at the MD Anderson RPPA Core. The paper’s data statement is request-from-the-corresponding-author. RPPA is not a mass-spectrometry deposit.

PRIDE keyword and EBI phrase searches for `CLDN4`, `claudin-4`, `claudin 4`, and `shCLDN4` returned no projects. OmicsDI restricted to `source:pride` or `source:massive` returned 0 for those terms. The PRIDE keyword `claudin` returned 30 projects; none is a CLDN4 knockdown and none is an RPPA. The token `Bitler` returned 9 PRIDE records that are bitter/bitterness stemming (cheese bitter peptides, tomato bitterness, denatonium), not this author.

**Accession: none.**

### What was digitized

Figshare supplements for MCT-21-0827. There is no RPPA matrix and no fold-change for DNA-PKcs, 53BP1, or XRCC1.

Table S1 is the TCGA ovarian transcriptome, CLDN4-high vs CLDN4-low (19,834 genes; 1,582 at q<0.05, matching the paper). It is not the knockdown proteome. Table S3’s compound list has no DNA-PKcs inhibitor.

### Scores

`kd_logic_score` is −1 only when the manuscript reports a significant protein decrease after CLDN4 knockdown, and 0 when the protein is not named. Folds were not invented. Stars on Figure 3 are not assigned to individual bars because the legend does not map them.

| Protein | Gene | KD protein | kd_logic_score | Table S1 log2(high/low) | q | In the FDR<0.05 set |
| --- | --- | --- | --- | --- | --- | --- |
| DNA-PKcs | PRKDC | not reported | **0** | **−0.27** (higher in CLDN4-low) | 0.00445 | yes |
| 53BP1 | TP53BP1 | down (Fig 3A–B, both shRNAs) | **−1** | −0.15 | 0.13 | no |
| XRCC1 | XRCC1 | down (Fig 3A, 3C, both shRNAs) | **−1** | +0.04 | 0.76 | no |

53BP1 is the NHEJ protein this RPPA actually supports as down after CLDN4 loss. XRCC1 is down as well and is an SSBR/BER scaffold. DNA-PKcs has no protein value. On the separate TCGA RNA table, PRKDC is the only one of the three inside FDR<0.05, and the transcript is higher when CLDN4 is low. That RNA direction is not the RPPA result.

Same-paper functional context, not re-measured here: homology-directed repair unchanged, NHEJ reduced, olaparib-induced 53BP1 foci reduced (Figure 4).

STING was not an RPPA call. In Table S1, STING1 is log2 +0.37, q=0.053 (outside the 1,582-gene set, slightly higher in CLDN4-high). CGAS, TBK1, IRF3, CCL5, CXCL10 and the other core NHEJ transcripts (XRCC4/5/6, LIG4, NHEJ1) are not FDR<0.05. This public table does not show a CLDN4-low STING/ISG increase.
