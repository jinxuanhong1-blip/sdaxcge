# TACSTD2 / CLDN4 in public lung-cancer scRNA-seq: epithelial restriction, not immune

**Slice:** `fable_tisch` · TISCH2 open NSCLC objects only · no LuCA (12 GB) · no fabricated numbers.

All statistics below were computed by
`scripts/fable_tisch/02_trop2_cldn4_tumor_vs_immune.py` from the four
TISCH2 matrices listed in [DATA_SOURCES.md](DATA_SOURCES.md). Expression
values are TISCH2 `log2(TPM/10 + 1)`. Cell-level Mann–Whitney *p* values
that underflow `float64` are written `p ≈ 0`. Per-patient Wilcoxon
*p* values are exact (they equal \(2^{-n}\) when every pair has the same
sign).

---

## English

### Question

Are **TACSTD2** (TROP2) and **CLDN4** expressed in lung-tumor epithelial /
malignant cells and largely absent from immune cells in public NSCLC
atlases? And is that contrast an infiltrate artifact (does tumor
TACSTD2/CLDN4 track the immune fraction)?

### Objects (all <2 GB or subsettable)

| Dataset | Malignant / proxy | Immune | Residual non-malignant epithelium |
|---|---:|---:|---|
| NSCLC_GSE131907 (Kim 2020; tumor-tissue subset) | 9,573 tumor-tissue Epithelial (proxy; TISCH2 has no Malignant call) | 53,114 | none left after the tumor-tissue restriction; compared instead to 3,265 **normal-lung** Epithelial cells |
| NSCLC_EMTAB6149 (Lambrechts 2018) | 5,728 Malignant | 30,694 | 1,346 Alveolar |
| NSCLC_GSE127465 (Zilionis 2019) | 3,995 Malignant | 26,055 | none annotated |
| NSCLC_GSE148071 (Wu 2021) | 48,118 Malignant | 20,077 | 8,147 Alveolar + Basal + Epithelial |

Controls: **EPCAM** (epithelial +) and **PTPRC** / CD45 (immune +).

### 1. Malignant vs pooled immune

| Dataset | Gene | mean malig | mean immune | % pos malig | % pos immune | log2FC | AUROC | cell-level *p* | paired Wilcoxon |
|---|---|---:|---:|---:|---:|---:|---:|---|---|
| GSE131907 | TACSTD2 | 1.534 | 0.051 | 81.4 | 4.2 | 4.90 | 0.895 | ≈ 0 | 9/9 patients malig > imm, *p* = 0.00391 |
| GSE131907 | CLDN4 | 1.417 | 0.040 | 80.3 | 3.4 | 5.16 | 0.891 | ≈ 0 | 9/9, *p* = 0.00391 |
| EMTAB6149 | TACSTD2 | 1.415 | 0.101 | 59.3 | 5.9 | 3.81 | 0.776 | ≈ 0 | no Patient column |
| EMTAB6149 | CLDN4 | 1.190 | 0.071 | 62.0 | 4.3 | 4.07 | 0.792 | ≈ 0 | no Patient column |
| GSE127465 | TACSTD2 | 1.208 | 0.016 | 58.7 | 0.9 | 6.22 | 0.790 | ≈ 0 | 7/7, *p* = 0.0156 |
| GSE127465 | CLDN4 | 0.896 | 0.007 | 50.3 | 0.5 | 6.91 | 0.749 | ≈ 0 | 7/7, *p* = 0.0156 |
| GSE148071 | TACSTD2 | 1.189 | 0.067 | 62.8 | 4.3 | 4.15 | 0.796 | ≈ 0 | 35/36, *p* = 2.48×10⁻⁷ |
| GSE148071 | CLDN4 | 1.116 | 0.057 | 61.3 | 3.9 | 4.30 | 0.790 | ≈ 0 | 36/36, *p* = 2.91×10⁻¹¹ |

Controls flip as expected: EPCAM is high in malignant / low in immune
(AUROC 0.76–0.90); PTPRC is the reverse (AUROC 0.15–0.22; paired
Wilcoxon 0/9, 0/7, 0/36 patients with malignant > immune).

### 2. Epithelial restriction (not tumor-vs-normal specificity)

TACSTD2 / CLDN4 mark **epithelium**, not malignancy *per se*.

- **GSE131907 tumor-tissue vs normal-lung Epithelial** (same lineage,
  different `Source`): TACSTD2 mean 1.534 vs 1.325 (% pos 81.4 vs 77.5),
  log2FC = 0.21, AUROC = 0.565, *p* = 5.17×10⁻²⁹. CLDN4 mean 1.417 vs
  1.351 (% pos 80.3 vs 81.3), log2FC = 0.069, AUROC = 0.525,
  *p* = 1.47×10⁻⁵. The cell-level *p* is driven by *n* (9,573 vs 3,265);
  the effect size is small. Both genes are already on in normal lung
  epithelium.
- **EMTAB6149 Malignant vs Alveolar**: TACSTD2 1.415 vs 1.395
  (% pos 59.3 vs 68.7), log2FC = 0.021, AUROC = 0.512, *p* = 0.151
  (not significant). CLDN4 is slightly *higher* in Alveolar
  (1.190 vs 1.313, log2FC = −0.142, *p* = 0.00267). Both remain far
  above immune (epi-vs-immune log2FC 3.79 / 4.21, AUROC 0.819 / 0.819).
- **GSE148071 Malignant vs leftover epithelium** (Alveolar+Basal+Epithelial,
  *n* = 8,147): TACSTD2 1.189 vs 1.014 (% pos 62.8 vs 62.3),
  log2FC = 0.229, AUROC = 0.55, *p* = 2.25×10⁻²⁸. CLDN4 is *higher* in
  leftover epithelium (1.116 vs 1.368, log2FC = −0.293). Epi-vs-immune
  log2FC 3.92 / 4.60.

So the malignant-vs-immune contrast is an **epithelial-vs-immune**
contrast. Residual lung epithelium (AT2/alveolar, basal) expresses
TACSTD2/CLDN4 at malignant-like levels.

### 3. Composition of the gene-positive population

Among cells with count > 0:

| Dataset | Gene | % of + cells in malignant | % in leftover epi | % in immune |
|---|---|---:|---:|---:|
| GSE131907 | TACSTD2 | 76.1 | 0 (no leftover epi in the tumor-tissue slice) | 21.6 |
| GSE131907 | CLDN4 | 80.2 | 0 | 18.7 |
| EMTAB6149 | TACSTD2 | 53.0 | 14.4 | 28.4 |
| EMTAB6149 | CLDN4 | 59.5 | 15.2 | 22.4 |
| GSE127465 | TACSTD2 | 90.3 | — | 8.9 |
| GSE127465 | CLDN4 | 93.7 | — | 6.1 |
| GSE148071 | TACSTD2 | 82.2 | 13.8 | 2.3 |
| GSE148071 | CLDN4 | 80.0 | 16.3 | 2.1 |

The non-zero immune share of TACSTD2+/CLDN4+ cells in GSE131907 and
EMTAB6149 is almost entirely low-level dropout / ambient RNA: immune
% pos is 0.5–5.9% versus 50–81% in malignant, and PTPRC+ cells are
97–98% immune. GSE127465 (higher neutrophil content, cleaner separation)
shows immune % pos ≤ 0.9%.

### 4. Versus immune fractions (per-patient)

Spearman correlation of **patient immune-cell fraction** vs **malignant
compartment mean** (patients with ≥20 malignant and ≥20 immune cells):

| Dataset | *n* patients | TACSTD2 ρ (*p*) | CLDN4 ρ (*p*) |
|---|---:|---|---|
| GSE131907 | 9 | −0.333 (0.381) | +0.283 (0.460) |
| GSE127465 | 7 | −0.107 (0.819) | −0.286 (0.535) |
| GSE148071 | 36 | −0.093 (0.591) | −0.055 (0.752) |

No dataset shows a significant association. Tumor TACSTD2/CLDN4 does
**not** track immune infiltrate fraction, so the malignant-vs-immune
contrast is not an artifact of sample composition.

### Methods (short)

- TISCH2 major-lineage `Malignant` when present; else tumor-tissue
  `Epithelial` (GSE131907 only), labelled as a proxy.
- Immune = TISCH2 major lineages in
  {CD8T, CD8Tex, CD4Tconv, Treg, Tprolif, NK, B, Plasma, Mono/Macro,
  DC, Mast, Neutrophils, …}.
- Cell-level: two-sided Mann–Whitney U, AUROC = *U*/(*n₁n₂*), BH-FDR
  across all gene × lineage tests inside a dataset
  (`stats_malignant_vs_immune.csv`).
- Patient-level: mean of malignant cells vs mean of immune cells,
  paired Wilcoxon signed-rank (minimum 5 patients, ≥20 cells/compartment).
- Immune-fraction: Spearman of `n_immune / n_cells` vs malignant mean.
- Epithelial restriction: leftover {Epithelial, Alveolar, AT1, AT2,
  Basal, Ciliated, Club} vs malignant vs immune; plus tumor-vs-normal
  Epithelial on the un-restricted GSE131907 object.

### Caveats

- GSE131907 has **no TISCH2 malignancy call**. Tumor-tissue Epithelial
  is a proxy (LUAD tumor epithelium is mostly cancer cells, but some
  residual AT2/club cells are mixed in).
- Cell-level tests are pseudoreplicated; the paired Wilcoxon is the
  confirmatory test and is available for 3/4 datasets.
- TISCH2 annotations are MAESTRO-transferred, not the original paper
  labels. GSE148071 `Tprolif` is called stromal in the malignancy
  column but is treated as immune by lineage name (634 cells).
- Ambient RNA can put a thin TACSTD2/CLDN4 tail on immune cells;
  % pos and AUROC, not the raw *p*, are the relevant effect sizes.
- No LuCA, no 12 GB download, no CELLxGENE objects were required.

### Reproduce

```bash
bash scripts/fable_tisch/run_all.sh data_fable_tisch
```

Outputs: `results/fable_tisch/<DATASET>/` and
`results/fable_tisch/cross_dataset_summary.csv`.

---

## 中文

### 问题

公开 NSCLC 单细胞图谱中，**TACSTD2（TROP2）** 与 **CLDN4** 是否主要在
肿瘤上皮 / 恶性细胞中表达、在免疫细胞中基本缺失？该差异是否只是免疫
浸润比例造成的假象？

### 对象（均 <2 GB 或可子集化；未下载 12 GB LuCA）

| 数据集 | 恶性 / 代理 | 免疫 | 残留非恶性上皮 |
|---|---:|---:|---|
| NSCLC_GSE131907（Kim 2020；仅肿瘤组织） | 9,573 肿瘤组织 Epithelial（代理；TISCH2 无 Malignant 标注） | 53,114 | 子集化后无残留上皮；另与 3,265 个**正常肺** Epithelial 比较 |
| NSCLC_EMTAB6149（Lambrechts 2018） | 5,728 Malignant | 30,694 | 1,346 Alveolar |
| NSCLC_GSE127465（Zilionis 2019） | 3,995 Malignant | 26,055 | 无 |
| NSCLC_GSE148071（Wu 2021） | 48,118 Malignant | 20,077 | 8,147 Alveolar + Basal + Epithelial |

对照基因：EPCAM（上皮阳性）、PTPRC/CD45（免疫阳性）。表达值为 TISCH2
`log2(TPM/10 + 1)`。

### 1. 恶性 vs 汇集免疫细胞

四个数据集方向一致：TACSTD2 / CLDN4 在恶性（或肿瘤上皮代理）中的
阳性率 50–81%，免疫细胞中 0.5–5.9%；log2FC 3.8–6.9；AUROC 0.75–0.90。
细胞水平 Mann–Whitney *p* 因样本量过大下溢为 0。

患者水平配对 Wilcoxon（同一患者恶性均值 vs 免疫均值）：

- GSE131907：9/9 患者恶性 > 免疫，*p* = 0.00391（TACSTD2 与 CLDN4 相同）
- GSE127465：7/7，*p* = 0.0156
- GSE148071：TACSTD2 35/36，*p* = 2.48×10⁻⁷；CLDN4 36/36，*p* = 2.91×10⁻¹¹
- EMTAB6149：无 Patient 列，无法做配对检验

对照成立：EPCAM 与靶基因同向；PTPRC 反向（配对 Wilcoxon 恶性 > 免疫
的患者数为 0/9、0/7、0/36）。

### 2. 上皮限制，而非肿瘤相对正常上皮的特异性

- **GSE131907 肿瘤组织 vs 正常肺 Epithelial**：TACSTD2 均值 1.534 vs
  1.325（阳性率 81.4% vs 77.5%），log2FC = 0.21，AUROC = 0.565。
  CLDN4 均值 1.417 vs 1.351（80.3% vs 81.3%），log2FC = 0.069，
  AUROC = 0.525。效应量小，正常肺上皮已经高表达。
- **EMTAB6149 恶性 vs 肺泡上皮**：TACSTD2 1.415 vs 1.395，
  log2FC = 0.021，*p* = 0.151（不显著）。CLDN4 在肺泡中略高
  （log2FC = −0.142，*p* = 0.00267）。两者相对免疫的 log2FC 仍为
  3.79 / 4.21。
- **GSE148071 恶性 vs 残留上皮**（*n* = 8,147）：TACSTD2 log2FC = 0.229
  （阳性率 62.8% vs 62.3%）；CLDN4 在残留上皮更高（log2FC = −0.293）。
  残留上皮 vs 免疫的 log2FC 为 3.92 / 4.60。

结论：恶性 vs 免疫的差异本质是**上皮 vs 免疫**。残留肺泡 / 基底上皮
的 TACSTD2/CLDN4 水平与恶性细胞相近。

### 3. 阳性细胞的区室构成

TACSTD2+ / CLDN4+ 细胞主要落在恶性 + 残留上皮。免疫区室占阳性细胞的
比例在 GSE148071 仅 2.1–2.3%，在 GSE127465 为 6–9%，在 GSE131907 /
EMTAB6149 较高（19–28%），但免疫细胞本身的阳性率只有 0.5–5.9%，
符合环境 RNA / dropout，而不是免疫细胞真表达。PTPRC+ 细胞 84–98%
落在免疫区室，方向相反。

### 4. 与免疫浸润比例

患者免疫细胞比例 vs 恶性区室基因表达均值的 Spearman 相关：

- GSE131907（*n* = 9）：TACSTD2 ρ = −0.333，*p* = 0.381；CLDN4 ρ = 0.283，*p* = 0.460
- GSE127465（*n* = 7）：ρ = −0.107 / −0.286，*p* = 0.819 / 0.535
- GSE148071（*n* = 36）：ρ = −0.093 / −0.055，*p* = 0.591 / 0.752

均不显著。肿瘤 TACSTD2/CLDN4 并不随免疫浸润比例变化，因此
恶性–免疫差异不是样本组成假象。

### 方法与限制

方法见英文 Methods。主要限制：GSE131907 无 TISCH2 恶性标注，使用
肿瘤组织 Epithelial 作为代理；细胞水平检验存在伪重复，以配对
Wilcoxon 为准；TISCH2 为 MAESTRO 迁移标注；未下载 LuCA。

复现：`bash scripts/fable_tisch/run_all.sh data_fable_tisch`。
数字来源：`results/fable_tisch/*/stats_summary.json` 与
`results/fable_tisch/cross_dataset_summary.csv`。
