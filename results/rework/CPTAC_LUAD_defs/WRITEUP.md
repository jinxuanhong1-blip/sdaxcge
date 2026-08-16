# CPTAC LUAD TACSTD2 protein vs MCP-counter / ESTIMATE / CIBERSORT — WES residual

Self-contained rework. Outputs only under `results/rework/CPTAC_LUAD_defs/`.
PR6 already reported the xCell slice. This folder does **not** replace PR6.
It asks whether that LUAD protein–immune anti-correlation is robust to
**other public immune definitions** after a **DNA purity residual**.

LSCC TACSTD2 protein vs xCell was **null** (PR23 / PR84). This folder is
**LUAD only**.

---

## English

### Question
In treatment-naive CPTAC lung adenocarcinoma, does **TACSTD2 protein**
associate with immune infiltration when immune is defined by public
**ESTIMATE**, public **CIBERSORT**, and **MCP-counter** (Becht 2016),
after residualizing on **WES_purity**?

PR6’s headline used a different definition: TACSTD2 protein vs freeze
**xCell immune score** ρ = −0.309 (p = 0.00100, primary-list q = 0.024).
That number is reproduced here as a check. It is **not** in this folder’s
FDR family.

This cohort has **no ICI response labels** (Gillette et al. *Cell* 2020,
PMID 32649874). Do not read these correlations as immunotherapy outcomes.

### Data
Open S3 freeze `data_freeze_v1.2_reorganized` (`cptac-pancancer-data`,
us-west-2). Protein + RNA + phenotype IDs align at **n = 110 tumors**.
Standalone `LUAD_{xcell,cibersort,mcpcounter,estimate}.txt` files are
HTTP **403**. ESTIMATE, CIBERSORT, xCell, and WES/WGS purity are taken
from `LUAD_phenotype.txt` (HTTP 200).

| Layer | n | Notes |
|---|---:|---|
| TACSTD2 protein `ENSG00000184292.7` | 110 | NA = 0 |
| CLDN4 protein `ENSG00000189143.9` | **79** | **NA = 31 (28.2%)**, TMT dropout |
| ESTIMATE ImmuneScore / CIBERSORT / xCell | 110 | RNA-derived freeze columns |
| MCP-counter (computed) | 110 | public RNA + Becht 2016 genes |
| WES_purity (primary residual) | 108 | DNA; 2 missing |
| WGS_purity (sensitivity) | 104 | DNA |
| ESTIMATE TumorPurity (sensitivity) | 110 | cosine of ESTIMATEScore; RNA-circular |

MCP-counter is **not** a freeze column. Scores are the mean of log2 RSEM
UQ-1500 marker genes (equivalent to log2 geometric mean of linear
expression). Official T-cell (15/15) and cytotoxic-lymphocyte (7/7)
markers are complete. Official “CD8 T cells” is a **single gene**
(`CD8B`) and was kept **exploratory**. NK cells miss `KIR3DS1`
(`ENSG00000275037`). Gene table:
`scripts/rework/CPTAC_LUAD_defs/mcp_counter_genes.tsv`.

CIBERSORT’s 22 freeze columns do **not** sum to 1 (median sum = 1.84,
range 1.31–3.78). They are not relative-mode fractions. `CIBERSORT_T_cell_CD8+`
is used as published.

### Method
Prespecified primary: TACSTD2 protein × 4 scores (ESTIMATE ImmuneScore,
MCP-counter T cells, MCP-counter cytotoxic lymphocytes, CIBERSORT CD8).
Primary residual = **WES_purity**. Partial Spearman = Pearson of
rank-residuals after OLS on ranked purity (df = n−3). BH-FDR across those
**four** WES partial tests only. CLDN4 is a comparator with its own 4-test
FDR (n = 79 / 77). xCell is a reproduction check only.

### Purity context (why residualizing is required)
| Feature | vs WES_purity ρ | p | n |
|---|---:|---:|---:|
| TACSTD2 protein | +0.074 | 0.446 | 108 |
| CLDN4 protein | +0.146 | 0.204 | 77 |
| ESTIMATE ImmuneScore | **−0.523** | 6.2×10⁻⁹ | 108 |
| MCP-counter T cells | **−0.500** | 3.5×10⁻⁸ | 108 |
| MCP-counter cytotoxic | −0.344 | 2.7×10⁻⁴ | 108 |
| CIBERSORT CD8 | −0.235 | 0.014 | 108 |
| xCell immune score | −0.402 | 1.6×10⁻⁵ | 108 |

TACSTD2 protein is **not** a DNA-purity marker. The immune scores are.
Residualizing immune on WES removes a large shared component by
construction. That is the test, not a bug.

### Reproduction of PR6 (not the rework claim)
| Score | n_marg | marg ρ | marg p | n_WES | partial ρ \| WES | partial p | 95% CI |
|---|---:|---:|---:|---:|---:|---:|---|
| xCell immune score | 110 | **−0.309** | **0.00100** | 108 | **−0.272** | **0.0047** | −0.44 to −0.09 |
| xCell CD8 | 110 | **−0.289** | **0.00219** | 108 | **−0.261** | **0.0067** | −0.43 to −0.07 |

These match PR6 to three decimals. Direction holds after WES residual.
This folder does **not** use them for FDR.

### Primary result (the rework)

**TACSTD2 protein is not significantly associated with ESTIMATE
ImmuneScore, MCP-counter T cells, MCP-counter cytotoxic lymphocytes, or
CIBERSORT CD8 after WES residual.** All four WES partial tests have
FDR > 0.36. Marginal ESTIMATE and MCP T cells are only nominal
(p ≈ 0.052) and do not survive multiplicity or purity adjustment.

The PR6 xCell finding is **real for xCell** and **not robust** to these
other public definitions. Do **not** write a deconvolution-general
“high TROP2 protein = immune-cold LUAD” claim from this freeze.

| Score | n_marg | marg ρ | marg p | marg FDR (4) | n_WES | partial ρ \| WES | partial p | FDR (4) | 95% CI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| ESTIMATE ImmuneScore | 110 | −0.185 | 0.052 | 0.105 | 108 | −0.129 | 0.184 | **0.360** | −0.31 to +0.06 |
| MCP-counter T cells | 110 | −0.185 | 0.053 | 0.105 | 108 | −0.134 | 0.170 | **0.360** | −0.32 to +0.06 |
| MCP-counter cytotoxic | 110 | −0.108 | 0.262 | 0.262 | 108 | −0.068 | 0.485 | **0.485** | −0.25 to +0.12 |
| CIBERSORT CD8 | 110 | −0.133 | 0.164 | 0.219 | 108 | −0.108 | 0.270 | **0.360** | −0.29 to +0.08 |

OLS residual Spearman (raw values on WES) agrees: all four p > 0.15.
WGS residual is the same story (all p > 0.24). ESTIMATE-purity residual
is **circular** with ImmuneScore (ImmuneScore vs ESTIMATE purity
ρ = −0.733) and is not an honest covariate for that score.

### CLDN4 protein (comparator, n = 79)
No primary WES partial reaches FDR < 0.05. Strongest marginal is
MCP-counter T cells ρ = −0.243 (p = 0.031); after WES ρ = −0.184
(p = 0.112, FDR = 0.447). CIBERSORT CD8 is null (partial ρ = −0.019,
p = 0.87). CLDN4 protein missingness is TMT dropout, not a join error.

### Exploratory (not the inferential claim)
These were not in the 4-test FDR family. Do not promote them after seeing
the numbers.

- Official MCP-counter **CD8 T cells** is `CD8B` alone. TACSTD2 protein vs
  that single gene is ρ = −0.289 (p = 0.00219), WES partial ρ = −0.263
  (p = 0.0063) — numerically close to xCell CD8. That is a one-gene RNA
  correlation, not a multi-gene MCP-counter score. It was prespecified
  exploratory for that reason.
- MCP-counter NK cells (8/9 genes) ρ = −0.277 (p = 0.0034); WES
  ρ = −0.245 (p = 0.011).
- MCP-counter B lineage ρ = −0.256 (p = 0.0069); WES ρ = −0.262
  (p = 0.0065).
- CIBERSORT leukocyte sum ρ = −0.192 (p = 0.045); WES ρ = −0.149
  (p = 0.125) — does not survive residual.
- ESTIMATE StromalScore, MCP monocytes / neutrophils / fibroblasts,
  CIBERSORT Treg / M2 / neutrophil: all null.

A post-hoc “lymphoid RNA” pattern (CD8B / NK / B) is **not** a
prespecified claim and is not the same as ESTIMATE / CIBERSORT / official
MCP T or cytotoxic scores.

### Honest interpretation
1. **PR6’s xCell number is reproducible** in this freeze (ρ = −0.309;
   WES residual ρ = −0.272, p = 0.0047).
2. **The same TACSTD2 protein is not significantly inverse** with public
   ESTIMATE ImmuneScore, public CIBERSORT CD8, or Becht MCP-counter T /
   cytotoxic scores after WES residual (all FDR ≥ 0.36). Direction is
   weakly negative and consistent with noise plus purity.
3. The LUAD protein–immune story is **definition-dependent**. xCell
   immune / CD8 show it; the three public families requested here do not
   at the prespecified primary endpoints.
4. LSCC TACSTD2 protein was already null vs xCell. Combined with this
   LUAD rework, a histology- and method-general “TROP2 protein = T-cell
   excluded” statement is **not supported**.
5. This is **not** an ICI cohort. Cross-layer (protein vs RNA scores)
   observational association only.

### Caveats
- No ICI labels; OS/PFS were not retested.
- Immune scores are RNA-derived; protein–score tests mix layers.
- MCP-counter was computed here because the freeze file is 403. Marker
  lists are the official Becht table, not a re-fit.
- Official MCP “CD8 T cells” is one gene; it was not a primary endpoint.
- CIBERSORT columns are not relative fractions (row sums ≠ 1).
- WES vs ESTIMATE ImmuneScore ρ = −0.52, so residualizing removes a
  large shared component.
- ESTIMATE purity is circular with ImmuneScore.
- Table-wise BH on 4 TACSTD2 WES partial tests only.
- CLDN4 protein missing 31/110.
- LUAD ≠ LSCC.

### How to rerun
```bash
pip install -r scripts/rework/CPTAC_LUAD_defs/requirements.txt
python3 scripts/rework/CPTAC_LUAD_defs/download.py
python3 scripts/rework/CPTAC_LUAD_defs/analyze.py
```
Matrices: `data/rework/CPTAC_LUAD_defs/` (not committed).
Tables/figures: `results/rework/CPTAC_LUAD_defs/`.
Methods: `results/rework/CPTAC_LUAD_defs/methods.md`.

---

## 中文

### 问题
PR6 报告 CPTAC **初治**肺腺癌中 **TACSTD2 蛋白**与冻存 **xCell immune
score** Spearman ρ = −0.309（p = 0.00100，主键 q = 0.024）。LSCC 的
TACSTD2 蛋白对 xCell **不显著**（PR23 / PR84）。本目录只重做 LUAD，问：
若免疫改用公开 **ESTIMATE**、公开 **CIBERSORT**、以及用公开 RNA 按
Becht 2016 计算的 **MCP-counter**，再对 **WES 纯度**做残差，关联是否还在？

本队列**没有 ICI 疗效标签**（Gillette 等 *Cell* 2020）。不能当成免疫治疗
结局。

### 数据与方法
开放 S3 冻存 `data_freeze_v1.2_reorganized`。肿瘤 110 例，蛋白 / RNA /
表型 ID 对齐。独立 `LUAD_mcpcounter.txt` 等为 HTTP 403。ESTIMATE /
CIBERSORT / xCell / 纯度取自 `LUAD_phenotype.txt`。MCP-counter 用公开
RNA + 官方 marker（T 细胞 15/15、细胞毒 7/7 齐全；官方 “CD8 T cells”
只有 `CD8B` 一基因，故放探索）。CIBERSORT 22 列行和中位数 1.84，**不是**
相对比例。主残差为 **WES_purity**（108/110）。推断用 TACSTD2 × 4 个主键
WES 偏相关的 BH-FDR。xCell 只作复现，不进本目录 FDR。

### 纯度背景
TACSTD2 蛋白与 WES 纯度无关（ρ = +0.074，p = 0.45）。ESTIMATE ImmuneScore
与 WES ρ = −0.523，MCP T 细胞 ρ = −0.500，xCell immune ρ = −0.402。
免疫评分本身就是纯度的反指标，残差会去掉大量共用成分。

### 主要结果（未编造）
- **PR6 的 xCell 数字可复现**：TACSTD2 蛋白 vs xCell immune ρ = −0.309
  （p = 0.00100）；WES 残差后 ρ = −0.272（p = 0.0047）。xCell CD8
  ρ = −0.289 → −0.261。这不是本目录的推断终点。
- **主键重做：ESTIMATE ImmuneScore、MCP-counter T 细胞、MCP-counter
  细胞毒、CIBERSORT CD8 在 WES 残差后均不显著**（偏相关 ρ = −0.129 /
  −0.134 / −0.068 / −0.108；p = 0.18 / 0.17 / 0.49 / 0.27；FDR 全部
  ≥ 0.36）。边缘 ESTIMATE / MCP T 仅名义 p ≈ 0.052，过不了多重检验，
  也过不了纯度校正。
- **不能**把 PR6 写成“对所有公开反卷积都成立的高 TROP2 蛋白 = 免疫冷”。
  xCell 有信号；这里指定的三个公开定义没有。
- CLDN4 蛋白（n = 79，NA = 31）主键 WES 偏相关全部 FDR > 0.44。
- 探索性：单基因 MCP CD8（`CD8B`）ρ = −0.289，WES 后 p = 0.006，和
  xCell CD8 接近，但这是单基因，**不是**主键。NK / B 谱系也有类似名义
  负相关；CIBERSORT 白细胞总和 WES 后 p = 0.13。这些都不能事后升格为主结论。
- 与 LSCC TACSTD2 蛋白对 xCell 的阴性合在一起，**不能**写成跨组织学、
  跨方法的“TROP2 蛋白 = T 细胞排斥”。

### 局限
无 ICI 标签；免疫评分来自 RNA；MCP-counter 为本地按官方基因表计算；
CIBERSORT 非相对比例；WES 与 ImmuneScore 高度负相关；ESTIMATE 纯度与
ImmuneScore 循环；多重检验仅为 4 项主检验；CLDN4 缺失 31 例；不可与
LSCC 未经预设交互而合并。

### 复现
见上文 English “How to rerun”。
