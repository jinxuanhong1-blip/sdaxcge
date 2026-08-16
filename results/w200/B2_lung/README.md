# B2 · DepMap/CCLE lung lines only · TACSTD2–CLDN4 ρ

**Honest verdict: the user-claimed ρ = 0.69 is not reproduced.**

On the pre-specified public slice — **DepMap Public 24Q4**, `OncotreeLineage == Lung`, complete RNA cases — Spearman **ρ = 0.607** (n = 214). That rounds to **0.61**, not 0.69. Pearson r = 0.562. A 5,000-resample bootstrap 95% CI is **[0.513, 0.687]**; **0.69 is outside that interval**.

We did not change lineage, histology, or correlation method to chase 0.69.

---

## 中文

### 一句话结论

用户给出的 **ρ = 0.69** 在公开 DepMap/CCLE **肺细胞系**上**不能复现**。预指定主分析（24Q4，Oncotree 谱系 = Lung，有 RNA 的细胞系）Spearman **ρ = 0.6075（n = 214）**，四舍五入为 **0.61**；Pearson r = 0.56。0.69 落在主分析 bootstrap 95% CI **之外**。

阳性共表达本身是真的（p = 5.6×10⁻²³），只是**点估计不是 0.69**。把 SCLC/神经内分泌系和 NSCLC 混在一起会把 ρ 拉低；仅 NSCLC 为 0.67，仍然不是 0.69。

### 用户数字 vs 本计算

| 项目 | 值 |
|---|---|
| 用户声称 | ρ = 0.69（未标明 Spearman/Pearson、版本、是否含 SCLC） |
| 预指定主分析 | 24Q4 Lung 细胞系，Spearman |
| 主分析点估计 | **0.6075 → 0.61** |
| 与 0.69 的差 | 0.083 |
| 是否在 2 位小数上相等 | **否** |
| 是否在 \|Δ\| ≤ 0.05 的“邻近”带 | **否** |
| bootstrap 95% CI | [0.513, 0.687] |
| 0.69 是否落在该 CI 内 | **否** |
| Pearson（同一队列） | 0.562 → 0.56（更远） |

### 数据

| 文件 | 来源 | 用途 |
|---|---|---|
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | [DepMap 24Q4 Figshare+](https://doi.org/10.25452/figshare.plus.27993248.v1) file 51065489 | log2(TPM+1)；列 `TACSTD2 (4070)`、`CLDN4 (1364)` |
| `Model.csv` | 同上 file 51065297 | `OncotreeLineage` / `OncotreePrimaryDisease` / `OncotreeSubtype` |
| `CCLE_expression.csv` + `sample_info.csv` | DepMap 22Q2 Figshare files 34989919 / 35020903 | 历史敏感性（旧 lineage=lung **含间皮瘤**） |
| `CCLE_RNAseq_genes_rpkm_20180929.gct.gz` | [Broad CCLE 2018](https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz) | 历史 RPKM；基因 `ENSG00000184292` / `ENSG00000189143` |

24Q4 表达矩阵 SHA256：`2a71dc94110efcc0221eae821bb93a9f03b54bea16f005818911a09d33383d56`（见 `download_manifest.json`）。门户网站有人机验证，本分析走 Figshare+，不爬 portal。

### 方法（预指定，不事后改）

1. 主队列：**OncotreeLineage == Lung** 且表达矩阵中 TACSTD2、CLDN4 均为数值。24Q4 中该谱系的 `ModelType` 全部是 Cell Line（260 个模型里 214 个有 RNA；46 个无 24Q4 RNA，已排除）。
2. 统计：`scipy.stats.spearmanr`（主）与 `pearsonr`；双侧；完整病例。Spearman 95% CI：对样本有放回重抽样 5000 次，种子 0，百分位法。
3. 敏感性（**报告但不用来“凑”0.69**）：NSCLC / NET-SCLC / LUAD / LUSC；22Q2 `lineage==lung`；2018 RPKM 的 `_LUNG` 后缀与 `Site_Primary==lung`；全谱系对照。
4. **没有**按表达阈值过滤、没有去掉 A549 等低表达系、没有改用蛋白组来贴近 0.69。

### 结果

**主分析（图 `fig_scatter_lung_cell_lines.png`）**

- n = 214；Spearman ρ = **0.607**（p = 5.6×10⁻²³）；Pearson r = 0.562。
- NSCLC（蓝，n = 143）多在双高象限；肺神经内分泌/SCLC（红，n = 60）大量堆在 TACSTD2 ≈ 0，CLDN4 从 0 到 9 不等。这就是全肺 ρ 低于 NSCLC-only 的原因。

**按组织学（同一 24Q4 矩阵）**

| 队列 | n | Spearman ρ | 四舍五入 | vs 0.69 |
|---|---:|---:|---:|---|
| Lung 全部（主） | 214 | 0.607 | 0.61 | 不匹配 |
| NSCLC | 143 | 0.667 | 0.67 | 邻近，不是 0.69；CI [0.56, 0.75] **含** 0.69 |
| LUAD | 80 | 0.756 | 0.76 | 不匹配（更高） |
| LUSC | 27 | 0.532 | 0.53 | 不匹配 |
| SCLC 亚型 | 59 | 0.348 | 0.35 | 不匹配 |
| NET（含多数 SCLC） | 60 | 0.370 | 0.37 | 不匹配 |
| 全谱系（对照，非本题） | 1673 | 0.721 | 0.72 | 邻近，但**不是肺** |

**历史 CCLE（命名敏感性，不替代主分析）**

| 队列 | n | Spearman ρ | 四舍五入 |
|---|---:|---:|---:|
| 22Q2 lineage=lung（含间皮瘤） | 207 | 0.664 | 0.66 |
| 22Q2 lung 去掉间皮瘤 | 187 | 0.618 | 0.62 |
| 22Q2 NSCLC | 135 | 0.653 | 0.65 |
| 2018 RPKM，ID 后缀 LUNG | 188 | 0.642 | 0.64 |
| 2018 RPKM，Site_Primary=lung | 181 | 0.674 | 0.67 |
| 2018 RPKM **全部组织**（对照） | 1019 | 0.679 | 0.68 |

没有任何一个**预指定的肺队列**在 2 位小数上等于 0.69。最接近的肺切片是 2018 `Site_Primary=lung`（0.67）和 24Q4 NSCLC（0.67）。2018 **全组织** RPKM 为 0.68，比“仅肺”更靠近 0.69——若有人把全 CCLE 相关系数标成“肺”，可能接近用户数字；**这不是对本题的确认**。

### 质控备忘

- 从 507 MB 矩阵抽出的 A549（ACH-000681）TACSTD2 = 0.202、CLDN4 = 2.205，与原始 CSV 逐列核对一致。该系在本 RNA 矩阵上 TACSTD2 极低，与常见蛋白 IHC 印象可以不一致；**未因此删点**。
- 工程株仅 2 个（A549_CRAF_KD、HCC-827-GR5）。去掉后 ρ = 0.603，可忽略。
- 旧版 22Q2 的 `lineage==lung` **包含间皮瘤**；24Q4 的 Oncotree Lung **不含**间皮瘤。混用注释会改变 n 和 ρ。

### 局限

1. 细胞系 ≠ 肿瘤；无基质、无免疫微环境。
2. mRNA ≠ TROP2/CLDN4 蛋白，不能当 ADC 靶点定量。
3. 用户未给出 DepMap 版本、是否 Spearman、是否排除 SCLC。不同合理定义给出 0.35–0.76，**没有一个是 0.69**。
4. 未计算全基因组共表达排名（那是 B1/B2“top coexp”另一题）；本题只回答这对基因的 ρ。

### 复现

```bash
pip install -r scripts/w200/B2_lung/requirements.txt
python3 scripts/w200/B2_lung/download.py
python3 scripts/w200/B2_lung/download_historical.py   # optional 22Q2 + 2018
python3 scripts/w200/B2_lung/analyze.py
```

输出目录：`results/w200/B2_lung/`。

---

## English

### TL;DR

The claimed **ρ = 0.69** for TACSTD2–CLDN4 in **DepMap/CCLE lung lines only** is **not reproduced**. The pre-specified analysis (DepMap Public 24Q4, `OncotreeLineage == Lung`, n = 214 complete RNA cases) gives Spearman **ρ = 0.6075** (rounds to **0.61**), Pearson r = 0.56. The bootstrap 95% CI is **[0.513, 0.687]**; **0.69 lies outside it**.

The two genes are positively correlated. The point estimate is not 0.69.

### What was pre-specified

- **Primary cohort:** lung lineage cell lines in 24Q4 (`OncotreeLineage == Lung`). All 260 lung models in `Model.csv` are `ModelType == Cell Line`; 214 have RNA.
- **Primary statistic:** Spearman ρ on log2(TPM+1). Pearson is reported, not used to pick a winner.
- **Match rule (set before looking at ρ):** equality at 2 decimal places vs 0.69. A looser |Δ| ≤ 0.05 band is labeled “nearby,” not a match.
- **No tuning:** we did not drop SCLC, drop low-expressors, switch releases, or switch to proteomics in order to land on 0.69.

### Main number

| | |
|---|---|
| User claim | 0.69 |
| 24Q4 lung lines Spearman | **0.6075 (n=214, p=5.6e-23)** |
| Rounds to | **0.61** |
| \|Δ\| vs 0.69 | 0.083 |
| Bootstrap 95% CI | [0.513, 0.687] |
| 0.69 in that CI? | **No** |
| Pearson, same n | 0.562 |

NSCLC-only is **0.667** (n=143; nearby; CI [0.56, 0.75] *does* contain 0.69). That is a different question than “lung lines only.” SCLC/NET is **0.37**. Mixing them is why the all-lung estimate sits at 0.61.

All-lineage 24Q4 ρ = **0.72**. CCLE 2018 RPKM *all tissues* ρ = **0.68**. Those are not lung-only results.

### Files

- `summary.json` — machine-readable verdict
- `correlations.csv` — every cohort, including historical
- `lung_cell_lines_expression.csv` — 214-line table used for the primary ρ
- `fig_scatter_lung_cell_lines.png` — primary scatter (NSCLC vs NET vs other)
- `fig_scatter_nsclc_cell_lines.png` — NSCLC-only scatter
- `download_manifest.json` — URLs and SHA256 for 24Q4
- two-gene extracts for 24Q4 / 22Q2 / 2018 so the full matrices need not be stored

### Bottom line

**Do not cite 0.69 as the DepMap/CCLE lung-line TACSTD2–CLDN4 correlation.** Cite **0.61 (n=214, 24Q4, all lung lines)** or, if the claim is restricted to NSCLC, **0.67 (n=143)**. Neither is 0.69.
