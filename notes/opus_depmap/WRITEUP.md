# TACSTD2 / CLDN4 in DepMap lung lines

**Release:** DepMap Public 24Q4 (figshare DOI `10.6084/m9.figshare.27993248`).  
**Open files only.** Every number below is written by `scripts/opus_depmap/` into `results/opus_depmap/tables/`. No statistic in this note was typed from memory.

Scripts, in order: `00_download_depmap.py` → `01_build_cohort.py` → `02_dependency.py` → `03_codependency.py` → `03b_family_and_sex.py` → `04_immune.py`.

---

## 中文

### 问题

在公开的 DepMap / Project Achilles CRISPR 数据与 CCLE RNA-seq 中，肺来源癌细胞系对 **TACSTD2**（TROP2）和 **CLDN4** 是否依赖？二者的 Chronos 谱是否与已知伴侣共依赖？它们的表达或依赖分数是否与 IFN / MHC-I / 肿瘤细胞免疫配体共变？

### 队列

`Model.csv` 中 `OncotreeLineage == Lung`，去掉 `Non-Cancerous`，得到 254 个肺癌模型。

| 子集 | n |
| --- | --- |
| 肺癌模型 | 254 |
| 有 CRISPR（Chronos gene effect） | 126（NSCLC 98，SCLC 25，其他 3） |
| 有表达 log2(TPM+1) | 208（NSCLC 143，SCLC 59，其他 6） |
| CRISPR 与表达都有 | 123（NSCLC 95） |
| 全谱系 CRISPR 对照集 | 1178 个模型 × 17916 个基因 |

来源：`cohort_summary.csv`。

Chronos 标尺：0 = 非必需对照中位数，−1 = 共同必需基因中位数。肺系中 726 个非必需对照基因的均值 SD = 0.116；1242 个共同必需对照均值 = −1.077（`02_dependency.py`）。DepMap 惯例：依赖概率 > 0.5 计为依赖系。

### 1. 依赖性（CRISPR gene effect）

| 基因 | 肺系 n | 均值 (SD) | 最小 | 依赖系（p>0.5） | 全癌依赖系 | 相对非必需对照的 z |
| --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | 126 | −0.0287 (0.121) | −0.408 | **0 / 126** | 1 / 1178 | −0.11 |
| CLDN4 | 126 | +0.0342 (0.132) | −0.294 | **0 / 126** | 0 / 1178 | +0.43 |

相对 0 的 Wilcoxon：TACSTD2 *p* = 0.00208，CLDN4 *p* = 0.000893。样本量下这个检验过强：效应落在非必需对照云团内（TACSTD2 肺均值在 17916 个基因中位于第 52 百分位，CLDN4 第 74 百分位；`dependency_genomewide_calibration.csv`）。肺 vs 非肺 Mann-Whitney：TACSTD2 *p* = 0.593，Cliff δ = −0.029；CLDN4 *p* = 0.969，δ = −0.002。

阳性对照行为正常：EEF2 肺均值 −2.28，125/126 系依赖；KRAS −0.70，69/126 依赖，且肺比非肺更负（*p* = 3.6×10⁻⁴，δ = −0.19）；NKX2-1 肺富集 *p* = 1.12×10⁻⁶。

自身表达能否预测自身 gene effect（Spearman，n = 123）：TACSTD2 *r* = −0.159（95% CI −0.327 至 0.019，*p* = 0.079）；CLDN4 *r* = +0.159（−0.019 至 0.327，*p* = 0.079）。高表达四分位并不构成依赖亚群。

**结论：** 在 24Q4 肺系中，TACSTD2 与 CLDN4 都不是 CRISPR 依赖。均值贴近 0，没有任何一系的依赖概率超过 0.5。

### 2. 共依赖

技术底噪（同一模型独立筛两次，肺系 n = 32）：TACSTD2 两次筛间 Pearson *r* = 0.430（*p* = 0.014）；CLDN4 *r* = 0.343（*p* = 0.054）。共同必需对照中位 *r* = 0.234，非必需对照中位 *r* = 0.134。两条查询基因的谱有弱但可测的重复性，远低于 EGFR / CTNNB1 一类强信号基因。

全基因组扫描（Pearson，BH FDR；置换 200 次只对两个靶基因）：

| 查询 | 队列 | n 系 | 最大 \|r\| | 置换零假设 95 分位 | 置换 *p* | FDR<5% 基因数 |
| --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | 全癌 | 1178 | 0.253 | 0.187 | 0.005 | 3328 |
| TACSTD2 | 肺 | 126 | 0.539 | 0.526 | 0.040 | 2 |
| CLDN4 | 全癌 | 1178 | 0.400 | 0.186 | 0.005 | 5630 |
| CLDN4 | 肺 | 126 | 0.506 | 0.526 | 0.095 | 165 |

全癌 FDR 计数被 n = 1178 放大：\|r\| ≈ 0.06 即可过 BH，不能当生物学信号读。真正可解释的靶向结果：

- **CLDN4–CLDN3**（预指定旁系）：全癌 *r* = 0.344，n = 1178，*p* = 4.90×10⁻³⁴；肺 *r* = 0.356，n = 126，*p* = 4.22×10⁻⁵，*q* = 3.94×10⁻⁴（`codependency_targeted_family.csv`）。
- TACSTD2 没有同等强度的伴侣。与 EPCAM 肺 *r* = 0.216，n = 126，*p* = 0.015；与 CLDN4 肺 *r* ≈ 0。
- 同一代码路径的阳性对照：CTNNB1–TCF7L2 全癌 *r* = 0.681；EGFR–GRB2 *r* = 0.458，并带回 SHC1 / GAB1；MYC–MAX *r* = 0.34。扫描本身没有坏。

供体性别不是肺 gene effect 的混杂：TACSTD2 男 91 / 女 33，Mann-Whitney *p* = 0.211；CLDN4 *p* = 0.088。与 7 个 Y 连锁基因 gene effect 的平均 *r* 分别为 −0.006 与 +0.023。

### 3. 表达及依赖 vs IFN / MHC-I / 免疫基因

预指定面板（全部存在于表达与 CRISPR 矩阵）：MHC-I 抗原呈递 16 个基因、IFN 信号 12 个、ISG 14 个、肿瘤细胞免疫配体 9 个，外加上皮阳性对照与淋巴系阴性对照。MHC-I 综合分 = 16 个 MHC-I 基因的列标准化均值；IFN 应答综合分 = 14 个 ISG 的同样处理。

**NSCLC vs SCLC 是强混杂。** 表达（Mann-Whitney）：

| 特征 | n NSCLC | n SCLC | 均值 NSCLC | 均值 SCLC | *p* | Cliff δ |
| --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 表达 | 143 | 59 | 5.10 | 1.09 | 7.18×10⁻¹⁴ | 0.67 |
| CLDN4 表达 | 143 | 59 | 6.09 | 5.25 | 0.00363 | 0.26 |
| MHC-I 综合分 | 143 | 59 | +0.24 | −0.58 | 2.31×10⁻¹² | 0.63 |
| IFN 应答综合分 | 143 | 59 | +0.16 | −0.38 | 4.80×10⁻¹⁰ | 0.56 |

SCLC 同时低表达 TACSTD2 和 MHC-I / ISG，所以**不校正谱系的肺系相关会被放大**。下面同时报告全肺、NSCLC-only 和偏相关（控制 NSCLC 指示变量）。

#### A. 表达 vs 表达（Spearman）

综合分：

| 查询 | 对象 | 队列 | n | *r* | 95% CI | *p* | *q* |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | MHC-I 综合分 | 肺 | 208 | 0.390 | 0.268, 0.500 | 5.66×10⁻⁹ | 1.91×10⁻⁸ |
| TACSTD2 | IFN 应答综合分 | 肺 | 208 | 0.489 | 0.378, 0.586 | 6.54×10⁻¹⁴ | 5.23×10⁻¹³ |
| CLDN4 | MHC-I 综合分 | 肺 | 208 | 0.197 | 0.063, 0.325 | 0.00429 | 0.00947 |
| CLDN4 | IFN 应答综合分 | 肺 | 208 | 0.246 | 0.114, 0.370 | 0.000341 | 0.00103 |
| TACSTD2 | MHC-I 综合分 | NSCLC | 143 | 0.145 | −0.020, 0.302 | 0.084 | 0.128 |
| TACSTD2 | IFN 应答综合分 | NSCLC | 143 | 0.341 | 0.187, 0.478 | 3.06×10⁻⁵ | 0.000140 |
| CLDN4 | MHC-I 综合分 | NSCLC | 143 | 0.093 | −0.072, 0.253 | 0.270 | 0.382 |
| CLDN4 | IFN 应答综合分 | NSCLC | 143 | 0.215 | 0.053, 0.367 | 0.00976 | 0.0265 |
| TACSTD2 | MHC-I 综合分 | 偏相关 NSCLC | 208 | 0.180 | — | 0.00929 | 0.0154 |
| TACSTD2 | IFN 应答综合分 | 偏相关 NSCLC | 208 | 0.342 | — | 4.59×10⁻⁷ | 2.26×10⁻⁶ |
| CLDN4 | MHC-I 综合分 | 偏相关 NSCLC | 208 | 0.107 | — | 0.126 | 0.201 |
| CLDN4 | IFN 应答综合分 | 偏相关 NSCLC | 208 | 0.172 | — | 0.0133 | 0.0314 |

单基因（肺，n = 208；*q* 为该查询×队列块内 BH）：

| 查询 | 基因 | *r* | 95% CI | *p* | *q* |
| --- | --- | --- | --- | --- | --- |
| TACSTD2 | CD274 | 0.528 | 0.422, 0.619 | 2.70×10⁻¹⁶ | 2.88×10⁻¹⁵ |
| TACSTD2 | IFIT3 | 0.509 | 0.401, 0.603 | 4.21×10⁻¹⁵ | 3.85×10⁻¹⁴ |
| TACSTD2 | IRF9 | 0.478 | 0.366, 0.577 | 2.70×10⁻¹³ | 1.92×10⁻¹² |
| TACSTD2 | OAS1 | 0.476 | 0.364, 0.575 | 3.68×10⁻¹³ | 2.35×10⁻¹² |
| TACSTD2 | MX1 | 0.429 | 0.311, 0.534 | 1.00×10⁻¹⁰ | 4.28×10⁻¹⁰ |
| TACSTD2 | B2M | 0.348 | 0.223, 0.462 | 2.56×10⁻⁷ | 5.84×10⁻⁷ |
| TACSTD2 | IRF1 | 0.367 | 0.244, 0.480 | 4.78×10⁻⁸ | 1.33×10⁻⁷ |
| TACSTD2 | HLA-A | 0.223 | 0.089, 0.348 | 0.00124 | 0.00158 |
| CLDN4 | CD274 | 0.337 | 0.211, 0.453 | 6.33×10⁻⁷ | 4.05×10⁻⁶ |
| CLDN4 | OAS1 | 0.326 | 0.199, 0.443 | 1.50×10⁻⁶ | 8.73×10⁻⁶ |
| CLDN4 | IRF1 | 0.271 | 0.140, 0.393 | 7.56×10⁻⁵ | 0.000302 |
| CLDN4 | B2M | 0.166 | 0.030, 0.295 | 0.0168 | 0.0269 |

上皮阳性对照按预期更强（CLDN4–CDH1 肺 *r* = 0.713；TACSTD2–KRT19 *r* = 0.718）。淋巴系阴性对照（CD3E / CD8A）在细胞系中接近不表达（均值 < 0.5 log2(TPM+1)），相关即使过 FDR 也不应解读为 T 细胞浸润。

NSCLC-only（n = 143）时，TACSTD2 与 MHC-I 综合分不再显著（*p* = 0.084），但与 IFN 应答综合分、IRF9（*r* = 0.430，*p* = 8.34×10⁻⁸）、MX1（*r* = 0.402）、CD274（*r* = 0.400）仍显著。CLDN4 与 MHC-I 综合分在 NSCLC 内不显著；与 IFN 应答综合分仍弱相关（*r* = 0.215，*p* = 0.00976，*q* = 0.0265）。

#### B. 依赖（Chronos）vs 免疫基因表达 — 本题的核心检验

128 个预指定检验 × 每个查询 × 队列。**FDR < 5% 的命中数 = 0**（肺 n = 123；NSCLC n = 95；偏相关同样为 0）。

| 查询 | 对象 | 队列 | n | *r* | 95% CI | *p* | *q* |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | MHC-I 综合分 | 肺 | 123 | −0.078 | −0.252, 0.100 | 0.389 | 0.867 |
| TACSTD2 | IFN 应答综合分 | 肺 | 123 | −0.052 | −0.227, 0.126 | 0.566 | 0.929 |
| CLDN4 | MHC-I 综合分 | 肺 | 123 | −0.000 | −0.177, 0.177 | 0.999 | 0.999 |
| CLDN4 | IFN 应答综合分 | 肺 | 123 | 0.129 | −0.049, 0.300 | 0.154 | 0.979 |
| TACSTD2 | MHC-I 综合分 | NSCLC | 95 | 0.021 | −0.182, 0.221 | 0.844 | 0.980 |
| CLDN4 | MHC-I 综合分 | NSCLC | 95 | −0.048 | −0.247, 0.155 | 0.643 | 0.931 |
| TACSTD2 | B2M 表达 | 肺 | 123 | −0.181 | −0.347, −0.004 | 0.0457 | 0.854 |
| TACSTD2 | CD274 表达 | 肺 | 123 | −0.140 | −0.309, 0.038 | 0.122 | 0.854 |

TACSTD2 vs B2M 的未校正 *p* = 0.046 在 128 次检验的 BH 之后消失。负 *r* 本意是“免疫基因越高，Chronos 越负（越依赖）”；这里看不到经多重检验后站得住的效应。这与第 1 节一致：查询基因的 Chronos 谱接近噪声，无法被免疫表达预测。

#### C. 依赖 vs 免疫基因依赖

肺系中唯一经靶向 FDR 站得住的中等效应仍是 **CLDN4–CLDN3**（*r* = 0.356，n = 126，*p* = 4.22×10⁻⁵，*q* = 0.00262）。CLDN4 与 ERAP1（*r* = −0.269，*p* = 0.00237，*q* = 0.0387）、HLA-B（*r* = 0.262，*q* = 0.0387）、JAK1（*r* = −0.261，*q* = 0.0387）的相关刚过 *q* = 0.05，效应小、n = 126，不作机制结论。TACSTD2 在肺系对免疫面板 gene effect 没有 *q* < 0.05 的伴侣。

### 总判断

1. **TACSTD2 与 CLDN4 都不是肺系 CRISPR 依赖。** 0/126 系依赖概率 > 0.5；均值在非必需对照范围内。
2. **CLDN4 的共依赖信号是旁系 CLDN3**，全癌与肺系都在。TACSTD2 没有同等伴侣。扫描管道被 EGFR / CTNNB1 / MYC 阳性对照验证。
3. **表达层面，TACSTD2（以及较弱的 CLDN4）与 IFN 应答 / 部分 MHC-I / CD274 共变。** 相当一部分来自 NSCLC vs SCLC。限制在 NSCLC 后，TACSTD2–IFN 与 TACSTD2–CD274 仍在，TACSTD2–MHC-I 综合分不再显著。
4. **依赖层面，Chronos 分数与 IFN / MHC-I / 免疫配体表达无经 FDR 校正的相关**（128 个预指定检验，0 个 *q* < 0.05）。细胞系 CRISPR 不支持“高 IFN/MHC-I 肺系更依赖 TACSTD2 或 CLDN4”。

### 限制

- 细胞系没有免疫微环境；“免疫基因”在这里是肿瘤细胞的内在转录程序，不是浸润。
- SCLC CRISPR n = 25，SCLC 单独相关只作描述。
- Chronos 对弱效应基因的两次筛间 *r* 只有 0.3–0.4，共依赖扫描的肺系顶命中（Y 连锁基因等）不可信；已用性别检验排除系统偏移。
- 24Q4 是撰写时 figshare 上可稳定拉取的公开整合发布；门户 25Q/26Q 被 Cloudflare 挡住，未使用。

---

## English

### Question

In public DepMap / Achilles CRISPR and CCLE RNA-seq, are lung cancer cell lines dependent on **TACSTD2** (TROP2) or **CLDN4**? Do their Chronos profiles co-depend with known partners? Do their expression or dependency scores track IFN, MHC-I or tumour-cell immune-ligand genes?

### Cohort

Lung Oncotree lineage in `Model.csv`, non-cancerous models dropped: 254 lung cancer models.

| subset | n |
| --- | --- |
| lung cancer models | 254 |
| with CRISPR (Chronos gene effect) | 126 (NSCLC 98, SCLC 25, other 3) |
| with expression, log2(TPM+1) | 208 (NSCLC 143, SCLC 59, other 6) |
| CRISPR and expression | 123 (NSCLC 95) |
| pan-lineage CRISPR matrix | 1178 models × 17,916 genes |

Source: `cohort_summary.csv`.

Chronos scale: 0 = median non-essential control, −1 = median common essential. In these lung lines the between-gene SD of 726 non-essential controls is 0.116; 1,242 common-essential controls average −1.077. DepMap convention: dependency probability > 0.5 counts as a dependent line.

### 1. Dependency

| gene | lung n | mean (SD) | min | dependent lines (p>0.5) | pan-cancer dependent | z vs non-essential controls |
| --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | 126 | −0.0287 (0.121) | −0.408 | **0 / 126** | 1 / 1178 | −0.11 |
| CLDN4 | 126 | +0.0342 (0.132) | −0.294 | **0 / 126** | 0 / 1178 | +0.43 |

Wilcoxon vs 0: TACSTD2 *p* = 0.00208, CLDN4 *p* = 0.000893. With n = 126 this test is overpowered; both locations sit inside the non-essential control cloud (TACSTD2 lung mean = 52nd percentile of 17,916 genes; CLDN4 = 74th; `dependency_genomewide_calibration.csv`). Lung vs non-lung Mann-Whitney: TACSTD2 *p* = 0.593, Cliff’s δ = −0.029; CLDN4 *p* = 0.969, δ = −0.002.

Positive controls behave: EEF2 lung mean −2.28, 125/126 lines dependent; KRAS −0.70, 69/126 dependent and more negative in lung than elsewhere (*p* = 3.6×10⁻⁴, δ = −0.19); NKX2-1 lung-enriched *p* = 1.12×10⁻⁶.

Does a gene’s own expression predict its own gene effect (Spearman, n = 123)? TACSTD2 *r* = −0.159 (95% CI −0.327 to 0.019, *p* = 0.079); CLDN4 *r* = +0.159 (−0.019 to 0.327, *p* = 0.079). The top expression quartile is not a dependent subset.

**Call:** neither TACSTD2 nor CLDN4 is a CRISPR dependency in 24Q4 lung lines. Means sit at 0; no line has dependency probability > 0.5.

### 2. Co-dependency

Technical floor (same model screened twice, lung n = 32): TACSTD2 between-screen Pearson *r* = 0.430 (*p* = 0.014); CLDN4 *r* = 0.343 (*p* = 0.054). Median *r* is 0.234 for common essentials and 0.134 for non-essential controls. The query profiles are weakly reproducible, far below EGFR / CTNNB1.

Genome-wide Pearson scan with BH FDR (200 label permutations on the two targets only):

| query | cohort | n lines | max \|r\| | permutation-null 95th pct | perm *p* | genes at FDR<5% |
| --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | pan-cancer | 1178 | 0.253 | 0.187 | 0.005 | 3328 |
| TACSTD2 | lung | 126 | 0.539 | 0.526 | 0.040 | 2 |
| CLDN4 | pan-cancer | 1178 | 0.400 | 0.186 | 0.005 | 5630 |
| CLDN4 | lung | 126 | 0.506 | 0.526 | 0.095 | 165 |

Pan-cancer FDR counts are inflated by n = 1178 (\|r\| ≈ 0.06 already clears BH) and should not be read as biology. The pre-specified, interpretable result:

- **CLDN4–CLDN3** (paralog): pan-cancer *r* = 0.344, n = 1178, *p* = 4.90×10⁻³⁴; lung *r* = 0.356, n = 126, *p* = 4.22×10⁻⁵, *q* = 3.94×10⁻⁴ (`codependency_targeted_family.csv`).
- TACSTD2 has no partner of that strength. vs EPCAM in lung: *r* = 0.216, n = 126, *p* = 0.015; vs CLDN4 in lung: *r* ≈ 0.
- Same code path, positive controls: CTNNB1–TCF7L2 pan-cancer *r* = 0.681; EGFR–GRB2 *r* = 0.458 plus SHC1 / GAB1; MYC–MAX *r* = 0.34. The scan is not broken.

Donor sex does not structure lung gene effect: TACSTD2 male 91 / female 33, Mann-Whitney *p* = 0.211; CLDN4 *p* = 0.088. Mean *r* vs seven Y-linked gene-effect profiles is −0.006 and +0.023.

### 3. Expression and dependency vs IFN / MHC-I / immune genes

Pre-specified panels (all present in both matrices): 16 MHC-I antigen-presentation genes, 12 IFN-signalling genes, 14 ISGs, 9 tumour-cell immune ligands, plus epithelial positive controls and lymphoid negative controls. MHC-I score = mean column-wise z of the 16 MHC-I genes; IFN-response score = the same for the 14 ISGs.

**NSCLC vs SCLC is a large confounder.** Expression (Mann-Whitney):

| feature | n NSCLC | n SCLC | mean NSCLC | mean SCLC | *p* | Cliff’s δ |
| --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 expression | 143 | 59 | 5.10 | 1.09 | 7.18×10⁻¹⁴ | 0.67 |
| CLDN4 expression | 143 | 59 | 6.09 | 5.25 | 0.00363 | 0.26 |
| MHC-I score | 143 | 59 | +0.24 | −0.58 | 2.31×10⁻¹² | 0.63 |
| IFN-response score | 143 | 59 | +0.16 | −0.38 | 4.80×10⁻¹⁰ | 0.56 |

SCLC is low for TACSTD2 and for MHC-I / ISG, so **unadjusted lung-wide correlations are inflated**. Lung, NSCLC-only, and partial Spearman (NSCLC indicator) are all reported.

#### A. Expression vs expression (Spearman)

Composite scores:

| query | partner | cohort | n | *r* | 95% CI | *p* | *q* |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | MHC-I score | lung | 208 | 0.390 | 0.268, 0.500 | 5.66×10⁻⁹ | 1.91×10⁻⁸ |
| TACSTD2 | IFN-response score | lung | 208 | 0.489 | 0.378, 0.586 | 6.54×10⁻¹⁴ | 5.23×10⁻¹³ |
| CLDN4 | MHC-I score | lung | 208 | 0.197 | 0.063, 0.325 | 0.00429 | 0.00947 |
| CLDN4 | IFN-response score | lung | 208 | 0.246 | 0.114, 0.370 | 0.000341 | 0.00103 |
| TACSTD2 | MHC-I score | NSCLC | 143 | 0.145 | −0.020, 0.302 | 0.084 | 0.128 |
| TACSTD2 | IFN-response score | NSCLC | 143 | 0.341 | 0.187, 0.478 | 3.06×10⁻⁵ | 0.000140 |
| CLDN4 | MHC-I score | NSCLC | 143 | 0.093 | −0.072, 0.253 | 0.270 | 0.382 |
| CLDN4 | IFN-response score | NSCLC | 143 | 0.215 | 0.053, 0.367 | 0.00976 | 0.0265 |
| TACSTD2 | MHC-I score | partial NSCLC | 208 | 0.180 | — | 0.00929 | 0.0154 |
| TACSTD2 | IFN-response score | partial NSCLC | 208 | 0.342 | — | 4.59×10⁻⁷ | 2.26×10⁻⁶ |
| CLDN4 | MHC-I score | partial NSCLC | 208 | 0.107 | — | 0.126 | 0.201 |
| CLDN4 | IFN-response score | partial NSCLC | 208 | 0.172 | — | 0.0133 | 0.0314 |

Per-gene (lung, n = 208; *q* is BH inside that query × cohort block):

| query | gene | *r* | 95% CI | *p* | *q* |
| --- | --- | --- | --- | --- | --- |
| TACSTD2 | CD274 | 0.528 | 0.422, 0.619 | 2.70×10⁻¹⁶ | 2.88×10⁻¹⁵ |
| TACSTD2 | IFIT3 | 0.509 | 0.401, 0.603 | 4.21×10⁻¹⁵ | 3.85×10⁻¹⁴ |
| TACSTD2 | IRF9 | 0.478 | 0.366, 0.577 | 2.70×10⁻¹³ | 1.92×10⁻¹² |
| TACSTD2 | OAS1 | 0.476 | 0.364, 0.575 | 3.68×10⁻¹³ | 2.35×10⁻¹² |
| TACSTD2 | MX1 | 0.429 | 0.311, 0.534 | 1.00×10⁻¹⁰ | 4.28×10⁻¹⁰ |
| TACSTD2 | B2M | 0.348 | 0.223, 0.462 | 2.56×10⁻⁷ | 5.84×10⁻⁷ |
| TACSTD2 | IRF1 | 0.367 | 0.244, 0.480 | 4.78×10⁻⁸ | 1.33×10⁻⁷ |
| TACSTD2 | HLA-A | 0.223 | 0.089, 0.348 | 0.00124 | 0.00158 |
| CLDN4 | CD274 | 0.337 | 0.211, 0.453 | 6.33×10⁻⁷ | 4.05×10⁻⁶ |
| CLDN4 | OAS1 | 0.326 | 0.199, 0.443 | 1.50×10⁻⁶ | 8.73×10⁻⁶ |
| CLDN4 | IRF1 | 0.271 | 0.140, 0.393 | 7.56×10⁻⁵ | 0.000302 |
| CLDN4 | B2M | 0.166 | 0.030, 0.295 | 0.0168 | 0.0269 |

Epithelial positive controls are stronger, as they should be (CLDN4–CDH1 lung *r* = 0.713; TACSTD2–KRT19 *r* = 0.718). Lymphoid negative controls (CD3E / CD8A) are near-absent in these lines (mean < 0.5 log2(TPM+1)); any FDR-passing correlation there is not T-cell infiltrate.

Inside NSCLC only (n = 143), TACSTD2 vs MHC-I score is no longer significant (*p* = 0.084), but TACSTD2 vs IFN-response score, IRF9 (*r* = 0.430, *p* = 8.34×10⁻⁸), MX1 (*r* = 0.402) and CD274 (*r* = 0.400) remain. CLDN4 vs MHC-I score is null in NSCLC; vs IFN-response score it stays weak (*r* = 0.215, *p* = 0.00976, *q* = 0.0265).

#### B. Dependency (Chronos) vs immune-gene expression — the test asked for

128 pre-specified tests per query × cohort. **Hits at FDR < 5% = 0** (lung n = 123; NSCLC n = 95; partial-NSCLC likewise 0).

| query | partner | cohort | n | *r* | 95% CI | *p* | *q* |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TACSTD2 | MHC-I score | lung | 123 | −0.078 | −0.252, 0.100 | 0.389 | 0.867 |
| TACSTD2 | IFN-response score | lung | 123 | −0.052 | −0.227, 0.126 | 0.566 | 0.929 |
| CLDN4 | MHC-I score | lung | 123 | −0.000 | −0.177, 0.177 | 0.999 | 0.999 |
| CLDN4 | IFN-response score | lung | 123 | 0.129 | −0.049, 0.300 | 0.154 | 0.979 |
| TACSTD2 | MHC-I score | NSCLC | 95 | 0.021 | −0.182, 0.221 | 0.844 | 0.980 |
| CLDN4 | MHC-I score | NSCLC | 95 | −0.048 | −0.247, 0.155 | 0.643 | 0.931 |
| TACSTD2 | B2M expression | lung | 123 | −0.181 | −0.347, −0.004 | 0.0457 | 0.854 |
| TACSTD2 | CD274 expression | lung | 123 | −0.140 | −0.309, 0.038 | 0.122 | 0.854 |

The unadjusted TACSTD2–B2M *p* = 0.046 does not survive BH over 128 tests. A negative *r* would have meant “higher immune expression, more negative (more dependent) Chronos”; nothing of that form remains after FDR. This matches section 1: the query Chronos profiles are near noise and are not predicted by immune expression.

#### C. Dependency vs immune-gene dependency

The only moderate, targeted-FDR-passing effect in lung is again **CLDN4–CLDN3** (*r* = 0.356, n = 126, *p* = 4.22×10⁻⁵, *q* = 0.00262). CLDN4 vs ERAP1 (*r* = −0.269, *p* = 0.00237, *q* = 0.0387), HLA-B (*r* = 0.262, *q* = 0.0387) and JAK1 (*r* = −0.261, *q* = 0.0387) just clear *q* = 0.05; effects are small and n = 126, so they are not treated as mechanism. TACSTD2 has no *q* < 0.05 partner among immune-panel gene effects in lung.

### Call

1. **Neither TACSTD2 nor CLDN4 is a CRISPR dependency in lung lines.** 0/126 lines with probability > 0.5; means sit in the non-essential range.
2. **CLDN4 co-depends with its paralog CLDN3** in both pan-cancer and lung. TACSTD2 has no partner of that strength. EGFR / CTNNB1 / MYC positive controls validate the scan.
3. **At the expression level, TACSTD2 (and more weakly CLDN4) co-varies with IFN-response genes, some MHC-I genes, and CD274.** A large share of the lung-wide signal is NSCLC vs SCLC. Inside NSCLC, TACSTD2–IFN and TACSTD2–CD274 remain; TACSTD2–MHC-I score does not.
4. **At the dependency level, Chronos scores do not track IFN / MHC-I / immune-ligand expression after FDR** (128 pre-specified tests, 0 with *q* < 0.05). Cell-line CRISPR does not support “high-IFN / high-MHC-I lung lines depend on TACSTD2 or CLDN4”.

### Limits

- Cell lines have no immune microenvironment; “immune genes” here are tumour-intrinsic programmes, not infiltrate.
- SCLC CRISPR n = 25; SCLC-only correlations are descriptive.
- Between-screen *r* for these weak-effect genes is only 0.3–0.4, so lung-only genome-wide top hits (including Y-linked genes) are not trusted; a sex check found no systematic shift.
- 24Q4 is the open integrated release that figshare served stably. Portal 25Q/26Q endpoints were behind Cloudflare and were not used.

### Files

| path | content |
| --- | --- |
| `results/opus_depmap/tables/cohort_*.csv` | cohort |
| `results/opus_depmap/tables/dependency_*.csv` | CRISPR dependency |
| `results/opus_depmap/tables/codependency_*.csv` | co-dependency, family, sex |
| `results/opus_depmap/tables/immune_*.csv` | IFN / MHC-I / ligand tests |
| `results/opus_depmap/figures/fig1_dependency_distributions.png` | gene-effect vs controls |
| `results/opus_depmap/figures/fig2_codependency.png` | co-dependency |
| `results/opus_depmap/figures/fig3_immune.png` | expression and dependency vs immune panel |
