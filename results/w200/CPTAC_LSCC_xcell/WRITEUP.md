# CPTAC LSCC protein TACSTD2/CLDN4 vs xCell CD8/immune — purity residual

Focused recompute. Outputs only under `results/w200/CPTAC_LSCC_xcell/`.
PR23 already analyzed this freeze; this folder reruns the **xCell CD8 / immune**
slice and adds a **DNA-purity residual**. It does not replace PR23.

---

## English

### Question
In treatment-naive CPTAC lung squamous cell carcinoma (LSCC / LUSC), do
**TACSTD2** and **CLDN4 protein** associate with freeze **xCell CD8** and
**xCell immune score** after residualizing on tumor purity?

This cohort has **no ICI response labels**. Do not read these correlations as
immunotherapy outcomes.

### Data
Open S3 freeze `data_freeze_v1.2_reorganized` (`cptac-pancancer-data`, us-west-2).
Protein tumor n=108. xCell and purity columns are in `LSCC_phenotype.txt`
(standalone `LSCC_xcell.txt` is HTTP 403). Satpathy et al. *Cell* 2021,
PMID 34358469; PDC000234.

| Layer | n | Notes |
|---|---:|---|
| TACSTD2 protein `ENSG00000184292.7` | 108 | NA=0 |
| CLDN4 protein `ENSG00000189143.9` | **78** | **NA=30 (27.8%)**, TMT dropout |
| xCell CD8 / immune | 108 | RNA-derived |
| WES_purity (primary residual) | 107 | DNA; 1 missing |
| WGS_purity (sensitivity) | 104 | DNA |
| ESTIMATE TumorPurity (sensitivity) | 108 | cosine of ESTIMATEScore; RNA-circular |

Marginal Spearman for the four primary pairs **reproduces PR23**
(TACSTD2–CD8 ρ=−0.080; TACSTD2–immune ρ=−0.133; CLDN4–CD8 ρ=−0.387;
CLDN4–immune ρ=−0.400).

### Method
Prespecified: 2 proteins × 2 xCell scores. Primary residual = **WES_purity**.
Partial Spearman = Pearson of rank-residuals after OLS on ranked purity
(df = n−3). A second “literal residual” Spearman uses OLS residuals of the
raw values. BH-FDR across the four WES partial tests only.

### Purity context (why residualizing is required)
| Feature | vs WES_purity ρ | p | n |
|---|---:|---:|---:|
| TACSTD2 protein | +0.066 | 0.497 | 107 |
| CLDN4 protein | **+0.305** | 0.0070 | 77 |
| xCell CD8 | −0.396 | 2.4×10⁻⁵ | 107 |
| xCell immune score | **−0.724** | 1.3×10⁻¹⁸ | 107 |

CLDN4 protein tracks DNA purity. xCell immune score is almost a purity
anti-marker. TACSTD2 protein is **not** a purity marker in this set.

### Primary result

**TACSTD2 protein is null** against xCell CD8 and xCell immune, before and
after WES residual. This does **not** support a Bessede-like
“high TROP2 = T-cell excluded” pattern in treatment-naive LSCC protein.

**CLDN4 protein is inverse** with both xCell scores marginally. After WES
residual the inverse **remains FDR&lt;0.05 but is smaller**. After WGS residual
(n=75) it **does not remain**. Do not treat a purity-independent CLDN4–immune
effect as settled.

| Protein | xCell | n_marg | marg ρ | marg p | n_WES | partial ρ \| WES | partial p | FDR (4 tests) | 95% CI |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| TACSTD2 | CD8 | 108 | −0.080 | 0.413 | 107 | −0.061 | 0.534 | 0.534 | −0.25 to +0.13 |
| TACSTD2 | immune | 108 | −0.133 | 0.169 | 107 | −0.129 | 0.188 | 0.251 | −0.31 to +0.06 |
| CLDN4 | CD8 | 78 | −0.387 | 4.7×10⁻⁴ | 77 | **−0.303** | 0.0078 | **0.023** | −0.49 to −0.08 |
| CLDN4 | immune | 78 | −0.400 | 2.9×10⁻⁴ | 77 | **−0.288** | 0.012 | **0.023** | −0.48 to −0.07 |

OLS residual Spearman (raw values) agrees in direction for CLDN4
(CD8 ρ=−0.289, p=0.011; immune ρ=−0.279, p=0.014). TACSTD2 OLS vs immune is
nominally p=0.046 but the rank-partial is p=0.188; that is method-sensitive
and is **not** a primary claim.

### Sensitivity (not the primary FDR family)
| Protein | xCell | partial ρ \| WGS (n) | p | partial ρ \| ESTIMATE (n) | p |
|---|---|---:|---:|---:|---:|
| TACSTD2 | CD8 | −0.033 (104) | 0.737 | −0.056 (108) | 0.565 |
| TACSTD2 | immune | −0.087 (104) | 0.381 | −0.121 (108) | 0.213 |
| CLDN4 | CD8 | −0.173 (75) | **0.141** | −0.315 (78) | 0.0053 |
| CLDN4 | immune | −0.049 (75) | **0.681** | −0.310 (78) | 0.0061 |

WGS residual nullifies CLDN4. ESTIMATE residual does not, but ESTIMATE purity
is built from stromal+immune RNA and is circular with xCell immune. Exploratory
CD8 subsets: CLDN4 vs naive / central-memory stay inverse after WES; effector
memory does not (partial p=0.23). TACSTD2 subsets are all null.

### Honest interpretation
1. **TACSTD2 protein is not associated with xCell CD8 or immune** in this
   treatment-naive LSCC proteome. Residualizing on purity does not create a
   signal that was not there.
2. **CLDN4 protein is the immune-cold partner at the protein layer**, matching
   PR23. Part of that inverse is purity (CLDN4 vs WES ρ=+0.305). After WES
   residual a smaller inverse remains (FDR=0.023). After WGS residual it is
   gone. The purity-independent claim is **covariate-sensitive**.
3. This is **not** an ICI cohort. Cross-layer (protein vs RNA-xCell)
   observational association only.

### Caveats
- No ICI labels; OS/PFS were not retested here.
- CLDN4 protein missing 30/108; those tumors are dropped.
- xCell is RNA deconvolution; protein–xCell mixes layers.
- WES vs xCell immune ρ=−0.72, so residualizing immune on WES removes a large
  shared component by construction. That is the test, not a bug.
- WGS is missing 4 tumors and 3 of the CLDN4-complete cases; the WGS-null
  result has less n.
- Table-wise BH on 4 primary tests only.
- LSCC ≠ LUAD.

### How to rerun
```bash
pip install -r scripts/w200/CPTAC_LSCC_xcell/requirements.txt
python3 scripts/w200/CPTAC_LSCC_xcell/download.py
python3 scripts/w200/CPTAC_LSCC_xcell/analyze.py
```
Matrices: `data/w200/CPTAC_LSCC_xcell/` (not committed).
Tables/figures: `results/w200/CPTAC_LSCC_xcell/`.

---

## 中文

### 问题
在 CPTAC **治疗初治**肺鳞癌（LSCC / LUSC）中，**TACSTD2** 与 **CLDN4 蛋白**
是否与冻存 **xCell CD8**、**xCell immune score** 相关？在用肿瘤纯度做残差
（partial Spearman）之后是否还在？本队列**没有 ICI 疗效标签**。

### 数据与方法
开放 S3 冻存 `data_freeze_v1.2_reorganized`。肿瘤蛋白 108 例。xCell / 纯度在
`LSCC_phenotype.txt`。TACSTD2 蛋白 108/108；**CLDN4 蛋白仅 78/108（NA=30）**。
预设终点：两蛋白 × xCell CD8 / immune。主残差协变量为 **WES_purity**（DNA，
107/108）；WGS 与 ESTIMATE 余弦纯度为敏感性分析。ESTIMATE 纯度由免疫/基质
RNA 变换而来，与 xCell immune **循环**，不能当作诚实残差。边缘 Spearman
与 PR23 一致。

### 主要结果（未编造）
- **TACSTD2 蛋白与 xCell CD8 / immune 均不显著**（边缘 ρ=−0.080 / −0.133；
  WES 偏相关 ρ=−0.061 / −0.129，FDR=0.53 / 0.25）。残差**不会**变出信号。
  **不能**写成 Bessede 式“高 TROP2 = T 细胞排斥”。
- **CLDN4 蛋白**边缘与 CD8 ρ=−0.387、immune ρ=−0.400（与 PR23 相同）。
  WES 残差后仍为负且 FDR=0.023，但变弱（ρ=−0.303 / −0.288）。
  **WGS 残差后不显著**（ρ=−0.173 / −0.049，p=0.14 / 0.68，n=75）。
  “独立于纯度的 CLDN4–免疫关联”**随协变量而变**，不能当成定论。
- 背景：CLDN4 蛋白与 WES 纯度 ρ=+0.305（p=0.007）；TACSTD2 与纯度无关
  （ρ=+0.066，p=0.50）；xCell immune 与 WES ρ=−0.724。

### 局限
无 ICI 标签；CLDN4 缺失 30 例；xCell 来自 RNA；WES 与 immune 高度负相关，
残差会去掉大量共用成分；WGS 例数更少；多重检验仅为 4 项主检验的 BH；
不可与 LUAD 未经预设交互而合并。

### 复现
见上文 English “How to rerun”。
