# NSCLC ICI exclusion scores playbook / NSCLC ICI 排斥评分操作手册

## English

### 1. What is implemented

| Output | Interpretation | Status |
|---|---|---|
| TIDE, Dysfunction, Exclusion, MDSC, CAF, TAM M2 | Author-defined immune-evasion model | Canonical through `run_tide.py` + authors' TIDEpy v1.3 |
| IPS, MHC, CP, EC, SC | Immunophenoscore (0–10) and components | Canonical public algorithm |
| IFNG6_mean, Higgs_IFNG4_mean, CTL5_mean | Inflamed/CTL marker means | Canonical public gene definitions; Higgs-4 has direct NSCLC durvalumab evidence |
| TIS18_public_mean | Unweighted public 18-gene summary | Proxy; **not** clinical TIS |
| CAF_collagen3_mean | Bulk fibroblast/collagen signal | Portable proxy; **not** TIDE CAF |
| HALLMARK_EMT_mean | Mean of the 200-gene Hallmark EMT set | Transparent alternative to ssGSEA |
| Thompson EMT/inflammation | Signed NSCLC EMT/inflammation model | Exact public formula; cohort-dependent z scores |
| hMENA–TGF-β CAF9, Pan-F-TBRS | Stromal exclusion gene sets | Exact ssGSEA bridge plus clearly named mean proxies |
| Bindea_cytotoxic_mean | Cytotoxic marker QC | Proxy; use ConsensusTME for Bessede replication |

All local scoring is dependency-free Python. See `SOURCES.md` before comparing
scores across studies.

### 2. Input contract and QC

Use a tab-separated genes-by-samples matrix:

```text
gene	S01	S02	S03
TACSTD2	8.2	4.1	6.0
CD8A	2.4	7.8	5.1
```

- First column: unique, current HGNC symbols. The scripts uppercase symbols.
- Remaining columns: finite numeric values from pretreatment tumor samples.
- For local means and IPS: use one consistently processed, log-scale matrix,
  preferably `log2(TPM + 1)` after sample/gene QC. Do not mix platforms.
- For TIDE: use whole transcriptome, follow TIDE's study-centering/reference
  instructions, and preserve the exact matrix submitted. TIDEpy may transform
  input internally; do not silently transform twice.
- Predefine duplicate-gene handling upstream. The scripts reject duplicates.
- Report per-score coverage. Do not interpret a score with `<80%` of genes;
  IPS additionally needs all four classes represented.
- Scores are relative research biomarkers, not clinical decision rules.

### 3. Run transparent signatures and IPS

From the repository root:

```bash
python3 methods/exclusion_scores/scripts/score_signatures.py \
  expression.tsv scores.tsv
```

Select only prespecified signatures by repeating `--signature`:

```bash
python3 methods/exclusion_scores/scripts/score_signatures.py \
  expression.tsv scores.tsv \
  --signature IFNG6_mean \
  --signature CAF_collagen3_mean \
  --signature HALLMARK_EMT_mean
```

The mean signatures preserve the input scale. IPS first z-scores all measured
genes within each sample, then reproduces the public factor/class aggregation.
Consequently, IPS requires a broad, comparable transcriptome—not a targeted
file containing only IPS genes.

Run the directly NSCLC-derived Thompson model separately because it requires
gene-wise z-scoring across the analyzed cohort:

```bash
python3 methods/exclusion_scores/scripts/score_thompson_emt.py \
  expression.tsv thompson_scores.tsv
```

It outputs signed 12-gene EMT, 27-gene inflammation, unweighted
`inflammation−EMT`, and the published fitted combination
`−0.60×EMT+0.19×inflammation`. The reported performance came from a small
retrospective cohort in which the model was internally fitted and evaluated;
do not reuse its reported cutoff as externally validated.

For exact hMENA–TGF-β CAF9 and Pan-F-TBRS ssGSEA:

```bash
Rscript methods/exclusion_scores/scripts/score_ssgsea.R \
  expression.tsv methods/exclusion_scores/gene_sets/signatures.tsv ssgsea.tsv
```

The output also records whether `TGFB1` exceeds the cohort 10th percentile,
the eligibility gate reported for the 2026 hMENA study. Because ssGSEA is
cohort- and implementation-sensitive, record the GSVA/Bioconductor version.

### 4. Run canonical TIDE

Install the authors' TIDEpy v1.3 in a pinned, isolated environment. Its model
objects and historical dependencies are the method; this repository does not
reimplement or silently approximate them.

```bash
python3 methods/exclusion_scores/scripts/run_tide.py \
  centered_whole_transcriptome.tsv tide.tsv
```

The bridge fixes cancer type to `NSCLC`, checks for whole-transcriptome input,
and verifies the expected result columns. Add `--pretreat` only for prior
immunotherapy, not prior chemotherapy or targeted treatment.

Record the TIDEpy commit/version, Python environment, cancer choice, pretreat
flag, input transformation, reference/centering method, and missing genes.
TIDEpy's own documentation says validation in NSCLC was limited and specifically
warns about unreferenced RNA-seq values.

### 5. Faithful Bessede-style TACSTD2 analysis

Bessede et al. analyzed POPLAR/OAK pretreatment tumors. Their immune comparison
was **TACSTD2-high versus TACSTD2-low ConsensusTME/Bindea cell estimates**; it
was not a TACSTD2-versus-TIDE-Exclusion test.

Faithful immune-estimate workflow in R:

```r
# expression: genes x samples, same normalized baseline matrix
bindea <- ConsensusTME::methodSignatures$Bindea
immune <- ConsensusTME::geneSetEnrichment(expression, bindea)
# Compare prespecified T, Tfh, cytotoxic and B-cell estimates between groups.
# Use two-sided Wilcoxon tests and BH-adjust across tested cell populations.
```

For a new cohort, avoid choosing a TACSTD2 cut point on the same outcome.
Primary analysis should keep TACSTD2 continuous (standardized); use a median
split only as a display/sensitivity analysis. The supplied survival script
expects:

```text
pfs_time pfs_event os_time os_event treatment TACSTD2 histology CD274 TLS
```

and runs adjusted continuous Cox models, a treatment-by-TACSTD2 interaction,
and a median-split sensitivity model:

```bash
Rscript methods/exclusion_scores/scripts/bessede_survival.R \
  clinical.tsv survival_results.tsv
```

A biomarker is treatment-predictive only if the interaction is supported (for
example, TACSTD2 has a different association under atezolizumab and docetaxel).
Significance in one arm and non-significance in another does not itself prove an
interaction. Check proportional hazards, nonlinearity, influential cases,
events per parameter, and cohort/trial stratification before inference.

### 6. TACSTD2 versus exclusion: prespecified extension

This answers a distinct mechanistic question: does higher epithelial TACSTD2
track a more excluded/stromal microenvironment?

There are two distinct “exclusion” extensions:

- **Transcriptomic:** canonical TIDE `Exclusion`, as described below.
- **Spatial:** a later Bessede-coauthored SITC poster defined a CD8 compartment
  IES and associated it with DDR1, not TACSTD2. No primary source reports a
  direct TACSTD2–IES test, and the exact IES formula/signature was not disclosed.

For exploratory spatial data, `score_spatial_ies.py` implements the declared
reconstruction
`[log2(stromal+epsilon)−log2(epithelial+epsilon)]/sqrt(2)`, where positive
values indicate stromal enrichment:

```bash
python3 methods/exclusion_scores/scripts/score_spatial_ies.py \
  compartment_cd8.tsv spatial_ies.tsv
```

Input columns are `sample`, `epithelial_cd8_density`, and
`stromal_cd8_density`. The default pseudocount is half the smallest positive
density and is written to output. Preserve density units and segmentation
rules. The output is labeled `reconstruction_not_author_exact`; never present
it as the poster's exact IES. Test continuous TACSTD2 against this IES by
Spearman correlation, then model IES with histology, PD-L1, tumor purity, batch,
and total CD8 density to distinguish redistribution from an immune desert.

1. Primary score: canonical TIDE `Exclusion`.
2. Sensitivity scores: TIDE `CAF`, `CAF_collagen3_mean`, `HALLMARK_EMT_mean`.
3. Positive immune-accessibility controls: `CTL5_mean`,
   `Bindea_cytotoxic_mean` (expected opposite direction).
4. Primary test: continuous Spearman correlation. Report rho, confidence
   interval (bootstrap in the final statistical workflow), exact/permutation
   two-sided P, n and a scatterplot. Adjust the prespecified score family by BH.
5. Sensitivity: median TACSTD2 groups; avoid outcome-optimized cut points.
6. Adjusted model: regress the score on continuous TACSTD2 plus histology,
   trial/batch, tumor purity, smoking, driver status and PD-L1/CD274 when
   available. Do not adjust for a mediator such as CAF unless estimating a
   direct effect.
7. Replicate in an independent ICI-treated NSCLC cohort. Bulk EMT/CAF can
   reflect stromal abundance rather than tumor-cell EMT; use pathology/purity
   or spatial data to resolve this.

Run the dependency-free screen:

```bash
python3 methods/exclusion_scores/scripts/test_tacstd2_exclusion.py \
  expression.tsv tide.tsv tacstd2_exclusion.json \
  --score Exclusion --score CAF --permutations 9999
```

For local scores, point the second argument to `scores.tsv`. The script reports
Spearman permutation tests and median-split effect sizes; it deliberately does
not claim covariate adjustment or survival inference.

### 7. Reproducibility checklist

- Lock cohort, endpoint, treatment line, biopsy timing and exclusions before
  viewing outcomes.
- Separate predictive analyses (treatment interaction) from prognostic ones.
- Preserve raw-to-score code, gene annotation, normalization and score coverage.
- Freeze score versions; never call `TIS18_public_mean` “TIS” or a local CAF
  mean “TIDE CAF.”
- Validate direction and scale on synthetic controls, then on a known dataset.
- Report all tested scores and multiplicity correction, including null results.
- Follow REMARK/TRIPOD as applicable; no score here is a companion diagnostic.

## 中文

### 1. 已实现内容

| 输出 | 含义 | 状态 |
|---|---|---|
| TIDE、Dysfunction、Exclusion、MDSC、CAF、TAM M2 | 作者定义的免疫逃逸模型 | 通过 `run_tide.py` 调用作者 TIDEpy v1.3，属于规范实现 |
| IPS、MHC、CP、EC、SC | 免疫表型评分及分量 | 公开算法的规范实现 |
| IFNG6_mean、Higgs_IFNG4_mean、CTL5_mean | 炎症/细胞毒性标志基因均值 | 公开定义；Higgs-4 有直接 NSCLC durvalumab 证据 |
| TIS18_public_mean | 公开 18 基因的无权重均值 | 代理指标，**不是**临床 TIS |
| CAF_collagen3_mean | 成纤维/胶原信号 | 可迁移代理指标，**不是** TIDE CAF |
| HALLMARK_EMT_mean | Hallmark EMT 200 基因均值 | 透明实现，不等同于 ssGSEA |
| Thompson EMT/炎症模型 | NSCLC 有符号 EMT/炎症评分 | 精确公开公式；依赖队列内 z 标准化 |
| hMENA–TGF-β CAF9、Pan-F-TBRS | 基质排斥签名 | ssGSEA 规范接口及明确标注的均值代理 |
| Bindea_cytotoxic_mean | 细胞毒性快速质控 | 代理指标；复现 Bessede 应使用 ConsensusTME |

本地评分脚本仅依赖 Python 标准库。跨队列比较前必须阅读
`SOURCES.md`。

### 2. 输入与质控

输入为制表符分隔的“基因 × 样本”矩阵，第一列为唯一的 HGNC 基因符号，
其余列为治疗前肿瘤样本的有限数值。推荐本地评分使用统一流程生成的
`log2(TPM + 1)`；不得混合平台。重复基因应在上游按预先规定的规则处理，
脚本会拒绝重复行。

TIDE 必须输入全转录组，并遵循其按研究队列中心化/参考样本的要求；保存
实际提交矩阵，避免重复 log 转换。每个评分都要报告覆盖率，低于 80%
不作解释；IPS 还要求四个类别均有基因。所有输出都是研究指标，不是临床
用药阈值。

### 3. 运行方法

本地签名与 IPS：

```bash
python3 methods/exclusion_scores/scripts/score_signatures.py \
  expression.tsv scores.tsv
```

规范 TIDE：

```bash
python3 methods/exclusion_scores/scripts/run_tide.py \
  centered_whole_transcriptome.tsv tide.tsv
```

应在隔离环境中固定作者 TIDEpy v1.3 及其依赖，并记录版本、癌种
`NSCLC`、既往免疫治疗参数、输入转换、中心化参考及缺失基因。

直接来源于 NSCLC 的 Thompson EMT/炎症模型必须单独运行，因为它按本队列
对每个基因做 z 标准化：

```bash
python3 methods/exclusion_scores/scripts/score_thompson_emt.py \
  expression.tsv thompson_scores.tsv
```

输出包括 12 基因有符号 EMT、27 基因炎症、`炎症−EMT` 以及论文组合
`−0.60×EMT+0.19×炎症`。原研究在同一小型回顾性队列拟合并评价，不能把
其阈值当作外部验证阈值。

hMENA–TGF-β CAF9 和 Pan-F-TBRS 的 ssGSEA：

```bash
Rscript methods/exclusion_scores/scripts/score_ssgsea.R \
  expression.tsv methods/exclusion_scores/gene_sets/signatures.tsv ssgsea.tsv
```

输出同时标记 `TGFB1` 是否高于队列第 10 百分位（2026 hMENA 研究的门控
条件）。应记录 GSVA/Bioconductor 版本。

### 4. Bessede 式 TACSTD2 分析

Bessede 等人的 2024 年研究使用 POPLAR/OAK 治疗前 RNA-seq。论文比较的
是 TACSTD2 高低组之间的 ConsensusTME/Bindea 细胞估计，并没有计算
TIDE Exclusion。严格复现应使用：

```r
bindea <- ConsensusTME::methodSignatures$Bindea
immune <- ConsensusTME::geneSetEnrichment(expression, bindea)
```

对预先指定的 T、Tfh、细胞毒性和 B 细胞评分做双侧 Wilcoxon 检验，并对
多个细胞类型做 BH 校正。新队列中，主要分析应使用连续标准化 TACSTD2；
中位数分组仅用于展示/敏感性分析，不应在同一结局上寻找“最佳”阈值。

生存及治疗交互模型：

```bash
Rscript methods/exclusion_scores/scripts/bessede_survival.R \
  clinical.tsv survival_results.tsv
```

输入列为 `pfs_time, pfs_event, os_time, os_event, treatment, TACSTD2,
histology, CD274, TLS`。只有治疗 × TACSTD2 交互得到支持时，才能称其为
治疗预测标志物；“一组显著、另一组不显著”并不等于交互显著。正式报告前
还需检查比例风险、非线性、异常点、事件数以及试验/队列分层。

### 5. TACSTD2 与排斥的扩展检验

该问题与论文复现不同：它检验 TACSTD2 是否伴随更强的免疫排斥/基质状态。

“排斥”还需区分转录组 TIDE Exclusion 与空间 CD8 IES。后者来自 Bessede
共同署名的 SITC 2024 poster，研究对象是 DDR1，并未直接检验 TACSTD2；
作者也未公开完整公式或 RNA 签名。探索性空间数据可运行：

```bash
python3 methods/exclusion_scores/scripts/score_spatial_ies.py \
  compartment_cd8.tsv spatial_ies.tsv
```

脚本明确采用重建公式
`[log2(基质密度+epsilon)−log2(上皮密度+epsilon)]/sqrt(2)`，正值表示 CD8
偏留于基质。输入列为 `sample、epithelial_cd8_density、
stromal_cd8_density`；输出会标记 `reconstruction_not_author_exact`，不得
写成作者原始 IES。调整模型应加入组织学、PD-L1、纯度、批次和总 CD8
密度，以区分真正的空间重分布与免疫荒漠。

1. 主要指标预先指定为规范 TIDE `Exclusion`。
2. 敏感性指标为 TIDE `CAF`、`CAF_collagen3_mean` 和
   `HALLMARK_EMT_mean`。
3. `CTL5_mean` 与 `Bindea_cytotoxic_mean` 作为方向相反的免疫可及性
   对照。
4. 主要统计量为连续值 Spearman rho；报告样本数、双侧置换 P 值、散点图，
   正式分析另加 bootstrap 置信区间，并对评分家族做 BH 校正。
5. 中位数分组只作敏感性分析，不使用结局优化阈值。
6. 调整模型可加入组织学、试验/批次、肿瘤纯度、吸烟、驱动基因及
   PD-L1/CD274。若目标是总效应，不要把可能的中介变量 CAF 纳入调整。
7. 在独立 ICI 治疗 NSCLC 队列复现。bulk EMT/CAF 可能反映基质含量，
   不一定是肿瘤细胞自身 EMT，应结合病理纯度或空间数据。

快速筛查命令：

```bash
python3 methods/exclusion_scores/scripts/test_tacstd2_exclusion.py \
  expression.tsv tide.tsv tacstd2_exclusion.json \
  --score Exclusion --score CAF --permutations 9999
```

该脚本只报告相关性置换检验和中位数组效应，不声称完成协变量调整、生存
推断或因果推断。

### 6. 最低复现要求

预先锁定队列、终点、治疗线次和活检时间；明确区分治疗预测与一般预后；
保存从原始数据到评分的代码、注释版本、归一化和覆盖率；冻结评分版本；
不得把代理 TIS/CAF 写成规范 TIS/TIDE CAF；报告全部检验和多重校正结果，
包括阴性结果。以上评分均不是伴随诊断。
