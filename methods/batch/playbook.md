# Multi-GEO TACSTD2 meta-analysis: batch-aware methods playbook

# 多 GEO TACSTD2 荟萃分析：批次感知方法手册

## 1. Scope and estimand / 范围与目标效应

**EN.** This protocol tests whether **TACSTD2** expression differs between biological response groups across independent Gene Expression Omnibus (GEO) studies. TACSTD2 is also commonly annotated as **EPCAM**; retain the original feature annotation and map both symbols to the current stable gene identifier (for example, the appropriate Ensembl gene ID for the organism) before analysis.

**中文。** 本方案用于检验多个独立 GEO 研究中，生物学应答组之间的 **TACSTD2** 表达是否存在差异。TACSTD2 也常被注释为 **EPCAM**；分析前应保留原始探针/基因注释，并将两个符号映射到对应物种的当前稳定基因标识符（例如 Ensembl gene ID）。

Define, before inspecting results:

- the biological contrast, such as responder versus non-responder or post-treatment versus baseline;
- the experimental unit (usually the patient, not a technical replicate);
- the primary endpoint and direction, defined so positive effects always have the same meaning;
- the eligible tissues, disease settings, organisms, and assay types;
- one primary TACSTD2 summary per independent cohort.

在查看结果之前预先定义：

- 生物学对比，例如应答者与非应答者，或治疗后与基线；
- 实验单位（通常是患者，而不是技术重复）；
- 主要终点及效应方向，使正效应始终代表同一含义；
- 可纳入的组织、疾病背景、物种和检测平台；
- 每个独立队列仅产生一个主要 TACSTD2 汇总效应。

The primary estimand is the **within-study adjusted response effect**, followed by a random-effects meta-analysis. It is not the coefficient obtained by concatenating all expression matrices and treating GEO accession as an ordinary batch.

主要目标效应是**研究内经协变量调整后的应答效应**，随后进行随机效应荟萃分析。不能把所有表达矩阵直接拼接，再把 GEO 编号当作普通批次来拟合一个系数。

## 2. Study inventory and eligibility / 研究清单与纳入标准

Create a frozen study manifest with one row per sample and these fields:

| Field / 字段 | Required content / 必需内容 |
|---|---|
| `study_id` | GEO series or biologically independent cohort / GEO 系列或生物学独立队列 |
| `sample_id` | Stable sample identifier / 稳定样本标识 |
| `subject_id` | Patient/animal identifier for paired or repeated data / 配对或重复数据的个体标识 |
| `response` | Harmonized biological response; preserve original label / 统一后的生物学应答，同时保留原标签 |
| `batch` | Library preparation, sequencing run, plate, processing date, or documented technical batch / 建库、测序批、芯片板、处理日期或其他已记录技术批次 |
| `platform` | GEO platform and assay type / GEO 平台及检测类型 |
| `tissue` | Tissue, compartment, and cell-selection method / 组织、区室及细胞筛选方法 |
| `timepoint` | Sampling time relative to intervention / 相对干预的采样时间 |
| `covariates` | Prespecified age, sex, stage, purity, treatment, etc. / 预先指定的年龄、性别、分期、纯度、治疗等 |
| `raw_available` | Raw counts/CEL/intensity availability / 原始计数、CEL 或强度数据可用性 |

**EN.** Include a cohort only when response labels, sample identity, and assay provenance can be reconciled. Detect overlapping cohorts by accession cross-references, sample titles, subject IDs, and publication methods. If two GEO accessions reuse subjects, retain one analysis or model the overlap; do not count them as independent studies.

**中文。** 仅当应答标签、样本身份和检测来源能够核对一致时纳入队列。应通过 GEO 编号交叉引用、样本标题、个体 ID 和论文方法识别重复队列。若两个 GEO 数据集重复使用同一批受试者，应只保留一个分析，或显式建模其相关性；不得把它们当作独立研究重复计数。

Record exclusions and reasons before effect estimation. Never exclude a study because its TACSTD2 result is null or opposite in direction.

在估计效应前记录排除项及理由。不得因某研究的 TACSTD2 结果不显著或方向相反而排除该研究。

## 3. Data preparation by assay / 按检测类型预处理

### 3.1 RNA-seq counts / RNA-seq 计数

**EN.**

1. Start from an integer gene-by-sample count matrix whenever possible.
2. Use one consistent gene annotation release within each study and document identifier mapping.
3. Remove failed libraries using prespecified sequencing QC, not TACSTD2 expression.
4. Filter low-expression genes without reference to response direction (for example, a CPM threshold in a minimum number of samples).
5. Estimate library-size normalization factors within each study (for example, TMM or median-ratio normalization).
6. Fit the response contrast to counts with an appropriate negative-binomial model, or use voom with precision weights. Include known batch variables in the design when estimable.

**中文。**

1. 尽可能从整数型“基因 × 样本”计数矩阵开始。
2. 每项研究内部使用一致的基因注释版本，并记录标识符映射。
3. 按预先规定的测序质控标准剔除失败文库，不得依据 TACSTD2 表达剔除样本。
4. 低表达过滤不得参考应答效应方向（例如要求至少若干样本达到指定 CPM）。
5. 在每项研究内部估计文库大小归一化因子（例如 TMM 或 median-ratio）。
6. 使用合适的负二项模型拟合应答对比，或使用带精度权重的 voom。若批次效应可识别，则在设计矩阵中纳入已知批次变量。

### 3.2 Microarray / 微阵列

**EN.**

1. Prefer raw files and process all samples from the same study together with the platform-appropriate background correction and normalization (for example, RMA for compatible Affymetrix arrays).
2. Apply platform-specific probe annotation. Remove probes with known cross-hybridization or ambiguous multi-gene mapping when such annotations are available.
3. Predefine TACSTD2 probe summarization. Prefer a uniquely mapping, reliably expressed probe; if several valid probes remain, use a documented gene-level summary or analyze probes and combine them without selecting the most significant one.
4. Fit a linear model within each study, including known technical and biological covariates.

**中文。**

1. 优先使用原始文件；同一研究的所有样本应一起采用平台适用的背景校正和归一化方法处理（例如对兼容的 Affymetrix 芯片使用 RMA）。
2. 使用平台特异的探针注释。若有可靠注释，应去除已知交叉杂交或映射到多个基因的探针。
3. 预先定义 TACSTD2 探针汇总规则。优先选择唯一映射且稳定表达的探针；若存在多个合格探针，应使用已记录的基因层面汇总方法，或分别分析后合并，不得挑选最显著的探针。
4. 在每项研究内部拟合线性模型，并纳入已知技术及生物学协变量。

Do not apply ComBat-seq to microarray intensities, TPM, FPKM, CPM, log-counts, variance-stabilized values, or any non-integer matrix.

不得将 ComBat-seq 用于微阵列强度、TPM、FPKM、CPM、对数计数、方差稳定化数据或任何非整数矩阵。

## 4. Decide whether batch is identifiable / 判断批次效应是否可识别

For every study, tabulate response by each candidate batch and inspect the design matrix:

```text
table(response, batch)
rank(model.matrix(~ response + batch + prespecified_covariates))
```

对每项研究，交叉列出应答与每个候选批次，并检查设计矩阵：

```text
table(response, batch)
rank(model.matrix(~ response + batch + prespecified_covariates))
```

Classify the design:

1. **Crossed or partially crossed / 完全或部分交叉：** both response groups occur in enough batches to estimate response and batch separately. Include batch in the model.
2. **Sparse but estimable / 稀疏但可估计：** some cells are empty, but the design matrix retains full rank. Use a parsimonious model and report sensitivity analyses.
3. **Fully confounded / 完全混杂：** all responders are in one batch and all non-responders in another, or an equivalent linear dependency exists. The response effect is not identifiable from these data.

分类如下：

1. **完全或部分交叉：** 足够多的批次同时包含两个应答组，可分别估计应答与批次。模型中应纳入批次。
2. **稀疏但可估计：** 某些组合为空，但设计矩阵仍满秩。采用简约模型并报告敏感性分析。
3. **完全混杂：** 所有应答者均处于一个批次、所有非应答者均处于另一个批次，或存在等价的线性依赖。此时数据无法识别应答效应。

**EN.** No batch-correction method—ComBat, ComBat-seq, SVA, RUV, Harmony, or otherwise—can reconstruct an unobserved response-versus-batch contrast. For a fully confounded study, obtain additional crossed samples, restrict inference to an estimable subset, or exclude that contrast with the reason “response–batch confounding.”

**中文。** 任何批次校正方法——包括 ComBat、ComBat-seq、SVA、RUV、Harmony 等——都无法重建未被观测到的“应答与批次”对比。对于完全混杂的研究，应补充跨批次样本、将推断限制在可估计子集，或以“应答–批次完全混杂”为由排除该对比。

## 5. Known batch: model it first / 已知批次：优先在模型中调整

The default inferential model is:

```text
expression ~ response + known_batch + prespecified_biological_covariates
```

默认推断模型为：

```text
expression ~ response + known_batch + prespecified_biological_covariates
```

For paired or repeated samples, include subject blocking or a suitable random/correlation structure. Do not treat repeated samples as independent.

对于配对或重复测量，应加入个体阻断项或合适的随机/相关结构，不得将重复样本视为相互独立。

**EN.** Direct covariate adjustment preserves uncertainty in the response coefficient and avoids making corrected values look more certain than they are. Batch-corrected matrices may be produced for PCA, heat maps, or clustering, but primary hypothesis tests should come from the fitted model on appropriately normalized, uncorrected observations.

**中文。** 直接在模型中调整协变量，可以保留应答系数的不确定性，也可避免“校正后数值看起来比实际更确定”。批次校正矩阵可用于 PCA、热图或聚类，但主要假设检验应来自对适当归一化、未经批次值替换的数据所拟合的模型。

## 6. ComBat-seq / ComBat-seq

Use ComBat-seq only when all of the following hold:

- the input is raw, non-negative integer RNA-seq counts;
- batches are known;
- response and batch are not fully confounded;
- each modeled batch has adequate replication;
- the intended downstream method can validly consume the adjusted integer counts.

仅在以下条件全部满足时使用 ComBat-seq：

- 输入为原始、非负整数 RNA-seq 计数；
- 批次标签已知；
- 应答与批次不存在完全混杂；
- 每个建模批次具有足够重复；
- 下游方法可以合理使用校正后的整数计数。

Specify the biological response and other protected biological covariates through the group/covariate design supported by the implementation. Never run an intercept-only batch correction when the response is the target of inference.

应通过实现所支持的分组/协变量设计明确指定需要保护的生物学应答及其他生物学协变量。当应答是推断目标时，绝不能进行仅含截距的批次校正。

**Recommended role / 推荐用途：**

- sensitivity analysis of clustering or visualization;
- a prespecified secondary analysis when a downstream count-based workflow requires one matrix;
- not the default route for the primary TACSTD2 response coefficient when batch can be included directly in the study-level model.

- 用于聚类或可视化的敏感性分析；
- 当下游计数流程必须使用一个矩阵时，作为预先指定的次要分析；
- 若可在研究层面模型中直接纳入批次，则不作为 TACSTD2 主要应答系数的默认分析路径。

After ComBat-seq, repeat library and sample QC, verify that protected response separation has not collapsed, and compare the TACSTD2 estimate with the model-adjusted primary analysis. A large discrepancy is a diagnostic signal, not a reason to select the preferred result.

使用 ComBat-seq 后，应重复文库及样本质控，确认受保护的应答差异未消失，并将 TACSTD2 效应与主要模型调整分析比较。明显不一致属于诊断信号，不能据此挑选更“理想”的结果。

## 7. SVA / 替代变量分析

Use surrogate variable analysis (SVA) for unmeasured, sample-level heterogeneity after accounting for known covariates.

在调整已知协变量后，可使用替代变量分析（SVA）处理未测量的样本层面异质性。

Define:

```r
mod  <- model.matrix(~ response + known_batch + covariates, metadata)
mod0 <- model.matrix(~ known_batch + covariates, metadata)
```

Then estimate the number of surrogate variables and the variables themselves using an approach appropriate to the data scale. For count data, use a count-aware SVA workflow or an appropriate normalized transformation as required by the chosen implementation. Refit the expression model with the estimated surrogate variables:

```text
expression ~ response + known_batch + covariates + SV1 + ... + SVk
```

其中：

```r
mod  <- model.matrix(~ response + known_batch + covariates, metadata)
mod0 <- model.matrix(~ known_batch + covariates, metadata)
```

随后根据数据尺度选择合适方法估计替代变量的数量及其取值。对计数数据，应采用计数感知的 SVA 流程，或按所选实现要求使用恰当的归一化变换。最后在表达模型中加入替代变量：

```text
expression ~ response + known_batch + covariates + SV1 + ... + SVk
```

The full model (`mod`) tells supervised SVA which biological signal must be protected; the null model (`mod0`) omits only the tested response term. Confirm the exact argument convention of the selected SVA function.

完整模型 `mod` 用于告诉监督式 SVA 哪些生物学信号必须保留；零模型 `mod0` 仅去除待检验的应答项。应核对所用 SVA 函数的具体参数约定。

Do not:

- estimate surrogate variables from TACSTD2 alone;
- choose the number of variables to maximize TACSTD2 significance;
- include so many variables that the response coefficient becomes unstable;
- interpret an SVA-adjusted association as causal control of an identified confounder.

不得：

- 仅依据 TACSTD2 估计替代变量；
- 以最大化 TACSTD2 显著性为目标选择替代变量个数；
- 加入过多替代变量而导致应答系数不稳定；
- 将 SVA 调整后的相关性解释为对某个已识别混杂因素的因果控制。

## 8. RUV / 非期望变异去除

RUV methods require credible negative-control genes or replicated samples:

- **RUVg:** uses negative-control genes;
- **RUVs:** uses replicate samples for which the biological factor of interest is constant;
- **RUVr:** uses residuals from a first-pass model and therefore needs extra care to avoid absorbing the target signal.

RUV 方法依赖可信的负对照基因或重复样本：

- **RUVg：** 使用负对照基因；
- **RUVs：** 使用目标生物学因素保持不变的重复样本；
- **RUVr：** 使用初始模型残差，因此尤其需要防止吸收目标信号。

Negative controls must be selected independently of the observed TACSTD2 response result. Suitable controls may come from validated spike-ins or externally justified stable genes. “Non-significant in this dataset” is not, by itself, proof that a gene is unaffected by response.

负对照必须独立于当前数据中观察到的 TACSTD2 应答结果进行选择。可使用经验证的外源 spike-in 或有外部证据支持的稳定基因。仅仅“在本数据中不显著”并不能证明某基因不受应答影响。

Select the unwanted-factor dimension `k` using prespecified diagnostics such as replicate agreement, relative log-expression behavior, PCA association with known technical factors, and stability of the genome-wide effect distribution. Report results over a small plausible range of `k`; do not optimize `k` for the TACSTD2 P value.

应依据预先规定的诊断标准选择非期望因子维度 `k`，例如重复样本一致性、相对对数表达行为、PCA 与已知技术因素的关联，以及全基因组效应分布的稳定性。应报告一个较小合理 `k` 范围内的结果，不得针对 TACSTD2 的 P 值优化 `k`。

Include estimated unwanted factors in the study-level model alongside the protected response term. If controls may respond biologically, prefer a different control set or treat RUV only as a sensitivity analysis.

将估计出的非期望因子与受保护的应答项一同纳入研究层面模型。若对照基因可能发生生物学应答，应更换对照集，或仅将 RUV 作为敏感性分析。

## 9. Why response must not be batch-corrected away / 为什么绝不能把应答信号当作批次消除

Batch correction estimates unwanted variation from patterns in the expression matrix. If biological response is correlated with batch and is not explicitly protected, the algorithm cannot know which part is technical. It may subtract true response signal, reverse its direction, compress its variance, or create an apparently clean but biased embedding.

批次校正根据表达矩阵中的模式估计非期望变异。如果生物学应答与批次相关，且未被明确保护，算法无法判断哪部分是技术效应。它可能减去真实应答信号、反转效应方向、压缩其方差，或生成看似整洁但有偏的低维图。

Use these non-negotiable safeguards:

1. Draw the causal/design relationships among response, treatment, tissue, timepoint, study, platform, and batch before correction.
2. Put the target response in every protected/full design used to estimate unwanted variation.
3. Never define “batch” from clusters that are themselves separated by response.
4. Never select a correction method because it yields the smallest TACSTD2 P value or the prettiest PCA.
5. Compare effect direction and magnitude before and after nuisance adjustment.
6. Plot known technical and biological labels on diagnostic embeddings.
7. Declare non-identifiability rather than forcing correction when response and batch are fully confounded.

必须遵守以下原则：

1. 校正前绘制应答、治疗、组织、时间点、研究、平台和批次之间的因果/设计关系。
2. 在所有用于估计非期望变异的保护/完整设计中加入目标应答项。
3. 不得根据本身由应答分离形成的聚类来定义“批次”。
4. 不得因为某方法产生最小的 TACSTD2 P 值或最“漂亮”的 PCA 而选择它。
5. 比较非期望变异调整前后的效应方向及大小。
6. 在诊断性降维图上同时标注已知技术因素与生物学标签。
7. 当应答与批次完全混杂时，应明确声明不可识别，而不是强行校正。

Loss of response separation in PCA is not evidence of successful correction. Conversely, preserved response separation is not proof of validity. The inferential checks are design identifiability, coefficient stability, calibrated uncertainty, and replication across independent studies.

PCA 中应答组不再分离并不代表校正成功；应答分离得到保留也不能证明方法有效。推断有效性应依据设计可识别性、系数稳定性、不确定性校准，以及独立研究间的重复验证。

## 10. Study-level TACSTD2 effect / 研究层面的 TACSTD2 效应

Estimate one prespecified effect and standard error per independent cohort.

每个独立队列估计一个预先指定的效应及其标准误。

Preferred effects:

- **Comparable log-expression scales:** adjusted mean difference or log2 fold change with its standard error.
- **Incompatible platforms/scales:** Hedges' \(g\), computed from the response contrast with small-sample correction.
- **Paired designs:** within-subject effect with a paired standard error.
- **Time-to-event or binary clinical outcomes:** use an outcome-appropriate regression coefficient; do not mix these coefficients with continuous-expression contrasts in one meta-analysis.

推荐效应量：

- **可比的对数表达尺度：** 经调整的均值差或 log2 倍数变化及其标准误；
- **平台/尺度不兼容：** 使用带小样本校正的 Hedges' \(g\)；
- **配对设计：** 个体内效应及配对标准误；
- **生存或二分类临床结局：** 使用与结局相匹配的回归系数；不得与连续表达对比混在同一个荟萃分析中。

Orient every effect so that a positive value has the same interpretation. Keep RNA-seq and microarray effects separate if their scales cannot be justified as commensurate; combine standardized effects only as a prespecified secondary analysis.

统一所有效应方向，使正值含义一致。若无法证明 RNA-seq 与微阵列效应尺度可比，应分别分析；仅可在预先指定的次要分析中合并标准化效应。

Report for every study:

- sample counts by response group;
- coefficient, standard error, confidence interval, and exact model;
- adjustment set and nuisance method;
- missing-data handling;
- probe or gene identifier used for TACSTD2;
- whether the effect was primary or sensitivity-only.

每项研究应报告：

- 各应答组样本数；
- 系数、标准误、置信区间及完整模型；
- 调整变量集及非期望变异处理方法；
- 缺失数据处理方式；
- TACSTD2 所用探针或基因标识；
- 该效应属于主要分析还是仅用于敏感性分析。

## 11. Random-effects meta-analysis / 随机效应荟萃分析

For study \(i\), let \(y_i\) be the TACSTD2 effect and \(s_i\) its standard error:

\[
y_i \sim N(\theta_i, s_i^2), \qquad
\theta_i \sim N(\mu, \tau^2).
\]

对研究 \(i\)，令 \(y_i\) 为 TACSTD2 效应，\(s_i\) 为其标准误：

\[
y_i \sim N(\theta_i, s_i^2), \qquad
\theta_i \sim N(\mu, \tau^2).
\]

Use a random-effects model because populations, tissues, platforms, treatments, and response definitions are expected to differ. Estimate \(\tau^2\) with a defensible method such as REML. Prefer Hartung–Knapp-type inference when the number of studies is small, while noting that no method makes inference reliable with very few studies.

由于人群、组织、平台、治疗和应答定义通常存在差异，应采用随机效应模型。可使用 REML 等合理方法估计 \(\tau^2\)。当研究数量较少时，优先考虑 Hartung–Knapp 类推断，但必须说明：研究极少时，任何方法都无法保证可靠推断。

Report:

- pooled effect with 95% confidence interval;
- prediction interval when estimable;
- \(\tau^2\) and \(I^2\), with uncertainty-aware interpretation;
- forest plot containing every eligible study;
- exact effect scale and positive direction.

报告内容包括：

- 合并效应及 95% 置信区间；
- 可估计时报告预测区间；
- \(\tau^2\) 和 \(I^2\)，并谨慎解释其不确定性；
- 包含所有合格研究的森林图；
- 明确效应尺度及正方向含义。

Do not use a fixed-effect model merely because a heterogeneity test is non-significant. With few studies, heterogeneity tests have low power.

不得仅因异质性检验不显著而采用固定效应模型。研究数量较少时，异质性检验的效能很低。

## 12. Prespecified sensitivity analyses / 预先指定的敏感性分析

Run, where data permit:

1. known-batch covariate adjustment only;
2. known batch plus SVA;
3. known batch plus a justified RUV method;
4. ComBat-seq as secondary analysis for eligible raw-count cohorts;
5. leave-one-study-out meta-analysis;
6. separate RNA-seq and microarray analyses;
7. separate tissue, treatment, response-definition, or timepoint strata;
8. alternative valid TACSTD2 probe summaries;
9. exclusion of studies at high risk of response-label or sample-overlap error.

在数据允许时开展：

1. 仅调整已知批次；
2. 已知批次加 SVA；
3. 已知批次加有充分依据的 RUV；
4. 对符合条件的原始计数队列，将 ComBat-seq 作为次要分析；
5. 逐一剔除研究的荟萃分析；
6. RNA-seq 与微阵列分别分析；
7. 按组织、治疗、应答定义或时间点分层；
8. 使用其他有效的 TACSTD2 探针汇总规则；
9. 排除应答标签错误或样本重复风险较高的研究。

Interpret sensitivity analyses by effect direction, magnitude, uncertainty, and heterogeneity—not by whether \(P < 0.05\). Label post hoc analyses clearly.

敏感性分析应依据效应方向、大小、不确定性和异质性解释，而不是仅判断 \(P < 0.05\)。所有事后分析必须明确标注。

## 13. Diagnostics and stopping rules / 诊断与停止规则

Before accepting a study-level estimate, verify:

- sample identities and response counts match the source metadata;
- normalization and expression distributions are plausible;
- TACSTD2 is detectably measured on that platform;
- the model matrix is full rank;
- no single sample dominates the coefficient;
- residual/PCA structure is not strongly associated with an unmodeled known technical factor;
- nuisance factors are not near-perfect proxies for response;
- correction has not produced impossible values for the downstream model.

接受研究层面效应前应确认：

- 样本身份及应答组数量与来源元数据一致；
- 归一化结果及表达分布合理；
- 该平台能够可靠检测 TACSTD2；
- 设计矩阵满秩；
- 不存在单一样本主导系数；
- 残差或 PCA 结构未与未建模的已知技术因素强相关；
- 非期望因子并非应答的近乎完美代理；
- 校正结果未产生不符合下游模型要求的数值。

Stop and mark the contrast non-estimable if the response coefficient is aliased, all response information is confined to one technical batch, sample identities cannot be reconciled, or TACSTD2 is not validly measured. More aggressive correction is not a remedy.

若应答系数发生别名/线性依赖、所有应答信息仅存在于一个技术批次、样本身份无法核对，或平台不能有效测量 TACSTD2，应停止分析并标记该对比不可估计。更激进的批次校正不能解决这些问题。

## 14. Reproducibility and reporting / 可重复性与报告

Archive:

- GEO accession and download date;
- raw-file checksums and retrieval commands;
- sample manifest and all label harmonization rules;
- annotation source and version;
- exact normalization, filtering, model, ComBat-seq, SVA, and RUV settings;
- software and package versions;
- study-level effects before meta-analysis;
- code-generated QC, forest, influence, and sensitivity outputs;
- protocol deviations with reasons.

归档内容：

- GEO 编号及下载日期；
- 原始文件校验和及下载命令；
- 样本清单及全部标签统一规则；
- 注释来源和版本；
- 归一化、过滤、模型、ComBat-seq、SVA 及 RUV 的完整参数；
- 软件与包版本；
- 荟萃分析前的研究层面效应；
- 由代码生成的质控、森林图、影响分析和敏感性分析结果；
- 方案偏离及原因。

The final report must distinguish:

1. raw-data preprocessing;
2. known-batch adjustment;
3. latent unwanted-variation adjustment;
4. visualization-only correction;
5. study-level biological effect estimation;
6. cross-study meta-analysis.

最终报告必须区分：

1. 原始数据预处理；
2. 已知批次调整；
3. 潜在非期望变异调整；
4. 仅用于可视化的校正；
5. 研究层面的生物学效应估计；
6. 跨研究荟萃分析。

This separation prevents a corrected expression matrix from being mistaken for evidence. The evidence is the set of identifiable, uncertainty-bearing, independently replicated study-level TACSTD2 effects.

这种区分可避免把“校正后的表达矩阵”误当作证据。真正的证据是：一组可识别、包含不确定性、并在独立研究中得到重复验证的 TACSTD2 研究层面效应。
