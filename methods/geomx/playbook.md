# GeoMx WTA Analysis Playbook — Compartment ROIs, Mixed Effects, and Survival
# GeoMx WTA 分析手册 —— 分区 ROI、混合效应模型与生存分析

**Scope / 适用范围:** GeoMx Digital Spatial Profiler (DSP) Whole Transcriptome Atlas (WTA) studies of the
"GSE271689 class": FFPE tissue microarrays or sections, **multiple ROIs per patient**, each ROI segmented
into **tumor (PanCK/CK+), leukocyte (CD45+) and macrophage (CD68+)** areas of illumination (AOIs),
with **patient-level clinical endpoints** (OS/PFS) and **compartment-specific target questions**
(e.g. is `TACSTD2`/`CLDN4` tumor-restricted?).

**中文:** 本手册面向 "GSE271689 类" 的 GeoMx DSP 全转录组 (WTA) 研究：FFPE 组织芯片或切片，**每位患者多个
ROI**，每个 ROI 再按形态学标记分割为 **肿瘤 (PanCK/CK+)、淋巴细胞 (CD45+)、巨噬细胞 (CD68+)** 三个照明区域
(AOI)，配套 **患者层面的临床终点** (OS/PFS)，并需要回答 **分区特异性的靶点问题**（例如 `TACSTD2`/`CLDN4`
是否为肿瘤特异表达）。

**Status / 状态:** methods-only. This directory contains no data and runs nothing by itself. All code under
`templates/` is a starting scaffold to be adapted, not a validated pipeline.
本目录仅包含方法学内容，不含数据；`templates/` 下的代码是需要按项目改写的脚手架，而非开箱即用的流水线。

---

## Table of contents / 目录

1. [Study anatomy and the unit-of-analysis rule / 研究结构与分析单元规则](#1)
2. [Annotation schema / 注释表规范](#2)
3. [Quality control / 质控](#3)
4. [Normalization and batch handling / 归一化与批次处理](#4)
5. [Design and mixed-effects models / 试验设计与混合效应模型](#5)
6. [Compartment contrasts: TACSTD2 / CLDN4 / 分区对比：TACSTD2 与 CLDN4](#6)
7. [From AOI to patient: aggregation and reliability / 从 AOI 到患者：聚合与可靠性](#7)
8. [Overall survival models / 总生存模型](#8)
9. [Signature scoring and deconvolution / signature 评分与解卷积](#9)
10. [Power and design planning / 效能与设计规划](#10)
11. [Pitfall table / 常见陷阱表](#11)
12. [Template index and run order / 模板索引与运行顺序](#12)
13. [Reporting / 报告要求](#13)

---

<a id="1"></a>
## 1. Study anatomy and the unit-of-analysis rule / 研究结构与分析单元规则

### 1.1 The nesting hierarchy / 嵌套层级

```
cohort (Yale / UQ / external validation)
└── slide / TMA block                      <- batch, hybridization run, sequencing lane
    └── patient                            <- the independent experimental unit
        └── core (TMA) or tissue region    <- sampling of the tumor
            └── ROI (typically 2-6/patient)
                └── segment / AOI: PanCK+ | CD45+ | CD68+
                    └── ~18,000 WTA targets
```

A WTA count matrix therefore has AOIs, not patients, in its columns. In a GSE271689-class design
(4 ROIs × 3 segments) one patient contributes up to **12 columns**. Every statement about statistical
significance must reflect the number of *patients*, not the number of *columns*.

**中文:** WTA 计数矩阵的列是 AOI 而不是患者。在 GSE271689 类设计中（4 个 ROI × 3 个分区），一位患者最多贡献
**12 列**。任何显著性结论都必须以 *患者数* 而非 *列数* 为准。

### 1.2 The unit-of-analysis rule / 分析单元规则

> **Rule.** Inference about a *patient-level* factor (treatment response, OS, stage, cohort) has
> `n = number of patients`. Inference about a *within-tissue* factor (segment type, tumor centre vs invasive
> margin, ROI-level immune density) can use the repeated AOIs, but only through a model that carries a
> patient random effect.

Treating AOIs as independent replicates ("pseudoreplication") inflates the effective sample size, deflates
standard errors, and produces false positives at a rate that grows with the intraclass correlation.
This is the single most common error in DSP papers and is explicitly called out in the GeoMx DSP Data
Analysis User Manual, the GeomxTools workflow, and the single-cell analogue literature
(Zimmerman et al. 2021 — see `references.md` R14).

**中文规则:** 涉及 *患者层面* 因素（疗效、OS、分期、队列）的推断，`n = 患者数`；涉及 *组织内部* 因素（分区类型、
瘤中心 vs 浸润边缘、ROI 免疫密度）的推断可以利用重复 AOI，但必须通过带患者随机效应的模型。把 AOI 当作独立重复
（伪重复, pseudoreplication）会虚增样本量、压低标准误，假阳性率随组内相关系数 (ICC) 上升而增加——这是 DSP
文献中最常见的错误。

### 1.3 Effective sample size / 有效样本量

For `m` AOIs per patient with intraclass correlation `ρ`, the design effect is `DE = 1 + (m − 1)ρ`, and the
effective sample size is `n_eff = n·m / DE`. With `m = 4` and a typical patient-level `ρ ≈ 0.5`, 4 ROIs buy
the information of **1.6 independent samples**, not 4. Report `ρ` per gene set or per target; it is the
number that tells you whether more ROIs or more patients is the better use of budget.

**中文:** 每位患者 `m` 个 AOI、组内相关 `ρ` 时，设计效应 `DE = 1 + (m − 1)ρ`，有效样本量 `n_eff = n·m / DE`。
当 `m = 4`、`ρ ≈ 0.5` 时，4 个 ROI 只相当于 **1.6 个独立样本**。请在论文中报告 `ρ`，它直接回答 "加患者还是加
ROI" 的问题。

---

<a id="2"></a>
## 2. Annotation schema / 注释表规范

Lock this table before touching counts. Half of all downstream modelling problems are annotation problems.
**中文:** 在碰计数矩阵之前先锁定这张表；下游一半的建模问题其实是注释问题。

| Column / 列 | Type | Required | Notes / 说明 |
|---|---|---|---|
| `dcc_file` | chr | yes | DCC filename, joins to the count matrix / 与计数矩阵连接的键 |
| `slide_name` | chr | yes | Batch proxy. Also record hyb kit lot, sequencing run / 批次代理变量，同时记录杂交试剂批号与测序批次 |
| `scan_name`, `roi_id` | chr | yes | ROI identity within slide / 玻片内 ROI 标识 |
| `patient_id` | chr | yes | **The clustering variable** / **聚类变量** |
| `core_id` | chr | if TMA | Multiple cores per patient are a distinct nesting level / TMA 中一名患者多个 core 时构成额外层级 |
| `segment` | factor | yes | `Tumor` / `CD45` / `CD68` (use fixed level order) / 固定因子水平顺序 |
| `roi_type` | factor | optional | e.g. `center` / `margin` / `TLS-adjacent` |
| `area_um2`, `nuclei` | num | yes | AOI size drivers; needed for QC, weighting, and sanity checks / 用于质控、加权与合理性检查 |
| `cohort` | factor | yes | Discovery vs validation; never pool blindly / 发现集 vs 验证集，切勿盲目合并 |
| `os_time`, `os_event` | num/int | yes | Patient-level, one row per patient in a separate clinical table / 患者层面，另存临床表 |
| `pfs_time`, `pfs_event` | num/int | optional | Same / 同上 |
| `age`, `sex`, `stage`, `histology`, `treatment_line`, `pdl1_tps`, `smoking` | mixed | yes | Adjustment covariates / 校正协变量 |
| `sample_date`, `block_age_years` | date/num | recommended | FFPE age is a real driver of RNA degradation / FFPE 存放年限显著影响 RNA 降解 |

**Two structural checks before modelling / 建模前的两项结构检查**

1. **Confounding check.** Cross-tabulate `slide_name × cohort`, `slide_name × outcome group`,
   `segment × slide_name`. If a slide is perfectly confounded with the comparison of interest, no statistical
   method can separate them — say so in the paper rather than adding the term and hoping.
   **中文：混杂检查。** 交叉列表检查玻片与队列/结局组/分区是否完全混杂；若完全混杂，任何方法都无法分离，应在文中
   明确说明，而不是加个协变量装作解决了。
2. **Balance check.** Tabulate AOIs per patient per segment. Unbalanced `m` (e.g. one patient with 12 AOIs,
   another with 2) is fine for mixed models but must be reported, and it changes patient-level aggregation
   (Section 7).
   **中文：平衡性检查。** 统计每位患者每个分区的 AOI 数；不平衡对混合模型可接受但必须报告，并影响患者层面聚合方式。

---

<a id="3"></a>
## 3. Quality control / 质控

### 3.1 Segment (AOI) QC / AOI 级质控

Recommended thresholds (GeomxTools defaults; the public vignette relaxes several of these for its demo
dataset, so do not copy vignette values blindly):

| Metric / 指标 | Threshold / 阈值 |
|---|---|
| Raw reads / 原始 reads | ≥ 1,000 |
| % trimmed, % stitched, % aligned | ≥ 80% each |
| Sequencing saturation / 测序饱和度 | ≥ 50% |
| Negative probe count geomean / 阴性探针几何均值 | ≥ 10 |
| NTC well count / 空白对照孔 | ≤ 1,000 |
| Nuclei / 细胞核数 | ≥ 100 (absolute floor 20) |
| AOI area / 面积 | ≥ 5,000 µm² (absolute floor 1,000) |

CD45+ and especially CD68+ AOIs are systematically smaller and lower-signal than PanCK+ AOIs. Applying one
nuclei/area threshold across segments silently deletes the immune compartment. **Set thresholds per segment
type and report the number dropped per segment.**

**中文:** CD45+ 尤其 CD68+ 的 AOI 系统性地更小、信号更低。对所有分区用同一个细胞核/面积阈值，等于悄悄删掉免疫分区。
**请按分区分别设阈值，并分区报告剔除数量。**

### 3.2 Probe QC / 探针质控

Remove probes failing the Grubbs outlier test globally or in ≥20% of AOIs, then collapse probes to gene
targets (geometric mean). WTA has one probe per target for most genes, so this step mainly matters for the
negative probe set.
**中文:** 用 Grubbs 检验剔除全局离群或在 ≥20% AOI 中离群的探针，然后按几何均值聚合为基因靶标。

### 3.3 Limit of quantification (LOQ) and gene filtering / 定量下限与基因过滤

```
LOQ_i = max( 2, geomean(NegProbe_i) × geoSD(NegProbe_i)^2 )
```

- **Segment filter:** drop AOIs detecting < 1% (conservative: < 5%) of targets above LOQ.
- **Gene filter:** keep genes detected above LOQ in ≥ 5–10% of AOIs.
- **Do the gene filter *within segment type*.** A gene detected in 60% of tumor AOIs but 0% of CD68 AOIs is
  informative; a pooled 5% filter keeps it, and then a naive tumor-vs-CD68 test reports a huge fold change
  that is really a detection artefact. Build one gene list per segment and analyse the intersection for
  cross-compartment contrasts, the union for within-compartment work — and state which you used.

**中文:**
- **AOI 过滤:** 剔除高于 LOQ 的检出率 < 1%（保守 < 5%）的 AOI。
- **基因过滤:** 保留在 ≥ 5–10% 的 AOI 中高于 LOQ 的基因。
- **基因过滤必须分区进行。** 某基因在 60% 肿瘤 AOI 中检出、在 CD68 AOI 中 0% 检出是有意义的信息；但用合并后的 5%
  阈值保留它，再做肿瘤 vs CD68 检验，得到的巨大 fold change 其实是检出率假象。跨分区对比用各分区基因表的**交集**，
  分区内分析用**并集**，并在文中写明用了哪一种。

### 3.4 What to report / 需报告内容

AOIs in → out at each step, split by segment; median genes detected per segment; median AOI area and nuclei
per segment; the LOQ formula and constants.
**中文:** 每一步的 AOI 进出数量（按分区拆分）、各分区检出基因中位数、面积与细胞核中位数、LOQ 公式与常数。

---

<a id="4"></a>
## 4. Normalization and batch handling / 归一化与批次处理

### 4.1 Choose deliberately, then show your work / 明确选择并展示依据

Q3 (75th percentile) scaling is the vendor default and remains the most published choice, but two lines of
evidence argue against accepting it uncritically:

- van Hijfte et al. (iScience 2023) showed Q3 fails standard distributional checks on GeoMx data and does not
  remove signal-to-noise-driven structure; rank-based (quantile) normalization performed best in their
  benchmark.
- Bhuva et al. (Genome Biology 2024) showed library size is confounded with biology in spatial data —
  in GeoMx, AOI area and cell content *are* biology, so any size-factor method partially removes the signal
  you came for.
- GeoDiff models background explicitly (Poisson background + negative-binomial threshold size factors)
  instead of scaling, which is the most principled option when many targets sit near LOQ — as they do in
  small CD68+ AOIs.

**Practical recommendation.** Run Q3, TMM and quantile (and GeoDiff if immune AOIs are background-dominated);
compare with RLE plots and PCA coloured by segment / slide / area / detection rate; pick one **a priori** as
primary and report the others as a sensitivity analysis. Report the primary choice in the abstract-level
methods, and put the comparison figure in the supplement.

**中文:** Q3（75 分位）是厂商默认、也是文献最常用的方法，但有两条证据提示不能不加检验地接受：van Hijfte 等
(iScience 2023) 显示 Q3 在 GeoMx 数据上不满足常规分布检验、无法消除信噪比驱动的结构，秩基（分位数）归一化表现最好；
Bhuva 等 (Genome Biol 2024) 指出空间数据中文库大小与生物学信息混杂——在 GeoMx 中 AOI 面积与细胞含量本身就是生物学，
任何 size factor 方法都会削掉一部分你要找的信号。GeoDiff 用显式背景建模（Poisson 背景 + 负二项阈值 size factor）
替代简单缩放，在大量靶标接近 LOQ 时（小的 CD68+ AOI 常见）最为稳妥。**实操建议：** 同时跑 Q3/TMM/quantile
（免疫 AOI 背景占优时加 GeoDiff），用 RLE 图与按分区/玻片/面积/检出率着色的 PCA 比较，**事先**指定一种为主分析，
其余作为敏感性分析报告。

### 4.2 Normalize per segment type or jointly? / 分区分别归一化还是合并归一化？

- **Cross-compartment contrasts (tumor vs CD45 vs CD68):** you must normalize jointly, otherwise the contrast
  is undefined. Accept that some of the tumor-vs-immune difference is compositional, and control for it as
  described in Section 6.3.
- **Within-compartment analyses (e.g. tumor-segment expression vs OS):** normalizing each segment type
  separately is defensible and often better, because AOI size distributions differ so strongly between
  compartments. Published DSP work has done exactly this to avoid scaling bias between large tumor AOIs and
  small immune AOIs.
- Never mix: do not compute a cross-compartment log-fold change from separately normalized matrices.

**中文:**
- **跨分区对比**（肿瘤 vs CD45 vs CD68）：必须合并归一化，否则对比无定义；跨分区差异中一部分来自组分差异，按 6.3 节
  的方法加以控制。
- **分区内分析**（如肿瘤分区表达与 OS 的关系）：按分区分别归一化是合理甚至更优的做法，因为各分区 AOI 面积分布差异
  极大；已发表 DSP 研究即采用该策略以避免大肿瘤 AOI 与小免疫 AOI 之间的缩放偏倚。
- 二者不可混用：绝不能用分别归一化的矩阵计算跨分区 log fold change。

### 4.3 Slide / batch effects / 玻片与批次效应

Use `standR::findNCGs()` to define negative control genes (least variable across slides within biology) and
`geomxBatchCorrection()` with RUV4; choose `k` by RLE/PCA stability, not by maximizing DE counts. For
inference, prefer keeping the batch term **in the model** (fixed effect if estimable, random if not) over
producing a corrected matrix; use batch-corrected values for visualization and clustering only.

**中文:** 用 `standR::findNCGs()` 选取阴性对照基因，配合 `geomxBatchCorrection()` 做 RUV4；`k` 依据 RLE/PCA
稳定性选择，不得以 "DE 基因数最多" 为标准。做统计推断时优先把批次项**留在模型里**（可估计则用固定效应，否则用随机
效应），而不是先生成校正后的矩阵；校正后的矩阵只用于可视化和聚类。

---

<a id="5"></a>
## 5. Design and mixed-effects models / 试验设计与混合效应模型

### 5.1 Two archetypes / 两种模型原型

| Question / 问题 | Structure / 结构 | Random effects / 随机效应 |
|---|---|---|
| Compare segments **co-existing** in the same tissue (tumor vs CD45 vs CD68) / 比较**同一组织内共存**的分区 | Paired, within patient / 患者内配对 | Random intercept **and random slope**: `(1 + segment \| patient)` |
| Compare **mutually exclusive** patient groups (responder vs non-responder, high vs low OS) within one segment / 在同一分区内比较**互斥**的患者组 | Between patient / 患者间 | Random intercept only: `(1 \| patient)` |
| Interaction: does the segment difference differ by group? / 交互：分区差异是否随患者组变化 | Mixed / 混合 | `(1 + segment \| patient)`, test `segment:group` |

The random-slope rule for co-existing structures is the GeoMx DSP Data Analysis User Manual's own guidance
and is reproduced in the GeomxTools workflow. Its purpose is to let the segment effect vary between patients
so that the *global* segment effect is tested against between-patient variability in that effect, not against
AOI noise.

**中文:** "共存结构用随机斜率" 的规则出自 GeoMx DSP 数据分析用户手册，GeomxTools 官方流程同样沿用。其作用是允许
分区效应在患者间变化，使 *全局* 分区效应针对 "该效应的患者间变异" 而非 "AOI 噪声" 进行检验。

### 5.2 Model menu / 模型菜单

| # | Method / 方法 | Use for / 适用 | R |
|---|---|---|---|
| A | limma-voom + `duplicateCorrelation(block = patient)` | **Primary genome-wide DE.** standR-recommended; borrows variance across genes, handles unbalanced designs, fast on 18k × 1k. / **全基因组 DE 首选** | `edgeR::voomLmFit(counts, design, block = patient, sample.weights = TRUE)` |
| B | `variancePartition::dream` | Crossed or multiple random effects (patient **and** slide), strongly unbalanced `m`, or when you need per-gene random effects rather than one consensus correlation | `dream(vobj, ~ segment + (1\|patient) + (1\|slide), meta)` |
| C | `lmerTest::lmer` on log2 normalized counts | Targeted genes, confirmatory analysis, random slopes, explicit variance components / 靶向基因、确证分析、随机斜率、方差成分 | `lmer(expr ~ segment + (1 + segment\|patient))` |
| D | Patient×segment pseudobulk + limma / t-test | Any **between-patient** contrast; the simplest defensible option and a mandatory sensitivity analysis / 任何**患者间**对比的最简稳健做法，且应作为必做敏感性分析 | aggregate then `lmFit` |
| E | `GeomxTools::mixedModelDE` | Reproducing the vendor/DSPDA analysis for comparability with prior publications / 复现厂商流程以便与既往文献比较 | `mixedModelDE(..., modelFormula = ~ testRegion + (1 + testRegion \| slide))` |

Report A (or B) as primary and D as sensitivity. If A and D disagree in direction for a headline gene, D wins
— the discrepancy means the A result was driven by within-patient replication.

**中文:** 以 A（或 B）为主分析、D 为敏感性分析。若两者在关键基因上方向不一致，以 D 为准——不一致说明 A 的结果由
患者内重复所驱动。

### 5.3 A worked contrast set / 一组可直接使用的对比

Group-means parameterization keeps contrasts readable:

```r
design <- model.matrix(~ 0 + segment_group, data = meta)   # segment_group = segment × outcome group
# e.g. levels: Tumor.R, Tumor.NR, CD45.R, CD45.NR, CD68.R, CD68.NR
contr <- limma::makeContrasts(
  Tumor_vs_CD45   = Tumor.R + Tumor.NR - CD45.R - CD45.NR,      # compartment specificity (pooled)
  Tumor_R_vs_NR   = Tumor.R - Tumor.NR,                         # outcome effect within tumor
  CD45_R_vs_NR    = CD45.R  - CD45.NR,                          # outcome effect within leukocyte
  CD68_R_vs_NR    = CD68.R  - CD68.NR,                          # outcome effect within macrophage
  Interaction     = (Tumor.R - Tumor.NR) - (CD45.R - CD45.NR),  # does the outcome effect differ by compartment
  levels = design)
```

Adjust p-values with Benjamini–Hochberg **within each contrast**, never across the stacked set of all
contrasts, and report the number of contrasts tested.

**中文:** 用 group-means 参数化写对比最清晰。BH 校正应在**每个对比内部**进行，不要把所有对比堆在一起校正，同时报告
共检验了多少个对比。

### 5.4 Diagnostics that must accompany any mixed model / 混合模型必须附带的诊断

1. `consensus.correlation` from `duplicateCorrelation` — report it. Values near 0 mean patient clustering is
   negligible; values > 0.3 mean pseudoreplication would have been severe. If > 0.5, refit
   `duplicateCorrelation` on residuals of the first fit and report both.
2. Singular fits / non-convergence rate for `lmer` across genes — if > 5–10% of genes give singular fits with
   a random slope, the design cannot support the slope; drop to a random intercept and say so.
3. Satterthwaite (default in `lmerTest`) or Kenward–Roger degrees of freedom for small `n`. Do not use naive
   Wald z-tests with ~20–60 patients.
4. Variance components: report `σ²_patient`, `σ²_slide`, `σ²_residual` and ICC for headline genes.

**中文:** ①报告 `duplicateCorrelation` 的 `consensus.correlation`（>0.5 时对残差再跑一次并同时报告）；②报告
`lmer` 跨基因的奇异拟合/不收敛比例（带随机斜率时超过 5–10% 说明设计撑不起斜率，应退回随机截距并说明）；③小样本用
Satterthwaite（`lmerTest` 默认）或 Kenward–Roger 自由度，20–60 例患者时不要用朴素 Wald z 检验；④报告关键基因的
`σ²_patient`、`σ²_slide`、`σ²_residual` 与 ICC。

---

<a id="6"></a>
## 6. Compartment contrasts: TACSTD2 / CLDN4 / 分区对比：TACSTD2 与 CLDN4

`TACSTD2` (TROP2) and `CLDN4` are ADC/theranostic targets. The question is rarely "is it differentially
expressed" — it is **"is it tumor-restricted, in every patient, in every region?"** That is three separate
quantities and each needs its own statistic.

**中文:** `TACSTD2` (TROP2) 与 `CLDN4` 是 ADC/诊疗一体化靶点。真正的问题不是 "是否差异表达"，而是
**"是否在每位患者、每个区域都局限于肿瘤"**——这是三个独立的量，各自需要不同的统计量。

### 6.1 Three quantities / 三个量

| Quantity / 量 | Statistic / 统计量 | Why it matters / 意义 |
|---|---|---|
| **Q1 Level** / 表达水平 | Tumor-segment log2 expression vs the AOI's LOQ; % of tumor AOIs above LOQ | An ADC target must be *present*, not merely *enriched* / ADC 靶点必须**存在**，而不只是**富集** |
| **Q2 Specificity** / 特异性 | Paired tumor − immune log2FC per patient; mixed-model estimate with 95% CI; % patients with LFC > 0 | Off-tumor expression drives toxicity and blurs the therapeutic window / 瘤外表达关系到毒性与治疗窗 |
| **Q3 Heterogeneity** / 异质性 | Between-ROI variance within patient (`σ²_within`), ICC, % ROIs concordant with the patient call | Determines whether one biopsy suffices to classify a patient / 决定单次活检能否代表患者 |

Literature anchors for sanity-checking your numbers: in localized anal cancer WTA-DSP, `TACSTD2` was
~3.5 log2FC higher in tumor than TME segments (adjusted p < 0.0001); single-cell references put `TACSTD2`
in 70–95% of malignant epithelial cells and < 2% of stromal/immune/endothelial cells; `CLDN4` in PDAC is
~16-fold higher in cancer than normal pancreas, with low spatial correlation to immune and fibroblast markers.
An observed tumor-vs-immune LFC of, say, 0.8 in your data therefore deserves a spillover investigation before
it is interpreted as biology.

**中文:** 用于核对量级的文献锚点：局部晚期肛管癌 WTA-DSP 中 `TACSTD2` 肿瘤分区较 TME 分区高约 3.5 log2FC
（校正 p < 0.0001）；单细胞数据显示 `TACSTD2` 在 70–95% 的恶性上皮细胞中表达、在基质/免疫/内皮细胞中 < 2%；
PDAC 中 `CLDN4` 较正常胰腺高约 16 倍，且与免疫、成纤维标记的空间相关性很低。若你的数据只得到 0.8 的
tumor-vs-immune LFC，应先排查渗漏（spillover）再谈生物学解释。

### 6.2 The primary model / 主模型

Paired within patient, all three segments, one gene at a time (or the whole panel with method A/B and then
extract the targets):

```r
lmerTest::lmer(expr ~ segment + (1 + segment | patient), data = df)   # co-existing segments -> random slope
# contrasts: Tumor - CD45, Tumor - CD68, CD45 - CD68  (emmeans, Kenward-Roger)
```

Add `+ (1 | slide)` when patients span slides; add `+ log10(area_um2)` **only** if you can argue that area is
technical here — for compartment contrasts area is largely biological and adjusting for it removes real
signal.

**中文:** 采用患者内配对、三分区同时建模；跨玻片时加 `(1 | slide)`。是否加 `log10(area_um2)` 需谨慎：在分区对比中
面积主要是生物学量而非技术量，盲目校正会删掉真实信号。

### 6.3 Spillover: the control that decides the paper / 渗漏：决定结论成败的对照

UV-based segmentation is imperfect. Immune AOIs adjacent to dense epithelium contain some tumor transcripts,
so a "tumor-restricted" gene will show non-zero counts in CD45/CD68 AOIs. Without a spillover control, you
cannot distinguish "TACSTD2 is expressed by macrophages" from "the CD68 AOI contains 8% tumor area".
Four controls, in increasing order of rigor:

1. **Epithelial spillover index.** Compute a per-AOI score from epithelial-restricted genes not related to
   your target (`EPCAM`, `KRT8`, `KRT18`, `KRT19`, `CDH1`). In immune AOIs this score estimates tumor
   contamination. Report the correlation between the index and your target's immune-AOI expression: if
   `r` is high, the immune signal is contamination.
2. **Adjust or stratify.** Add the spillover index as a covariate in the immune-AOI model, or restrict the
   analysis to the immune AOIs in the lowest tertile of the index and show the contrast is unchanged.
3. **Deconvolution.** `SpatialDecon::runspatialdecon()` with `is_pure_tumor = (segment == "Tumor")` derives
   tumor profiles from the PanCK+ AOIs and appends them to the safeTME matrix, so the tumor-derived component
   of the immune AOI is modelled rather than attributed to immune cells. Report the estimated tumor fraction
   in CD45/CD68 AOIs — it is a QC metric for the segmentation itself.
4. **Orthogonal validation.** IHC/RNAscope on the same blocks. Correlate the DSP tumor-segment value with
   H-score (published TACSTD2 DSP-vs-IHC correlations are moderate, r ≈ 0.4, which is itself worth stating).

Whatever you do, phrase the conclusion accordingly: the tumor-vs-immune LFC is a **lower bound** on true
compartment specificity, because contamination is unidirectional (tumor into immune AOIs more than the
reverse, given typical tumor:immune area ratios).

**中文:** 紫外分割并不完美，紧邻致密上皮的免疫 AOI 会混入肿瘤转录本，因此 "肿瘤特异" 基因在 CD45/CD68 AOI 中也会
有非零计数。没有渗漏对照，就无法区分 "巨噬细胞表达 TACSTD2" 与 "该 CD68 AOI 含 8% 肿瘤面积"。四级对照（严格度递增）：
①**上皮渗漏指数**：用与靶点无关的上皮限定基因（`EPCAM`、`KRT8/18/19`、`CDH1`）计算每个 AOI 的评分，在免疫 AOI 中
即为肿瘤污染估计；报告该指数与靶点在免疫 AOI 中表达的相关性，`r` 高即提示为污染。②**校正或分层**：把渗漏指数作为
协变量，或仅取指数最低三分位的免疫 AOI 重复分析并显示结论不变。③**解卷积**：`SpatialDecon::runspatialdecon()` 设
`is_pure_tumor = (segment == "Tumor")`，从 PanCK+ AOI 推导肿瘤表达谱并并入 safeTME 矩阵，使免疫 AOI 中的肿瘤来源
成分被显式建模；报告 CD45/CD68 AOI 的肿瘤成分估计值，它本身就是分割质量的质控指标。④**正交验证**：同一蜡块的
IHC/RNAscope，与 DSP 肿瘤分区值做相关（已发表的 TACSTD2 DSP-vs-IHC 相关约 r ≈ 0.4，这个数字本身值得写出来）。
无论采用哪一级，结论都应表述为：tumor-vs-immune LFC 是真实分区特异性的**下界**，因为污染是单向占优的。

### 6.4 Reporting template for a target gene / 靶基因报告模板

> `TACSTD2` was detected above LOQ in **X/Y (Z%)** tumor AOIs, **a/b (c%)** CD45 AOIs and **d/e (f%)** CD68
> AOIs. In a linear mixed model with a patient random intercept and random segment slope, expression was
> **L log2 units** higher in tumor than CD45 AOIs (95% CI …, Kenward–Roger p = …) and **M** higher than CD68
> AOIs (…). The difference was directionally consistent in **k/n** patients. The epithelial spillover index
> explained **r²** of the residual TACSTD2 signal in immune AOIs; restricting to the lowest spillover tertile
> changed the estimate to **L′**. Within-patient between-ROI variance accounted for **ICC** of total variance,
> implying **q** ROIs are needed for a patient-level estimate with reliability ≥ 0.8.

**中文:** 报告应至少包含：各分区高于 LOQ 的 AOI 比例；混合模型给出的分区差异与置信区间及所用自由度方法；方向一致的
患者比例；渗漏指数的解释力与低渗漏子集的敏感性分析；ICC 及达到可靠性 0.8 所需 ROI 数。

---

<a id="7"></a>
## 7. From AOI to patient: aggregation and reliability / 从 AOI 到患者：聚合与可靠性

Survival analysis needs one number per patient per segment. How you get there is a real statistical choice,
not bookkeeping.

**中文:** 生存分析需要 "每位患者每个分区一个数值"。如何得到这个数值是一个实质性的统计选择，而非简单汇总。

| Method / 方法 | Formula | When / 适用 | Caveat / 注意 |
|---|---|---|---|
| Arithmetic mean of log2 / log2 均值 | `mean(log2 expr)` | Balanced `m`, low outlier risk | Sensitive to one bad AOI / 对单个异常 AOI 敏感 |
| Median / 中位数 | `median(log2 expr)` | Default when `m ≥ 3` / `m ≥ 3` 时的默认选择 | Loses efficiency at small `m` |
| Area/nuclei-weighted mean / 加权均值 | `Σ wᵢ xᵢ / Σ wᵢ`, `w = nuclei` | When AOI size varies widely | Weights re-introduce a technical axis / 权重会引入技术维度 |
| **LMM BLUP (recommended)** / 混合模型 BLUP（推荐） | patient random intercept from `lmer(expr ~ 1 + (1\|patient))` | Unbalanced `m`, `m` varying 1–8 | Shrinks patients with few AOIs toward the mean — this is a feature (measurement-error correction), but it must be disclosed / 会把 AOI 少的患者向总体均值收缩，这是特性（测量误差校正），但必须披露 |

**Reliability.** With ICC `ρ` and `m` AOIs, the reliability of the patient mean is the Spearman–Brown
quantity `R = mρ / (1 + (m−1)ρ)`. Unreliable exposures attenuate hazard ratios toward 1 — a null OS result
with `R = 0.4` is uninformative, not negative. Compute `R` and report it next to every OS hazard ratio.

**中文:** **可靠性。** ICC 为 `ρ`、每人 `m` 个 AOI 时，患者均值的可靠性 `R = mρ / (1 + (m−1)ρ)`。暴露变量不可靠会
把风险比向 1 衰减——`R = 0.4` 时的阴性 OS 结果是 "信息不足"，不是 "阴性"。请在每个 OS 风险比旁一并报告 `R`。

**Derived exposures worth prespecifying / 值得事先指定的衍生暴露:**

- `tumor_score` = patient-level tumor-segment expression (the ADC-target relevant quantity).
- `immune_score` = CD45 (or CD68) segment expression.
- `specificity_score` = `tumor_score − immune_score` (compartment contrast as a patient-level biomarker).
- `heterogeneity_score` = within-patient SD across ROIs (an ADC-relevant "will one biopsy do?" measure).

These four answer different clinical questions and should be tested as four prespecified exposures, not as a
fishing exercise.
**中文:** 上述四个衍生变量回答不同的临床问题，应作为四个事先指定的暴露分别检验，而非事后挑选。

---

<a id="8"></a>
## 8. Overall survival models / 总生存模型

### 8.1 Default: patient-level Cox / 默认：患者层面 Cox 模型

```r
coxph(Surv(os_time, os_event) ~ scale(tumor_score) + age + sex + stage + histology + strata(cohort),
      data = pdat)
```

- **One row per patient.** This is the primary analysis. `n` = patients, events = deaths.
- **Continuous exposure, scaled per SD**, so the HR reads "per 1 SD increase". Do not dichotomize for the
  primary analysis (Altman 1994; Polley & Dignam 2021).
- **Events per variable ≥ 10.** With 60 patients and 30 deaths you can afford ~3 covariates. State the count.
- **Stratify, don't adjust, for cohort/slide** when baseline hazards plausibly differ and you do not need the
  cohort effect estimate.
- **Check PH** with scaled Schoenfeld residuals (`cox.zph`); if violated, use a time-varying coefficient or
  report a restricted mean survival time difference instead of a single HR.
- **Check functional form** with restricted cubic splines (3 knots) before assuming log-linearity.

**中文:** 主分析为**每位患者一行**的 Cox 模型；暴露按 SD 标准化后以连续变量进入，HR 解释为 "每增加 1 个 SD"；
主分析不做二分（Altman 1994；Polley & Dignam 2021）。遵守 "每变量至少 10 个事件" 并报告事件数。队列/玻片用
`strata()` 分层而非校正。用 `cox.zph` 检验比例风险，违背时改用时依系数或 RMST 差值。用限制性立方样条（3 节点）
先检查函数形式，再假定对数线性。

### 8.2 If you model at AOI level / 若在 AOI 层面建模

Sometimes reviewers ask for it, or you want to use all AOIs without aggregating. Two valid options, with
different targets of inference:

```r
# Marginal / population-average: robust sandwich SE, patient as cluster
coxph(Surv(os_time, os_event) ~ expr + age + cluster(patient), data = aoi_dat)

# Conditional / cluster-specific: shared gamma or Gaussian frailty
coxme::coxme(Surv(os_time, os_event) ~ expr + age + (1 | patient), data = aoi_dat)
```

Both repeat the same event time for every AOI of a patient, which is a duplication of the outcome, not extra
information about it. They control the standard errors but they do **not** create statistical power, and the
frailty variance is poorly identified when every cluster has exactly one event. **Use them as sensitivity
analyses; keep the patient-level model as primary.** Report which target of inference you mean: marginal
(population-average) for the robust-SE model, conditional (given the patient's frailty) for `coxme`.

**中文:** 有时审稿人会要求或你希望不聚合直接用全部 AOI。两种可行方案的推断目标不同：稳健三明治方差
（`cluster(patient)`）给出边际/人群平均效应；共享脆弱模型（`coxme`）给出条件效应。二者都把同一事件时间在患者的每个
AOI 上重复，属于结局的复制而非新增信息：它们能修正标准误，但**不会**带来统计效能，且每个 cluster 只有一个事件时
脆弱方差几乎不可识别。**只作敏感性分析，主分析仍用患者层面模型**，并写明推断目标是边际还是条件效应。

### 8.3 Cut points, if you must / 若必须使用切点

1. Prefer prespecified cut points: median, tertiles, or a clinically anchored threshold.
2. If data-driven (`maxstat`), you must (a) report the number of candidate thresholds searched,
   (b) apply a multiplicity-corrected p-value (Lausen–Schumacher / `maxstat` correction), and
   (c) bootstrap-correct the HR for optimism. An uncorrected min-p HR is biased away from 1, badly, in
   small studies.
3. Best practice, as done in the Nature Genetics NSCLC spatial-signature study: **fix the cut point in the
   training cohort and apply it unchanged to validation cohorts**, and state whether validation tests were
   one- or two-sided.

**中文:** ①优先使用事先指定的切点（中位数、三分位或临床阈值）；②若用数据驱动（`maxstat`），必须报告搜索的候选切点
数、使用多重性校正 p 值（Lausen–Schumacher 校正）、并对 HR 做 bootstrap 乐观偏倚校正——小样本中未校正的最小 p 值
HR 会严重偏离 1；③最佳实践（如 Nature Genetics NSCLC 空间 signature 研究）：**在训练集固定切点、原封不动地用于验证
集**，并说明验证集检验是单侧还是双侧。

### 8.4 Genome-wide OS screening / 全基因组 OS 筛查

If you scan all ~18k genes per segment for OS association:
- Fit per-gene Cox models on patient-level aggregates (fast, honest) and BH-adjust within segment.
- Report the number of genes tested per segment and the FDR threshold.
- Expect very few survivors at 18k × 3 segments with n < 100 patients. Prespecified targets (TACSTD2, CLDN4)
  and signature scores are the analyses that will replicate; frame the genome-wide scan as hypothesis-generating.
- Any multi-gene score built on the same patients whose OS you then test must be validated externally, or at
  minimum with nested cross-validation and a bootstrap optimism-corrected C-index. Reporting the apparent
  C-index of a LASSO score fit on the same data is not a result.

**中文:** 若对每个分区扫描全部 ~1.8 万基因：用患者层面聚合值逐基因拟合 Cox，在分区内做 BH 校正；报告每个分区检验的
基因数与 FDR 阈值。n < 100 时 18k × 3 的扫描通常几无幸存者——能重复的是事先指定的靶点（TACSTD2、CLDN4）与
signature 评分，全基因组扫描应定位为假设生成。任何在同一批患者上构建、又在同一批患者上检验 OS 的多基因评分，必须
外部验证，至少要用嵌套交叉验证与 bootstrap 乐观校正的 C-index；把 LASSO 在训练集上的表观 C-index 当作结果是不成立的。

### 8.5 Other endpoint issues / 其他终点问题

- **Immortal time.** If tissue is collected during treatment or at progression, landmark the analysis at a
  fixed time after treatment start and exclude patients with events before the landmark.
- **Competing risks.** OS has none by definition; for cancer-specific death or time-to-next-treatment use
  Fine–Gray subdistribution models and report both cause-specific and subdistribution HRs.
- **Discrimination.** Report time-dependent AUC (e.g. at 12/24/36 months) and Harrell's C with bootstrap CIs,
  not just a log-rank p.
- **Kaplan–Meier plots** must carry numbers at risk and the censoring pattern.

**中文:** **永生时间偏倚**：若组织在治疗中或进展时采集，需设定 landmark 时点并排除该时点前发生事件的患者；
**竞争风险**：OS 本身无竞争风险，肿瘤特异死亡或至下一线治疗时间应使用 Fine–Gray 模型并同时报告病因特异 HR 与
亚分布 HR；**区分度**：报告时依 AUC（12/24/36 个月）与 Harrell C 及 bootstrap 置信区间，而不只是 log-rank p 值；
**KM 曲线**必须带风险人数表与删失情况。

---

<a id="9"></a>
## 9. Signature scoring and deconvolution / signature 评分与解卷积

- **Score within segment type.** A rank-based score (`singscore`) computed on a pooled tumor+immune matrix
  ranks compartment identity, not biology. Compute scores separately per segment, or include segment in the
  model.
- **singscore vs GSVA:** `singscore` is rank-based and single-sample stable (useful when AOIs are added
  later); `GSVA` is population-relative and shifts when the cohort changes. For a locked biomarker, prefer
  single-sample methods.
- **Deconvolution of immune AOIs:** `SpatialDecon` with `safeTME` (built to avoid cancer-expressed genes) and
  `is_pure_tumor` set from the PanCK+ AOIs; supply `cell_counts = nuclei` to get results on an absolute cell
  scale. Deconvolving a CD68+ AOI into 14 immune types is over-reach — the AOI was already gated on CD68;
  interpret the output as relative composition within the gated population.
- **CIBERSORTx/LM22** on the immune compartment is an acceptable cross-check (used in the NSCLC spatial
  multi-omics study) but LM22 was built for bulk PBMC-like mixtures; treat differences between it and
  safeTME as method variance, not biology.

**中文:** **分区内评分**：在肿瘤+免疫合并矩阵上算秩基评分（`singscore`）排的是分区身份而非生物学，应分区计算或把
分区纳入模型。**singscore vs GSVA**：前者秩基、单样本稳定（后续追加 AOI 不改变已有值），后者相对于队列、队列一变
数值就变；锁定型生物标志物优先用单样本方法。**免疫 AOI 解卷积**：用 `SpatialDecon` + `safeTME`（构建时已规避癌细胞
表达基因），`is_pure_tumor` 取自 PanCK+ AOI，`cell_counts = nuclei` 可得到绝对细胞尺度；把已按 CD68 门控的 AOI 再
解卷积成 14 种免疫细胞属于过度解读，应解释为门控群体内的相对组成。**CIBERSORTx/LM22** 可作交叉验证（NSCLC 空间多组学
研究即如此），但 LM22 面向类 PBMC 的 bulk 混合物，与 safeTME 的差异应视为方法学差异而非生物学差异。

---

<a id="10"></a>
## 10. Power and design planning / 效能与设计规划

### 10.1 The trade-off / 权衡

Variance of a patient-level mean from `m` AOIs: `Var = σ²_between + σ²_within / m`. Adding AOIs shrinks only
the second term, so returns collapse after `m ≈ 3–4`. Adding patients shrinks the whole thing linearly.
**Under a fixed budget, patients beat ROIs beyond 3–4 ROIs per patient** — except when the scientific question
*is* within-patient heterogeneity (Section 6, Q3), where ROIs are the point.

**中文:** 患者均值的方差为 `Var = σ²_between + σ²_within / m`；增加 AOI 只压缩第二项，`m ≈ 3–4` 之后收益急剧递减，
而增加患者数是线性收益。**预算固定时，超过每人 3–4 个 ROI 后应优先增加患者**——除非研究问题本身就是患者内异质性
（见 6.1 的 Q3），那时 ROI 才是重点。

### 10.2 Simulate, don't formula-hunt / 用模拟而非套公式

Closed-form power for a three-level mixed model with unbalanced clusters and 18k tests is not worth deriving.
`templates/R/07_power_simulation.R` simulates: patients → ROIs → segments → counts, with configurable
`σ²_patient`, `σ²_roi`, `σ²_resid`, dropout, and effect size, then reports empirical power for
(a) the segment contrast, (b) a between-group contrast, and (c) an OS hazard ratio at a given event rate and
exposure reliability.

**中文:** 三层嵌套、簇不平衡、上万次检验的情形不值得推导解析公式。`templates/R/07_power_simulation.R` 直接模拟
患者→ROI→分区→计数的生成过程（可配置各层方差、脱落率与效应量），给出三类经验效能：分区对比、组间对比、以及在给定
事件率与暴露可靠性下的 OS 风险比。

### 10.3 Anchors from the literature / 文献锚点

The Nature Genetics NSCLC study (GSE271689-class Yale cohort) profiled **4 ROIs per tumor** across three
segments, with 131 patients in the transcriptomic arm, and validated signatures in two external cohorts.
That is roughly the scale at which compartment-specific OS signatures have been shown to validate externally.
A 20-patient single-cohort DSP study can characterize compartment biology; it cannot establish a prognostic
biomarker.

**中文:** Nature Genetics 的 NSCLC 研究（GSE271689 类的 Yale 队列）为**每个肿瘤 4 个 ROI** × 三分区，转录组部分 131
例患者，并在两个外部队列验证。这大致就是分区特异 OS signature 能够外部验证的规模量级。20 例单队列的 DSP 研究可以刻画
分区生物学，但不足以确立预后标志物。

---

<a id="11"></a>
## 11. Pitfall table / 常见陷阱表

| # | Pitfall / 陷阱 | Why it breaks / 后果 | Fix / 处理 |
|---|---|---|---|
| 1 | t-test / naive limma across AOIs / 直接对 AOI 做 t 检验 | Pseudoreplication; FPR grows with ICC / 伪重复，假阳性随 ICC 上升 | Patient random effect or pseudobulk (§5.2) |
| 2 | One QC threshold across segments / 各分区共用一个质控阈值 | Silently deletes CD68 compartment / 悄悄删掉 CD68 分区 | Per-segment thresholds, report drops (§3.1) |
| 3 | Pooled gene-detection filter / 合并计算检出率过滤 | Detection artefacts become "fold changes" / 检出率假象被当成差异表达 | Per-segment gene lists (§3.3) |
| 4 | Cross-compartment LFC from separately normalized matrices / 用分别归一化的矩阵算跨分区 LFC | Undefined contrast / 对比无定义 | Joint normalization for cross-compartment (§4.2) |
| 5 | "TACSTD2 is expressed in macrophages" without spillover control / 无渗漏对照就下 "巨噬细胞表达" 结论 | Segmentation contamination misread as biology / 分割污染被当成生物学 | Spillover index + SpatialDecon (§6.3) |
| 6 | Random slope on a design that can't support it / 设计撑不起随机斜率 | Singular fits, unstable p-values / 奇异拟合、p 值不稳 | Check convergence rate, fall back to intercept (§5.4) |
| 7 | AOI-level Cox as primary / 以 AOI 层面 Cox 为主分析 | Event duplication; n inflated / 事件复制、n 虚增 | Patient-level primary, clustered/frailty as sensitivity (§8.2) |
| 8 | Min-p cut point, uncorrected / 未校正的最小 p 值切点 | HR biased away from 1 / HR 严重偏离 1 | Prespecify, or correct + bootstrap (§8.3) |
| 9 | Apparent C-index of a score fit on the same data / 同一数据拟合又报告表观 C-index | Optimism, not performance / 是乐观偏倚而非性能 | External or nested CV + optimism correction (§8.4) |
| 10 | Slide confounded with cohort, "adjusted for" anyway / 玻片与队列完全混杂却仍 "校正" | Non-identifiable; estimate is arbitrary / 不可识别，估计任意 | Declare the confounding (§2) |
| 11 | Adjusting compartment contrasts for AOI area / 对分区对比校正 AOI 面积 | Removes biological signal / 删掉真实生物学信号 | Only adjust when area is demonstrably technical (§6.2) |
| 12 | BH across all contrasts stacked together / 把所有对比堆在一起做 BH | Wrong error rate family / 错误的检验族 | BH within contrast, report contrast count (§5.3) |
| 13 | Null OS result with unreliable exposure / 暴露不可靠时报告阴性 OS | Attenuation mistaken for absence / 把衰减当成无关联 | Report reliability `R` alongside HR (§7) |
| 14 | Ignoring FFPE block age / 忽略蜡块年限 | Degradation correlates with era, which correlates with treatment / 降解与年代相关，年代又与治疗相关 | Record and check as covariate (§2) |

---

<a id="12"></a>
## 12. Template index and run order / 模板索引与运行顺序

```
methods/geomx/
├── playbook.md                       # this file
├── references.md                     # annotated 2024-2026 + foundational citations
├── public_data.md                    # GSE271689 / GSE292098 only; no restricted data
├── README.md
├── checklists/reporting_checklist.md
└── templates/
    ├── config/
    │   ├── study_config.example.yml  # paths, segment names, contrasts, endpoints, thresholds
    │   └── annotation_template.csv   # annotation table with required columns
    ├── R/
    │   ├── 00_setup.R                # packages, config loader, session info
    │   ├── utils_geomx.R             # LOQ, per-segment filtering, aggregation, ICC
    │   ├── 01_qc_dcc_to_spe.R        # DCC/PKC -> GeoMxSet -> QC -> SpatialExperiment
    │   ├── 02_normalization_batch.R  # Q3 / TMM / quantile comparison, RUV4, diagnostics
    │   ├── 03_de_compartment.R       # limma-voom + duplicateCorrelation; dream; lmer
    │   ├── 04_target_genes.R         # TACSTD2 / CLDN4: level, specificity, heterogeneity, spillover
    │   ├── 05_patient_aggregation.R  # AOI -> patient scores, ICC, reliability
    │   ├── 06_survival_os.R          # Cox (patient-level primary), clustered/frailty, validation
    │   ├── 07_power_simulation.R     # nested simulation for design planning
    │   └── 08_fetch_public_geo.R     # public GSE271689 / GSE292098 DCC only
    └── methods_text/
        ├── methods_paragraph_en.md
        └── methods_paragraph_zh.md
```

Run order: `08` (public DCC) → `00 → 01 → 02 → 03 → 04 → 05 → 06`, with `07` used at design time.
Each script reads `config/study_config.yml`, writes to `results/`, and is safe to run standalone once
the previous step's outputs exist. If public phenotype has no OS table, stop after `04`.

**中文:** 运行顺序为 `08`（公开 DCC）→ `00 → 01 → 02 → 03 → 04 → 05 → 06`，`07` 在设计阶段使用。
每个脚本读取 `config/study_config.yml`，输出到 `results/`，在上一步产物存在时可独立运行。若公开表型
没有 OS 表，做到 `04` 即止。

---

<a id="13"></a>
## 13. Reporting / 报告要求

See `checklists/reporting_checklist.md` for the itemized list (bilingual) and
`templates/methods_text/` for fill-in methods paragraphs. The four things most often missing from DSP
papers, and the four a reviewer should ask for:

1. **The clustering variable and how it entered every model** (patient random effect? pseudobulk? nothing?).
2. **Per-segment QC and filtering counts**, so the reader knows the immune compartment survived.
3. **The normalization choice with a comparison**, not just "Q3 as per the manufacturer".
4. **Patients and events**, stated separately from AOIs, wherever an n appears.

**中文:** 逐条清单见 `checklists/reporting_checklist.md`，方法学段落模板见 `templates/methods_text/`。
DSP 论文最常缺失、也最应被审稿人追问的四点：①聚类变量是什么、如何进入每一个模型（患者随机效应？pseudobulk？还是没处理？）；
②分区拆分的质控与过滤计数，以证明免疫分区确实存活；③归一化方法的选择依据与比较，而不是一句 "按厂商建议用 Q3"；
④凡出现 n 的地方，患者数与事件数必须与 AOI 数分开陈述。

---

## References / 参考文献

Full annotated list with DOIs and PMIDs in [`references.md`](references.md). Anchor citations for this
playbook: standR (Liu et al., *Nucleic Acids Res* 2024), library-size confounding (Bhuva et al.,
*Genome Biol* 2024), GeoMx normalization bias (van Hijfte et al., *iScience* 2023), SpatialDecon
(Danaher et al., *Nat Commun* 2022), dream (Hoffman & Roussos, *Bioinformatics* 2021), pseudoreplication
(Zimmerman et al., *Nat Commun* 2021), spatial NSCLC signatures (Aung et al., *Nat Genet* 2025),
continuous-biomarker statistics (Polley & Dignam, *J Nucl Med* 2021), optimal-cutpoint bias
(Altman et al., *JNCI* 1994), smiDE contamination/spatial-correlation modelling (Vasconcelos et al.,
*Genome Biol* 2026).
