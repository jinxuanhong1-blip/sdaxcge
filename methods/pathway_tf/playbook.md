# Pathway and TF activity playbook / 通路与转录因子活性分析手册

## English

### Scope

This workflow estimates pathway and transcription-factor (TF) activity after
`CLDN4` or `TACSTD2` perturbation and tests the same programs in tumors.

- PROGENy: pathway activity from footprint genes.
- DoRothEA: TF activity from signed TF-target regulons.
- decoupleR: the inference framework; report ULM and MLM separately.
- GSVA or ssGSEA: sample-level enrichment for Hallmark EMT, Hallmark IFN-alpha,
  Hallmark IFN-gamma, and a versioned custom junction signature.

There is no MSigDB Hallmark "junction" set. Do not relabel a GO or custom
junction set as Hallmark. Freeze its gene list, direction, source, species, and
retrieval date in the GMT description or an adjacent manifest.

### Files and input contract

Copy and edit the files in `templates/`:

1. `config.template.yml` -> `config.yml`.
2. `metadata.template.tsv` -> one row per biological sample.
3. `contrasts.template.tsv` -> the prespecified models.
4. `custom_signatures.template.gmt` -> replace every placeholder with approved,
   species-matched gene symbols.
5. `signature_manifest.template.tsv` -> provenance for every tested signature.

The expression RDS must contain one numeric matrix with genes in rows and
samples in columns. Use unique HGNC symbols for human or MGI symbols for mouse.
Columns must exactly equal `metadata$sample_id`, in the same order. Supply
normalized log expression (for example, log2-CPM or a variance-stabilized
assay), not raw counts, TPM mixed across cohorts, z-scored genes, differential
expression statistics, or a prefiltered "significant genes" matrix.

Keep all adequately expressed genes for footprint scoring. Resolve duplicate
symbols before scoring (prefer an annotation-based rule or the row with highest
mean expression), and record the rule. Never silently average mixed gene-ID
types.

### Design before scoring

For perturbations:

- Keep biological replicates as separate samples. Technical replicates may be
  collapsed before scoring.
- Encode perturbation direction explicitly (knockout/knockdown versus
  overexpression), reagent, dose, time, cell line, batch, and matched control.
- Compare within target and time point. Include batch and cell line in the
  model; use donor as a block for paired primary samples.
- The experimental unit, not each cell in a single-cell experiment, determines
  replication. For single-cell data, pseudobulk by sample and relevant cell
  type first.

For tumors:

- Use continuous `CLDN4` or `TACSTD2` log expression as the primary exposure.
  A median high/low split loses information and is visualization-only.
- Adjust at minimum for tumor purity, cohort, and known processing batch.
  Add cancer type, stage, or major molecular subtype when the cohort spans
  them. Do not adjust for a variable that is plausibly downstream of the target
  unless the estimand requires it.
- Analyze repeated specimens with patient blocking or a mixed model. Analyze
  cancer types separately before a prespecified meta-analysis when effects may
  be heterogeneous.
- Because the exposure gene can occur in a signature, repeat tumor analyses
  after removing `CLDN4`/`TACSTD2` from each gene set. This leave-one-gene-out
  sensitivity check guards against a tautological association.

### Reproducible setup

Use a project-local R library (for example, `renv`) and record package versions.
Install Bioconductor packages with `BiocManager`; obtain `msigdbr` from CRAN.
The templates require a current GSVA release with parameter objects
(`ssgseaParam`/`gsvaParam`).

```r
install.packages(c("BiocManager", "msigdbr", "renv", "yaml"))
BiocManager::install(c("decoupleR", "GSVA"))
renv::snapshot()
```

Run from the repository root:

```bash
Rscript methods/pathway_tf/R/score_decoupler.R methods/pathway_tf/config.yml
Rscript methods/pathway_tf/R/score_gsva.R methods/pathway_tf/config.yml
```

Each script writes scores, an RDS, and `sessionInfo()` under
`analysis.output_dir`. The GSVA script also writes per-set expression overlap.

### Prespecified scoring choices

1. Run PROGENy with the top 500 footprint genes per pathway.
2. Run DoRothEA confidence levels A-C. Use A-B as a stricter sensitivity
   analysis; do not select confidence levels after seeing associations.
3. Run ULM and MLM. Treat agreement in direction as robustness; do not average
   methods before testing. A positive score means inferred activation relative
   to other samples under that method, not measured protein activation.
4. Use ssGSEA as the default rank-based enrichment method. GSVA is an acceptable
   prespecified alternative. Do not choose whichever gives the smaller
   p-value.
5. Score Hallmark EMT, IFN-alpha response, and IFN-gamma response separately.
   Do not merge the two IFN sets. Score directional junction sets separately
   (`JUNCTION_UP`, `JUNCTION_DOWN`) rather than mixing opposing genes.

Fit each score as an outcome using the formula in
`contrasts.template.tsv`. For a perturbation, the tested coefficient is
perturbed minus matched control after covariate adjustment. For tumors, it is
the score change per unit of target log expression. If a paired/block design is
present, use `limma::duplicateCorrelation`, a fixed patient effect when
estimable, or an appropriate mixed model.

Standardize each score across samples only when effect sizes must be compared
across methods; retain and archive unscaled scores. Report the standardized
beta, 95% confidence interval, raw p-value, and Benjamini-Hochberg FDR. Correct
within a clearly declared family (for example, all tested PROGENy pathways for
one contrast); also provide a global correction as a sensitivity analysis when
many method/signature/contrast combinations are screened.

### Quality control and interpretation

Before inference:

- Confirm perturbation efficacy from independent evidence and expression of the
  intended target.
- Report library/assay QC, sample outliers, gene-set overlap, and the number of
  regulon targets detected per sample.
- Inspect score distributions and batch effects without removing samples based
  on the desired biological result.
- Flag sets below the configured minimum size and do not interpret them.
- Keep PROGENy, TF, and gene-set enrichment results as complementary evidence.
  Correlated methods are not independent validation.

A result is prioritized when its direction is reproducible across biological
replicates, robust to the prespecified covariates and sensitivity analyses, and
consistent between perturbation and tumors. A tumor association alone does not
establish that CLDN4 or TACSTD2 caused the activity change. Validate key claims
with orthogonal assays or an independent cohort.

Minimum report:

- organism, annotation release, expression transform, filtering, and duplicate
  handling;
- perturbation construct/direction, controls, dose, time, sample size, and
  model formula;
- PROGENy footprint size; DoRothEA confidence levels; decoupleR methods and
  versions;
- MSigDB release, Hallmark set names, custom-junction provenance, and overlap;
- effect, interval, p-value, FDR family, sensitivity analyses, and all exclusions.

Primary references and APIs:
[PROGENy](https://saezlab.github.io/progeny/),
[DoRothEA](https://saezlab.github.io/dorothea/),
[decoupleR](https://saezlab.github.io/decoupleR/),
[GSVA](https://bioconductor.org/packages/GSVA/), and
[MSigDB Hallmark](https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=H).

---

## 中文

### 范围

本流程用于评估 `CLDN4` 或 `TACSTD2` 扰动后的通路与转录因子（TF）活性，
并在肿瘤样本中检验相同生物学程序。

- PROGENy：根据响应基因足迹推断通路活性。
- DoRothEA：根据有方向的 TF-靶基因调控网络推断 TF 活性。
- decoupleR：推断框架；ULM 与 MLM 应分别报告。
- GSVA 或 ssGSEA：计算 Hallmark EMT、Hallmark IFN-alpha、Hallmark
  IFN-gamma 及自定义细胞连接（junction）签名的样本级富集分数。

MSigDB Hallmark 中没有“junction”基因集。不得把 GO 或自定义连接签名标为
Hallmark。应冻结并记录其基因列表、上下调方向、来源、物种和获取日期。

### 文件与输入约定

复制并修改 `templates/` 中的文件：

1. `config.template.yml` 改为 `config.yml`；
2. `metadata.template.tsv`：每个生物学样本一行；
3. `contrasts.template.tsv`：预先指定统计模型；
4. `custom_signatures.template.gmt`：把全部占位符替换为经审核且物种匹配的
   基因符号；
5. `signature_manifest.template.tsv`：记录每个签名的来源与版本。

表达量 RDS 必须保存“基因 × 样本”的数值矩阵。人样本使用唯一 HGNC symbol，
小鼠使用 MGI symbol。矩阵列名必须与 `metadata$sample_id` 完全一致且顺序相同。
输入应为标准化后的 log 表达量（如 log2-CPM 或方差稳定化 assay），不能使用
原始 counts、跨队列混合 TPM、按基因 z-score、差异分析统计量或仅含“显著基因”
的矩阵。

足迹分析应保留所有充分表达的基因。评分前按预先规定的注释规则解决重复 symbol
（或保留平均表达最高的一行）并记录规则；不可静默混用不同基因 ID。

### 评分前的研究设计

扰动实验：

- 生物学重复必须保留为独立样本；技术重复可在评分前合并。
- 明确记录敲除/敲低或过表达方向、试剂、剂量、时间、细胞系、批次及匹配对照。
- 在同一靶点和时间点内比较，模型中纳入批次与细胞系；配对原代样本以供体分块。
- 单细胞实验的重复单位是样本而非单个细胞；先按样本和相关细胞类型做
  pseudobulk。

肿瘤队列：

- 主要暴露变量使用连续的 `CLDN4` 或 `TACSTD2` log 表达量。中位数高低分组会
  损失信息，仅用于可视化。
- 至少校正肿瘤纯度、队列和已知处理批次。跨癌种时纳入癌种、分期或主要分子
  亚型。除非目标估计量要求，否则不要校正可能处于靶点下游的变量。
- 重复取样需按患者分块或采用混合模型。若预期癌种间异质，应先分癌种分析，再
  进行预先指定的 meta-analysis。
- 暴露基因可能属于被评分签名，因此应从每个基因集中删除 `CLDN4`/`TACSTD2`
  后重复肿瘤分析，作为 leave-one-gene-out 敏感性检验。

### 可复现运行

使用项目级 R 环境（如 `renv`）并保存版本。Bioconductor 包用
`BiocManager` 安装，`msigdbr` 从 CRAN 安装。模板要求支持
`ssgseaParam`/`gsvaParam` 参数对象的新版 GSVA。

```r
install.packages(c("BiocManager", "msigdbr", "renv", "yaml"))
BiocManager::install(c("decoupleR", "GSVA"))
renv::snapshot()
```

在仓库根目录运行：

```bash
Rscript methods/pathway_tf/R/score_decoupler.R methods/pathway_tf/config.yml
Rscript methods/pathway_tf/R/score_gsva.R methods/pathway_tf/config.yml
```

脚本在 `analysis.output_dir` 下输出分数、RDS 和 `sessionInfo()`；GSVA 脚本还
输出每个基因集与表达矩阵的重叠基因数。

### 预先指定的评分策略

1. PROGENy 每条通路使用排名前 500 的足迹基因。
2. DoRothEA 主分析使用 A-C 置信等级；A-B 作为更严格的敏感性分析。不得根据
   结果选择等级。
3. 同时运行 ULM 与 MLM，方向一致视为稳健性证据，但检验前不得取平均。正分数
   表示该方法下的推断活化，并不等同于实测蛋白活化。
4. 默认使用基于排序的 ssGSEA；也可预先指定 GSVA。不得根据较小 p 值事后选择。
5. Hallmark EMT、IFN-alpha 和 IFN-gamma 分别评分；两类 IFN 不合并。
   `JUNCTION_UP` 与 `JUNCTION_DOWN` 分开评分，不能混合方向相反的基因。

以每个活性分数为结局，采用 `contrasts.template.tsv` 中的公式。扰动分析检验
校正协变量后的“扰动减对照”；肿瘤分析检验靶基因 log 表达每增加一个单位对应的
分数变化。配对/分块设计可采用 `limma::duplicateCorrelation`、可估计的患者
固定效应或合适的混合模型。

仅在跨方法比较效应量时对各分数做样本间标准化，并保留原始分数。报告标准化
beta、95% 置信区间、原始 p 值和 Benjamini-Hochberg FDR。需明确多重检验家族
（例如某一 contrast 下全部 PROGENy 通路）；若筛选很多方法、签名与 contrast
组合，还应提供全局校正的敏感性结果。

### 质控与解释

推断前应：

- 用独立证据确认扰动有效，并确认目标基因表达变化；
- 报告 assay 质控、异常样本、基因集重叠及每个 regulon 检出的靶基因数；
- 检查分数分布和批次效应，不得依据期望生物学结果删除样本；
- 低于最小基因数的集合应标记并停止解释；
- PROGENy、TF 活性和基因集富集是互补证据；相关方法不构成独立验证。

优先结果应在生物学重复间方向一致，对预设协变量和敏感性分析稳健，并在扰动
实验与肿瘤中一致。单独的肿瘤相关性不能证明 CLDN4 或 TACSTD2 导致活性变化。
关键结论需用正交实验或独立队列验证。

最低报告内容包括：物种与注释版本、表达变换与过滤、重复 ID 处理；扰动方向、
对照、剂量、时间、样本量及模型公式；PROGENy 足迹大小、DoRothEA 置信等级、
decoupleR 方法与版本；MSigDB 版本、Hallmark 名称、自定义 junction 来源与
重叠数；效应值、区间、p 值、FDR 家族、敏感性分析及所有排除项。
