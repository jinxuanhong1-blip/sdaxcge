# Small-cohort ICI biomarker ML playbook / 小样本 ICI 生物标志物机器学习手册

This playbook is for exploratory prediction of immune-checkpoint inhibitor (ICI)
response from small GEO transcriptomic cohorts. It is not a recipe for a
clinically deployable test. The primary goals are honest uncertainty,
leakage-resistant evaluation, and a prespecified route to independent
validation.

本手册适用于使用 GEO 小样本转录组队列探索性预测免疫检查点抑制剂（ICI）疗效。
它不是临床可用检测的开发捷径。首要目标是如实表达不确定性、避免数据泄漏，并预先规定
独立验证路径。

## 1. Start with the estimand / 从研究问题开始

Before touching expression data, write down:

- **Population:** cancer type, treatment line, ICI drug/regimen, and eligibility.
- **Endpoint:** e.g. RECIST objective response at a fixed assessment, durable
  clinical benefit at a prespecified time, or time-to-event outcome. Do not mix
  endpoint definitions silently.
- **Sampling time:** pretreatment tissue only unless the intended clinical use
  explicitly includes on-treatment samples.
- **Unit of analysis:** patient, not biopsy or array. Multiple samples from one
  patient must remain in the same fold.
- **Intended use:** prognostic enrichment, treatment-specific prediction, or a
  biological hypothesis generator. A single-arm ICI cohort cannot distinguish
  a predictive treatment-interaction biomarker from a general prognostic marker.

在处理表达矩阵前，先写清：

- **目标人群：** 癌种、治疗线次、ICI 药物/方案和纳入标准。
- **终点：** 例如固定评估时点的 RECIST 客观缓解、预设时长的持久临床获益，或生存终点；
  不可无说明地混用定义。
- **取样时点：** 除非目标用途明确要求治疗中样本，否则仅使用治疗前组织。
- **分析单位：** 患者而非活检或芯片；同一患者的多个样本必须始终位于同一折。
- **用途：** 预后富集、治疗特异性预测，还是机制假设生成。单臂 ICI 队列不能区分
  “治疗交互型预测标志物”和一般预后标志物。

Create a signed or version-controlled analysis specification before inspecting
outcome associations. Define the primary metric and all exclusions there.

在查看结局关联前，建立带版本记录的分析方案，并在其中指定主要指标和所有排除规则。

## 2. Data audit and preprocessing / 数据审计与预处理

1. Reconstruct patient/sample identifiers and detect repeated measures.
2. Verify labels against source metadata; record ambiguous or missing labels
   rather than guessing.
3. Record platform, study, processing date, site, tissue source, purity,
   histology, treatment, and other plausible technical or clinical confounders.
4. Normalize from the least processed data available using a method appropriate
   to the platform. Freeze probe-to-gene mapping and duplicate-gene aggregation.
5. Apply outcome-independent sample and gene quality-control rules.
6. Keep an immutable raw-to-analysis manifest and report every exclusion.

1. 重建患者与样本标识，识别重复测量。
2. 对照原始元数据核验标签；模糊或缺失标签应记录，不可猜测。
3. 记录平台、研究、处理日期、中心、组织来源、纯度、组织学类型、治疗方案，以及可能的
   技术或临床混杂因素。
4. 从可获得的最原始数据开始，采用适合平台的标准化方法；冻结探针到基因的映射及重复
   基因聚合规则。
5. 采用与结局无关的样本和基因质控规则。
6. 保留不可变的“原始数据→分析数据”清单，并报告每项排除。

### Batch is a confounder, not merely noise / 批次是混杂因素，不只是噪声

Cross-tabulate outcome against study, platform, site, processing batch, tissue
source, and collection period. Plot embeddings colored separately by outcome
and batch. If response is nearly determined by batch, no statistical correction
can identify a biological response signal from these data.

将结局分别与研究、平台、中心、处理批次、组织来源和采集时期交叉制表；在降维图上分别
按结局和批次着色。如果疗效几乎由批次决定，则任何统计校正都无法从该数据中识别生物学
疗效信号。

Batch handling must occur **inside each training fold**:

- Estimate scaling, filtering thresholds, missing-value imputation, PCA, and
  batch-adjustment parameters on training samples only.
- Apply frozen parameters to the held-out fold. Never run ComBat, feature
  selection, or normalization jointly on all samples before CV.
- Prefer designs that preserve identifiability: grouped folds by study/site,
  leave-one-study-out evaluation, or models including prespecified batch
  covariates when future batch labels will be available.
- Do not use batch correction that requires outcome labels from the validation
  set. Report performance by batch/study as well as pooled performance.

批次处理必须在**每个训练折内部**完成：

- 仅用训练样本估计缩放、过滤阈值、缺失值填补、PCA 和批次校正参数。
- 将冻结后的参数应用于留出折。禁止在交叉验证前对全体样本共同运行 ComBat、特征选择
  或标准化。
- 优先采用可识别的设计：按研究/中心分组划折、逐研究留一验证，或在未来可获得批次
  标签时纳入预先指定的批次协变量。
- 不得使用需要验证集结局标签的批次校正。除汇总性能外，还应分批次/研究报告性能。

## 3. Feature policy / 特征策略

For small cohorts, use a small, prespecified biological feature space (for
example, pathway scores defined without this cohort's outcomes) whenever
possible. If gene-level modeling is unavoidable:

- perform variance/missingness filtering within the training fold;
- perform any univariate screening within the training fold;
- use regularization and tune it only in inner CV;
- keep clinical covariates and benchmark against a clinical-only model;
- record gene identifiers and transformations so an external sample can be
  processed without refitting.

对小队列，应尽可能使用小规模、预先指定的生物学特征空间（例如不利用本队列结局定义的
通路评分）。若必须进行基因层面建模：

- 方差/缺失率过滤只能在训练折内完成；
- 任一单变量筛选只能在训练折内完成；
- 使用正则化，并仅在内层交叉验证中调参；
- 保留临床协变量，并与“仅临床”模型比较；
- 记录基因标识和变换，使外部样本无需重新拟合即可处理。

## 4. Nested cross-validation / 嵌套交叉验证

Use nested CV to separate model selection from performance estimation.

**Outer loop:** estimates generalization. Use stratified grouped folds, grouping
by patient and, where scientifically appropriate, study/site/batch. Choose the
number of folds based on the minority class; every test fold must contain both
classes. With very few events, repeated outer splits may describe instability
but do not manufacture independent information.

**Inner loop:** on each outer-training set, repeat the entire trainable workflow:
feature filtering/selection, preprocessing, hyperparameter search, and any
model choice. Select one configuration without seeing the outer test fold.

**Prediction ledger:** save exactly one out-of-fold prediction per patient per
outer repeat, together with fold, study, batch, observed label, and model
version. Compute metrics from this ledger, not from training predictions.

嵌套交叉验证用于分离“模型选择”和“性能估计”。

**外层循环：**估计泛化性能。采用按患者分组的分层折；在科学上合适时，还应按研究、
中心或批次分组。折数受少数类别样本数限制，每个测试折必须包含两类。事件极少时，重复
外层划分可描述不稳定性，但不会创造新的独立信息。

**内层循环：**在每个外层训练集中，重新执行全部可学习步骤：特征过滤/选择、预处理、
超参数搜索及模型选择。不得查看外层测试折。

**预测台账：**每次外层重复中，为每名患者保存且仅保存一个折外预测，同时记录折、
研究、批次、真实标签和模型版本。所有指标均从该台账计算，不得使用训练集预测。

Report ROC AUC and PR AUC with uncertainty, but emphasize proper scoring rules
(Brier score and log loss) and threshold-specific sensitivity/specificity when
clinically relevant. Obtain confidence intervals by patient-level bootstrap of
the full out-of-fold ledger (or by repeating the entire nested procedure);
acknowledge that CV predictions are correlated. Never choose a threshold on an
outer test fold.

报告 ROC AUC、PR AUC 及其不确定性，同时重视恰当评分规则（Brier 分数、对数损失）和
临床相关阈值下的敏感度/特异度。置信区间可对完整折外预测台账进行患者级 bootstrap
（或重复完整嵌套流程）；需说明 CV 预测彼此相关。禁止在外层测试折上选择阈值。

## 5. Why a 20-gene signature should not be trained on n=16
/ 为什么不应在 n=16 上训练 20 基因签名

This is not solved by choosing a penalized model. Sixteen patients provide at
most 16 independent outcome observations, often fewer effective observations
after class imbalance, missingness, and batch clustering. Selecting 20 genes
from thousands adds a large hidden search space. Many signatures will separate
the training data by chance, and leave-one-out CV remains highly variable while
reusing almost the same samples for selection. Coefficients, selected genes,
decision thresholds, apparent AUC, and calibration will all be unstable.

正则化模型不能解决这一问题。16 名患者最多只提供 16 个独立结局观测；考虑类别不平衡、
缺失和批次聚类后，有效信息更少。从数千基因中选择 20 个基因会引入巨大的隐性搜索空间。
大量签名会偶然分开训练数据；留一法仍会因反复使用几乎相同的样本进行筛选而高度不稳定。
系数、入选基因、决策阈值、表观 AUC 和校准都会不稳定。

For n=16, defensible outputs are descriptive plots, effect sizes with wide
intervals, assessment of a **fully locked pre-existing** score, and hypotheses
for a larger study. Do not derive and claim validation of a new multigene
signature in the same 16 patients. A sample-size target should instead be
derived from expected event prevalence, desired precision/calibration, number
of candidate parameters, and planned validation—not from a universal
events-per-variable slogan.

对于 n=16，较合理的产出包括描述性图、带宽置信区间的效应量、对**完全锁定的既有**
评分进行评估，以及为更大研究提出假设。不可在同 16 名患者中同时开发并宣称验证新的
多基因签名。样本量应根据预期事件率、目标精度/校准、候选参数数量和验证设计推导，而非
套用固定的“每变量事件数”口号。

## 6. Calibration and decision use / 校准与决策用途

Discrimination is not calibration. For every evaluation set:

- show a calibration plot with uncertainty and report calibration intercept
  (ideal 0), slope (ideal 1), Brier score, and log loss;
- avoid many-bin reliability plots in tiny samples; use a smooth curve with a
  bootstrap band and display the prediction distribution;
- if calibration is trained, fit Platt scaling or isotonic regression only
  within the inner training process. Isotonic regression is usually too flexible
  for very small samples;
- lock prevalence assumptions and any recalibration rule before external
  validation;
- use decision-curve analysis only with clinically justified threshold ranges,
  and show uncertainty. It cannot rescue an unvalidated model.

区分度不等于校准度。对每个评估集：

- 展示带不确定性的校准图，并报告校准截距（理想值 0）、斜率（理想值 1）、Brier 分数
  和对数损失；
- 极小样本不宜使用多分箱可靠性图；应采用带 bootstrap 区间的平滑曲线，并显示预测
  分布；
- 若训练校准器，Platt scaling 或保序回归只能在内层训练过程中拟合；保序回归通常对
  极小样本过于灵活；
- 在外部验证前锁定患病率假设及任何再校准规则；
- 决策曲线分析只应使用具有临床依据的阈值范围并展示不确定性；它不能挽救未经验证的模型。

## 7. Locked external validation: OAK/POPLAR / 锁定的外部验证：OAK/POPLAR

OAK and POPLAR are randomized atezolizumab-versus-docetaxel NSCLC trials and
are valuable for separating prognostic effects from treatment-by-biomarker
interaction. Molecular/clinical data access may be controlled (for example
through EGA); verify the current data-use terms and exact accession from the
trial data provider rather than copying an unverified identifier.

OAK 和 POPLAR 是阿替利珠单抗对比多西他赛的随机 NSCLC 试验，可用于区分一般预后效应
与“治疗×生物标志物”交互效应。分子和临床数据可能受控访问（例如通过 EGA）；应向试验
数据提供方核实现行数据使用条款和准确登录号，不要照搬未经核验的编号。

Prespecify this sequence before obtaining outcomes:

1. **Access and governance:** obtain data-use approval; define permitted
   endpoints, linkage, publication, and model-export rules.
2. **Transportability audit:** compare eligibility, histology, line of therapy,
   specimen timing, assay/platform, tissue handling, endpoint definitions, and
   missingness with the GEO development cohort.
3. **Freeze the artifact:** container/environment, gene mapping, normalization
   references, feature list, coefficients, intercept, calibration mapping,
   missing-feature policy, and decision thresholds. Hash and timestamp it.
4. **Blind processing:** process expression and clinical covariates without
   outcome access. Do not center against the external cohort or substitute genes
   based on their outcome association.
5. **Primary validation:** evaluate the locked model in the prespecified
   atezolizumab population, reporting discrimination, calibration, clinical
   operating points, and confidence intervals.
6. **Predictive claim:** in the randomized intention-to-treat population, fit a
   prespecified outcome model containing treatment, the continuous locked
   biomarker, and their interaction (plus prespecified stratification factors).
   For survival endpoints, check proportional-hazards assumptions and report
   absolute risk at clinically relevant times where possible. A significant
   result in the ICI arm alone does not establish treatment prediction.
7. **Trial roles:** preferably use one trial for locked validation and the other
   for confirmation. If pooling is necessary, preserve trial strata and report
   trial-specific estimates and heterogeneity.
8. **One allowed update:** if calibration transport fails, report the original
   locked result first. Any intercept-only recalibration is a model update and
   requires a subsequent untouched validation set.

在获得结局前预先规定以下顺序：

1. **访问与治理：**取得数据使用批准；明确允许使用的终点、数据链接、发表和模型导出规则。
2. **可迁移性审计：**比较 GEO 开发队列与外部队列的入组标准、组织学类型、治疗线次、
   取样时点、检测平台、组织处理、终点定义和缺失情况。
3. **冻结模型工件：**固定容器/环境、基因映射、标准化参考、特征列表、系数、截距、
   校准映射、缺失特征策略和决策阈值；生成哈希并加时间戳。
4. **盲态处理：**在无法访问结局的情况下处理表达和临床协变量。不得利用外部队列重新
   中心化，也不得按结局关联替换基因。
5. **主要验证：**在预设的阿替利珠单抗人群中评价锁定模型，报告区分度、校准度、临床
   操作点和置信区间。
6. **预测性主张：**在随机化意向治疗人群中，拟合预先指定的结局模型，包含治疗、连续型
   锁定标志物及其交互项（以及预设分层因素）。对于生存终点，检查比例风险假设，并尽可能
   报告临床相关时点的绝对风险。仅在 ICI 组显著不能证明治疗预测作用。
7. **试验分工：**优先用一个试验作锁定验证、另一个作确认。若必须合并，应保留试验分层，
   并报告各试验效应及异质性。
8. **仅允许一次更新：**若校准迁移失败，先报告原始锁定结果。任何仅截距再校准都属于
   模型更新，需在后续未使用的数据集中再次验证。

Register the analysis, publish a flow diagram and missing-data table, and report
results according to TRIPOD (and applicable extensions). Release the locked
code and model card where data-use agreements permit.

注册分析方案，发表样本流程图和缺失数据表，并依据 TRIPOD（及适用扩展）报告结果。
在数据使用协议允许时公开锁定代码和模型卡。

## 8. Minimal go/no-go checklist / 最小继续或停止清单

Proceed to model development only if all are true:

- the endpoint and patient-level split unit are unambiguous;
- each class is large enough to populate every inner and outer validation fold;
- outcome is not aliased with study/platform/batch;
- all learned preprocessing can be fitted within folds and applied prospectively;
- the candidate feature space and primary analysis are prespecified;
- an external validation population and locked transport pipeline are plausible.

仅当以下条件全部满足时才进入模型开发：

- 终点和患者级划分单位明确；
- 两类别样本数足以填充每个内层和外层验证折；
- 结局未与研究/平台/批次完全混淆；
- 所有可学习预处理均可在折内拟合，并可前瞻应用；
- 候选特征空间和主要分析已预先指定；
- 存在可行的外部验证人群和锁定迁移流程。

The companion `elastic_net_template.py` intentionally refuses tiny datasets.
Its guardrails are demonstration defaults, not a sample-size justification.
Passing them does not imply that a cohort is adequate.

配套的 `elastic_net_template.py` 会主动拒绝极小数据集。其保护阈值仅为演示默认值，
不是样本量论证；通过阈值不代表队列足够。
