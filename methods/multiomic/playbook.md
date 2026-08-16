# TACSTD2/CLDN4 junction multi-omic integration playbook (2024–2026)

## 0. Decision question / 决策问题

**English.** Test one prespecified claim: a **TACSTD2/CLDN4 junction program** is associated with (i) immune exclusion and (ii) resistance to immune-checkpoint inhibition (ICI). Keep three objects distinct:

1. **Junction event:** directly observed split reads or an equivalently specific validated assay for a TACSTD2–CLDN4 junction.
2. **Junction program:** a frozen, direction-signed feature set learned without using immune or ICI outcomes.
3. **Proxy phenotype:** TACSTD2/CLDN4 expression, protein, or epithelial-state scores. A proxy may support mechanism but cannot establish the junction event.

Primary estimands are:

- immune exclusion: change in a prespecified exclusion endpoint per 1-SD higher program score;
- ICI resistance: log hazard ratio for PFS/OS, log odds ratio for non-response, or response-score difference per 1-SD higher program score;
- interaction, where powered: program × ICI treatment, preferred over a prognostic association in ICI-only cohorts.

**中文。** 预先规定并检验一个命题：**TACSTD2/CLDN4 连接（junction）程序**是否与（i）免疫排斥及（ii）免疫检查点抑制剂（ICI）耐药相关。必须区分：

1. **连接事件：** split reads 直接证据，或具有同等特异性的验证实验；
2. **连接程序：** 不使用免疫或 ICI 结局训练、特征方向已冻结的基因/蛋白集合；
3. **替代表型：** TACSTD2/CLDN4 表达、蛋白或上皮状态评分，只能提供机制支持，不能证明连接事件。

主要效应量为：程序评分每升高 1 个标准差时，免疫排斥终点的变化；PFS/OS 的 log(HR)、不应答的 log(OR) 或应答评分差异；样本量允许时优先检验“程序 × ICI 治疗”交互，而非仅在 ICI 队列中检验预后相关性。

## 1. Freeze the analysis contract / 冻结分析契约

Before inspecting outcomes, register:

- biological unit (patient, specimen, region, cell, or spot) and one row per unit;
- exact junction coordinates/build, strandedness, minimum unique split-read support, blacklist and mapping-QC rules;
- frozen program features, signs, weights, missing-feature policy, and score version;
- one primary exclusion endpoint and one primary ICI endpoint per cohort;
- covariates: at minimum tumor type, stage/line, tissue site, assay batch, purity or malignant-cell fraction, and clinically relevant treatment factors;
- repeated-sample handling, spatial neighborhood radius, scRNA pseudobulk definition, and patient-level resampling;
- compatible analysis strata for any meta-analysis (Section 8).

在查看结局前锁定：分析单位、连接坐标与参考基因组、split-read/QC 阈值、冻结的程序特征与方向/权重、每队列的主要终点、协变量、重复测量处理、空间邻域、单细胞 pseudobulk 定义，以及可合并的分析层（第 8 节）。

Do not let the same patients enter discovery and validation. Hash stable patient identifiers and check overlap across controlled and open releases. / 发现集和验证集不得包含同一患者；对稳定患者 ID 做哈希并检查受控与开放数据之间的重复。

## 2. Evidence tiers and data governance / 证据等级与数据治理

| Tier | Evidence / 证据 | Permitted claim / 可支持结论 |
|---|---|---|
| A | Direct junction reads plus immune/ICI endpoint in the same patient | Junction-associated exclusion/resistance |
| B | Frozen junction-program score plus endpoint in an independent cohort | Program-associated exclusion/resistance |
| C | TACSTD2/CLDN4 RNA/protein/spatial proxy only | Mechanistic consistency, not junction validation |
| D | In vitro perturbation or model system | Causal mechanism in that model only |

**Open data.** Store accession, release/version, license, reference build, download date, checksums, and processing provenance. Derived patient-level data must follow the source terms.

**Controlled data.** Keep raw reads, germline-risk fields, dates, and direct/quasi-identifiers inside the approved enclave. Export only disclosure-reviewed sufficient statistics: `estimate`, `se`, effective patient `n`, endpoint/contrast labels, QC counts, and non-identifying score summaries. Never move row-level controlled and open records into a common unsecured workspace.

**开放数据。** 记录 accession、版本、许可、参考基因组、下载日期、checksum 和处理流程。  
**受控数据。** 原始 reads、潜在胚系信息、日期及直接/间接标识符必须留在获批环境中；仅导出经过披露审查的效应量、标准误、有效患者数、终点/对比标签、QC 数量和不可识别的评分摘要。

Every result row must carry `access_class` (`open`, `controlled`, `mixed_summary`) and `evidence_tier`. Combining summary statistics does not change the most restrictive source obligations.

## 3. Modality-specific modules / 各模态分析模块

### 3.1 Bulk RNA-seq

- Re-align or use junction-aware outputs on one declared genome/transcriptome build. Require unique anchoring on both sides and inspect paralog/read-through artifacts.
- Quantify the event separately from gene expression. Report detection limit and fraction of samples with adequate junction coverage.
- Compute the frozen program on normalized expression; do not substitute the mean of TACSTD2 and CLDN4.
- Model exclusion with an orthogonal endpoint where possible (pathology, deconvolution, T-cell–inflamed signature). Adjust for purity and tumor type; report within-tumor estimates before pan-cancer estimates.
- For ICI, use treatment-aware models. In nonrandomized cohorts, distinguish predictive treatment interaction from prognostic association and adjust only prespecified baseline confounders.

### 3.2 scRNA-seq

- Call the junction only when read structure and depth support it; absence in sparse 3′ data is **not** evidence of absence.
- Identify malignant cells using CNV/genotype plus markers; avoid circular labeling from TACSTD2/CLDN4 program genes.
- Score cells with rank-based or control-gene-adjusted methods, then aggregate to patient × cell state (median, detection fraction, or pseudobulk). The inferential unit is the patient, not the cell.
- Primary exclusion readouts: malignant-program score versus patient-level CD8/T-cell abundance, dysfunctional state, or tumor–immune mixing. Fit mixed models only when patient-level uncertainty and clustering are preserved.

### 3.3 Protein

- Predefine analyte identity: targeted MS peptide, validated antibody, CITE-seq tag, RPPA, or IHC. TACSTD2 (TROP2) and CLDN4 abundance is Tier C unless the assay specifically measures a junction product.
- Normalize within platform/batch; retain lower-limit-of-detection flags. Use H-score/positive-cell fraction for IHC and patient-level summaries for single-cell protein.
- Test concordance with the frozen RNA program and immune endpoint separately. Do not convert protein abundance into “junction-positive.”

### 3.4 Spatial

- Segment tumor and immune compartments without using the tested program. Score malignant spots/cells and define exclusion before analysis.
- Suggested primary endpoint: distance from each CD8 cell/spot to the nearest malignant boundary, or observed tumor–CD8 contact versus a label-permutation null. Alternatives must be labeled secondary.
- Summarize region results to patient-level estimates; bootstrap/permutations must resample patients at the outer level. Adjust for section area, tumor content, platform, and region-selection scheme.
- Serial sections are linked evidence, not independent replicates. Registration uncertainty and spot mixing must be reported.

## 4. Frozen scoring and cross-dataset z-scores / 冻结评分与跨数据集 z 分数

Let feature \(x_{ij}\) be feature \(j\) in sample/unit \(i\), with frozen sign \(s_j\) and optional weight \(w_j\).

1. Normalize features using a modality-appropriate pipeline.
2. Standardize each feature **within dataset × modality × biological compartment**, using reference samples defined without outcomes:
   \[
   z_{ij}=(x_{ij}-\mathrm{median}_j)/(1.4826\,\mathrm{MAD}_j)
   \]
   Use mean/SD only when prespecified. Winsorize only by a frozen rule.
3. Score \(S_i=\sum_j w_js_jz_{ij}/\sum_j|w_j|\). Require a frozen minimum feature coverage (recommended ≥70%); otherwise mark missing.
4. Orient every score so higher values mean stronger junction program, then z-standardize the **final score within dataset** for effect estimates per 1 SD.

Cross-dataset z-scores align scale, not biology. They do not remove platform, composition, tissue, or batch differences and must not be pooled at row level. / 跨数据集 z 分数只统一量纲，不消除平台、组织组成、癌种或批次差异，禁止把不同数据集的样本行直接拼接后推断。

Report both robust-z and rank-percentile sensitivity analyses when zero inflation or saturation is substantial. For scRNA/spatial, standardize after patient-level aggregation for patient-level inference.

## 5. MOFA and consensus scoring / MOFA 与共识评分

### MOFA as a triangulation model

Use MOFA/MOFA2 only in participants with sufficiently connected views; record the missing-view pattern. Inputs are one matrix per view with the same patient key. Scale features within view, cap highly redundant features, and prevent a large transcriptomic view from dominating by feature selection or view weighting.

- Fit without immune/ICI outcomes.
- Choose factor number and convergence rules before outcome testing; run multiple seeds.
- Label the “junction-program factor” only after testing loading enrichment for the frozen signed program and factor stability under leave-one-dataset/view-out analysis.
- Align factor sign to the frozen program. A post hoc factor selected because it best predicts response is exploratory.
- Validate factor scores out of sample or project held-out samples with fixed loadings. Report variance explained by view and missing-view sensitivity.

MOFA is not required for cohorts with one view or weak cross-view overlap. Those cohorts contribute modality-specific estimates instead.

### Consensus score

For each patient with at least two eligible view scores \(Z_{iv}\):

\[
C_i=\frac{\sum_v q_{iv}r_v Z_{iv}}{\sum_v q_{iv}r_v}
\]

where \(q_{iv}\) is prespecified sample QC (0–1) and \(r_v\) is a frozen view-reliability weight estimated without outcomes. Also report:

- unweighted median of available oriented view z-scores;
- number and identity of views;
- sign concordance and maximum pairwise disagreement;
- leave-one-view-out scores.

Do not impute a direct junction event from a consensus of proxies. Primary analysis should require a frozen set of views; an “available-view” score is sensitivity analysis because its meaning changes with missingness.

**中文要点。** MOFA 必须在不使用免疫/ICI 结局的条件下拟合；因子命名依据冻结程序的 loading 富集和稳定性，而不是依据应答预测最好。共识分数按样本 QC 与预先确定的模态可靠性加权，同时报告未加权中位数、模态数量、方向一致性和 leave-one-view-out。多个替代指标的共识不能推断直接连接事件。

## 6. Endpoint models / 终点模型

| Question | Preferred unit/model | Exported effect |
|---|---|---|
| Exclusion, continuous | patient-level linear/robust model | standardized beta, SE |
| Exclusion, binary | logistic model | log OR, SE |
| ICI response | logistic model with fixed response definition | log OR for non-response, SE |
| PFS/OS | Cox model after PH check | log HR for progression/death, SE |
| Treatment prediction | treatment × score model | interaction coefficient, SE |
| Spatial contact/distance | patient-level permutation/mixed model | standardized beta or log ratio, SE |

All effects must be oriented so positive means **more exclusion or more resistance**. Keep different effect scales separate. Use patient-clustered or patient-bootstrap SEs for repeated regions/cells. Report unadjusted and prespecified adjusted estimates, but designate one as primary.

## 7. Result-table plug-in contract / 结果表接入规范

Other agents should emit UTF-8 CSV/TSV files with one row per estimate. Required columns:

| Column | Meaning |
|---|---|
| `dataset_id`, `analysis_id` | stable dataset and unique analysis identifiers |
| `modality` | `bulk_rna`, `scrna`, `protein`, `spatial`, `multiomic` |
| `access_class`, `evidence_tier` | governance and A–D evidence level |
| `tumor_type`, `treatment`, `endpoint` | explicit population and endpoint |
| `biological_unit` | usually `patient`; flag cell/spot/region analyses |
| `score_name`, `score_version` | frozen score identity |
| `contrast` | e.g. `per_1sd_higher_program` |
| `effect_type` | `beta`, `log_or`, `log_hr`, `interaction_beta`, `log_ratio` |
| `estimate`, `se`, `n_patients` | oriented estimate, standard error, independent patient count |
| `direction` | must be `higher_is_more_exclusion_or_resistance` |
| `adjustment_set` | stable label such as `primary_v1` |
| `analysis_tier` | `primary`, `sensitivity`, or `exploratory` |
| `meta_group` | frozen compatibility group; blank means no meta-analysis |
| `qc_status` | `pass`, `warn`, or `fail` |
| `junction_measure` | `direct`, `program`, `proxy`, or `perturbation` |
| `endpoint_definition`, `time_origin` | frozen endpoint details; use `not_applicable` explicitly |
| `population` | stable cohort/eligibility stratum label |
| `participant_set_id` | non-identifying overlap-audit label |

Recommended: `n_events`, `n_features_used`, `view_count`, `genome_build`, `cohort_version`, `estimate_ci_low`, `estimate_ci_high`, `notes`, and provenance URI/checksum.

`combine_results.py` validates and harmonizes these rows. It computes two-sided p-values and within-compatible-group inverse-variance fixed-effect summaries, plus heterogeneity diagnostics. It deliberately does not fabricate a standard error, convert effect types, or merge incompatible groups.

## 8. Meta-analysis gate: when not to combine / Meta 分析门控：何时不能合并

Meta-analyze only rows sharing all prespecified compatibility fields:

- estimand and effect scale;
- endpoint definition and time origin/window;
- score version, contrast, and direction;
- biological unit and patient-level uncertainty;
- population/tumor stratum, treatment context, and primary adjustment set;
- evidence target (direct junction, program, or proxy);
- non-overlapping participants.

**Do not meta-analyze**:

- log(HR), log(OR), correlation, AUC, and standardized beta together;
- ICI response with untreated prognosis;
- direct junction evidence with expression/protein proxies as if equivalent;
- cell/spot-level naive SEs with patient-level estimates;
- pan-cancer and tumor-specific effects when composition defines the contrast;
- different response criteria, survival time origins, treatment lines, or score versions;
- overlapping public/controlled releases;
- observational treatment associations with randomized treatment interactions;
- estimates selected because they were significant.

These are “apples and oranges.” Present them as a stratified evidence matrix and test qualitative concordance. A random-effects model does not repair incompatible estimands. With compatible rows, report fixed-effect synthesis by default for a common estimand; add REML/Hartung–Knapp only when justified and sufficiently populated. With fewer than three studies, emphasize individual CIs; heterogeneity estimates are unstable.

**中文。** 只有效应尺度、终点定义、评分版本、分析单位、癌种/治疗场景、调整集、证据对象和患者来源均兼容时才合并。随机效应模型不能修复“苹果和橘子”的目标量不一致。不能合并时，应使用分层证据矩阵、方向一致性和各自置信区间，而非一个总效应。

## 9. Multiplicity, robustness, and falsification / 多重检验、稳健性与证伪

- One primary exclusion and one primary ICI test; control FDR within clearly named secondary families.
- Negative controls: housekeeping/random matched programs, immune-rich nonmalignant compartment scores, and outcomes not plausibly linked to ICI mechanism.
- Sensitivities: leave-one-dataset/view/tumor-out; direct-event-only; high-coverage-only; purity/malignant-fraction adjustment; rank score; alternate exclusion assay; overlap removal.
- Check nonlinearity, proportional hazards, influential cohorts, batch association, and score–purity correlation.
- For spatial/scRNA, permute at patient level or within valid exchangeability blocks.
- Separate “association,” “predictive interaction,” and “mechanistic perturbation” in conclusions.

## 10. Execution and reporting checklist / 执行与报告清单

1. Freeze hypothesis, score, endpoints, covariates, and meta groups. / 冻结假设、评分、终点、协变量和合并组。
2. Audit identity overlap and access restrictions. / 审核患者重复及访问限制。
3. Run modality modules independently; export contract-compliant tables. / 各模态独立分析并导出标准表。
4. Orient effects and validate patient-level uncertainty. / 统一效应方向并验证患者层级不确定性。
5. Fit MOFA/consensus only where views overlap; keep modality estimates. / 仅在模态重叠充分时拟合 MOFA/共识。
6. Run the combiner; resolve `fail` rows, review `warn` rows, and inspect excluded groups.
7. Report a tiered forest/evidence matrix, heterogeneity, missingness, QC, and sensitivity analyses.
8. Conclude “supported” only if independent Tier A/B evidence agrees across at least two suitable datasets and no dominant falsification test fails; otherwise use “suggestive,” “inconsistent,” or “not testable.”

Suggested command:

```bash
python3 methods/multiomic/combine_results.py \
  --input results/bulk.csv results/scrna.tsv results/protein.csv results/spatial.csv \
  --outdir combined_multiomic
```

The schematic in `schematic.mmd` maps this workflow. / `schematic.mmd` 给出完整流程图。
