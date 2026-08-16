# Hierarchical meta-analysis of TACSTD2 / CLDN4 in open ICI GEO cohorts
# TACSTD2 / CLDN4 开放 ICI GEO 队列的分层 / 贝叶斯合并

> **Methods only · 仅方法。** This folder is a template. It does not contain a
> finding, a recommended biomarker, or a pooled estimate from real patients.
> It reads `results/*.csv` **if those files already exist** from sibling
> per-cohort analyses. If they do not, the script prints the input contract and
> exits.
>
> **本目录是方法模板，不含结论。** 只在兄弟分析已经写出 `results/*.csv` 时读取；
> 否则打印输入约定并退出，不编造效应值。

**Open cohorts only · 仅开放队列。** The likelihood is restricted to public
GEO / ArrayExpress / Zenodo-style matrices. EGA, dbGaP, and DAC-controlled
trial transcriptomes (CheckMate, KEYNOTE, OAK/POPLAR, IMpower, SU2C-MARK, …)
are **catalogued as missing**, never downloaded, and enter the analysis only
through the MNAR sensitivity module.
**似然只含公开 GEO / ArrayExpress 等矩阵。** EGA、dbGaP、DAC 受控试验转录组
**只作为缺失队列登记**，不下载，只进入 MNAR 敏感性分析。

**Genes · 基因:** `TACSTD2` (TROP2) and `CLDN4`, run separately (one hierarchical
model per gene × endpoint × effect family).
**基因：** `TACSTD2`（TROP2）与 `CLDN4`，按「基因 × 终点 × 效应族」分别拟合。

---

## Table of contents · 目录

0. [Scope and open-only gate · 范围与开放门控](#0-scope-and-open-only-gate--范围与开放门控)
1. [Estimand · 估计目标](#1-estimand--估计目标)
2. [Input contract · 输入约定](#2-input-contract--输入约定)
3. [Hierarchical model · 分层模型](#3-hierarchical-model--分层模型)
4. [Shrinkage · 收缩](#4-shrinkage--收缩)
5. [Frequentist companions · 频率派对照](#5-frequentist-companions--频率派对照)
6. [Publication and selection bias · 发表 / 选择偏倚](#6-publication-and-selection-bias--发表--选择偏倚)
7. [Missing EGA / dbGaP trials · 缺失的 EGA / dbGaP 试验](#7-missing-ega--dbgap-trials--缺失的-ega--dbgap-试验)
8. [Prior sensitivity · 先验敏感性](#8-prior-sensitivity--先验敏感性)
9. [How to run · 如何运行](#9-how-to-run--如何运行)
10. [What not to claim · 不可声称](#10-what-not-to-claim--不可声称)
11. [References · 参考文献](#11-references--参考文献)

---

## 0. Scope and open-only gate · 范围与开放门控

**EN.** Sibling agents in this repository have already analysed many small
human ICI GEO series for `TACSTD2` and `CLDN4` (response / DCB / MPR / PFS).
Those per-cohort estimates are the **level-1 data** of this template. This
step does not re-download expression matrices and does not invent accessions.

A real, already-seen open set (catalogued without effect sizes in
[`example/known_open_ici_cohorts.csv`](example/known_open_ici_cohorts.csv)):

| Accession | Why it is in scope | Typical open endpoint |
| --- | --- | --- |
| GSE126044 | open bulk RNA-seq, NSCLC anti-PD-1 | responder / non-responder |
| GSE135222 | open bulk RNA-seq, NSCLC anti-PD-1/PD-L1 | PFS (DCB often PFS ≥ 180 d) |
| GSE166449 | open bulk RNA-seq, lung ICI | responder / non-responder |
| GSE207422 | open bulk RNA-seq, neoadjuvant PD-1 + chemo | MPR / RECIST |
| GSE136961, GSE93157 | open, but **genes absent** from the panel | not estimable; do not impute |

Controlled resources that stay in the **missing** registry (list only; no
download attempted in this repo): `EGAS00001005013` (OAK/POPLAR),
CheckMate-153 `EGAS00001007508`, SU2C-MARK `phs002822`, CheckMate 017/057/227/9LA
and IMpower150 trial transcriptomes when no open processed matrix exists,
and GSE205335 **raw** (`EGAD00001008703`).

`access_policy: open_only` (default) drops any input row whose `access` token
is EGA / dbGaP / controlled, or whose `cohort_id` starts with `EGA` / `phs`.
A GEO accession without an access column is treated as open. The drop is
written to `input_audit.csv` so the count is reconcilable.

**中文.** 本仓库已有多个小样本人源 ICI GEO 系列的 `TACSTD2` / `CLDN4` 单队列分析。
那些逐队列估计是本模板的**第一层数据**。本步骤不重新下载表达矩阵，不编造登录号。

上表中的开放队列属于范围；面板上没有这两个基因的系列（GSE136961、GSE93157）
记为不可估计，**不补值**。EGA / dbGaP 受控试验只进缺失登记。默认
`access_policy: open_only` 会丢掉受控行，并写入 `input_audit.csv`。

---

## 1. Estimand · 估计目标

**EN.** Two effect families are supported and **must not be pooled together**:

1. **`smd`** (default for GEO response analyses). Hedges' *g* (or Cohen's *d*)
   of gene expression, responder minus non-responder. Already on a standardized
   scale. This is what `notes/gpt_meta` and several GEO slices produce.
2. **`log_ratio`**. log hazard ratio or log odds ratio per 1 SD of expression
   (or another declared unit). Survival / logistic association.

**Orientation `harm`:** after harmonisation, a **positive** value means higher
`TACSTD2` or `CLDN4` is associated with a **worse** ICI outcome. Responder-minus-
nonresponder SMDs and response log-ORs are sign-flipped. Overall survival log-HRs
already have this sign and are left alone.

**One model per (gene, endpoint, effect family, adjustment set).** DCB, ORR, MPR,
PFS and OS are different estimands. Univariable and multivariable rows are
different estimands. Do not dump them into one forest plot.

**Exchangeability.** Cohorts are treated as exchangeable draws from a
between-cohort distribution. That is a modelling choice, not a fact. Lung
anti-PD-1 monotherapy pretreatment bulk RNA-seq is closer to exchangeable than
a mix of neoadjuvant chemo-IO, blood, spatial AOIs, and melanoma. If the
registry says the mix is wild, either subset or stop.

**中文.** 支持两种效应族，**禁止混池**：

1. **`smd`**：表达量的 Hedges' *g*（应答者减非应答者）。GEO 应答分析的默认。
2. **`log_ratio`**：每 1 SD 表达的 log HR / log OR。

**方向 `harm`：** 正值 = 高表达与更差 ICI 结局相关。应答 SMD / 应答 log-OR 会改号。
**每个（基因、终点、效应族、调整集）单独一个模型。** 队列被当作可交换的随机效应；
若适应证 / 治疗线 / 组织类型混杂过重，应先子集再合并，或停止合并。

---

## 2. Input contract · 输入约定

**EN.** Search `--results-dir` (default `results/`) for `*.csv`. Column names
are matched case-insensitively after stripping punctuation. Minimum fields:

| Column | Aliases (examples) | Role |
| --- | --- | --- |
| `cohort_id` | cohort, gse, accession | required |
| `gene` | gene_symbol, symbol | required; `TACSTD2` or `CLDN4` |
| `endpoint` | outcome | required; `DCB`, `ORR`, `OS`, `PFS`, `MPR`, … |
| `unit` | effect_unit, scale | `hedges_g` / `smd` / `per_sd` / `per_log2` / `high_vs_low` |
| `yi` | estimate, hedges_g, log_hr, beta | required\* |
| `sei` | se, stderr | required\* |
| `access` | access_type, availability | `open` / `geo` / `ega` / `dbgap` / `controlled` |

\* Instead of `(yi, sei)` you may supply `hr`/`or` + CI (log_ratio family only),
`yi` + CI, or `yi` + two-sided `pvalue` (least preferred; recorded as
`se_source=pvalue`).

Recommended: `ni`, `events`, `responders`, `x_sd`, `model_adjust`,
`cancer_type`, `ici_target`, `tissue_timing`, `sample_group`, `platform`,
`study`. `x_sd` is required to convert a per-log2 slope onto the per-SD scale.
`sample_group` marks overlapping patient sets (pre- vs on-treatment from one
trial); default policy keeps the pretreatment row.

Every input row is written to `input_audit.csv` with `kept` and `reason`.

**中文.** 在 `--results-dir`（默认 `results/`）下找 `*.csv`。列名忽略大小写与标点。
最低字段见上表。没有 `(yi, sei)` 时可用 CI 或 p 值反推标准误。`access` 用于开放
门控。`sample_group` 标记重叠样本。每一行的去留原因写入 `input_audit.csv`。

If `results/` is absent, nothing is written. That is the correct behaviour
until the per-cohort step exists.
若 `results/` 不存在，不写任何输出——这是在单队列步骤完成之前的正确行为。

---

## 3. Hierarchical model · 分层模型

**EN.** Normal–normal random-effects model

\[
y_i \mid \theta_i \sim \mathcal N(\theta_i, V_{ii}), \qquad
\theta_i \mid \mu,\tau \sim \mathcal N(\mu, \tau^2),
\]

with a Gaussian prior on \(\mu\) and a half-normal / half-Cauchy / uniform /
exponential prior on \(\tau\). \(V\) may be non-diagonal when overlapping
cohorts are kept (`overlap_policy: correlate`).

Conditional on \(\tau\) the model is Gaussian-conjugate, so
\(p(\mu \mid y, \tau)\), \(p(\theta_i \mid y, \tau)\) and the posterior
predictive for a new cohort are closed form. The \(\tau\) posterior is obtained
by one-dimensional quadrature on a grid (default 801 points). The reported
posterior of \(\mu\) is a finite mixture of those Gaussians. **No MCMC**, so
prior-sensitivity differences are not sampler noise.

Implementation: [`templates/hier_meta.py`](templates/hier_meta.py)
(`fit_bayes`). The R companion [`templates/hier_meta.R`](templates/hier_meta.R)
fits `metafor::rma.uni` (FE / DL / PM / REML+Hartung–Knapp) on the same
harmonised rows and prints a `brms` formula; it does not replace the Python
primary estimator.

**中文.** 正态–正态随机效应。对 \(\tau\) 做一维数值积分，对 \(\mu\) 与各队列
\(\theta_i\) 用共轭高斯，得到有限混合后验。无 MCMC。Python 为主要实现；R 脚本
只提供 metafor 对照。

Outputs: `pooled_estimates.csv`, `mu_posterior.csv`, `tau_posterior.csv`,
`heterogeneity.csv`, `run_manifest.json`.

---

## 4. Shrinkage · 收缩

**EN.** The posterior mean of \(\theta_i\) is a precision-weighted average of
the cohort estimate and \(\mu\). The reported shrinkage weight is

\[
B_i = \mathbb E\bigl[s_i^2 / (s_i^2 + \tau^2) \bigm| y\bigr].
\]

\(B_i \to 1\) means the cohort is almost replaced by the grand mean (typical of
n ≈ 16–30 GEO series). \(B_i \to 0\) means the cohort is left almost unpooled.
With K small, a spectacular single-cohort z-score that shrinks onto a near-null
\(\mu\) is the expected behaviour, not a bug. That is the point of partial
pooling on this problem.

Also written: leave-one-out refits (`leave_one_out.csv`) and a forest of
observed vs partially pooled intervals (`fig_forest_shrinkage.png` when
matplotlib is available).

**中文.** \(B_i\) 接近 1 表示该小队列几乎被总均数取代——这正是 n≈16–30 的 GEO
系列上分层模型要做的事。单队列显著、收缩后接近零，是预期行为。另有留一法与森林图。

---

## 5. Frequentist companions · 频率派对照

**EN.** The same \(y, V\) are passed through:

- common-effect (inverse-variance)
- DerSimonian–Laird, Paule–Mandel, REML \(\tau^2\)
- Hartung–Knapp–Sidik–Jonkman intervals (t reference; recommended at small K)
- Cochran Q, \(I^2\), Higgins–Thompson prediction interval

These are **cross-checks**, not a menu from which to pick the most significant.
When they disagree, believe the Bayesian prior-sensitivity table and the
prediction interval, not the smallest p-value.

**中文.** 固定效应、DL / PM / REML、Hartung–Knapp、Q / \(I^2\)、预测区间。
用作对照，不按「最显著」挑选。分歧时看先验敏感性表与预测区间。

---

## 6. Publication and selection bias · 发表 / 选择偏倚

**EN.** Open GEO is not a random sample of ICI biology. Small significant
series are more likely to be deposited with usable response labels; large
negative randomized trials sit behind EGA. The diagnostics below are
**sensitivity tools**. None of them "correct" the estimate.

| Tool | File / row | What it assumes | When it is weak |
| --- | --- | --- | --- |
| Egger regression | `smallstudy_bias.csv` | funnel asymmetry = small-study effect | K < 10, almost always here |
| Egger `inv_n` | same, `egger_predictor: inv_n` | better for log-OR (SE is a function of the estimate) | needs `ni` |
| Begg rank | same | rank correlation of standardised effect vs variance | very low power |
| Trim-and-fill | same | missing studies are exact mirrors | a *what-if*, never a corrected truth |
| PET / PEESE | same | effect is linear in SE or SE² | identified from the same small K |
| Step selection model | `selection_model_*` | non-significant cohorts observed with relative probability \(\omega\) | weakly identified; report the \(\omega\) grid, not just the MLE |
| Gene-panel empirical p | `gene_null_calibration.csv` | optional CSV of identically pooled negative-control genes | needs a prespecified null panel |

Do not treat a non-significant Egger p as evidence of no bias. With K = 3–8
the test cannot see the bias that the missing-EGA module is designed to bound.

**中文.** 开放 GEO 不是 ICI 生物学的随机样本。上表全部是敏感性工具，**没有一个
能「校正」估计**。K=3–8 时 Egger 不显著 ≠ 无偏倚；真正要看的是第 7 节的缺失试验
界限。

---

## 7. Missing EGA / dbGaP trials · 缺失的 EGA / dbGaP 试验

**EN.** This is the publication-bias problem that actually matters for ICI
transcriptomes. The large randomized trials that would dominate a well-designed
meta-analysis are usually **not in the likelihood**.

Set `missingness.registry` to a CSV with one row per **known** cohort, open or
not (`status` = `included` / `missing` / `gene_absent` / …; `access` =
`open` / `ega` / `dbgap`). The synthetic example is
[`example/cohort_registry.csv`](example/cohort_registry.csv). A catalog of real
accessions **without invented effects** is
[`example/known_open_ici_cohorts.csv`](example/known_open_ici_cohorts.csv).

Three complementary bounds:

1. **Pattern-mixture delta sweep** (`missingness_delta_sweep.csv`). Unobserved
   cohorts are given effects centred at \(\mu + \delta\). \(\delta = 0\) is MAR
   and barely moves \(\mu\). The **tipping point** is the smallest \(|\delta|\)
   at which the credible interval covers 0.
2. **Manski-style envelope.** All unobserved cohorts are placed at a low bound,
   then a high bound (default: observed min/max ± 2\(\tau\)). The envelope is
   the union of those credible intervals.
3. **Fail-safe K.** Smallest number of additional *null* cohorts of typical
   open-GEO precision that would push \(P(\mu\) has the observed sign\()\) below
   0.95. This replaces Rosenthal's discredited fail-safe N.

Standard errors for unobserved cohorts: `1/sqrt(events)` for survival,
`1/sqrt(n p (1-p))` for binary response, else the median observed SE. These
are order-of-magnitude devices for a sensitivity analysis, not design-based
variances.

**Do not download the EGA files to "fill in" \(\delta\).** If a DAC is later
granted, re-run the per-cohort step under that DAC and add the rows as
`access=open` only after the data-use agreement allows it. Until then the
honest sentence is: *the open-GEO pool is consistent with X; it would take
unobserved trials shifted by \(\delta\) to overturn it.*

**中文.** ICI 转录组真正的发表偏倚是：**本应占权重的大随机试验通常不在似然里。**
用登记表列出每一个已知队列（开放或受控）。三种界限：δ 扫描与 tipping point、
Manski 包络、推翻结论所需的「零效应开放精度」队列数 K。不要为了填 δ 去下载
EGA。在获得 DAC 之前，诚实表述是：开放 GEO 合并与 X 相容；要使结论翻转，未观测
试验需要偏移 δ。

---

## 8. Prior sensitivity · 先验敏感性

**EN.** With K small, \(\tau\) is barely identified and the half-normal scale
is doing real work. `prior_sensitivity.csv` refits the same data under:

- half-normal \(\tau\) scales 0.25 / 0.5 / 1.0
- half-Cauchy 0.5
- uniform \(\tau\) on [0, 2]
- skeptical \(N(0, 0.35^2)\) and flat \(\mu\)
- optional empirical log-normal on \(\tau^2\) (Turner / Rhodes tables) —
  **disabled until the table cell is filled in**

If the sign of \(\mu\) or whether the interval covers 0 flips across this
list, the data cannot support a stable claim. Report the table, not a single
prior.

A ROPE (region of practical equivalence) is declared on the effect scale
(default 0.2 for SMD, `log(1.1)` for log-ratio). `prob_in_rope` is the
posterior probability that the pooled effect is scientifically small, which
is a different question from \(P(\mu>0)\).

**中文.** K 小时 \(\tau\) 先验在真正起作用。必须并列报告多套先验。若正负号或
区间是否含 0 随先验翻转，数据撑不起稳定结论。ROPE 内概率问的是「效应是否小」，
与 \(P(\mu>0)\) 不是同一个问题。

---

## 9. How to run · 如何运行

```bash
# 1. No results yet: contract only, exit 0
python3 methods/hier_meta/templates/hier_meta.py --results-dir results

# 2. Smoke-test on the bundled SYNTHETIC example (numbers mean nothing)
python3 methods/hier_meta/templates/hier_meta.py --demo \
    --genes TACSTD2,CLDN4 \
    --outdir methods/hier_meta/out

# 3. Real per-cohort CSVs, both genes, open cohorts only
python3 methods/hier_meta/templates/hier_meta.py \
    --results-dir results --recursive \
    --config methods/hier_meta/config/config.json \
    --genes TACSTD2,CLDN4 \
    --endpoint DCB --effect-family smd \
    --outdir methods/hier_meta/out

# 4. Frequentist companion (optional; requires metafor)
Rscript methods/hier_meta/templates/hier_meta.R --demo \
    --gene TACSTD2 --endpoint DCB --outdir methods/hier_meta/out/TACSTD2

# 5. Self-test
python3 methods/hier_meta/templates/selftest.py
```

Dependencies: Python 3 + `numpy`, `scipy` (`matplotlib` optional for figures).
R + `metafor` optional.

When real CSVs arrive, keep one endpoint and one effect family per invocation.
A survival (`log_ratio`, `OS`) run is a **separate** call from a DCB SMD run.

**中文.** 无 `results/` 时只打印约定。`--demo` 跑捆绑的合成例子（数字无科学含义）。
真实分析时两个基因分别输出到子目录；每个终点、每个效应族单独调用一次。

---

## 10. What not to claim · 不可声称

**EN.**

- Do not claim a cross-cancer or cross-drug TACSTD2 / CLDN4 ICI biomarker from
  this template alone.
- Do not treat a pooled Hedges' *g* as a treatment-effect HR, or mix SMD with
  log-HR in one forest.
- Do not treat trim-and-fill, PET-PEESE, or a selection-model MLE as the
  "bias-corrected" truth.
- Do not download or analyse EGA / dbGaP files under this playbook.
- Do not impute effects for panel-missing series (GSE136961, GSE93157).
- Do not read a cohort-level meta-regression slope as a patient-level
  interaction (aggregation bias).
- Do not hide the missing-trial envelope. If the tipping-point \(\delta\) is
  small, the open-GEO result is fragile and must be said to be fragile.
- The `--demo` numbers are simulated. Copying them into a results table is
  fabricating data.

**中文.**

- 不能单凭本模板声称跨癌种 / 跨药物生物标志物。
- 不能把 Hedges' *g* 当成治疗 HR，不能把 SMD 与 log-HR 画进同一张森林图。
- 不能把 trim-and-fill / PET-PEESE / 选择模型 MLE 当成「校正后的真值」。
- 本手册下不下载、不分析 EGA / dbGaP。
- 面板缺基因的系列不补值。
- 队列水平元回归不是患者水平交互。
- 必须报告缺失试验包络；tipping point 很小就要写「脆弱」。
- `--demo` 数字是合成的，写入结果表等于造假。

---

## 11. References · 参考文献

- Gelman A, Hill J. *Data Analysis Using Regression and Multilevel/Hierarchical
  Models.* Cambridge University Press, 2007. (partial pooling)
- Raudenbush SW, Bryk AS. *Hierarchical Linear Models.* 2nd ed. Sage, 2002.
- DerSimonian R, Laird N. Meta-analysis in clinical trials. *Control Clin Trials.*
  1986;7:177-188.
- Paule RC, Mandel J. Consensus values and weighting factors. *J Res NBS.*
  1982;87:377-385.
- Viechtbauer W. Conducting meta-analyses in R with the metafor package.
  *J Stat Softw.* 2010;36(3):1-48.
- Hartung J, Knapp G. A refined method for the meta-analysis of controlled
  clinical trials with binary endpoints. *Stat Med.* 2001;20:3875-3889.
- Higgins JPT, Thompson SG, Spiegelhalter DJ. A re-evaluation of
  random-effects meta-analysis. *J R Stat Soc A.* 2009;172:137-159.
- Egger M, et al. Bias in meta-analysis detected by a simple, graphical test.
  *BMJ.* 1997;315:629-634.
- Sterne JAC, et al. Recommendations for examining and interpreting funnel
  plot asymmetry. *BMJ.* 2011;343:d4002.
- Duval S, Tweedie R. Trim and fill. *Biometrics.* 2000;56:455-463.
- Stanley TD, Doucouliagos H. Meta-regression approximations to reduce
  publication-selection bias. *Res Synth Methods.* 2014;5:60-78.
- Vevea JL, Hedges LV. A general linear model for estimating effect size in
  the presence of publication bias. *Psychometrika.* 1995;60:419-435.
- Copas J, Shi JQ. Meta-analysis, funnel plots and sensitivity analysis.
  *Biostatistics.* 2000;1:247-262.
- Manski CF. *Partial Identification of Probability Distributions.* Springer, 2003.
- Turner RM, et al. Predicting the extent of heterogeneity in meta-analysis,
  using empirical data from the Cochrane Database of Systematic Reviews.
  *Int J Epidemiol.* 2012;41:818-827.
- Rhodes KM, et al. Predictive distributions were developed for the extent of
  heterogeneity in meta-analyses of continuous outcome data. *J Clin Epidemiol.*
  2015;68:52-60.
- Hedges LV. Distribution theory for Glass's estimator of effect size and
  related estimators. *J Educ Stat.* 1981;6:107-128.

Repo-internal catalogs used to name open vs controlled resources (no effect
sizes copied): `notes/ici_catalog.tsv` (PR #5),
`results/gpt_checkmate/ega_catalog.tsv` (PR #12),
`example/known_open_ici_cohorts.csv`.
