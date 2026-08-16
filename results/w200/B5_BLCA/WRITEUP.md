# B5_BLCA — CLDN4-high and ICI objective response in urothelial carcinoma

**Verdict: the claim is NOT supported. Conclusion label: `inconsistent`.**

The urothelial arm of claim B5 asserts that CLDN4-high tumours respond to immune
checkpoint inhibition with an odds ratio of **0.42**. In three public urothelial
ICI cohorts (n = 406 response-evaluable patients, 91 objective responders) the
pooled odds ratio is **1.31 (95% CI 0.82–2.10)** — the *opposite* direction, and
statistically incompatible with 0.42 (z = 4.74, p = 2.1 × 10⁻⁶).

We did not tune any specification toward 0.42 (`did_we_tune_to_0.42: false`).

---

## 1. Question and design

**Question.** In patients with urothelial/bladder carcinoma treated with anti-PD-1
or anti-PD-L1 therapy, do tumours with high CLDN4 expression have lower odds of
RECIST objective response than tumours with low CLDN4 expression, and is the
effect size 0.42?

**Design.** Retrospective association study across three independent public
cohorts, with a single prespecified primary test per cohort and one random-effects
pooled estimate. The analysis plan (`analysis_plan.md`) was committed to git
**before** any CLDN4-versus-response statistic was computed; the commit ordering is
auditable in the branch history.

**Status of the claim.** No published CLDN4 × ICI meta-analysis could be located
(PubMed `"claudin 4" AND ("checkpoint inhibitor" OR "anti-PD-1")` returns 0
records). OR = 0.42 is therefore treated as an **unsourced assertion under test**,
not as a published result being reproduced. Claim B5 refers to 11 cohorts; only
**3** eligible public urothelial ICI cohorts with both RNA and RECIST exist as far
as we could establish, and we do not name cohorts we cannot enumerate.

## 2. Data provenance and accession verification

| Cohort | Source | Treatment | Expression scale |
|---|---|---|---|
| IMvigor210 | `IMvigor210CoreBiologies` `cds.RData` (Mariathasan 2018, Nature, doi:10.1038/nature25501) | atezolizumab | counts → log2(TPM+1) |
| BACI | GEO **GSE176307** (Robertson 2021) | atezolizumab / pembrolizumab / other ICI | salmon TPM → log2(TPM+1) |
| Snyder | bhklab PredictIO, **Zenodo 7058399** (Snyder 2017, PLoS Med; harmonised by Bareche 2022, PMID 36055464) | atezolizumab | log2(TPM+0.001) as deposited |

Every file is recorded with URL, byte count and SHA-256 in `download_manifest.json`
and `catalog.tsv`.

**IMvigor210 required special handling and we state it openly.** The official
distribution point `http://research-pub.gene.com/IMvigor210CoreBiologies/` now
returns **HTTP 404**, so `cds.RData` was taken from a third-party GitHub mirror. A
mirror is not self-validating, so it was accepted only after passing three checks
(`validation.json`):

1. Sample count is **348**, matching the published RNA-seq cohort.
2. The RECIST distribution reproduces the published values **exactly**: CR 25,
   PR 43, SD 63, PD 167, non-evaluable 50.
3. The immune-phenotype distribution reproduces the published values **exactly**:
   desert 76, excluded 134, inflamed 74, unknown 64.

Additionally, an **independently harmonised copy of the same trial** (bhklab
PredictIO, processed from raw data by a different group with a different pipeline)
gives an identical RECIST tally, and its CLDN4 values agree with ours at
**Spearman ρ = 0.991** across all 348 samples. The mirror is therefore accepted as
faithful.

CLDN4 was resolved as **ENSG00000189143** where Ensembl identifiers were available
and by exact symbol match otherwise; it matched exactly one row in every cohort.

## 3. Sample size

| Cohort | Available | CLDN4 present | Response evaluable (analysed) | Responders | Non-responders | ORR | Excluded, non-evaluable |
|---|---|---|---|---|---|---|---|
| IMvigor210 | 348 | 348 | **298** | 68 | 230 | 22.8 % | 50 |
| BACI | 89 | 89 | **87** | 16 | 71 | 18.4 % | 2 |
| Snyder | 25 | 25 | **21** | 7 | 14 | 33.3 % | 4 |
| **Total pooled** | 462 | 462 | **406** | **91** | 315 | 22.4 % | 56 |

`n` is at the patient level throughout. 2 BACI patients had 2 RNA-seq runs each;
replicate runs were averaged within patient and never counted as independent `n`. The PredictIO copy of IMvigor210 (348 samples) is **excluded from pooling** to
prevent double counting and is used only for the concordance check in §6 S7.

## 4. Methods and reproducibility

**Endpoint.** ORR from RECIST best overall response: responder = {CR, PR},
non-responder = {SD, PD}. Non-evaluable records were excluded from the primary
analysis and counted in §3.

**Exposure.** CLDN4-high = expression strictly above the cohort-specific median,
computed within the analysis set and **without reference to response**; ties go to
the low group. No filter anywhere in the pipeline references CLDN4 or response.

**Tests.** Two-sided Fisher exact test per cohort; odds ratio of *response* for
CLDN4-high versus CLDN4-low as the conditional maximum-likelihood estimate with an
exact 95 % CI, so **OR < 1 is the direction the claim asserts**. Pooling used Woolf
log ORs (Haldane–Anscombe 0.5 correction only if a cell were zero; no cell was) in
a DerSimonian–Laird random-effects model as the prespecified primary, with
Hartung–Knapp limits because k = 3, plus fixed-effect inverse-variance and
Mantel–Haenszel.

**Reproducibility.** Seed 20260816 in `repro/seed.txt`; environment in
`repro/pip-freeze.txt`; `scripts/w200/B5_BLCA/{download,prepare_cohorts,analyze}.py`
plus `extract_imvigor210.R` regenerate every number and figure from the raw
downloads. Expression scale is declared per cohort and carried in an
`expression_scale` column of `per_sample_expression.csv`.

## 5. Multiple testing

The primary analysis is **one** prespecified pooled test, so no correction applies
to it. The exploratory gene family (CLDN1, CLDN2, CLDN3, CLDN7, CLDN18, TACSTD2,
EPCAM, CDH1, OCLN, TJP1, CD274) is **m = 11** tests corrected by
**Benjamini–Hochberg** at threshold **q < 0.05**; **0 of 11** genes pass, and
adjusted values are reported as q-values, not as "FDR". CLDN4 is the primary
hypothesis and is deliberately excluded from that family. The sensitivity analyses
in §6 are robustness checks on a single hypothesis, not independent confirmatory
tests, and no alpha is claimed for them.

## 6. Results

### Primary

| Cohort | CLDN4-high ORR | CLDN4-low ORR | OR (95 % CI) | Fisher p |
|---|---|---|---|---|
| IMvigor210 | 26.2 % (39/149) | 19.5 % (29/149) | 1.47 (0.82–2.64) | 0.214 |
| BACI | 18.6 % (8/43) | 18.2 % (8/44) | 1.03 (0.30–3.53) | 1.000 |
| Snyder | 30.0 % (3/10) | 36.4 % (4/11) | 0.76 (0.08–6.52) | 1.000 |

**Pooled (n = 406).** All four pooling methods agree:

| Model | OR | 95 % CI | p |
|---|---|---|---|
| Random effects (DerSimonian–Laird, primary) | **1.31** | 0.82–2.10 | 0.258 |
| Random effects, Hartung–Knapp | 1.31 | 0.71–2.43 | 0.199 |
| Fixed effect, inverse variance | 1.31 | 0.82–2.10 | 0.258 |
| Mantel–Haenszel | 1.31 | 0.82–2.10 | — |

Heterogeneity is absent: Q = 0.71 on 2 df, p = 0.700, I² = 0 %, τ² = 0.

**Comparison with the claimed value.** The claim is not merely unconfirmed; it is
excluded. 0.42 lies **outside** the pooled 95 % CI, and formally testing
H₀: OR = 0.42 gives **z = 4.74, p = 2.1 × 10⁻⁶**. Only 1 of 3 cohorts even points
in the claimed direction, and that one (Snyder, n = 21) has an interval spanning
0.08–6.52. No cohort-level test reached p < 0.05. Label: **`DOES_NOT_MATCH`**.

**This is not simply a power failure.** Two independent lines of evidence:

- *Simulation power.* Had the true OR been 0.42, this design would have detected it
  with probability **0.70** at α = 0.05 (20 000 simulations at the observed group
  sizes and CLDN4-low response rates). Power is moderate rather than high, so a
  small true effect remains possible — but the observed point estimate is on the
  *other side of 1*, and 0.42 is directly rejected.
- *Positive control.* The same median-split Fisher machinery, on the same patients,
  recovers the CD8 T-effector/IFN-γ association that is published for IMvigor210:
  pooled **OR 1.76 (1.09–2.83), p = 0.021**, and continuously in IMvigor210
  **OR 1.46 per SD, p = 0.0069**. The pipeline detects a real effect of this
  magnitude when one is present.
- *Negative control.* Housekeeping genes are null as expected: pooled
  **OR 0.88 (0.55–1.41), p = 0.592**. The pipeline does not manufacture signal.

### Sensitivity — every specification agrees

| Analysis | Result |
|---|---|
| S1 cutoff: median / tertile / quartile / upper-quartile / 60-40 | pooled OR 1.31 / 1.52 / 1.25 / 1.47 / 1.33 — never near 0.42 |
| S2 continuous, OR per SD | IMvigor210 1.15 (0.87–1.52) p = 0.338; BACI 0.87 (0.52–1.47) p = 0.608; Snyder 0.96 (0.38–2.43) p = 0.930 |
| S3 Mann–Whitney of CLDN4 by response | p = 0.239 / 0.865 / 0.799 |
| S4 ITT, non-evaluable counted as non-responders | pooled OR 1.30 (0.82–2.06) p = 0.270 |
| S5 IMvigor210 bladder-site biopsies only (n = 168) | OR 1.14 (0.53–2.43) p = 0.859 |
| S6 DESeq size-factor normalisation instead of TPM | OR 1.47 (0.82–2.64) p = 0.214 |
| S7 independent processing pipeline (PredictIO) | CLDN4 ρ = 0.991 (n = 348); OR 1.47, identical conclusion |
| S8 multivariable adjustment (n = 205) | dichotomous adjusted OR 1.03 (0.50–2.13) p = 0.938; continuous adjusted OR 0.88 (0.67–1.16) p = 0.352 |

The continuous estimates straddle 1 in both directions while the dichotomised ones
sit slightly above 1; the honest reading of that pattern is **no association of any
consistent direction**, not a positive association.

### Exploratory (not tests of the claim)

No claudin, junction or checkpoint gene in the m = 11 family passes BH q < 0.05.

One finding deserves explicit mention because it is the most claim-favourable
result in this slice, and it still does not rescue the claim. In IMvigor210, CLDN4
is **weakly negatively correlated with the CD8 T-effector score (ρ = −0.133,
p = 0.013, n = 348)**, and median CLDN4 decreases monotonically across immune
phenotypes — desert 4.67, excluded 4.47, inflamed 4.17 (Kruskal–Wallis H = 3.56,
p = 0.169, not significant). So a faint immune-exclusion correlate in the expected
direction is present at the level of the *mechanistic proxy*, yet it does **not**
translate into reduced objective response at the level of the *clinical endpoint*.
Under this project's evidence grading that is Tier C mechanistic consistency at
best; it is not evidence for OR = 0.42.

## 7. Limitations

- **Three cohorts, not eleven.** The urothelial arm alone cannot confirm or refute a
  pooled 11-cohort number. It can only say that the urothelial evidence is
  incompatible with OR = 0.42.
- **Power is moderate (0.70 for OR = 0.42), not high.** A modest true effect cannot
  be excluded; what is excluded is an effect of the claimed magnitude.
- **Snyder contributes 21 evaluable patients**, so its interval (0.08–6.52) is
  nearly uninformative on its own.
- **BACI is real-world**, with mixed ICI agents and investigator-assessed rather
  than centrally reviewed response.
- **IMvigor210 came from a mirror** because the official source now 404s; accepted
  only after exact reproduction of published RECIST and immune-phenotype
  distributions and ρ = 0.991 agreement with an independent harmonisation.
- **Expression scales differ across cohorts.** The within-cohort median split makes
  the contrast scale-free but does not harmonise platforms; a cross-cohort absolute
  CLDN4 threshold is not established here.
- **Association only.** No causal, mechanistic or predictive-utility claim is
  implied, and no patient-selection recommendation follows.
- **Biopsy site varies** within IMvigor210 (195 bladder, 67 kidney, 26 ureter, 26
  lymph node, and others); S5 restricts to bladder and is unchanged.

## 8. Deviations from the prespecified plan

Two clarifications, neither affecting the primary result:

1. The plan said "cohort-specific median" without stating the reference set. We
   computed it within the **analysis set** (response-evaluable patients with CLDN4),
   which keeps the split response-blind and the groups balanced. Computing it over
   all RNA samples changes no conclusion.
2. The plan did not prespecify a formal test against 0.42 or a power simulation.
   Both were added because "we failed to reject 1" is a weaker and less honest
   statement than "0.42 is rejected at p = 2 × 10⁻⁶ with 70 % power". Both are
   reported with their assumptions.

Nothing was removed from the plan, and no cohort, cutoff or covariate was changed
after seeing a result.

## 9. Post-hoc open-cohort expansion (requested after the primary was locked)

Two further **open** BLCA ICI RNA+CLDN4 cohorts were added after commit `b1e1159`,
because the follow-up request was to report OR/n/p for every open BLCA ICI dataset
with CLDN4. They are **not** part of the locked primary and were not used to chase
0.42.

| Cohort | Source | n | responders | OR (95% CI) | p | Endpoint |
|---|---|---|---|---|---|---|
| UC-GENOME | cBioPortal `blca_bcan_hcrn_2022` (Damrauer 2022, PMID 36333289) | 89 | 34 | 0.71 (0.27–1.82) | 0.515 | IO best response CR/PR vs SD/PD |
| GSE111636 | GEO GSE111636 (pembrolizumab, HTA-2.0) | 11 | 6 | 6.44 (0.33–490) | 0.242 | depositor binary responder/progressor |

UC-GENOME is independent of BACI/GSE176307: the UC-GENOME paper used GSE176307 as
an external validation set, and the patient identifiers do not overlap.

Post-hoc 4-cohort pool (prespecified 3 + UC-GENOME; GSE111636 held out because
its endpoint is not RECIST): **OR 1.14 (0.75–1.72), p = 0.540, I² = 0%**. Testing
H₀: OR = 0.42 still rejects at **z = 4.73, p = 2.3 × 10⁻⁶**. Adding the only
open cohort that points in the claimed direction moves the estimate from 1.31
toward 1, not toward 0.42.

Full OR/n/p table: `or_n_p.md` / `or_n_p.csv`.

---

# 中文 / Chinese

# B5_BLCA — 尿路上皮癌中 CLDN4 高表达与 ICI 客观缓解

**结论：该主张不成立。结论标签：`inconsistent`（不一致）。**

主张 B5 的尿路上皮癌部分断言，CLDN4 高表达肿瘤对免疫检查点抑制治疗的缓解比值比为
**0.42**。在三个公开的尿路上皮癌 ICI 队列中（n = 406 例可评估缓解的患者，91 例客观
缓解），合并比值比为 **1.31（95% CI 0.82–2.10）**，方向*相反*，且与 0.42 在统计学
上不相容（z = 4.74，p = 2.1 × 10⁻⁶）。

我们没有为接近 0.42 而调整任何分析设定（`did_we_tune_to_0.42: false`）。

## 1. 问题与设计

**问题.** 在接受抗 PD-1 或抗 PD-L1 治疗的尿路上皮癌/膀胱癌患者中，CLDN4 高表达肿瘤
获得 RECIST 客观缓解的比值是否低于 CLDN4 低表达肿瘤，效应量是否为 0.42？

**设计.** 跨三个独立公开队列的回顾性关联研究，每个队列一个预先设定的主要检验，并给出
一个随机效应合并估计。分析方案（`analysis_plan.md`）在计算任何 CLDN4 与缓解相关的
统计量**之前**已提交到 git；提交顺序可在分支历史中审计。

**主张的来源状态.** 未能检索到任何已发表的 CLDN4 × ICI 荟萃分析（PubMed
`"claudin 4" AND ("checkpoint inhibitor" OR "anti-PD-1")` 返回 0 条记录）。因此
OR = 0.42 被视为**无来源的待检验断言**，而非正在复现的已发表结果。主张 B5 提到 11 个
队列；据我们所能确认，同时具备 RNA 与 RECIST 的合格公开尿路上皮癌 ICI 队列只有 **3**
个，我们不会命名无法列举的队列。

## 2. 数据来源与编号核验

| 队列 | 来源 | 治疗 | 表达尺度 |
|---|---|---|---|
| IMvigor210 | `IMvigor210CoreBiologies` `cds.RData`（Mariathasan 2018, Nature, doi:10.1038/nature25501） | atezolizumab | counts → log2(TPM+1) |
| BACI | GEO **GSE176307**（Robertson 2021） | atezolizumab / pembrolizumab / 其他 ICI | salmon TPM → log2(TPM+1) |
| Snyder | bhklab PredictIO，**Zenodo 7058399**（Snyder 2017, PLoS Med；由 Bareche 2022 统一处理，PMID 36055464） | atezolizumab | 存档为 log2(TPM+0.001) |

每个文件的 URL、字节数与 SHA-256 均记录于 `download_manifest.json` 与 `catalog.tsv`。

**IMvigor210 需要特殊处理，我们公开说明。** 官方分发地址
`http://research-pub.gene.com/IMvigor210CoreBiologies/` 现返回 **HTTP 404**，因此
`cds.RData` 取自第三方 GitHub 镜像。镜像本身不能自证，故仅在通过以下三项检查后采用
（`validation.json`）：

1. 样本数为 **348**，与已发表的 RNA-seq 队列一致。
2. RECIST 分布**完全**复现已发表数值：CR 25、PR 43、SD 63、PD 167、不可评估 50。
3. 免疫表型分布**完全**复现已发表数值：desert 76、excluded 134、inflamed 74、
   未知 64。

此外，同一试验的**独立统一处理副本**（bhklab PredictIO，由另一团队用不同流程从原始
数据处理）给出完全相同的 RECIST 计数，其 CLDN4 值与我们的结果在全部 348 个样本上达到
**Spearman ρ = 0.991**。因此认定该镜像忠实可信。

CLDN4 在有 Ensembl 编号时按 **ENSG00000189143** 解析，否则按符号精确匹配；在每个队列
中均恰好匹配一行。

## 3. 样本量

| 队列 | 可用 | 有 CLDN4 | 缓解可评估（纳入分析） | 缓解者 | 非缓解者 | ORR | 因不可评估而排除 |
|---|---|---|---|---|---|---|---|
| IMvigor210 | 348 | 348 | **298** | 68 | 230 | 22.8 % | 50 |
| BACI | 89 | 89 | **87** | 16 | 71 | 18.4 % | 2 |
| Snyder | 25 | 25 | **21** | 7 | 14 | 33.3 % | 4 |
| **合并总计** | 462 | 462 | **406** | **91** | 315 | 22.4 % | 56 |

全文 `n` 均为患者层面。BACI 中有 2 例患者各有 2 次 RNA-seq 检测；重复检测在患者内取
平均，绝不作为独立 `n` 计数。IMvigor210 的 PredictIO 副本（348 个样本）**不参与合并**，
以避免重复计数，仅用于 §6 S7 的一致性检查。

## 4. 方法与可重复性

**终点.** 由 RECIST 最佳总体缓解定义 ORR：缓解者 = {CR, PR}，非缓解者 = {SD, PD}。
不可评估记录从主要分析中排除，并在 §3 中计数。

**暴露.** CLDN4 高 = 表达严格高于队列特异中位数，在分析集内计算且**不参考缓解信息**；
等于中位数者归入低表达组。流程中任何过滤步骤均不引用 CLDN4 或缓解。

**检验.** 每队列采用双侧 Fisher 精确检验；CLDN4 高对低的*缓解*比值比采用条件极大似然
估计及精确 95 % CI，因此 **OR < 1 即为该主张所断言的方向**。合并采用 Woolf 对数比值比
（仅在出现零格时施加 Haldane–Anscombe 0.5 校正；本研究无零格），以 DerSimonian–Laird
随机效应模型作为预设主要方法；因 k = 3 另报 Hartung–Knapp 界限，并报固定效应逆方差与
Mantel–Haenszel。

**可重复性.** 种子 20260816 见 `repro/seed.txt`；环境见 `repro/pip-freeze.txt`；
`scripts/w200/B5_BLCA/{download,prepare_cohorts,analyze}.py` 与
`extract_imvigor210.R` 可从原始下载重新生成全部数值与图形。表达尺度按队列声明，并写入
`per_sample_expression.csv` 的 `expression_scale` 列。

## 5. 多重检验

主要分析为**一个**预设的合并检验，因此不施加校正。探索性基因族（CLDN1、CLDN2、CLDN3、
CLDN7、CLDN18、TACSTD2、EPCAM、CDH1、OCLN、TJP1、CD274）共 **m = 11** 个检验，采用
**Benjamini–Hochberg** 校正，阈值 **q < 0.05**；**11 个中有 0 个**通过，校正后数值报作
q 值，而非 "FDR"。CLDN4 为主要假设，特意不纳入该基因族。§6 的敏感性分析是对单一假设的
稳健性检查，而非独立的验证性检验，不为其主张任何 alpha。

## 6. 结果

### 主要结果

| 队列 | CLDN4 高 ORR | CLDN4 低 ORR | OR（95 % CI） | Fisher p |
|---|---|---|---|---|
| IMvigor210 | 26.2 %（39/149） | 19.5 %（29/149） | 1.47（0.82–2.64） | 0.214 |
| BACI | 18.6 %（8/43） | 18.2 %（8/44） | 1.03（0.30–3.53） | 1.000 |
| Snyder | 30.0 %（3/10） | 36.4 %（4/11） | 0.76（0.08–6.52） | 1.000 |

**合并（n = 406）.** 四种合并方法结果一致：

| 模型 | OR | 95 % CI | p |
|---|---|---|---|
| 随机效应（DerSimonian–Laird，主要） | **1.31** | 0.82–2.10 | 0.258 |
| 随机效应，Hartung–Knapp | 1.31 | 0.71–2.43 | 0.199 |
| 固定效应，逆方差 | 1.31 | 0.82–2.10 | 0.258 |
| Mantel–Haenszel | 1.31 | 0.82–2.10 | — |

无异质性：Q = 0.71（df = 2），p = 0.700，I² = 0 %，τ² = 0。

**与所称数值的比较.** 该主张不仅未被证实，而且被排除。0.42 落在合并 95 % CI **之外**，
对 H₀: OR = 0.42 进行正式检验得 **z = 4.74，p = 2.1 × 10⁻⁶**。3 个队列中仅 1 个指向所
称方向，而该队列（Snyder，n = 21）的区间跨度为 0.08–6.52。没有任何队列层面的检验达到
p < 0.05。标签：**`DOES_NOT_MATCH`**。

**这并非单纯的检验效能不足.** 有两条独立证据：

- *模拟效能.* 若真实 OR 为 0.42，本设计在 α = 0.05 下检出的概率为 **0.70**（按观察到的
  分组规模与 CLDN4 低表达组缓解率进行 20 000 次模拟）。效能为中等而非很高，因此不能排除
  较小的真实效应；但观察到的点估计位于 1 的*另一侧*，且 0.42 被直接拒绝。
- *阳性对照.* 同一套中位数分割 Fisher 流程，在同一批患者上，重现了 IMvigor210 已发表的
  CD8 T 效应/IFN-γ 关联：合并 **OR 1.76（1.09–2.83），p = 0.021**；在 IMvigor210 中按
  连续变量为 **每 SD OR 1.46，p = 0.0069**。当效应真实存在时，本流程能够检出该量级的
  效应。
- *阴性对照.* 管家基因如预期为无效应：合并 **OR 0.88（0.55–1.41），p = 0.592**。本流程
  不会制造出信号。

### 敏感性分析 — 所有设定结果一致

| 分析 | 结果 |
|---|---|
| S1 切点：中位数 / 三分位 / 四分位 / 上四分位 / 60-40 | 合并 OR 1.31 / 1.52 / 1.25 / 1.47 / 1.33 — 从未接近 0.42 |
| S2 连续变量，每 SD 的 OR | IMvigor210 1.15（0.87–1.52）p = 0.338；BACI 0.87（0.52–1.47）p = 0.608；Snyder 0.96（0.38–2.43）p = 0.930 |
| S3 按缓解分组的 CLDN4 Mann–Whitney 检验 | p = 0.239 / 0.865 / 0.799 |
| S4 ITT，不可评估计为非缓解 | 合并 OR 1.30（0.82–2.06）p = 0.270 |
| S5 仅 IMvigor210 膀胱部位活检（n = 168） | OR 1.14（0.53–2.43）p = 0.859 |
| S6 以 DESeq 大小因子归一化替代 TPM | OR 1.47（0.82–2.64）p = 0.214 |
| S7 独立处理流程（PredictIO） | CLDN4 ρ = 0.991（n = 348）；OR 1.47，结论相同 |
| S8 多因素校正（n = 205） | 二分类校正后 OR 1.03（0.50–2.13）p = 0.938；连续校正后 OR 0.88（0.67–1.16）p = 0.352 |

连续估计在 1 的两侧摆动，而二分类估计略高于 1；对该模式的诚实解读是**不存在任何方向
一致的关联**，而非存在正向关联。

### 探索性分析（并非对该主张的检验）

m = 11 的基因族中没有任何 claudin、连接或检查点基因通过 BH q < 0.05。

有一项结果值得明确指出，因为它是本切片中最有利于该主张的结果，但仍不足以挽救该主张。
在 IMvigor210 中，CLDN4 与 CD8 T 效应评分呈**弱负相关（ρ = −0.133，p = 0.013，
n = 348）**，且 CLDN4 中位数在免疫表型间单调下降——desert 4.67、excluded 4.47、
inflamed 4.17（Kruskal–Wallis H = 3.56，p = 0.169，不显著）。因此在*机制替代指标*层面
存在方向符合预期的微弱免疫排斥相关性，但它**未**转化为*临床终点*层面客观缓解的降低。
按本项目的证据分级，这最多属于 C 级机制一致性；它不构成 OR = 0.42 的证据。

## 7. 局限性

- **三个队列，而非十一个.** 单凭尿路上皮癌部分无法证实或否证一个 11 队列的合并数值，
  只能说明尿路上皮癌证据与 OR = 0.42 不相容。
- **效能为中等（对 OR = 0.42 为 0.70），并不高.** 不能排除较小的真实效应；被排除的是所称
  量级的效应。
- **Snyder 仅贡献 21 例可评估患者**，其区间（0.08–6.52）单独看几乎不提供信息。
- **BACI 为真实世界队列**，ICI 药物混杂，且缓解由研究者评估而非中心化复核。
- **IMvigor210 取自镜像**，因官方来源现已 404；仅在完全复现已发表 RECIST 与免疫表型分布
  并与独立统一处理结果达到 ρ = 0.991 一致后采用。
- **各队列表达尺度不同.** 队列内中位数分割使对比不依赖尺度，但并未统一平台；本研究未
  建立跨队列的 CLDN4 绝对阈值。
- **仅为关联.** 不暗示任何因果、机制或预测效用的主张，也不导出任何患者筛选建议。
- **IMvigor210 内活检部位不一**（膀胱 195、肾 67、输尿管 26、淋巴结 26 及其他）；S5 限定
  为膀胱部位，结论不变。

## 8. 与预设方案的偏离

两处澄清，均不影响主要结果：

1. 方案中写的是"队列特异中位数"，未说明参照集合。我们在**分析集**（有 CLDN4 且缓解可
   评估的患者）内计算，以保持分割对缓解盲法且分组均衡。若在全部 RNA 样本上计算，结论
   不变。
2. 方案未预设针对 0.42 的正式检验与效能模拟。两者均已加入，因为"未能拒绝 1"是比
   "0.42 在 70 % 效能下以 p = 2 × 10⁻⁶ 被拒绝"更弱、更不诚实的表述。两者均连同其假设
   一并报告。

方案中没有任何内容被删除，也没有在看到结果后更改任何队列、切点或协变量。

## 9. 事后开放队列扩展（主要分析锁定后应要求加入）

在提交 `b1e1159` 之后，又加入了两个**开放**的膀胱/尿路上皮癌 ICI RNA+CLDN4 队列，因后续
要求报告每一个开放队列的 OR/n/p。它们**不属于**已锁定的主要分析，也未被用于追逐 0.42。

| 队列 | 来源 | n | 缓解者 | OR（95% CI） | p | 终点 |
|---|---|---|---|---|---|---|
| UC-GENOME | cBioPortal `blca_bcan_hcrn_2022`（Damrauer 2022，PMID 36333289） | 89 | 34 | 0.71（0.27–1.82） | 0.515 | IO 最佳缓解 CR/PR 对 SD/PD |
| GSE111636 | GEO GSE111636（pembrolizumab，HTA-2.0） | 11 | 6 | 6.44（0.33–490） | 0.242 | 提交者二分的缓解/进展 |

UC-GENOME 与 BACI/GSE176307 相互独立：UC-GENOME 论文将 GSE176307 用作外部验证集，且
患者编号无重叠。

事后 4 队列合并（预设 3 个 + UC-GENOME；因终点不是 RECIST，GSE111636 未纳入）：
**OR 1.14（0.75–1.72），p = 0.540，I² = 0%**。对 H₀: OR = 0.42 仍以
**z = 4.73，p = 2.3 × 10⁻⁶** 拒绝。加入唯一指向所称方向的开放队列，只是把估计从 1.31
拉向 1，而不是拉向 0.42。

完整 OR/n/p 表见 `or_n_p.md` / `or_n_p.csv`。
