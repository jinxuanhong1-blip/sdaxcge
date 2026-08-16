# TACSTD2 / CLDN4 与 ICI 疗效：开放 NSCLC 队列 Meta 分析

## 中文

### 结论

在可估计的 3 个开放 NSCLC 队列中，没有证据表明治疗前 **TACSTD2** 或 **CLDN4** 转录本表达与 ICI 应答稳定相关。效应值定义为“应答者减非应答者”的 Hedges' g；正值表示应答者表达更高。

| 基因 | 固定效应 g（95% CI） | P | 随机效应 g（95% CI） | P | I² |
|---|---:|---:|---:|---:|---:|
| TACSTD2 | -0.046 (-0.587, 0.495) | 0.867 | -0.050 (-0.610, 0.509) | 0.860 | 6.2% |
| CLDN4 | -0.120 (-0.662, 0.422) | 0.664 | -0.134 (-0.732, 0.465) | 0.661 | 17.2% |

两种模型的置信区间都跨越 0，而且效应估计接近 0。结果不支持把这两个基因单独作为跨队列 ICI 应答标志物；但小样本量和终点不完全一致意味着也不能据此排除较小效应。

### 队列与终点

- **GSE126044**：16 例治疗前 NSCLC；GEO 直接标注 responder 5 例、non-responder 11 例。
- **GSE135222**：27 例治疗前 NSCLC；GEO 仅开放 PFS 状态和时间，没有直接 RECIST 标签。本分析按已发表的 6 个月持久临床获益规则，以 PFS ≥180 天为应答（7 例），其余为非应答（20 例）。这是一个预先明确但近似的应答终点。
- **GSE136961**：21 例、9 DCB / 12 NDB；其 Oncomine 395 基因靶向面板的开放 TPM 矩阵中不含 TACSTD2 或 CLDN4，因此两个效应都不可估计，未进入合并估计。没有补值。
- **GSE166449**：22 例治疗前 NSCLC；样本标题直接标注 responder 7 例、non-responder 15 例。

各队列的效应、样本量、均值、置信区间和 P 值见 `results/gpt_meta/cohort_effects.tsv`。GSE126044 中两个基因方向均为应答者较低，但置信区间跨 0；另外两个可估计队列的方向较弱且相反。因此，合并结果接近零不是由一致的强效应产生。

### 方法

1. 从 NCBI GEO FTP 下载 family SOFT 元数据和作者提交的表达矩阵；URL、字节数及 SHA-256 记录于 `results/gpt_meta/source_manifest.tsv`。
2. GSE126044 原始计数转换为 log2(CPM+1)；GSE135222 TPM 转换为 log2(TPM+1)；GSE166449 使用提交矩阵中的既有对数尺度值。
3. 每个队列、每个基因计算应答者与非应答者间的 Hedges' g（小样本校正标准化均差）及近似 95% CI。
4. 使用逆方差固定效应模型及 DerSimonian–Laird 随机效应模型分别合并。随机效应的 `tau²` 使用非负截断。
5. 全流程仅用 Python 标准库，可运行：

```bash
python3 scripts/gpt_meta/meta_analyze.py
```

### 限制

- 实际可合并的只有 3 个小队列（每个基因总计 65 例），随机效应方差和异质性估计不稳定。
- GSE135222 的开放 GEO 元数据没有直接 RECIST 最佳疗效，PFS 180 天二分法是持久获益代理终点；这与其他队列的直接应答/DCB 标签不完全相同。
- 平台涵盖全转录组 RNA-seq 和靶向面板；标准化均差减少了量纲差异，但不能消除平台、组织处理和患者构成差异。
- 分析是未调整的单基因比较，没有控制组织学、肿瘤纯度、PD-L1、治疗药物或其他混杂因素。
- 每基因只有 3 个可估计队列，不适合可靠评估发表偏倚。

## English

### Conclusion

Across the three estimable open NSCLC cohorts, there was no evidence that pretreatment **TACSTD2** or **CLDN4** transcript expression was consistently associated with ICI response. Effects are Hedges' g for responders minus nonresponders; positive values indicate higher expression in responders.

| Gene | Fixed-effect g (95% CI) | P | Random-effect g (95% CI) | P | I² |
|---|---:|---:|---:|---:|---:|
| TACSTD2 | -0.046 (-0.587, 0.495) | 0.867 | -0.050 (-0.610, 0.509) | 0.860 | 6.2% |
| CLDN4 | -0.120 (-0.662, 0.422) | 0.664 | -0.134 (-0.732, 0.465) | 0.661 | 17.2% |

Both models produced estimates near zero with confidence intervals crossing zero. These data do not support either gene alone as a cross-cohort ICI-response marker. The small samples and non-identical endpoints nevertheless leave modest effects unresolved.

### Cohorts and outcomes

- **GSE126044:** 16 pretreatment NSCLC samples; GEO directly labels 5 responders and 11 nonresponders.
- **GSE135222:** 27 pretreatment NSCLC samples. GEO exposes PFS status/time but no direct RECIST label. Following the published six-month durable-benefit rule, this analysis classifies PFS ≥180 days as response (7) and the remainder as nonresponse (20). This is an explicit but approximate response endpoint.
- **GSE136961:** 21 samples, 9 DCB / 12 NDB. Neither target gene is present in the deposited Oncomine 395-gene TPM matrix, so neither effect is estimable or pooled. No values were imputed.
- **GSE166449:** 22 pretreatment NSCLC samples; GEO sample titles directly identify 7 responders and 15 nonresponders.

Per-cohort effects, group sizes, means, confidence intervals, and P values are in `results/gpt_meta/cohort_effects.tsv`. Both genes were lower in responders in GSE126044, although the confidence intervals crossed zero; the other two estimable cohorts had weaker effects in the opposite direction. The near-null pooled estimates therefore do not represent a set of concordant strong effects.

### Methods

1. Family SOFT metadata and submitter-provided expression matrices were downloaded from NCBI GEO FTP. URLs, byte counts, and SHA-256 checksums are in `results/gpt_meta/source_manifest.tsv`.
2. GSE126044 counts were transformed to log2(CPM+1), GSE135222 TPM to log2(TPM+1), and GSE166449 was analyzed on the deposited log-scale values.
3. For each cohort and gene, Hedges' g (small-sample-corrected standardized mean difference) and its approximate 95% CI compared responders with nonresponders.
4. Effects were pooled with inverse-variance fixed-effect and DerSimonian–Laird random-effects models. Random-effects `tau²` was truncated at zero.
5. The complete pipeline uses only the Python standard library:

```bash
python3 scripts/gpt_meta/meta_analyze.py
```

### Limitations

- Only three small cohorts were estimable (65 samples per gene in total), making heterogeneity and random-effects variance estimates imprecise.
- GSE135222 lacks an open direct RECIST best-response field; dichotomized 180-day PFS is a durable-benefit proxy and is not identical to the direct response/DCB labels in the other cohorts.
- Platforms include whole-transcriptome RNA-seq and a targeted panel. Standardized mean differences reduce scale incompatibility but do not remove platform, tissue-processing, or case-mix differences.
- This is an unadjusted single-gene analysis and does not control for histology, tumor purity, PD-L1, drug, or other confounders.
- Three estimable cohorts per gene are insufficient for meaningful publication-bias assessment.

### Reproducible outputs

- `results/gpt_meta/sample_expression.tsv`: analysis-scale values and labels
- `results/gpt_meta/cohort_effects.tsv`: cohort-level effects and non-estimable rows
- `results/gpt_meta/meta_summary.tsv`: fixed/random pooled estimates
- `results/gpt_meta/forest_plot.svg`: forest plot
- `results/gpt_meta/source_manifest.tsv`: source provenance and checksums
