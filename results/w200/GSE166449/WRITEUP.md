# GSE166449: TACSTD2 and CLDN4 vs pembrolizumab response

# GSE166449：TACSTD2 与 CLDN4 对 pembrolizumab 应答

**Verdict / 结论:** this 22-patient single-arm LUAD pembrolizumab cohort provides
**no persuasive evidence** that pretreatment `TACSTD2`, `CLDN4`, their
equal-weight two-gene score, or a bulk T-cell proxy is associated with
response. `TACSTD2` is numerically **higher** in responders, opposite a simple
high-TROP2 resistance story, but the interval includes both directions.

**结论：** 这个 22 例单臂肺腺癌 pembrolizumab 队列**不能提供有说服力的证据**
表明治疗前 `TACSTD2`、`CLDN4`、等权双基因评分或 bulk T 细胞代理指标与应答
相关。`TACSTD2` 在应答者中数值上**更高**，与“高 TROP2 即耐药”的简单假说
相反，但置信区间同时包含两个方向。

---

## 1. Cohort / 队列

| Item | Value |
|------|--------|
| GEO | [GSE166449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166449) |
| Paper | Lee et al., *Cell* 2021, PMID [33857424](https://pubmed.ncbi.nlm.nih.gov/33857424/), DOI [10.1016/j.cell.2021.03.030](https://doi.org/10.1016/j.cell.2021.03.030) |
| Disease | advanced lung adenocarcinoma |
| Treatment | pembrolizumab (paper); GEO titles say only “immunotherapy” |
| Biopsy | pretreatment |
| Labels | 7 `Immunotherapy_Responder` vs 15 `Immunotherapy_nonResponder` |
| Response rule in the paper | RECIST 1.1; CR/PR = responder, SD/PD = non-responder |
| Public covariates | none (no PD-L1, line, survival, or patient-level RECIST category) |

GEO itself does not name the drug. The pembrolizumab assignment comes from the
Cell 2021 STAR Methods description of the Samsung Medical Center cohort.

GEO 本身未写明药物。pembrolizumab 来自 Cell 2021 STAR Methods 对三星医疗中心
队列的描述。

## 2. Methods (pre-specified) / 方法（预先指定）

- **Scale.** Use the deposited matrix as `log2(TPM + 1)`. Do not re-log.
- **Primary test.** Exact two-sided Mann–Whitney U for `TACSTD2` and `CLDN4`.
- **Effect size.** Mean difference with Welch 95% CI; Hedges’ g; AUC =
  P(R > NR); rank-biserial = 2·AUC − 1; univariable OR per 1 SD.
- **Secondary, descriptive.** Median-split Fisher exact OR (high vs low).
  This loses information and is not used to claim a cutoff.
- **Exploratory score.** Mean of within-cohort z-scored `TACSTD2` and `CLDN4`.
  Exact label-permutation p for the mean difference.
- **Context, not a hunt.** Bulk T-cell proxy = mean z of `CD3D`, `CD3E`,
  `CD8A`. Spearman vs each gene; same response tests as a power/context check.
- **Multiplicity.** BH q across the two primary gene tests only.
- **Not done.** Cutpoint search, multivariable fitting, gene-set fishing, or
  pooling with other GEO series.

未因结果更接近预期假说而改检验、改分组或改基因。

## 3. Results / 结果

### 3.1 Primary genes vs response

| Feature | Median R | Median NR | Mean R−NR (95% CI) | Exact MW p | BH q | AUC (bootstrap 95% CI) | High/low OR (95% CI) | Fisher p |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| TACSTD2 | 3.346 | 2.436 | 0.520 (-1.204, 2.244) | 0.407 | 0.814 | 0.619 (0.333, 0.867) | 3.750 (0.540, 26.046) | 0.361 |
| CLDN4 | 1.310 | 1.511 | 0.136 (-0.733, 1.006) | 0.945 | 0.945 | 0.514 (0.238, 0.781) | 0.656 (0.108, 4.003) | 1.000 |

`TACSTD2` median difference R−NR = 0.911.
`CLDN4` median difference R−NR = -0.201.

### 3.2 Exploratory two-gene score

Equal-weight z-mean: AUC 0.610
(bootstrap 95% CI 0.333–0.857);
exact label-permutation p=0.546. Not validated.

### 3.3 T-cell proxy and co-expression

The T-cell proxy does not separate response groups: AUC
0.667
(bootstrap 95% CI 0.400–0.886);
exact MW p=0.237;
exact permutation p=0.139.
This is a context check: the cohort is too small to recover even a conventional
immune-infiltration marker.

| Pair | Spearman rho | p |
|---|---:|---:|
| TACSTD2 vs CLDN4 | 0.475 | 0.026 |
| TACSTD2 vs T-cell proxy | -0.010 | 0.966 |
| CLDN4 vs T-cell proxy | -0.164 | 0.465 |

`TACSTD2` and `CLDN4` are moderately correlated, so they are not independent
signals. Neither gene shows an immune-cold correlation here (`rho ≈ 0`).

## 4. Honest reading / 诚实解读

1. **No association.** Neither primary gene, the two-gene score, nor the T-cell
   proxy reaches even an unadjusted 0.05 threshold.
2. **Direction.** `TACSTD2` is higher in responders. That is the opposite of
   “high TROP2 predicts pembrolizumab resistance” and also opposite the
   GSE126044 / GSE207422 direction reported elsewhere in this repository.
   `CLDN4` is essentially random and does not support a CLDN4-high → worse
   response claim in this cohort.
3. **Not a predictive biomarker test.** There is no control arm, so prognostic
   and treatment-interaction effects cannot be separated.
4. **Power.** With 7 vs 15, only large effects are detectable. A null p-value
   is not proof of no effect. The wide AUC intervals (crossing 0.5) are the
   honest result.
5. **Bulk RNA-seq.** Expression mixes tumor, stroma, and purity. No covariates
   are available for adjustment. The T-cell proxy is not a cell count.

## 5. Files / 文件

- `tables/sample_level_expression.csv`
- `tables/association_statistics.csv`
- `tables/gene_correlation.csv`
- `figures/expression_by_response.*`
- `figures/tacstd2_cldn4_scatter.*`
- `figures/tcell_proxy_scatter.*`
- `summary.json`
- `file_manifest.json`
- `analyze.py`, `download.py`

Reproduce:

```bash
python3 results/w200/GSE166449/download.py
python3 results/w200/GSE166449/analyze.py
# or: python3 scripts/w200/GSE166449.py
```
