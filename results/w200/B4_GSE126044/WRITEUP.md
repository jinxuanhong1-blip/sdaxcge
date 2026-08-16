# B4 exact recompute: GSE126044 CLDN4 and TJ score, R vs NR

# B4 精确重算：GSE126044 CLDN4 与紧密连接（TJ）评分，应答者 vs 非应答者

**Verdict / 结论:** the claimed **“NR higher TJ, p = 0.019” is not recovered.**
Direction is consistently **higher in non-responders**. The primary two-sided
Mann–Whitney p-values are **0.115 (CLDN4)** and **0.115 (core TJ mean-z)**,
not 0.019. The closest a-priori TJ module in this repo is `CUSTOM_TJ_CORE`
at **p = 0.069**. No pre-specified definition reaches p = 0.019.

**结论：** 声称的 **“NR 更高 TJ，p = 0.019”无法复现。** 方向一致为
**非应答者更高**。主要双侧 Mann–Whitney p 值为 **CLDN4 = 0.115**、
**核心 TJ mean-z = 0.115**，不是 0.019。本仓库先验 TJ 模块中最接近的是
`CUSTOM_TJ_CORE`（**p = 0.069**）。没有任何预先指定的定义达到 p = 0.019。

---

## 1. What was claimed / 原声称

Batch B4 / w200 cited GSE126044 as **NR higher TJ, p = 0.019**.
This folder recomputes that contrast from the raw count matrix, honestly,
for **CLDN4 alone** and for **tight-junction signature scores**.

批次 B4 / w200 引用 GSE126044 为 **NR 更高 TJ，p = 0.019**。
本目录从原始计数矩阵诚实重算 **CLDN4 单基因** 与 **紧密连接签名评分**。

## 2. Data / 数据

| Item | Value |
|------|--------|
| GEO | [GSE126044](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126044) (Cho et al., NSCLC, pre-treatment anti-PD-1) |
| Matrix | `GSE126044_counts.txt.gz` — integer gene × sample counts, HGNC symbols (18,747 genes × 16 samples) |
| Labels | GEO `patient_response`: 5 responder / 11 non-responder |
| Sample type | 11 fresh, 5 FFPE. **All 5 FFPE are non-responders** (response fully confounded with FFPE) |
| Source in this repo | same count file as `fable_geo_2019_2021` / `fable_ici_bulk` |

Per-sample values: `tables/per_sample_scores.csv`.
Response × sample-type: `tables/response_by_sample_type.csv`.

## 3. Methods (pre-specified) / 方法（预先指定）

- **Normalisation.** `log2 CPM = log2(1e6 × count / library_size + 1)`.
  Same as the existing GSE126044 pipelines in this repo, so CLDN4 numbers
  are comparable.
- **Primary test.** Two-sided Mann–Whitney U (Wilcoxon rank-sum),
  5 responders vs 11 non-responders. Welch t-test is secondary.
- **Effect size.** AUC = P(R > NR) from the U statistic; rank-biserial
  = 2·AUC − 1; Cohen’s d (R − NR).
- **TJ score (primary).** Unweighted **mean of per-gene z-scores**
  (z across the 16 samples on log2 CPM) over a curated structural TJ set
  (claudins, OCLN/MARVELD2/3, TJP1–3, JAMs, cingulin, apical polarity
  scaffolds; 40 genes present with non-zero variance). Gene list:
  `data/geneset_core_tj.txt`.
- **A-priori sensitivity sets (not p-hunted).**
  - core TJ without CLDN4
  - core TJ as mean log2 CPM (no z-score)
  - `CUSTOM_TJ_CORE` (26 genes) and `CUSTOM_CLAUDIN_PAR_FOCAL` (5 genes)
    from this repo’s `hunt-tj-gsea` module
  - KEGG Tight Junction `hsa04530` (current REST dump, 163 genes used)
- **Fresh-only sensitivity.** Drop the 5 FFPE samples (all NR) and
  retest CLDN4 / core TJ / CUSTOM_TJ_CORE on 5 R vs 6 fresh NR.
- Script: `scripts/w200/B4_GSE126044_recompute.py`.

Nothing was selected because it was closer to 0.019.

未因更接近 0.019 而挑选基因集、检验或子集。

## 4. Results / 结果

### 4.1 All samples (primary; n = 16 = 5 R + 11 NR)

| Feature | n genes | Median R | Median NR | Direction | AUC (R>NR) | MW p | Welch p |
|---------|---------|----------|-----------|-----------|------------|------|---------|
| **CLDN4 log2 CPM** | 1 | 2.54 | 3.77 | NR higher | 0.236 | **0.1149** | 0.166 |
| TACSTD2 log2 CPM (reference) | 1 | 6.17 | 6.29 | NR higher | 0.364 | 0.4409 | 0.322 |
| **core TJ (mean-z)** | 40 | −0.076 | 0.157 | NR higher | 0.236 | **0.1149** | 0.152 |
| core TJ, no CLDN4 | 39 | −0.068 | 0.150 | NR higher | 0.255 | 0.1451 | 0.153 |
| core TJ, mean expr | 40 | 3.36 | 3.54 | NR higher | 0.291 | 0.2212 | 0.212 |
| CUSTOM_TJ_CORE (mean-z) | 24 | −0.122 | 0.170 | NR higher | 0.200 | 0.0687 | 0.164 |
| CUSTOM_CLAUDIN_PAR_FOCAL | 5 | −0.018 | 0.188 | NR higher | 0.327 | 0.3196 | 0.203 |
| KEGG hsa04530 (mean-z) | 163 | −0.001 | 0.092 | NR higher | 0.364 | 0.4409 | 0.349 |
| KEGG hsa04530, mean expr | 163 | 4.81 | 4.86 | NR higher | 0.400 | 0.5833 | 0.358 |

CLDN4 medians and p = 0.1149 **exactly match** two independent prior
analyses in this repo (`fable_ici_bulk`, `fable_geo_2019_2021`).
That number is real. It is not 0.019.

CLDN4 中位数与 p = 0.1149 **精确吻合**本仓库中两项独立分析
（`fable_ici_bulk`、`fable_geo_2019_2021`）。该数字是真实的，但不是 0.019。

### 4.2 Fresh biopsies only (sensitivity; 5 R vs 6 NR)

All five FFPE libraries are NR, so sample type and response are
linearly dependent. Restricting to fresh tissue does **not** create
a significant result:

| Feature | MW p (fresh only) |
|---------|-------------------|
| CLDN4 | 0.329 |
| core TJ mean-z | 0.178 |
| CUSTOM_TJ_CORE | 0.126 |

### 4.3 Could a one-sided test rescue p = 0.019? / 单侧检验能否得到 0.019？

No. A one-sided Mann–Whitney in the claimed direction (NR higher)
is approximately half the two-sided p:

- CLDN4 / core TJ: ~0.057
- CUSTOM_TJ_CORE: ~0.034

Still not 0.019.

不能。按“NR 更高”方向的单侧 Mann–Whitney 约为双侧 p 的一半
（CLDN4/核心 TJ ≈ 0.057；CUSTOM_TJ_CORE ≈ 0.034），仍不是 0.019。

## 5. Honest reading / 诚实解读

1. **Direction.** Every pre-specified CLDN4 / TJ contrast is higher in
   non-responders. That part of the B4 claim is consistent with the data.
2. **Significance.** With n = 16 (5 vs 11), none of the pre-specified
   tests is significant at α = 0.05, and **none is p = 0.019**.
   Calling this “NR higher TJ p = 0.019” overstates the evidence.
3. **CLDN4 vs TJ score.** The primary 40-gene TJ mean-z has the **same**
   Mann–Whitney U (13) and p (0.1149) as CLDN4 alone. The signature does
   not add a stronger response association than the single gene.
4. **Confounding.** FFPE vs fresh is completely aligned with NR
   (5/5 FFPE = NR). Any analysis that ignores sample type can mix a
   technical effect into the response contrast.
5. **Where 0.019 might have come from.** This recompute does not recover
   that exact p under any a-priori TJ definition, scoring method
   (mean-z vs mean expression), or the fresh-only subset. Possible
   sources of a stray 0.019 (ssGSEA/GSVA, a different gene set, a
   different test, a different sample filter, or a different dataset
   mis-labelled as GSE126044) were **not** reverse-engineered.
   Hunting until a gene set hits 0.019 would not be honest.

1. **方向。** 所有预先指定的 CLDN4 / TJ 对比均为非应答者更高。B4 声称的方向与数据一致。
2. **显著性。** n = 16（5 vs 11）下，预先指定检验均未达 α = 0.05，**也没有 p = 0.019**。
   写成“NR 更高 TJ p = 0.019”夸大了证据。
3. **CLDN4 与 TJ 评分。** 40 基因 TJ mean-z 的 U（13）和 p（0.1149）与 CLDN4 单基因相同。
   签名并不比单基因更强。
4. **混杂。** FFPE 与应答完全纠缠（5/5 FFPE = NR）。忽略样本类型会把技术效应混进应答对比。
5. **0.019 可能从何而来。** 在先验 TJ 定义、评分方法和 fresh-only 子集下均未得到该 p 值。
   未反向工程 ssGSEA/GSVA、其他基因集或其他过滤。为凑 0.019 去搜基因集不诚实。

## 6. Files / 文件

```
results/w200/B4_GSE126044/
  WRITEUP.md
  summary.json
  data/GSE126044_counts.txt.gz
  data/GSE126044_clinical.csv
  data/kegg_hsa04530.txt
  data/geneset_*.txt
  tables/R_vs_NR_stats.csv
  tables/per_sample_scores.csv
  tables/geneset_membership.csv
  tables/response_by_sample_type.csv
  figures/CLDN4_R_vs_NR.png
  figures/coreTJ_R_vs_NR.png
  figures/customTJ_R_vs_NR.png
  figures/keggTJ_R_vs_NR.png
scripts/w200/B4_GSE126044_recompute.py
```
