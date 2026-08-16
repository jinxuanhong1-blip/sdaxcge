# TACSTD2 / CLDN4 vs ICI outcomes in GEO 2019–2021 human lung series
# GEO 2019–2021 人类肺癌 ICI 系列中 TACSTD2 / CLDN4 与免疫治疗结局的关系

> Parallel slice. All outputs live only under `notes/fable_geo_2019_2021/`,
> `scripts/fable_geo_2019_2021/`, `results/fable_geo_2019_2021/`.
> All GEO accessions are real NCBI identifiers; none were invented.
> 本平行切片的全部产物仅位于上述三个目录；所有 GEO 编号均为 NCBI 真实编号，无任何虚构 ID。

---

## English

### 1. Objective
Exhaustively identify human **lung-cancer immune-checkpoint-inhibitor (ICI)**
GEO series published **2019–2021**, verify them, download the open processed
data (< 2 GB each), and test whether tumor expression of **TACSTD2** and
**CLDN4** is associated with ICI outcomes.

### 2. Search & verification
An NCBI E-utilities search of the `gds` database (GSE, *Homo sapiens*, PDAT
2019–2021, lung × ICI terms) returned **90 real series**. Rule-based triage
plus manual review of metadata (`results/tables/triage_all_candidates.csv`)
yielded **4 analyzable** cohorts with processed expression matrices, per-sample
ICI outcomes, and measurable target genes, plus **1** genuine NSCLC anti-PD-1
cohort (GSE136961) that had to be excluded because its Oncomine 395-gene immune
panel does not measure TACSTD2/CLDN4. The remaining 85 were single-cell,
methylation, cell-line/in-vitro, non-lung, or non-ICI studies. See
`verification_log.md` for the full audit trail.

### 3. Datasets analyzed
| GSE | Cohort | n | Outcome model |
|---|---|---|---|
| GSE126044 | NSCLC tumor biopsy, anti-PD-1 | 16 (5 R / 11 NR) | responder vs non-responder |
| GSE135222 | NSCLC tumor, anti-PD-1/PD-L1 | 27 | progression-free survival (PFS) |
| GSE182328 | advanced/limited lung tumor, ICI-treated | 44 | Akkermansia group (prognostic **surrogate**) |
| GSE111414 | PBMC CD8+ T cells, nivolumab | 20 | responder vs non-responder (epithelial-null control) |

Normalisation: counts → log2 CPM; TPM → log2(TPM+1). Tests: Mann–Whitney U with
AUC for binary outcomes; Cox regression (per log2 unit) + Kaplan–Meier
median-split log-rank for PFS.

### 4. Results
Full table: `results/tables/combined_TACSTD2_CLDN4_vs_ICI_outcomes.csv`.
Figures: `results/figures/`.

- **GSE126044 (response):** both genes trend **higher in non-responders**
  (TACSTD2 AUC 0.36, p = 0.44; CLDN4 AUC 0.24, p = 0.11). Neither significant.
- **GSE135222 (PFS):** no association (TACSTD2 Cox HR 1.04, p = 0.78, log-rank
  p = 0.17; CLDN4 HR 1.05, p = 0.61, log-rank p = 0.91).
- **GSE182328 (Akkermansia surrogate):** CLDN4 higher in the Akkermansia-negative
  (worse-prognosis) group (AUC 0.32, nominal **p = 0.045**, uncorrected);
  TACSTD2 not significant (p = 0.43).
- **GSE111414 (control):** TACSTD2/CLDN4 ≈ 0 in CD8+ T cells, as expected;
  not interpretable.

### 5. Interpretation
Across the verified 2019–2021 GEO human lung ICI cohorts, **neither TACSTD2 nor
CLDN4 shows a statistically robust association with ICI response or PFS**. There
is a **consistent, non-significant trend** that higher TACSTD2/CLDN4 marks
non-responder / poorer-prognosis tumors (both genes higher in GSE126044
non-responders; CLDN4 higher in GSE182328 Akkermansia-negative tumors, nominal
p = 0.045). This is **hypothesis-generating only**.

### 6. Limitations
Small n (16–44); underpowered; no multiple-testing correction; GSE182328 uses a
microbiome surrogate rather than a deposited clinical endpoint; results are
exploratory and not a substitute for a pooled, harmonised meta-analysis.

### 7. Reproducibility
Run in order: `01_search_geo.py` → `02_fetch_metadata.py` → `03_download.py`
→ `04_build_clinical.py` → `05_analyze.py` → `06_triage.py` (Python 3,
pandas/numpy/scipy/lifelines/matplotlib).

### 8. Leftover / no-skip pass (do not skip for size or tissue)
A second pass revisited all **53 leftover** lung+ICI series (blood, large scRNA,
RAW.tar-only). Full audit: `results/noskip/GEO_2019_2021/` (catalog, gene scan,
computes, bilingual writeup).

The only leftover series with a processed target-gene matrix **and** a deposited
patient ICI-related survival endpoint is **GSE190266** (n=70, CLDN4 present,
TACSTD2 absent): 6-month PFS Cox HR = 0.917, p = 0.1246; log-rank p = 0.2105.
**Not significant.** Blood leftovers GSE152590 / GSE141479 were downloaded;
target gene symbols are absent. GSE176021 has response labels but no processed
expression matrix — not computed. **No statistics were fabricated.**

See `results/noskip/GEO_2019_2021/WRITEUP.md`.

---

## 中文

### 1. 目标
穷尽检索 **2019–2021 年**发表、涉及**人类肺癌免疫检查点抑制剂（ICI）**的 GEO
系列，进行核验，下载开放的已处理数据（每个 < 2 GB），并检验肿瘤中 **TACSTD2** 与
**CLDN4** 的表达是否与 ICI 治疗结局相关。

### 2. 检索与核验
使用 NCBI E-utilities 对 `gds` 数据库检索（GSE、人类、发表日期 2019–2021、肺癌 ×
ICI 关键词），共返回 **90 个真实系列**。基于规则的分诊加人工元数据审阅
（`results/tables/triage_all_candidates.csv`）筛出 **4 个可分析**队列（具备已处理
表达矩阵、样本级 ICI 结局、且目标基因可测量），另有 **1 个**真实的 NSCLC 抗 PD-1
队列（GSE136961）因其 Oncomine 395 基因免疫面板不测 TACSTD2/CLDN4 而被排除。其余
85 个为单细胞、甲基化、细胞系/体外、非肺癌或非 ICI 研究。完整审计见
`verification_log.md`。

### 3. 纳入分析的数据集
| GSE | 队列 | n | 结局模型 |
|---|---|---|---|
| GSE126044 | NSCLC 肿瘤活检，抗 PD-1 | 16（5 应答 / 11 非应答） | 应答 vs 非应答 |
| GSE135222 | NSCLC 肿瘤，抗 PD-1/PD-L1 | 27 | 无进展生存（PFS） |
| GSE182328 | 晚期/局限期肺肿瘤，ICI 治疗 | 44 | Akkermansia 分组（预后**替代指标**） |
| GSE111414 | 外周血 CD8+ T 细胞，纳武利尤单抗 | 20 | 应答 vs 非应答（上皮基因阴性对照） |

标准化：计数 → log2 CPM；TPM → log2(TPM+1)。检验：二分类用 Mann–Whitney U 及 AUC；
PFS 用 Cox 回归（每 log2 单位）与 Kaplan–Meier 中位数分组 log-rank。

### 4. 结果
完整表：`results/tables/combined_TACSTD2_CLDN4_vs_ICI_outcomes.csv`；图：`results/figures/`。

- **GSE126044（应答）：** 两个基因均**在非应答者中偏高**（TACSTD2 AUC 0.36，p =
  0.44；CLDN4 AUC 0.24，p = 0.11），均不显著。
- **GSE135222（PFS）：** 无相关（TACSTD2 HR 1.04，p = 0.78，log-rank 0.17；CLDN4
  HR 1.05，p = 0.61，log-rank 0.91）。
- **GSE182328（Akkermansia 替代）：** CLDN4 在 Akkermansia 阴性（预后较差）组更高
  （AUC 0.32，名义 **p = 0.045**，未校正）；TACSTD2 不显著（p = 0.43）。
- **GSE111414（对照）：** CD8+ T 细胞中 TACSTD2/CLDN4 ≈ 0，符合预期，不可解释。

### 5. 结论解读
在已核验的 2019–2021 GEO 人类肺癌 ICI 队列中，**TACSTD2 与 CLDN4 均未显示与 ICI
应答或 PFS 的稳健统计学关联**。存在一个**一致但不显著的趋势**：TACSTD2/CLDN4 越高，
越倾向于非应答/预后较差的肿瘤（GSE126044 两基因在非应答者更高；GSE182328 中 CLDN4
在 Akkermansia 阴性肿瘤更高，名义 p = 0.045）。该结论**仅具假设生成价值**。

### 6. 局限
样本量小（16–44），检验效能不足；未做多重比较校正；GSE182328 使用微生物组替代指标而
非 GEO 中登记的临床终点；结果为探索性，不能替代经过统一整合的合并 meta 分析。

### 7. 可复现
按序运行：`01_search_geo.py` → `02_fetch_metadata.py` → `03_download.py` →
`04_build_clinical.py` → `05_analyze.py` → `06_triage.py`（Python 3，依赖
pandas/numpy/scipy/lifelines/matplotlib）。

### 8. 遗留 / 不跳过体积与组织类型
第二轮复查全部 **53 个遗留** 肺癌+ICI 系列（血液、大体积 scRNA、仅有 RAW.tar）。
完整审计见 `results/noskip/GEO_2019_2021/`。

遗留中唯一同时具备已处理目标基因矩阵与患者 ICI 相关生存终点的是 **GSE190266**
（n=70，有 CLDN4、无 TACSTD2）：6 个月 PFS Cox HR=0.917，p=0.1246；log-rank
p=0.2105。**不显著。** 血液系列 GSE152590 / GSE141479 已下载，目标基因符号不在
沉积矩阵中。GSE176021 有应答标签但无已处理表达矩阵，故未计算。**未编造任何统计量。**
