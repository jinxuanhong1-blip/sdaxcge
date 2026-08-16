# TACSTD2 / CLDN4 vs ICI response in new GEO lung-cancer series (2024–2026)

> Parallel slice. All outputs live under `notes/fable_geo_2026/`,
> `scripts/fable_geo_2026/`, and `results/w200/GEO_2026/`.
> 中文在后 / English first, 中文 second.
> Search and source checks refreshed on **2026-08-16 UTC**.

---

## English

### Objective
Search GEO (2024–2026) for **new human lung-cancer immune-checkpoint-inhibitor
(ICI) series** that are **not** part of the prior slices
(GSE126044, GSE135222, GSE136961, GSE166449, GSE93157, GSE207422, GSE205335),
download the open processed data (< 2 GB each), and test whether tumor
**TACSTD2 (TROP2)** and **CLDN4** expression is associated with ICI response
where per-sample response labels exist.

### Method (reproducible pipeline, `scripts/fable_geo_2026/`)
1. `01_geo_search.py` — NCBI E-utilities `esearch`/`esummary` on `db=gds`:
   lung-cancer terms **AND** ICI terms **AND** `Homo sapiens` **AND**
   `gse[Entry Type]` **AND** publication date 2024/01/01–2026/12/31.
   → **156 hits → 155 candidate GSEs** after removing the 7 excluded accessions.
2. `02_triage.py` — keyword scoring (patient cohort vs cell-line/mouse studies).
3. `03_list_supp.py` — GEO FTP supplementary listing + file sizes (enforce
   "open processed" and "< 2 GB").
4. `04`/`06_*` — series-matrix sample-characteristic scans; **17 of 155** series
   carry a response/outcome field (`results/w200/GEO_2026/label_scan_all.tsv`).
5. `05_download.py` — download qualifying processed matrices.
6. `07_analyze.py` — TACSTD2/CLDN4 vs response (Mann–Whitney U, Cliff's delta,
   Spearman, log-rank).

No accessions were invented — every GSE / GSM / URL is captured from live
E-utilities and GEO FTP responses stored under `results/w200/GEO_2026/`.

### Why these genes need tumor tissue
TACSTD2 (TROP2) and CLDN4 are **epithelial / tumor-cell** genes. They are only
meaningfully measured in tumor tissue (bulk RNA-seq, spatial WTA/CTA, or tumor
scRNA), **not** in blood/PBMC assays. This drove dataset selection.

### Datasets analyzed (tumor tissue **with** ICI response labels)

| GSE | Cohort / regimen | Assay (processed file) | TACSTD2 | CLDN4 | Response label |
|-----|------------------|------------------------|:-------:|:-----:|----------------|
| **GSE261345** | ES-SCLC, chemo-immunotherapy (CANTABRICO) | GeoMx DSP, Cancer Transcriptome Atlas (normalized counts) | ✓ | ✗ (not on CTA panel) | best RECIST + progression/death dates |
| **GSE261348** | ES-SCLC, chemo-immunotherapy (IMfirst) | GeoMx DSP, CTA (normalized counts) | ✓ | ✗ | best RECIST + progression/death dates |
| **GSE233203** | NSCLC (EGFR-TKI resistant), atezolizumab+bevacizumab+chemo (ABCP) | 10x scRNA-seq (7 patients) | ✓ | ✓ | Response / Non-response |

### Results

Effect direction is expressed as **Cliff's delta** (responder vs non-responder;
positive = higher in responders). P-values are two-sided Mann–Whitney U
(group comparison) or Spearman/log-rank (survival). Full table:
`results/w200/GEO_2026/analysis_results.tsv`.

**GSE261345 (ES-SCLC, 26 patients, 121 tumor ROIs)**
- TACSTD2, responder (CR/PR, n=17) vs non-responder (SD/PD, n=9): medians 7.24 vs
  7.36 log2(Q3+1); Cliff's δ = −0.24; **p = 0.33** (ns).
- TACSTD2 vs PFS days: Spearman ρ = 0.15 (p = 0.48); median-split PFS log-rank
  p = 0.59 (ns).

**GSE261348 (ES-SCLC, 32 patients, 175 tumor ROIs)**
- TACSTD2, responder (n=21) vs non-responder (n=9): medians 7.64 vs 7.73;
  Cliff's δ = −0.04; **p = 0.89** (ns).
- TACSTD2 vs PFS days: Spearman ρ = −0.26 (p = 0.15); log-rank p = 0.64 (ns).

**GSE233203 (NSCLC, 3 responders vs 4 non-responders, whole-sample scRNA pseudobulk)**
- TACSTD2: responder median 5.74 vs non-responder 3.70 log2(CPM+1);
  **Cliff's δ = 0.83** (large, higher in responders); p = 0.11 (ns, underpowered).
- CLDN4: responder 5.54 vs non-responder 3.18; Cliff's δ = 0.50; p = 0.40 (ns).
- The markers are strongly concordant across the seven samples (Spearman
  ρ = 0.89, p = 0.0068). This is not an independent outcome test and may
  reflect their shared epithelial program and/or tumor-cell abundance.

### Interpretation
- In the two **ES-SCLC** GeoMx cohorts, tumor **TACSTD2 shows no association**
  with RECIST response or PFS (all p > 0.3, small/near-zero effect sizes).
  CLDN4 is not on the CTA panel, so it could not be tested there.
- In the small **NSCLC** scRNA cohort, both **TACSTD2 and CLDN4 trend higher in
  responders**; the TACSTD2 effect size is large (δ = 0.83) but does **not**
  reach significance with only 3 vs 4 patients.
- Overall: **no statistically significant TACSTD2/CLDN4–response association**
  in the datasets that had usable tumor expression + ICI response labels. The
  NSCLC trend is hypothesis-generating and needs a larger labeled cohort.

### Datasets found but not usable for this specific question (documented, not hidden)
- **GSE309652** (NSCLC, anti-PD-(L)1, clean R/NR, n=72): NanoString **metabolism**
  panel (768 genes) — TACSTD2/CLDN4/EPCAM absent. Downloaded and checked.
- **GSE329813** (NSCLC neoadjuvant chemo-immuno spatial, n=127): expression
  present (TACSTD2 ✓, CLDN4 ✗) but **no per-ROI response/MPR label** in GEO
  metadata (only `batch`, `tissue`). Not used (would require inventing labels).
- **GSE253564** (durvalumab±RT, bulk lung tumor FPKM, full transcriptome, both
  genes present): **no per-sample response label** in GEO (treatment arm only).
- **GSE292421** (pan-cancer incl. NSCLC bulk FPKM, both genes present): only a
  TME immunophenotype (inflamed/excluded/desert), **not response**; and its
  FPKM column IDs (sequencing-chip barcodes) have **no crosswalk** to the GSM
  metadata, so labels cannot be joined reliably.
- Response-labeled **blood/PBMC** series (GSE266219, GSE295969, GSE306542,
  GSE295601, GSE310370) — epithelial genes not meaningfully expressed.
- Cell-line / mechanistic (GSE252437, GSE255144, GSE271377, GSE253718) and the
  pan-cancer scRNA atlas GSE218989 — no usable per-patient lung ICI response.

### Caveats
- Sample sizes are small (especially GSE233203, n=7); p-values are exploratory
  and not multiple-testing corrected.
- DSP values are Q3-normalized counts from mixed "Full ROI" segments;
  patient-level means were used to avoid ROI pseudoreplication.
- GSE233203 values are pseudobulked over **all captured cells**, because GEO
  provides count matrices but no validated cell-type labels. They mix per-cell
  expression with epithelial/tumor-cell abundance and must not be interpreted
  as tumor-cell-intrinsic expression.
- Spearman tests against PFS ignore right-censoring and are descriptive. The
  accompanying log-rank tests use event indicators but rely on an
  information-losing median split. Neither survival analysis is confirmatory.
- GSE261345 biopsy sites include metastatic tissue (e.g. skin) as well as lung;
  all are ES-SCLC patients on chemo-immunotherapy.

### Reproduce
```bash
python3 scripts/fable_geo_2026/01_geo_search.py
python3 scripts/fable_geo_2026/02_triage.py
python3 scripts/fable_geo_2026/03_list_supp.py
python3 scripts/fable_geo_2026/06_label_scan_all.py
python3 scripts/fable_geo_2026/05_download.py
python3 scripts/fable_geo_2026/07_analyze.py
```
Dependencies: `pandas scipy numpy matplotlib openpyxl lifelines`.
Large downloads are reproducible and are not committed (see
`results/w200/GEO_2026/data/.gitignore`).

---

## 中文

> 检索与数据源核查更新于 **2026-08-16 UTC**。

### 目标
在 GEO 数据库中检索 **2024–2026 年新发布的、人类肺癌免疫检查点抑制剂（ICI）相关
系列**，且**不属于**既有切片已覆盖的编号（GSE126044、GSE135222、GSE136961、
GSE166449、GSE93157、GSE207422、GSE205335）；下载开放的已处理数据（每个文件
< 2 GB）；在**存在样本级疗效标签**的数据集中，检验肿瘤 **TACSTD2（TROP2）** 与
**CLDN4** 表达是否与 ICI 疗效相关。

### 方法（可复现流程，见 `scripts/fable_geo_2026/`）
1. `01_geo_search.py` — 通过 NCBI E-utilities（`db=gds`）检索：肺癌关键词 **且**
   ICI 关键词 **且** `Homo sapiens` **且** `gse[Entry Type]` **且** 发表日期
   2024/01/01–2026/12/31。→ **命中 156 条 → 剔除 7 个排除编号后得 155 个候选 GSE**。
2. `02_triage.py` — 关键词打分（区分患者队列与细胞系/小鼠机制研究）。
3. `03_list_supp.py` — 读取 GEO FTP 附件清单与文件大小（确保"开放已处理"且
   "< 2 GB"）。
4. `04`/`06_*` — 扫描 series-matrix 样本特征；155 个中有 **17 个**含疗效/结局字段
   （`results/w200/GEO_2026/label_scan_all.tsv`）。
5. `05_download.py` — 下载符合条件的已处理矩阵。
6. `07_analyze.py` — TACSTD2/CLDN4 与疗效关联分析（Mann–Whitney U、Cliff's
   delta、Spearman、log-rank）。

**未杜撰任何编号**：所有 GSE / GSM / URL 均来自实时 E-utilities 与 GEO FTP 响应，
并保存在 `results/w200/GEO_2026/` 下。

### 为什么这两个基因需要肿瘤组织
TACSTD2（TROP2）与 CLDN4 是**上皮/肿瘤细胞**基因，只有在肿瘤组织
（bulk RNA-seq、空间 WTA/CTA、或肿瘤 scRNA）中才能有意义地测得，**不能**用
血液/PBMC 数据评估。这决定了数据集的筛选。

### 纳入分析的数据集（肿瘤组织 **且** 具备 ICI 疗效标签）

| GSE | 队列/方案 | 平台（已处理文件） | TACSTD2 | CLDN4 | 疗效标签 |
|-----|-----------|--------------------|:-------:|:-----:|----------|
| **GSE261345** | 广泛期小细胞肺癌，化免联合（CANTABRICO） | GeoMx DSP，肿瘤转录组图谱 CTA（归一化计数） | ✓ | ✗（CTA 面板不含） | 最佳 RECIST + 进展/死亡日期 |
| **GSE261348** | 广泛期小细胞肺癌，化免联合（IMfirst） | GeoMx DSP，CTA（归一化计数） | ✓ | ✗ | 最佳 RECIST + 进展/死亡日期 |
| **GSE233203** | 非小细胞肺癌（EGFR-TKI 耐药后），阿替利珠+贝伐+化疗（ABCP） | 10x 单细胞 RNA-seq（7 例） | ✓ | ✓ | 缓解/未缓解 |

### 结果

效应方向以 **Cliff's delta**（缓解 vs 未缓解，正值＝缓解者更高）表示；P 值为
双侧 Mann–Whitney U（组间）或 Spearman/log-rank（生存）。完整结果见
`results/w200/GEO_2026/analysis_results.tsv`。

**GSE261345（广泛期 SCLC，26 例，121 个肿瘤 ROI）**
- TACSTD2：缓解者（CR/PR，n=17）vs 未缓解者（SD/PD，n=9）中位数 7.24 vs 7.36；
  Cliff's δ = −0.24；**p = 0.33**（不显著）。
- TACSTD2 与 PFS 天数：Spearman ρ = 0.15（p = 0.48）；中位切分 log-rank
  p = 0.59（不显著）。

**GSE261348（广泛期 SCLC，32 例，175 个肿瘤 ROI）**
- TACSTD2：缓解者（n=21）vs 未缓解者（n=9）中位数 7.64 vs 7.73；
  Cliff's δ = −0.04；**p = 0.89**（不显著）。
- TACSTD2 与 PFS 天数：Spearman ρ = −0.26（p = 0.15）；log-rank p = 0.64（不显著）。

**GSE233203（NSCLC，缓解 3 例 vs 未缓解 4 例，全样本 scRNA 拟 bulk）**
- TACSTD2：缓解者中位 5.74 vs 未缓解 3.70 log2(CPM+1)；**Cliff's δ = 0.83**
  （效应大，缓解者更高）；p = 0.11（不显著，样本量不足）。
- CLDN4：缓解者 5.54 vs 未缓解 3.18；Cliff's δ = 0.50；p = 0.40（不显著）。
- 两个标志物在 7 个样本间高度一致（Spearman ρ = 0.89，p = 0.0068）。这并非独立的
  疗效检验，可能反映二者共享的上皮程序和/或肿瘤细胞丰度。

### 解读
- 在两个 **广泛期 SCLC** GeoMx 队列中，肿瘤 **TACSTD2 与 RECIST 疗效或 PFS 均无
  关联**（p 全部 > 0.3，效应量极小）。CTA 面板不含 CLDN4，无法在此检验。
- 在小样本 **NSCLC** scRNA 队列中，**TACSTD2 与 CLDN4 在缓解者中均有升高趋势**；
  TACSTD2 效应量大（δ = 0.83），但在仅 3 vs 4 例下**未达统计显著**。
- 总体：在具备可用肿瘤表达 + ICI 疗效标签的数据集中，**未发现 TACSTD2/CLDN4 与
  疗效的统计学显著关联**。NSCLC 中的趋势属于产生假设的线索，需更大的带标签队列
  验证。

### 检索到但不适用于本问题的数据集（如实记录，非隐瞒）
- **GSE309652**（NSCLC，抗 PD-(L)1，清晰 R/NR，n=72）：NanoString **代谢**面板
  （768 基因），不含 TACSTD2/CLDN4/EPCAM。已下载并核实。
- **GSE329813**（NSCLC 新辅助化免空间转录组，n=127）：表达存在（TACSTD2 ✓、
  CLDN4 ✗），但 GEO 元数据**无每个 ROI 的疗效/MPR 标签**（仅 `batch`、`tissue`）。
  未使用（否则需杜撰标签）。
- **GSE253564**（durvalumab±放疗，肺肿瘤 bulk FPKM，全转录组，两基因均在）：
  GEO 中**无样本级疗效标签**（仅治疗分组 Arm1/Arm2）。
- **GSE292421**（泛癌含 NSCLC，bulk FPKM，两基因均在）：仅有 TME 免疫表型
  （inflamed/excluded/desert），**并非疗效**；且其 FPKM 列名（测序芯片条码）与
  GSM 元数据**无对应关系**，无法可靠连接标签。
- 具疗效标签的**血液/PBMC**系列（GSE266219、GSE295969、GSE306542、GSE295601、
  GSE310370）——上皮基因几乎不表达。
- 细胞系/机制研究（GSE252437、GSE255144、GSE271377、GSE253718）与泛癌 scRNA
  图谱 GSE218989——无可用的患者级肺癌 ICI 疗效。

### 注意事项
- 样本量偏小（尤其 GSE233203，n=7）；P 值为探索性，未做多重检验校正。
- DSP 数值为混合 "Full ROI" 片段的 Q3 归一化计数；采用患者级均值以避免 ROI 伪重复。
- GSE233203 仅提供计数矩阵、无经验证的细胞类型标签，因此数值是在**全部捕获细胞**
  上拟 bulk；它混合了单细胞表达与上皮/肿瘤细胞丰度，不能解读为肿瘤细胞内在表达。
- TACSTD2 与 PFS 的 Spearman 分析忽略右删失，仅为描述性；配套 log-rank 使用事件
  指标但依赖损失信息的中位数切分。两种生存分析均非确证性。
- GSE261345 活检部位含转移灶（如皮肤）与肺；均为接受化免联合的广泛期 SCLC 患者。

### 复现
```bash
python3 scripts/fable_geo_2026/01_geo_search.py
python3 scripts/fable_geo_2026/02_triage.py
python3 scripts/fable_geo_2026/03_list_supp.py
python3 scripts/fable_geo_2026/06_label_scan_all.py
python3 scripts/fable_geo_2026/05_download.py
python3 scripts/fable_geo_2026/07_analyze.py
```
依赖：`pandas scipy numpy matplotlib openpyxl lifelines`。大文件可复现且不纳入
版本库（见 `results/w200/GEO_2026/data/.gitignore`）。
