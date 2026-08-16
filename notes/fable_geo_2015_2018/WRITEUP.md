# GEO 2015–2018 Human Lung ICI / PD-1 / PD-L1 / CTLA-4 Series — Exhaustive Mining & TACSTD2/CLDN4 Analysis

> Parallel slice. All outputs live under `notes/fable_geo_2015_2018/`,
> `scripts/fable_geo_2015_2018/`, and `results/fable_geo_2015_2018/`.
> **No invented IDs**: every accession below was returned by NCBI E-utilities and
> re-verified against its authoritative GEO SOFT record.

Bilingual write-up: **English first, 中文在后。**

---

## 1. Objective (EN)

Exhaustively enumerate **GEO Series (GSE) released 2015–2018** that concern
**human lung cancer** and **immune-checkpoint biology / therapy**
(ICI, PD-1/PDCD1, PD-L1/CD274, CTLA-4). For every accession: verify it, download
all **open processed** files **< 2 GB**, and — **where treatment-response labels
exist** — test whether **TACSTD2 (Trop-2)** and **CLDN4 (Claudin-4)** expression
associates with response.

## 2. Pipeline & reproducibility (EN)

Five ordered scripts (pure Python + NCBI public endpoints, no API key):

| step | script | what it does | key output |
|---|---|---|---|
| 1 | `01_search.py` | Two complementary GEO `gds` queries (lung × ICI × human × GSE × 2015–2018 PDAT), unioned | `search_candidates.tsv` (47 GSE) |
| 2 | `02_verify.py` | Fetch SOFT-brief + FTP `suppl/`+`matrix/` listings for **every** accession | `verified_series.tsv`, `suppl_files.tsv` |
| 3 | `03_download.py` | Download all open processed files **< 2 GB** for in-scope series | `download_manifest.tsv`, `data/` |
| 4 | `04_build_clinical.py` | Parse per-sample characteristics from series matrices; flag response/survival labels | `clinical/*.tsv`, `clinical_label_index.tsv` |
| 5 | `05_target_gene_analysis.py` | Gene-availability audit + TACSTD2/CLDN4 vs response/survival | `analysis/*` |
| 6 | `06_curate.py` | Join everything into one master table | `master_summary.tsv` |

Run in order from the repo root; step 3 re-creates the (git-ignored) `data/`.

## 3. Search & verification results (EN)

- **47 real GSE** matched and were verified. Organism: **30 human, 16 mouse,
  1 human+mouse**. Public-date years: 3 (2015), 12 (2016), 11 (2017), 21 (2018) —
  all inside the 2015–2018 window.
- **23 series are in scope** (`is_human` AND `is_lung` AND `is_ici`); the other 24
  are retained in `verified_series.tsv` with the reason they fall out (mouse model,
  or human-but-non-lung such as melanoma/RCC/AML, or a false "PD" hit — e.g.
  **GSE87879** matched only via *"PD 0332991"* = palbociclib, a CDK4/6 inhibitor,
  and is correctly flagged `is_ici = False`).
- **Downloads:** **79 files, ~3.18 GB, 0 failures.** Exactly **one** file exceeded
  the limit and was skipped: `GSE72094_RAW.tar` (2.1 GB) — its expression is still
  obtained from the series matrix.

The 23 in-scope accessions:

```
GSE72094 GSE81089 GSE84789 GSE84797 GSE90728 GSE90729 GSE91061 GSE93157
GSE99254 GSE99531 GSE100860 GSE101929 GSE102286 GSE106420 GSE108819 GSE109010
GSE109020 GSE110390 GSE111360 GSE113972 GSE115305 GSE121682 GSE124199
```

Full per-series evidence (dates, PubMed, type, sample count, labels, gene
measurability, download size) is in
[`results/fable_geo_2015_2018/master_summary.tsv`](../../results/fable_geo_2015_2018/master_summary.tsv).

## 4. The central finding: labels vs. measured genes (EN)

The task's conditional — *"TACSTD2/CLDN4 vs response **if labels exist**"* — exposes
a real structural gap in the 2015–2018 human-lung ICI corpus:

| Category | Series | Response/outcome labels? | TACSTD2 / CLDN4 measured? |
|---|---|---|---|
| **Lung + ICI treatment-response** | **GSE93157** (nivolumab/pembrolizumab; nCounter 775-gene panel) | **Yes** — `best.resp` (RECIST CR/PR/SD/PD), binary `response`, `pfs`, `drug` | **No** — immune panel; only CD274 present |
| | **GSE110390** (durvalumab; 21-gene IFN-γ panel) | No per-sample labels in GEO | **No** — only CD274 |
| **Lung + genome-wide, but survival (not ICI)** | **GSE72094** (LUAD array, n=442) | Overall survival (`vital_status`, `survival_time_in_days`) | **Yes** |
| | **GSE81089** (NSCLC RNA-seq, n=199) | Overall survival (`dead`, dates) | **Yes** |
| **Genome-wide + ICI response, but melanoma** | **GSE91061** (nivolumab, Riaz 2017) | **Yes** — RECIST `response`, pre/on-treatment | **Yes** |

**In short:** within 2015–2018, the human **lung** series that carry genuine ICI
**response** labels use targeted immune panels that **do not measure TACSTD2 or
CLDN4**, while the lung series that **do** measure those epithelial genes
genome-wide are surgical cohorts with **overall-survival** (not ICI-response)
endpoints. We therefore report exactly what the data support, and add a
matched-design **melanoma ICI benchmark** (GSE91061) so the intended
response analysis is still demonstrated end-to-end — clearly labelled as
non-lung. Gene measurability per processed file:
[`analysis/gene_availability_audit.tsv`](../../results/fable_geo_2015_2018/analysis/gene_availability_audit.tsv).

## 5. Quantitative results (EN)

### 5a. GSE91061 — melanoma, nivolumab, RECIST response (ICI benchmark)
Pre-treatment tumors, log2(FPKM+1), Mann-Whitney U (PR/CR responders vs PD
non-responders):

| gene | n(PR/CR) | n(PD) | median resp | median non-resp | U | p |
|---|---|---|---|---|---|---|
| TACSTD2 | 10 | 23 | 0.425 | 0.716 | 97 | **0.49** |
| CLDN4 | 10 | 23 | 0.246 | 0.155 | 133 | **0.49** |
| CD274 (PD-L1) | 10 | 23 | 1.832 | 1.578 | 129 | 0.60 |

No significant TACSTD2/CLDN4 difference by response in this cohort.
Plot: `analysis/GSE91061_response_boxplots.png`.

### 5b. GSE72094 — LUAD (n=398 analysable), overall survival
Median-split log-rank; Spearman vs CD274:

| gene | log-rank p (OS) | ρ vs CD274 | p(ρ) |
|---|---|---|---|
| TACSTD2 | 0.41 | −0.016 | 0.75 |
| CLDN4 | 0.58 | −0.014 | 0.78 |

### 5c. GSE81089 — NSCLC (n=196 analysable), overall survival
| gene | log-rank p (OS) | ρ vs CD274 | p(ρ) |
|---|---|---|---|
| TACSTD2 | 0.33 | −0.005 | 0.95 |
| CLDN4 | 0.73 | −0.112 | 0.12 |

Neither TACSTD2 nor CLDN4 shows a significant overall-survival split, and both are
essentially uncorrelated with PD-L1/CD274 in these lung cohorts. Plots:
`analysis/GSE72094_survival_km.png`, `analysis/GSE81089_survival_km.png`.
Full numbers: `analysis/analysis_results.json`.

## 6. Interpretation & limitations (EN)

- **Honest negatives.** With no in-window human-lung dataset providing both ICI
  response labels *and* TACSTD2/CLDN4 measurements, no lung-specific
  TACSTD2/CLDN4-vs-ICI-response claim can be made from 2015–2018 GEO. The
  melanoma benchmark and the two lung survival cohorts are all null for these
  genes; this is a reportable result, not a pipeline failure.
- **Classification caveats.** Text-based scope flags are inclusive: **GSE91061** is
  melanoma (its abstract merely mentions NSCLC) and **GSE93157** is a mixed
  lung/HNSCC/melanoma cohort (35/65 samples are lung). Both are treated with the
  correct tissue label in the analysis. **GSE101929/GSE102286** (African- vs
  European-American NSCLC) and **GSE109010/GSE109020** (SWI/SNF) match on PD-L1 /
  "immunotherapy" mentions but are not ICI-treatment studies.
- **Scope of "processed".** Downloads target open processed supplementary +
  series-matrix files < 2 GB; raw SRA reads were not fetched. Sorted-cell series
  (e.g. GSE99531, GSE99254 T cells) technically contain TACSTD2/CLDN4 rows but at
  ~0 (epithelial genes in lymphocytes) and lack per-sample response labels.

## 7. Output inventory (EN)

- `results/fable_geo_2015_2018/master_summary.tsv` — 47 verified series, joined evidence.
- `.../verified_series.tsv|json`, `.../search_candidates.tsv`, `.../suppl_files.tsv` — search + verification.
- `.../download_manifest.tsv|json` — every file, size, status, URL (0 failures, 1 skipped ≥2 GB).
- `.../clinical/*.tsv`, `.../clinical_label_index.tsv` — per-sample metadata + label flags.
- `.../analysis/analysis_results.json`, `.../analysis/gene_availability_audit.tsv`, 3 PNGs.
- `.../data/` — the downloaded files (git-ignored; reproducible via `03_download.py`).

---

# 中文版本

## 1. 目标
穷尽式检索 **2015–2018 年发布的 GEO 系列 (GSE)**，聚焦 **人类肺癌** 与
**免疫检查点 (ICI，PD-1/PDCD1、PD-L1/CD274、CTLA-4)** 相关研究。对每个编号：
逐一**核实**，下载全部 **公开的已处理** 且 **小于 2 GB** 的文件；并在**存在治疗
应答标签时**，检验 **TACSTD2 (Trop-2)** 与 **CLDN4 (Claudin-4)** 的表达是否与应答相关。
**绝不编造编号**——所有编号均来自 NCBI E-utilities，并经 GEO SOFT 原始记录复核。

## 2. 流程与可复现性
六个有序脚本（纯 Python + NCBI 公共接口，无需 API key），见上表；从仓库根目录
依次运行即可，第 3 步会重新生成被 git 忽略的 `data/`。

## 3. 检索与核实结果
- 共匹配并核实 **47 个真实 GSE**：**人类 30、小鼠 16、人鼠混合 1**；公开年份
  分布为 2015 年 3 个、2016 年 12 个、2017 年 11 个、2018 年 21 个，全部落在窗口内。
- **23 个属于范围内**（人类 且 肺 且 ICI）；其余 24 个仍保留在 `verified_series.tsv`
  中并注明剔除原因（小鼠模型、或人类但非肺癌如黑色素瘤/肾癌/白血病、或 "PD" 误命中
  ——例如 **GSE87879** 仅因 *"PD 0332991"*（palbociclib，CDK4/6 抑制剂）被命中，已
  正确标记 `is_ici = False`）。
- **下载：79 个文件，约 3.18 GB，0 失败。** 仅 **1 个** 超过限制被跳过：
  `GSE72094_RAW.tar`（2.1 GB）——其表达数据仍可从 series matrix 获得。

## 4. 核心发现：有标签的没测基因，测了基因的没标签
任务的条件句“**若存在标签**才做 TACSTD2/CLDN4 vs 应答”恰好暴露了 2015–2018
人肺 ICI 语料的结构性缺口（详见第 4 节英文表）：

- **有真正 ICI 应答标签的肺癌系列**（GSE93157 nivolumab/pembrolizumab，nCounter
  775 基因 panel；GSE110390 durvalumab 21 基因 panel）都使用**靶向免疫 panel**，
  **未测 TACSTD2 / CLDN4**（仅有 CD274）。
- **全基因组测到这两个上皮基因的肺癌系列**（GSE72094 芯片 442 例、GSE81089 RNA-seq
  199 例）是**手术切除队列**，只有**总生存 (OS)**，并非 ICI 应答。
- 唯一**同时**具备全基因组表达与 ICI 应答标签的数据是 **GSE91061**（nivolumab，
  Riaz 2017），但属于**黑色素瘤**（非肺）。

因此我们**严格按数据允许的范围报告**，并额外以 GSE91061 作为**同设计的黑色素瘤
ICI 基准**，完整演示预期的应答分析（明确标注为非肺）。

## 5. 定量结果
- **GSE91061（黑色素瘤，nivolumab，RECIST 应答）**：治疗前肿瘤，log2(FPKM+1)，
  Mann-Whitney（PR/CR vs PD）：TACSTD2 p=0.49、CLDN4 p=0.49、CD274 p=0.60——均无显著差异。
- **GSE72094（肺腺癌，n=398，总生存）**：中位数分组 log-rank：TACSTD2 p=0.41、
  CLDN4 p=0.58；与 CD274 的 Spearman ρ 均≈0。
- **GSE81089（NSCLC，n=196，总生存）**：TACSTD2 p=0.33、CLDN4 p=0.73；与 CD274 ρ≈0。

三个数据集中，TACSTD2/CLDN4 既未显著区分生存/应答，也与 PD-L1/CD274 基本无关。
图见 `analysis/*.png`，完整数值见 `analysis/analysis_results.json`。

## 6. 解读与局限
- **诚实的阴性结果**：窗口内不存在“肺 + ICI 应答标签 + 测到 TACSTD2/CLDN4”的数据集，
  故无法就 2015–2018 GEO 得出肺癌特异的 TACSTD2/CLDN4 与 ICI 应答关系；上述结果均为阴性，
  这是可报告的结论，而非流程失败。
- **分类说明**：基于文本的范围标记较为宽松——GSE91061 实为黑色素瘤（摘要仅提及 NSCLC）；
  GSE93157 为肺/头颈/黑色素瘤混合队列（65 例中 35 例为肺）；GSE101929/GSE102286
  （非裔 vs 欧裔 NSCLC）、GSE109010/GSE109020（SWI/SNF）仅因 PD-L1 / “免疫治疗”字样命中，
  并非 ICI 治疗研究。分析中均按正确组织学标签处理。
- **“已处理”范围**：仅下载公开的已处理补充文件与 series matrix（<2 GB），未抓取原始 SRA。
  分选细胞系列（如 GSE99531、GSE99254 的 T 细胞）虽含 TACSTD2/CLDN4 行但表达≈0（上皮基因在
  淋巴细胞中几乎不表达），且无逐样本应答标签。

## 7. 产物清单
见上文英文第 7 节；所有文件均位于
`notes/fable_geo_2015_2018/`、`scripts/fable_geo_2015_2018/`、`results/fable_geo_2015_2018/`。
