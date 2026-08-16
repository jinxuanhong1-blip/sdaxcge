# TACSTD2 (TROP2) & CLDN4 vs ICI Response — Open Human Bulk Lung Cohorts
# TACSTD2 (TROP2) 与 CLDN4 表达 vs 免疫检查点抑制剂（ICI）应答 —— 公开人类 bulk 肺癌队列

> Parallel slice. All outputs are under `notes/fable_ici_bulk/`, `scripts/fable_ici_bulk/`,
> `results/fable_ici_bulk/`. Real data, real statistics only. No FASTQ; processed matrices
> only (all < 2 GB).
>
> 平行任务分片。所有产物仅位于 `notes/fable_ici_bulk/`、`scripts/fable_ici_bulk/`、
> `results/fable_ici_bulk/`。仅使用真实数据与真实统计；不下载 FASTQ，仅用处理后的表达矩阵
> （均 < 2 GB）。

---

## 1. Scope / 研究范围

**EN.** We asked whether pre-/on-treatment tumor expression of the epithelial markers
**TACSTD2 (TROP2)** and **CLDN4** is associated with response to anti-PD-1/PD-L1 immune
checkpoint inhibitors (ICI) in **open, human, bulk lung** cohorts. Six GEO series were
requested: GSE126044, GSE135222, GSE136961, GSE166449, GSE93157, GSE207422 (bulk arm only).
Each GEO page was verified individually (platform, organism, assay type, response metadata,
and whether the two target genes are actually measured).

**中文.** 我们评估在**公开的人类 bulk 肺癌**队列中，上皮标志物 **TACSTD2 (TROP2)** 与
**CLDN4** 的肿瘤表达是否与抗 PD-1/PD-L1 免疫检查点抑制剂（ICI）疗效相关。任务指定 6 个
GEO 数据集：GSE126044、GSE135222、GSE136961、GSE166449、GSE93157、GSE207422（仅取 bulk 部分）。
每个 GEO 页面均逐一核实（平台、物种、检测类型、应答标注，以及两个目标基因是否被真实检测）。

---

## 2. Dataset eligibility / 数据集可用性

| GSE | Cohort / 队列 | Platform / 平台 | Response label / 应答标注 | TACSTD2 & CLDN4 measured? | Verdict / 结论 |
|-----|------|----------|-----------|----|------|
| **GSE126044** | 16 NSCLC, pre-tx anti-PD-1 | Illumina HiSeq 2500, bulk RNA-seq (raw counts) | responder / non-responder | **Yes** | **Usable / 可用** |
| **GSE135222** | 27 NSCLC, anti-PD-(L)1 | Illumina HiSeq 2500, bulk RNA-seq (TPM) | PFS event + PFS time (survival) | **Yes** | **Usable (survival) / 可用（生存）** |
| GSE136961 | 21 NSCLC, anti-PD-1 | Ion Torrent S5 XL, **Oncomine Immune Response** targeted panel (395 immune genes) | DCB / NDB | **No** | **Unusable / 不可用** |
| **GSE166449** | 22 lung cancer, pre-tx immunotherapy | Illumina HiSeq 2000, bulk RNA-seq (TPM) | responder / non-responder | **Yes** | **Usable / 可用** |
| GSE93157 | 65 pts (NSCLC/melanoma/HNSCC), anti-PD-1 | **NanoString PanCancer Immune 730** panel | RECIST + PFS | **No** | **Unusable / 不可用** |
| **GSE207422** | NSCLC neoadjuvant anti-PD-1 + chemo; **24 bulk** RNA-seq | Illumina NovaSeq 6000, bulk (log2TPM) | pathologic response MPR / NMPR | **Yes** | **Usable (bulk only) / 可用（仅 bulk）** |

**Why two are unusable / 两个数据集不可用的原因.**
GSE136961 (Oncomine Immune Response, 395 genes) and GSE93157 (NanoString PanCancer Immune,
730 genes) are **targeted immune-oncology panels**. A `grep` of their processed matrices for
`TACSTD|CLDN|TROP` returned **nothing** — the epithelial markers TACSTD2 and CLDN4 are simply
not on these panels, so the requested analysis cannot be performed on them regardless of
sample size. This is documented, not a failure to download (both files were downloaded and
inspected). See `notes/fable_ici_bulk/geo_verification.md`.

GSE136961（Oncomine 免疫应答，395 基因）与 GSE93157（NanoString PanCancer 免疫，730 基因）
均为**靶向免疫肿瘤 panel**。对其处理矩阵检索 `TACSTD|CLDN|TROP` 无任何命中——这两个 panel
根本不包含上皮标志物 TACSTD2 与 CLDN4，故无论样本量如何都无法进行本分析。此为“记录为不可用”，
而非下载失败（两个文件均已下载并检查）。

GSE207422 also contains a 175 MB scRNA-seq matrix; per the task ("bulk if present") we used
**only the bulk RNA-seq arm** (`GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz`, 5.4 MB) and
ignored the single-cell data.

GSE207422 另含 175 MB 单细胞矩阵；按任务要求（“如有 bulk 则用 bulk”），我们**仅使用其 bulk 部分**
（`GSE207422_NSCLC_bulk_RNAseq_log2TPM.txt.gz`，5.4 MB），忽略单细胞数据。

---

## 3. Methods / 方法

- **Downloads / 下载.** Processed supplementary matrices + `series_matrix` metadata from the
  NCBI GEO FTP. No FASTQ. Files staged in `/tmp/ici_bulk_data` (not committed; raw GEO data).
- **Gene IDs / 基因编号.** GSE135222 uses Ensembl IDs; TACSTD2 = `ENSG00000184292`,
  CLDN4 = `ENSG00000189143`. Others use HGNC symbols.
- **Normalization / 归一化.**
  - GSE126044 raw counts → **log2 CPM** (per-sample library-size scaling, `log2(CPM+1)`).
  - GSE135222 / GSE166449 TPM → `log2(TPM+1)`.
  - GSE207422 values are already log2TPM → used as provided.
- **Binary comparison / 二分类比较** (GSE126044, GSE166449, GSE207422):
  two-sided **Mann–Whitney U** (Wilcoxon rank-sum) between responders and non-responders;
  effect size reported as **AUC** = P(responder expr > non-responder expr) and rank-biserial
  correlation. Response definitions: responder/non-responder (GSE126044, GSE166449);
  **MPR / MPR(pCR) = responder vs NMPR = non-responder** (GSE207422 pathologic response).
- **Survival / 生存** (GSE135222): **Cox proportional-hazards** on continuous `log2(TPM+1)`
  expression for progression-free survival (event = PFS event, duration = PFS time), reporting
  HR per log2 unit with 95% CI; plus **median-split log-rank** test and Kaplan–Meier curves.
- **Multiple testing / 多重检验.** Benjamini–Hochberg FDR across all 8 primary tests
  (6 Mann–Whitney + 2 Cox).
- **Reproducibility / 可复现.** `scripts/fable_ici_bulk/analyze.py` (+ `geo_utils.py`).
  Per-sample values in `results/fable_ici_bulk/tables/per_sample_expression.csv`.

---

## 4. Results / 结果

### 4.1 Binary response (Mann–Whitney U) / 二分类应答

| Dataset | Gene | n resp / non | Median resp | Median non | AUC (resp>non) | p | FDR |
|---------|------|--------------|-------------|------------|----------------|---|-----|
| GSE126044 | TACSTD2 | 5 / 11 | 6.17 | 6.29 | 0.364 | 0.441 | 0.705 |
| GSE126044 | CLDN4 | 5 / 11 | 2.54 | 3.77 | 0.236 | 0.115 | 0.705 |
| GSE166449 | TACSTD2 | 7 / 15 | 2.12 | 1.78 | 0.619 | 0.407 | 0.705 |
| GSE166449 | CLDN4 | 7 / 15 | 1.21 | 1.33 | 0.514 | 0.945 | 0.945 |
| GSE207422 | TACSTD2 | 9 / 15 | 5.00 | 6.22 | 0.378 | 0.340 | 0.705 |
| GSE207422 | CLDN4 | 9 / 15 | 5.85 | 6.27 | 0.356 | 0.257 | 0.705 |

(Expression units: GSE126044 log2CPM; GSE166449 log2(TPM+1); GSE207422 log2TPM.)

### 4.2 Survival — GSE135222 PFS / 生存

| Gene | n (events) | Cox HR /log2 unit (95% CI) | Cox p | Median-split log-rank p |
|------|-----------|-----------------------------|-------|--------------------------|
| TACSTD2 | 27 (21) | 1.03 (0.82–1.31) | 0.776 | 0.430 |
| CLDN4 | 27 (21) | 1.05 (0.87–1.27) | 0.608 | 0.704 |

Full machine-readable tables: `results/fable_ici_bulk/tables/{binary_response_stats,survival_stats,per_sample_expression}.csv`.
Figures (boxplots + Kaplan–Meier): `results/fable_ici_bulk/figures/`.

### 4.3 Bottom line / 核心结论

**EN.** Across the four usable open bulk lung ICI cohorts, **neither TACSTD2 (TROP2) nor
CLDN4 tumor expression was significantly associated with ICI response or PFS.** All 8 primary
tests had raw p > 0.05, and every Benjamini–Hochberg FDR was > 0.70. Effect directions were
**inconsistent across cohorts** (e.g., TACSTD2 AUC 0.62 in GSE166449 but 0.36–0.38 in
GSE126044/GSE207422), which is not what a real, reproducible predictive biomarker looks like.
The Cox HRs for PFS were ~1.0 with CIs spanning 1, i.e. no PFS signal. Cohorts are small
(n = 16–27), so the study is underpowered to detect a small effect; but there is **no evidence**
here that either gene predicts ICI benefit in bulk lung tumors.

**中文.** 在 4 个可用的公开 bulk 肺癌 ICI 队列中，**TACSTD2 (TROP2) 与 CLDN4 的肿瘤表达
均与 ICI 应答或 PFS 无显著相关**。全部 8 项主要检验的原始 p 值均 > 0.05，BH 校正 FDR 均 > 0.70。
效应方向在不同队列间**不一致**（例如 TACSTD2 在 GSE166449 的 AUC 为 0.62，但在
GSE126044/GSE207422 为 0.36–0.38），这与一个真实、可复现的预测性生物标志物应有的表现不符。
GSE135222 中 PFS 的 Cox HR ≈ 1.0 且置信区间跨越 1，即无 PFS 信号。各队列样本量较小
（n = 16–27），检验效能不足以发现小效应；但现有证据**不支持**这两个基因可预测 bulk 肺肿瘤的
ICI 获益。

---

## 5. Caveats / 局限性

- Small cohorts (n = 16–27) → limited power; a modest true effect could be missed.
  样本量小（n = 16–27），效能有限，可能漏检中等真实效应。
- Response definitions differ (RECIST responder/non-responder; PFS; pathologic MPR/NMPR).
  各队列应答定义不同（RECIST 应答；PFS；病理 MPR/NMPR），不可直接混合。
- Mixed treatment lines/settings: GSE207422 is **neoadjuvant anti-PD-1 + chemotherapy**
  (combination), not ICI monotherapy; GSE126044/GSE166449/GSE135222 are advanced-stage
  anti-PD-(L)1. 治疗线与设置混杂：GSE207422 为**新辅助 抗 PD-1 + 化疗**联合方案，非单药 ICI。
- Bulk expression mixes tumor + stroma + immune compartments.
  bulk 表达为肿瘤 + 间质 + 免疫成分的混合信号。
- Two of the six requested series (GSE136961, GSE93157) do not measure the target genes and
  are documented as unusable, as instructed.
  6 个数据集中有 2 个（GSE136961、GSE93157）未检测目标基因，按要求记录为不可用。
- No cross-cohort meta-analysis / batch correction was performed given heterogeneous
  endpoints and units. 因终点与单位异质，未做跨队列 meta 分析或批次校正。

---

## 6. Files / 文件清单

```
notes/fable_ici_bulk/
  WRITEUP.md              <- this file / 本文件
  geo_verification.md     <- per-GEO verification + eligibility / 逐一核实与可用性
scripts/fable_ici_bulk/
  geo_utils.py            <- series_matrix metadata parser
  analyze.py              <- download-aware analysis (stats + figures)
results/fable_ici_bulk/
  tables/binary_response_stats.csv
  tables/survival_stats.csv
  tables/per_sample_expression.csv
  figures/GSE126044_{TACSTD2,CLDN4}.png
  figures/GSE166449_{TACSTD2,CLDN4}.png
  figures/GSE207422_{TACSTD2,CLDN4}.png
  figures/GSE135222_{TACSTD2,CLDN4}_KM.png
```

Reproduce / 复现:
```bash
pip3 install pandas scipy statsmodels matplotlib lifelines openpyxl
# download the 6 processed matrices + series_matrix files into /tmp/ici_bulk_data
python3 scripts/fable_ici_bulk/analyze.py
```
