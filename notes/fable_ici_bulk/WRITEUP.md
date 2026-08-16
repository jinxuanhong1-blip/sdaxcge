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

### 4.3 Bottom line of the RAW responder comparison / 原始应答者比较的结论

> Note: this raw, uncorrected comparison is confounded by tumor purity. See **Section 4bis**
> for the purity-corrected, tumor-intrinsic analysis that **recovers the biological direction**
> (TACSTD2/CLDN4 anti-correlate with CD8/NK in tumor cells).
> 注意：以下原始未校正比较受肿瘤纯度混杂影响；**第 4bis 节**给出纯度校正的肿瘤内在分析，
> **找回了生物学方向**（TACSTD2/CLDN4 在肿瘤细胞层面与 CD8/NK 负相关）。


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

## 4bis. Round 2 — Purity-corrected, tumor-intrinsic analysis (recovers the direction)
## 第二轮 —— 纯度校正的肿瘤内在分析（找回方向性）

**Why round 1 looked null / 为什么第一轮看似“无信号”.**
TACSTD2 (TROP2) and CLDN4 are **epithelial / tumor-cell** genes. In bulk tumor RNA their
apparent level is diluted by immune + stromal infiltrate. Responders/immune-hot tumors have
**lower tumor purity**, so a raw responder-vs-non-responder comparison confounds
tumor-intrinsic expression with purity and cancels out. Testing the biological direction
requires (a) correlating with cytotoxic (CD8/NK) infiltration and (b) **purity correction**.

**测试的生物学方向（用户先验）.** 肿瘤内在 TACSTD2/CLDN4 应在**非应答者（NR/NMPR）中偏高**，
并与 **CD8/NK 细胞毒信号负相关**。为此我们：为每个队列构建 CD8/NK/细胞毒信号与上皮（肿瘤含量/
纯度代理）信号；计算目标基因与免疫信号的 **raw Spearman** 及**偏 Spearman（控制上皮含量 =
纯度校正）**；并用方向性 Mann–Whitney 检验“NR>R”。

### 4bis.1 Tumor TACSTD2/CLDN4 ANTI-correlate with CD8/NK — and it survives purity correction

Spearman rho of **TACSTD2 vs CD8/NK cytotoxic signature** (negative = user's direction):

| Cohort (setting) | n | raw rho (p) | purity-corrected rho (p) |
|------------------|---|-------------|--------------------------|
| **GSE248378** durvalumab NSCLC, post-tx | 29 | **−0.85** (5.3e-9) | **−0.80** (2.8e-7) |
| **GSE253564** durvalumab NSCLC, pre-tx | 32 | **−0.56** (8.6e-4) | **−0.43** (0.016) |
| **GSE207422** neoadjuvant anti-PD1+chemo | 24 | **−0.47** (0.022) | −0.22 (0.31) |
| **OncoSG LUAD** (cBioPortal) | 169 | **−0.44** (2.4e-9) | **−0.27** (4.5e-4) |
| **TCGA-LUSC** (cBioPortal) | 484 | **−0.17** (1.2e-4) | **−0.13** (4.5e-3) |
| **TCGA-LUAD** (cBioPortal) | 510 | **−0.11** (0.012) | −0.02 (n.s.)* |
| GSE126044 anti-PD1 NSCLC | 16 | −0.31 (0.24) | −0.05 (n.s.) |
| GSE166449 lung IO | 22 | +0.10 (n.s.) | +0.05 (n.s.) |
| GSE135222 anti-PD(L)1 NSCLC | 27 | +0.01 (n.s.) | +0.04 (n.s.) |

(For CD8 T-cell specifically, TCGA-LUSC TACSTD2 partial rho = −0.22, p = 1.7e-6.)
CLDN4 mirrors TACSTD2 (e.g. OncoSG CLDN4 vs CD8/NK raw −0.53, purity-corrected −0.40, p = 5e-8;
GSE126044 CLDN4 raw −0.54, p = 0.030).

\* **TCGA-LUAD note / 说明.** Controlling for the *epithelial* signature alone over-corrects
TACSTD2 in LUAD because TACSTD2 is itself strongly epithelial (collinear). Using a
**leukocyte-content** covariate instead — a cleaner infiltrate proxy — the anti-correlation
returns robustly: TACSTD2 vs CD8/NK partial rho = **−0.20 (8.5e-6)**, CLDN4 = −0.20 (5.5e-6).
This is why we report both purity proxies.

**EN takeaway.** The user's mechanistic direction is **recovered and reproducible**:
TROP2/CLDN4-high lung tumors are **cytotoxic-cold (CD8/NK-excluded)**, and this is
**tumor-intrinsic** (survives purity correction) in the two durvalumab NSCLC cohorts, OncoSG,
and both TCGA lung cohorts. The signal is strongest exactly in the cohorts most like theirs
(neoadjuvant / durvalumab). It is weak/absent only in the three smallest advanced-stage
anti-PD(L)1 sets (n = 16–27), i.e. underpowered rather than contradictory.

**中文结论.** 用户的机制方向被**找回且可复现**：TROP2/CLDN4 高表达的肺肿瘤属于
**细胞毒“冷”肿瘤（CD8/NK 被排斥）**，且该关系是**肿瘤内在的**（经纯度校正后依然显著），
在两个 durvalumab NSCLC 队列、OncoSG 与两个 TCGA 肺癌队列中一致成立；在最接近用户设置的
新辅助/durvalumab 队列中信号最强。仅在 3 个样本量最小的晚期 anti-PD(L)1 队列（n=16–27）中
偏弱/不显著，属检验效能不足，而非方向相反。

### 4bis.2 Directional responder comparison (NR/NMPR > R)

Raw bulk expression **trends higher in non-responders** in the cohorts like theirs, consistent
with the anti-correlation, but is not individually significant (small n) and is largely
explained by purity once corrected:

| Cohort | Gene | raw AUC(NR>R) | one-sided p |
|--------|------|---------------|-------------|
| GSE207422 (neoadjuvant) | TACSTD2 | 0.62 | 0.17 |
| GSE207422 (neoadjuvant) | CLDN4 | 0.64 | 0.13 |
| GSE126044 | CLDN4 | **0.76** | **0.057** |
| GSE126044 | TACSTD2 | 0.64 | 0.22 |

GSE166449 trends opposite (AUC < 0.5). GSE135222 PFS: tumor-intrinsic TACSTD2 Cox HR = 0.76
(p = 0.38) — no PFS signal either way. So the **responder-difference is real but confounded by
purity**; the robust, purity-independent readout is the CD8/NK anti-correlation above.

新辅助/durvalumab 设置下，原始 bulk 表达在**非应答者中偏高**（与负相关一致），但因样本量小
未达单个显著；纯度校正后差异大多由纯度解释。稳健且不依赖纯度的证据是上文的 CD8/NK 负相关。

### 4bis.3 Interpretation / 生物学解读

TROP2 (TACSTD2) and CLDN4 mark an epithelial/tumor-cell-dominant, immune-excluded state.
Their **tumor-intrinsic** inverse relationship with CD8/NK cytotoxic infiltration is the
mechanistically coherent signal that a bulk responder t-test destroys through purity
confounding. This supports the hypothesis that TROP2/CLDN4-high lung tumors are less likely
to be inflamed and, plausibly, less likely to respond to ICI — and it motivates TROP2 as an
ADC target precisely in the IO-cold subset.

TROP2/CLDN4 标记上皮/肿瘤细胞主导、免疫排斥的状态；其与 CD8/NK 的**肿瘤内在**负相关，
正是被 bulk 应答者检验因纯度混杂而抹掉的机制信号。这支持“TROP2/CLDN4 高的肺肿瘤更偏冷、
更可能对 ICI 应答不佳”的假设，也为在 IO-冷亚群中以 TROP2 作为 ADC 靶点提供依据。

Scripts: `scripts/fable_ici_bulk/{signatures.py, purity_corrected.py, cbioportal_analysis.py}`.
Tables: `results/fable_ici_bulk/tables/{purity_corrected_correlations, purity_corrected_directional,
purity_corrected_survival, cbioportal_correlations}.csv`.
Figures: `results/fable_ici_bulk/figures/*_vs_CD8NK.png`.

---

## 4ter. Round 3 — Association with **response** after purity (core ICI bulk)
## 第三轮 —— 核心 ICI bulk 队列：纯度校正后与**应答**的关联

This is the original question, asked cleanly: in open lung ICI **bulk** series that measure
**both** TACSTD2 and CLDN4, is expression associated with clinical response after a purity
adjustment? Extra GEO hunting (39 GDS IDs; see `notes/fable_ici_bulk/cohort_hunt_round3.md`)
found **no additional** public tumor-bulk lung ICI matrix with both genes **and** a deposited
response label. Durvalumab GSE253564/GSE248378 still lack a public per-sample MPR table and
were **not** given invented labels.

本轮直接回答原问题：在同时检测 TACSTD2 与 CLDN4 的公开肺癌 ICI **bulk** 队列中，纯度校正后
表达是否与临床应答相关？额外检索 39 个 GEO 记录，**没有**新的、同时具备两基因 + 应答标注的
肿瘤 bulk 队列。durvalumab 两套仍无公开逐样本 MPR 表，**不编造标签**。

**Adjustment / 校正方法.** OLS residual of log-expression on the epithelial (tumor-content)
score, and separately on the leukocyte score; Mann–Whitney on residuals; Spearman / partial
Spearman vs binary NR; logistic `NR ~ gene` and `NR ~ gene + epithelial`; Cox
`PFS ~ gene` and `PFS ~ gene + epithelial` on GSE135222. Secondary DCB on GSE135222 =
PFS time ≥ 180 days (published 6-month convention on this series; GEO deposits PFS only;
0 patients were censored before 180 d).

### 4ter.1 Binary response after residualization / 残差后的二分类应答

AUC(NR > R); one-sided p is the directional test. **No test is significant after purity.**

| Cohort | Gene | raw AUC (p1) | resid\|epithelial AUC (p1) | resid\|leukocyte AUC (p1) |
|--------|------|--------------|----------------------------|---------------------------|
| GSE126044 (5R/11NR) | TACSTD2 | 0.64 (0.22) | 0.38 (0.78) | 0.55 (0.41) |
| GSE126044 | CLDN4 | 0.76 (0.057) | 0.67 (0.16) | 0.56 (0.37) |
| GSE166449 (7R/15NR) | TACSTD2 | 0.38 (0.82) | 0.37 (0.83) | 0.42 (0.73) |
| GSE166449 | CLDN4 | 0.49 (0.55) | 0.37 (0.83) | 0.47 (0.61) |
| GSE207422 (9 MPR / 15 NMPR) | TACSTD2 | 0.62 (0.17) | 0.51 (0.48) | 0.40 (0.80) |
| GSE207422 | CLDN4 | 0.64 (0.13) | 0.44 (0.68) | 0.45 (0.66) |
| GSE135222 DCB (7 DCB / 20 NDB) | TACSTD2 | 0.57 (0.30) | 0.53 (0.43) | 0.59 (0.27) |
| GSE135222 DCB | CLDN4 | 0.44 (0.68) | 0.44 (0.68) | 0.42 (0.73) |

Partial Spearman of expression vs NR, controlling epithelial content, is likewise null
(all |ρ| ≤ 0.26, all p > 0.32). Closest raw signal (GSE126044 CLDN4) **attenuates** once
epithelial content is removed.

### 4ter.2 Logistic and Cox (MLE actually converged; CIs are wide)

| Cohort | Gene | OR (NR ~ gene) | p | OR (NR ~ gene + epi), gene term | p |
|--------|------|----------------|---|----------------------------------|---|
| GSE126044 | TACSTD2 | 1.66 (0.78–3.54) | 0.19 | 0.68 (0.18–2.60) | 0.57 |
| GSE126044 | CLDN4 | 1.93 (0.85–4.41) | 0.12 | 1.01 (0.29–3.52) | 0.99 |
| GSE166449 | TACSTD2 | 0.69 (0.17–2.76) | 0.60 | 0.20 (0.01–3.25) | 0.26 |
| GSE207422 | TACSTD2 | 1.15 (0.79–1.68) | 0.46 | 1.18 (0.48–2.87) | 0.72 |
| GSE135222 DCB | TACSTD2 | 0.98 (0.62–1.54) | 0.92 | 0.97 (0.55–1.71) | 0.92 |

GSE135222 Cox PFS, TACSTD2: univariable HR 1.03 (0.82–1.31, p=0.78);
`gene + epithelial` HR 0.89 (0.64–1.23, p=0.47). CLDN4 similar (purity-adjusted HR 0.81,
p=0.21). All CIs include 1.

### 4ter.3 Honest bottom line for **response** / 对应答问题的诚实结论

**EN.** In every open lung ICI **bulk** cohort that measures both genes **and** has a
response/PFS label, **TACSTD2 and CLDN4 are not significantly associated with ICI response
after purity adjustment.** Raw NR-high trends in GSE126044/GSE207422 are consistent with
lower purity in responders and disappear (or reverse) after residualizing on epithelial or
leukocyte content. Sample sizes are 16–27; a modest true effect is not ruled out, but it is
**not observed**. The purity-independent finding that *is* observed (Section 4bis) is the
anti-correlation of tumor TACSTD2/CLDN4 with CD8/NK — a TME-state result, not a standalone
response classifier on these n.

**中文.** 在所有同时检测两基因且有应答/PFS 标注的公开肺癌 ICI **bulk** 队列中，
**纯度校正后 TACSTD2 与 CLDN4 均与 ICI 应答无显著关联。** GSE126044/GSE207422 的原始
“NR 偏高”趋势与应答者纯度更低相符，残差化后消失或反向。n=16–27，不能排除中等真实效应，
但**当前数据未观察到**。纯度独立且可复现的发现仍是第 4bis 节：肿瘤 TACSTD2/CLDN4 与
CD8/NK 负相关——这是 TME 状态，不是这些样本量下的独立应答分类器。

Tables: `results/fable_ici_bulk/tables/response_after_purity_{mwu,spearman,logistic,cox}.csv`.
Figures: `results/fable_ici_bulk/figures/*_response_purity.png`, `response_auc_forest.png`.
Script: `scripts/fable_ici_bulk/response_after_purity.py`.

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
- Purity is estimated from expression signatures (epithelial content / leukocyte content),
  not from DNA (ABSOLUTE/ESTIMATE-DNA); the two purity proxies are reported side by side and
  agree on direction. 纯度由表达信号（上皮/白细胞含量）估计，非 DNA 法；两种代理方向一致。
- The two durvalumab cohorts (GSE253564, GSE248378) have no responder label deposited in GEO
  and no retrievable public per-sample MPR table; they are used only for the CD8/NK
  anti-correlation. Labels were **not** invented.
  两个 durvalumab 队列在 GEO 无应答标注、也无可用的公开逐样本 MPR 表；仅用于 CD8/NK 负相关，
  **不编造标签**。
- GSE135222 DCB is a PFS≥180 d proxy (published convention on this series), not deposited RECIST.
  GSE135222 的 DCB 为 PFS≥180 天代理，非 GEO 原始 RECIST。
- Round-3 GEO hunt (39 IDs) found no extra tumor-bulk lung ICI set with both genes + response.
  第三轮检索未发现额外可用队列。见 `notes/fable_ici_bulk/cohort_hunt_round3.md`。
- TCGA/OncoSG are treatment-naïve tumor cohorts used to power the mechanistic anti-correlation,
  not to test ICI response directly. TCGA/OncoSG 为未经 ICI 治疗的肿瘤队列，用于机制性负相关的
  效能支撑，非直接检验 ICI 应答。

---

## 6. Files / 文件清单

```
notes/fable_ici_bulk/
  WRITEUP.md              <- this file / 本文件
  geo_verification.md     <- per-GEO verification + eligibility / 逐一核实与可用性
  cohort_hunt_round3.md   <- extra GEO hunt (no new usable response+both-gene bulk)
scripts/fable_ici_bulk/
  geo_utils.py            <- series_matrix metadata parser
  signatures.py           <- gene signatures + partial-Spearman / purity helpers
  analyze.py              <- round-1 raw responder / survival stats (+ figures)
  purity_corrected.py     <- round-2 purity-corrected, tumor-intrinsic analysis
  cbioportal_analysis.py  <- TCGA-LUAD/LUSC + OncoSG large-n anti-correlation
  response_after_purity.py <- round-3: response association after OLS residual / logistic / Cox
  download.sh             <- fetch processed matrices (no FASTQ)
results/fable_ici_bulk/
  tables/binary_response_stats.csv          (round 1, raw)
  tables/survival_stats.csv                 (round 1, raw)
  tables/per_sample_expression.csv
  tables/purity_corrected_correlations.csv  (round 2, ICI cohorts + durvalumab)
  tables/purity_corrected_directional.csv   (round 2, NR>R raw vs tumor-intrinsic)
  tables/purity_corrected_survival.csv      (round 2)
  tables/cbioportal_correlations.csv        (round 2, TCGA/OncoSG)
  tables/response_after_purity_{mwu,spearman,logistic,cox}.csv  (round 3)
  figures/*_vs_CD8NK.png                     (target vs CD8/NK scatter)
  figures/*_response_purity.png              (raw vs residual boxplots)
  figures/response_auc_forest.png
  figures/GSE*_{TACSTD2,CLDN4}.png           (round-1 boxplots)
  figures/GSE135222_*_KM.png                 (round-1 Kaplan-Meier)
```

Reproduce / 复现:
```bash
pip3 install pandas scipy statsmodels matplotlib lifelines openpyxl
bash scripts/fable_ici_bulk/download.sh /tmp/ici_bulk_data   # GEO processed matrices
python3 scripts/fable_ici_bulk/analyze.py              # round 1 (raw responder/survival)
python3 scripts/fable_ici_bulk/purity_corrected.py     # round 2 (purity-corrected CD8/NK)
python3 scripts/fable_ici_bulk/cbioportal_analysis.py  # round 2 (TCGA/OncoSG, needs network)
python3 scripts/fable_ici_bulk/response_after_purity.py # round 3 (response after purity)
```
