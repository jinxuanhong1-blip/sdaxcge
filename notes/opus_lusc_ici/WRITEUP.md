# TACSTD2 / CLDN4 in public LUSC ICI transcriptomes

**Slice paths only:** `notes/opus_lusc_ici/`, `scripts/opus_lusc_ici/`, `results/opus_lusc_ici/`.  
**Rule:** LUSC/LSCC only for the primary claims. Histology is split, not mixed. Every number below is computed from the public matrices; nothing is filled in.

Reproduce: `bash scripts/opus_lusc_ici/run_all.sh` (Python 3, packages listed in `scripts/opus_lusc_ici/requirements.txt`). Raw downloads stay in `/tmp/opus_lusc_ici_cache` and are not committed.

---

## English

### Question

Do **TACSTD2** (TROP2) and **CLDN4** associate with immune-checkpoint inhibitor (ICI) **response** or with **immune programmes** in **lung squamous-cell carcinoma (LUSC / LSCC)**? Public ICI RNA-seq is almost always NSCLC-mixed. This slice keeps LUSC separate from LUAD / non-squamous disease.

### What was included, and what was not

Seven GEO series with tumour RNA and an ICI context were downloaded. Two further series were inspected and **excluded**:

| Accession | Why used / unused |
|---|---|
| GSE190265 (Dijon France-3) | Pathology histology + anti-PD-1 monotherapy + PFS. **Included.** |
| GSE190266 (Dijon France-4) | Pathology histology + anti-PD-1 + 6-month PFS. **Included for CLDN4 only** — the deposited 16,383-gene matrix has **no TACSTD2 row**. |
| GSE207422 | Neoadjuvant anti-PD-1 + chemotherapy; bulk log2 TPM + MPR. **Included** (pre-treatment bulk only). |
| GSE283829 | Advanced ICI; RECIST best response; pathology SqCC. **Included.** |
| GSE253564 | Neoadjuvant durvalumab ± SBRT, pre-treatment FPKM. Pathology squamous. **Immune-context only** (no public response label). |
| GSE135222 (SMC) | Anti-PD-1/PD-L1 + PFS. **No histology in GEO.** Used only in the sensitivity set after a high-confidence transcriptomic LUSC call. |
| GSE126044 (Yonsei) | Anti-PD-1 + responder/non-responder. **No histology in GEO.** Same rule as SMC. |
| GSE93157 (Prat NanoString) | Has “SQUAMOUS LUNG CANCER” labels, but the 770-gene panel **lacks TACSTD2 and CLDN4**. Unused. |
| GSE136961 | Oncomine Immune Response panel (395 genes). Not used (panel, incomplete raw data). |

TCGA-LUSC (UCSC Xena GDC STAR TPM, n = 501 primary tumours) is a **treatment-naive surgical** reference. It is used to (i) train / calibrate the histology score and (ii) replicate the **immune-axis** correlations. It is **not** an ICI-response cohort.

### Histology split (not mixed)

Pathology labels in GEO are never overwritten.

For SMC and Yonsei a squamous-minus-adenocarcinoma marker score was computed on within-sample percentile ranks over a 17,226-gene universe shared with TCGA. High-confidence inferred calls use gates **Δ ≥ 0.20 → LUSC** and **Δ ≤ −0.20 → non-LUSC**; the interval is left **indeterminate** and is excluded from LUSC analyses.

Hold-out on pathology-annotated ICI samples (classifier not used for their labels):

| Rule | n called | accuracy | sensitivity (LUSC) | specificity | indeterminate |
|---|---:|---:|---:|---:|---:|
| High-confidence gates | 121 | 0.983 | 0.976 | 0.988 | 38 |
| Midpoint forced call | 159 | 0.881 | 0.962 | 0.840 | 0 |

Primary LUSC counts (pathology, pre-treatment):

| Cohort | LUSC n | with binary endpoint | benefit / no-benefit | TACSTD2 | CLDN4 |
|---|---:|---:|---|---:|---:|
| FRANCE3 | 11 | 11 | 5 / 6 | 11 | 11 |
| FRANCE4 | 11 | 10 | 1 / 9 | 0 | 11 |
| NEOCHEMO | 12 | 12 | 4 / 8 | 12 | 12 |
| PLA | 9 | 9 | 2 / 7 | 9 | 9 |
| DURVART | 10 | 0 | — | 10 | 10 |
| **Total pathology LUSC** | **53** | **42** | **12 / 30** | **42** | **53** |

High-confidence inferred LUSC added in the sensitivity set only: SMC 11 (1 benefit / 10 no-benefit), Yonsei 7 (2 / 5).

### Endpoints

- **Advanced monotherapy (France-3/4, SMC):** durable benefit = progression-free at the **6-month landmark**. Patients censored before 6 months are left missing.
- **Neoadjuvant chemo-ICI (GSE207422):** benefit = **MPR** (including pCR) vs NMPR.
- **PLA (GSE283829):** benefit = **CR** vs SD/PD (the GEO field `disease stage` stores RECIST best response).
- **DURVART:** no public response field.

These endpoints are **not interchangeable**. The meta-analysis therefore pools a **direction of association** (Hedges’ *g* of expression in benefit vs no-benefit), not a single clinical definition.

A cohort enters a two-group test only if **both arms have ≥ 2 samples with the gene measured**. France-4 LUSC (1 benefit) is therefore **dropped from the binary response tests** (it remains in the CLDN4 PFS Cox model).

### Primary result — ICI response in pathology LUSC

Random-effects (DerSimonian–Laird) meta-analysis of Hedges’ *g* (benefit − no-benefit). Permutation *p* values randomise benefit labels **within cohort** (10,000 permutations).

| Gene | k cohorts | n (benefit / no) | *g* | 95% CI | *p* | *I*² |
|---|---:|---|---:|---|---:|---:|
| TACSTD2 | 3 | 32 (11 / 21) | −0.29 | −1.05 to 0.47 | 0.46 | 0% |
| CLDN4 | 3 | 32 (11 / 21) | +0.03 | −0.97 to 1.03 | 0.96 | 39% |

Housekeeping controls on the same 32 samples: ACTB *g* = −0.002 (*p* = 1.00), GAPDH *g* = −0.18 (*p* = 0.65), PPIA *g* = +0.03 (*p* = 0.94).

Per-cohort TACSTD2 (pathology LUSC): France-3 *g* = −0.99 (MWU *p* = 0.18, permutation *p* = 0.11); neoadjuvant *g* = +0.01; PLA *g* = +0.28. No cohort is individually significant. Leave-one-cohort-out does not create a significant pooled effect.

Adding high-confidence inferred LUSC (sensitivity set) does not change the conclusion: TACSTD2 *g* = −0.20 (k = 4, n = 39, *p* = 0.58); CLDN4 *g* = −0.02 (k = 4, n = 39, *p* = 0.95).

**Pathology non-LUSC contrast** (not the primary claim): TACSTD2 *g* = +0.17 (k = 3, n = 38, *p* = 0.66); CLDN4 *g* = +0.46 (k = 4, n = 86, *p* = 0.17). Mixing LUSC with LUAD would have hidden that the (already null) LUSC estimate is not the same object as the non-LUSC estimate.

### PFS (Cox, per +1 SD within the analysed subset)

Pathology LUSC, France-3 (n = 11, 7 events): TACSTD2 HR = 4.08 (95% CI 0.82–20.3, *p* = 0.086); CLDN4 HR = 1.74 (0.69–4.40, *p* = 0.24).  
France-4 LUSC CLDN4 (n = 11, 9 events): HR = 1.03 (0.51–2.08, *p* = 0.93).  
Pooled within-cohort z-scored CLDN4 (n = 22, 16 events): HR = 1.24 (0.77–2.01, *p* = 0.37).

These intervals are wide. They are compatible with no effect and with a harmful effect of higher TACSTD2; they are **not** evidence of benefit.

### Immune context

Signature scores are the mean within-cohort z-score of the genes present (gene lists and sources: `scripts/opus_lusc_ici/signatures.json`). Spearman ρ is meta-analysed with Fisher *z* (random effects).

**Pathology LUSC ICI (primary set):**

| Gene | Signature | k | n | ρ | 95% CI | *p* | BH *q* |
|---|---|---:|---:|---:|---|---:|---:|
| TACSTD2 | EPITHELIAL | 4 | 42 | +0.72 | 0.50 to 0.85 | 6.8×10⁻⁷ | 1.8×10⁻⁵ |
| CLDN4 | EPITHELIAL | 5 | 53 | +0.69 | 0.48 to 0.82 | 2.1×10⁻⁷ | 1.1×10⁻⁵ |
| TACSTD2 | CHECKPOINT | 4 | 42 | −0.41 | −0.66 to −0.07 | 0.018 | 0.11 |
| TACSTD2 | CD8_T_CELL | 4 | 42 | −0.32 | −0.60 to 0.03 | 0.068 | 0.20 |
| TACSTD2 | CYT | 4 | 42 | −0.32 | −0.59 to 0.03 | 0.073 | 0.21 |
| TACSTD2 | IFNG_6 | 4 | 42 | +0.04 | −0.31 to 0.38 | 0.82 | 0.91 |

After a partial Spearman that removes the epithelial score, the ICI LUSC immune associations are smaller and are not claimed as independent of tumour-cell content.

**TCGA-LUSC (n = 501, treatment-naive) — immune-axis replication, not ICI response:**

| Gene | Signature | ρ | *p* | BH *q* | partial ρ (epithelial) |
|---|---|---:|---:|---:|---:|
| TACSTD2 | CD8_T_CELL | −0.285 | 7.8×10⁻¹¹ | 7.0×10⁻¹⁰ | −0.231 |
| TACSTD2 | TIS_18 | −0.247 | 2.0×10⁻⁸ | 7.0×10⁻⁸ | −0.202 |
| TACSTD2 | IFNG_6 | −0.247 | 2.1×10⁻⁸ | 7.0×10⁻⁸ | −0.206 |
| TACSTD2 | CYT | −0.205 | 3.6×10⁻⁶ | 9.2×10⁻⁶ | −0.153 |
| TACSTD2 | EPITHELIAL | +0.246 | 2.3×10⁻⁸ | 7.0×10⁻⁸ | — |
| CLDN4 | EPITHELIAL | +0.515 | 2.5×10⁻³⁵ | 4.5×10⁻³⁴ | — |
| CLDN4 | CD8_T_CELL | −0.115 | 0.010 | 0.020 | +0.040 (NS after epithelial) |

So: in LUSC, both genes track an **epithelial / tumour-cell** programme. TACSTD2 additionally shows a **modest inverse** correlation with CD8 / IFN-γ / cytolytic scores in TCGA-LUSC that partly survives epithelial adjustment. The same direction is visible but **not significant after BH** in the small ICI LUSC set. CLDN4’s weak inverse immune correlations in TCGA **collapse** after epithelial adjustment.

### Verification (what was actually run)

1. **Histology hold-out** on pathology ICI samples (table above). Primary analyses use pathology labels only.
2. **Permutation** of benefit labels within cohort (10,000).
3. **Housekeeping genes** (ACTB, GAPDH, PPIA) as response negative controls — all null.
4. **Leave-one-cohort-out** random-effects meta — still null.
5. **Sensitivity** adding high-confidence inferred LUSC — still null.
6. **Histology contrast** (pathology non-LUSC) reported separately, never pooled with LUSC for the primary claim.
7. **Partial Spearman** vs epithelial score.
8. **TCGA-LUSC** replication of the immune axis only.
9. **Gene coverage:** France-4 has no TACSTD2; that cohort is omitted from every TACSTD2 test rather than imputed.
10. **Processed size:** results tree is well under 2 GB (see `results/opus_lusc_ici/manifest/processed_size.json`).

### Interpretation (kept inside the data)

Public LUSC-on-ICI RNA-seq is **small**. With 11 benefit events across three pathology-confirmed cohorts, this slice **does not support** TACSTD2 or CLDN4 as ICI-response biomarkers in LUSC. The point estimates are near zero and the confidence intervals are wide.

The **immune** result is more stable: both genes are epithelial-linked; TACSTD2 is weakly anti-correlated with T-cell / IFN-γ programmes in a 501-tumour LUSC reference. That is a tumour-intrinsic / purity-adjacent association, not an ICI-outcome claim.

ADC-relevant surface targets (TROP2, claudin-4) can be highly expressed in LUSC without that expression marking ICI benefit in the public record assembled here.

### Limitations

- n is small; several cohorts have 2 benefit events.
- Endpoints differ (6-month PFS vs MPR vs CR).
- Platforms differ (TPM, log2 TPM, counts→logCPM, FPKM); tests are within-cohort.
- France-4 is a reduced gene panel.
- SMC / Yonsei histology is inferred; they are sensitivity-only.
- No PD-L1 IHC or TMB was uniformly available.
- Bulk RNA cannot assign TACSTD2/CLDN4 to tumour vs epithelium vs contaminant.

### Files

- Scripts: `scripts/opus_lusc_ici/`
- Tables: `results/opus_lusc_ici/tables/`
- Figures: `results/opus_lusc_ici/figures/`
- Download manifest (URL, bytes, SHA-256): `results/opus_lusc_ici/manifest/raw_downloads.json`

---

## 中文

### 问题

在**肺鳞癌（LUSC / LSCC）**中，**TACSTD2**（TROP2）和 **CLDN4** 是否与免疫检查点抑制剂（ICI）**疗效**或**免疫程序**相关？公开 ICI 转录组几乎都是 NSCLC 混杂队列。本切片把 LUSC 与腺癌 / 非鳞癌**拆开**，不混合统计。

### 纳入与排除

下载了 7 个带肿瘤 RNA 且有 ICI 背景的 GEO 系列。另检查并**排除**：

| 登录号 | 使用 / 不用的原因 |
|---|---|
| GSE190265（Dijon France-3） | 病理组织学 + 抗 PD-1 单药 + PFS。**纳入。** |
| GSE190266（Dijon France-4） | 病理组织学 + 抗 PD-1 + 6 个月 PFS。沉积矩阵 16,383 基因**没有 TACSTD2**，仅用于 CLDN4。 |
| GSE207422 | 新辅助抗 PD-1 + 化疗；bulk log2 TPM + MPR。**纳入**（仅治疗前 bulk）。 |
| GSE283829 | 晚期 ICI；RECIST 最佳疗效；病理 SqCC。**纳入。** |
| GSE253564 | 新辅助 durvalumab ± SBRT。有鳞癌标注，**无公开疗效字段**，只做免疫背景。 |
| GSE135222（SMC） | 抗 PD-1/PD-L1 + PFS。**GEO 无组织学。**仅在高置信转录组 LUSC 判定后进入敏感性集。 |
| GSE126044（延世） | 抗 PD-1 + 应答/无应答。**GEO 无组织学。**规则同 SMC。 |
| GSE93157（Prat NanoString） | 有“SQUAMOUS LUNG CANCER”标签，但 770 基因面板**不含 TACSTD2、CLDN4**。不用。 |
| GSE136961 | Oncomine 免疫面板。不用。 |

TCGA-LUSC（Xena GDC STAR TPM，n = 501 原发瘤）是**未接受 ICI 的手术标本**，只用于（i）校准组织学评分，（ii）重复**免疫轴**相关。它**不是** ICI 疗效队列。

### 组织学拆分（不混合）

GEO 病理标签从不覆盖。

SMC / 延世：在与 TCGA 共享的 17,226 基因宇宙上计算样本内百分位，再做鳞癌标记均值 − 腺癌标记均值。高置信门控为 **Δ ≥ 0.20 → LUSC**、**Δ ≤ −0.20 → 非 LUSC**；中间为**不确定**，不进入 LUSC 分析。

在**有病理标注**的 ICI 样本上做外推检查（这些样本的标签仍用病理，不用分类器）：

| 规则 | 已判定 n | 准确率 | LUSC 灵敏度 | 特异度 | 不确定 |
|---|---:|---:|---:|---:|---:|
| 高置信门控 | 121 | 0.983 | 0.976 | 0.988 | 38 |
| 中点强制分类 | 159 | 0.881 | 0.962 | 0.840 | 0 |

主要分析集（病理确诊 LUSC，治疗前）：

| 队列 | LUSC n | 有二分类终点 | 获益 / 未获益 | TACSTD2 | CLDN4 |
|---|---:|---:|---|---:|---:|
| FRANCE3 | 11 | 11 | 5 / 6 | 11 | 11 |
| FRANCE4 | 11 | 10 | 1 / 9 | 0 | 11 |
| NEOCHEMO | 12 | 12 | 4 / 8 | 12 | 12 |
| PLA | 9 | 9 | 2 / 7 | 9 | 9 |
| DURVART | 10 | 0 | — | 10 | 10 |
| **病理 LUSC 合计** | **53** | **42** | **12 / 30** | **42** | **53** |

敏感性集额外加入高置信推断 LUSC：SMC 11（1 / 10），延世 7（2 / 5）。

### 终点

- **晚期单药（France-3/4、SMC）：** 6 个月里程碑无进展 = 获益；6 个月前删失则缺失。
- **新辅助（GSE207422）：** MPR（含 pCR）vs NMPR。
- **PLA：** CR vs SD/PD（GEO 字段 `disease stage` 实为 RECIST）。
- **DURVART：** 无公开疗效。

终点**不可互换**。荟萃分析合并的是表达量在获益 vs 未获益之间的 **Hedges’ *g* 方向**，不是同一临床定义。

两组检验要求**每臂至少 2 例测到该基因**。France-4 LUSC 仅 1 例获益，因此**不进入二分类疗效检验**（仍进入 CLDN4 的 PFS Cox）。

### 主要结果 — 病理 LUSC 的 ICI 疗效

随机效应（DerSimonian–Laird）合并 Hedges’ *g*（获益 − 未获益）。置换 *p* 在队列内打乱获益标签（10,000 次）。

| 基因 | k | n（获益 / 未） | *g* | 95% CI | *p* | *I*² |
|---|---:|---|---:|---|---:|---:|
| TACSTD2 | 3 | 32（11 / 21） | −0.29 | −1.05 ~ 0.47 | 0.46 | 0% |
| CLDN4 | 3 | 32（11 / 21） | +0.03 | −0.97 ~ 1.03 | 0.96 | 39% |

同 32 例管家基因：ACTB *g* = −0.002（*p* = 1.00），GAPDH *g* = −0.18（*p* = 0.65），PPIA *g* = +0.03（*p* = 0.94）。

分队列 TACSTD2：France-3 *g* = −0.99（MWU *p* = 0.18，置换 *p* = 0.11）；新辅助 *g* = +0.01；PLA *g* = +0.28。无一队列单独显著。留一队列法也不能得到显著合并效应。

加入高置信推断 LUSC 后结论不变：TACSTD2 *g* = −0.20（k = 4，n = 39，*p* = 0.58）；CLDN4 *g* = −0.02（k = 4，n = 39，*p* = 0.95）。

**病理非 LUSC 对照**（非主结论）：TACSTD2 *g* = +0.17（k = 3，n = 38，*p* = 0.66）；CLDN4 *g* = +0.46（k = 4，n = 86，*p* = 0.17）。若把 LUSC 与 LUAD 混在一起，会把两个不同对象当成一个。

### PFS（Cox，分析子集内每 +1 SD）

病理 LUSC，France-3（n = 11，7 事件）：TACSTD2 HR = 4.08（95% CI 0.82–20.3，*p* = 0.086）；CLDN4 HR = 1.74（0.69–4.40，*p* = 0.24）。  
France-4 LUSC CLDN4（n = 11，9 事件）：HR = 1.03（0.51–2.08，*p* = 0.93）。  
队列内 z 分数合并 CLDN4（n = 22，16 事件）：HR = 1.24（0.77–2.01，*p* = 0.37）。

区间很宽，与无效或与高表达有害都相容，**不能**解读为获益证据。

### 免疫背景

签名分为队列内基因 z 分数均值（基因列表与来源见 `signatures.json`）。Spearman ρ 经 Fisher *z* 随机效应合并。

**病理 LUSC ICI（主分析集）：**

| 基因 | 签名 | k | n | ρ | 95% CI | *p* | BH *q* |
|---|---|---:|---:|---:|---|---:|---:|
| TACSTD2 | EPITHELIAL | 4 | 42 | +0.72 | 0.50 ~ 0.85 | 6.8×10⁻⁷ | 1.8×10⁻⁵ |
| CLDN4 | EPITHELIAL | 5 | 53 | +0.69 | 0.48 ~ 0.82 | 2.1×10⁻⁷ | 1.1×10⁻⁵ |
| TACSTD2 | CHECKPOINT | 4 | 42 | −0.41 | −0.66 ~ −0.07 | 0.018 | 0.11 |
| TACSTD2 | CD8_T_CELL | 4 | 42 | −0.32 | −0.60 ~ 0.03 | 0.068 | 0.20 |
| TACSTD2 | CYT | 4 | 42 | −0.32 | −0.59 ~ 0.03 | 0.073 | 0.21 |
| TACSTD2 | IFNG_6 | 4 | 42 | +0.04 | −0.31 ~ 0.38 | 0.82 | 0.91 |

用上皮评分做偏 Spearman 后，ICI LUSC 的免疫相关进一步减弱，本文不声称其独立于肿瘤细胞含量。

**TCGA-LUSC（n = 501，非 ICI）— 只重复免疫轴：**

| 基因 | 签名 | ρ | *p* | BH *q* | 偏 ρ（上皮） |
|---|---|---:|---:|---:|---:|
| TACSTD2 | CD8_T_CELL | −0.285 | 7.8×10⁻¹¹ | 7.0×10⁻¹⁰ | −0.231 |
| TACSTD2 | TIS_18 | −0.247 | 2.0×10⁻⁸ | 7.0×10⁻⁸ | −0.202 |
| TACSTD2 | IFNG_6 | −0.247 | 2.1×10⁻⁸ | 7.0×10⁻⁸ | −0.206 |
| TACSTD2 | CYT | −0.205 | 3.6×10⁻⁶ | 9.2×10⁻⁶ | −0.153 |
| TACSTD2 | EPITHELIAL | +0.246 | 2.3×10⁻⁸ | 7.0×10⁻⁸ | — |
| CLDN4 | EPITHELIAL | +0.515 | 2.5×10⁻³⁵ | 4.5×10⁻³⁴ | — |
| CLDN4 | CD8_T_CELL | −0.115 | 0.010 | 0.020 | +0.040（校正后不显著） |

结论：在 LUSC 中两基因都贴**上皮 / 肿瘤细胞**程序。TACSTD2 在 501 例 TCGA-LUSC 中与 CD8 / IFN-γ / 细胞毒评分呈**弱负相关**，部分经得起上皮校正。ICI LUSC 小样本方向相同，但 **BH 后不显著**。CLDN4 的弱免疫负相关在上皮校正后**消失**。

### 核实（实际做了的）

1. 病理 ICI 样本上的组织学外推检查（见上表）。主分析只用病理标签。  
2. 队列内置换获益标签（10,000 次）。  
3. 管家基因作为疗效负对照 — 全部无效。  
4. 留一队列随机效应 — 仍无效。  
5. 加入高置信推断 LUSC 的敏感性分析 — 仍无效。  
6. 病理非 LUSC 对照单独报告，不与 LUSC 主结论合并。  
7. 对上皮评分的偏 Spearman。  
8. 仅在 TCGA-LUSC 重复免疫轴。  
9. 基因覆盖：France-4 无 TACSTD2，该队列从所有 TACSTD2 检验中删除，不填补。  
10. 处理后体积远小于 2 GB（见 `processed_size.json`）。

### 解释（不超出数据）

公开的“LUSC + ICI + RNA”样本**很少**。病理确诊队列里只有 11 例获益，本切片**不支持**把 TACSTD2 或 CLDN4 当作 LUSC 的 ICI 疗效标志物。点估计接近零，置信区间很宽。

免疫结果更稳：两基因都与上皮程序同向；TACSTD2 在 501 例 LUSC 参考中与 T 细胞 / IFN-γ 程序弱负相关。这是肿瘤内在 / 纯度邻近的关联，**不是** ICI 结局结论。

TROP2、claudin-4 作为 ADC 靶点可以在 LUSC 中高表达，但在本切片汇集的公开记录里，这种表达并不标记 ICI 获益。

### 限制

- 样本量小，部分队列仅 2 例获益。  
- 终点不同（6 个月 PFS / MPR / CR）。  
- 平台不同；统计均在队列内进行。  
- France-4 为缩减基因面板。  
- SMC / 延世组织学为推断，仅敏感性。  
- 无统一 PD-L1 / TMB。  
- Bulk RNA 不能把 TACSTD2/CLDN4 定位到肿瘤细胞。

### 文件

- 脚本：`scripts/opus_lusc_ici/`  
- 表格：`results/opus_lusc_ici/tables/`  
- 图：`results/opus_lusc_ici/figures/`  
- 下载清单（URL、字节、SHA-256）：`results/opus_lusc_ici/manifest/raw_downloads.json`
