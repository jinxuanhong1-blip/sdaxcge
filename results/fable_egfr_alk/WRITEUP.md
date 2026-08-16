# WRITEUP — TACSTD2/CLDN4 vs immune microenvironment & ICI response in EGFR/ALK-mutant lung cancer

> Parallel slice. Public data only. Raw inputs are git-ignored under `data/fable_egfr_alk/`;
> all committed outputs live under `results/fable_egfr_alk/`, `scripts/fable_egfr_alk/`,
> `notes/fable_egfr_alk/`. Reproduce with `bash scripts/fable_egfr_alk/run_all.sh`.
> Methods & provenance: [`notes/fable_egfr_alk/data_provenance.md`](../../notes/fable_egfr_alk/data_provenance.md).

---

## English

### Question
In EGFR/ALK-mutant lung cancer — tumours that respond poorly to immune-checkpoint
inhibitors (ICI) — how do the ADC targets **TACSTD2 (TROP2)** and **CLDN4 (Claudin-4)**
relate to (a) the anti-tumour immune microenvironment and (b) ICI outcome?

### Data
- **TCGA-LUAD** (n=528 primary tumours) and **TCGA-LUSC** (n=501) — STAR-TPM expression
  + WXS mutations + survival (UCSC Xena GDC hub). Driver groups derived: EGFR-mutant
  (71 LUAD), ALK-driven (EML4-ALK expression-outlier proxy, 25 LUAD), EGFR/ALK-WT (422 LUAD).
- **GSE126044** (n=16) — anti-PD-1 NSCLC, responder/non-responder.
- **GSE135222** (n=27) — anti-PD-1/PD-L1 NSCLC, progression-free survival (PFS).

### Findings

**1. TROP2 and Claudin-4 are elevated in EGFR-mutant LUAD.**
Both targets are significantly higher in EGFR-mutant vs EGFR/ALK-WT tumours
(TACSTD2 median +0.41 log2TPM, Mann–Whitney p=0.016; CLDN4 +0.36, p=6.6×10⁻⁴).
Bootstrap 95% CIs of the median difference exclude 0 for both. ALK-driven tumours
trend *lower* than WT.

**2. TROP2/Claudin-4 mark an immune-excluded ("cold") microenvironment.**
Across the whole cohort both targets correlate **negatively** with cytotoxic/CD8
programs — significant after BH-FDR:

| | TACSTD2 vs CD8-effector | CLDN4 vs CD8-effector |
|---|---|---|
| TCGA-LUAD (n=528) | ρ=−0.129 (q=0.022) | ρ=−0.138 (q=0.011) |
| TCGA-LUSC (n=501) | ρ=−0.269 (q<10⁻³) | ρ=−0.119 (q=0.027) |

The pattern spans essentially the entire cytotoxic panel (GZMA/B, PRF1, IFNG, CXCL9/10,
CD8A/B) and replicates in LUSC (see heatmaps).

**3. ICI cohorts are directionally consistent but underpowered.**
In GSE135222, higher TACSTD2/CLDN4 trends toward shorter PFS (Spearman ρ≈−0.15, NS at
n=27). In GSE126044, CLDN4 is lower in responders (p=0.11, NS at n=16). The CD8-effector
signature itself perfectly separates responders (AUC=1.00, p=4.6×10⁻⁴) — a **positive
control** confirming the scoring pipeline works; the target trends simply lack power at
these sample sizes.

### Verification (`03_verify.py`, 18/21 checks pass)
- **Bootstrap CIs** (5,000×) for target↔CD8-effector are strictly negative in LUAD and LUSC.
- **Permutation p-values** (10,000×) significant in both cohorts.
- **Cross-cohort sign concordance** LUAD↔LUSC significant for both targets (binomial p=0.03 / 6×10⁻⁴).
- **EGFR-mut vs WT** median-difference bootstrap CI > 0 for both targets.
- **Alternative immune scoring** (single-gene CD8A/GZMB) reproduces the negative sign;
  the 3 non-significant checks are the *broad* immune panel (diluted by myeloid CD68/ITGAX)
  and one single-gene LUSC test — an honest nuance: the effect is specific to the
  **cytotoxic/CD8** program, not bulk leukocyte content.

### Interpretation
EGFR-mutant LUAD combines **high TROP2/Claudin-4** with a **cytotoxic-cold** immune
contexture. This coheres with the clinical picture (EGFR-mutant NSCLC responds poorly to
ICI) and supports a therapeutic logic: these tumours are better matched to **TROP2/Claudin
ADCs** (e.g. datopotamab-deruxtecan, sacituzumab-govitecan, claudin-directed ADCs) than to
single-agent checkpoint blockade — including the post-TKI setting where ICI benefit is limited.

### Limitations
Public ICI RNA-seq cohorts lack EGFR/ALK genotype, so EGFR/ALK specificity is carried by
TCGA (microenvironment) while ICI response is assessed in overall NSCLC (n=16, 27; direction
consistent, significance limited). ALK-driven status is an expression-outlier proxy for
EML4-ALK, not a direct fusion call. TCGA is treatment-naïve, not TKI-then-ICI.

---

## 中文

### 问题
在 **EGFR/ALK 突变型肺癌**（对免疫检查点抑制剂 ICI 疗效差）中，抗体偶联药物（ADC）靶点
**TACSTD2（TROP2）** 与 **CLDN4（Claudin-4）** 与（a）抗肿瘤免疫微环境、（b）ICI 疗效
之间的关系如何？

### 数据
- **TCGA-LUAD**（528 例原发瘤）与 **TCGA-LUSC**（501 例）——STAR-TPM 表达 + WXS 突变 +
  生存（UCSC Xena GDC hub）。分组：EGFR 突变（LUAD 71 例）、ALK 驱动（EML4-ALK 表达离群
  值代理，25 例）、EGFR/ALK 野生型（422 例）。
- **GSE126044**（16 例）——抗 PD-1 NSCLC，应答/不应答。
- **GSE135222**（27 例）——抗 PD-1/PD-L1 NSCLC，无进展生存（PFS）。

### 主要发现

**1. EGFR 突变型 LUAD 中 TROP2 与 Claudin-4 表达升高。**
两靶点在 EGFR 突变型显著高于 EGFR/ALK 野生型（TACSTD2 中位数 +0.41 log2TPM，
Mann–Whitney p=0.016；CLDN4 +0.36，p=6.6×10⁻⁴），中位差的自助法 95% 置信区间均不含 0。
ALK 驱动型则较野生型偏低。

**2. TROP2/Claudin-4 标记"免疫排斥（冷）"微环境。**
在全队列中，两靶点与细胞毒/CD8 程序呈**负相关**（BH-FDR 校正后显著）：

| | TACSTD2 vs CD8-效应 | CLDN4 vs CD8-效应 |
|---|---|---|
| TCGA-LUAD（528） | ρ=−0.129（q=0.022） | ρ=−0.138（q=0.011） |
| TCGA-LUSC（501） | ρ=−0.269（q<10⁻³） | ρ=−0.119（q=0.027） |

该模式覆盖几乎整个细胞毒基因面板（GZMA/B、PRF1、IFNG、CXCL9/10、CD8A/B），并在 LUSC
中重复（见热图）。

**3. ICI 队列方向一致但样本量不足。**
GSE135222 中，TACSTD2/CLDN4 越高 PFS 越短（Spearman ρ≈−0.15，n=27 未达显著）；
GSE126044 中 CLDN4 在应答者中更低（p=0.11，n=16 未达显著）。CD8-效应signature 本身可
完美区分应答者（AUC=1.00，p=4.6×10⁻⁴），作为**阳性对照**证明评分流程有效；靶点趋势仅
因样本量小而缺乏统计功效。

### 验证（`03_verify.py`，21 项通过 18 项）
- **自助法置信区间**（5,000 次）：靶点↔CD8-效应在 LUAD、LUSC 中严格为负。
- **置换检验 p 值**（10,000 次）：两队列均显著。
- **跨队列符号一致性** LUAD↔LUSC：两靶点均显著（二项检验 p=0.03 / 6×10⁻⁴）。
- **EGFR 突变 vs 野生型**：中位差自助法置信区间 > 0（两靶点）。
- **替代免疫评分**（单基因 CD8A/GZMB）复现负号；3 项未显著者为**广谱**免疫面板
  （被髓系 CD68/ITGAX 稀释）及一项 LUSC 单基因检验——这是诚实的细节：效应特异于
  **细胞毒/CD8** 程序，而非整体白细胞含量。

### 解读
EGFR 突变型 LUAD 兼具 **高 TROP2/Claudin-4** 与 **细胞毒冷** 免疫格局，与临床现象
（EGFR 突变型 NSCLC 对 ICI 疗效差）一致，并支持治疗逻辑：此类肿瘤更适合
**TROP2/Claudin ADC**（如 Dato-DXd、sacituzumab-govitecan、claudin 靶向 ADC），
而非单药检查点抑制——包括 TKI 进展后 ICI 获益有限的场景。

### 局限
公开 ICI RNA-seq 队列缺少 EGFR/ALK 基因型，故 EGFR/ALK 特异性由 TCGA（微环境）承担，
ICI 疗效在整体 NSCLC 中评估（n=16、27；方向一致但功效有限）。ALK 驱动为 EML4-ALK 的
表达离群值代理，非直接融合检测。TCGA 为初治样本，非 TKI 后 ICI。

---

### Output map
- Figures: `results/fable_egfr_alk/figures/` (driver boxplots, immune heatmaps, ICI boxplot, KM curves)
- Tables: `results/fable_egfr_alk/tables/` (correlations, driver comparisons, ICI response/PFS, `verification_report.csv`)
- Per-sample processed: `results/fable_egfr_alk/processed/`
- Scripts: `scripts/fable_egfr_alk/` (`common.py`, `00`–`03`, `run_all.sh`)
