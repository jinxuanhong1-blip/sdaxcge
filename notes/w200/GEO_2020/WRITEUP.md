# GEO 2020 leftover: TACSTD2 / CLDN4 vs lung ICI outcomes

# 2020 年 GEO 遗漏切片：肺癌 ICI 结局与 TACSTD2 / CLDN4

> Wave `w200`. All outputs live under `notes/w200/GEO_2020/`,
> `scripts/w200/GEO_2020/`, `results/w200/GEO_2020/`.
> Every GSE accession is a real NCBI identifier. None were invented.
> 本切片全部产物仅位于上述三个目录。所有 GEO 编号均为 NCBI 真实编号。

---

## English

### Honest answer

**There is no leftover 2020 GEO series that can test tumor TACSTD2 or CLDN4
against a per-patient lung ICI outcome.**

A calendar-2020-only NCBI `gds` search (GSE, *Homo sapiens*, lung × ICI terms,
PDAT 2020/01/01–2020/12/31, plus per-drug queries) returned **28 real series**.
That set is identical to the 2020-dated subset of the earlier 2019–2021 wave.
Two of the 28 were already finished:

| GSE | Why it is not a leftover |
|---|---|
| **GSE126044** | Tumor RNA-seq, anti-PD-1, 5 responders / 11 non-responders, both genes present. Already analyzed (TACSTD2 AUC 0.36, p = 0.44; CLDN4 AUC 0.24, p = 0.11; both higher in non-responders; neither significant). |
| **GSE136961** | NSCLC anti-PD-1 DCB/NDB, but the Oncomine Immune Response 395-gene panel does not measure TACSTD2 or CLDN4. Already documented. |

The 2019–2021 wave had parked five 2020 series as `LUNG_ICI_other`. This slice
opened those records (series-matrix headers, characteristics, supplementary
listings, platform maps) instead of trusting the keyword triage.

| GSE | What it actually is | Can it answer the question? |
|---|---|---|
| **GSE141479** | 74 Clariom D arrays, **sorted CD8+ PBMC**, NSCLC, nivolumab pre/post (41 patients). Both genes are on GPL23126. GEO has **no** responder / RECIST / PFS / OS field. | **No.** Wrong compartment; no clinical endpoint in GEO. |
| **GSE154286** | 133 NSCLC biopsies, SAKK19/09, **chemotherapy** diagnosis vs progression. Custom NanoString **201-gene** panel. TACSTD2 and CLDN4 are absent (panel IDs committed). | **No.** Not ICI; genes not measured. |
| **GSE150972** | One EGFR-mutant patient who benefited from anti-PD-1. Empty expression table; RAW only. | **No.** n = 1, no processed matrix. |
| **GSE99995** | 12 LUAD tumors stratified by IFN-γ / PD-L1. No ICI treatment. | **No.** Not an ICI leftover. Not used as a surrogate. |
| **GSE124885** | Pediatric autoimmune lung / post-infection CD4 study, n = 4. | **No.** Not lung-cancer ICI. |

Single-cell, methylation companion arms of GSE126044, cell-line, and
non-lung hits are listed in `results/w200/GEO_2020/tables/triage_2020.csv`.
They were not re-analyzed.

### What was measured on the one leftover that has the genes

GSE141479 (Hatae / Chamoto et al., *JCI Insight* 2020, PMID 31855576).
Probes: TACSTD2 `TC0100014340.hg.1`, CLDN4 `TC0700007993.hg.1`.
Values are linear Clariom D normalized intensities.

| Gene | n (pre / post) | median pre | median post | Mann–Whitney p | paired Wilcoxon p (n=33) | probe-mean percentile among 135,750 probes |
|---|---|---|---|---|---|---|
| TACSTD2 | 41 / 33 | 15.97 | 16.92 | 0.40 | 0.49 | 23.8th (array median-of-means = 31.2) |
| CLDN4 | 41 / 33 | 10.03 | 10.13 | 0.62 | 0.76 | 3.7th |

Both probes sit in the **low tail** of the array, which is what an epithelial
gene should do in sorted CD8+ T cells. The pre vs post contrast is **not** an
ICI outcome test. GEO does not deposit response labels; this slice does not
invent them from the paper.

### What this does not change

The only 2020 GEO tumor cohort that can speak to TACSTD2 / CLDN4 vs ICI
response remains **GSE126044**, already reported, underpowered (n = 16), and
non-significant. This leftover pass did not find a second 2020 cohort. It
closed the parked `LUNG_ICI_other` list instead of leaving it as a maybe.

### Reproducibility

Scripts in `scripts/w200/GEO_2020/`. Master tables:
`tables/triage_2020.csv`, `tables/leftover_probe.csv`,
`tables/leftover_verdict.csv`,
`tables/combined_TACSTD2_CLDN4_leftover_2020.csv`.
Full audit: `notes/w200/GEO_2020/verification_log.md`.

---

## 中文

### 诚实结论

**2020 年 GEO 中没有遗漏的、可用肿瘤 TACSTD2 / CLDN4 检验肺癌 ICI 结局的系列。**

仅限日历年 2020 的 NCBI `gds` 检索（GSE、人类、肺癌 × ICI、PDAT
2020/01/01–2020/12/31，并加药名检索）得到 **28 个真实系列**，与 2019–2021
波次中日期为 2020 的子集一致。其中 2 个早已做完：

- **GSE126044**：肿瘤 RNA-seq，抗 PD-1，5 应答 / 11 非应答，两基因均在。已分析（TACSTD2 AUC 0.36，p = 0.44；CLDN4 AUC 0.24，p = 0.11；均在非应答者偏高，均不显著）。
- **GSE136961**：NSCLC 抗 PD-1 DCB/NDB，但 Oncomine 395 免疫 panel **不含** TACSTD2 / CLDN4。

2019–2021 波次把 5 个 2020 系列标成 `LUNG_ICI_other`。本切片打开了这些记录，而不是停留在关键词分诊。

- **GSE141479**：74 张 Clariom D，**分选的外周血 CD8+**，纳武利尤单抗治疗前后。两基因在平台上。GEO **没有**应答 / RECIST / PFS / OS。不能回答问题。
- **GSE154286**：SAKK19/09 **化疗**前后活检，定制 NanoString 201 基因，两基因不在 panel 上。
- **GSE150972**：1 例 EGFR 突变、抗 PD-1 获益；无处理后表达矩阵。
- **GSE99995**：IFN-γ / PD-L1 分层的肺腺癌，不是 ICI 治疗队列。不作替代终点。
- **GSE124885**：儿童自身免疫肺 / 感染后，不是肺癌 ICI。

### 唯一测到两基因的遗漏系列

GSE141479。TACSTD2 / CLDN4 探针强度落在整张芯片的低尾（分别为第 23.8 与 3.7 百分位），符合上皮基因在 CD8+ T 细胞中的预期。治疗前后比较 **不是** ICI 临床结局检验。GEO 未提供应答标签，本切片不从论文里编造样本级标签。

### 这不改变什么

2020 年 GEO 里唯一能谈肿瘤 TACSTD2 / CLDN4 与 ICI 应答的队列仍是 **GSE126044**，已经报告过，n = 16，不显著。这次遗漏排查没有找到第二个 2020 队列，只是把此前挂起的 `LUNG_ICI_other` 名单关掉。
