# GEO 2025–2026 lung IO scRNA: malignant TACSTD2 / CLDN4 vs response and T/NK

Bilingual notes follow the repo convention. Numbers will be filled from `combined_results.json` after the analysis run. Accessions are live GEO IDs only.

## 中文摘要

对 NCBI GEO 做了一次现场 E-utilities 检索（2025–2026，人肺 + IO + scRNA）。72 条 GSE 中 49 条已在既有候选表（PR #35）中列出。**14 条 2026 年新 GSE 里，没有“肺原发 + ICI 治疗 + 带缓解标签 + 开放处理后矩阵”的系列。**

本切片因此做了两件事：

1. **GSE291670**（2025-03-31，已出现在候选表但从未做过本终点）：6 例 NSCLC，新辅助安罗替尼 + 卡瑞利珠单抗，样本名即为 MPR / Non-MPR。分析恶性细胞 TACSTD2、CLDN4 与病理缓解、以及与 T/NK 比例的关系。
2. **GSE325414**（2026-06-10，**新** NSCLC scRNA）：脉冲电场消融 treat-and-resect，**不是 ICI**，无缓解标签。仅检验恶性 TACSTD2/CLDN4 与 T/NK（及 mTLS）。

## English summary

A live GEO search found **no new 2026 lung-primary ICI scRNA series with response labels** beyond accessions already listed. This slice analyzes:

- **GSE291670** (listed as a candidate, never analyzed for this endpoint): malignant TACSTD2/CLDN4 vs MPR and vs T/NK.
- **GSE325414** (new 2026 NSCLC scRNA, not ICI): malignant TACSTD2/CLDN4 vs T/NK.

See `00_selection.md` for the full accession-level exclusion table. Results: `GSE291670_results.json`, `GSE325414_results.json`, `combined_results.json`, `tacstd2_cldn4_summary.png`.
