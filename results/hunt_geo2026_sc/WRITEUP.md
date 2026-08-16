# GEO 2025–2026 lung IO scRNA: malignant TACSTD2 / CLDN4 vs response and T/NK

Accessions are live GEO IDs from `geo_esearch_raw.json` / sample records. No IDs were invented. Full include/exclude table: `00_selection.md`.

## 中文摘要

现场检索 NCBI GEO（人肺 + IO + scRNA，PDAT 2025–2026）：**72** 条 GSE，其中 **49** 条已在 PR #35 候选表中列出。**14** 条 2026 年新 GSE **没有一条**同时满足：肺原发、ICI 治疗、带缓解标签、开放处理后矩阵。

最接近的 2026 年 ICB scRNA 是 **GSE295600**（RACIN，纳武利尤单抗+伊匹木单抗+低剂量放疗），但 GEO 样本的 cancer type 为卵巢/前列腺/结肠/胰腺/乳腺，**无肺原发**。

本切片因此分析：

1. **GSE291670**（2025-03-31）：曾出现在候选表，但从未做过「恶性 TACSTD2/CLDN4 × MPR / T-NK」终点。6 例 NSCLC，新辅助安罗替尼 + 卡瑞利珠单抗；样本名即为 MPR-1/2/3 与 Non-MPR-1/2/3。
2. **GSE325414**（2026-06-10，**新**）：NSCLC BD-Rhapsody scRNA，脉冲电场消融 treat-and-resect，**不是 ICI**，无缓解标签。只用作者细胞类型标签检验恶性 TACSTD2/CLDN4 与 T/NK（及 mTLS）。

### GSE291670 结果（患者水平 n=3 vs 3，检验力不足）

| 终点 | MPR 均值 | non-MPR 均值 | Mann–Whitney p |
|------|----------|--------------|----------------|
| 恶性 TACSTD2 | 0.234 | 0.075 | **0.40** |
| 恶性 CLDN4 | 0.169 | 0.105 | **0.40** |

方向：MPR 更高，但 **n=3 vs 3，不能称为显著**。Spearman（n=6）：TACSTD2 vs T/NK 比例 ρ=−0.37（p=0.47）；CLDN4 ρ=−0.66（p=0.16）。细胞水平 p 值（TACSTD2 1.7×10⁻⁴⁵）是伪重复，**不是患者水平结论**。

QC 后 28,847 细胞；标记基因 Leiden 注释：恶性 14,750，T/NK 1,210。

### GSE325414 结果（无 ICI 标签）

156,467 细胞，25 供者，62 样本。作者标签：恶性 43,328，T/NK 44,290。表达用作者 `nCount_RNA` 做 log1p(CP10k)，避免 12 基因子集再 normalize 的膨胀。

| 相关 | n | ρ | p |
|------|---|---|---|
| 样本水平 TACSTD2 vs T/NK 比例 | 62 | −0.22 | 0.086 |
| 样本水平 CLDN4 vs T/NK 比例 | 62 | −0.05 | 0.71 |
| 供者水平 TACSTD2 vs T/NK 比例 | 25 | −0.07 | 0.73 |
| 供者水平 CLDN4 vs T/NK 比例 | 25 | −0.12 | 0.57 |
| 恶性 TACSTD2：mTLS High vs Low | 25 vs 16 | — | 0.14 |
| 恶性 CLDN4：mTLS High vs Low | 25 vs 16 | — | 0.067 |

**结论：** 新 2026 年公开 GEO 里没有可用的「肺 ICI scRNA + 缓解标签」增量队列。GSE291670 患者水平为无效/检验力不足（方向是 MPR 更高，与「TROP2 高 → 冷肿瘤/耐药」相反）。GSE325414 恶性 TACSTD2 与 T/NK 呈弱负相关趋势（样本 p=0.086），供者水平消失。

## English summary

Live GEO search: **72** lung+IO+scRNA GSE (2025–2026); **49** already listed; **14** new 2026 series. **None of the new 2026 series is a lung-primary ICI-treated scRNA cohort with response labels and open processed matrices.** Closest new 2026 ICB scRNA (GSE295600) has no lung-primary patients on the GEO sample records.

### GSE291670 (listed candidate, first analysis of this endpoint)

6 NSCLC tumors, neoadjuvant anlotinib + camrelizumab. Sample titles = MPR vs Non-MPR (verbatim GEO). After QC: 28,847 cells; 14,750 malignant; 1,210 T/NK (marker-score Leiden).

Patient-level means (log-norm): TACSTD2 0.234 (MPR) vs 0.075 (non-MPR); CLDN4 0.169 vs 0.105. **Mann–Whitney p=0.40 for both (n=3 vs 3).** Direction is MPR-higher, opposite a “TROP2-high = immune-cold / ICI-resistant” claim, and underpowered. Spearman vs T/NK fraction (n=6): TACSTD2 ρ=−0.37 (p=0.47); CLDN4 ρ=−0.66 (p=0.16). Cell-level p-values treat cells as independent and are **not** a patient-level result.

### GSE325414 (new 2026 NSCLC scRNA; not ICI)

PEF ablation treat-and-resect. Author labels: 43,328 malignant, 44,290 T/NK of 156,467 cells / 25 donors. Expression: log1p(CP10k) using author `nCount_RNA`.

Sample-level TACSTD2 vs T/NK fraction: ρ=−0.22, p=0.086 (n=62). Donor-level ρ=−0.07, p=0.73 (n=25). CLDN4 vs T/NK null. mTLS High vs Low: TACSTD2 p=0.14, CLDN4 p=0.067.

### Honest takeaway

No new 2026 public GEO lung-ICI scRNA response cohort. The only unused labeled lung ICI scRNA with malignant cells and open 10x matrices in this hunt (GSE291670) is a 3-vs-3 null. The new 2026 NSCLC atlas (GSE325414) does not carry ICI labels; TACSTD2 vs T/NK is a weak, non-significant negative trend that does not survive donor aggregation.

## Files

- `geo_esearch_raw.json` / `geo_esearch_2026_broad.json` — live E-utilities dumps
- `candidates_new.tsv` / `broad_2026_lung_scrna.tsv` / `already_listed.txt`
- `GSE291670_per_sample_summary.csv`, `GSE291670_lineage_counts.csv`, `GSE291670_results.json`
- `GSE325414_per_sample_summary.csv`, `GSE325414_per_donor_summary.csv`, `GSE325414_results.json`
- `combined_results.json`, `tacstd2_cldn4_summary.png`
- scripts: `scripts/hunt_geo2026_sc/`
