# GEO 2021 leftover lung ICI × TACSTD2 / CLDN4

> Outputs live only under `results/w200/GEO_2021/`. Every GSE accession is a real NCBI identifier.

---

## English

### 1. What “leftover” means
A 2021-only NCBI `gds` search (human GSE, lung × ICI terms, PDAT 2021-01-01 to 2021-12-31, plus per-drug queries) returned **38 real series** (`candidates_2021.json`). Two were already analyzed in the 2019–2021 slice (GSE111414, GSE182328). The leftover question is whether any remaining 2021 series can still test **tumor TACSTD2 / CLDN4 versus ICI outcome**.

### 2. Leftover triage
Full table: `leftover_triage.csv`.

| decision | n | meaning |
|---|---|---|
| already_analyzed | 2 | GSE111414, GSE182328 |
| analyze | 2 | GSE190265, GSE190266 |
| side_audit | 1 | GSE146100 (one patient, three nodules) |
| exclude | 33 | wrong disease, wrong assay, T-cell/TIL only, cell line, no ICI endpoint, or missing target genes |

No additional leftover bulk-tumor + ICI-outcome cohort appeared. GSE180347 is lung but Nanostring-without-targets and has no ICI endpoint. GSE181820 has both genes but no deposited response/survival. T-cell series (GSE173351/176021/176022/179994/144945/184053) cannot measure tumor epithelial TACSTD2/CLDN4.

### 3. The only leftover analyzable cohorts
GSE190270 subseries (Limagne / Ghiringhelli; PMID 35051357). Advanced NSCLC, anti-PD-1. Tumor TPM from GEO.

| cohort | n | genes present | outcome |
|---|---|---|---|
| GSE190265 France3 | 43 aligned suppl. records (34 GSM) | TACSTD2, CLDN4, CXCL10, CD274, OPTN, TLR9 | uncapped PFS |
| GSE190266 France4 | 70 | CLDN4, CXCL10, CD274; **not TACSTD2** | PFS already capped at 6 months |

Honest data limits:
- France4’s matrix has **16,383 genes and ends at MTMR14**. That is an Excel-style 2^14 column cap. TACSTD2 sorts after T and is absent. Inventing a France4 TACSTD2 result would be fabrication.
- France4 Cox models use **capped time** and are biased; DCB (still progression-free at the 6-month cap) is the cleaner binary readout there.
- France3 histology is labeled for only 25/43 overlapping GSM records.
- No RECIST field is deposited. DCB = PFS ≥ 6 months; NDB = event before 6 months.

### 4. Primary leftover result: TACSTD2 / CLDN4
`survival_results.csv`. Expression = log2(TPM+1). Univariate unless noted.

| cohort | gene | n DCB/NDB | AUC | Wilcoxon p (BH) | Cox HR (p) | log-rank p |
|---|---|---|---|---|---|---|
| GSE190265 | TACSTD2 | 14/29 | 0.46 | 0.67 (0.93) | 1.08 (0.39) | 0.83 |
| GSE190265 | CLDN4 | 14/29 | 0.49 | 0.93 (0.93) | 1.08 (0.58) | 0.94 |
| GSE190266 | TACSTD2 | — | — | **absent** | — | — |
| GSE190266 | CLDN4 | 17/52 | 0.66 | 0.052 (0.15) | 0.92 (0.12) | 0.21 |

Histology-adjusted Cox p-values remain >0.15. Spearman vs uncapped France3 PFS is null (TACSTD2 ρ=−0.03; CLDN4 ρ=0.05).

**Verdict:** leftover 2021 GEO does **not** support TACSTD2 or CLDN4 as ICI PFS/DCB markers. The single unadjusted CLDN4 Wilcoxon in France4 (p=0.052) fails BH correction, fails log-rank, and is absent in France3.

### 5. Control check (not the leftover claim)
The same files *can* detect ICI-related signal. Paper-claimed CXCL10 and deposited CD274 associate with DCB (`control_results.csv`): France3 CD274 Wilcoxon p=0.0028 (BH q=0.017); France4 CXCL10 p=0.0064 (q=0.019). That makes the TACSTD2/CLDN4 null more credible, not less. These controls were not the leftover hypothesis and are not a new biomarker claim.

### 6. GSE146100 side audit
One woman, three post-pembrolizumab nodules. Co-detection of TACSTD2+CLDN4: W1 2.17%, W2 4.68%, W3 14.81%. Not a consistent response association (nonresponding W1 is lowest). Cells are epithelial-like, not proven malignant. See `GSE146100_note.md`.

### 7. Reproduce
```bash
python3 search_2021.py
python3 analyze_leftovers.py
python3 analyze.py GSE146100_NormData.txt.gz
```

---

## 中文

### 1. “leftover” 指什么
对 NCBI `gds` 做 **2021 年**限定检索（人类 GSE、肺癌 × ICI、PDAT 2021-01-01 至 2021-12-31，并按药物补检索）得到 **38 个真实系列**。其中 GSE111414、GSE182328 已在 2019–2021 切片分析过。本 leftover 问题是：剩下的 2021 系列里，是否还有队列能检验**肿瘤 TACSTD2 / CLDN4 与 ICI 结局**。

### 2. 分诊
完整表：`leftover_triage.csv`。2 个已分析，2 个可分析（GSE190265/266），1 个单病例旁证（GSE146100），33 个排除（病种不对、T 细胞/TIL、细胞系、无 ICI 终点、或测不到目标基因）。没有新的“肿瘤 bulk + ICI 结局” leftover 队列。

### 3. 唯一可分析 leftover
GSE190270 子系列（PMID 35051357）。GSE190266 的 TPM 只有 16,383 个基因、止于 MTMR14，属 Excel 式 16384 列上限，**TACSTD2 不在文件中**；其 PFS 已被截断到 6 个月。GSE190265 补充文件有 43 条与表达对齐的记录，GEO 只挂 34 个 GSM。无 RECIST 字段。

### 4. 主结果
TACSTD2 在 France3 为无效（AUC 0.46，Wilcoxon p=0.67，Cox p=0.39）。CLDN4 在 France3 无效；在 France4 未校正 Wilcoxon p=0.052，BH q=0.15，log-rank p=0.21。**2021 leftover GEO 不支持这两个基因作为 ICI PFS/DCB 标志物。**

### 5. 对照
同一批文件能检出 ICI 相关信号：France3 的 CD274、France4 的 CXCL10 与 DCB 相关。这使 TACSTD2/CLDN4 的阴性更可信，而不是新的 leftover 结论。

### 6. GSE146100
一例三病灶，TACSTD2+/CLDN4+ 细胞比例并非按应答整齐分开。不能当作患者级证据。
