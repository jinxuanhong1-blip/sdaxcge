# Leftover GEO 2019–2021 lung ICI series (no-skip pass)
# 2019–2021 GEO 肺癌 ICI 遗留系列（不因体积/组织类型跳过）

> Outputs live under `results/noskip/GEO_2019_2021/`.
> Every accession is a real NCBI GSE. **No fabricated statistics.**
> 全部编号均为真实 GEO 编号。**无任何编造统计量。**

---

## English

### What this pass did
The first pass analyzed 4 series and skipped leftovers for tissue (blood/PBMC/T cell)
or practical size/format (large scRNA, RAW.tar-only). This pass **revisited all 53
leftover lung+ICI candidates**, inventoried supplementary files, downloaded
processed matrices **including blood and files up to 958 MB**, scanned for
TACSTD2/CLDN4, and computed **only** what the deposited files support.

### Leftover ICI-outcome compute (the one that actually exists)

| GSE | n | genes | deposited endpoint | result |
|---|---|---|---|---|
| **GSE190266** | 70 lung-tumor biopsies, chemo ± anti-PD-1/PD-L1 | **CLDN4 present; TACSTD2 absent** from the 16,383-gene TPM matrix | 6-month truncated PFS (`pfs_time`, `pfs_evt`) | Cox HR = **0.917** per log2(TPM+1), p = **0.1246**; KM median-split log-rank p = **0.2105**. 52 events / 18 censored. **Not significant.** |

This is the only leftover 2019–2021 series that has (i) a processed matrix containing
at least one target gene and (ii) a deposited patient-level ICI-related survival
endpoint. Direction (HR < 1) is opposite the first-pass non-significant trend of
“higher CLDN4 → worse”; it is **not** a confirmation.

### Other leftover computes (honest labels — not ICI endpoints)

| GSE | what was computed | result | why it is not an ICI-outcome test |
|---|---|---|---|
| GSE190265 | TACSTD2/CLDN4 vs squamous vs non-squamous (n=25 with histology) | TACSTD2 p=0.31; CLDN4 p=0.12 | No response/PFS in GEO |
| GSE181820 | both genes vs GEO groups A/B/C (n=22) | Kruskal p=0.0305 / 0.0369; group C lower | Group meaning **not deposited**; not ICI response |
| GSE146100 | scRNA per-nodule means (1 patient, 3 nodules, pembrolizumab) | W3 (NR) highest; W2 (R) mid | n=3, **no inferential test** |
| GSE145896 | PD-1+ TIL CD28+ vs CD28− (n=8) | TACSTD2 p=0.69; CLDN4 p=1.00 | Subset RNA, not patient response |
| GSE131933 | Tumor-cell TACSTD2/CLDN4 at explant D0 (T1/T3) | descriptive means only | Ex-vivo durvalumab; T1D7 (958 MB) downloaded, **genes absent** in that file |
| GSE184053 | TACSTD2 in Temra | mean 0.41, mostly 0 | No ICI endpoint; T-cell compartment |

### Blood / large leftovers revisited (do-not-skip)

| GSE | action | honest finding |
|---|---|---|
| GSE152590 | downloaded blood CD8 TPM xlsx | TACSTD2 and CLDN4 **absent** (24530 genes) |
| GSE141479 | downloaded 57 MB blood CD8 series_matrix | probe IDs `TC*.hg.1`; **no gene symbols**; no response/PFS |
| GSE179994 | downloaded 441 MB T-cell RDS | TACSTD2 string present; **no patient-level stats extracted**; GEO has pre/on-treatment only, no response/PFS |
| GSE180347 | series_matrix only (144 LUAD) | PD-L1 group present; **RAW.tar only — no processed matrix** |
| GSE176021 | series_matrix | **response status is deposited** but suppl is annotation RDS + RAW.tar — **no processed gene matrix**, so **not computed** |
| GSE173351 | series_matrix | RAW.tar only |

### Interpretation (leftover + first pass, still honest)
Across 2019–2021 GEO human lung ICI series that actually deposit processed
TACSTD2/CLDN4 **and** a patient ICI endpoint, the evidence remains **negative /
underpowered**:

- First pass: GSE126044, GSE135222 — non-significant; GSE182328 CLDN4 vs Akkermansia
  surrogate nominal p=0.045 (not a deposited clinical endpoint).
- This pass: GSE190266 CLDN4 vs 6-month PFS — **not significant** (p=0.12 / 0.21).

No leftover series changes that conclusion. Several attractive leftovers
(GSE176021 response labels; GSE180347 n=144) **cannot be computed** because GEO
did not deposit a processed expression matrix.

### Reproducibility
`07_noskip_inventory.py` → `08_noskip_download.py` → `09_noskip_detect.py` →
`10_noskip_analyze.py` → `11_noskip_catalog.py`.
Large raw downloads are not stored in git; they are re-fetched by script 08.

---

## 中文

### 本轮做了什么
首轮分析了 4 个系列，其余因组织类型（血液/T 细胞）或体积/格式被搁置。本轮
**复查全部 53 个遗留肺癌+ICI 候选**，清点补充文件，下载已处理矩阵
（**包括血液及最大 958 MB 的文件**），扫描 TACSTD2/CLDN4，并**只计算沉积文件
能够支持的统计**。

### 真正具有 ICI 结局的遗留计算

| GSE | n | 基因 | 沉积终点 | 结果 |
|---|---|---|---|---|
| **GSE190266** | 70 例肺肿瘤活检，化疗 ± 抗 PD-1/PD-L1 | **有 CLDN4；TPM 矩阵无 TACSTD2** | 截断于 6 个月的 PFS | Cox HR=**0.917**，p=**0.1246**；KM log-rank p=**0.2105**。**不显著。** |

这是遗留系列中唯一同时具备「已处理矩阵含目标基因」和「样本级 ICI 相关生存终点」的数据集。
方向（HR<1）与首轮「CLDN4 越高越差」的不显著趋势相反，**不能视为验证**。

### 其他计算（明确标注：不是 ICI 终点）
GSE190265 仅为组织学；GSE181820 的 A/B/C 分组含义未在 GEO 中给出（不可当作 ICI 应答）；
GSE146100 为 1 名患者 3 个结节的描述性均值，**未做推断检验**；GSE145896 / GSE131933 /
GSE184053 分别为 TIL 亚群、体外药敏、T 细胞隔室。

### 血液/大文件复查
GSE152590、GSE141479：已下载，目标基因符号不在沉积矩阵中。
GSE179994：441 MB RDS 已下载，未提取、未编造患者级统计。
GSE176021：GEO **有应答标签**，但**没有已处理表达矩阵**，故未计算。
GSE180347 / GSE173351：仅有 RAW.tar。

### 结论
把遗留系列算进来后，2019–2021 GEO 中「有已处理 TACSTD2/CLDN4 + 患者 ICI 终点」的证据
仍然是**阴性/效能不足**。GSE190266 的 CLDN4–PFS 不显著。没有编造任何缺失统计量。
