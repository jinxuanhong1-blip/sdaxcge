# TACSTD2 (Trop-2) & CLDN4 (Claudin-4) vs. ICI outcomes in GEO 2022–2023 human lung series

*Parallel slice `fable_geo_2022_2023`. All outputs live under
`notes/fable_geo_2022_2023/`, `scripts/fable_geo_2022_2023/`, and
`results/fable_geo_2022_2023/`.*

---

## English

### 1. Objective
Perform an **exhaustive** search of GEO for **human lung-cancer immune
checkpoint inhibitor (ICI) series published in 2022–2023**, verify them,
download the open *processed* data that is < 2 GB, and test whether tumour
expression of **TACSTD2 (Trop-2)** and **CLDN4 (Claudin-4)** — two epithelial
antigens that are also antibody–drug-conjugate targets — is associated with
**ICI treatment outcome**.

### 2. Pipeline (fully reproducible)
| Stage | Script | What it does |
|---|---|---|
| 01 | `scripts/fable_geo_2022_2023/01_search_geo.py` | Exhaustive NCBI E-utilities search of `db=gds` (GSE series only, `Homo sapiens`, `2022/01/01–2023/12/31` `[PDAT]`), union of a broad lung×ICI query plus per-drug queries. |
| 02 | `02_verify_candidates.py` | Text relevance flags (human/lung/ICI/expression/outcome) + supplementary-file listing & byte sizes from the NCBI FTP HTTPS mirror; flags open, processed download < 2 GB. |
| 03 | `03_fetch_series_metadata.py` | Downloads every `*_series_matrix.txt.gz` header, parses per-sample `!Sample_characteristics_ch*`, detects outcome fields and whether expression is embedded (array) vs. metadata-only (RNA-seq). |
| 04 | `04_download_and_parse.py` | Downloads the processed expression matrices and parses the full per-sample metadata table. |
| 05 | `05_analyze_markers.py` | Extracts TACSTD2/CLDN4 per sample, runs Mann–Whitney U (responder vs. non-responder), Cliff's delta, Spearman correlation with residual tumour, and writes boxplots/scatterplots. |
| 06 | `06_build_catalog.py` | Assembles the master verified catalog and `catalog_summary.md`. |

### 3. Search & verification result (exhaustive catalog)
- **94** unique GSE series returned by the exhaustive query.
- **42** pass text relevance (human + lung + ICI + expression profiling).
- **30** carry per-sample outcome annotation in GEO metadata.
- **41** have an open, processed supplementary download < 2 GB.

Full tables:
`results/fable_geo_2022_2023/geo_lung_ici_catalog.csv` (master),
`geo_search_candidates.csv`, `geo_verified.csv`, `geo_suppl_files.csv`,
`series_characteristics_summary.csv`; per-series metadata JSON under
`results/fable_geo_2022_2023/series_matrix_meta/`.

**Critical verification finding — why most series cannot answer the question.**
The two largest, cleanest NSCLC ICI *survival* cohorts, **GSE162520** (n=92,
PD-1/PD-L1, OS/PFS) and **GSE161537** (n=82, immunotherapy, RECIST/OS/PFS), both
distribute expression as a **~2,560-gene targeted immune panel that does *not*
contain TACSTD2 or CLDN4** (it carries CLDN3, EPCAM, CD274, PDCD1, etc., but not
our two genes). Most other high-`n` series are **single-cell / sorted-immune**
studies (GSE160903, GSE185206, GSE235500/GSE235603), **cell-line** killing
assays (GSE214992), **xenografts** (GSE190731), or **platelet** RNA
(GSE216297) — none of which provide bulk tumour epithelial expression paired
with a per-patient ICI outcome. After verification, exactly **one lung cohort**
provides open, whole-transcriptome bulk tumour expression *and* a per-patient
ICI outcome and contains both genes.

### 4. Cohorts analysed
| Cohort | Tissue | n (analysed) | Expression | ICI regimen | Outcome |
|---|---|---|---|---|---|
| **GSE207422** (primary, lung) | NSCLC baseline (pre-treatment) tumour | 24 | bulk RNA-seq, log2 TPM | neoadjuvant anti-PD-1 (toripalimab/sintilimab/camrelizumab) + platinum chemo | Major pathologic response (MPR incl. pCR, n=9) vs. non-MPR (NMPR, n=15); residual-tumour fraction; RECIST |
| GSE243238 (cross-check, **non-lung**) | acral melanoma | 14 | bulk RNA-seq, raw counts → log2 CPM | ICI | clinical benefit yes (n=3) / no (n=11) |

### 5. Results
Association tests are in `results/fable_geo_2022_2023/analysis/marker_outcome_tests.csv`; figures in the same folder.

**GSE207422 (NSCLC, baseline tumour):**
| Gene | Outcome | Responder median | Non-responder median | p | effect |
|---|---|---|---|---|---|
| TACSTD2 | MPR vs NMPR | 5.00 | 6.22 (log2 TPM) | 0.34 | Cliff's δ −0.24 |
| CLDN4 | MPR vs NMPR | 5.85 | 6.27 | 0.26 | Cliff's δ −0.29 |
| TACSTD2 | Spearman vs residual-tumour % | — | — | 0.23 | ρ = +0.25 |
| CLDN4 | Spearman vs residual-tumour % | — | — | 0.31 | ρ = +0.22 |

**GSE243238 (acral melanoma, cross-check):**
| Gene | Outcome | Benefit median | No-benefit median | p | effect |
|---|---|---|---|---|---|
| TACSTD2 | clinical benefit | 0.40 | 3.55 (log2 CPM) | 0.29 | Cliff's δ −0.46 |
| CLDN4 | clinical benefit | 0.04 | 0.64 | 0.085 | Cliff's δ −0.70 |

**Consistent direction, not statistically significant.** In every test the
*better* ICI outcome is associated with *lower* baseline TACSTD2 and CLDN4
(responders/benefit have lower expression; higher marker levels track with more
residual tumour). The direction is identical across the lung cohort and the
independent melanoma cohort, and it is biologically coherent — high Trop-2 /
Claudin-4 mark a more epithelial, potentially less-immunogenic tumour phenotype.
However **no test reaches p < 0.05**; the melanoma CLDN4 signal (p≈0.085) is the
strongest and is limited by only 3 benefit cases.

### 6. Interpretation & limitations
- The hypothesis (high TACSTD2/CLDN4 → worse ICI response) is **directionally
  supported but not confirmed**. Power is the binding constraint: the only
  eligible lung cohort has 24 baseline samples, and the melanoma cross-check has
  a 3-vs-11 imbalance.
- The most informative lung ICI *survival* cohorts of 2022–2023 (GSE162520,
  GSE161537) are **unusable for these two genes** because they used a targeted
  panel — a concrete data-availability gap worth flagging for future ADC-target
  biomarker work.
- GSE207422 is a **neoadjuvant chemo-immunotherapy** setting, so the response
  signal is a combination effect, not ICI monotherapy.
- All expression values are the authors' processed matrices (log2 TPM / log2
  CPM); no re-alignment was performed.

### 7. Conclusion
Across an exhaustive 2022–2023 GEO screen, only one open lung ICI cohort with
whole-transcriptome bulk tumour data (GSE207422) can address TACSTD2/CLDN4, and
one non-lung ICI cohort (GSE243238) serves as a cross-check. Both show the same
non-significant trend: **lower baseline Trop-2 and Claudin-4 expression
associates with better ICI outcome.** This is a hypothesis-generating signal
that warrants testing in larger whole-transcriptome ICI cohorts.

---

## 中文

### 1. 目标
在 GEO 中**穷尽式**检索 **2022–2023 年发表的人类肺癌免疫检查点抑制剂（ICI）系列**，
进行核实，下载体积 < 2 GB 的公开**已处理**数据，并检验肿瘤中
**TACSTD2（Trop-2）** 与 **CLDN4（Claudin-4）**（两种上皮抗原，同时也是抗体偶联药物 ADC 的靶点）
的表达是否与 **ICI 治疗结局**相关。

### 2. 流程（完全可复现）
| 阶段 | 脚本 | 功能 |
|---|---|---|
| 01 | `01_search_geo.py` | 通过 NCBI E-utilities 检索 `db=gds`（仅 GSE 系列、`Homo sapiens`、发表日期 2022/01/01–2023/12/31），取"肺癌 × ICI"广义查询与各药物查询的并集。 |
| 02 | `02_verify_candidates.py` | 文本相关性标记（人类/肺/ICI/表达/结局）+ 从 NCBI FTP HTTPS 镜像列出补充文件及字节大小；标记 < 2 GB 的公开已处理下载。 |
| 03 | `03_fetch_series_metadata.py` | 下载每个 `series_matrix` 头部，解析每样本 `!Sample_characteristics_ch*`，识别结局字段以及表达是否内嵌（芯片）或仅元数据（RNA-seq）。 |
| 04 | `04_download_and_parse.py` | 下载已处理表达矩阵并解析完整的每样本元数据表。 |
| 05 | `05_analyze_markers.py` | 提取每样本 TACSTD2/CLDN4，进行 Mann–Whitney U 检验（应答 vs 非应答）、Cliff's delta、与残余肿瘤的 Spearman 相关，并绘制箱线图/散点图。 |
| 06 | `06_build_catalog.py` | 汇总核实后的主目录与 `catalog_summary.md`。 |

### 3. 检索与核实结果（穷尽式目录）
- 穷尽查询返回 **94** 个唯一 GSE 系列。
- **42** 个通过文本相关性（人类 + 肺 + ICI + 表达谱）。
- **30** 个在 GEO 元数据中带有每样本结局注释。
- **41** 个具有 < 2 GB 的公开已处理补充下载。

完整表格见
`results/fable_geo_2022_2023/geo_lung_ici_catalog.csv`（主表）等。

**关键核实发现——为何多数系列无法回答该问题。**
两个规模最大、最规范的 NSCLC ICI *生存*队列 **GSE162520**（n=92，PD-1/PD-L1，OS/PFS）与
**GSE161537**（n=82，免疫治疗，RECIST/OS/PFS），其表达数据均为
**约 2,560 个基因的靶向免疫面板，其中并不包含 TACSTD2 或 CLDN4**（含 CLDN3、EPCAM、
CD274、PDCD1 等，但没有我们关注的两个基因）。其余高样本量系列多为
**单细胞/分选免疫细胞**研究（GSE160903、GSE185206、GSE235500/GSE235603）、
**细胞系**杀伤实验（GSE214992）、**异种移植**（GSE190731）或**血小板** RNA（GSE216297），
均无法提供与每患者 ICI 结局配对的肿瘤上皮批量表达。核实后，恰好**只有一个肺癌队列**
同时具备公开的全转录组批量肿瘤表达、每患者 ICI 结局且包含这两个基因。

### 4. 纳入分析的队列
| 队列 | 组织 | n（分析） | 表达 | ICI 方案 | 结局 |
|---|---|---|---|---|---|
| **GSE207422**（主要，肺） | NSCLC 基线（治疗前）肿瘤 | 24 | 批量 RNA-seq，log2 TPM | 新辅助抗 PD-1（特瑞普利/信迪利/卡瑞利珠）+ 铂类化疗 | 主要病理缓解（MPR 含 pCR，n=9）vs 非 MPR（NMPR，n=15）；残余肿瘤比例；RECIST |
| GSE243238（交叉验证，**非肺**） | 肢端黑色素瘤 | 14 | 批量 RNA-seq，原始计数 → log2 CPM | ICI | 临床获益 是（n=3）/ 否（n=11） |

### 5. 结果
关联检验见 `results/fable_geo_2022_2023/analysis/marker_outcome_tests.csv`，图见同目录。

**GSE207422（NSCLC，基线肿瘤）：**
| 基因 | 结局 | 应答中位 | 非应答中位 | p | 效应量 |
|---|---|---|---|---|---|
| TACSTD2 | MPR vs NMPR | 5.00 | 6.22（log2 TPM） | 0.34 | Cliff's δ −0.24 |
| CLDN4 | MPR vs NMPR | 5.85 | 6.27 | 0.26 | Cliff's δ −0.29 |
| TACSTD2 | 与残余肿瘤 % 的 Spearman | — | — | 0.23 | ρ = +0.25 |
| CLDN4 | 与残余肿瘤 % 的 Spearman | — | — | 0.31 | ρ = +0.22 |

**GSE243238（肢端黑色素瘤，交叉验证）：**
| 基因 | 结局 | 获益中位 | 无获益中位 | p | 效应量 |
|---|---|---|---|---|---|
| TACSTD2 | 临床获益 | 0.40 | 3.55（log2 CPM） | 0.29 | Cliff's δ −0.46 |
| CLDN4 | 临床获益 | 0.04 | 0.64 | 0.085 | Cliff's δ −0.70 |

**方向一致，但未达统计学显著。** 在所有检验中，**更好的** ICI 结局都与**更低的**基线
TACSTD2 与 CLDN4 相关（应答/获益者表达更低；标志物越高，残余肿瘤越多）。该方向在肺癌
队列与独立的黑色素瘤队列中完全一致，且符合生物学逻辑——高 Trop-2/Claudin-4 提示更偏
上皮、可能免疫原性更低的肿瘤表型。但**没有一项检验达到 p < 0.05**；黑色素瘤 CLDN4 信号
（p≈0.085）最强，但仅有 3 例获益病例，限制明显。

### 6. 解读与局限
- 假设（高 TACSTD2/CLDN4 → ICI 应答更差）**方向上获支持但未获证实**。样本量是核心制约：
  唯一合格的肺癌队列仅 24 例基线样本，黑色素瘤交叉验证为 3 vs 11 的不平衡。
- 2022–2023 年最有信息量的肺癌 ICI *生存*队列（GSE162520、GSE161537）**因采用靶向面板
  而无法用于这两个基因**——这是一个值得关注的数据可得性缺口，尤其对未来 ADC 靶点的
  生物标志物研究。
- GSE207422 为**新辅助化免联合**背景，应答信号是联合效应，而非 ICI 单药。
- 所有表达值均为作者提供的已处理矩阵（log2 TPM / log2 CPM），未重新比对。

### 7. 结论
在 2022–2023 年 GEO 的穷尽式筛查中，只有一个具备全转录组批量肿瘤数据的公开肺癌 ICI 队列
（GSE207422）可用于评估 TACSTD2/CLDN4，另有一个非肺 ICI 队列（GSE243238）作为交叉验证。
二者呈现相同的（非显著）趋势：**更低的基线 Trop-2 与 Claudin-4 表达与更好的 ICI 结局相关。**
这是一个可产生假设的信号，值得在更大规模的全转录组 ICI 队列中进一步检验。
