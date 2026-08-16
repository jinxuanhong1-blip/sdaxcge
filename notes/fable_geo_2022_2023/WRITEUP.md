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
| 07 | `07_leftover_audit.py` | Deep metadata + gene-presence audit of every leftover / borderline series. |
| 08 | `08_leftover_analyze.py` | TACSTD2/CLDN4 tests on leftover series that have a real outcome join. |

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
with a per-patient ICI outcome. The first pass therefore analysed **GSE207422**
as the only whole-transcriptome lung tumour + ICI-outcome + both-genes cohort.

A leftover pass then finished every remaining relevant series and recovered
**GSE221733** (missed because the ICI regex required `immunotherap\b`, which
does not match the word *immunotherapy*). See
`notes/fable_geo_2022_2023/leftover_verdicts.md`.

### 4. Cohorts analysed
| Cohort | Tissue | n (analysed) | Expression | ICI regimen | Outcome |
|---|---|---|---|---|---|
| **GSE207422** (primary, lung) | NSCLC baseline (pre-treatment) tumour | 24 | bulk RNA-seq, log2 TPM | neoadjuvant anti-PD-1 (toripalimab/sintilimab/camrelizumab) + platinum chemo | Major pathologic response (MPR incl. pCR, n=9) vs. non-MPR (NMPR, n=15); residual-tumour fraction; RECIST |
| GSE243238 (cross-check, **non-lung**) | acral melanoma | 14 | bulk RNA-seq, raw counts → log2 CPM | ICI | clinical benefit yes (n=3) / no (n=11) |
| **GSE221733** (leftover, recovered) | NSCLC GeoMx DSP PanCK+ AOIs | 34 patients | NanoString CTA (TACSTD2 only; CLDN4 absent) | immunotherapy | Responder n=15 vs Non-responder n=19; OS (follow-up + vital status) |
| **GSE248378** (leftover) | NSCLC post-ICI resected tumour (non-MPR set) | 29 | bulk RNA-seq FPKM | neoadjuvant durvalumab +/- SBRT | recurrence 9 vs 20 (Nat Commun source-data join); DFS |
| GSE193049 (leftover, **BALF not tumour**) | bronchoalveolar lavage | 7 | bulk RNA-seq | PD-1 blockade | responder n=4 vs non-responder n=3 |

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

**Leftover tumour tests** (`analysis/leftover_marker_outcome_tests.csv`):

**GSE221733 (NSCLC DSP, TACSTD2 only, patient-level PanCK+):**
| Gene | Outcome | Responder median | Non-responder median | p | effect |
|---|---|---|---|---|---|
| TACSTD2 | ICI response | 9.268 | 9.378 (CTA norm.) | 0.5324 | Cliff's δ −0.13 |
| TACSTD2 | OS log-rank (median split) | high n=17 | low n=17 | 0.1545 | — |

**GSE248378 (NSCLC post-ICI non-MPR FPKM; paper n=29):**
| Gene | Outcome | No-recurrence median | Recurrence median | p | effect |
|---|---|---|---|---|---|
| TACSTD2 | recurrence | 41.98 | 80.01 (FPKM) | 0.0562 | Cliff's δ −0.456 |
| CLDN4 | recurrence | 38.97 | 77.37 | 0.3108 | Cliff's δ −0.244 |
| TACSTD2 | DFS log-rank (median split) | high n=15 | low n=14 | 0.1220 | — |

Sensitivity excluding arm-discordant `45-M-PO`: TACSTD2 recurrence p=0.0766, CLDN4 p=0.3759.

**GSE193049 (BALF, not tumour, n=7):** TACSTD2 p=0.6286 (δ=+0.333); CLDN4 p=0.4000 (δ=+0.50) — opposite direction to the tumour series; underpowered.

**Direction in tumour tissue, still not p < 0.05.** In every *tumour* test the
better ICI outcome is associated with *lower* TACSTD2 (and, where measured,
CLDN4). The leftover GSE248378 TACSTD2–recurrence test is the closest
(p=0.0562) and uses the paper’s own 9-vs-20 non-MPR contrast. No leftover
test was invented: GSE248378 has no MPR contrast in GEO (29/29 are non-MPR);
GSE221733 has no CLDN4; GSE193049 is BALF.

### 6. Interpretation & limitations
- The hypothesis (high TACSTD2/CLDN4 → worse ICI response) is **directionally
  supported in tumour tissue but not confirmed**. No tumour test reached
  p < 0.05. The strongest leftover signal is GSE248378 TACSTD2 vs recurrence
  (p=0.0562) in post-treatment non-MPR tumours — not a baseline predictor.
- GSE221733 was recovered only after leftover audit; CLDN4 is absent from the
  CTA panel, so that series cannot test Claudin-4.
- GSE162520 / GSE161537 remain unusable for these two genes (targeted panel).
- GSE207422 and GSE248378 are neoadjuvant (chemo-IO or durvalumab +/- SBRT),
  not ICI monotherapy in metastatic disease.
- All values are author-processed matrices; no re-alignment.

### 7. Conclusion
After finishing every leftover 2022–2023 GEO lung ICI series, three lung
datasets can be joined to an ICI outcome without inventing labels: GSE207422
(baseline bulk, both genes), GSE221733 (DSP PanCK+, TACSTD2 only), and
GSE248378 (post-Rx non-MPR bulk, both genes). Tumour tests share the same
non-significant direction: **lower Trop-2 (and CLDN4 when present) tracks with
better ICI outcome.** BALF (GSE193049) does not. This remains
hypothesis-generating.

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
| 07 | `07_leftover_audit.py` | 对全部剩余/边界系列做元数据与基因存在性深核。 |
| 08 | `08_leftover_analyze.py` | 对能真实连接结局的剩余系列做 TACSTD2/CLDN4 检验。 |

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
均无法提供与每患者 ICI 结局配对的肿瘤上皮批量表达。首轮因此仅分析 **GSE207422**。
剩余轮次补完所有相关系列，并找回因 ICI 正则 `immunotherap\b` 匹配不到英文单词
*immunotherapy* 而漏掉的 **GSE221733**。详见
`notes/fable_geo_2022_2023/leftover_verdicts.md`。

### 4. 纳入分析的队列
| 队列 | 组织 | n（分析） | 表达 | ICI 方案 | 结局 |
|---|---|---|---|---|---|
| **GSE207422**（主要，肺） | NSCLC 基线（治疗前）肿瘤 | 24 | 批量 RNA-seq，log2 TPM | 新辅助抗 PD-1（特瑞普利/信迪利/卡瑞利珠）+ 铂类化疗 | 主要病理缓解（MPR 含 pCR，n=9）vs 非 MPR（NMPR，n=15）；残余肿瘤比例；RECIST |
| GSE243238（交叉验证，**非肺**） | 肢端黑色素瘤 | 14 | 批量 RNA-seq，原始计数 → log2 CPM | ICI | 临床获益 是（n=3）/ 否（n=11） |
| **GSE221733**（剩余，找回） | NSCLC GeoMx DSP PanCK+ | 34 例患者 | NanoString CTA（仅 TACSTD2；无 CLDN4） | 免疫治疗 | 应答 15 vs 非应答 19；OS |
| **GSE248378**（剩余） | NSCLC ICI 后切除肿瘤（非 MPR 集） | 29 | 批量 RNA-seq FPKM | 新辅助 durvalumab ± SBRT | 复发 9 vs 20（源数据连接）；DFS |
| GSE193049（剩余，**BALF 非肿瘤**） | 支气管肺泡灌洗液 | 7 | 批量 RNA-seq | PD-1 阻断 | 应答 4 vs 非应答 3 |

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

**剩余肿瘤检验**（`analysis/leftover_marker_outcome_tests.csv`）：

**GSE221733（NSCLC DSP，仅 TACSTD2，患者水平 PanCK+）：**
| 基因 | 结局 | 应答中位 | 非应答中位 | p | 效应量 |
|---|---|---|---|---|---|
| TACSTD2 | ICI 应答 | 9.268 | 9.378 | 0.5324 | Cliff's δ −0.13 |
| TACSTD2 | OS log-rank（中位分组） | 高 n=17 | 低 n=17 | 0.1545 | — |

**GSE248378（NSCLC ICI 后非 MPR FPKM；论文 n=29）：**
| 基因 | 结局 | 无复发中位 | 复发中位 | p | 效应量 |
|---|---|---|---|---|---|
| TACSTD2 | 复发 | 41.98 | 80.01（FPKM） | 0.0562 | Cliff's δ −0.456 |
| CLDN4 | 复发 | 38.97 | 77.37 | 0.3108 | Cliff's δ −0.244 |
| TACSTD2 | DFS log-rank（中位分组） | 高 n=15 | 低 n=14 | 0.1220 | — |

剔除臂别不一致的 `45-M-PO` 后：TACSTD2 复发 p=0.0766，CLDN4 p=0.3759。

**GSE193049（BALF，非肿瘤，n=7）：** TACSTD2 p=0.6286（δ=+0.333）；CLDN4 p=0.4000（δ=+0.50）——与肿瘤系列方向相反，效能不足。

**肿瘤组织方向一致，仍未达 p < 0.05。** 所有*肿瘤*检验中，更好的 ICI 结局都与更低的
TACSTD2（以及可测时的 CLDN4）相关。剩余系列中最接近显著的是 GSE248378 的
TACSTD2–复发（p=0.0562）。未编造任何统计：GSE248378 在 GEO 中无 MPR 对照
（29/29 均为非 MPR）；GSE221733 无 CLDN4；GSE193049 是 BALF。

### 6. 解读与局限
- 假设（高 TACSTD2/CLDN4 → ICI 应答更差）在**肿瘤组织中方向获支持但未证实**。
  没有一项肿瘤检验达到 p < 0.05。最强剩余信号是 GSE248378 治疗后非 MPR 肿瘤的
  TACSTD2 与复发（p=0.0562），不是基线预测因子。
- GSE221733 仅在剩余核验中找回；CTA 面板无 CLDN4。
- GSE162520 / GSE161537 仍因靶向面板无法用于这两个基因。
- GSE207422 与 GSE248378 均为新辅助背景，而非晚期单药 ICI。
- 全部为作者已处理矩阵，未重新比对。

### 7. 结论
补完 2022–2023 全部剩余 GEO 肺癌 ICI 系列后，可在不编造标签的前提下连接到 ICI
结局的肺癌数据有三个：GSE207422（基线批量，两基因）、GSE221733（DSP PanCK+，
仅 TACSTD2）、GSE248378（治疗后非 MPR 批量，两基因）。肿瘤检验方向一致：
**更低的 Trop-2（及可测时的 CLDN4）与更好的 ICI 结局相关。** BALF（GSE193049）
则不然。仍为产生假设的信号。
