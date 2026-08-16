# TACSTD2 (TROP2) / CLDN4 vs ICI response in Chinese & East-Asian NSCLC cohorts — 中国及东亚 NSCLC 免疫治疗队列中 TACSTD2/CLDN4 与疗效的关联

Parallel slice `fable_china_ici` — 2026-08-16.
Outputs: `notes/fable_china_ici/`, `scripts/fable_china_ici/`, `results/fable_china_ici/`.

---

## English

### 1. Objective

Find publicly available (GEO / GSA-open) processed expression datasets from Chinese or East-Asian NSCLC patients treated with immune-checkpoint inhibitors (ICIs) — prioritizing CameL-, ORIENT- and RATIONALE-related trials — and test whether tumor **TACSTD2 (TROP2)** or **CLDN4** expression associates with ICI response. Separately, catalog GSA-Human / OMIX controlled-access datasets. Only open processed matrices were analyzed (no HRA downloads). Total processed downloads were kept < 2 GB (core ~12 MB; optional scRNA ~1.82 GB).

### 2. Landmark-trial data availability (key negative finding)

None of the major Chinese registrational ICI trials has open molecular data:

| Trial | Drug | Data status |
|---|---|---|
| **ORIENT-11** (NCT03607539) | sintilimab | Tumor whole-transcriptome RNA-seq exists (MHC-II biomarker paper, PMID 34358724) but is only requestable via **Vivli** (Lilly), with genomic data excluded from standard requests. A related processed set, **OMIX012889** (SYSUCC, n=171 expression profiling), is controlled and not yet HGRAC-registered — currently not downloadable. |
| **ORIENT-3** | sintilimab | RNA-seq of 110 biomarker-evaluable patients with response labels was analyzed in Cancer Communications 2022 (doi 10.1002/cac2.12385); no public deposition. |
| **CameL** (NCT03134872) | camrelizumab | Data on request to corresponding author only; explicitly subject to HGRAC (中国人类遗传资源管理) restrictions. |
| **RATIONALE-307** (NCT03594747) | tislelizumab | EdgeSeq Precision IO panel GEP exists; requests via BeOne Medicines (ClinicalTrials@beonemed.com). No repository deposition. |

Consequence: response-annotated *processed* data from Chinese ICI-treated NSCLC are only available from investigator-initiated / academic cohorts, mostly neoadjuvant.

### 3. Open datasets analyzed

| Dataset | Center | Setting | Data | Response label | n used |
|---|---|---|---|---|---|
| **GSE207422** | Shanghai Pulmonary Hospital (CN) | neoadjuvant anti-PD-1 (toripalimab / camrelizumab / sintilimab) + chemo, resectable NSCLC | bulk RNA-seq log2TPM (pre-treatment biopsies) + scRNA-seq | MPR(incl. pCR) vs NMPR | 24 bulk (9/15); 12 scRNA post-tx |
| **GSE241934** (NEOTIDE/CTONG2104) | Peking University (CN) | neoadjuvant sintilimab + chemo (EGFR-mut IIT, n=11) and real-world EGFR-WT immunochemo (n=34) | scRNA-seq of post-treatment resected tumors, author cell-type annotations | MPR/pCR vs non-MPR | 11 + 24 patients with ≥20 epithelial cells |
| **GSE126044** | Yonsei (KR) | anti-PD-1 monotherapy, advanced NSCLC | bulk RNA-seq counts → log2CPM | responder vs non-responder | 16 (5/11) |
| **GSE135222** | KR | anti-PD-1/PD-L1, advanced NSCLC | bulk RNA-seq (normalized) | PFS event/time; DCB = PFS ≥ 180 d | 27 (7/20) |
| **GSE260770** | Guangzhou Medical University (CN) | sintilimab for high-risk GGO lesions | **plasma exosomal** mRNA FPKM (exploratory; not tumor) | lesion response R vs NR | 50 (26/24) |
| **GSE329813** | Peking Union Medical College (CN) | NCT05383716 neoadjuvant pembrolizumab + platinum-doublet | GeoMx DSP Q3-norm (1812-gene panel; TACSTD2 present, **CLDN4 absent**) | MPR vs NMPR (11/11) | 22 patients (tumor-bed ROI means) |

Open but not used for these genes: GSE179994 (PKU, T cells only), GSE136961 (KR, 395-gene panel without TACSTD2/CLDN4), GSE225620 (CN, whole blood), GSE309446 (CN, TCR-seq), **GSE243013** (PKU, 234-pt immune-only scRNA, 6.6 GB), GSE291670 (Tianjin, n=6 scRNA), GSE229353 (Tianjin, CD45+ only), GSE146100 (GLCI, 1 patient), GSE337519 (Shanghai Chest, 1 patient). Full catalog: `results/fable_china_ici/tables/dataset_catalog.csv`.

### 4. Results

Statistics: two-sided Mann-Whitney U; AUC oriented so that >0.5 = higher expression in responders. Stouffer combination of signed z from independent pre-treatment tumor cohorts. Full numbers in `results/fable_china_ici/tables/` (`bulk_stats.csv`, `bulk_stats_extended.csv`, `stouffer_meta.csv`, scRNA tables); figures in `results/fable_china_ici/figures/`.

**Pre-treatment tumor bulk RNA (primary question):**

| Dataset | Gene | Median resp. vs non-resp. | p | AUC (resp-high) |
|---|---|---|---|---|
| GSE207422 (CN, MPR vs NMPR) | TACSTD2 | 5.00 vs 6.22 log2TPM | 0.34 | 0.38 |
| GSE207422 | CLDN4 | 5.85 vs 6.27 | 0.26 | 0.36 |
| GSE207422 (same 24, RECIST CR/PR vs SD) | TACSTD2 | 6.05 vs 6.22 | 0.76 | 0.55 |
| GSE207422 RECIST | CLDN4 | 5.87 vs 6.49 | 0.32 | 0.36 |
| GSE126044 (KR, R vs NR) | TACSTD2 | 6.17 vs 6.29 log2CPM | 0.44 | 0.36 |
| GSE126044 | CLDN4 | 2.54 vs 3.77 | 0.11 | 0.24 |
| GSE135222 (KR, DCB vs NDB) | TACSTD2 | 7.34 vs 7.87 | 0.61 | 0.43 |
| GSE135222 | CLDN4 | 7.69 vs 7.13 | 0.69 | 0.56 |

GSE135222 survival: median-split logrank p = 0.43 (TACSTD2) / 0.70 (CLDN4); Cox per-unit HR 1.03 / 1.05, both n.s.

**Stouffer meta (3 independent pre-tx tumor cohorts, N=67):** TACSTD2 Z = −1.29, p = 0.20; CLDN4 Z = −1.33, p = 0.18. Direction is consistently *lower* expression in responders, but not significant. Forest: `figures/forest_auc_pretreatment.png`.

**Post-treatment GeoMx (GSE329813, Chinese, n=22; TACSTD2 only):** patient-mean tumor-bed TACSTD2 is **significantly higher in NMPR** (median 4.80 vs 3.94; MW p = 0.002; responder-high AUC = 0.11). ROI-level sensitivity p = 0.029. **This is the only nominally significant tumor-tissue result**, but it is post-resection after neoadjuvant chemo-IO: MPR beds have less residual epithelium, so TACSTD2 (an epithelial/TROP2 marker) is expected to drop with tumor clearance. Treat as residual-tumor content, not a pre-treatment predictive biomarker.

**Post-treatment epithelial cells (scRNA, descriptive):** GSE207422 TACSTD2 lower in MPR (p = 0.15, AUC 0.22); GSE241934 IIT slightly higher in MPR (AUC 0.71, p = 0.32) and real-world null (AUC 0.61 / 0.43) — inconsistent, same residual-tumor confound.

**Exploratory plasma exosome (GSE260770):** exosomal TACSTD2 marginally *higher in non-responders* (p = 0.041), but both group medians are 0 (very sparse detection) — hypothesis-generating only.

### 5. Interpretation

- **No open East-Asian *pre-treatment* cohort supports TACSTD2/CLDN4 as a positive predictor of ICI benefit.** Direction is consistently *higher in non-responders* (responder-high AUC 0.24–0.43) but Stouffer meta is n.s. (p ≈ 0.18–0.20, N=67).
- The only significant tumor-tissue signal is **post-treatment GeoMx TACSTD2 higher in NMPR (GSE329813, p=0.002)**. This is most parsimoniously residual epithelium after MPR, not a predictive biomarker. It does, however, confirm that TROP2 remains abundant in ICI-refractory residual tumor — relevant to TROP2-ADC sequencing after IO.
- For TROP2-ADC ± IO in Chinese NSCLC (sacituzumab tirumotecan/SKB264 + pembrolizumab, Dato-DXd + IO), open data remain too small to claim a predictive cutoff.
- **Power is the binding constraint.** Adequately powered Chinese RNA-seq (ORIENT-11 n=397; ORIENT-3 n=110) exists but is locked (Section 2). GSE243013 (n=234 immune scRNA) is open but immune-only and 6.6 GB.

### 6. Caveats

- Neoadjuvant cohorts (GSE207422, GSE241934, GSE329813) combine ICI with chemotherapy; pathologic response is not a pure ICI-response readout. Advanced-line Korean cohorts are ICI monotherapy but not Chinese.
- GSE329813 / GSE207422 / GSE241934 post-treatment assays are confounded by residual-tumor fraction. GSE329813 GeoMx has no CLDN4. GSE207422 scRNA has no per-cell type annotation on GEO (EPCAM+/PTPRC− proxy).
- GSE135222 DCB threshold (PFS ≥ 180 d) is our derivation from the provided PFS fields, not an author label.
- No multiple-testing correction; all results are descriptive for a two-gene targeted question. Stouffer meta assumes independence and does not weight by n.

### 7. Controlled-access list (GSA-Human / OMIX / sponsor-held)

For requests (see `dataset_catalog.csv` for DAC contacts):

- **HRA001033** — Shanghai Pulmonary Hospital; raw scRNA (15 pts) + bulk RNA (24 pts) of GSE207422; DAC zhangpeng1121@tongji.edu.cn.
- **HRA006493** — same group; >1M cells, paired pre/post PD-1 ± chemo or anti-VEGFA; PRJCA022724.
- **HRA003360** — Zhejiang University 2nd Affiliated Hospital; neoadjuvant anti-PD-1 NSCLC; PRJCA009311.
- **HRA005191** — Peking University (Zemin Zhang lab); post-treatment scRNA+TCR, response/resistance; likely raw counterpart of GSE241934.
- **HRA004391** — National Cancer Center/CAMS; LUAD neoadjuvant immunotherapy scRNA.
- **OMIX012889** — SYSUCC; ORIENT-11 expression profiling + clinical (n=171); controlled, not yet HGRAC-registered.
- **Sponsor-held (no accession):** ORIENT-11 (Vivli), ORIENT-3 (Innovent), CameL (Hengrui), RATIONALE-307 (BeOne Medicines).

### 8. Reproduction

```bash
bash scripts/fable_china_ici/01_download_data.sh core   # ~12 MB; `all` adds scRNA (~1.8 GB)
python3 scripts/fable_china_ici/02_bulk_analysis.py     # bulk + KM/Cox
python3 scripts/fable_china_ici/05_combined_analysis.py # GSE329813 + RECIST + Stouffer/forest
python3 scripts/fable_china_ici/03_scrna_gse207422.py   # needs `all` download
python3 scripts/fable_china_ici/04_scrna_gse241934.py   # needs `all` download
```

Python deps: pandas, numpy, scipy, matplotlib, openpyxl, lifelines. Data cached in `$DATA_DIR` (default `/tmp/fable_china_ici_data`); raw downloads are not committed to the repo.

---

## 中文

### 1. 目标

检索中国/东亚 NSCLC 免疫检查点抑制剂（ICI）治疗队列中**公开可下载**（GEO / GSA 开放级）的处理后表达数据——优先关注 CameL、ORIENT、RATIONALE 相关试验——并分析肿瘤 **TACSTD2（TROP2）** 与 **CLDN4** 表达是否与 ICI 疗效相关；同时整理 GSA-Human/OMIX **受控访问**数据集清单。仅分析开放处理后矩阵（不下载 HRA）。处理后数据下载总量控制在 2 GB 以内（核心约 12 MB；可选单细胞约 1.82 GB）。

### 2. 注册试验数据可及性（重要阴性结论）

中国主要 ICI 注册试验均**无开放分子数据**：

- **ORIENT-11**（信迪利单抗）：肿瘤全转录组 RNA-seq 存在（MHC-II 生物标志物论文，PMID 34358724），但仅可通过礼来 **Vivli** 平台申请，且基因组学数据被排除在常规申请之外；相关处理后数据 **OMIX012889**（中山大学肿瘤防治中心，n=171）为受控访问且尚未完成人遗办（HGRAC）备案，目前无法下载。
- **ORIENT-3**（信迪利单抗）：110 例带疗效标注的 RNA-seq 已用于发表（Cancer Communications 2022，doi 10.1002/cac2.12385），但未存入公共库。
- **CameL**（卡瑞利珠单抗，NCT03134872）：数据仅可向通讯作者申请，明确受人类遗传资源管理规定约束。
- **RATIONALE-307**（替雷利珠单抗）：EdgeSeq Precision IO 面板表达数据存在，需向百济/BeOne（ClinicalTrials@beonemed.com）申请，无公共库存档。

因此，带疗效标注的中国 ICI-NSCLC 开放数据目前只来自研究者发起的（多为新辅助）学术队列。

### 3. 纳入分析的开放数据集

- **GSE207422**（上海市肺科医院）：可切除 NSCLC 新辅助抗 PD-1（特瑞普利/卡瑞利珠/信迪利单抗）+化疗；治疗前活检 bulk RNA-seq log2TPM（n=24，MPR 9 / NMPR 15）+ 单细胞（多为术后样本）。
- **GSE241934**（北京大学，NEOTIDE/CTONG2104）：新辅助信迪利单抗+化疗（EGFR 突变 IIT 队列 n=11）及真实世界 EGFR 野生型免疫化疗队列（n=34）的术后单细胞数据，含作者细胞类型注释与 MPR 标注。
- **GSE126044**（韩国延世）：晚期 NSCLC 抗 PD-1 单药，n=16（R 5 / NR 11），bulk counts。
- **GSE135222**（韩国）：晚期 NSCLC 抗 PD-1/PD-L1，n=27，含 PFS 事件与时间（DCB 定义为 PFS ≥ 180 天）。
- **GSE260770**（广州医科大学附一院）：多原发肺癌高危 GGO 病灶信迪利单抗治疗，**血浆外泌体** mRNA（探索性，非肿瘤组织），n=50。
- **GSE329813**（北京协和医学院）：NCT05383716 新辅助帕博利珠单抗+铂类双药，GeoMx DSP（1812 基因；含 TACSTD2，**无 CLDN4**），22 例术后肿瘤床（MPR 11 / NMPR 11）。

开放但未用于本问题：GSE179994（仅 T 细胞）、GSE136961（免疫面板无目标基因）、GSE225620（全血）、GSE309446（仅 TCR）、**GSE243013**（北大 234 例免疫细胞 scRNA，6.6 GB）、GSE291670（天津，n=6）、GSE229353（仅 CD45+）、GSE146100 / GSE337519（各 1 例）。完整清单见 `dataset_catalog.csv`。

### 4. 结果

统计：双侧 Mann-Whitney U；AUC>0.5 表示应答者表达更高。完整数值见 `results/fable_china_ici/tables/`，图见 `results/fable_china_ici/figures/`。

**治疗前肿瘤 bulk（核心问题）**：三个队列方向一致——**TACSTD2/CLDN4 在无应答者中更高**（应答者高的 AUC 0.24–0.43），但均未达显著。GSE207422 的 RECIST（CR/PR vs SD）同样阴性。GSE135222 生存分析阴性。**Stouffer 合并（N=67）**：TACSTD2 Z=−1.29，p=0.20；CLDN4 Z=−1.33，p=0.18。

**术后 GeoMx（GSE329813，中国，n=22）**：肿瘤床患者均值 TACSTD2 在 NMPR 中显著更高（p=0.002，AUC=0.11）。这是唯一达到名义显著的肿瘤组织结果，但最合理的解释是 **MPR 后残留上皮减少**（TROP2 为上皮标志），不能当作治疗前预测标志物。它确实说明 ICI 难治残留灶仍高表达 TROP2，对 TROP2-ADC 序贯有提示意义。

**术后上皮单细胞（描述性）**：GSE207422 / GSE241934 方向不一致、均不显著，同样受残留肿瘤比例混杂。

**外泌体探索（GSE260770）**：外泌体 TACSTD2 在无应答者中略高（p=0.041），两组中位数均为 0，仅作提示。

### 5. 解读与局限

- 现有开放**治疗前**东亚数据**不支持** TACSTD2/CLDN4 高表达预测 ICI 获益；趋势为无应答者偏高，但合并检验不显著。术后 GeoMx 的显著结果受残留肿瘤混杂。
- 对 TROP2-ADC±IO（如芦康沙妥珠单抗 SKB264+帕博利珠单抗）而言，开放数据只能提示难治残留灶仍表达靶点，不能给出预测阈值。
- 主要局限：各队列 16–27 例；新辅助为 ICI+化疗复合疗效；GSE329813 无 CLDN4；GSE207422 单细胞无逐细胞注释；GSE135222 的 DCB 为派生；未做多重检验校正。
- 破局点仍在受控数据：ORIENT-11（n=397）与 ORIENT-3（n=110）。GSE243013（n=234）虽开放但是免疫细胞且超过 2 GB。

### 6. 受控访问清单

GSA-Human：**HRA001033**（上海肺科，GSE207422 原始数据）、**HRA006493**（同组，>100 万细胞配对治疗前后）、**HRA003360**（浙大二院）、**HRA005191**（北大张泽民实验室，疑为 GSE241934 原始数据）、**HRA004391**（医科院肿瘤医院）；OMIX：**OMIX012889**（ORIENT-11 相关，未备案不可下载）；申办方持有（无库编号）：ORIENT-11（Vivli）、ORIENT-3、CameL、RATIONALE-307。各 DAC 联系方式见 `dataset_catalog.csv`。

### 7. 复现

见英文第 8 节命令；依赖 pandas / numpy / scipy / matplotlib / openpyxl / lifelines；数据缓存于 `$DATA_DIR`（默认 `/tmp/fable_china_ici_data`），原始下载文件不入库。
