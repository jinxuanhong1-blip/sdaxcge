# TACSTD2 (TROP2) / CLDN4 vs ICI response in Chinese & East-Asian NSCLC cohorts — 中国及东亚 NSCLC 免疫治疗队列中 TACSTD2/CLDN4 与疗效的关联

Parallel slice `fable_china_ici` — 2026-08-16.
Outputs: `notes/fable_china_ici/`, `scripts/fable_china_ici/`, `results/fable_china_ici/`.

---

## English

### 1. Objective

Find publicly available (GEO / GSA-open) processed expression datasets from Chinese or East-Asian NSCLC patients treated with immune-checkpoint inhibitors (ICIs) — prioritizing CameL-, ORIENT- and RATIONALE-related trials — and test whether tumor **TACSTD2 (TROP2)** or **CLDN4** expression associates with ICI response. Separately, catalog GSA-Human / OMIX controlled-access datasets. Total processed downloads were kept < 2 GB (actual ≈ 1.84 GB, of which 1.82 GB is optional scRNA).

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

Open but not analyzable for these genes: GSE179994 (PKU, anti-PD-1, sorted T cells only), GSE136961 (KR, 395-gene immune panel without TACSTD2/CLDN4), GSE225620 (CN, whole blood), GSE309446 (CN, TCR-seq only). Full catalog incl. all notes: `results/fable_china_ici/tables/dataset_catalog.csv`.

### 4. Results

Statistics: two-sided Mann-Whitney U; AUC oriented so that >0.5 = higher expression in responders. Full numbers in `results/fable_china_ici/tables/bulk_stats.csv`, `gse207422_scrna_stats.csv`, `gse241934_scrna_stats.csv`; figures in `results/fable_china_ici/figures/`.

**Pre-treatment tumor bulk RNA (primary question):**

| Dataset | Gene | Median resp. vs non-resp. | p | AUC (resp-high) |
|---|---|---|---|---|
| GSE207422 (CN, MPR vs NMPR) | TACSTD2 | 5.00 vs 6.22 log2TPM | 0.34 | 0.38 |
| GSE207422 | CLDN4 | 5.85 vs 6.27 | 0.26 | 0.36 |
| GSE126044 (KR, R vs NR) | TACSTD2 | 6.17 vs 6.29 log2CPM | 0.44 | 0.36 |
| GSE126044 | CLDN4 | 2.54 vs 3.77 | 0.11 | 0.24 |
| GSE135222 (KR, DCB vs NDB) | TACSTD2 | 7.34 vs 7.87 | 0.61 | 0.43 |
| GSE135222 | CLDN4 | 7.69 vs 7.13 | 0.69 | 0.56 |

GSE135222 survival: median-split logrank p = 0.43 (TACSTD2) / 0.70 (CLDN4); Cox per-unit HR 1.03 / 1.05, both n.s.

**Post-treatment epithelial cells (scRNA, descriptive):** GSE207422 TACSTD2 lower in MPR (p = 0.15, AUC 0.22); GSE241934 IIT slightly higher in MPR (AUC 0.71, p = 0.32) and real-world null (AUC 0.61 / 0.43) — inconsistent, confounded by residual-tumor content after therapy.

**Exploratory plasma exosome (GSE260770):** exosomal TACSTD2 marginally *higher in non-responders* (p = 0.041), but both group medians are 0 (very sparse detection) — hypothesis-generating only.

### 5. Interpretation

- **No open East-Asian cohort supports TACSTD2/CLDN4 as a positive predictor of ICI benefit.** The consistent (though individually non-significant) direction across all three pre-treatment tumor cohorts is *higher TACSTD2/CLDN4 in non-responders* (responder-high AUC 0.24–0.43). This aligns with reports that TROP2-high NSCLC tends toward an immune-excluded phenotype with poorer ICI outcomes.
- For the TROP2-ADC ± IO combination question in Chinese NSCLC (e.g. sacituzumab tirumotecan/SKB264 + pembrolizumab, Dato-DXd + IO), this direction is, if anything, convenient: ICI non-responders may retain or enrich TROP2 target expression, but current open data are far too small to be more than suggestive.
- **Power is the binding constraint**: 16–27 patients per cohort detects only large effects. The adequately powered Chinese datasets (ORIENT-11 n=397 RNA-seq; ORIENT-3 n=110) exist but are locked (Section 2).

### 6. Caveats

- Neoadjuvant cohorts (GSE207422, GSE241934) combine ICI with chemotherapy; pathologic response is not a pure ICI-response readout. Advanced-line Korean cohorts are ICI monotherapy but not Chinese.
- GSE207422/GSE241934 scRNA epithelium is post-treatment (resected residual tumor) in most samples; MPR cases contribute few malignant cells and "Epi" includes normal epithelium. GSE207422 has no per-cell type annotation on GEO, so an EPCAM+/PTPRC− proxy was used.
- GSE135222 DCB threshold (PFS ≥ 180 d) is our derivation from the provided PFS fields, not an author label.
- No multiple-testing correction; all results are descriptive for a two-gene targeted question.

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
bash scripts/fable_china_ici/01_download_data.sh core   # ~10 MB; `all` adds scRNA (~1.8 GB)
python3 scripts/fable_china_ici/02_bulk_analysis.py     # bulk + KM/Cox
python3 scripts/fable_china_ici/03_scrna_gse207422.py   # needs `all` download
python3 scripts/fable_china_ici/04_scrna_gse241934.py   # needs `all` download
```

Python deps: pandas, numpy, scipy, matplotlib, openpyxl, lifelines. Data cached in `$DATA_DIR` (default `/tmp/fable_china_ici_data`); raw downloads are not committed to the repo.

---

## 中文

### 1. 目标

检索中国/东亚 NSCLC 免疫检查点抑制剂（ICI）治疗队列中**公开可下载**（GEO / GSA 开放级）的处理后表达数据——优先关注 CameL、ORIENT、RATIONALE 相关试验——并分析肿瘤 **TACSTD2（TROP2）** 与 **CLDN4** 表达是否与 ICI 疗效相关；同时整理 GSA-Human/OMIX **受控访问**数据集清单。处理后数据下载总量控制在 2 GB 以内（实际约 1.84 GB，其中 1.82 GB 为可选的单细胞矩阵）。

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

开放但不适用于本问题的：GSE179994（北大，仅分选 T 细胞）、GSE136961（韩国，395 基因免疫面板，不含目标基因）、GSE225620（全血）、GSE309446（仅 TCR-seq）。完整清单见 `results/fable_china_ici/tables/dataset_catalog.csv`。

### 4. 结果

统计：双侧 Mann-Whitney U；AUC>0.5 表示应答者表达更高。完整数值见 `results/fable_china_ici/tables/`，图见 `results/fable_china_ici/figures/`。

**治疗前肿瘤 bulk（核心问题）**：三个队列方向一致——**TACSTD2/CLDN4 在无应答者中更高**（应答者高的 AUC 0.24–0.43），但均未达显著（p 0.11–0.69）。GSE135222 生存分析亦为阴性（logrank p=0.43/0.70；Cox HR 1.03/1.05）。

**术后上皮单细胞（描述性）**：GSE207422 中 MPR 者上皮 TACSTD2 偏低（p=0.15）；GSE241934 两个队列方向不一致、均不显著——受术后残留肿瘤比例混杂。

**外泌体探索（GSE260770）**：外泌体 TACSTD2 在无应答者中略高（p=0.041），但两组中位数均为 0（检出极稀疏），仅作提示。

### 5. 解读与局限

- 现有开放东亚数据**不支持** TACSTD2/CLDN4 高表达预测 ICI 获益；一致的趋势反而是无应答者偏高，与 TROP2 高表达偏"免疫排斥"表型的报道一致。对 TROP2-ADC±IO 联用（如芦康沙妥珠单抗 SKB264+帕博利珠单抗）而言，这一方向提示 ICI 无应答人群可能保留靶点表达，但样本量太小，仅为线索。
- 主要局限：各队列 16–27 例，效能不足；新辅助队列为 ICI+化疗的复合疗效；单细胞为术后样本且 GSE207422 缺少逐细胞注释（用 EPCAM+/PTPRC− 代理）；GSE135222 的 DCB 为我们基于 PFS 派生；未做多重检验校正。
- 破局点在受控数据：ORIENT-11（n=397 RNA-seq）与 ORIENT-3（n=110）样本量充足，需走 Vivli/申办方申请通道。

### 6. 受控访问清单

GSA-Human：**HRA001033**（上海肺科，GSE207422 原始数据）、**HRA006493**（同组，>100 万细胞配对治疗前后）、**HRA003360**（浙大二院）、**HRA005191**（北大张泽民实验室，疑为 GSE241934 原始数据）、**HRA004391**（医科院肿瘤医院）；OMIX：**OMIX012889**（ORIENT-11 相关，未备案不可下载）；申办方持有（无库编号）：ORIENT-11（Vivli）、ORIENT-3、CameL、RATIONALE-307。各 DAC 联系方式见 `dataset_catalog.csv`。

### 7. 复现

见英文第 8 节命令；依赖 pandas / numpy / scipy / matplotlib / openpyxl / lifelines；数据缓存于 `$DATA_DIR`（默认 `/tmp/fable_china_ici_data`），原始下载文件不入库。
