# TROP2 / CLDN4 lung IHC–mIF scoring tables — public-data slice

**Scope.** Public TROP2 or CLDN4 immunohistochemistry (IHC) or multiplex immunofluorescence (mIF) scoring tables in lung cancer, with or without immune-checkpoint inhibitor (ICI) annotation. Outputs only under `notes/grok_ihc/`, `scripts/grok_ihc/`, `results/grok_ihc/`.

**检索范围。** 公开的肺（±ICI）TROP2 或 CLDN4 IHC / mIF 评分表（期刊增补、Figshare、Zenodo、PMC、HPA）。产出仅落在上述三个目录。

---

## English

### Headline

**Bessede mIF is not public.** The 5-plex TROP2 mIHF images and quantitative scores from Bessede et al. (*Clin Cancer Res* 2024; 10.1158/1078-0432.CCR-23-2566) are withheld. The paper states that immunofluorescence datasets are not publicly available because they could compromise research-participant consent; reuse requires a request to A. Italiano and ethics-committee approval. Transcriptomes sit on EGA `EGAS00001005013` (controlled). Extra clinical data go through Vivli/Roche.

What *is* public on AACR Figshare are **clinical summary tables only** (S1–S4). Supplementary Table S3 describes the n=50 BIP mIHF cohort (age, sex, PS, PD-L1, RECIST, DCB) and does **not** contain per-patient TROP2 membrane / intracellular / nuclear scores.

No other public **patient-level** TROP2 mIF scoring table in lung ± ICI was found. No public **patient-level** CLDN4 IHC/mIF table with ICI outcomes was found.

### What was downloaded (open)

| Source | Marker / assay | Lung ± ICI | Table type | Local path |
|---|---|---|---|---|
| Bessede 2024 Figshare S1–S4 | TROP2 context; **not** mIF scores | ICI (OAK/POPLAR RNA; BIP mIHF n=50 clinical) | Aggregated clinical / Cox | `results/grok_ihc/raw/bessede_ccr2024/`, `extracted/bessede_ccr2024/` |
| Hashimoto 2025 *Sci Rep* + MOESM1 | TROP2 IHC (intensity 0–3 × proportion 1–4; OE = product 12) | **Yes — Nivo-Ipi n=110** | Aggregated Table 1–2 + Suppl A1 | `raw/hashimoto2025_scirep/`, `extracted/hashimoto2025_scirep/` |
| Dum 2022 *Pathobiology* Table 1 | TROP2 IHC TMA, 4-tier | Lung ADC/SCC/SCLC/NET; no ICI | Category % (not IPD) | `extracted/dum2022_trop2_tma/` |
| Kuo 2025 *PLoS ONE* + S1–S7 TIF | Trop-2 IHC SP295 H-score 0–300 | NSCLC n=107 + 158; PD-L1 in set 1; IPD request-only | Summary stats + figures | `raw/kuo2025_plosone/` |
| Pak 2012 *WJSO* | TROP2 IHC total score 0–12; OE >4 | Resected AdC 100 / SCC 64; no ICI | Aggregated Tables 1–2 | `raw/pak2012_wjso/` |
| Omori 2021 supplements | TROP2 IHC 0–3; high = 3 | Lung; 5 ICI cases mentioned (no change) | Aggregated high vs low | `raw/omori2021_trop2_lung/` |
| HPA `pathology.tsv` v23 | TACSTD2 + CLDN4 ordinal IHC | Lung cancer n=10 / 11; no ICI | Cancer-type counts | `extracted/hpa/` |
| Gyulai 2023 *POR* | CLDN4 H-score 0–300 | Rare lung ACC+MEC n=35; no ICI | Aggregated medians | `raw/gyulai2023_cldn_rare_lung/` |
| Moldvay 2016 *POR* | CLDN1–7 IHC 0–5 (incl. CLDN4) | Stage I NSCLC n=137; no ICI | Aggregated comparisons | `raw/moldvay2016_cldn_nsclc/` |

Scripts: `scripts/grok_ihc/download_open_tables.py`, `retry_and_extract.py`, `build_catalog.py`. Full inventory: `results/grok_ihc/inventory/source_catalog.csv`.

### Restricted / missing (named so they are not mistaken for open)

- **Bessede mIF / mIHF scores and images — not public** (see headline).
- **Zenodo 10.5281/zenodo.18543127** (also 18494664): AI TROP-2 membrane/cytoplasm H-scores + TMA images for ~1,142 resected NSCLCs. Record is visible; **files are restricted**. `cohort_data.csv` download returned a login page.
- **Kuo 2025 individual H-scores**: exist internally; PLOS data policy points to `datarequest@gilead.com`. Supplements are figures, not xlsx IPD.
- **Dum 2022 raw TMA IPD**: “available upon reasonable request.”
- **Inamura 2017 Oncotarget** (ADC 172/270 high, SqCC 150/201, HGNET 21/115): CC-BY paper; publisher and PMC PDF were 403 / interstitial at fetch. No IPD table on Figshare.
- **Jung 2014 CLDN4 lung ADC**: open paper exists; PMC fetch did not yield a usable score table. No Figshare/Zenodo IPD.
- **TROPION-Lung01 QCS / NMR** and **EVOKE-02 Trop-2 H-score** analyses: conference/press summaries only; no public patient-level scoring tables.

### Scoring systems actually used (lung)

- **H-score 0–300** (Kuo SP295; Gyulai CLDN4): \(\sum\) intensity(1–3) × % cells.
- **Product score 0–12** (Pak; Hashimoto uses intensity × proportion with OE at 12): proportion 0–4 × intensity 0–3.
- **4-tier TMA** (Dum): negative / weak / moderate / strong from intensity + % rules.
- **3-bin intensity** (Inamura): 0 / 1 / 2 with a high vs no-low cut.
- **Ordinal 0–5** (Moldvay CLDN panel).
- **HPA**: High / Medium / Low / Not detected (tiny lung n).
- **Bessede mIHF** (not public): total / membrane / intracellular TROP2; only intracellular associated with worse PFS (median 1.4 vs 13.4 mo) in the paper’s figures.

### Usable public numbers (not IPD)

**TROP2 + ICI (aggregated).** Hashimoto 2025: 46/110 (41.8%) TROP2-overexpressed (product = 12) on Nivo-Ipi; overexpression associated with worse PFS/OS, stronger in PD-L1 <1% and non-adenocarcinoma. Suppl Table A1 is the only downloaded open table that crosses TROP2 IHC with ICI-era PD-L1/histology.

**TROP2 IHC prevalence (no ICI IPD).** Dum 2022 analyzable cores: lung ADC 181 — 6.1% neg / 4.4% weak / 33.7% mod / 55.8% strong; lung SCC 72 — 0 / 12.5 / 6.9 / 80.6; SCLC 11 — 27.3 / 63.6 / 9.1 / 0. Kuo: set 1 median H-score 68 (0–251); set 2 medians 150 (ADC) and 170 (SCC). Pak: TROP2 OE 23% AdC vs 64.1% SCC. HPA TACSTD2 lung: 0 high / 2 med / 3 low / 5 not detected (n=10; antibody-dependent and not aligned with Dum/Kuo prevalence).

**CLDN4.** No ICI-linked public score table. HPA lung: 0 high / 9 med / 2 low / 0 not detected (n=11). Gyulai: CLDN4 H-scores in 35 rare salivary-type lung tumors (grade 2 vs 1 median 60 vs 1); no IPD. Moldvay: CLDN4 IHC 0–5 differs ADC vs lepidic ADC (p=0.001) in 137 stage I cases; no IPD. Diagnostic CLDN4 literature (carcinoma vs mesothelioma) is binary positivity, not H-score IPD.

### Bottom line

Open tables that can be reused immediately are **aggregated**. The only open ICI-era TROP2 **protein** table is Hashimoto 2025 (plus Bessede S3 clinical metadata without scores). Bessede mIF remains closed. The largest advertised patient-level TROP2 H-score file (Zenodo ADC/AI cohort) is restricted. CLDN4 lung IPD with ICI does not appear in public supplements or Figshare/Zenodo as of this slice.

---

## 中文

### 结论先行

**Bessede 的 mIF 不公开。** Bessede 等（*Clin Cancer Res* 2024）5 色 mIHF（PanCK / TROP2 / CD8 / PD-L1 / DAPI）的图像与定量评分因知情同意不能公开；须向通讯作者 A. Italiano 申请并经伦理委员会批准。转录组在 EGA `EGAS00001005013`（受控）。额外临床数据走 Vivli / Roche。

AACR Figshare 上公开的 S1–S4 **只是临床汇总表**。S3 是 BIP 研究 mIHF 亚组（n=50）的基线特征，**没有**逐例 TROP2 膜 / 胞内 / 核评分。

未发现其他公开的、带逐例患者 TROP2 mIF 评分的肺 ± ICI 表。也未发现带 ICI 结局的逐例 CLDN4 IHC/mIF 评分表。

### 已下载的开放表

- **Bessede Figshare S1–S4**：OAK/POPLAR 特征与多因素（TACSTD2 **RNA**）、BIP mIHF/血浆临床表。见 `results/grok_ihc/extracted/bessede_ccr2024/`。
- **Hashimoto 2025**（Nivo-Ipi，n=110）：TROP2 IHC，强度 0–3 × 比例 1–4，过表达定义为乘积 12（46/110 = 41.8%）。增补 Table A1 已下载。这是本切片中**唯一**同时具备 TROP2 **蛋白**评分与 ICI 的开放表（仍为汇总，非 IPD）。
- **Dum 2022 TMA Table 1**：肺腺癌 / 鳞癌 / SCLC / 肺 NET 的四档阳性率（非 IPD）。已转录至 `extracted/dum2022_trop2_tma/`。
- **Kuo 2025**：SP295 膜 H-score 汇总（中位 68 与 150/170）；增补为图（TIF），个例数据需向 Gilead 申请。
- **Pak 2012、Omori 2021**：开放汇总表；Omori 提到 5 例 ICI 前后 TROP2 无变化。
- **HPA pathology.tsv v23**：TACSTD2 / CLDN4 肺癌序数 IHC 计数（样本很少）。
- **Gyulai 2023、Moldvay 2016**：CLDN4（及 claudin 组）H-score / 0–5 分汇总；无 ICI、无 IPD。

完整目录：`results/grok_ihc/inventory/source_catalog.csv`。

### 明确不开放 / 未拿到

- Bessede mIF 图像与评分：**不公开**。
- Zenodo `10.5281/zenodo.18543127`：约 1142 例切除 NSCLC 的 AI TROP-2 H-score + TMA 图，**文件受限**。
- Kuo 个例 H-score、Dum 原始 TMA IPD：需申请。
- TROPION-Lung01 QCS/NMR、EVOKE-02 Trop-2 H-score：仅有会议摘要，无公开评分表。
- Inamura 2017、Jung 2014：论文开放，但本次未能稳定下到可用 PDF/IPD。

### 可直接引用的公开数字（均为汇总）

- **ICI + TROP2 蛋白**：Hashimoto，过表达 41.8%，与 Nivo-Ipi 后更差的 PFS/OS 相关，在 PD-L1 <1% 与非腺癌中更明显。
- **TROP2 流行率（无 ICI IPD）**：Dum 肺腺癌 181 例可评，强阳 55.8%；肺鳞癌 72 例，强阳 80.6%；SCLC 几乎无强阳。Kuo 任一水平阳性约 82–90%。Pak 过表达腺癌 23%、鳞癌 64.1%。
- **CLDN4**：HPA 肺癌 11 例中 9 例中等、2 例低、0 例高。Gyulai 罕见唾液腺型肺癌 H-score 中位数（2 级 vs 1 级：60 vs 1）。Moldvay I 期 NSCLC 中 CLDN4 在腺癌与贴壁为主腺癌之间有差异。均无 ICI 个例表。

### 一句话

能立刻复用的开放资源几乎全是**汇总表**；ICI 场景下唯一开放的 TROP2 **蛋白**表是 Hashimoto 2025。Bessede mIF 不公开。最大的逐例 TROP2 H-score 文件（Zenodo）受限。带 ICI 的 CLDN4 肺评分 IPD 在公开增补 / Figshare / Zenodo 中未出现。
