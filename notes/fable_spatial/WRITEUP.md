# Spatial lung slice — TACSTD2/CLDN4 vs immune neighborhoods / 肺癌空间转录组切片 — TACSTD2/CLDN4 与免疫邻域

Parallel slice `fable_spatial`. Outputs live only under `notes/fable_spatial/`, `scripts/fable_spatial/`, `results/fable_spatial/`. Raw/processed inputs (~1.1 GB, under the 2 GB budget) are downloaded to `/workspace/data/fable_spatial/` and are **not** committed.

并行切片 `fable_spatial`。所有产出仅位于 `notes/fable_spatial/`、`scripts/fable_spatial/`、`results/fable_spatial/`。原始/处理后输入数据（约 1.1 GB，低于 2 GB 预算）下载到 `/workspace/data/fable_spatial/`，**不**入库。

---

## 1. Dataset verification / 数据集核验

| Dataset | Platform | Content | TACSTD2+CLDN4 in panel? | Decision |
|---|---|---|---|---|
| [GSE271689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271689) | GeoMx DSP **WTA** (~18k genes), 586 AOI DCC files | Advanced NSCLC, first-line PD-1-based ICI; Yale discovery + Greek validation cohorts; paper: [Nat Genet 2025, s41588-025-02351-7](https://www.nature.com/articles/s41588-025-02351-7) | Yes (WTA; both genes present in all compartment matrices) | **Included** |
| GSE221733 (UQ arm of same study) | GeoMx DSP **CTA** targeted panel (1,827 genes) | UQ validation cohort | **No — TACSTD2 present, CLDN4 absent** | **Skipped** per instruction "skip targeted panels unless gene list includes both" |
| [E-MTAB-13530](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530) | 10x Visium | 8 NSCLC tumours (20 sections), 8 matched non-involved lung (16 sections), 2 healthy donors (4 sections); `filtered_feature_bc_matrix.h5` + `spatial.tar` per section (~622 MB total) | Yes (whole transcriptome) | **Included** (all 40 sections pass QC) |
| [GSE189487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189487) | 10x Visium | LUAD progression atlas, 6 sections: TD5/TD8 = AIS, TD3/TD6 = MIA, TD1/TD2 = IAC (PMID 36434043); RAW tar 220 MB | Yes | **Included** |

**Raw-data consistency check for GSE271689 / GSE271689 原始数据一致性核验.** We re-derived the patient-level tumor-compartment (PanCK/CK) matrix from the 586 GEO `.dcc` files (probe→gene via NanoString `Hs_R_NGS_WTA_v1.0.pkc`, MD5 `46ed2c63…` matching the Zenodo 12752405 copy; Q3 normalization per AOI; mean of CK AOIs per patient `spotid`). Against the published patient-level matrix (Nat Genet source data Fig. 6b): mean per-patient log-expression Pearson r = **0.94** (range 0.73–0.99, 37 shared patients, 10,069 shared genes); across patients Spearman rho = **0.963 (TACSTD2)** and **0.960 (CLDN4)**. GEO 数据与发表矩阵高度一致，验证通过。Files: `gse271689_dcc_vs_published_*.csv`.

Clinical outcomes: Yale per-patient OS/PFS come embedded in the paper's source-data sheets (Fig. 6b/6f); Greek OS (`OS_1` days, `Death`) in Fig. 6d/6h sheets. GEO metadata itself carries only spotid/segment/treatment. 临床随访（OS）取自论文 source data 工作表，GEO 元数据本身不含随访。

## 2. Methods / 方法

### GeoMx (GSE271689)

- Yale (n=37 first-line ICI patients; tumor CK compartment, Q3-normalized, log2(x+1)); Greek (n=61; sheets are batch-corrected and contain negative values → rank-based inverse-normal transform).
- Stromal immune context ("immune neighborhood" at AOI level): mean z-score of available T-effector genes in the paired stromal/leukocyte compartment (Yale: CD8A, GZMB, NKG7, CXCL9, CXCL13, CD2, TRAC; Greek: 12 genes).
- Stats: Spearman correlation tumor gene vs stromal score; Cox PH (lifelines) per SD, univariable and adjusted for the stromal T-effector score; Kaplan–Meier median split + log-rank.
- 方法：肿瘤区室基因表达与配对基质区室 T 效应细胞评分做 Spearman 相关；Cox 比例风险（每 SD），单变量及校正基质免疫评分；中位数分组 KM + log-rank。

### Visium (E-MTAB-13530, GSE189487)

Per section: spots ≥250 UMI; CP10K + log1p; epithelial score (EPCAM, KRT7/8/18/19/17, CDH1, ELF3, SFN — deliberately excludes TACSTD2/CLDN4 and close paralogs to avoid circularity); broad immune score (23 lineage markers) and T-effector score (14 genes). Neighborhood = mean immune score of the ≤6 nearest spatial neighbor spots (self excluded, KD-tree on pixel coordinates, distance-capped). In epithelial spots (top 50% epithelial score): **partial Spearman** of TACSTD2/CLDN4 vs neighborhood immune score, controlling for the spot's own epithelial and immune scores (removes the trivial "more epithelium = fewer immune neighbors" and self-composition effects). Per-section naive permutation p (B=1000) is reported for transparency; **primary inference is cross-section**: Wilcoxon signed-rank on per-section rho + signed Stouffer combination, stratified by tissue class.

每张切片：spot ≥250 UMI，CP10K+log1p 标准化；上皮评分（不含 TACSTD2/CLDN4 及近缘旁系，避免循环论证）；广义免疫评分与 T 效应评分。邻域 = 最近 ≤6 个相邻 spot（不含自身）的免疫评分均值。在上皮 spot（上皮评分前 50%）中做**偏 Spearman**（控制自身上皮/免疫评分）。主要推断在切片层面（Wilcoxon 符号秩 + 带符号 Stouffer 合并），避免 spot 层面空间自相关导致的假阳性。

## 3. Results / 结果

### 3.1 GeoMx ICI cohorts: OS association (GSE271689)

Cox OS, HR per SD of tumor-compartment expression (`gse271689_os_cox.csv`):

| Cohort | Gene | Univariable HR (95% CI), p | Adjusted for stromal T-eff, p |
|---|---|---|---|
| Yale (n=37, 17 events) | TACSTD2 | **1.69 (1.01–2.81), p=0.045** | **2.01 (1.14–3.56), p=0.016** |
| Yale | CLDN4 | 1.50 (0.96–2.36), p=0.076 | 1.57 (0.99–2.49), p=0.057 |
| Greek (n=61, 33 events) | TACSTD2 | 0.86 (0.61–1.20), p=0.37 | 0.88 (0.62–1.23), p=0.45 |
| Greek | CLDN4 | 0.86 (0.61–1.22), p=0.40 | 0.85 (0.59–1.20), p=0.35 |

KM median split (`gse271689_os_km_median_split.csv`, figure `fig_gse271689_yale_km.png`): Yale CLDN4-high log-rank p=0.024 (median OS 606 d vs not reached); TACSTD2-high p=0.12; Greek both null.

Tumor-compartment TACSTD2/CLDN4 were **not** correlated with stromal T-effector activity (Yale rho≈0.14, p≈0.43; Greek rho≈0.01, ns; `gse271689_tumor_vs_stromal_immune_corr.csv`) — i.e., their prognostic signal is not a proxy for immune infiltration.

**Interpretation / 解读:** 在 Yale 发现队列中，肿瘤区室高 TACSTD2（TROP2）与 ICI 治疗后较差 OS 相关，且独立于基质 T 效应活性；CLDN4 呈同向边缘显著。但在 Greek 验证队列（批次校正数据、更晚期人群）**未复现**（HR 甚至 <1，均不显著）。该信号只能视为“发现集阳性、验证失败”，不可作外推结论。This is a discovery-only signal that failed external validation; treat accordingly.

### 3.2 Visium: TACSTD2/CLDN4 vs immune neighborhoods

Cross-section meta (`emtab13530_meta.csv`, `gse189487_meta.csv`, figure `fig_visium_neighborhood_rhos.png`); median per-section partial rho [sections positive / total], Wilcoxon p:

| Dataset / tissue | Gene | Broad immune neighborhood | T-effector neighborhood |
|---|---|---|---|
| E-MTAB-13530 tumour (20 sect.) | TACSTD2 | −0.010 [5/20], p=0.097 (Stouffer p=0.0027, negative) | +0.007 [11/20], p=0.35 |
| E-MTAB-13530 tumour | CLDN4 | −0.013 [9/20], p=0.55 | **+0.019 [14/20], p=0.014** (Stouffer p=1.5e−6) |
| E-MTAB-13530 non-involved (16) | TACSTD2 | **+0.039 [12/16], p=0.034** | **+0.033 [13/16], p=4.3e−4** |
| E-MTAB-13530 non-involved | CLDN4 | +0.019 [10/16], p=0.40 | +0.018 [11/16], p=0.058 |
| E-MTAB-13530 donor (4) | both | positive (Stouffer p<0.005) | positive |
| GSE189487 LUAD (6) | TACSTD2 | −0.016 [1/6], p=0.0625 (Stouffer p=0.0077, negative) | −0.001 [3/6], p=1.0 |
| GSE189487 LUAD | CLDN4 | −0.063 [1/6], p=0.0625 (Stouffer p=6e−6, negative) | +0.003 [4/6], p=0.84 |

(For n=6 sections, p=0.0625 is the smallest attainable Wilcoxon value; 5/6 or 6/6 sections were negative.)

**Effect sizes are uniformly small (|median partial rho| ≤ 0.06).** Direction flips by tissue and immune axis: within **tumours**, TACSTD2/CLDN4-high epithelial spots sit in slightly **immune-poorer broad neighborhoods** (myeloid/B/T lineages combined; consistent across both Visium datasets), while CLDN4 shows a weak **positive** coupling to T-effector-rich neighborhoods in E-MTAB tumours. In **non-malignant lung** (non-involved, donor), associations are weakly positive throughout. 结论：肿瘤内高 TACSTD2/CLDN4 的上皮 spot 邻域广义免疫评分略低（两套 Visium 数据方向一致，但效应量极小），不构成强“免疫排斥”表型；CLDN4 与 T 效应邻域呈弱正相关。正常肺组织中两基因与免疫邻域弱正相关。

### 3.3 Expression context / 表达背景

- Paired tumour vs non-involved epithelium (E-MTAB-13530, n=8 patients; `emtab13530_tumour_vs_noninvolved.csv`, figure `fig_emtab13530_expression_by_tissue.png`): TACSTD2 mean log1p 1.35 vs 0.51, **Wilcoxon p=0.0078**; CLDN4 0.98 vs 0.47, p=0.055. TROP2 在肿瘤上皮显著上调，支持其作为 ADC 靶点的组织学选择性。
- GSE189487 stage trend AIS→MIA→IAC (n=6 sections): TACSTD2 rho=0.48, p=0.34; CLDN4 rho=0.00 — no significant stage trend at this n. 无显著分期趋势（样本量小）。

## 4. Limitations / 局限

1. Visium spots are multicellular mixtures; partial correlation mitigates but cannot remove composition confounding. Visium spot 为多细胞混合物，偏相关只能部分校正。
2. Within-section permutation p-values ignore spatial autocorrelation (anti-conservative); conclusions rest on cross-section replication tests. 切片内置换检验未考虑空间自相关，推断以跨切片检验为准。
3. Greek GeoMx sheets are batch-corrected (negative values); rank-based transforms were used, but attenuation of a true signal cannot be excluded. Greek 队列矩阵经批次校正，可能稀释真实信号。
4. Small cohorts (Yale n=37, 17 events); median-split KM is descriptive; no formal multiple-testing correction across the 2 genes × 2 axes × tissue strata — nominal p-values reported. 样本量小，未做全局多重检验校正，均为名义 p 值。
5. UQ CTA cohort skipped (CLDN4 not in panel), so ICI-OS replication relies on the Greek cohort only. UQ 靶向面板缺 CLDN4 被跳过，OS 验证仅剩 Greek 队列。

## 5. Reproduction / 复现

```bash
pip install scanpy h5py lifelines statsmodels pandas scipy openpyxl matplotlib
# data download (GEO tars, E-MTAB h5+spatial, PKC from Zenodo, paper source data) as in scripts
python3 scripts/fable_spatial/run_gse271689_geomx.py
python3 scripts/fable_spatial/run_gse271689_dcc_verify.py
python3 scripts/fable_spatial/run_emtab13530_visium.py
python3 scripts/fable_spatial/run_gse189487_visium.py
python3 scripts/fable_spatial/make_figures.py
```

Scripts: `scripts/fable_spatial/` (`common.py` holds marker sets, loaders, neighborhood + stats). Result tables/figures: `results/fable_spatial/`.
