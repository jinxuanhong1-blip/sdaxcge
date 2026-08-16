# Spatial transcriptomics methods playbook · 空间转录组方法手册

> Scope · 适用范围: **10x Visium** and **NanoString GeoMx DSP-WTA** (whole
> transcriptome). Targeted imaging panels (**CosMx-1K, Xenium IO**) are out of
> scope **unless both `TACSTD2` and `CLDN4` are on that panel** — see
> [Panel gate](#0-scope--panel-gate--适用范围与面板门控).
>
> 适用于 **10x Visium** 与 **NanoString GeoMx DSP-WTA**（全转录组）。靶向成像面板
> （**CosMx-1K、Xenium IO**）默认不在范围内，**除非该面板同时包含 `TACSTD2` 与
> `CLDN4`**，详见[面板门控](#0-scope--panel-gate--适用范围与面板门控)。

**Central biological question · 核心生物学问题:** characterise the niche /
neighbourhood of **`TACSTD2` (TROP2) / `CLDN4`-high tumour spots** relative to
**T-cell, B-cell, and tertiary-lymphoid-structure (TLS) spots**, and connect
spatial structure to immune-checkpoint-inhibitor (ICI) overall survival when
clinical labels exist.
刻画 **`TACSTD2` (TROP2) / `CLDN4` 高表达肿瘤 spot** 相对 **T 细胞、B 细胞、三级
淋巴结构（TLS）spot** 的生态位/邻域，并在有临床标签时将空间结构与免疫检查点抑制剂
（ICI）总生存关联。

**Non-negotiable guardrail · 不可逾越的底线:** **no invented accessions and no
invented statistics.** Every number in an output comes from a tool run on real
data. Stages with missing inputs are skipped with a clear message, never faked.
**不编造登录号、不编造统计量。** 输出中的每个数值都来自在真实数据上运行的工具。缺少
输入的阶段会明确提示并跳过，绝不伪造。

---

## Table of contents · 目录

0. [Scope & panel gate · 适用范围与面板门控](#0-scope--panel-gate--适用范围与面板门控)
1. [Datasets · 数据集](#1-datasets--数据集)
2. [Environments · 环境](#2-environments--环境)
3. [Visium track · Visium 流程](#3-visium-track--visium-流程)
4. [GeoMx WTA track · GeoMx WTA 流程](#4-geomx-wta-track--geomx-wta-流程)
5. [ICI overall survival · ICI 总生存](#5-ici-overall-survival--ici-总生存)
6. [Reproducibility & guardrails · 复现与底线](#6-reproducibility--guardrails--复现与底线)
7. [References · 参考文献](#7-references--参考文献)

---

## 0. Scope & panel gate · 适用范围与面板门控

**EN.** This playbook targets whole-transcriptome spatial assays. Visium and
GeoMx-WTA measure the entire transcriptome, so `TACSTD2` and `CLDN4` are always
present and the niche analysis is always defined. Targeted imaging panels only
carry a curated gene list; the standard CosMx 1,000-plex and Xenium immuno-
oncology (IO) panels are **not guaranteed** to include both genes. Therefore:

- **Run** the niche analysis on Visium / GeoMx-WTA unconditionally.
- **Run** it on a targeted panel **only if both `TACSTD2` and `CLDN4` are on that
  panel.** If either is missing, the templates raise and stop — no proxy-gene
  substitution, no silent partial analysis.

This is enforced in code by `assert_panel(...)` in
[`templates/_utils.py`](templates/_utils.py) and by `panel_gate` in
[`config/config.yaml`](config/config.yaml).

**中文.** 本手册面向全转录组空间检测。Visium 与 GeoMx-WTA 覆盖全转录组，故 `TACSTD2`
与 `CLDN4` 恒存在，生态位分析恒有定义。靶向成像面板仅含精选基因列表；标准 CosMx
1000-plex 与 Xenium 免疫肿瘤（IO）面板**不保证**同时含这两个基因。因此：

- 对 Visium / GeoMx-WTA **无条件**运行生态位分析；
- 对靶向面板**仅当面板同时含 `TACSTD2` 与 `CLDN4` 时**运行；若缺任一，模板报错停止，
  **不使用替代基因、不做静默的部分分析**。

代码中由 [`templates/_utils.py`](templates/_utils.py) 的 `assert_panel(...)` 与
[`config/config.yaml`](config/config.yaml) 的 `panel_gate` 强制执行。

---

## 1. Datasets · 数据集

Both are real public accessions; details verified from the repositories. See
[`demo/README.md`](demo/README.md) for download and required-file specifics.
两者均为真实公开登录号，细节已从数据库核对。下载与所需文件见 [`demo/README.md`](demo/README.md)。

| Accession · 登录号 | Platform · 平台 | Content · 内容 | Role · 用途 |
| --- | --- | --- | --- |
| **GSE271689** | GeoMx DSP-**WTA** | NSCLC, PD-1 immunotherapy, 586 AOIs; `RAW.tar` ≈ 36 MB | Primary demo (< 2 GB) — end-to-end GeoMx + ICI OS · 主演示（<2GB），端到端 GeoMx + ICI OS |
| **E-MTAB-13530** | 10x **Visium** | Human NSCLC lesions + non-involved lung, 36 sections | Secondary demo — Visium track (full raw may exceed 2 GB; use one section) · 次演示；完整原始可能>2GB，取单切片 |

Why GSE271689 is the driven demo · 为何以 GSE271689 为主演示: its processed GeoMx
archive is ~36 MB (**far under the 2 GB budget**), it is genuinely GeoMx-**WTA**,
and it carries the immunotherapy-outcome labels the ICI OS template needs. The
source study reports its own survival statistics; **this playbook does not reuse
those numbers as if the templates produced them** — you compute your own on the
data you assemble.
GSE271689 的处理档约 36 MB（**远低于 2GB**），确为 GeoMx-**WTA**，且带 ICI OS 模板所需
的免疫治疗结局标签。原研究报告了其自身生存统计量；**本手册不把这些数字当作模板产出复用**，
你应在自建数据上计算自己的结果。

---

## 2. Environments · 环境

Two separate environments (Python for Visium, R/Bioconductor for GeoMx). Tool
generation is 2024–2026.
两套环境（Visium 用 Python，GeoMx 用 R/Bioconductor），工具为 2024–2026 代际。

```bash
# Visium (Python) — templates 01-05 + demo helpers
python -m venv .venv && source .venv/bin/activate
pip install -r methods/spatial/environment/requirements-visium.txt

# GeoMx (R/Bioconductor) — templates 06-09
Rscript methods/spatial/environment/install-geomx.R
```

Key versions · 关键版本 (as of 2026-08): `squidpy` 1.8.2 (Jun 2026, also 1.7.0
line), `cell2location` 0.1.5 (Sep 2025) / 0.1.4 (Oct 2024), `SpatialDE`
(SpatialDE2 API), `scanpy` ≥ 1.10; `GeomxTools`, `SpatialDecon`, `standR` on
Bioconductor ≥ 3.18. GPU strongly recommended for cell2location.
GPU 强烈建议用于 cell2location。

---

## 3. Visium track · Visium 流程

Pipeline: **QC → normalization → deconvolution → niche/neighbourhood → SVG.**
All commands assume you have activated the Python env and are at the repo root.
流程：**质控 → 归一化 → 去卷积 → 生态位/邻域 → 空间可变基因。** 命令均假设已激活 Python
环境且位于仓库根目录。

### 3.1 QC — [`templates/01_visium_qc.py`](templates/01_visium_qc.py)

**EN.** Reads a spaceranger `outs` directory or an `.h5ad`, computes per-spot
`total_counts`, `n_genes_by_counts`, and `pct_counts_mito`, writes the
distributions **before** filtering (so thresholds are auditable), then applies
the `config.yaml:qc.visium` gates and filters low-prevalence genes. Raw counts
are preserved in `layers['counts']` for cell2location and SpatialDE2.
**中文.** 读取 spaceranger `outs` 或 `.h5ad`，计算每个 spot 的总计数、基因数、线粒体
占比，**先**输出过滤前分布（阈值可审计），再按 `config.yaml:qc.visium` 过滤，并去除低
流行度基因；原始 counts 保存在 `layers['counts']` 供 cell2location 与 SpatialDE2 使用。

```bash
python methods/spatial/templates/01_visium_qc.py \
    --input <outs_dir_or.h5ad> --sample S1
```

Inspect `qc_distributions.png` and adjust thresholds to the tissue — the config
values are starting points, not universal truths.
查看 `qc_distributions.png` 并按组织调整阈值——配置值是起点而非定论。

### 3.2 Normalization + clustering — [`templates/02_visium_normalization.py`](templates/02_visium_normalization.py)

**EN.** Default is log-normalization + top-2000 HVGs; analytic Pearson residuals
are available (`normalization.visium.method: pearson_residuals`) and generally
give cleaner count-based clustering. Produces PCA/UMAP/Leiden and per-spot
marker-set scores (`score_tumor_epithelial`, `score_t_cell`, `score_b_cell`,
`score_tls`, …) that the niche step consumes.
**中文.** 默认对数归一化 + 2000 高变基因；可选解析 Pearson 残差，通常聚类更干净。产出
PCA/UMAP/Leiden 及每个 spot 的标记评分，供生态位步骤使用。

```bash
python methods/spatial/templates/02_visium_normalization.py \
    --input <S1_qc.h5ad> --sample S1
```

### 3.3 Deconvolution — [`templates/03_visium_deconvolution_cell2location.py`](templates/03_visium_deconvolution_cell2location.py)

**EN.** cell2location, two stages: (A) estimate reference signatures from an
annotated scRNA-seq `.h5ad`; (B) map cell-type abundances onto Visium spots
(**raw counts required**). This mirrors the E-MTAB-13530 study design. The
per-spot dominant type is stored as `obs['cell_type']`, which the niche step can
use instead of marker scores.
**中文.** cell2location 两阶段：(A) 从带注释 scRNA-seq 估计参考特征；(B) 将细胞类型丰度
映射到 Visium spot（**需原始 counts**）。与 E-MTAB-13530 设计一致。每个 spot 的主导
类型存为 `obs['cell_type']`，可替代标记评分供生态位步骤使用。

```bash
python methods/spatial/templates/03_visium_deconvolution_cell2location.py ref \
    --reference <annotated_scrna.h5ad> --label-key cell_type
python methods/spatial/templates/03_visium_deconvolution_cell2location.py map \
    --input <S1_qc.h5ad> --sample S1
```

If no scRNA reference is available, skip this and drive the niche step from
marker scores (RCTD/stereoscope/destVI are alternatives).
若无 scRNA 参考，可跳过，用标记评分驱动生态位步骤（也可用 RCTD/stereoscope/destVI）。

### 3.4 Niche / neighbourhood — [`templates/04_visium_niche_neighborhood.py`](templates/04_visium_niche_neighborhood.py)

**EN.** The core analysis. Each spot is labelled one of
`{tumor_TACSTD2_CLDN4, T_cell, B_cell, TLS, other}`. A spot earns the tumour
label **only if both `TACSTD2` and `CLDN4` are detected** and the tumour lineage
score wins; otherwise it takes the top of T/B/TLS (if > 0) or `other`. Then:

1. Build the Visium hex-grid graph (`squidpy.gr.spatial_neighbors`, 6 neighbours).
2. **Neighbourhood enrichment** (`nhood_enrichment`, permutation z-scores):
   quantifies whether tumour spots are preferentially adjacent to — or excluded
   from — TLS/B/T spots. Saved as a labelled z-score matrix CSV.
3. **Co-occurrence vs distance** (`co_occurrence`).
4. **Compositional niches**: cluster spots by their neighbour-label fractions
   (KMeans, `k = niche.n_niches`) and export the mean composition per niche, so
   you can read off which niche is the tumour-adjacent TLS/immune interface.

If you ran deconvolution, pass `--cluster-key cell_type` to use abundance-derived
labels instead of marker scores.

**中文.** 核心分析。每个 spot 被标为
`{tumor_TACSTD2_CLDN4, T_cell, B_cell, TLS, other}` 之一。**仅当 `TACSTD2` 与
`CLDN4` 均检出**且肿瘤谱系评分最高时才标为肿瘤；否则取 T/B/TLS 中最高者（>0）或
`other`。随后：

1. 构建 Visium 六边形图（6 邻居）；
2. **邻域富集**（置换 z 分数）：量化肿瘤 spot 是否偏好与 TLS/B/T 相邻或被排斥，输出带
   标签的 z 分数矩阵 CSV；
3. **共现-距离曲线**；
4. **组成型生态位**：按邻居标签比例对 spot 聚类（KMeans），导出每个生态位的平均组成，
   据此读出哪个生态位是肿瘤邻接的 TLS/免疫界面。

如已做去卷积，用 `--cluster-key cell_type` 以丰度标签替代标记评分。

```bash
python methods/spatial/templates/04_visium_niche_neighborhood.py \
    --input <S1_norm_or_c2l.h5ad> --sample S1
# outputs: *_nhood_enrichment_zscore.csv, *_niche_composition.csv, *_niche.h5ad
```

Interpretation · 解读: a **positive** tumour↔TLS z-score with a co-occurrence
peak at short distance indicates a tumour-adjacent lymphoid interface; a
**negative** tumour↔T/B z-score indicates immune exclusion around
`TACSTD2/CLDN4`-high tumour.
肿瘤↔TLS 的 **正** z 分数且短距离共现峰，提示肿瘤邻接淋巴界面；肿瘤↔T/B 的 **负** z
分数，提示 `TACSTD2/CLDN4` 高肿瘤周围的免疫排斥。

### 3.5 Spatially variable genes — [`templates/05_visium_svg_spatialde2.py`](templates/05_visium_svg_spatialde2.py)

**EN.** Primary: **SpatialDE2** omnibus test (on counts) plus optional
expression-based tissue segmentation. Fallback (dependency-light): squidpy
Moran's I autocorrelation on HVGs. Reports FDR-significant genes at
`svg.fdr_alpha`.
**中文.** 主用 **SpatialDE2** omnibus 检验（基于 counts）与可选组织分区；备选 squidpy
Moran's I。按 `svg.fdr_alpha` 报告 FDR 显著基因。

```bash
python methods/spatial/templates/05_visium_svg_spatialde2.py \
    --input <S1_norm.h5ad> --sample S1 --method spatialde2   # or: --method moran
```

---

## 4. GeoMx WTA track · GeoMx WTA 流程

GeoMx AOIs are multi-cellular segments collected per **compartment** (e.g.
tumour / immune / stroma, or PanCK+ vs PanCK-) inside ROIs on patient slides.
The pipeline: **QC + Q3 normalization → compartment models → deconvolution.**
GeoMx 的 AOI 是按**区室**（如 肿瘤/免疫/基质，或 PanCK± ）在患者切片 ROI 内采集的多
细胞片段。流程：**质控 + Q3 归一化 → 区室模型 → 去卷积。**

### 4.1 QC + normalization — [`templates/06_geomx_qc_normalization.R`](templates/06_geomx_qc_normalization.R)

**EN.** Builds a `NanoStringGeoMxSet` from DCC + PKC + annotation, applies
segment (AOI) QC (reads, % trimmed/aligned/saturation, nuclei, area) and probe
QC, aggregates probes to gene targets, computes the **limit of quantification
(LOQ)** from negative probes, filters genes by detection rate, and applies **Q3
(upper-quartile) normalization** — the GeoMx standard.
**中文.** 由 DCC + PKC + 注释构建 `NanoStringGeoMxSet`，做 AOI 质控（读数、修剪/比对/
饱和度、核数、面积）与探针质控，探针聚合到基因，由阴性探针计算**定量下限（LOQ）**，按
检出率过滤基因，并做 **Q3（上四分位）归一化**——GeoMx 标准。

```bash
Rscript methods/spatial/templates/06_geomx_qc_normalization.R \
    --dcc methods/spatial/demo/data/GSE271689 \
    --pkc <path/Hs_R_NGS_WTA_v1.0.pkc> \
    --annotation <path/geomx_annotation.csv> \
    --out methods/spatial/demo/out/geomx
```

> The WTA `.pkc` and the AOI annotation sheet are **not** in the GEO archive; you
> supply them (see [`demo/README.md`](demo/README.md)).
> WTA `.pkc` 与 AOI 注释表**不在** GEO 档案中，需自行提供（见 `demo/README.md`）。

### 4.2 Compartment models — [`templates/07_geomx_compartment_models.R`](templates/07_geomx_compartment_models.R)

**EN.** Differential expression between compartments must respect the nested
design (many AOIs per ROI, many ROIs per patient). Use a **linear mixed model**
with a random intercept per patient/slide — the GeoMx-recommended `mixedModelDE`
— or the `standR` + limma-voom route with `duplicateCorrelation(block=patient)`.
Pick with `--engine {mixed|standr}`.
**中文.** 区室间差异表达须考虑嵌套设计（每 ROI 多 AOI、每患者多 ROI）。采用带患者/切片
随机截距的**线性混合模型**（GeoMx 推荐的 `mixedModelDE`），或 `standR` + limma-voom
（`duplicateCorrelation(block=患者)`）。以 `--engine` 选择。

```bash
Rscript methods/spatial/templates/07_geomx_compartment_models.R \
    --rds methods/spatial/demo/out/geomx/geomx_target_qnorm.rds \
    --compartment segment --contrast Tumor,Immune --subject patient \
    --engine mixed --out methods/spatial/demo/out/geomx
```

For the `TACSTD2/CLDN4` question in GeoMx, contrast the tumour (PanCK+/epithelial)
compartment against immune/TLS compartments and read the tumour-epithelial DE
signature; the resulting per-AOI signature score is a natural input to the ICI
OS template.
GeoMx 中的 `TACSTD2/CLDN4` 问题：将肿瘤（PanCK+/上皮）区室与免疫/TLS 区室对比，读取
肿瘤上皮 DE 特征；所得每 AOI 特征评分是 ICI OS 模板的自然输入。

### 4.3 Deconvolution — [`templates/08_geomx_deconvolution_spatialdecon.R`](templates/08_geomx_deconvolution_spatialdecon.R)

**EN.** `SpatialDecon` estimates cell-type abundances per AOI from a
gene×cell-type profile matrix, using the negative-probe background as the noise
model. Supply a lung/TME profile matrix (SafeTME-style or study-matched
scRNA-derived).
**中文.** `SpatialDecon` 以基因×细胞类型特征矩阵、阴性探针背景为噪声模型，估计每个 AOI
的细胞类型丰度。需提供肺/TME 特征矩阵。

```bash
Rscript methods/spatial/templates/08_geomx_deconvolution_spatialdecon.R \
    --rds methods/spatial/demo/out/geomx/geomx_target_qnorm.rds \
    --profile <cell_profile_matrix.csv> \
    --out methods/spatial/demo/out/geomx
```

### 4.4 One-command guarded demo · 一键守卫式演示

[`demo/run_demo.sh`](demo/run_demo.sh) downloads the GSE271689 DCC files, checks
every prerequisite (Rscript, `GeomxTools`, PKC, annotation), and runs 06→07 only
when all are present. If anything is missing it prints exactly what to supply and
exits **without fabricating output**.
[`demo/run_demo.sh`](demo/run_demo.sh) 下载 GSE271689 DCC，检查全部前置条件，齐备时
才运行 06→07；有缺失则明确提示并退出，**不伪造输出**。

---

## 5. ICI overall survival · ICI 总生存

[`templates/09_ici_os_survival.R`](templates/09_ici_os_survival.R)

**EN.** **Runs only if survival labels exist.** Given a table that joins a
spatial feature (compartment DE signature score, deconvolved cell-type fraction,
or `TACSTD2/CLDN4` niche fraction) to clinical OS, it performs Kaplan–Meier +
log-rank by group, univariable Cox per feature, and a multivariable Cox
(leading feature + covariates). If the `time`/`event` columns are absent it
prints a notice and exits 0 — nothing imputed, nothing invented. Aggregate AOIs
to the patient level (`--id patient`) to avoid pseudo-replication.

**中文.** **仅当存在生存标签时运行。** 给定将空间特征（区室 DE 特征评分、去卷积细胞比例、
或 `TACSTD2/CLDN4` 生态位比例）与临床 OS 关联的表格，执行 KM + log-rank、逐特征单因素
Cox 及多因素 Cox（主特征 + 协变量）。若缺 `time`/`event` 列，则提示并以 0 退出——不插补、
不编造。用 `--id patient` 聚合到患者级，避免伪重复。

```bash
Rscript methods/spatial/templates/09_ici_os_survival.R \
    --table <clinical_features.csv> \
    --time os_months --event os_event --group signature_high \
    --covariates age,sex,stage --id patient \
    --out methods/spatial/demo/out/survival
```

GSE271689 suits this (GeoMx-WTA NSCLC on PD-1 therapy with outcomes). The source
study's published hazard ratios are **their** results; do not copy them into your
report — compute your own on the features and labels you assemble.
GSE271689 契合此步（PD-1 治疗的 GeoMx-WTA NSCLC 含结局）。原研究发表的风险比是**其**
结果，勿抄入你的报告——请在自建特征与标签上计算自己的结果。

---

## 6. Reproducibility & guardrails · 复现与底线

- **No invented accessions / statistics · 不编造登录号/统计量.** Only
  `GSE271689` and `E-MTAB-13530` are referenced; all numbers come from tool runs
  on real data. 仅引用这两个登录号；所有数值来自真实数据的工具运行。
- **Fail loud, never fake · 显式失败，绝不伪造.** Missing inputs → skip with a
  message (see `run_demo.sh`, template 09). 缺输入即带提示跳过。
- **Panel gate enforced in code · 面板门控在代码中强制** (`_utils.assert_panel`,
  `config.panel_gate`).
- **Determinism · 确定性:** `project.random_seed` threads through Leiden, KMeans,
  and permutation tests. 随机种子贯穿 Leiden、KMeans、置换检验。
- **Data hygiene · 数据卫生:** downloads live under `demo/data/` and are
  git-ignored; only code is committed. 下载数据在 `demo/data/` 且被 git 忽略，仅提交代码。
- **Patient-level inference · 患者级推断:** aggregate AOIs/spots to patients
  before survival and cross-sample tests to avoid pseudo-replication.
  生存与跨样本检验前先聚合到患者，避免伪重复。
- **Config-first · 配置优先:** tune everything in
  [`config/config.yaml`](config/config.yaml); templates read from it.
  所有参数在配置文件中调整，模板从中读取。

---

## 7. References · 参考文献

- **E-MTAB-13530** — 10x Visium of human NSCLC. Study: *Single-cell and spatial
  transcriptomics analysis of non-small cell lung cancer*, Nat Commun (2024),
  s41467-024-48700-8; data at ArrayExpress/BioStudies.
- **GSE271689** — GeoMx DSP-WTA, NSCLC, PD-1 immunotherapy; NCBI GEO.
- **squidpy** — Palla et al., Nat Methods (2022); v1.7–1.8.2 (2025–2026).
- **cell2location** — Kleshchevnikov et al., Nat Biotechnol (2022);
  v0.1.4 (2024) / v0.1.5 (2025).
- **SpatialDE2** — Kats, Vento-Tormo, Stegle (2021), bioRxiv
  10.1101/2021.10.27.466045.
- **GeomxTools / GeoMxWorkflows / SpatialDecon / standR** — Bioconductor.
- **Giotto** — Dries et al., Genome Biol (2021) — alternative Visium/imaging
  ecosystem (R).
- Survival modelling — `survival` / `survminer` (Cox PH, Kaplan–Meier).

> Tool versions confirmed against upstream release channels as of 2026-08.
> 工具版本截至 2026-08 已对照上游发布渠道核对。
