# TROP2 / CLDN4 multiplex IF & IMC — methods playbook

# TROP2 / CLDN4 多重免疫荧光与成像质谱流式 — 方法手册

**Scope / 范围:** methods only under `methods/mif/`. No results, no reconstructed images, no claimed re-analysis of Bessede.

**范围说明:** 仅方法，全部放在 `methods/mif/`。无结果、无重建图像、不声称复现 Bessede 分析。

**Last checked / 检索日期:** 2026-08-16.

---

## 0. One-sentence honesty / 一句话诚实声明

**There is no public TROP2+CLDN4 multiplex IF or IMC image set to segment in this repo. Bessede mIF (images and cell-level tables) is not public. If a user-supplied cell table exists, spatial stats can start from that table. If it does not, stop — do not invent pixels or cells.**

**本仓库没有可分割的公开 TROP2+CLDN4 多重 IF / IMC 图像。Bessede 多重免疫荧光（图像与单细胞表）不公开。若使用者提供细胞表，可从表做空间统计；若没有，停止——不要编造像素或细胞。**

---

## 1. What is public vs what is not / 什么公开、什么不公开

### 1.1 Bessede et al. (the paper people usually mean) / Bessede 等人（通常所指的那篇）

Bessede A, Peyraud F, Besse B, et al. *TROP2 is associated with primary resistance to immune checkpoint inhibition in patients with advanced non–small cell lung cancer.* Clin Cancer Res. 2024;30(4):779–785. doi:[10.1158/1078-0432.CCR-23-2566](https://doi.org/10.1158/1078-0432.CCR-23-2566). PMID: 38048058.

| Item / 项目 | Public? / 是否公开 | Notes / 说明 |
|---|---|---|
| Paper + methods text / 论文与方法正文 | Yes (CC BY-NC-ND 4.0) | Enough to know *how they stained and scored*, not enough to re-segment. 足以知道他们如何染色与打分，不足以重新分割。 |
| Supplementary patient tables (S3/S4) / 补充患者表 | Journal supplements | Clinical characteristics of the BIP mIHF subset, **not** a cell table. BIP mIHF 亚组临床特征，**不是**细胞表。 |
| RNA-seq (POPLAR/OAK, BIP) / 转录组 | Restricted (EGA / Vivli) | Transcripts, not images. 转录组，不是图像。 |
| **mIF / mIHF images** | **No** | Explicit data-availability statement (quoted below). 数据可用性声明原文见下。 |
| **inForm / FlowJo cell-level intensity table** | **No** | Not deposited on GEO, Zenodo, Figshare, or EGA as a reusable cell table (checked 2026-08-16). 未作为可复用细胞表存放（2026-08-16 检索）。 |
| CLDN4 channel / CLDN4 通道 | **Never in this panel** | Bessede 5-plex is **PanCK, TROP2, CD8, PD-L1, DAPI**. No CLDN4. 五色面板无 CLDN4。 |

Bessede data-availability statement (verbatim from the paper):

> The immunofluorescence datasets are not publicly available due to information that could compromise research participant consent. According to French/European regulations, any reuse of the data must be approved by the ethics committee. Each request for access to the immunofluorescence dataset (including the images) will be granted after a request is sent to the corresponding author (A. Italiano) and approval by the ethics committee.

**Do not scrape figure panels and treat them as analysis-grade images.** Figure 2A is a representative field, RGB, compressed, and not a spectrally unmixed component TIFF.

**不要把论文插图当分析级图像。** 图 2A 是代表性视野，RGB、压缩，不是光谱拆分后的 component TIFF。

### 1.2 Bessede panel is not a TROP2/CLDN4 panel / Bessede 面板不是 TROP2/CLDN4 面板

What they actually ran (Ventana Discovery XT + Akoya Opal + PhenoImager HT + inForm v2.6.0):

- CD8 — clone C8/144B (Dako/Agilent)
- PanCK — AE1/AE3/PK26 (Ventana)
- PD-L1 — QR-1 (Diagomics)
- TROP2 — EPR2043 (Abcam)
- DAPI (Akoya spectral)

They segmented cells on **DAPI + fluorescent membrane signal**, extracted mean intensities, GaussNorm-normalized (`flowStats`), and thresholded in FlowJo. They scored **total / membrane / intracellular TROP2** separately. Intracellular (not membrane) TROP2 associated with worse ICI outcome in *n* = 50.

A TROP2/CLDN4 multiplex study is a **different experiment**. CLDN4 (tight-junction claudin-4) was not measured. Related biology (TROP2 interacting with claudins / immune exclusion) is discussed elsewhere (e.g. Zhao et al., JITC 2026, TROP2/claudin-7 in TNBC) but that work used HALO on confocal IF, not a public TROP2+CLDN4 IMC/mIF cell table.

TROP2/CLDN4 多重实验是**另一项实验**。CLDN4 未被测量。相关生物学（TROP2 与 claudin / 免疫排斥）见其他工作（如 Zhao 等 JITC 2026，TNBC 中 TROP2/claudin-7），使用 HALO 分析共聚焦 IF，同样**没有**公开的 TROP2+CLDN4 IMC/mIF 细胞表。

### 1.3 Public IMC / mIF that is *not* a substitute / 公开 IMC / mIF **不能**替代

Checked 2026-08-16; none of these is a TROP2+CLDN4 multiplex:

| Dataset / 数据集 | Why it does not replace Bessede or TROP2/CLDN4 / 为何不能替代 |
|---|---|
| Jackson, Fischer et al. 2020 breast IMC (`imcdatasets`) | 35-plex; **no TROP2, no CLDN4**. 35 色；无 TROP2、无 CLDN4。 |
| Other `imcdatasets` (Damond pancreas, Hoch-Schulz melanoma, …) | Standard immune/epithelial panels; TROP2/CLDN4 not part of the curated objects. 标准免疫/上皮面板，整理对象不含 TROP2/CLDN4。 |
| QuPath demo `LuCa-7color` | Teaching image (PD-L1/CD8/FoxP3/CD68/PD-1/CK). Not TROP2/CLDN4, not Bessede. 教学图，不是 TROP2/CLDN4，不是 Bessede。 |
| Zenodo TNBC / metastasis IMC tables | Cohort-specific panels; no identified public table with both TROP2 and CLDN4 columns. 队列特异面板；未发现同时含 TROP2 与 CLDN4 列的公开表。 |

**Conclusion / 结论:** image-first QuPath work in this folder is a **protocol for when you have authorized images**. Table-first spatial stats run **only if** you drop a cell table into `methods/mif/data/`. This repo ships neither.

**本文件夹的 QuPath 流程是“有授权图像时怎么做”的方案。** 表优先空间统计**仅在**把细胞表放入 `methods/mif/data/` 后运行。本仓库两者都不附带。

---

## 2. Two legal entry points / 两条合法入口

```
A. You have images (mIF OME-TIFF / IMC OME-TIFF) after ethics + MTA
   → QuPath (or steinbock) segmentation → cell table → spatial stats
   有伦理与 MTA 后的图像 → QuPath（或 steinbock）分割 → 细胞表 → 空间统计

B. You have only a cell table (x, y, marker intensities, optional phenotype)
   → skip segmentation → spatial stats in python/analyze_cell_table.py
   仅有细胞表 → 跳过分割 → python/analyze_cell_table.py
```

If neither A nor B exists: **methods documentation only**. The analyzer exits with a documented skip.

若 A、B 都没有：**只保留方法文档**。分析脚本以记录在案的 skip 退出。

Put user tables here (gitignored contents except README):

将用户表放在此处（除 README 外内容被 gitignore）：

`methods/mif/data/cells.csv`  
or / 或 `methods/mif/data/cells.parquet`

---

## 3. Recommended panels (if you generate new data) / 若自行产生新数据时的推荐面板

This is **not** Bessede’s panel. It is what a TROP2/CLDN4 multiplex needs if the scientific question is co-expression + immune exclusion.

**这不是 Bessede 面板。** 若科学问题是共表达 + 免疫排斥，TROP2/CLDN4 多重需要这些。

### 3.1 Minimal mIF (5–6 plex) / 最小 mIF（5–6 色）

| Channel / 通道 | Role / 角色 | Compartment to quantify / 定量亚细胞区室 |
|---|---|---|
| DAPI or spectral DAPI | Nuclei / 核 | nucleus |
| PanCK (or E-cadherin) | Tumor epithelium / 肿瘤上皮 | cytoplasm / membrane |
| TROP2 | Target / 靶点 | **membrane, cytoplasm, nucleus separately**（必须分区室；Bessede 的关键发现） |
| CLDN4 | Tight junction / 紧密连接 | **membrane** (expected) / 预期膜 |
| CD8 | Cytotoxic T cells / 细胞毒 T | whole cell |
| Optional: PD-L1 or CD3 | Context / 背景 | membrane / whole cell |

TROP2 is cleaved into ECD (membrane/cytoplasm) and ICD (can accumulate in nucleus). If you only export one “cell mean TROP2”, you cannot test the Bessede-style intracellular hypothesis.

TROP2 裂解为 ECD（膜/胞质）与 ICD（可在核内蓄积）。若只导出一个“细胞平均 TROP2”，无法检验 Bessede 式胞内假说。

### 3.2 IMC notes / IMC 注意

- Nucleus channel is usually **Ir191 / Ir193**, not DAPI. 核通道通常是 **Ir191 / Ir193**，不是 DAPI。
- Pixel size is typically **~1 µm**. Do not run a 0.5 µm fluorescence StarDist model without testing. 像素约 **1 µm**。未经测试不要直接套用 0.5 µm 荧光 StarDist 模型。
- Prefer **steinbock** (Bodenmiller) for MCD → OME-TIFF, spillover compensation, and DeepCell/Mesmer segmentation; QuPath is excellent for QC and export. MCD → OME-TIFF、溢出校正、DeepCell/Mesmer 分割优先 **steinbock**；QuPath 适合质控与导出。
- Metal-tagged TROP2 and CLDN4 clones must be validated (positive/negative cell lines, tonsil/skin controls, isotope purity). 金属标记 TROP2、CLDN4 克隆必须验证。
- No public IMC panel identified that already includes both markers. 未发现已同时包含两标记的公开 IMC 面板。

---

## 4. QuPath: project setup and cell segmentation / QuPath：项目与细胞分割

Scripts live in `methods/mif/qupath/`. They do **nothing** without images you add locally.

脚本在 `methods/mif/qupath/`。没有本地图像时**不会做任何事**。

Official references / 官方参考:

- Multiplex tutorial: <https://qupath.readthedocs.io/en/stable/docs/tutorials/multiplex_analysis.html>
- StarDist: <https://qupath.readthedocs.io/en/stable/docs/deep/stardist.html>
- InstanSeg (newer alternative): <https://qupath.readthedocs.io/en/stable/docs/deep/instanseg.html>
- Bankhead et al. 2017, *Sci Rep* (cite QuPath). 引用 QuPath 请用此文。
- Schmidt et al. 2018, MICCAI (cite StarDist). 引用 StarDist 请用此文。

### 4.1 Project / 项目

1. Create a QuPath project. Add **unmixed component** OME-TIFFs (Akoya) or IMC OME-TIFFs — not RGB snapshots. 建项目。加入**拆分后的 component** OME-TIFF 或 IMC OME-TIFF，不要 RGB 截图。
2. Set image type to **Fluorescence** (use this also for IMC). 图像类型设为 **Fluorescence**（IMC 也用这个）。
3. Rename channels to short marker names (`DAPI`, `TROP2`, `CLDN4`, `PanCK`, `CD8`, …). Classifiers reuse these names. 通道名改为短标记名。分类器会复用这些名字。
4. Annotate analyzable tumor (exclude folds, necrosis, dust, edge flare). Pathologist review if the claim is clinical. 标注可分析肿瘤区（排除折叠、坏死、灰尘、边缘光晕）。若结论偏临床，需病理复核。

Groovy helper: `qupath/00_set_channels.groovy`.

### 4.2 Nucleus / cell segmentation / 核与细胞分割

**Order of preference / 优先顺序**

1. **InstanSeg** (QuPath extension) if available — whole-cell, multiplex-aware. 若可用，优先（全细胞、对多重更友好）。
2. **StarDist** fluorescence model `dsb2018_heavy_augment.pb` on DAPI (mIF) or Ir193 (IMC), then cell expansion. 荧光模型用于 DAPI（mIF）或 Ir193（IMC），再膨胀出细胞。
3. Built-in **Cell detection** (watershed) only as a fallback. 内置分水岭仅作后备。

StarDist template: `qupath/01_stardist_cells.groovy`.

Critical parameters you must tune on *your* images (do not copy-paste as truth):

必须在**自己的图像**上调参（不要把默认值当真理）：

| Parameter / 参数 | mIF starting point / mIF 起点 | IMC starting point / IMC 起点 |
|---|---|---|
| Detection channel / 检测通道 | `DAPI` | `Ir193` (or `DNA1`) |
| `pixelSize` | ~0.5 µm | try 1.0 µm first / 先试 1.0 µm |
| Probability threshold / 概率阈 | 0.5, then inspect / 再目视 | often 0.4–0.6 |
| Cell expansion / 细胞膨胀 | 2–5 µm | 2–4 µm (large pixels) |
| Percentile norm / 百分位归一 | 1–99 | 1–99; IMC dynamic range is ugly / IMC 动态范围很差 |

Expand nuclei to cells so **membrane** measurements exist. TROP2 and CLDN4 are membrane proteins; nuclear-only ROIs will systematically under-call them.

必须把核膨胀为细胞，否则没有**膜**测量。TROP2 与 CLDN4 是膜蛋白；仅核 ROI 会系统性低估。

QC the segmentation before any threshold:

分割质控先于任何阈值：

- Split / merged nuclei in dense tumor. 密集团块中的核分裂/融合。
- Stromal thin nuclei vs tumor. 间质细长核 vs 肿瘤。
- Empty “cells” in necrosis. 坏死区空细胞。
- IMC: laser-ablation stripes and hot pixels. IMC：剥蚀条纹与热像素。

Script `qupath/01_stardist_cells.groovy` can drop detections with impossible nuclear area or eccentricity (edit the gates).

脚本可按核面积/偏心率丢掉不可能的检测（请改阈值）。

### 4.3 Measurements / 测量

After detection, QuPath writes per-compartment intensities for every named channel:

检测后 QuPath 会按区室写出每个命名通道的强度：

- `Nucleus: <Marker> mean` / `max` / `std dev`
- `Cytoplasm: <Marker> mean`
- `Cell: <Marker> mean`
- `Membrane: <Marker> mean` (if membrane is modeled; otherwise use cytoplasm ring as proxy) 若未建模膜，可用胞质环近似

For TROP2, export **all three** of nucleus / cytoplasm (or membrane) / cell. Bessede’s intracellular score is not the same as cell-mean.

TROP2 请导出核 / 胞质（或膜）/ 全细胞。Bessede 的胞内分数 ≠ 细胞均值。

For CLDN4, the biologically expected signal is **membranous**. A useful derived feature (compute after export if QuPath did not):

CLDN4 生物学预期为**膜阳性**。导出后可算：

```
cldn4_membrane_ratio = membrane_CLDN4 / max(cytoplasm_CLDN4, epsilon)
```

High ratio + high membrane intensity → junction-like staining. High cytoplasm / low membrane → mis-segmentation or internalized / non-junction protein.

高比值 + 高膜强度 → 连接样染色。高胞质 / 低膜 → 分割错误或内化 / 非连接蛋白。

### 4.4 Phenotyping in QuPath / 在 QuPath 中表型

Follow the official multiplex tutorial: **one classifier per marker**, then combine.

按官方多重教程：**每个标记一个分类器**，再合并。

1. `Classify → Object classification → Create single measurement classifier`  
   - Tumor: `Cell: PanCK mean` (or cytoplasm).  
   - TROP2 membrane: `Membrane: TROP2 mean` or `Cytoplasm: TROP2 mean`.  
   - TROP2 intracellular: `Nucleus: TROP2 mean` (only if bleed-through from membrane is controlled). 仅在膜串扰可控时使用。  
   - CLDN4: `Membrane: CLDN4 mean`.  
   - CD8: `Cell: CD8 mean` (or cytoplasm).
2. Or train a random-trees object classifier per channel (`Ignore*` vs marker class). 或按通道训练随机树（`Ignore*` vs 标记类）。
3. `Load object classifier` and combine → classes such as `PanCK: TROP2: CLDN4`. 合并后类别形如 `PanCK: TROP2: CLDN4`。

Do **not** use a global intensity cutoff across slides without a transfer rule (percentile, mixture model, or control-tissue anchor). Slide-to-slide Opal/IMC intensity is not a common scale.

**不要**在没有转移规则（百分位、混合模型或对照组织锚点）时用一张片子的强度阈套所有片子。Opal/IMC 强度不是同一标尺。

Export: `qupath/02_export_cell_table.groovy` → TSV/CSV with centroids + measurements + class.

### 4.5 IMC-specific QuPath path / IMC 在 QuPath 中的路径

1. Convert MCD → OME-TIFF (`imctools` / steinbock / Standard BioTools).  
2. Confirm voxel size metadata (~1 µm). 确认体素元数据。  
3. Set Ir193 (or DNA) as the StarDist channel in `01_stardist_cells.groovy`.  
4. Compensate spillover *before* thresholds (steinbock `compensate`). 阈值前做溢出校正。  
5. Hot-pixel / dead-pixel QC; IMC is sparse and noisy. 热像素/死像素质控。

If steinbock already produced a cell table, **skip QuPath segmentation** and go to §5.

若 steinbock 已产出细胞表，**跳过 QuPath 分割**，直接 §5。

---

## 5. Table-first analysis (no images required) / 仅细胞表分析（不需要图像）

### 5.1 Expected schema / 期望表结构

Dictionary: `schema/cell_table_columns.md`. Header-only example: `schema/cell_table.header.csv`.

Minimum useful columns / 最低可用列:

| Column (aliases accepted) / 列名（接受别名） | Required? / 必需？ | Meaning / 含义 |
|---|---|---|
| `cell_id` | yes | Unique within image. 图像内唯一。 |
| `image_id` (or `roi_id`, `sample_id`) | yes | Spatial graph is **per image**. 空间图必须**按图**计算。 |
| `x`, `y` | yes for spatial / 空间分析必需 | Centroids in µm preferred; if pixels, pass `--pixel-size-um`. 质心优先 µm；若是像素，传 `--pixel-size-um`。 |
| `trop2` or compartment columns | for TROP2 questions | Intensity. 强度。 |
| `cldn4` or compartment columns | for CLDN4 questions | Intensity. 强度。 |
| `panck` / `cd8` / `pd_l1` | optional | Lineage / immune. 谱系 / 免疫。 |
| `phenotype` or `class` | optional | If missing, script can threshold. 缺失时脚本可阈值。 |
| `area_um2`, `nucleus_area` | QC | Drop debris / clumps. 去掉碎片/团块。 |

Compartment aliases understood by the script / 脚本认识的区室别名:

`trop2_nucleus`, `trop2_cytoplasm`, `trop2_membrane`, `trop2_cell`, and the same for `cldn4`.

### 5.2 Run / 运行

```bash
# Honest default: looks for methods/mif/data/cells.csv or cells.parquet
# 默认：查找 methods/mif/data/cells.csv 或 cells.parquet
python methods/mif/python/analyze_cell_table.py

# If you have a table somewhere else / 表在别处
python methods/mif/python/analyze_cell_table.py --table /path/to/cells.csv --out methods/mif/out

# Parser smoke test only (synthetic geometry, NOT biology)
# 仅测试解析器（合成几何，不是生物学）
python methods/mif/python/analyze_cell_table.py --synthetic-smoke
```

If no table is found, the script prints the public-data limitation and writes `methods/mif/out/SKIPPED_NO_TABLE.md`. Exit code 0 (skip, not a crash).

若无表，脚本打印公开数据限制并写入 `methods/mif/out/SKIPPED_NO_TABLE.md`。退出码 0（跳过，不是崩溃）。

### 5.3 What the script computes when a table exists / 有表时脚本计算什么

Per image, then optionally aggregated / 先按图，再可选汇总:

1. **QC** — missing coords, duplicate IDs, intensity quantiles, area filters. 缺失坐标、重复 ID、强度分位数、面积过滤。
2. **Phenotype** — use `phenotype`/`class` if present; else per-image thresholds (default: image-wise 75th percentile, override with `--q`). 有表型列则用；否则按图阈值（默认图内 75 分位，`--q` 可改）。
3. **Co-expression** — fractions of tumor cells that are TROP2+, CLDN4+, double-positive. 肿瘤细胞中 TROP2+、CLDN4+、双阳比例。
4. **Compartment scores** (if columns exist) — intracellular TROP2, membrane CLDN4 ratio. 若有区室列：胞内 TROP2、膜 CLDN4 比值。
5. **Spatial stats** (requires x, y) — see §6. 空间统计见 §6。

It does **not** invent survival *p*-values or Bessede cutoffs. Survival needs a clinical table joined on `sample_id`; that join is out of scope until you provide both tables.

脚本**不会**编造生存 *p* 值或 Bessede 截断。生存分析需要按 `sample_id` 连接的临床表；两表都提供之前不做。

---

## 6. Spatial statistics (from a cell table) / 空间统计（基于细胞表）

All distances in **micrometres**. Compute **within one image / ROI**. Never pool (x, y) across slides.

距离一律用**微米**。在**单张图 / ROI 内**计算。绝不要把不同切片的 (x, y) 拼成一张点图。

### 6.1 Questions that match TROP2/CLDN4 biology / 与 TROP2/CLDN4 生物学匹配的问题

| Question / 问题 | Metric / 指标 | Needs / 需要 |
|---|---|---|
| Do TROP2 and CLDN4 sit on the same tumor cells? 是否在同一肿瘤细胞上？ | Double-positive fraction; Spearman of intensities on PanCK+ cells | intensities, not just x,y |
| Is CLDN4 membranous where TROP2 is membranous? 两者是否同为膜阳性？ | Joint membrane ratios | compartment columns |
| Are CD8 cells excluded from TROP2+CLDN4+ tumor cores? CD8 是否被排斥在双阳肿瘤核心外？ | NN distance CD8 → tumor subset; G-function; mixing score | x,y + phenotypes |
| Do double-positive cells cluster? 双阳细胞是否聚集？ | Ripley *K* / pair correlation vs CSR | x,y + phenotype |

Without CD8 (or another T-cell marker) you **cannot** test immune exclusion. A TROP2+CLDN4-only table can still test co-expression and clustering.

没有 CD8（或其他 T 细胞标记）**不能**检验免疫排斥。仅 TROP2+CLDN4 的表仍可检验共表达与聚集。

### 6.2 Metrics implemented / 已实现指标

See `python/analyze_cell_table.py`.

**Nearest-neighbor (NN) distance / 最近邻距离**

For each CD8 cell, distance to nearest `TROP2+CLDN4+` tumor cell (and, as a control, to nearest any-tumor cell). Summarize as median NN per image. Compare the two medians: if CD8 is farther from double-positive tumor than from tumor in general, that is *consistent with* local exclusion — not proof.

对每个 CD8，算到最近 `TROP2+CLDN4+` 肿瘤细胞的距离（对照：到任意肿瘤）。按图取中位 NN。若 CD8 离双阳肿瘤比离一般肿瘤更远，**符合**局部排斥，不是证明。

**G-function (empirical) / 经验 G 函数**

Fraction of CD8 cells that have ≥1 target cell within radius *r* (default *r* ∈ {10, 20, 30, 50} µm). Permute CD8 labels among non-tumor cells (or all cells) *n* times to get a null envelope.

半径 *r* 内至少 1 个靶细胞的 CD8 比例。将 CD8 标签在非肿瘤（或全部）细胞中置换 *n* 次得到零假设包络。

**Mixing / isolation score / 混合 / 隔离分数**

Among *k* nearest neighbors of each double-positive tumor cell, fraction that are CD8. Low mixing + large NN = exclusion-like. Sensitive to density: always report tumor and CD8 densities (cells / mm²).

每个双阳肿瘤细胞的 *k* 近邻中 CD8 比例。低混合 + 大 NN = 排斥样。对密度敏感：务必同时报告肿瘤与 CD8 密度（细胞 / mm²）。

**Ripley *K* (border-corrected, translation) / Ripley *K*（平移边界校正）**

Clustering of the double-positive point pattern versus complete spatial randomness. Implemented as a simple translation-edge correction; for publication-grade *K* use `spatstat` in R (`Kest`, `Lest`, `pcf`).

双阳点格局相对完全空间随机的聚集。此处为简单平移校正；发表级 *K* 请用 R `spatstat`（`Kest`, `Lest`, `pcf`）。

### 6.3 What not to do / 不要做的事

- Do not compute Ripley *K* on a whole TMA punch treated as one rectangle if the tissue is a crescent — use the tissue mask or at least the convex hull of cells. 组织呈新月形时不要把整个 TMA 孔当矩形算 *K*——用组织掩膜，至少用细胞凸包。
- Do not interpret a 20 µm G-function on IMC (1 µm pixels, ~cell-sized) the same as on 0.25 µm mIF. 不要把 IMC（1 µm 像素）上的 20 µm G 函数与 0.25 µm mIF 同等解释。
- Do not threshold TROP2 on all cells including stroma and call it “TROP2-high tumor”. Restrict to PanCK+ (or equivalent). 不要在包括间质的全部细胞上阈值 TROP2 再称为“TROP2 高肿瘤”。限制在 PanCK+。
- Do not copy Bessede’s maxstat PFS cutoff onto another cohort. maxstat is anti-conservative without validation. 不要把 Bessede 的 maxstat PFS 截断抄到另一队列。未经验证的 maxstat 偏乐观。

### 6.4 R extras (optional, not required here) / 可选 R 补充

If you already live in R: `spatstat`, `imcRtools` (`colPair`, `aggregateNeighbors`), `spicyR`, `Statial`. These need the same cell table. This playbook’s runnable path is Python so the folder stays dependency-light.

若已在 R 中：`spatstat`、`imcRtools`、`spicyR`、`Statial`。需要同一张细胞表。本手册可运行路径用 Python，以保持依赖轻量。

---

## 7. Intensity normalization and cutoffs / 强度归一与截断

Bessede: `GaussNorm` (`flowStats`) then FlowJo thresholds. That is a **per-batch cytometry-style** norm, not a spatial method.

Bessede：`GaussNorm`（`flowStats`）再 FlowJo 阈值。这是**批次内细胞术风格**归一，不是空间方法。

Practical rules if you only have a table / 只有表时的可行规则:

1. Prefer **within-image ranks or z-scores** for co-expression plots. 共表达图优先**图内秩或 z**。
2. If multiple batches, normalize **within batch** before any global cutoff. 多批次先**批内**再全局截断。
3. Publish the exact measurement name (`Nucleus: TROP2 mean` vs `Cell: TROP2 mean`). 发表时写明测量名。
4. If you must dichotomize: pre-specify percentile or control-tissue threshold; lock it before outcome tests. 若必须二分：预先规定百分位或对照组织阈值；在结局检验前锁定。
5. Membrane vs intracellular TROP2 must use **different** measurements. 膜 vs 胞内 TROP2 必须用**不同**测量。

---

## 8. Quality-control checklist / 质控清单

**Images (if you have them) / 若有图像**

- [ ] Unmixed components, not RGB JPEG. 拆分通道，不是 RGB JPEG。
- [ ] Channel names match the panel. 通道名与面板一致。
- [ ] Tumor annotation excludes artifact. 肿瘤标注排除伪影。
- [ ] Segmentation reviewed on ≥3 fields (tumor, stroma, edge). 至少 3 个视野复核分割。
- [ ] TROP2/CLDN4 look membranous on known-positive epithelium. 已知阳性上皮上呈膜染色。
- [ ] No obvious spectral bleed of TROP2 into DAPI (would fake “nuclear TROP2”). TROP2 无明显串到 DAPI（会假造“核 TROP2”）。

**Table / 表**

- [ ] One row = one cell; IDs unique per image. 一行一细胞；图内 ID 唯一。
- [ ] x, y in a documented unit. x, y 单位已写明。
- [ ] Marker columns distinguished from QC columns. 标记列与质控列分开。
- [ ] Sample-level metadata in a **separate** table (do not repeat survival 10⁵ times). 样本元数据**另表**（不要把生存重复 10⁵ 次）。

**Stats / 统计**

- [ ] Spatial metrics per image, then a patient-level summary (mean / median of images). 空间指标按图，再汇总到患者。
- [ ] Densities reported with every proximity metric. 每个邻近指标都带密度。
- [ ] Permutation or analytic null for clustering claims. 聚集结论要有置换或解析零假设。

---

## 9. What this folder will never do / 本文件夹永远不会做的事

- Download or reconstruct Bessede mIF. 不下载、不重建 Bessede mIF。
- Treat paper figures as raw data. 不把论文图当原始数据。
- Pretend Jackson-Fischer IMC or LuCa-7color is TROP2/CLDN4. 不把公开 IMC/教学图假装成 TROP2/CLDN4。
- Report “replication” of Bessede intracellular-TROP2 survival without their images **and** clinical table. 没有他们的图像**和**临床表，不报告胞内 TROP2 生存“复现”。
- Ship a fake cell table labeled as BIP/Bessede. 不把假细胞表标成 BIP/Bessede。

A `--synthetic-smoke` flag exists only to test that the parser and KD-tree code run. The points are random. They are not a biological result.

`--synthetic-smoke` 只用于测试解析器与 KD 树能跑。点是随机的，不是生物学结果。

---

## 10. How to request Bessede mIF (if you actually need it) / 若确实需要，如何申请 Bessede mIF

1. Email the corresponding author (A. Italiano, as stated in the paper) requesting the immunofluorescence dataset **including images** and, separately, the inForm/FlowJo cell table if it exists. 按论文联系通讯作者，申请免疫荧光数据**含图像**，并另询 inForm/FlowJo 细胞表是否存在。
2. Expect ethics-committee approval under French/EU rules; budget time. 按法国/欧盟规则预留伦理审批时间。
3. Ask explicitly for: unmixed component TIFFs, pixel size, channel map, tumor annotations, cell table with membrane/nucleus/cytoplasm TROP2, and the exact FlowJo gates. 明确索取：拆分 TIFF、像素尺寸、通道图、肿瘤标注、含膜/核/胞质 TROP2 的细胞表、以及 FlowJo 门控。
4. Even after approval you still **do not have CLDN4**. A TROP2/CLDN4 question needs new staining or a different cohort. 批准后你仍然**没有 CLDN4**。TROP2/CLDN4 问题需要新染色或另一队列。

---

## 11. File map / 文件地图

```
methods/mif/
  playbook.md                 ← this file / 本文件
  data/README.md              ← drop cells.csv here / 把 cells.csv 放这里
  schema/cell_table_columns.md
  schema/cell_table.header.csv
  qupath/00_set_channels.groovy
  qupath/01_stardist_cells.groovy
  qupath/02_export_cell_table.groovy
  python/analyze_cell_table.py
  python/requirements.txt
```

---

## 12. References / 参考文献

1. Bessede A, et al. Clin Cancer Res. 2024;30:779–785. doi:10.1158/1078-0432.CCR-23-2566. **mIF not public.**
2. Bankhead P, et al. QuPath: Open source software for digital pathology image analysis. Sci Rep. 2017;7:16878.
3. Schmidt U, et al. Cell detection with star-convex polygons. MICCAI 2018.
4. QuPath multiplexed analysis tutorial. <https://qupath.readthedocs.io/en/stable/docs/tutorials/multiplex_analysis.html>
5. QuPath StarDist. <https://qupath.readthedocs.io/en/stable/docs/deep/stardist.html>
6. Windhager J, et al. An end-to-end workflow for multiplexed image processing and analysis. Nat Protoc. 2023 (steinbock / IMC).
7. Jackson HW, Fischer JR, et al. The single-cell pathology landscape of breast cancer. Nature. 2020;578:615–620. **Public IMC; no TROP2/CLDN4.**
8. Zhao et al. TROP2/claudin program mediates immune exclusion… J Immunother Cancer. 2026;14:e012265. **CLDN7 confocal + HALO; not a public TROP2/CLDN4 IMC table.**
9. Baddeley A, Rubak E, Turner R. *Spatial Point Patterns: Methodology and Applications with R.* CRC, 2015 (`spatstat`).
10. Bressan D, et al. Protocol for whole-slide image analysis of human multiplexed tumor tissues using QuPath and R. STAR Protoc. 2024.
