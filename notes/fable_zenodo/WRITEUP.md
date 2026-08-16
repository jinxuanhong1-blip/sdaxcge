# Open lung-ICI / TACSTD2 / CLDN4 matrices from Zenodo — search, download & analysis

> Bilingual write-up. English first, 中文在后。
> All artifacts live under `notes/fable_zenodo/`, `scripts/fable_zenodo/`,
> `results/fable_zenodo/`.

---

## English

### 1. Goal

Search Zenodo / figshare / OSF for **open, processed matrices** relevant to
**lung immune-checkpoint-inhibitor (ICI) therapy** and to the epithelial /
antibody-drug-conjugate targets **TACSTD2 (TROP2)** and **CLDN4**. Download the
open, processed files under 2 GB, skip restricted Zenodo records, verify DOIs,
and analyze what is feasible.

### 2. Search

`scripts/fable_zenodo/search_repositories.py` queried the Zenodo REST API,
figshare v2 search, and OSF v2 nodes with 8 query strings. It returned **217
unique records** (119 Zenodo, 98 figshare; OSF matched nothing on title). After
title-based relevance scoring (`filter_candidates.py`), 40 relevant records
exposed a downloadable matrix-like file. Full details:
`results/fable_zenodo/search_results.json`, `shortlist.json`, and
`notes/fable_zenodo/search_summary.md`.

### 3. DOI verification

Every candidate DOI was resolved through `https://doi.org/<doi>`; the two used
datasets plus two additional on-topic candidates all return HTTP 200. Machine
record: `results/fable_zenodo/doi_verification.json`.

### 4. Downloads (open, CC-BY-4.0, < 2 GB)

Manifest with sizes + MD5: `results/fable_zenodo/download_manifest.json`.

| Zenodo DOI | Dataset | Files pulled | Type |
|------------|---------|--------------|------|
| [10.5281/zenodo.10731914](https://doi.org/10.5281/zenodo.10731914) | scRNA-seq of lung adenocarcinoma (LUAD) patients & mice | `cell_matrix_sce_8LUADs.csv.gz` (58 MB), `cell_annotation_sce_8LUADs.csv` (9.6 MB) | scRNA-seq counts, 28,457 human cells |
| [10.5281/zenodo.8041882](https://doi.org/10.5281/zenodo.8041882) | Acquired resistance to **anti-PD1** therapy in NSCLC | `sce-clustering-allsamples-original_final.RDS` (160 MB), `AK_List.csv` | Imaging Mass Cytometry, 41 markers × 99,659 cells |

**Restricted Zenodo records were skipped** (4 relevant ones, listed in the search
summary and `doi_verification.json`).

### 5. Analysis

#### 5a. TACSTD2 & CLDN4 in human LUAD scRNA-seq (Zenodo 10731914)

Script: `scripts/fable_zenodo/analyze_luad.py`. The matrix is genes × cells (raw
UMI); we normalized per cell to CP10K using `nCount_RNA` and summarized by the
`Cell_type.unimodel` annotation (28,457 cells across 7 compartments).

Both genes are **strongly and specifically epithelial (tumor-cell) restricted**:

| Gene | Epithelial % positive | Epithelial mean CP10K | Next-highest compartment |
|------|----------------------|----------------------|--------------------------|
| **TACSTD2** | 73.5 % | 10.18 | Myeloid 5.3 % / 0.18 |
| **CLDN4** | 80.0 % | 6.09 | Fibroblast 2.5 % / 0.08 |

Immune (T/NK, B, myeloid, mast) and stromal (endothelial, fibroblast) cells are
essentially negative (< ~5 % positive, mean CP10K < 0.2). This matches known
biology: TACSTD2/TROP2 and CLDN4 are epithelial surface antigens, which is why
TROP2 is the target of the ADC sacituzumab govitecan and CLDN4 is an emerging
solid-tumor target. In an ICI context, this tumor-cell-restricted expression is
what makes them attractive, largely immune-cell-sparing, targets.

Tables: `results/fable_zenodo/luad_expression_by_celltype.csv`,
`luad_expression_global.csv`. Figure: `luad_tacstd2_cldn4_dotplot.png`.

#### 5b. Anti-PD1 NSCLC micro-environment (Zenodo 8041882)

Script: `scripts/fable_zenodo/analyze_ici_imc.py`. The RDS is an S4
`SingleCellExperiment`; no R is installed, so it was parsed directly in Python
with the `rdata` package. It turned out to be an **Imaging Mass Cytometry (IMC)**
object — a **41-protein panel over 99,659 single cells** (7 patients), annotated
with `TherapyStatus` (pre/post anti-PD1), `ImmuneStatus` (T-cell high/low) and a
cluster label. **TACSTD2/CLDN4 are transcriptomic markers and are not on this
targeted protein panel, so they cannot be measured here**; instead the dataset
characterizes the ICI-treated tumor immune micro-environment.

Comparing mean arcsinh expression **post vs pre anti-PD1** (46,638 pre / 53,021
post cells):

- **Up post-therapy:** CD15 (+0.48), CD45RO memory-T (+0.41), CD134/OX40 (+0.32),
  CD73 (+0.25), Carbonic Anhydrase IX / hypoxia (+0.23), Granzyme B (+0.13).
- **Down post-therapy:** CD163 (−0.15), CD14 (−0.11), Siglec-1/CD169 (−0.11) —
  i.e. a drop in M2-like / suppressive myeloid markers.

The direction is consistent with the paper's theme (an immunosuppressive shift on
acquired anti-PD1 resistance, with CD73/OX40/memory-T remodeling). Tables:
`results/fable_zenodo/ici_imc_marker_pre_post.csv`, `ici_imc_panel.csv`,
`ici_imc_summary.json`. Figure: `ici_imc_pre_post_markers.png`.

### 6. Reproduce

```bash
pip install pandas scipy h5py matplotlib rdata
python3 scripts/fable_zenodo/search_repositories.py
python3 scripts/fable_zenodo/filter_candidates.py
python3 scripts/fable_zenodo/download_selected.py
python3 scripts/fable_zenodo/analyze_luad.py
python3 scripts/fable_zenodo/analyze_ici_imc.py
```

### 7. Caveats

- Relevance ranking is title-based (descriptions were not fetched for all 217
  hits), so some borderline datasets may be missed; the two chosen datasets were
  manually confirmed.
- CP10K normalization uses the provided `nCount_RNA`; no additional QC filtering
  was re-applied beyond the authors' processing.
- The IMC pre/post comparison is a pooled cell-level mean (not patient-level
  mixed model), so treat magnitudes as descriptive, not inferential.
- No raw/restricted files were downloaded; only open CC-BY processed matrices.

---

## 中文

### 1. 目标

在 Zenodo / figshare / OSF 上检索与**肺癌免疫检查点抑制剂（ICI）治疗**以及上皮/
抗体偶联药物（ADC）靶点 **TACSTD2（TROP2）** 和 **CLDN4** 相关的**开放、已处理矩阵**。
下载 2 GB 以内的开放处理文件，跳过受限 Zenodo 记录，核验 DOI，并在可行范围内分析。

### 2. 检索

`scripts/fable_zenodo/search_repositories.py` 用 8 组关键词查询了 Zenodo REST API、
figshare v2 搜索与 OSF v2 节点，共得到 **217 条去重记录**（Zenodo 119、figshare 98；
OSF 标题无匹配）。经基于标题的相关性打分（`filter_candidates.py`）后，有 40 条相关
记录带有可下载的矩阵类文件。详见 `results/fable_zenodo/search_results.json`、
`shortlist.json` 及 `notes/fable_zenodo/search_summary.md`。

### 3. DOI 核验

所有候选 DOI 均通过 `https://doi.org/<doi>` 解析；最终使用的两个数据集及另外两个
切题候选均返回 HTTP 200。机器可读记录见 `results/fable_zenodo/doi_verification.json`。

### 4. 下载（开放、CC-BY-4.0、< 2 GB）

含大小与 MD5 的清单：`results/fable_zenodo/download_manifest.json`。

| Zenodo DOI | 数据集 | 下载文件 | 类型 |
|------------|--------|----------|------|
| [10.5281/zenodo.10731914](https://doi.org/10.5281/zenodo.10731914) | 肺腺癌（LUAD）患者与小鼠 scRNA-seq | `cell_matrix_sce_8LUADs.csv.gz`（58 MB）、`cell_annotation_sce_8LUADs.csv`（9.6 MB） | scRNA-seq 计数，28,457 个人类细胞 |
| [10.5281/zenodo.8041882](https://doi.org/10.5281/zenodo.8041882) | NSCLC 对 **anti-PD1** 治疗的获得性耐药 | `sce-clustering-allsamples-original_final.RDS`（160 MB）、`AK_List.csv` | 成像质谱流式（IMC），41 标记 × 99,659 细胞 |

**受限 Zenodo 记录已跳过**（共 4 条切题记录，列于检索小结与 `doi_verification.json`）。

### 5. 分析

#### 5a. 人类 LUAD scRNA-seq 中的 TACSTD2 与 CLDN4（Zenodo 10731914）

脚本：`scripts/fable_zenodo/analyze_luad.py`。矩阵为「基因 × 细胞」原始 UMI；我们用
`nCount_RNA` 将每个细胞归一化为 CP10K，并按 `Cell_type.unimodel` 注释汇总（28,457 个
细胞、7 类群体）。

两个基因都**高度且特异地局限于上皮（肿瘤）细胞**：

| 基因 | 上皮阳性率 | 上皮平均 CP10K | 次高群体 |
|------|-----------|----------------|----------|
| **TACSTD2** | 73.5% | 10.18 | 髓系 5.3% / 0.18 |
| **CLDN4** | 80.0% | 6.09 | 成纤维 2.5% / 0.08 |

免疫细胞（T/NK、B、髓系、肥大细胞）与基质细胞（内皮、成纤维）基本为阴性（阳性率
< 约 5%，平均 CP10K < 0.2）。这与已知生物学一致：TACSTD2/TROP2 与 CLDN4 是上皮表面
抗原——TROP2 正是 ADC 药物 sacituzumab govitecan 的靶点，CLDN4 也是新兴的实体瘤靶点。
在 ICI 语境下，这种肿瘤细胞特异表达使它们成为对免疫细胞影响较小的理想靶点。

表格：`results/fable_zenodo/luad_expression_by_celltype.csv`、
`luad_expression_global.csv`；图：`luad_tacstd2_cldn4_dotplot.png`。

#### 5b. Anti-PD1 NSCLC 微环境（Zenodo 8041882）

脚本：`scripts/fable_zenodo/analyze_ici_imc.py`。该 RDS 为 S4 `SingleCellExperiment`；
环境未装 R，故用 Python 的 `rdata` 包直接解析。它实际是一个**成像质谱流式（IMC）**对象
——**41 个蛋白标记 × 99,659 个单细胞**（7 名患者），带有 `TherapyStatus`（治疗前/后
anti-PD1）、`ImmuneStatus`（T 细胞高/低）与聚类标签。**TACSTD2/CLDN4 属于转录组标记，
不在该靶向蛋白面板中，因此无法在此测量**；该数据集刻画的是 ICI 治疗后的肿瘤免疫微环境。

比较 **anti-PD1 治疗后 vs 治疗前**的平均 arcsinh 表达（治疗前 46,638 / 治疗后 53,021 细胞）：

- **治疗后升高**：CD15（+0.48）、CD45RO 记忆 T（+0.41）、CD134/OX40（+0.32）、
  CD73（+0.25）、碳酸酐酶 IX / 缺氧（+0.23）、颗粒酶 B（+0.13）。
- **治疗后降低**：CD163（−0.15）、CD14（−0.11）、Siglec-1/CD169（−0.11）——即 M2 样/
  抑制性髓系标记下降。

方向与原文主题一致（获得性 anti-PD1 耐药中的免疫抑制性重塑，伴 CD73/OX40/记忆 T 变化）。
表格：`results/fable_zenodo/ici_imc_marker_pre_post.csv`、`ici_imc_panel.csv`、
`ici_imc_summary.json`；图：`ici_imc_pre_post_markers.png`。

### 6. 复现

```bash
pip install pandas scipy h5py matplotlib rdata
python3 scripts/fable_zenodo/search_repositories.py
python3 scripts/fable_zenodo/filter_candidates.py
python3 scripts/fable_zenodo/download_selected.py
python3 scripts/fable_zenodo/analyze_luad.py
python3 scripts/fable_zenodo/analyze_ici_imc.py
```

### 7. 注意事项

- 相关性排序基于标题（未对全部 217 条抓取描述），可能漏掉个别边缘数据集；最终选用的
  两个数据集已人工确认。
- CP10K 归一化使用作者提供的 `nCount_RNA`；除作者已有处理外未再做额外 QC 过滤。
- IMC 治疗前/后比较为汇总的细胞级均值（非患者级混合模型），量级仅供描述性参考。
- 未下载任何原始/受限文件，仅使用开放的 CC-BY 处理矩阵。
