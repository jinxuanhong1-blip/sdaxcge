# Lung ICI scRNA-seq Methods Playbook (Seurat v5 / Scanpy / OSCA)
# 肺癌免疫检查点抑制剂 (ICI) 单细胞 RNA 测序方法手册（Seurat v5 / Scanpy / OSCA）

> **Scope / 范围**: This is a **methods-only** playbook. It documents *how* to run a 2025–2026 single-cell RNA-seq (scRNA-seq) workflow for lung cancer immune-checkpoint-inhibitor (ICI) cohorts. It contains **no datasets, no results, and no invented numbers, gene lists, effect sizes, or citations of specific findings**. Every place where you must supply cohort-specific values (sample paths, batch keys, gene signatures, thresholds) is marked `# TODO`.
>
> **范围**：这是一份**仅方法学 (methods-only)** 的手册，记录如何在 2025–2026 年针对肺癌 ICI 队列运行单细胞 RNA 测序流程。手册**不包含任何数据集、结果，也不编造任何数字、基因列表、效应量或对具体研究结论的引用**。凡需要你根据自己队列填写的位置（样本路径、批次键、基因签名、阈值）均以 `# TODO` 标注。

---

## 0. How to read this document / 如何阅读本文档

**EN** — Each section states the *decision* first (what to choose and why), then the *mechanics* (the exact API calls). Runnable, parameterized templates live next to this file:

- `templates/seurat_v5_template.R` — Seurat v5 end-to-end (layers, BPCells, sketching, `IntegrateLayers`, Azimuth).
- `templates/osca_bioc_template.R` — OSCA / Bioconductor end-to-end (`SingleCellExperiment`, `scran`/`scater`/`batchelor`, `SingleR`, `AUCell`).
- `templates/scanpy_template.py` — Scanpy end-to-end (AnnData layers, HVG, PCA/neighbors/UMAP/Leiden, Harmony, `score_genes`).
- `templates/celltypist_scvi_template.py` — CellTypist annotation + scVI/scANVI integration.
- `signatures/tacstd2_cldn4_modules.md` + `signatures/modules.template.json` — how to *build and validate* module gene sets (you supply the genes).

**中文** — 每节先给出**决策**（选什么、为什么），再给出**操作**（具体 API 调用）。可运行、可参数化的模板与本文件放在同一目录（见上）。

---

## 1. Data model: layers & assays / 数据模型：layers 与 assays

### 1.1 Seurat v5 — the `Assay5` object and layers

**EN** — Seurat v5 replaces the old `@slots` (`counts`/`data`/`scale.data`) with named **layers** inside an `Assay5` object. Key ideas:

- A single assay can hold many layers, e.g. `counts`, `data` (log-normalized), `scale.data`.
- When you `split` an object by a metadata column (e.g. sample/patient), Seurat creates *per-group layers* named `counts.<group>`, `data.<group>`. This is the mechanism that powers layer-aware integration (Section 5).
- `JoinLayers()` collapses per-group layers back into single `counts`/`data` layers — **required before differential expression / pseudobulk**.
- Multi-modal data (CITE-seq ADT, HTO) live as **separate assays** (`obj[["ADT"]]`), each with its own layers.

```r
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample_id)   # per-sample layers
Layers(obj[["RNA"]])                                      # inspect layers
obj <- JoinLayers(obj)                                    # collapse before DE
```

**中文** — Seurat v5 用命名 **layers** 取代了旧版的 `@slots`（`counts`/`data`/`scale.data`）。要点：

- 一个 assay 可容纳多个 layer，如 `counts`、`data`（log 归一化）、`scale.data`。
- 当你按某个 metadata 列（如 sample/patient）`split` 对象时，Seurat 会创建**按组分开的 layer**，命名为 `counts.<组名>`、`data.<组名>`。这是支撑“分层感知整合”（第 5 节）的机制。
- `JoinLayers()` 会把分组 layer 合并回单一的 `counts`/`data`——**在差异表达 / pseudobulk 之前必须执行**。
- 多模态数据（CITE-seq 的 ADT、HTO）作为**独立 assay** 存放（`obj[["ADT"]]`），各自拥有自己的 layers。

### 1.2 Scanpy / AnnData — `.X`, `.layers`, `.obs`, `.var`, `.obsm`

**EN** — In AnnData, the analogous structure is:

| AnnData field | Holds | Seurat v5 analogue |
|---|---|---|
| `adata.X` | the "active" matrix (convention: log-normalized) | `data` layer |
| `adata.layers["counts"]` | raw counts (keep a copy!) | `counts` layer |
| `adata.obs` | per-cell metadata | `obj@meta.data` |
| `adata.var` | per-gene metadata (HVG flags, etc.) | feature meta |
| `adata.obsm["X_pca"]`, `["X_umap"]`, `["X_scVI"]` | embeddings | `Reductions` |
| `adata.raw` | frozen pre-HVG state (optional) | — |

Convention that avoids the most common bugs: **always keep raw counts in `adata.layers["counts"]`** and treat `adata.X` as log-normalized after preprocessing. Tools such as scVI read counts from a named layer.

**中文** — AnnData 的对应结构见上表。避免常见错误的约定：**始终把原始 counts 保存在 `adata.layers["counts"]`**，并在预处理后把 `adata.X` 作为 log 归一化矩阵；scVI 等工具会从指定 layer 读取 counts。

### 1.3 OSCA / Bioconductor — `SingleCellExperiment`

**EN** — OSCA uses `SingleCellExperiment` (SCE). Matrices live in `assays(sce)` (`counts`, `logcounts`); cell metadata in `colData`; gene metadata in `rowData`; embeddings in `reducedDims(sce)` (`PCA`, `UMAP`, `corrected`). Alternative experiments (ADT) sit in `altExp(sce, "ADT")`. This is the interchange format for `scran`, `scater`, `batchelor`, `SingleR`, and `AUCell`.

**中文** — OSCA 使用 `SingleCellExperiment`（SCE）。矩阵在 `assays(sce)`（`counts`、`logcounts`）；细胞元数据在 `colData`；基因元数据在 `rowData`；降维在 `reducedDims(sce)`；多模态在 `altExp(sce, "ADT")`。它是 `scran`/`scater`/`batchelor`/`SingleR`/`AUCell` 的通用交换格式。

### 1.4 Interoperability / 互操作

**EN** — Convert deliberately, and always carry raw counts across the boundary:
- Seurat ⇄ SCE: `Seurat::as.SingleCellExperiment()` / `Seurat::as.Seurat()`.
- AnnData ⇄ Seurat/SCE: `sceasy`, `zellkonverter::readH5AD()/writeH5AD()`, or `anndata`/`schard` in R. Prefer round-tripping the **counts** layer; recompute normalization/HVG/PCA in the destination tool rather than trusting a converted `scale.data`/`X`.

**中文** — 转换要有意为之，并始终把原始 counts 带过边界：Seurat ⇄ SCE 用 `as.SingleCellExperiment()`/`as.Seurat()`；AnnData ⇄ R 用 `zellkonverter`/`sceasy`。优先传递 **counts** layer，并在目标工具中重新计算归一化/HVG/PCA，而不要信任转换后的 `scale.data`/`X`。

---

## 2. QC & preprocessing / 质控与预处理

**EN** — Order of operations (identical logic across the three ecosystems):

1. **Ambient RNA / empty drops**: `DropletUtils::emptyDrops` (OSCA) on raw unfiltered matrices; optionally `SoupX` / `CellBender` for ambient correction *upstream*.
2. **Per-cell QC metrics**: library size, detected genes, `% mitochondrial` (`^MT-`), `% ribosomal` (`^RP[SL]`), `% hemoglobin` (`^HB[^(P)]`). For lung/tumor samples, mitochondrial fraction is frequently elevated in stressed tumor/epithelial cells — **set the MT cutoff per-cohort from the observed distribution, do not hard-code a universal 10%/20%**.
3. **Doublets**: `scDblFinder` (OSCA/R) or `scrublet`/`sc.pp.scrublet` (Scanpy). Run **per sample/lane**, never on the merged object.
4. **Filter**, then **normalize**: log-normalization (`NormalizeData` / `sc.pp.normalize_total`+`log1p` / `scran::computeSumFactors`+`logNormCounts`) or `SCTransform`. Deconvolution size factors (`scran`) are the OSCA default and handle composition differences well.

> **No invented thresholds**: numeric QC cutoffs below are `# TODO`. Derive them from your own per-sample distributions (e.g. MAD-based outlier detection via `scater::isOutlier`).

**中文** — 操作顺序（三套生态逻辑一致）：

1. **环境 RNA / 空液滴**：在未过滤的原始矩阵上用 `DropletUtils::emptyDrops`（OSCA）；上游可选 `SoupX`/`CellBender` 做环境 RNA 校正。
2. **单细胞 QC 指标**：文库大小、检出基因数、线粒体比例（`^MT-`）、核糖体比例（`^RP[SL]`）、血红蛋白比例（`^HB[^(P)]`）。肺/肿瘤样本中，应激的肿瘤/上皮细胞线粒体比例常偏高——**MT 阈值应按队列从实际分布确定，不要写死通用的 10%/20%**。
3. **双细胞**：`scDblFinder`（R）或 `scrublet`（Scanpy）。**按样本/lane 逐个运行**，绝不要在合并后的对象上运行。
4. **过滤**后**归一化**：log 归一化或 `SCTransform`；`scran` 的反卷积 size factor 是 OSCA 默认，能较好处理组成差异。

> **不编造阈值**：下方数值 QC 阈值均为 `# TODO`，请从你自己的分布推导（如 `scater::isOutlier` 的 MAD 离群检测）。

---

## 3. Scaling to large cohorts: BPCells & sketching / 大队列扩展：BPCells 与 sketching

### 3.1 BPCells (on-disk) / 磁盘存储

**EN** — Lung ICI atlases routinely reach 10⁶–10⁷ cells. **BPCells** (Greenleaf Lab) stores the counts matrix on disk in a bit-packed, streaming C++ format, so a `Seurat` object can reference a matrix far larger than RAM. Workflow:

- Ingest each sample's counts once and write a BPCells directory: `BPCells::write_matrix_dir()`.
- Wrap it as a layer: the `counts` layer of the `RNA` assay is a BPCells matrix pointer, not an in-memory `dgCMatrix`.
- Normalization and feature statistics stream over disk; only the **sketch** (Section 3.2) is loaded into memory.

```r
library(BPCells)
mat <- open_matrix_10x_hdf5(path = "sample.h5")   # or open_matrix_dir()
write_matrix_dir(mat = mat, dir = "bpcells/sample")  # one-time conversion
mat.disk <- open_matrix_dir("bpcells/sample")
obj <- CreateSeuratObject(counts = mat.disk)          # counts stay on disk
```

Scanpy analogues for out-of-core work: **`AnnData` in `backed="r"` mode**, **Dask-backed AnnData** (`anndata` experimental / `rapids-singlecell`), or TileDB-SOMA (`cellxgene-census` / `tiledbsoma`). Choose one and document it — do not mix.

**中文** — 肺癌 ICI 图谱常达 10⁶–10⁷ 细胞。**BPCells**（Greenleaf 实验室）以位打包的流式 C++ 格式把 counts 矩阵存于磁盘，使 `Seurat` 对象可引用远超内存的矩阵。流程：用 `write_matrix_dir()` 一次性转换每个样本；把磁盘矩阵作为 `RNA` assay 的 `counts` layer；归一化与特征统计以流式方式在磁盘上完成，只有 **sketch**（3.2 节）载入内存。Scanpy 的对应方案：`AnnData` 的 `backed="r"` 模式、Dask 支持的 AnnData（`rapids-singlecell`）或 TileDB-SOMA（`tiledbsoma`/`cellxgene-census`）。选定一种并记录，不要混用。

### 3.2 Sketching (leverage-score subsampling) / Sketch（杠杆分数子采样）

**EN** — **Sketching** selects a small, information-rich subset ("atoms") of cells for in-memory analysis, then projects results back to all cells. Unlike uniform downsampling, **leverage-score** sampling oversamples rare populations, which is essential in lung TME where rare states (e.g. tissue-resident memory T cells, specific DC subsets, rare epithelial/tumor states) drive ICI biology.

Seurat v5:

```r
obj <- SketchData(obj, ncells = 50000, method = "LeverageScore", sketched.assay = "sketch")
DefaultAssay(obj) <- "sketch"        # fast, in-memory analysis here
# ... normalize/HVG/PCA/cluster on the sketch ...
obj <- ProjectData(                  # extend labels/embeddings to the full set
  obj, assay = "RNA", sketched.assay = "sketch",
  sketched.reduction = "pca", full.reduction = "pca.full",
  umap.model = "umap", dims = 1:50,
  refdata = list(cluster_full = "seurat_clusters"))
```

Guidance: `ncells` per layer is a memory/resolution trade-off (`# TODO`, e.g. 5k per sample for many samples, or 50k total). Cluster/integrate on the sketch, then `ProjectData`/`ProjectIntegration` to the full dataset; switch `DefaultAssay` back to `RNA` for the on-disk full data. For Scanpy, the analogue is `sc.pp.subsample`/geometric sketching (`geosketch`) followed by kNN label transfer (`sc.tl.ingest` or a trained scANVI model) — geometric sketching similarly preserves rare cells.

**中文** — **Sketch** 选取一小部分信息丰富的细胞（"atoms"）在内存中分析，再把结果投影回全部细胞。与均匀降采样不同，**杠杆分数**采样会过采样稀有群体——这在肺 TME 中至关重要，因为稀有状态（组织驻留记忆 T 细胞、特定 DC 亚群、稀有上皮/肿瘤状态）驱动 ICI 生物学。用 `SketchData()` 生成 sketch，把 `DefaultAssay` 切到 `sketch` 做快速内存分析，再用 `ProjectData()`/`ProjectIntegration()` 投影回全集，之后把 `DefaultAssay` 切回磁盘上的 `RNA`。`ncells` 是内存/分辨率的权衡（`# TODO`）。Scanpy 对应：`geosketch` 几何 sketch + kNN 标签迁移（`sc.tl.ingest` 或训练好的 scANVI 模型）。

---

## 4. Dimensionality reduction & clustering / 降维与聚类

**EN** — Standard path in all three: HVG selection → scale/PCA → neighbor graph → UMAP + graph clustering.

- **HVG**: `FindVariableFeatures` (Seurat, `vst`), `sc.pp.highly_variable_genes` (Scanpy; use `flavor="seurat_v3"` on counts, or `"pearson_residuals"`), `scran::modelGeneVar` + `getTopHVGs` (OSCA). Consider **batch-aware HVGs** (`batch_key=` in Scanpy; `block=` in `modelGeneVar`) so batch genes don't dominate.
- **PCA / denoise**: `RunPCA`; `sc.pp.pca`; `scran::denoisePCA` (chooses PCs by retained biological variance).
- **Neighbors + clustering**: `FindNeighbors`+`FindClusters` (Louvain/Leiden); `sc.pp.neighbors`+`sc.tl.leiden`; `scran::buildSNNGraph`+`igraph`/`clusterCells`.
- **Leiden note (2025)**: in current Scanpy, set `sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False)` — the default flavor is changing and the `igraph` backend is faster and reproducible. Sweep `resolution` (`# TODO`) rather than trusting one value.

**中文** — 三套流程一致：HVG 选择 → scale/PCA → 近邻图 → UMAP + 图聚类。HVG 建议**批次感知**（Scanpy 的 `batch_key=`、`modelGeneVar` 的 `block=`），避免批次基因主导。PCA 去噪可用 `scran::denoisePCA`。聚类用 Louvain/Leiden。**Leiden 提示（2025）**：当前 Scanpy 建议显式设置 `flavor="igraph", n_iterations=2, directed=False`，更快且可复现。对 `resolution` 做扫描（`# TODO`），不要只信一个值。

---

## 5. Batch integration APIs / 批次整合 API

**EN** — In lung ICI cohorts the dominant batch axis is usually **patient/sample** (and sometimes chemistry/site/timepoint pre/post-treatment). Integrate to remove technical batch **without erasing the treatment/response signal you care about** — so *never* integrate over the response variable, and always compare an unintegrated embedding side-by-side. Assess with kBET / iLISI / cLISI / `scib-metrics` (`# TODO`).

### 5.1 Seurat v5 `IntegrateLayers` (verified API)

**EN** — One function, five methods, operating on the split layers from Section 1.1. It returns a **new dimensional reduction**; the expression matrix is not modified.

```r
obj <- IntegrateLayers(obj, method = CCAIntegration,   orig.reduction = "pca", new.reduction = "integrated.cca")
obj <- IntegrateLayers(obj, method = RPCAIntegration,  orig.reduction = "pca", new.reduction = "integrated.rpca")
obj <- IntegrateLayers(obj, method = HarmonyIntegration,orig.reduction = "pca", new.reduction = "harmony")
obj <- IntegrateLayers(obj, method = FastMNNIntegration,orig.reduction = "pca", new.reduction = "integrated.mnn")
obj <- IntegrateLayers(obj, method = scVIIntegration,  new.reduction = "integrated.scvi",
                       conda_env = "/path/to/scvi-env")   # TODO: env with scvi-tools
# then cluster/UMAP on the chosen reduction, e.g.:
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- JoinLayers(obj)   # before DE / pseudobulk
```

Method picking: **RPCA** for speed/robustness on large or low-overlap datasets; **CCA** for maximal alignment of shared types; **Harmony** for many samples and fast iteration; **FastMNN** as a well-characterized graph/MNN option; **scVI** when you want a probabilistic latent space reusable for scANVI label transfer. For SCT-normalized data pass `normalization.method = "SCT"`.

**中文** — 一个函数、五种方法，作用于第 1.1 节 split 出来的 layer，返回**新的降维**，不修改表达矩阵。方法选择：大数据/低重叠用 **RPCA**（快而稳）；共享类型最大对齐用 **CCA**；样本多、需快速迭代用 **Harmony**；成熟的图/MNN 方案用 **FastMNN**；需要可复用于 scANVI 标签迁移的概率隐空间用 **scVI**。SCT 数据请传 `normalization.method = "SCT"`。整合后在所选 reduction 上聚类/UMAP，DE/pseudobulk 前先 `JoinLayers()`。**关键原则**：整合是为去除技术批次，绝不要在响应变量（response/treatment）上做整合，并始终与未整合的嵌入并排对比（用 `scib-metrics` 等评估，`# TODO`）。

### 5.2 Scanpy integration / Scanpy 整合

| Method | Call | Output |
|---|---|---|
| Harmony | `sce.pp.harmony_integrate(adata, key="sample")` | `adata.obsm["X_pca_harmony"]` |
| BBKNN | `sce.pp.bbknn(adata, batch_key="sample")` | modifies neighbor graph |
| Scanorama | `sce.pp.scanorama_integrate(adata, key="sample")` | `adata.obsm["X_scanorama"]` |
| scVI/scANVI | `scvi-tools` (see `templates/celltypist_scvi_template.py`) | `adata.obsm["X_scVI"]` |

(`sce` = `scanpy.external`.) Then `sc.pp.neighbors(adata, use_rep="X_pca_harmony")` → UMAP/Leiden on the corrected embedding.

**中文** — Scanpy 整合选项见上表（`sce` 即 `scanpy.external`）。之后用 `sc.pp.neighbors(use_rep=...)` 在校正后的嵌入上做 UMAP/Leiden。

### 5.3 OSCA / batchelor / OSCA 整合

**EN** — `batchelor::fastMNN()` (returns a `corrected` reducedDim), `batchelor::rescaleBatches()` / `regressBatches()` for lightweight correction, or `harmony::RunHarmony()` on the PCA. Then `buildSNNGraph` on the corrected embedding.

**中文** — 用 `batchelor::fastMNN()`（产出 `corrected` 降维）、`rescaleBatches()`/`regressBatches()` 轻量校正，或对 PCA 用 `harmony::RunHarmony()`；随后在校正嵌入上 `buildSNNGraph`。

---

## 6. Cell-type annotation / 细胞类型注释

**EN** — Use a *layered* strategy: (a) automated reference/classifier labels, (b) marker-based curation, (c) manual review of ambiguous clusters. Never rely on a single automated call for a publication-grade lung atlas.

### 6.1 CellTypist (Python) — verified API

**EN** — Logistic-regression classifiers. Input must be **log1p-normalized to 10,000 counts per cell** (CellTypist's expected scale). Relevant built-in models: `Immune_All_Low.pkl` / `Immune_All_High.pkl` (immune sub-populations) and `Human_Lung_Atlas.pkl` (lung). `majority_voting=True` refines per-cell labels over clusters.

```python
import celltypist
from celltypist import models
models.download_models(model=["Immune_All_Low.pkl", "Human_Lung_Atlas.pkl"])  # TODO pick models
# adata.X must be log1p CPM-10k; keep raw in adata.layers["counts"]
pred = celltypist.annotate(adata, model="Immune_All_Low.pkl", majority_voting=True)
adata = pred.to_adata()   # adds predicted_labels, majority_voting, conf_score
```

**中文** — 逻辑回归分类器。输入必须是**每细胞归一化到 10,000 counts 后再 log1p** 的数据（CellTypist 期望的尺度）。相关内置模型：`Immune_All_Low.pkl`/`Immune_All_High.pkl`（免疫亚群）与 `Human_Lung_Atlas.pkl`（肺）。`majority_voting=True` 会在聚类层面细化每细胞标签。用 `celltypist.annotate` 后 `to_adata()` 得到 `predicted_labels`、`majority_voting`、`conf_score`。

### 6.2 Azimuth (R) — verified reference

**EN** — Reference-mapping onto curated atlases. The lung reference is **`lungref` (v2.0.0)**, which corresponds to the Human Lung Cell Atlas (HLCA) v2. `RunAzimuth` starts from raw counts and adds predicted labels at multiple annotation levels plus a mapping/prediction score.

```r
library(Azimuth); library(SeuratData)
InstallData("lungref")                    # Azimuth Reference: lung, v2.0.0 (HLCA v2)
obj <- RunAzimuth(obj, reference = "lungref")
# adds predicted.ann_level_* , predicted.ann_finest_level, mapping.score, prediction.score.*
```

Treat `mapping.score`/`prediction.score` as confidence; low-confidence cells (`# TODO` threshold) should fall through to manual review.

**中文** — 参考映射到精选图谱。肺参考为 **`lungref`（v2.0.0）**，对应人类肺细胞图谱 HLCA v2。`RunAzimuth` 从原始 counts 开始，添加多层级预测标签及映射/预测分数。把 `mapping.score`/`prediction.score` 视为置信度，低置信细胞（阈值 `# TODO`）应进入人工复核。

### 6.3 SingleR (OSCA) / SingleR（OSCA）

**EN** — `SingleR::SingleR()` labels cells against a reference SCE by correlation of marker genes; pairs naturally with `celldex` bulk references or a custom lung reference SCE. Good cross-check against CellTypist/Azimuth.

**中文** — `SingleR::SingleR()` 通过标记基因相关性对参考 SCE 进行标注；可搭配 `celldex` bulk 参考或自建肺参考 SCE，是对 CellTypist/Azimuth 的良好交叉验证。

### 6.4 Marker-based confirmation / 基于标记的确认

**EN** — Confirm major lung TME compartments with canonical, literature-established markers (supply the list yourself, `# TODO`): T/NK, B/plasma, myeloid (macrophage/monocyte/DC), mast, epithelial (incl. AT1/AT2/club/basal and malignant), endothelial, fibroblast/smooth-muscle. Use `FindAllMarkers` / `sc.tl.rank_genes_groups` / `scran::scoreMarkers`, and dotplots/feature plots to reconcile automated labels with markers.

**中文** — 用文献确立的经典标记确认肺 TME 主要区室（标记列表由你提供，`# TODO`）：T/NK、B/浆细胞、髓系（巨噬/单核/DC）、肥大细胞、上皮（含 AT1/AT2/club/basal 及恶性）、内皮、成纤维/平滑肌。用 `FindAllMarkers`/`rank_genes_groups`/`scoreMarkers` 及点图/特征图，将自动标签与标记对齐。

---

## 7. Scoring the TACSTD2 / CLDN4 modules / TACSTD2 与 CLDN4 模块打分

> **Why these genes / 为什么关注这两个基因**: **TACSTD2** (TROP2) and **CLDN4** (Claudin-4) are epithelial cell-surface molecules and clinically relevant antibody-drug-conjugate (ADC) targets in lung/NSCLC. In the HLCA v2 marker panels, **TACSTD2** is a canonical epithelial/airway-epithelium marker; **CLDN4** is a tight-junction epithelial marker. Module scoring quantifies coordinated expression of a curated gene program per cell, which is more robust than a single gene. **This playbook does not provide a signature gene list — you must curate it** (Section 7.4).

### 7.1 What "module score" means / “模块打分”的含义

**EN** — A module score summarizes the average expression of a gene set per cell, corrected for a background of control genes matched by expression level (so highly expressed housekeeping genes don't inflate the score). The two mainstream algorithms:
- **Seurat `AddModuleScore`** / **Scanpy `sc.tl.score_genes`**: mean of the module genes minus mean of expression-binned random control genes. (Same idea in both; identical enough to compare.)
- **UCell (R/Python)** / **AUCell (R)**: rank-based (Mann-Whitney U / AUC) scoring — insensitive to normalization/depth and well-suited to sparse data; good for cross-sample robustness.

**中文** — 模块分数概括每个细胞中一组基因的平均表达，并以按表达量匹配的对照基因作背景校正（避免高表达管家基因抬高分数）。两类主流算法：`AddModuleScore`/`sc.tl.score_genes`（模块均值减去按表达分箱的随机对照均值）；UCell/AUCell（基于秩的 AUC 打分，对归一化/测序深度不敏感，适合稀疏数据与跨样本稳健性）。

### 7.2 Seurat v5 (R)

```r
# modules is a named list you curate; see signatures/tacstd2_cldn4_modules.md
modules <- list(
  TACSTD2_module = tacstd2_genes,   # TODO: curated character vector, includes "TACSTD2"
  CLDN4_module   = cldn4_genes      # TODO: curated character vector, includes "CLDN4"
)
DefaultAssay(obj) <- "RNA"; obj <- JoinLayers(obj)   # score on joined log-normalized data
obj <- AddModuleScore(obj, features = modules, name = names(modules), seed = 42, ctrl = 100)
# -> meta.data columns TACSTD2_module1, CLDN4_module2  (Seurat appends an index)

# Rank-based alternative (recommended for robustness):
# library(UCell); obj <- AddModuleScore_UCell(obj, features = modules)
```

### 7.3 Scanpy (Python)

```python
# gene lists you curate; keep only genes present in adata.var_names
tacstd2_genes = [g for g in TACSTD2_MODULE if g in adata.var_names]  # TODO curate; includes "TACSTD2"
cldn4_genes   = [g for g in CLDN4_MODULE   if g in adata.var_names]  # TODO curate; includes "CLDN4"
# adata.X should be log1p-normalized
sc.tl.score_genes(adata, tacstd2_genes, score_name="TACSTD2_module", ctrl_size=100, random_state=0)
sc.tl.score_genes(adata, cldn4_genes,   score_name="CLDN4_module",   ctrl_size=100, random_state=0)
# Rank-based alternative:
# import ucell  # or scanpy-compatible UCell port; AUCell available in R
```

### 7.4 Building a defensible module (no invented data) / 构建可辩护的模块（不编造数据）

**EN** — Because fabricated signatures are not acceptable, build the gene set by one of these transparent routes and **record the provenance** in `signatures/tacstd2_cldn4_modules.md`:
1. **Anchor + curated program**: start from the anchor gene (`TACSTD2` or `CLDN4`) and add genes from a **named, citable database** (e.g. MSigDB Hallmark/Reactome epithelial or tight-junction gene sets, or the published HLCA epithelial marker panels). Cite the exact set/version.
2. **Data-driven co-expression**: compute genes co-expressed with the anchor **in your own data** (e.g. Spearman correlation on the epithelial subset, or a `hdWGCNA`/`Hotspot` module containing the anchor). This produces a cohort-specific module with a reproducible recipe rather than an assumed list.
3. **Validate**: check the module localizes to the expected compartment (epithelial/malignant), correlates with the anchor gene, is stable across samples (score per-sample), and is not merely tracking total counts (correlate score vs. `nCount`; regress out if needed).

**中文** — 由于不接受编造签名，请用以下**可追溯**方式之一构建基因集，并在 `signatures/tacstd2_cldn4_modules.md` 中**记录来源**：
1. **锚基因 + 精选程序**：从锚基因（`TACSTD2` 或 `CLDN4`）出发，补充来自**具名可引用数据库**的基因（如 MSigDB Hallmark/Reactome 的上皮或紧密连接基因集、已发表的 HLCA 上皮标记面板），并注明确切集合与版本。
2. **数据驱动共表达**：在**你自己的数据**中计算与锚基因共表达的基因（如在上皮子集上做 Spearman 相关，或用 `hdWGCNA`/`Hotspot` 得到包含锚基因的模块），得到有可复现配方的队列特异模块，而非假设的列表。
3. **验证**：确认模块定位到预期区室（上皮/恶性）、与锚基因相关、跨样本稳定（按样本打分），且并非只反映总 counts（将分数与 `nCount` 相关，必要时回归剔除）。

### 7.5 Downstream / 下游

**EN** — Once scored, compare module scores across annotated compartments and (if you have treatment metadata) across pre/post-ICI or responder/non-responder groups **using per-sample aggregation (pseudobulk) and appropriate mixed models** — do not run per-cell tests treating cells as independent replicates. That statistical design is out of scope here; this playbook stops at producing the per-cell scores.

**中文** — 打分后，在注释区室间、以及（若有治疗元数据）ICI 前后 / 响应与否分组间比较模块分数时，**应按样本聚合（pseudobulk）并用合适的混合模型**，不要把细胞当作独立重复做每细胞检验。该统计设计不在本手册范围内；本手册止步于生成每细胞分数。

---

## 8. Reproducibility & environment / 可复现性与环境

**EN** — Pin versions and record them with the run:
- R: Seurat v5.x, SeuratObject v5.x, BPCells, Signac (if ATAC), Azimuth + SeuratData (`lungref` 2.0.0), scran/scater/batchelor/SingleR/AUCell (Bioconductor release matched to your R version), UCell. Capture `sessionInfo()`.
- Python: scanpy ≥1.10, anndata, scvi-tools, celltypist ≥1.7, harmonypy, scanorama, bbknn, geosketch, scib-metrics, leidenalg/igraph. Capture `sc.logging.print_versions()` and freeze with a lockfile (conda/`uv`/`pip freeze`).
- Set all random seeds (`SketchData(seed=)`, `random_state=`, `seed.use=`). Persist intermediate objects (`.rds`/`.h5ad`/BPCells dirs). Keep raw counts immutable.

**中文** — 固定并随运行记录版本：R 端（Seurat v5.x、BPCells、Azimuth+SeuratData 的 `lungref` 2.0.0、scran/scater/batchelor/SingleR/AUCell、UCell，记录 `sessionInfo()`）；Python 端（scanpy ≥1.10、scvi-tools、celltypist ≥1.7、harmonypy、scanorama、bbknn、geosketch、scib-metrics、leidenalg/igraph，记录 `print_versions()` 并用 lockfile 冻结）。设定所有随机种子；持久化中间对象；保持原始 counts 不可变。

---

## 9. End-to-end order of operations / 端到端操作顺序

```
ingest per sample ─▶ QC (empty drops, per-cell metrics, doublets) ─▶ normalize
   │                                                                    │
   ├── (large data) BPCells on-disk counts ──▶ SketchData ─────────────┤
   │                                                                    ▼
   └────────────────────────────────────────────────▶ HVG ─▶ PCA ─▶ IntegrateLayers / Harmony / scVI
                                                                        │
                          UMAP + Leiden on integrated reduction ◀───────┤
                                    │                                    │
        ProjectData/ProjectIntegration to full data (if sketched)       │
                                    │                                    ▼
             annotation: CellTypist / Azimuth(lungref) / SingleR + marker curation
                                    │
             JoinLayers ─▶ module scoring (TACSTD2 / CLDN4) ─▶ per-cell scores (STOP: analysis is out of scope)
```

**EN** — This is the canonical flow. Every dataset-specific value (paths, batch keys, thresholds, gene sets) is a `# TODO` you fill from your own cohort. The playbook and templates never assume or fabricate them.

**中文** — 这是标准流程。所有数据集特定的值（路径、批次键、阈值、基因集）都是需你按队列填写的 `# TODO`。本手册与模板绝不假设或编造这些值。
