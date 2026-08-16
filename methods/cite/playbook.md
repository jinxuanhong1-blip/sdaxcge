# CITE-seq Methods Playbook — GSE154826-class (RNA + immune ADT panel)

> **METHODS ONLY.** This document is a reproducible analysis protocol. It contains
> no results, no numbers, and no biological findings — only the methodology,
> parameters, and decision rules for analyzing a CITE-seq dataset of the
> GSE154826 class.
>
> Bilingual: **English** first, **中文** second (mirrored content). See the
> Chinese section starting at [中文版](#中文版-cite-seq-方法手册).

---

## 0. Dataset class and central premise

**Class.** GSE154826 (Leader et al., *Cancer Cell* 2021; Merad Lab) is a human
non-small-cell lung cancer (NSCLC) CITE-seq resource. A "GSE154826-class" dataset
has, per cell:

- **RNA** — 10x 3′ gene expression (GEX).
- **ADT** — antibody-derived tags for an **immune-focused surface-marker panel**
  (CD3/CD4/CD8/CD14/CD19/CD56/PD-1, etc.) plus **isotype controls**.
- **HTO** — cell-hashing tags (in GSE154826 raised against B2M/CD298) for sample
  multiplexing.
- Optionally **TCR** (VDJ).
- Tissues: tumor, non-involved lung (nLung), PBMC — and, in GSE154826, **spiked
  mouse cells** used as a technical control.

**Central analytical premise (read before touching the data).**
`TACSTD2` (which encodes the **TROP2** protein) and `CLDN4` (Claudin-4) are
**epithelial / carcinoma** RNA markers. The ADT panel in this dataset class is an
**immune panel**, so a **TROP2 antibody is usually absent** — i.e. TROP2 protein
is *not measured*. Consequences that drive every downstream choice:

1. **Epithelial identity is RNA-only.** You cannot orthogonally confirm
   `TACSTD2`/`CLDN4` at the protein level from the ADT panel.
2. **The two modalities are asymmetric.** ADT is informative for immune lineages;
   epithelial/tumor structure lives almost entirely in RNA.
3. **No method can rescue an unmeasured protein.** Neither WNN nor totalVI can
   "denoise into existence" a TROP2 ADT that was never in the panel. totalVI
   denoises *measured* proteins only.

Treat these three points as invariants throughout the pipeline.

---

## 1. Inputs and data model

Collect, per library/sample:

| Item | Purpose |
| --- | --- |
| RNA raw counts (genes × cells) | expression, HVG, clustering |
| ADT raw counts (proteins × cells), incl. isotype controls | protein signal + background estimation |
| HTO raw counts | demultiplexing |
| **Unfiltered / raw droplet matrices** (empty droplets kept) | ambient-RNA and ADT-background estimation |
| Sample / patient / tissue / batch metadata | integration, pseudobulk |
| Panel definition (exact ADT names) | modality-aware annotation and the TROP2 caveat |

**Recommended structure** — MuData with two modalities, raw counts preserved in
layers:

```python
import mudata as md, anndata as ad
mdata = md.MuData({"rna": rna_adata, "prot": prot_adata})
mdata["rna"].layers["counts"] = mdata["rna"].X.copy()   # keep integer counts
mdata["prot"].layers["counts"] = mdata["prot"].X.copy()  # totalVI needs raw ADT
```

For a Seurat/WNN track, use an object with `RNA` and `ADT` assays (plus `HTO`).

---

## 2. Demultiplexing (hashing)

1. Normalize HTO (CLR across cells).
2. Assign singlet / doublet / negative with `HTODemux` (Seurat),
   `hashsolo` (scanpy), or `GMM-Demux`.
3. Drop inter-sample **doublets** and **negatives**; carry the sample label.
4. If the library spiked mouse cells (GSE154826), remove them by aligning to a
   barnyard (human+mouse) reference or by a mouse-gene fraction threshold before
   any human-only analysis.

---

## 3. Quality control (per modality, then joint)

**RNA per cell:** `n_genes_by_counts`, `total_counts`, `pct_counts_mt`,
`pct_counts_ribo`, `pct_counts_hb`. Set thresholds per tissue (tumor cells
tolerate different mito/gene ranges than lymphocytes) — prefer **MAD-based**
outlier cutoffs over hard global numbers.

**ADT per cell:** total ADT counts, number of proteins detected, and
**isotype-control level** (high isotype → sticky/aggregating cell → drop).

**Joint rule:** a cell is retained only if it passes **both** RNA and ADT QC and
was classified a hashing singlet.

**Doublets (expression-based):** run `scDblFinder` / `DoubletFinder` /
`scrublet` in addition to HTO doublet calls; cross-modality doublets (e.g. an
immune + epithelial mixture) are a known source of spurious epithelial RNA in
immune barcodes — see §9.

---

## 4. Ambient / background correction (do not skip)

**RNA ambient contamination.** Lysed tumor epithelium releases `TACSTD2`/`CLDN4`
transcripts into the droplet suspension, so these genes can appear at low levels
in *bona fide* immune cells purely as ambient signal. Estimate and remove it with
**CellBender `remove-background`**, **SoupX**, or **DecontX**, using the
unfiltered droplet matrix. This step is a prerequisite for interpreting epithelial
RNA in non-epithelial clusters (§9).

**ADT background.** Antibody background is large and structured. Two accepted
routes:

- **DSB normalization** (Mulè et al. 2022): uses empty droplets + isotype controls
  to define and subtract background; feed DSB values into the WNN track.
- **totalVI's protein background mixture** (§6B): the model learns a per-protein
  background component from **raw** ADT counts — do **not** pre-DSB the counts you
  give totalVI.

---

## 5. Normalization

- **RNA:** `log1p` on CP10k, or SCTransform / analytic Pearson residuals for HVG
  selection (2,000–3,000 HVGs typical). Keep raw counts for totalVI and for
  pseudobulk DE.
- **ADT (WNN track):** DSB (preferred when empty droplets exist) or CLR; state the
  CLR margin (across cells is the common choice).
- **ADT (totalVI track):** none — totalVI consumes raw ADT counts.

---

## 6. Multimodal integration — run both tracks

Running WNN and totalVI in parallel and comparing them is the recommended QC: they
should agree on immune structure; disagreement flags modality-weight or
background-model problems.

### 6A. Weighted Nearest Neighbor (WNN, Seurat v4; Hao et al. 2021)

WNN learns a **cell-specific** weight for each modality and builds a joint graph.

```r
library(Seurat)
# assumes obj has SCT-processed RNA and DSB/CLR-normalized ADT
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, reduction.name = "pca")
DefaultAssay(obj) <- "ADT"
VariableFeatures(obj) <- rownames(obj[["ADT"]])          # small panel: use all
obj <- ScaleData(obj) |> RunPCA(reduction.name = "apca")

obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("pca", "apca"),
  dims.list      = list(1:30, 1:min(20, nrow(obj[["ADT"]]) - 1)),
  modality.weight.name = c("RNA.weight", "ADT.weight")
)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap")
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3)  # Leiden
```

**Modality-weight caveat.** Inspect `RNA.weight` / `ADT.weight`. Because there is
no TROP2 (or any epithelial) ADT, epithelial cells will be driven almost entirely
by **RNA weight**; WNN cannot use protein to sharpen epithelial boundaries. This is
expected, not a bug — but it means epithelial clustering quality is an RNA-only
property here.

### 6B. totalVI (scvi-tools; Gayoso et al. 2021)

totalVI is a joint generative model that (a) yields a shared latent space,
(b) denoises ADT by separating foreground from background, and (c) corrects batch.

```python
import scvi

scvi.model.TOTALVI.setup_mudata(
    mdata,
    rna_layer="counts",
    protein_layer="counts",
    batch_key="sample",                 # patient/sample; strong batch here
    modalities={"rna_layer": "rna", "protein_layer": "prot"},
)
model = scvi.model.TOTALVI(mdata, n_latent=20, gene_likelihood="nb")
model.train(max_epochs=400, early_stopping=True)

# shared latent space for neighbors / UMAP / clustering
mdata.obsm["X_totalVI"] = model.get_latent_representation()

# denoised (background-removed) protein — for MEASURED proteins only
denoised_rna, denoised_prot = model.get_normalized_expression(
    n_samples=25, return_mean=True,
    include_protein_background=False,   # True only to compare against raw ADT
)
```

Then, in scanpy: `sc.pp.neighbors(mdata, use_rep="X_totalVI")` → `sc.tl.umap` →
`sc.tl.leiden`. Differential expression via `model.differential_expression(...)`
gives RNA and protein Bayes factors in one framework.

**totalVI notes for this dataset class.**

- Give totalVI **raw** ADT counts (no prior DSB/CLR).
- If different libraries used **different ADT panels**, set a protein missing in a
  batch to all-zeros for that batch; totalVI detects and masks it. Optionally use
  `panel_key`.
- `get_normalized_expression` denoises **only measured** proteins. It will never
  produce a TROP2 value, because TROP2 is not in the panel — reaffirming §0.3.

---

## 7. Batch / sample integration

GSE154826-class data has strong patient/sample effects (cells tend to cluster by
patient). Handle explicitly:

- **totalVI:** pass `batch_key="sample"` (or patient); optionally add
  `categorical_covariate_keys` for tissue.
- **WNN:** batch-correct the RNA PCA (Harmony or Seurat RPCA) *before*
  `FindMultiModalNeighbors`; rescale ADT per batch if needed.
- Verify integration with kBET / iLISI and, negatively, that biological structure
  (immune lineages) is preserved.

---

## 8. Annotation strategy (modality-aware)

- **Immune lineages** — annotate with **ADT + RNA** (e.g. CD3/CD4/CD8 for T,
  CD19/CD20 for B, CD14/CD16 for myeloid, CD56 for NK). Protein is often cleaner
  than RNA for these surface markers.
- **Epithelial / tumor** — annotate with **RNA only**: `EPCAM`, `TACSTD2`,
  `CLDN4`, keratins (`KRT8/18/19`), etc. **Label these clusters explicitly as
  "RNA-defined epithelial; no ADT confirmation (TROP2 not in panel)."**
- Optionally project onto a reference (Azimuth, scANVI, or a totalVI reference) —
  but reference protein features must intersect this panel.

---

## 9. Reconciling RNA `TACSTD2`/`CLDN4` with absent TROP2 protein (core block)

This is the analytical crux of the dataset class. Follow this decision procedure
whenever epithelial RNA markers appear:

1. **Document the panel.** Confirm and record in the methods that the ADT panel
   contains **no TROP2 (`TACSTD2`) antibody**. State it as a measurement boundary.
2. **Never claim protein-level TROP2.** Report `TACSTD2` strictly as *RNA*
   evidence. Do not describe cells as "TROP2⁺" from this assay.
3. **Separate true epithelial cells from ambient/doublet artifacts.** A cluster is
   credibly epithelial only if, *after* ambient correction (§4), it shows a
   **coordinated epithelial program** (multiple keratins + `EPCAM` + `CLDN4`),
   **low immune ADT**, and passes doublet filtering (§3). Isolated low-level
   `TACSTD2` in an otherwise immune, high-CD45 cluster is most likely ambient
   contamination or a doublet, not epithelial identity.
4. **If protein-level TROP2 is genuinely required**, it must come from a *different*
   assay — flow cytometry, IHC/multiplex-IF, spatial proteomics, or a re-run with
   a TROP2 ADT added to the panel. Record this as an explicit limitation of the
   current data.
5. **Do not over-interpret totalVI protein denoising.** Denoised protein exists
   only for measured ADTs; the absence of a TROP2 channel is a hard limit, not a
   background-noise problem to be modeled away.

---

## 10. Differential expression

- **Preferred:** **pseudobulk** aggregated to patient/sample, tested with
  DESeq2 or edgeR — this respects the dominant unit of replication and avoids
  pseudoreplication across cells of one patient.
- **totalVI DE** (`model.differential_expression`) for a unified RNA+protein view
  on measured features.
- **MAST** for single-cell-level tests when needed; always include patient as a
  covariate / random effect.

---

## 11. Reproducibility

- Pin versions: `scanpy`, `scvi-tools`, `mudata`, `anndata`, `Seurat`,
  `SeuratObject`, DSB, CellBender. Record them in an environment file.
- Set and log all seeds (`scvi.settings.seed`, `set.seed`); note that GPU
  training is not bit-deterministic — record hardware.
- Log every parameter (HVG count, PCA dims, `n_latent`, epochs, thresholds).
- Persist intermediate objects (`.h5mu` / `.rds`) and the exact ADT panel table.

---

## 12. Methods-paragraph template (fill in, do not fabricate)

> Paired transcriptome (10x 3′ GEX) and surface-protein (ADT, immune-focused
> panel with isotype controls) measurements were demultiplexed by cell hashing
> (`<tool>`). Ambient RNA was removed with `<CellBender/SoupX>` and QC applied per
> modality (`<thresholds>`). ADT background was handled by `<DSB / totalVI>`.
> Multimodal integration used `<WNN / totalVI (n_latent=…, batch_key=…)>`, and
> clusters were annotated using ADT for immune lineages and RNA for epithelial
> identity. **The ADT panel did not include a TROP2 (TACSTD2) antibody; epithelial
> markers TACSTD2/CLDN4 are therefore reported at the RNA level only and were not
> confirmed at the protein level.** Differential expression used `<pseudobulk /
> totalVI>`. Software versions and seeds are in `<env file>`.

---

## 13. References

- Leader A.M., Grout J.A., et al. Single-cell analysis of human NSCLC lesions
  refines tumor classification and patient stratification (CITE-seq; **GSE154826**).
  *Cancer Cell*, 2021 (bioRxiv 2020.07.16.207605).
- Stoeckius M., et al. Simultaneous epitope and transcriptome measurement in
  single cells (**CITE-seq**). *Nat. Methods*, 2017.
- Stoeckius M., et al. Cell Hashing with barcoded antibodies. *Genome Biol.*, 2018.
- Hao Y., et al. Integrated analysis of multimodal single-cell data (**WNN**,
  Seurat v4). *Cell*, 2021.
- Gayoso A., Steier Z., et al. Joint probabilistic modeling of paired
  transcriptome and proteome (**totalVI**). *Nat. Methods*, 2021.
- Mulè M.P., et al. Normalizing and denoising protein expression data with
  **DSB**. *Nat. Commun.*, 2022.
- Fleming S.J., et al. **CellBender** remove-background. *Nat. Methods*, 2023.
- Young M.D., Behjati S. **SoupX**. *GigaScience*, 2020.
- Germain P.-L., et al. **scDblFinder**. *F1000Research*, 2021.

---
---

# 中文版 CITE-seq 方法手册

> **仅方法（METHODS ONLY）。** 本文档是一份可复现的分析协议，**不包含任何结果、
> 数字或生物学结论**，只提供针对 GSE154826 类 CITE-seq 数据的分析方法、参数与决策
> 规则。内容与上方英文部分一一对应。

---

## 0. 数据集类别与核心前提

**类别。** GSE154826（Leader 等，*Cancer Cell* 2021；Merad 实验室）是一套人非小细胞
肺癌（NSCLC）CITE-seq 资源。"GSE154826 类"数据集在单细胞层面包含：

- **RNA** —— 10x 3′ 基因表达（GEX）。
- **ADT** —— 抗体衍生标签，对应一个**以免疫细胞为主的表面标志物 panel**
  （CD3/CD4/CD8/CD14/CD19/CD56/PD-1 等），并含**同型对照（isotype control）**。
- **HTO** —— 细胞哈希标签（GSE154826 中针对 B2M/CD298），用于样本混合复用。
- 可选 **TCR**（VDJ）。
- 组织类型：肿瘤、非受累肺（nLung）、PBMC；GSE154826 中还有作为技术对照的
  **掺入小鼠细胞（spiked mouse cells）**。

**核心分析前提（动数据前务必阅读）。**
`TACSTD2`（编码 **TROP2** 蛋白）与 `CLDN4`（Claudin-4）是**上皮 / 癌细胞** RNA 标志物。
本类数据集的 ADT 是**免疫 panel**，因此**通常不含 TROP2 抗体**，即 TROP2 蛋白**未被
测量**。由此产生贯穿全流程的三条不变量：

1. **上皮身份只有 RNA 证据。** 无法用 ADT 在蛋白层面交叉验证 `TACSTD2`/`CLDN4`。
2. **两种模态不对称。** ADT 对免疫谱系信息量大；上皮/肿瘤结构几乎完全存在于 RNA。
3. **任何方法都无法"补出"未测量的蛋白。** WNN 与 totalVI 都不能凭空产生从未在 panel
   中的 TROP2 ADT；totalVI 只对**已测量**的蛋白做去噪。

整个流程中，请把以上三点当作不可违背的约束。

---

## 1. 输入与数据模型

按文库/样本收集：

| 项目 | 用途 |
| --- | --- |
| RNA 原始计数（基因 × 细胞） | 表达、HVG、聚类 |
| ADT 原始计数（蛋白 × 细胞），含同型对照 | 蛋白信号 + 背景估计 |
| HTO 原始计数 | 去复用（demultiplex） |
| **未过滤/原始液滴矩阵**（保留空液滴） | 环境 RNA 与 ADT 背景估计 |
| 样本/患者/组织/批次元数据 | 整合、伪散装（pseudobulk） |
| panel 定义（精确 ADT 名称） | 模态感知注释与 TROP2 注意事项 |

**推荐结构** —— 使用含两个模态的 MuData，并在 layers 中保留原始计数：

```python
import mudata as md, anndata as ad
mdata = md.MuData({"rna": rna_adata, "prot": prot_adata})
mdata["rna"].layers["counts"] = mdata["rna"].X.copy()   # 保留整数计数
mdata["prot"].layers["counts"] = mdata["prot"].X.copy()  # totalVI 需要原始 ADT
```

若走 Seurat/WNN 路线，使用含 `RNA`、`ADT`（及 `HTO`）assay 的对象。

---

## 2. 去复用（哈希）

1. 归一化 HTO（跨细胞 CLR）。
2. 用 `HTODemux`（Seurat）、`hashsolo`（scanpy）或 `GMM-Demux` 判定
   单细胞（singlet）/ 双细胞（doublet）/ 阴性（negative）。
3. 删除跨样本的**双细胞**与**阴性**；保留样本标签。
4. 若文库掺入小鼠细胞（GSE154826），在任何仅人类分析之前，通过人+鼠联合参考
   （barnyard）比对或小鼠基因占比阈值将其去除。

---

## 3. 质量控制（先分模态，再联合）

**RNA 每细胞：** `n_genes_by_counts`、`total_counts`、`pct_counts_mt`、
`pct_counts_ribo`、`pct_counts_hb`。阈值按组织设定（肿瘤细胞与淋巴细胞可接受的
线粒体/基因范围不同）——优先使用**基于 MAD** 的离群阈值，而非全局硬性数字。

**ADT 每细胞：** ADT 总计数、检出蛋白数，以及**同型对照水平**（同型偏高 →
黏附/聚集细胞 → 剔除）。

**联合规则：** 仅当细胞**同时**通过 RNA 与 ADT 质控且被判为哈希单细胞时才保留。

**双细胞（基于表达）：** 在 HTO 双细胞判定之外，另跑
`scDblFinder` / `DoubletFinder` / `scrublet`；跨模态双细胞（如免疫+上皮混合）是免疫
条形码中出现虚假上皮 RNA 的已知来源——见 §9。

---

## 4. 环境 / 背景校正（不可跳过）

**RNA 环境污染。** 裂解的肿瘤上皮会把 `TACSTD2`/`CLDN4` 转录本释放进液滴悬液，
因此这些基因可能仅因环境信号就在**真正的免疫细胞**中低水平出现。请使用未过滤液滴矩阵，
以 **CellBender `remove-background`**、**SoupX** 或 **DecontX** 估计并去除。此步骤是
在非上皮簇中解释上皮 RNA（§9）的前提。

**ADT 背景。** 抗体背景大且有结构。两条公认路线：

- **DSB 归一化**（Mulè 等 2022）：用空液滴 + 同型对照定义并扣除背景；将 DSB 值送入
  WNN 路线。
- **totalVI 蛋白背景混合模型**（§6B）：模型从**原始** ADT 计数中学习每个蛋白的背景
  分量——因此送入 totalVI 的计数**不要**预先做 DSB。

---

## 5. 归一化

- **RNA：** CP10k 后 `log1p`；或用 SCTransform / 解析 Pearson 残差做 HVG 选择
  （常取 2000–3000 个 HVG）。为 totalVI 与伪散装 DE 保留原始计数。
- **ADT（WNN 路线）：** DSB（有空液滴时优先）或 CLR；注明 CLR 的边（margin，常按
  跨细胞）。
- **ADT（totalVI 路线）：** 不归一化 —— totalVI 直接吃原始 ADT 计数。

---

## 6. 多模态整合 —— 两条路线都跑

并行运行 WNN 与 totalVI 并相互比较，是推荐的质控手段：两者应在免疫结构上一致；不一致
往往提示模态权重或背景模型存在问题。

### 6A. 加权最近邻（WNN，Seurat v4；Hao 等 2021）

WNN 为每个细胞学习**各模态的权重**并构建联合图。

```r
library(Seurat)
# 假设 obj 已含 SCT 处理的 RNA 与 DSB/CLR 归一化的 ADT
DefaultAssay(obj) <- "SCT"
obj <- RunPCA(obj, reduction.name = "pca")
DefaultAssay(obj) <- "ADT"
VariableFeatures(obj) <- rownames(obj[["ADT"]])          # panel 小：全部使用
obj <- ScaleData(obj) |> RunPCA(reduction.name = "apca")

obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("pca", "apca"),
  dims.list      = list(1:30, 1:min(20, nrow(obj[["ADT"]]) - 1)),
  modality.weight.name = c("RNA.weight", "ADT.weight")
)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap")
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3)  # Leiden
```

**模态权重注意事项。** 检查 `RNA.weight` / `ADT.weight`。由于没有 TROP2（也没有任何
上皮）ADT，上皮细胞几乎完全由 **RNA 权重**驱动；WNN 无法借蛋白来锐化上皮边界。这是
预期行为而非 bug——但意味着此处上皮聚类质量本质上是一个 RNA-only 的属性。

### 6B. totalVI（scvi-tools；Gayoso 等 2021）

totalVI 是联合生成模型，可（a）给出共享隐空间，（b）通过区分前景/背景对 ADT 去噪，
（c）校正批次。

```python
import scvi

scvi.model.TOTALVI.setup_mudata(
    mdata,
    rna_layer="counts",
    protein_layer="counts",
    batch_key="sample",                 # 患者/样本；此处批次效应强
    modalities={"rna_layer": "rna", "protein_layer": "prot"},
)
model = scvi.model.TOTALVI(mdata, n_latent=20, gene_likelihood="nb")
model.train(max_epochs=400, early_stopping=True)

# 用于近邻 / UMAP / 聚类的共享隐空间
mdata.obsm["X_totalVI"] = model.get_latent_representation()

# 去噪（去背景）蛋白 —— 仅针对“已测量”的蛋白
denoised_rna, denoised_prot = model.get_normalized_expression(
    n_samples=25, return_mean=True,
    include_protein_background=False,   # 仅在与原始 ADT 对比时设为 True
)
```

随后在 scanpy 中：`sc.pp.neighbors(mdata, use_rep="X_totalVI")` → `sc.tl.umap` →
`sc.tl.leiden`。用 `model.differential_expression(...)` 可在同一框架内得到 RNA 与蛋白
的贝叶斯因子。

**本类数据集的 totalVI 要点。**

- 给 totalVI 的是**原始** ADT 计数（不要预先 DSB/CLR）。
- 若不同文库使用了**不同 ADT panel**，将某批次缺失的蛋白在该批次全置 0；totalVI 会
  自动检测并屏蔽。可选用 `panel_key`。
- `get_normalized_expression` **只**对已测量蛋白去噪。它永远不会产出 TROP2 值，因为
  TROP2 不在 panel 中——再次印证 §0 第 3 条。

---

## 7. 批次 / 样本整合

GSE154826 类数据的患者/样本效应很强（细胞倾向按患者聚类）。需显式处理：

- **totalVI：** 传入 `batch_key="sample"`（或患者）；可另加 `categorical_covariate_keys`
  表示组织。
- **WNN：** 在 `FindMultiModalNeighbors` **之前**对 RNA PCA 做批次校正（Harmony 或
  Seurat RPCA）；必要时按批次重标定 ADT。
- 用 kBET / iLISI 验证整合效果，并反向确认生物学结构（免疫谱系）得以保留。

---

## 8. 注释策略（模态感知）

- **免疫谱系** —— 用 **ADT + RNA** 注释（如 T 细胞 CD3/CD4/CD8，B 细胞 CD19/CD20，
  髓系 CD14/CD16，NK CD56）。对这些表面标志物，蛋白往往比 RNA 更干净。
- **上皮 / 肿瘤** —— **仅用 RNA** 注释：`EPCAM`、`TACSTD2`、`CLDN4`、角蛋白
  （`KRT8/18/19`）等。**必须把这些簇明确标注为"RNA 定义的上皮；无 ADT 确认
  （panel 中无 TROP2）"。**
- 可选投影到参考（Azimuth、scANVI 或 totalVI 参考）——但参考的蛋白特征须与本 panel
  交集。

---

## 9. 调和 RNA `TACSTD2`/`CLDN4` 与缺失的 TROP2 蛋白（核心环节）

这是本类数据集的分析关键。每当出现上皮 RNA 标志物时，按以下决策流程处理：

1. **记录 panel。** 在方法中确认并写明 ADT panel **不含 TROP2（`TACSTD2`）抗体**，
   将其作为一条测量边界陈述。
2. **切勿声称蛋白层面的 TROP2。** `TACSTD2` 严格作为 *RNA* 证据报告；不要基于本实验
   把细胞描述为"TROP2⁺"。
3. **区分真上皮细胞与环境/双细胞假象。** 只有当一个簇在**环境校正之后**（§4）呈现
   **协调一致的上皮程序**（多个角蛋白 + `EPCAM` + `CLDN4`）、**低免疫 ADT**，并通过
   双细胞过滤（§3）时，才可信地判为上皮。若在一个高 CD45、整体免疫的簇中孤立地出现
   低水平 `TACSTD2`，最可能是环境污染或双细胞，而非上皮身份。
4. **若确需蛋白层面的 TROP2**，必须来自*另一种*实验——流式、IHC/多重免疫荧光、空间
   蛋白组，或在 panel 中加入 TROP2 ADT 重跑。请将此作为当前数据的明确局限记录。
5. **不要过度解读 totalVI 的蛋白去噪。** 去噪蛋白只对已测量 ADT 存在；缺少 TROP2 通道
   是硬性限制，而非可被模型消除的背景噪声问题。

---

## 10. 差异表达

- **首选：** 聚合到患者/样本层面的**伪散装（pseudobulk）**，用 DESeq2 或 edgeR 检验
  ——这尊重主要的重复单元，避免同一患者众多细胞造成的伪重复。
- **totalVI DE**（`model.differential_expression`）：在已测量特征上统一给出 RNA+蛋白
  视图。
- 需要单细胞级检验时用 **MAST**；始终把患者作为协变量/随机效应纳入。

---

## 11. 可复现性

- 固定版本：`scanpy`、`scvi-tools`、`mudata`、`anndata`、`Seurat`、`SeuratObject`、
  DSB、CellBender，并记录于环境文件。
- 设置并记录所有随机种子（`scvi.settings.seed`、`set.seed`）；注意 GPU 训练非比特级
  确定，需记录硬件。
- 记录每个参数（HVG 数、PCA 维度、`n_latent`、epochs、阈值）。
- 持久化中间对象（`.h5mu` / `.rds`）与精确的 ADT panel 列表。

---

## 12. 方法段落模板（填空，勿编造）

> 配对的转录组（10x 3′ GEX）与表面蛋白（ADT，以免疫为主、含同型对照的 panel）测量，
> 经细胞哈希（`<工具>`）去复用。用 `<CellBender/SoupX>` 去除环境 RNA，并按模态施加质控
> （`<阈值>`）。ADT 背景以 `<DSB / totalVI>` 处理。多模态整合采用
> `<WNN / totalVI（n_latent=…, batch_key=…）>`，并用 ADT 注释免疫谱系、用 RNA 判定
> 上皮身份。**ADT panel 不含 TROP2（TACSTD2）抗体；因此上皮标志物 TACSTD2/CLDN4 仅在
> RNA 层面报告，未在蛋白层面确认。** 差异表达采用 `<pseudobulk / totalVI>`。软件版本
> 与随机种子见 `<环境文件>`。

---

## 13. 参考文献

- Leader A.M., Grout J.A. 等。人 NSCLC 病灶单细胞分析，精化肿瘤分类与患者分层
  （CITE-seq；**GSE154826**）。*Cancer Cell*，2021（bioRxiv 2020.07.16.207605）。
- Stoeckius M. 等。单细胞中同时测量表位与转录组（**CITE-seq**）。*Nat. Methods*，2017。
- Stoeckius M. 等。基于条码抗体的 Cell Hashing。*Genome Biol.*，2018。
- Hao Y. 等。多模态单细胞数据的整合分析（**WNN**，Seurat v4）。*Cell*，2021。
- Gayoso A.、Steier Z. 等。配对转录组与蛋白组的联合概率建模（**totalVI**）。
  *Nat. Methods*，2021。
- Mulè M.P. 等。用 **DSB** 归一化并去噪蛋白表达数据。*Nat. Commun.*，2022。
- Fleming S.J. 等。**CellBender** remove-background。*Nat. Methods*，2023。
- Young M.D.、Behjati S.。**SoupX**。*GigaScience*，2020。
- Germain P.-L. 等。**scDblFinder**。*F1000Research*，2021。
