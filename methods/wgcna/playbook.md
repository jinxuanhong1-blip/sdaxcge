# WGCNA / hdWGCNA playbook: TACSTD2–CLDN4 junction module / TACSTD2–CLDN4 连接模块共表达网络手册

> **Scope / 范围.** Methods only. This document specifies how to *discover* a TACSTD2–CLDN4 junction co-expression module in bulk RNA-seq (WGCNA) and scRNA-seq (hdWGCNA), then *apply* the module eigengene (ME) to immune-checkpoint inhibitor (ICI) response. It does not report results, invent accessions, or claim a validated biomarker.
>
> 仅限方法。本文规定如何在 bulk RNA-seq（WGCNA）和 scRNA-seq（hdWGCNA）中*发现* TACSTD2–CLDN4 连接共表达模块，再将模块特征基因（module eigengene, ME）*应用于*免疫检查点抑制剂（ICI）疗效关联。不报告结果、不编造登录号、不宣称已验证的生物标志物。

Related playbooks (when present): `methods/ml/playbook.md` (small-cohort prediction), `methods/causal/playbook.md` (TACSTD2 → CD8 → ICI), `methods/survival/` (time-to-event), `methods/batch/` (ComBat/sva), `methods/scrna/` (Seurat/Scanpy).

配套手册（若已存在）：`methods/ml/playbook.md`（小样本预测）、`methods/causal/playbook.md`（TACSTD2 → CD8 → ICI）、`methods/survival/`（生存）、`methods/batch/`（批次）、`methods/scrna/`（单细胞预处理）。

---

## 1. Estimand / 研究问题

Write these down *before* inspecting module–trait *p*-values.

在查看模块–性状 *p* 值之前先写清以下内容。

- **Population:** lung cancer (NSCLC/SCLC; LUAD/LUSC stated separately), ICI regimen (PD-1 / PD-L1 / CTLA-4 / chemo-IO), treatment line, specimen timing (pretreatment unless on-treatment is the stated use).
- **Exposure / module:** a *junction module* — a co-expression module in which **TACSTD2** (TROP2; ENSG00000184292; UniProt P09758) and **CLDN4** (ENSG00000189143) both have high module membership (MM / kME), *or* a module enriched for tight-junction / epithelial-barrier genes that contains both seeds. Mouse orthologs: `Tacstd2`, `Cldn4`.
- **Unit of analysis:** patient (donor). Multiple aliquots, regions, or cells from one patient stay together.
- **Primary ICI endpoint (pick one and freeze it):** RECIST objective response, durable clinical benefit at a fixed landmark, MPR/pCR (neoadjuvant), or a time-to-event endpoint (PFS/OS). Do not mix definitions after seeing results.
- **Secondary biology:** ME versus immune signatures (T-cell inflamed GEP, IFNG, CD8, CYT, TLS, TGF-β exclusion) and versus TACSTD2 / CLDN4 single-gene expression.
- **What this is *not*:** a treatment-predictive claim. A single-arm ICI cohort can only support a *prognostic* association of the ME with outcome on ICI. Predictive (treatment-interaction) claims require a comparator arm (see `methods/causal/` and the OAK/POPLAR note in `methods/ml/`).

- **目标人群：** 肺癌（NSCLC/SCLC；LUAD/LUSC 分开写）、ICI 方案（PD-1 / PD-L1 / CTLA-4 / 化疗联合免疫）、治疗线次、取样时点（除非明确分析治疗中样本，否则仅用治疗前）。
- **暴露 / 模块：** *连接模块*——TACSTD2 与 CLDN4 均具有高模块隶属度（MM / kME）的共表达模块，*或*富集紧密连接 / 上皮屏障基因且同时包含这两个种子基因的模块。小鼠同源基因：`Tacstd2`、`Cldn4`。
- **分析单位：** 患者（供体）。同一患者的多份样本、多区域或细胞必须绑定。
- **主要 ICI 终点（选一个并冻结）：** RECIST 客观缓解、固定时点的持久临床获益、新辅助 MPR/pCR，或生存终点（PFS/OS）。不可在看结果后再混用定义。
- **次要生物学：** ME 与免疫签名（T 细胞炎症 GEP、IFNG、CD8、CYT、TLS、TGF-β 排斥）以及与单基因 TACSTD2 / CLDN4 的关系。
- **这不是什么：** 治疗预测性主张。单臂 ICI 队列只能支持 ME 在 ICI 治疗下的*预后*关联。预测性（治疗交互）主张需要对照臂。

**Prior art to test, not assume:** Bessede et al., *Clin Cancer Res* 2024 (PMID 38048058) — high TACSTD2 with atezolizumab primary resistance and less T-cell infiltration (OAK/POPLAR; RNA often EGA-controlled). Co-expression of TACSTD2 with CLDN4 is a hypothesis, not a given. The two genes may land in different modules; that is a reportable negative, not a reason to force a merge.

**待检验、不可预设的既往证据：** Bessede 等（PMID 38048058）——高 TACSTD2 与阿替利珠单抗原发耐药及 T 细胞浸润减少相关（OAK/POPLAR；RNA 常为 EGA 受控）。TACSTD2 与 CLDN4 共表达是假设，不是事实。二者可能落入不同模块；这是可报告的阴性结果，不是强行合并的理由。

---

## 2. Sample-size limits (hard gates) / 样本量硬门槛

WGCNA estimates *thousands* of pairwise correlations. Those correlations are noisy when *n* is small. Official WGCNA FAQ (Langfelder & Horvath; table updated December 2017): **do not attempt WGCNA on fewer than 15 samples**; **≥20 is preferred**. Finer results (hub ranking, module splitting) degrade first.

WGCNA 要估计*数千对*相关，小样本下噪声很大。官方 FAQ（Langfelder & Horvath；2017 年 12 月更新）：**少于 15 个样本不要做 WGCNA**；**最好 ≥20**。更细的结果（hub 排序、模块拆分）最先不可靠。

Companion guard: `methods/wgcna/sample_size_gates.R`. Passing a gate is permission to *run a method*, not evidence that the cohort is adequate.

配套门槛脚本：`methods/wgcna/sample_size_gates.R`。通过门槛只表示*允许运行该方法*，不表示队列足够。

### 2.1 Network construction (bulk) / 网络构建（bulk）

| Patient *n* (unique donors after QC) | Network construction | Soft-threshold if scale-free *R*² < 0.80 | Claims allowed |
| --- | --- | --- | --- |
| *n* < 15 | **Forbidden.** Do not call `blockwiseModules` / `pickSoftThreshold`. | — | Locked junction *score* + pairwise TACSTD2–CLDN4 only (Section 7). |
| 15 ≤ *n* < 20 | Exploratory signed (or signed-hybrid) network only. | Signed power **18**; unsigned/hybrid **9**. | Module membership of the two seeds. No hub-gene list as a “signature.” No ME–ICI model with covariates. |
| 20 ≤ *n* < 30 | Allowed, conservative. | Signed **16**; unsigned/hybrid **8**. | Module identity + unadjusted ME–trait. Label **exploratory**. |
| 30 ≤ *n* < 40 | Preferred for ME–trait. | Signed **14**; unsigned/hybrid **7**. | Standard module–trait with BH-FDR. Still no new multigene classifier. |
| *n* ≥ 40 | Standard. | Signed **12**; unsigned/hybrid **6** (or data-driven *β* if *R*² ≥ 0.80 and mean connectivity is not huge). | Discovery cohort for a *locked* module to project elsewhere. |

Use the FAQ power *only* when the scale-free fit fails for reasonable powers (<30 signed; <15 unsigned/hybrid) *or* when a strong global driver keeps mean connectivity in the hundreds. Prefer the smallest *β* that reaches *R*² ≥ 0.80 with a descending mean-connectivity curve.

仅当合理幂次下无标度拟合失败（signed <30；unsigned/hybrid <15），或强全局驱动使平均连接度仍达数百时，才改用 FAQ 幂次表。优先选择使 *R*² ≥ 0.80 且平均连接度下降的最小 *β*。

**Do not discover the module and test ME versus ICI in the same *n* = 16 cohort.** Open ICI bulk sets in this project are typically *n* ≈ 16–28 (e.g. GSE126044, GSE166449, GSE135222, GSE207422 bulk). Those sizes sit at or below the exploratory band. Default design:

**禁止在同一个 n=16 队列里既发现模块又检验 ME 与 ICI 的关联。** 本项目开放 ICI bulk 队列通常约 16–28 例。默认设计：

1. **Discover** the junction module in a large, treatment-naive transcriptome (TCGA LUAD and/or LUSC, run separately; optional CPTAC protein network). No ICI labels are required for discovery.
2. **Lock** the gene list, the ME sign convention, and the scoring recipe (Section 6.3).
3. **Apply** the locked score / projected ME in each ICI cohort. That application still obeys the *association* gates in Section 2.2.

1. 在大型、未接受 ICI 的转录组中**发现**连接模块（TCGA LUAD / LUSC 分开跑；可选 CPTAC 蛋白网络）。发现阶段不需要 ICI 标签。
2. **锁定**基因列表、ME 符号约定和评分配方（第 6.3 节）。
3. 在各 ICI 队列中**应用**锁定评分 / 投影 ME。应用阶段仍遵守第 2.2 节的关联门槛。

Consensus WGCNA across ICI studies is allowed only if *each* contributing study has *n* ≥ 20 after QC, or if studies are concatenated *after* a documented, outcome-blind batch correction fitted without ICI labels (see `methods/batch/`). Concatenating *n* = 16 + 22 to “make 38” without batch handling is not a valid *n* = 38 network.

仅当每个参与研究质控后 *n* ≥ 20，或在无结局标签的批次校正后再合并时，才允许跨 ICI 研究做共识 WGCNA。把 n=16 与 n=22 硬拼成“38”而不处理批次，不能当作 n=38 网络。

### 2.2 Module eigengene versus ICI / 模块特征基因与 ICI

These gates apply to **locked** scores as well as to de novo MEs.

以下门槛同时适用于**锁定**评分和从头计算的 ME。

| Analysis | Minimum *n* | Rule |
| --- | --- | --- |
| Spearman / bicor (ME vs continuous trait or gene) | *n* ≥ 4 complete pairs | Project convention: *n* < 4 → NA. Report *n*, *ρ*, *p*, BH-FDR. |
| Two-group ME contrast (response vs no response) | ≥5 per arm to *run*; ≥10 per arm to *report* as a primary contrast | Mann–Whitney U + rank-biserial *r*. 5–9 per arm: exploratory supplement only. <5 in either arm: do not test. |
| Logistic ME + covariates | ≥10 events per fitted coefficient | Typical open ICI bulks cannot support this. |
| Cox PH (ME vs PFS/OS) | ≥10 events | Check PH. Median-split log-rank is descriptive, not a new cutoff discovery. |
| New multigene “junction signature” trained on ICI labels | Not allowed at these *n* | See `methods/ml/`: do not train a 20-gene signature on *n* = 16. |

Binary or ordinal traits: if using `bicor`, set `maxPOutliers = 0.05` (or 0.10) and `robustY = FALSE` (Langfelder & Horvath 2011; WGCNA FAQ). Otherwise bicor can treat one class as “outliers.”

二分类或有序性状：若用 `bicor`，设 `maxPOutliers = 0.05`（或 0.10）且 `robustY = FALSE`。否则 bicor 可能把其中一个类别当成“离群值”。

### 2.3 hdWGCNA (scRNA / spatial) / 单细胞 hdWGCNA

hdWGCNA (Morabito et al., *Cell Reports Methods* 2023) builds metacells because WGCNA is sensitive to dropout. Metacells are **not** extra biological replicates.

hdWGCNA（Morabito 等，2023）用 metacell 降低 dropout。Metacell **不是**额外的生物学重复。

| Gate | Default | Fail action |
| --- | --- | --- |
| Grouping | Always `group.by = c(cell_type, Sample)` (or equivalent). Never mix donors or lineages in one metacell. | Abort. |
| Neighbor size *k* | 20–75; start at 25. Smaller *k* for small datasets. | Tune; do not go below 10. |
| `max_shared` | 10 (tutorial default); keep overlap modest. | If almost all cells are shared, raise *k* or lower `max_shared`. |
| `min_cells` per sample × cell type | ≥50 (prefer ≥ *k* × 2). Exclude rare types (Morabito: underrepresented types yield bad metacells). | Drop that sample × type; do not lower `min_cells` to “keep everyone.” |
| Metacells for *network construction* | ≥20–30 metacells in the analysis group (e.g. epithelial) | Do not construct. Project a locked module (Section 6.3) or use the Section 7 score. |
| Epithelial / tumor cells per *donor* (junction module) | Enough to pass `min_cells` in most ICI donors | If many ICI donors fail, discover in a larger atlas and project. |
| Donor *n* for ME–ICI | Same table as Section 2.2, counting **donors**, not cells or metacells | Cells are pseudoreplicates. |

**Compartment:** construct the junction-module network in **epithelial / malignant** cells only. TACSTD2 and CLDN4 are epithelial products; an immune-cell network will not define a junction module. Immune, stromal, and endothelial compartments may be run as *separate* hdWGCNA experiments (different `wgcna_name`) if a question about those lineages is prespecified.

**区室：** 连接模块网络只在**上皮 / 恶性**细胞中构建。TACSTD2 与 CLDN4 是上皮产物；免疫细胞网络不能定义连接模块。免疫、基质、内皮可作为*独立* hdWGCNA 实验（不同 `wgcna_name`），但须预先提出问题。

Project file-size rule: skip any single file >2 GB and record it (e.g. some GSE131907 matrices). Do not invent a substitute matrix.

项目文件大小规则：单个文件 >2 GB 则跳过并记录。不得编造替代矩阵。

---

## 3. Inputs and eligibility / 输入与入选

Do not invent GEO / ArrayExpress / PRIDE / CPTAC / EGA accessions. Verify each accession on its official page. Use open processed matrices. Skip FASTQ / SRA / raw MS. Controlled data (EGA / dbGaP / GSA-Human, including OAK/POPLAR RNA): catalog only unless access is granted.

不要编造登录号。在官方页面核验每一条。只用开放的处理后矩阵。跳过 FASTQ / SRA / 原始质谱。受控数据（含 OAK/POPLAR RNA）仅编目，除非已获授权。

**Feature-space check (mandatory):** TACSTD2 and CLDN4 (or mouse orthologs) must be present. Targeted panels that lack both genes (examples previously flagged in this project: some nCounter / Oncomine immune panels) are ineligible for expression analysis. Report coverage; do not impute missing genes.

**特征空间检查（强制）：** 必须能测到 TACSTD2 与 CLDN4（或小鼠同源基因）。缺少二者的靶向 panel 不得做表达分析。报告覆盖率；不得填补缺失基因。

**Discovery-friendly (large *n*, usually no ICI labels):** TCGA LUAD and LUSC RNA-seq (separate networks); optional CPTAC LUAD / LSCC protein (treatment-naive; protein WGCNA is allowed but is not an ICI test). Large scRNA atlases only if the matrix is open and ≤2 GB, or if a processed Seurat/AnnData object is already local.

**发现用（大样本，通常无 ICI 标签）：** TCGA LUAD、LUSC RNA-seq（分开建网）；可选 CPTAC LUAD / LSCC 蛋白（初治；可做蛋白 WGCNA，但不是 ICI 检验）。大型单细胞图谱仅在矩阵开放且 ≤2 GB、或本地已有处理后对象时使用。

**Application (ICI labels):** open bulk / scRNA cohorts that pass the feature-space check and have response or survival labels (e.g. pretreatment anti-PD-1 RNA with RECIST; neoadjuvant PD-1+chemo with MPR; scRNA with MPR/RECIST). Mouse ICI RNA uses `Tacstd2` / `Cldn4` and the same gates.

**应用（有 ICI 标签）：** 通过特征空间检查且具有疗效或生存标签的开放 bulk / scRNA 队列。小鼠 ICI RNA 使用 `Tacstd2` / `Cldn4` 及同一套门槛。

**Batch / purity:** record study, platform, site, histology, specimen site, tumor purity (or ESTIMATE / IHC), PD-L1, TMB, driver (EGFR/ALK/STK11/KEAP1 when present), and treatment. Cross-tabulate response × batch before any ME test. If response is aliased with batch, stop (see `methods/ml/` and `methods/batch/`).

**批次 / 纯度：** 记录研究、平台、中心、组织学、取材部位、肿瘤纯度、PD-L1、TMB、驱动基因和治疗。在任何 ME 检验前交叉表“疗效 × 批次”。若疗效与批次完全混淆，则停止。

---

## 4. Bulk WGCNA procedure / Bulk WGCNA 流程

### 4.1 Expression matrix / 表达矩阵

1. One row per gene symbol (or ENSG, then map). Collapse duplicate symbols by **maximum mean** (project convention). Keep TACSTD2 and CLDN4 even if they fail a variance filter.
2. One column per **patient**. Prefer pretreatment primary tumor. If a patient has multiple tumors, prespecify (e.g. one index lesion) — do not average after seeing outcome.
3. Transform: counts → `log2(CPM + 1)` or DESeq2 VST; TPM → `log2(TPM + 1)`. Do not mix scales inside one network.
4. Filter genes by median absolute deviation or variance (e.g. top 5,000–10,000, or MAD > 0), **without** using ICI labels or other traits. WGCNA FAQ: do not pre-filter by differential expression versus the trait you will later correlate.
5. Sample QC: `goodSamplesGenes`; hierarchical clustering of samples; remove only extreme outliers with a prespecified height cut, documented in the manifest. Recheck *n* against Section 2.1 after removal.

1. 每基因一行（符号或 ENSG）。重复符号按**最大均值**合并。即使方差低也保留 TACSTD2 与 CLDN4。
2. 每**患者**一列。优先治疗前原发灶。同一患者多灶须预先规定选取规则，不可看结局后再平均。
3. 变换：counts → `log2(CPM + 1)` 或 DESeq2 VST；TPM → `log2(TPM + 1)`。同一网络内不混用尺度。
4. 按 MAD 或方差过滤基因（例如 top 5,000–10,000），**不得**使用 ICI 标签。FAQ：不要按稍后要关联的性状做差异表达预过滤。
5. 样本质控：`goodSamplesGenes`；样本层次聚类；仅按预设高度剔除极端离群并写入清单。剔除后按第 2.1 节重新核对 *n*。

### 4.2 Network and modules / 网络与模块

Default (preferred): **signed** or **signed hybrid** network, `corType = "bicor"`, `maxPOutliers = 0.05`, `pearsonFallback = "individual"`, `TOMType` matching the network, `minModuleSize = 30`, `deepSplit = 2`, `mergeCutHeight = 0.25`, `maxBlockSize` large enough to avoid arbitrary blocks (e.g. 20,000 if memory allows).

默认（优先）：**signed** 或 **signed hybrid**，`corType = "bicor"`，`maxPOutliers = 0.05`，`pearsonFallback = "individual"`，TOM 类型与网络一致，`minModuleSize = 30`，`deepSplit = 2`，`mergeCutHeight = 0.25`，`maxBlockSize` 尽量大以免人为分块。

```r
# Illustrative skeleton — not a result.
sft <- pickSoftThreshold(
  datExpr,
  powerVector = c(1:10, seq(12, 30, 2)),
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(maxPOutliers = 0.05)
)
# If sft$fitIndices R^2 never exceeds 0.80 at beta < 30, use Section 2.1 FAQ power.
net <- blockwiseModules(
  datExpr,
  power = beta,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  maxPOutliers = 0.05,
  pearsonFallback = "individual",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  maxBlockSize = 20000
)
MEs <- net$MEs
kME <- signedKME(datExpr, MEs, corFnc = "bicor",
                 corOptions = "use = 'p'; maxPOutliers = 0.05")
```

Enable `enableWGCNAThreads()` only on a machine you control; set a seed for any stochastic block split.

仅在可控机器上启用 `enableWGCNAThreads()`；对随机分块设种子。

### 4.3 Calling the junction module / 判定连接模块

Prespecified rules, in order:

预先规定的判定顺序：

1. **Same-module rule:** TACSTD2 and CLDN4 are both assigned to a non-grey module *M*. Require MM (absolute kME for *M*) ≥ 0.60 for both. If both ≥ 0.80, they are *hub-class* members of *M*. *M* is the primary junction module.
2. **Split-module rule:** the two seeds are in different non-grey modules. Report both modules, both MM values, and the ME–ME correlation. Do **not** merge modules post hoc to force co-membership. Optional secondary score: mean of the two locked MEs, declared *before* ICI tests.
3. **Grey / absent rule:** either seed is grey or was filtered. Stop calling a de novo junction module. Fall back to Section 7.
4. **Enrichment (descriptive):** test *M* for GO / Reactome tight-junction and epithelial-barrier terms on the background of genes that entered the network. Enrichment supports interpretation; it does **not** override rules 1–3.

1. **同模块规则：** 两基因同属非灰色模块 *M*，且对 *M* 的 |kME| 均 ≥ 0.60。若均 ≥ 0.80，则为 hub 级成员。*M* 为主要连接模块。
2. **拆分规则：** 两基因落在不同非灰色模块。报告两模块、MM 及 ME–ME 相关。**禁止**事后合并以强行同属。可选次要评分：两个锁定 ME 的均值，须在 ICI 检验前声明。
3. **灰色 / 缺失规则：** 任一种子为灰色或被过滤。停止称其为发现的连接模块，改走第 7 节。
4. **富集（描述性）：** 在进入网络的基因背景下检验紧密连接 / 上皮屏障条目。富集只支持解释，**不能**覆盖规则 1–3。

Hub genes for reporting are the top kME genes *inside the locked module*, not a new ICI-supervised list. Do not feed ICI *p*-values into hub selection.

报告用 hub 是锁定模块内 kME 最高的基因，不是经 ICI 监督重排的列表。

---

## 5. hdWGCNA procedure (scRNA) / hdWGCNA 流程（单细胞）

Prerequisites: a processed Seurat object (v5: `UpdateSeuratObject` if needed) with normalization, variable genes, scaling, PCA, batch-aware reduction (e.g. Harmony), UMAP, and **cell-type / malignant labels**. Do not subset the object after `SetupForWGCNA` (hdWGCNA does not support that). QC and annotation belong in `methods/scrna/`; this playbook starts from a frozen annotation.

前提：已处理的 Seurat 对象（必要时 `UpdateSeuratObject`），含标准化、高变基因、缩放、PCA、批次感知降维（如 Harmony）、UMAP 以及**细胞类型 / 恶性标注**。`SetupForWGCNA` 之后不要再 subset。质控与注释见 `methods/scrna/`；本手册从冻结注释开始。

### 5.1 Setup and metacells / 设置与 metacell

```r
# Epithelial / malignant cells only for the junction experiment.
epi <- subset(seurat_obj, subset = compartment == "epithelial")  # prespecified label
epi <- SetupForWGCNA(
  epi,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "junction_epi"
)
# After setup, verify TACSTD2 and CLDN4 are in GetWGCNAGenes(epi).
# If a seed fails the fraction cut, add it via gene_select = "custom"
# (union of fraction genes + c("TACSTD2","CLDN4")) — document the exception.

epi <- MetacellsByGroups(
  seurat_obj = epi,
  group.by = c("cell_type", "Sample"),
  reduction = "harmony",
  k = 25,
  max_shared = 10,
  min_cells = 50,
  ident.group = "cell_type"
)
epi <- NormalizeMetacells(epi)
```

If `GetWGCNAGenes` drops a seed, switch to `gene_select = "custom"` with the union of fraction-passing genes and the two seeds. Never add ICI-DE genes at this step.

若分数阈值丢掉种子基因，改用 `gene_select = "custom"`，取“过阈值基因 ∪ 两个种子”。此步不得加入 ICI 差异基因。

Count metacells in the epithelial group. If <20, **do not** `ConstructNetwork`. Use projection (Section 6.3) or Section 7.

统计上皮组 metacell 数。若 <20，**不要** `ConstructNetwork`，改用投影或第 7 节。

### 5.2 Network, MEs, connectivity / 网络、ME、连接度

Follow the official hdWGCNA single-cell tutorial (`TestSoftPowers` → `ConstructNetwork` → `ModuleEigengenes` → `ModuleConnectivity` → `GetHubGenes`). Use `networkType = "signed"`. Harmonize MEs by sample when `ModuleEigengenes(..., group.by.vars = "Sample")` is appropriate (batch across donors).

按官方单细胞教程：`TestSoftPowers` → `ConstructNetwork` → `ModuleEigengenes` → `ModuleConnectivity` → `GetHubGenes`。`networkType = "signed"`。跨供体批次明显时，`ModuleEigengenes(..., group.by.vars = "Sample")` 对 ME 做 Harmony。

Apply the same junction-module calling rules (Section 4.3) to hdWGCNA module assignments and kME.

连接模块判定规则与第 4.3 节相同。

**Donor-level ME for ICI:** average the (harmonized) module eigengene across epithelial cells — or across epithelial metacells — *within each donor*. That donor vector is the only object that may enter Section 2.2 tests. Do not treat cells as independent ICI outcomes.

**供体水平 ME：** 在每个供体内，对外皮细胞（或上皮 metacell）的（和谐化）ME 取平均。只有这条供体向量可以进入第 2.2 节检验。不可把细胞当作独立的 ICI 结局。

Optional: `ProjectModules` from a large atlas network onto a smaller ICI scRNA object when the ICI object fails metacell gates. Projection uses the locked gene weights; it is not a second discovery.

可选：当 ICI 对象未过 metacell 门槛时，用 `ProjectModules` 把大型图谱网络投影到较小的 ICI scRNA。投影使用锁定权重，不是第二次发现。

Spatial transcriptomics: hdWGCNA spatial tutorials apply the same sample-size logic to **spots aggregated by sample and niche**, not to raw sparse spots. GeoMx AOIs are closer to bulk; use Section 4 if each AOI class has enough *donors*.

空间转录组：按**样本 × niche** 聚合 spot，而不是对原始稀疏 spot 建网。GeoMx AOI 更接近 bulk；若每类 AOI 的*供体*数足够，走第 4 节。

---

## 6. Correlating the eigengene with ICI response / 将特征基因与 ICI 疗效关联

### 6.1 Prespecify the primary module / 预先指定主模块

The primary ME is:

主 ME 定义为：

- the locked discovery module that satisfied Section 4.3 rule 1, **or**
- the locked Section 7 score if discovery failed or ICI *n* forbids de novo WGCNA.

- 满足第 4.3 节规则 1 的锁定发现模块，**或**
- 发现失败、或 ICI 样本量禁止从头 WGCNA 时的第 7 节锁定评分。

Do not scan all module–trait *p*-values and then name the smallest one “the junction module.” Other modules may be reported as a BH-FDR–controlled secondary heatmap.

禁止先扫遍所有模块–性状 *p* 值，再把最小者命名为“连接模块”。其他模块可作为 BH-FDR 控制的次要热图。

### 6.2 Tests / 检验

For each ICI cohort, on the donor-level primary ME (sign-aligned so that higher ME means higher TACSTD2–CLDN4 program; flip the ME if `cor(ME, TACSTD2) < 0`):

在每个 ICI 队列中，对供体水平主 ME（符号对齐：ME 越高表示 TACSTD2–CLDN4 程序越强；若 `cor(ME, TACSTD2) < 0` 则翻转）：

1. **Binary response:** Mann–Whitney U, rank-biserial *r*, exact *n* per arm. Optional ROC AUC of the *locked* ME is descriptive — no threshold tuning on the same cohort (`methods/ml/`).
2. **Continuous immune scores:** Spearman *ρ* (pairwise complete). Same signature definitions as the project library (mean *z* of genes present; store `used/total`).
3. **Survival:** Cox PH per SD of ME; KM by median only as a display. No data-driven cutoff search.
4. **Single-gene comparators:** repeat 1–3 for TACSTD2 and for CLDN4 alone. The module must beat, or at least be reported beside, the seeds; a module that is just TACSTD2 is not a new finding.
5. **Multiple testing:** BH-FDR within the prespecified family (e.g. all modules × one primary endpoint, or primary ME × a small locked signature set). Do not FDR-pool exploratory splits with the primary.

1. **二分类疗效：** Mann–Whitney U、秩双列 *r*、每臂确切 *n*。锁定 ME 的 ROC AUC 仅作描述，不在同一队列上调阈值。
2. **连续免疫评分：** Spearman *ρ*（成对完整）。签名定义与项目库一致（现存基因的均 *z*；记录 `used/total`）。
3. **生存：** 按 ME 每 1 SD 的 Cox PH；中位数 KM 仅作展示，不做数据驱动切点搜索。
4. **单基因对照：** 对 TACSTD2、CLDN4 单独重复 1–3。模块须优于或至少与种子并列报告；若模块只是 TACSTD2 的化身，则不是新发现。
5. **多重检验：** 在预设家族内做 BH-FDR。不要把探索性拆分与主分析混在一个 FDR 池里。

Covariate-adjusted models (purity, histology, PD-L1) only when Section 2.2 logistic/Cox *n* is met. Otherwise report stratified descriptives (e.g. ME by histology) without a 6-parameter GLM on 16 patients.

仅当第 2.2 节的 logistic/Cox 样本量满足时才做协变量调整。否则只做分层描述，不在 16 例上拟合 6 参数 GLM。

### 6.3 Locking and projection / 锁定与投影

A lock file (JSON or YAML) must contain: discovery accession(s), *n*, *β*, network type, module color/label, gene symbols + kME, ME sign flip, scoring formula, software versions, and a content hash. Application cohorts may:

锁定文件（JSON 或 YAML）须包含：发现队列登录号、*n*、*β*、网络类型、模块颜色/标签、基因符号 + kME、ME 符号翻转、评分公式、软件版本和内容哈希。应用队列可以：

- compute a **fixed weighted score** \(\sum_g w_g z_g\) with \(w_g\) = discovery kME (or 1/|*M*| for an unweighted mean *z*), using only genes present; or
- compute the first principal component of the locked genes (orientation matched to TACSTD2); or
- use hdWGCNA `ProjectModules` / WGCNA `moduleEigengenes` on the locked genes.

- 用发现集 kME 作权重计算**固定加权分** \(\sum_g w_g z_g\)（或等权均 *z*），只用现存基因；或
- 对锁定基因做第一主成分（方向与 TACSTD2 对齐）；或
- 对锁定基因调用 hdWGCNA `ProjectModules` / WGCNA `moduleEigengenes`。

Do not re-estimate *β*, re-cut the dendrogram, or drop genes because they look weak versus ICI in the application set.

不得在应用集上重估 *β*、重切树、或因相对 ICI “看起来弱”而删基因。

---

## 7. Fallback when gates fail / 门槛失败时的回退

If *n* < 15, metacell counts fail, a seed is missing, or scale-free topology is uninterpretable **and** the FAQ power would be used on a cohort that also has a strong unmodeled batch, **do not build a network**.

若 *n* < 15、metacell 不足、种子缺失，或无标度拓扑无法解释且 FAQ 幂次还要建在未处理的强批次上，则**不要建网**。

Use a **locked, outcome-blind junction score**:

改用**锁定、与结局无关的连接评分**：

1. **Core pair (always):** TACSTD2, CLDN4. Report their Spearman *ρ* and the mean *z* of the pair (`used/total` = 2/2 or 1/2).
2. **Prespecified panel:** protein-coding genes in a frozen ontology snapshot — e.g. Reactome “Tight junction interactions” and/or GO “tight junction” (GO:0070160) **plus** TACSTD2 if it is not already in those lists. Download the list once, store it under `methods/wgcna/` (gene symbols + ontology release), and never edit it after ICI tests begin.
3. Score = mean *z* of genes present in that cohort (project convention). No weights learned from ICI labels.

1. **核心对（始终）：** TACSTD2、CLDN4。报告 Spearman *ρ* 及二者均 *z*（`used/total` = 2/2 或 1/2）。
2. **预设 panel：** 冻结本体快照中的蛋白编码基因，例如 Reactome “Tight junction interactions” 和/或 GO “tight junction”（GO:0070160），若 TACSTD2 不在其中则**并入**。下载一次，存入 `methods/wgcna/`（符号 + 本体版本），ICI 检验开始后不得再改。
3. 评分 = 该队列现存基因的均 *z*。权重不得从 ICI 标签学习。

This fallback is the correct primary analysis for typical *n* = 16 ICI bulks. WGCNA on those cohorts is optional exploration and must be labeled as such.

对典型 n=16 的 ICI bulk，该回退才是正确的主分析。在这些队列上做 WGCNA 只是可选探索，必须标明。

---

## 8. Diagnostics (methods, not results) / 诊断（方法，不是结果）

Produce and archive, for every network that was *allowed* to run:

对每一个*被允许*运行的网络，生成并归档：

- sample dendrogram (colored by batch, histology, response — response coloring is QC, not the primary test);
- scale-free *R*² and mean connectivity versus *β*;
- module size table, including grey fraction;
- kME table with TACSTD2 and CLDN4 highlighted;
- MM vs gene-significance scatter for the primary endpoint *only if* Section 2.2 allows the test;
- donor-level ME boxplots by response with *n* annotated;
- metacell UMAP (hdWGCNA) colored by sample and cell type;
- a one-row gate report from `sample_size_gates.R`.

- 样本树（按批次、组织学、疗效着色——疗效着色是质控，不是主检验）；
- *R*² 与平均连接度对 *β* 的曲线；
- 模块大小表（含灰色比例）；
- 突出 TACSTD2、CLDN4 的 kME 表；
- 仅当第 2.2 节允许时，主终点的 MM–基因显著性散点；
- 按疗效分组的供体 ME 箱线，标注 *n*；
- hdWGCNA metacell UMAP（样本、细胞类型着色）；
- `sample_size_gates.R` 的一行门槛报告。

If diagnostics show a two-block sample tree (batch or a global driver) and *R*² never reaches 0.80, stop and fix batch (`methods/batch/`) or fall back to Section 7. Do not “fix” topology by trait-based gene filters.

若诊断显示样本树两块分裂（批次或全局驱动）且 *R*² 达不到 0.80，应停止并处理批次，或回退第 7 节。禁止用基于性状的基因过滤来“修好”拓扑。

---

## 9. Outputs this methods slice should produce / 本方法切片应产出的文件

When an analysis *run* uses this playbook, write under that run’s `results/` (not into this methods folder):

当某次分析*运行*采用本手册时，写入该次运行的 `results/`（不要写进本 methods 目录）：

| Artifact | Content |
| --- | --- |
| `wgcna_gate_report.tsv` | One row per cohort: *n*, events, metacells, decision (`forbid` / `explore` / `allow`), reason |
| `module_membership.tsv` | Gene, module, kME, is_seed |
| `junction_lock.json` | Section 6.3 lock (discovery only) |
| `me_ici_associations.tsv` | Cohort, endpoint, *n*, *n* per arm, stat, *p*, FDR, note `exploratory`/`primary` |
| `coexpression_TACSTD2_CLDN4.tsv` | Pairwise *ρ* (always, including fallback) |

No fabricated statistics. If a test was forbidden, the association table gets a row with `stat = NA` and the gate reason — not a silent omission.

禁止编造统计量。若检验被禁止，关联表仍写一行 `stat = NA` 及门槛原因，而不是 silently 省略。

---

## 10. Software / 软件

| Tool | Role | Reference |
| --- | --- | --- |
| WGCNA (R) | Bulk networks, ME, kME, bicor | Langfelder & Horvath, *BMC Bioinformatics* 2008; FAQ (2017 power table) |
| hdWGCNA (R) | Metacells, sc/spatial networks, projection | Morabito et al., *Cell Reports Methods* 2023; https://smorabit.github.io/hdWGCNA/ |
| Seurat v5 | scRNA container | Hao et al.; see `methods/scrna/` |
| DESeq2 VST or log-CPM / log-TPM | Bulk transform | Outcome-blind |

Record exact package versions in the lock file. Do not mix unsigned tutorial defaults with this signed SOP without saying so.

在锁定文件中记录确切版本。不得在未说明的情况下把 unsigned 教程默认值与本 signed SOP 混用。

---

## 11. What this playbook does not do / 本手册明确不做的事

- Does not download or analyze data by itself; it is not a results writeup.
- Does not invent accessions, ICI labels, or *p*-values.
- Does not treat CPTAC treatment-naive protein as an ICI-response cohort.
- Does not treat cells, metacells, or spots as independent ICI outcomes.
- Does not train an ICI-supervised junction signature on *n* < 40 (and even at *n* ≥ 40, prediction belongs in `methods/ml/` with nested CV).
- Does not use EGA/dbGaP matrices without access.
- Does not run WGCNA on targeted panels that lack TACSTD2 and CLDN4.

- 本身不下载、不分析数据；不是结果文稿。
- 不编造登录号、ICI 标签或 *p* 值。
- 不把 CPTAC 初治蛋白当作 ICI 疗效队列。
- 不把细胞、metacell 或 spot 当作独立 ICI 结局。
- 不在 *n* < 40 上训练 ICI 监督的连接签名（即便 *n* ≥ 40，预测也应走 `methods/ml/` 的嵌套 CV）。
- 无权限时不使用 EGA/dbGaP 矩阵。
- 不在缺少 TACSTD2 与 CLDN4 的靶向 panel 上跑 WGCNA。
