# 高级去卷积方法手册 / Advanced Deconvolution Playbook

**Scope: METHODS ONLY — 仅方法学**

> 中文 — 本文件只描述**方法、协议与决策规则**。此处**不包含任何数据、结果、效应量、p 值或结论**。文中出现的所有数字均为**协议参数**（阈值、迭代次数、随机种子等），不是任何数据分析的产出。任何结果必须写入独立的 `results/` 目录，并回引本文件的章节编号。
>
> EN — This file specifies **methods, protocols, and decision rules only**. It contains **no data, no results, no effect sizes, no p-values, and no conclusions**. Every number below is a **protocol parameter** (threshold, iteration count, random seed), never an output of a data analysis. Results must live in a separate `results/` tree and cite the section numbers here.

| | |
|---|---|
| Path | `methods/deconv_advanced/playbook.md` |
| Coverage window | 2024–2026 method landscape |
| Methods in scope | BayesPrism / InstaPrism, CIBERSORTx, MuSiC / MuSiC2, EcoTyper, InstantDL (scope-corrected, see §1.5) |
| Applied question | Correlating `TACSTD2` (TROP2) with CD8 T-cell and TLS estimates in bulk ICI cohorts (§7) |
| Companion protocol | BayesPrism / EcoTyper vs **both** `TACSTD2` and `CLDN4`: `bayesprism_ecotyper_vs_tacstd2_cldn4.md` (§13) |
| Comparison protocol | §8 (methods) and §B6 (BayesPrism vs EcoTyper, two seeds) |
| Languages | 中文 + English, interleaved per section; code/tables shared |

---

## 目录 / Table of Contents

- [§0 前置约定 / Conventions and Non-Goals](#0-前置约定--conventions-and-non-goals)
- [§1 方法全景 / Method Landscape](#1-方法全景--method-landscape-20242026)
- [§2 输入数据与协调 / Inputs and Harmonization](#2-输入数据与协调--inputs-and-harmonization)
- [§3 参考构建 / Reference Construction](#3-参考构建--reference-construction)
- [§4 各方法运行规程 / Per-Method Run Procedures](#4-各方法运行规程--per-method-run-procedures)
- [§5 输出统一化 / Output Harmonization](#5-输出统一化--output-harmonization)
- [§6 TLS 量化模块 / TLS Quantification Module](#6-tls-量化模块--tls-quantification-module)
- [§7 TACSTD2 × CD8/TLS 关联协议 / The TACSTD2 Correlation Protocol](#7-tacstd2--cd8tls-关联协议--the-tacstd2-correlation-protocol)
- [§8 方法比较协议 / Comparison Protocol](#8-方法比较协议--comparison-protocol)
- [§9 可复现性与报告 / Reproducibility and Reporting](#9-可复现性与报告--reproducibility-and-reporting)
- [§10 目录结构 / Directory Layout](#10-目录结构--directory-layout)
- [§11 陷阱清单 / Pitfall Register](#11-陷阱清单--pitfall-register)
- [§12 参考文献 / References](#12-参考文献--references)
- [§13 BayesPrism / EcoTyper × TACSTD2 / CLDN4](#13-bayesprism--ecotyper--tacstd2--cldn4)

---

## §0 前置约定 / Conventions and Non-Goals

中文 — 本手册服务于一个具体的分析形态：**手上只有 bulk RNA-seq（多为 ICI 治疗队列，常为 FFPE、样本量中等），需要推断细胞组成，再把某个上皮基因（`TACSTD2`）与免疫读数（CD8、TLS）关联起来**。这个形态有三个结构性难点，全书围绕它们展开：

1. **细胞比例是成分数据（compositional）**：所有比例之和为 1，任一细胞类型的比例天然与肿瘤纯度反向绑定。未经处理的 Pearson/Spearman 会把"闭合效应"读成生物学信号。见 §5.3、§7.5。
2. **`TACSTD2` 是上皮/肿瘤区室基因**：bulk 中它的表达同时编码"每个肿瘤细胞表达多少"和"样本里有多少肿瘤细胞"。不拆开这两者，任何与免疫比例的相关都无法解释。见 §7.2、§7.4。
3. **方法之间不一致是常态**：2024–2026 的多项基准研究一致显示不存在普适最优方法，且 T 细胞亚群之间的溢出（spillover）是共性弱点。因此**主结论必须在方法间稳定**，这一点被写成正式终点，见 §8.4。

非目标 / Non-goals：本手册不涉及空间转录组的原生分析（仅在 §8.1 作为金标准来源出现）、不涉及 scRNA-seq 自身的聚类注释流程（仅作为参考输入，见 §3）、不给出任何队列层面的生物学解释。

EN — This playbook targets one concrete analysis shape: **you have bulk RNA-seq (typically ICI-treated cohorts, often FFPE, moderate n), you need cell composition, and you then want to relate an epithelial gene (`TACSTD2`) to immune readouts (CD8, TLS)**. Three structural difficulties organize everything that follows:

1. **Cell fractions are compositional data.** They sum to one, so every immune fraction is mechanically anti-correlated with tumor purity. An unadjusted Pearson/Spearman coefficient reads the closure effect as biology. See §5.3 and §7.5.
2. **`TACSTD2` is an epithelial/malignant-compartment gene.** Its bulk expression confounds "how much each tumor cell expresses" with "how many tumor cells are present." Until those are separated, no correlation with an immune fraction is interpretable. See §7.2 and §7.4.
3. **Between-method disagreement is the normal case.** Benchmarks published across 2024–2026 consistently report that no method dominates and that spillover among T-cell subsets is a shared weakness. Therefore **conclusion stability across methods is a formal endpoint**, not a robustness afterthought. See §8.4.

Non-goals: native spatial-transcriptomics analysis (appears only as a ground-truth source in §8.1), scRNA-seq clustering/annotation itself (consumed as a reference input, §3), and any cohort-level biological interpretation.

### §0.1 术语与记号 / Notation

| Symbol | Meaning |
|---|---|
| $B_{g \times n}$ | bulk expression matrix, genes × samples |
| $R_{g \times k}$ | reference signature / profile matrix, genes × cell types |
| $\theta_{k \times n}$ | inferred cell-type fraction matrix (columns sum to 1) |
| $Z_{g \times k \times n}$ | cell-type-specific expression posterior (BayesPrism) or imputed GEP (CIBERSORTx HiRes) |
| $\pi_n$ | tumor purity / malignant fraction of sample $n$ |
| $E^{\text{mal}}_{g,n}$ | expression of gene $g$ **within the malignant compartment** of sample $n$ |
| $\mathrm{clr}(\cdot)$ | centered log-ratio transform |
| CS / CE | EcoTyper cell state / carcinoma ecotype |

### §0.2 强制性规则 / Hard Rules

中文 — 以下规则在本项目中不可协商，违反即视为分析无效：

- **R0-1** 任何比例读数在进入线性模型前必须完成成分变换（§5.3），或改用 beta/Dirichlet 回归。
- **R0-2** 任何 `TACSTD2` 与免疫读数的相关必须同时报告"未校正"和"纯度校正"两个版本（§7.3）。只报告其中一个视为选择性报告。
- **R0-3** 主分析方法必须在看到任何关联结果**之前**冻结（§8.0 预注册）。
- **R0-4** 每个关联分析必须附带阳性对照与阴性对照（§7.8）。
- **R0-5** 跨队列合并必须使用随机效应元分析或含队列随机截距的混合模型，禁止简单池化（§7.6）。

EN — Non-negotiable in this project; violating any of them invalidates the analysis:

- **R0-1** Any fraction entering a linear model must first be compositionally transformed (§5.3), or else modeled with beta/Dirichlet regression.
- **R0-2** Every `TACSTD2`–immune association must report both the unadjusted and the purity-adjusted estimate (§7.3). Reporting only one counts as selective reporting.
- **R0-3** The primary method must be frozen **before** any association result is seen (§8.0 pre-registration).
- **R0-4** Every association analysis ships with its positive and negative controls (§7.8).
- **R0-5** Cross-cohort aggregation uses random-effects meta-analysis or a mixed model with a cohort random intercept; naive pooling is prohibited (§7.6).

---

## §1 方法全景 / Method Landscape (2024–2026)

中文 — 下面五张"方法卡"给出每个工具的模型假设、输入输出契约、以及在本项目中的**指定角色**。角色分配比工具本身重要：BayesPrism/InstaPrism 与 CIBERSORTx HiRes 承担区室解析（§7.4 的核心），MuSiC/MuSiC2 承担快速比例估计与条件失配控制，EcoTyper 承担细胞状态层面的表型，InstantDL 承担影像侧正交验证。

EN — Five method cards below give each tool's modeling assumptions, input/output contract, and its **assigned role** in this project. The role assignment matters more than the tool: BayesPrism/InstaPrism and CIBERSORTx HiRes carry compartment resolution (the core of §7.4), MuSiC/MuSiC2 carry fast fraction estimation and condition-mismatch control, EcoTyper carries cell-state-level phenotyping, and InstantDL carries orthogonal image-side validation.

### §1.1 BayesPrism / InstaPrism

| Field | Specification |
|---|---|
| Model | Bayesian generative model; joint posterior over cell-type fractions $\theta$ and cell-type-specific expression $Z$, with a scRNA-seq reference as prior |
| Key property | Robust to reference–bulk misspecification (within-cell-type heterogeneity, platform/technical shift), as reported by multiple independent 2023–2026 benchmarks |
| Malignant handling | Explicit `key` argument marks the malignant cell type; its profile is updated per sample rather than held fixed — this is why it can serve §7.4 |
| Inference | BayesPrism: Gibbs sampling (slow, hours at cohort scale). InstaPrism: derandomized fixed-point algorithm, effectively equivalent output, seconds-scale |
| Input units | Raw counts for both reference and bulk (do **not** feed TPM) |
| Output | $\theta$ (type and state level); $Z$ via `get.exp()`; InstaPrism compresses $Z$ into a 2D scaling matrix |
| Role here | **Primary candidate** for both fractions and compartment-resolved `TACSTD2` |
| Cost | InstaPrism makes cohort-scale + bootstrap + sensitivity sweeps feasible; plain BayesPrism generally does not |

中文 — 实践要点：(a) InstaPrism 提供预编译的癌种特异参考集，但**若你的癌种与参考癌种不匹配，自建参考优先**（§3）；(b) `outlier.cut` / `outlier.fraction` 用于剔除在 bulk 中异常高表达、会主导似然的基因（线粒体、核糖体、血红蛋白），必须在配置中显式写死而非依赖默认；(c) `update = TRUE`（参考更新模式）与不更新模式会给出不同的 $Z$，在 §8.3 中作为一个稳健性轴。

EN — Practical points: (a) InstaPrism ships curated per-cancer-type references, but **build your own when your tumor type does not match** (§3); (b) `outlier.cut` / `outlier.fraction` remove genes whose bulk expression dominates the likelihood (mitochondrial, ribosomal, hemoglobin) and must be pinned explicitly in config rather than left to defaults; (c) reference-update mode (`update = TRUE`) versus no-update yields different $Z$, and this becomes one robustness axis in §8.3.

### §1.2 CIBERSORTx

| Field | Specification |
|---|---|
| Model | $\nu$-SVR deconvolution against a signature matrix, with two batch-correction modes |
| Mode: Fractions | Relative fractions by default; absolute mode rescales by median signature expression |
| Mode: High-Resolution (HiRes) | Imputes gene-level, **cell-type-specific expression per sample** — the input EcoTyper consumes and a second route to §7.4 |
| Mode: Group GEP | Cell-type expression per group rather than per sample; use when per-sample HiRes is too sparse |
| Batch correction | **B-mode** for bulk-vs-signature platform mismatch (e.g. array signature, RNA-seq bulk); **S-mode** for scRNA-seq-derived signatures (droplet dropout) |
| Input units | Non-log, linear space; TPM/CPM for mixture; counts or CPM for the single-cell reference sample file |
| Access | Web portal or Docker image with a username/token; license is restricted (non-commercial by default) — record the license terms in `env/LICENSES.md` |
| Role here | **Required** if EcoTyper recovery is run; secondary compartment-resolution route; the historical reference method that reviewers expect |

中文 — 实践要点：(a) **B-mode 与 S-mode 不可同时开**，且选错模式是该工具最常见的误用来源；(b) HiRes 输出对低丰度细胞类型会产生大量缺失/不可估计条目，`TACSTD2` 在恶性区室通常丰度足够，但 CD8 区室的低表达基因会大面积缺失——这决定了 §7.4 只在恶性区室做 `TACSTD2`；(c) 相对比例与绝对比例不可混用于同一模型，必须在配置中固定一种并在 §8 中作为敏感性轴。

EN — Practical points: (a) **B-mode and S-mode are mutually exclusive**, and picking the wrong one is this tool's single most common misuse; (b) HiRes output carries many missing/non-estimable entries for low-abundance cell types — `TACSTD2` is typically well-powered in the malignant compartment, but low-expression genes in the CD8 compartment drop out broadly, which is precisely why §7.4 restricts `TACSTD2` imputation to the malignant compartment; (c) relative and absolute fractions must never be mixed inside one model; pin one in config and treat the other as a sensitivity axis in §8.

### §1.3 MuSiC / MuSiC2

| Field | Specification |
|---|---|
| Model (MuSiC) | Weighted non-negative least squares; weights favor genes with cross-subject and cross-cell consistency; multi-subject reference is required by design |
| Collinearity | Tree-guided recursive estimation: cluster related cell types, estimate cluster proportions, then recurse within cluster using low-within-cluster-variance genes |
| Model (MuSiC2) | Iterative two-step procedure for **condition mismatch** between reference and bulk: deconvolve → detect and remove cell-type-specific DE genes → re-deconvolve, until convergence |
| Why MuSiC2 matters here | ICI cohorts are treated/diseased tissue; most public scRNA-seq references are treatment-naive or from different disease states. MuSiC2 directly targets that mismatch |
| Input units | Counts; reference as `SingleCellExperiment` with `clusters` and `samples` columns |
| Output | $\theta$ only (no per-sample cell-type expression posterior) |
| Role here | Fast fraction estimator; the condition-mismatch control arm; **not** usable for §7.4 |

中文 — 实践要点：(a) MuSiC 的多受试者加权在参考只有 1–2 个供体时失效，参考至少需要 ≥5 个供体（§3.2）；(b) 树引导递归的分组必须**预先按免疫学层级指定**（T/NK 一支、髓系一支、B/浆细胞一支、上皮/恶性一支），不要用数据驱动的树，否则不同队列会得到不同的树结构而破坏可比性；(c) MuSiC2 的 `control` / `case` 划分在 ICI 场景中应对应"参考条件"与"待测条件"，而不是"响应/不响应"——把响应状态放进去会把结局信息泄漏进比例估计。

EN — Practical points: (a) MuSiC's multi-subject weighting degenerates when the reference has 1–2 donors; require ≥5 donors (§3.2); (b) fix the tree-guided grouping **a priori along immunological lineage** (T/NK, myeloid, B/plasma, epithelial/malignant) rather than learning it from data, otherwise different cohorts get different tree structures and lose comparability; (c) in MuSiC2 the `control`/`case` split must encode "reference condition" versus "query condition," **not** responder versus non-responder — putting outcome status into the split leaks outcome information into the fraction estimates.

### §1.4 EcoTyper

| Field | Specification |
|---|---|
| Model | NMF over cell-type-purified expression (typically CIBERSORTx HiRes output) to define **cell states (CS)** within each cell type, then co-occurrence modeling across states to define **ecotypes / cellular communities (CE)** |
| Two workflows | **Discovery** (derive new CS/CE from your own bulk, scRNA-seq, or sorted profiles) and **Recovery** (project previously defined CS/CE onto new bulk, scRNA-seq, or spatial data) |
| Pre-built models | Pan-carcinoma CS/CE set (CE1–CE10) and a DLBCL-specific set |
| Critical preprocessing | Genes log2-transformed and **scaled to unit variance within each tumor type / dataset** before recovery; skipping this silently invalidates recovery |
| Significance | Recovery is assessed per cell state with a permutation test; a z-score threshold (commonly 1.65) declares a state recovered |
| Output | Per-sample CS abundance and CS assignment; per-sample CE abundance |
| Role here | Cell-**state** resolution layer: distinguishes CD8 states and B/plasma states that a fraction-only view collapses; supplies a multicellular-community readout to sit alongside the TLS score (§6) |

中文 — 实践要点：(a) **Recovery 优先于 Discovery**。在中等样本量的 ICI 队列上做 de novo discovery 会得到队列特异的状态定义，跨队列无法比较；只有当样本量充足且明确要建立新分型时才做 discovery，并必须在独立队列上 recovery 验证。(b) EcoTyper 的状态标签（如 CD8 的第 N 个状态）**在不同 discovery 运行之间没有稳定语义**，引用时必须绑定到具体模型版本。(c) 未通过 permutation 阈值的状态必须整体从下游剔除，而不是当作"低丰度"继续使用。

EN — Practical points: (a) **Prefer recovery over discovery.** De novo discovery on a moderate-sized ICI cohort yields cohort-specific state definitions that cannot be compared across cohorts; run discovery only with ample n and an explicit goal of defining new subtypes, and then validate by recovery in an independent cohort. (b) EcoTyper state labels (e.g. "CD8 state N") **carry no stable meaning across discovery runs**; always bind a label to a specific model version when citing it. (c) States failing the permutation threshold must be dropped wholesale from downstream analysis, not carried forward as "low abundance."

### §1.5 InstantDL — 范围更正 / Scope Correction

> 中文 — **重要**：`InstantDL` **不是** bulk 转录组去卷积工具。它是一个生物医学**影像**深度学习流水线（语义分割、实例分割、逐像素回归、图像分类），由 Waibel、Shetab Boushehri 与 Marr 发表于 *BMC Bioinformatics* (2021)。把它当作 BayesPrism/CIBERSORTx 的同类进行"细胞比例"比较是范畴错误，本手册不这样做。
>
> EN — **Important**: `InstantDL` is **not** a bulk-transcriptome deconvolution tool. It is a biomedical **imaging** deep-learning pipeline (semantic segmentation, instance segmentation, pixel-wise regression, image classification) from Waibel, Shetab Boushehri, and Marr, *BMC Bioinformatics* (2021). Treating it as a peer of BayesPrism/CIBERSORTx in a cell-fraction comparison is a category error, and this playbook does not do that.

中文 — 因此 InstantDL 在本项目中被重新指派为**两个合法角色**：

**角色 A — 影像侧正交金标准（推荐，见 §8.1）。** 对配套的 H&E 或多重免疫荧光（mIF）切片，用 InstantDL 的实例分割 + 分类任务产出 (i) CD8⁺ 细胞密度，(ii) TLS 区域的检测与面积占比。这些是**与转录组完全正交**的测量，用来校准去卷积输出——正是当前基准研究最缺的东西。其内置的不确定性估计（分类头带 dropout）可直接给出每张切片的置信权重，用于 §8.2 的加权指标。

**角色 B — 深度学习去卷积的位置由替代工具填补。** 如果最初的意图是"纳入一个深度学习类去卷积方法"，那么正确的候选是 **Scaden**（模拟数据训练的前馈集成）、**TAPE**（自编码器，可输出细胞类型特异表达）或 **DISSECT** 一类半监督方法，而不是 InstantDL。2024–2026 的基准中，Scaden 在**溢出（spillover）**指标上反复表现为最低之一，这正是 CD8 与 CD4/NK/Treg 混淆的关键指标，因此若要加一个 DL 臂，Scaden 是与本项目问题最匹配的选择。

EN — InstantDL is therefore reassigned to **two legitimate roles**:

**Role A — orthogonal image-side ground truth (recommended, see §8.1).** On matched H&E or multiplex immunofluorescence (mIF) slides, use InstantDL's instance-segmentation and classification tasks to produce (i) CD8⁺ cell density and (ii) TLS region detection and area fraction. These are measurements **fully orthogonal to the transcriptome**, used to calibrate deconvolution output — exactly what current benchmarks most lack. Its built-in uncertainty estimation (dropout in the classification head) yields a per-slide confidence weight usable for the weighted metrics in §8.2.

**Role B — the deep-learning deconvolution slot is filled by substitutes.** If the original intent was "include a deep-learning deconvolution method," the correct candidates are **Scaden** (feed-forward ensemble trained on simulated mixtures), **TAPE** (autoencoder, also emits cell-type-specific expression), or a semi-supervised method such as **DISSECT** — not InstantDL. Across 2024–2026 benchmarks Scaden repeatedly ranks among the lowest on **spillover**, the metric that governs CD8-versus-CD4/NK/Treg confusion, so Scaden is the DL arm best matched to this project's question.

> 记录要求 / Documentation requirement: 在 `results/` 的任何输出中提到 InstantDL 时，必须附上本节的范围说明。/ Any mention of InstantDL in `results/` must carry this scope note.

### §1.6 方法—角色矩阵 / Method–Role Matrix

| Method | Fractions $\theta$ | Cell-type-specific expression $Z$ | Cell states | Handles condition mismatch | Assigned role |
|---|---|---|---|---|---|
| InstaPrism / BayesPrism | yes | yes (native posterior) | via cell-state labels | partially (robust prior) | Primary fractions + primary §7.4 route |
| CIBERSORTx | yes (rel/abs) | yes (HiRes imputation) | no (feeds EcoTyper) | via B-/S-mode batch correction | Secondary fractions + EcoTyper input + secondary §7.4 route |
| MuSiC | yes | no | no | no | Fast comparator arm |
| MuSiC2 | yes | no | no | **yes (explicit target)** | Condition-mismatch control arm |
| EcoTyper | state abundance | consumes $Z$ | **yes (CS/CE)** | n/a | Cell-state + ecotype layer |
| Scaden (DL substitute) | yes | no | no | via simulation design | Optional DL arm, low-spillover comparator |
| InstantDL | **no** | **no** | no | n/a | Image-side orthogonal ground truth (§8.1) |

---

## §2 输入数据与协调 / Inputs and Harmonization

### §2.1 Bulk ICI 队列清单 / Bulk ICI Cohort Register

中文 — 建立 `data/cohorts/registry.tsv`，每个队列一行，字段固定如下。**登记必须在任何分析前完成**，因为队列的可及性和平台差异决定了后续统计模型的形式（§7.6）。

EN — Maintain `data/cohorts/registry.tsv`, one row per cohort, with the fixed schema below. **Registration must be complete before any analysis**, because accessibility and platform differences determine the form of the statistical model in §7.6.

| Column | Content |
|---|---|
| `cohort_id` | stable short key used in all filenames |
| `tumor_type` | OncoTree code |
| `agent_class` | anti-PD-1 / anti-PD-L1 / anti-CTLA-4 / combination |
| `timepoint` | pre-treatment / on-treatment / post-progression |
| `assay` | RNA-seq (poly-A / rRNA-depleted / exome-capture) or microarray or NanoString |
| `fixation` | FFPE / fresh-frozen |
| `quantification` | raw counts / TPM / FPKM / normalized intensities |
| `n_samples` | integer |
| `has_response` | RECIST availability |
| `has_survival` | PFS / OS availability |
| `has_purity_dna` | DNA-derived purity availability (WES/SNP) |
| `has_matched_imaging` | H&E / mIF availability (gates §8.1 Role A) |
| `access` | open / dbGaP / EGA + accession |

中文 — 常用的公开或受控访问 bulk ICI 队列（作为**候选清单**，不预设纳入）：泌尿上皮癌 atezolizumab（IMvigor210，随 `IMvigor210CoreBiologies` R 包分发）；黑色素瘤抗 PD-1（Riaz、Hugo、Liu、Gide）与抗 CTLA-4（Van Allen）；ccRCC（Braun/CheckMate 009-010-025、Miao）；胃癌 pembrolizumab（Kim）；NSCLC（POPLAR/OAK，受控访问）；GBM（Zhao）。整合层可用 CRI iAtlas、TIDE、ICBatlas 等资源，但**必须记录其内部已做的归一化**，否则会与 §2.3 冲突。

EN — Commonly used public or controlled-access bulk ICI cohorts (a **candidate list**, not a pre-committed inclusion set): urothelial atezolizumab (IMvigor210, distributed with the `IMvigor210CoreBiologies` R package); melanoma anti-PD-1 (Riaz, Hugo, Liu, Gide) and anti-CTLA-4 (Van Allen); ccRCC (Braun / CheckMate 009-010-025, Miao); gastric pembrolizumab (Kim); NSCLC (POPLAR/OAK, controlled access); GBM (Zhao). Aggregators such as CRI iAtlas, TIDE, and ICBatlas may be used, but **you must record the normalization they already applied**, or it will collide with §2.3.

### §2.2 单位契约 / Units Contract

中文 — 单位错误是本流程中最高频的静默失败。下表是硬性契约，写进 `config/units.yaml` 并在运行前由断言检查。

EN — Unit errors are the highest-frequency silent failure here. The table below is a hard contract; encode it in `config/units.yaml` and assert it before every run.

| Consumer | Bulk input | Reference input | Log? |
|---|---|---|---|
| BayesPrism / InstaPrism | raw counts | raw counts (cell × gene) | no |
| CIBERSORTx fractions | TPM or CPM, linear | scRNA-seq counts or CPM ("refsample" format) | no |
| CIBERSORTx HiRes | same mixture file as fractions run | signature matrix from the fractions run | no |
| MuSiC / MuSiC2 | raw counts | counts in `SingleCellExperiment` | no |
| EcoTyper recovery | TPM/FPKM or array intensities | pre-built model | **log2, then unit-variance scaling per dataset** |
| ssGSEA / GSVA (TLS score) | log2(TPM+1) | gene sets | yes |
| singscore (TLS score) | rank-based, any monotone transform | gene sets | irrelevant (ranks) |

```yaml
# config/units.yaml (excerpt)
bulk:
  counts_matrix: data/processed/{cohort_id}/counts.tsv        # integers, no log
  tpm_matrix:    data/processed/{cohort_id}/tpm.tsv           # linear, colSums == 1e6
  log_matrix:    data/processed/{cohort_id}/log2tpm1.tsv      # log2(tpm + 1)
assertions:
  counts_are_integer: true
  tpm_colsum_tolerance: 1.0        # |colSum - 1e6| must be below this
  no_negative_values: true
  gene_id_space: "ENSG_unversioned"
```

### §2.3 基因标识与归一化 / Gene Identifiers and Normalization

中文 — 规程：

1. **统一到无版本号的 Ensembl gene ID**，在最后一步再映射到 HGNC symbol 供人类阅读。symbol 在不同年份的注释间漂移，直接用 symbol 做跨队列匹配会静默丢基因。
2. **多对一映射**：同一 symbol 对应多个 ENSG 时，保留表达量总和最高者，并把决策写入 `data/processed/{cohort_id}/gene_map_decisions.tsv`。
3. **`TACSTD2` 特别注意**：该基因的编码区位于单个外显子（intronless）。后果有二：(a) 在**单核**（snRNA-seq）参考中，intronless 基因的 pre-mRNA/内含子读段策略与多外显子基因不同，可能系统性偏移其参考谱；(b) 在有基因组 DNA 残留的文库中，单外显子基因更容易被 DNA 污染抬高。因此若参考来自 snRNA-seq，必须在 §8.3 中把"参考来源 = scRNA vs snRNA"作为一条敏感性轴，并在注释中确认外显子结构。
4. **不做跨队列的 batch 校正后再去卷积**。去卷积在每个队列内独立运行；批次效应在**下游统计层**通过队列随机效应处理（§7.6）。对 bulk 矩阵先做 ComBat 再去卷积，会破坏去卷积方法所假设的线性混合尺度。
5. **FFPE 队列**：额外记录 3′ 偏倚指标与 rRNA 残留比例，作为 §7.6 模型的候选协变量与 §8.3 的分层轴。

EN — Procedure:

1. **Normalize to unversioned Ensembl gene IDs**, mapping to HGNC symbols only at the final human-readable step. Symbols drift between annotation releases, and matching cohorts on symbols silently drops genes.
2. **Many-to-one mappings**: when one symbol maps to several ENSGs, keep the highest-total-expression entry and log the decision to `data/processed/{cohort_id}/gene_map_decisions.tsv`.
3. **Special handling for `TACSTD2`**: its coding sequence sits in a single exon (intronless). Two consequences: (a) in **single-nucleus** (snRNA-seq) references, intron/pre-mRNA counting strategies treat intronless genes differently from multi-exon genes and can systematically shift their reference profile; (b) in libraries with residual genomic DNA, single-exon genes are more easily inflated by DNA contamination. If the reference derives from snRNA-seq, carry "reference source = scRNA vs snRNA" as a sensitivity axis in §8.3, and confirm the exon structure in your annotation build.
4. **Do not batch-correct across cohorts before deconvolution.** Run deconvolution independently within each cohort and handle batch at the **downstream statistical layer** via cohort random effects (§7.6). Running ComBat on the bulk matrix first destroys the linear mixture scale that deconvolution methods assume.
5. **FFPE cohorts**: additionally record 3′ bias metrics and residual rRNA fraction as candidate covariates for §7.6 and as a stratification axis in §8.3.

### §2.4 样本级 QC 门 / Sample-Level QC Gates

中文 — 在去卷积前应用，失败样本记录原因后剔除，剔除清单进入 `data/processed/{cohort_id}/excluded.tsv`：

EN — Applied before deconvolution; failures are removed with a recorded reason into `data/processed/{cohort_id}/excluded.tsv`:

| Gate | Criterion (pin exact values in config) |
|---|---|
| Library size | below a pre-specified floor of mapped reads |
| Detected genes | below a pre-specified floor |
| rRNA / mitochondrial fraction | above a pre-specified ceiling |
| Hemoglobin fraction | above a pre-specified ceiling (blood contamination distorts immune fractions) |
| Sex-check | `XIST` / `RPS4Y1` consistency with recorded sex (**required**, because one TLS signature contains a Y-linked gene, §6.2) |
| Duplicate / relatedness | expression-correlation outliers flagged and adjudicated |
| Purity availability | recorded, not a gate; drives the §7.3 adjustment strategy |

---

## §3 参考构建 / Reference Construction

中文 — 参考决定了去卷积能回答什么问题。2024–2026 基准的一致教训是：**参考构建的贡献常常大于算法选择的贡献**。因此参考本身被当作一个受控变量处理，而不是一个背景设定。

EN — The reference determines which questions deconvolution can answer at all. The consistent lesson from 2024–2026 benchmarks is that **reference construction often contributes more than algorithm choice**. The reference is therefore treated as a controlled variable, not as background configuration.

### §3.1 细胞类型粒度阶梯 / Granularity Ladder

中文 — 预先定义三级粒度，**所有方法在所有三级上都要跑**，因为溢出行为随粒度剧烈变化：

EN — Define three levels a priori and **run every method at all three**, because spillover behavior changes sharply with granularity:

| Level | Cell types |
|---|---|
| **L1 (coarse)** | malignant/epithelial, T/NK, B/plasma, myeloid, fibroblast, endothelial |
| **L2 (intermediate)** | malignant, CD8 T, CD4 T (conv), Treg, NK, B, plasma, monocyte/macrophage, cDC, pDC, mast, fibroblast, endothelial |
| **L3 (fine)** | L2 with CD8 split (naive / effector-memory / exhausted-terminal / tissue-resident), CD4 split (naive / memory / Tfh / Treg), macrophage polarization axes, fibroblast subsets (incl. matrix CAF) |

中文 — 规则：**主分析终点定义在 L2 的 `CD8 T`**（§7.9）。L1 用作稳定性下界（如果 L1 都不稳定，问题不在粒度），L3 用作解释性上界，但 L3 的 CD8 亚群受溢出影响最大，不得作为主终点。`Tfh` 在 L3 出现，是与 TLS 生物学最相关的 T 细胞状态（§6），但仍属探索性。

EN — Rule: **the primary endpoint is defined on L2 `CD8 T`** (§7.9). L1 serves as a stability floor (if L1 is unstable, granularity is not the problem); L3 serves as an interpretive ceiling, but L3 CD8 subsets suffer the most spillover and may not carry the primary endpoint. `Tfh` appears at L3 and is the T-cell state most relevant to TLS biology (§6), yet remains exploratory.

### §3.2 参考 scRNA-seq 的最低要求 / Minimum Reference Requirements

| Requirement | Specification |
|---|---|
| Donors | ≥5 independent donors (hard requirement for MuSiC's cross-subject weighting) |
| Tissue match | Same tumor type and, where possible, same site (primary vs metastasis) |
| Treatment context | Record whether the reference is treatment-naive; if it is and the bulk is on-treatment, MuSiC2 becomes mandatory as a comparator arm (§1.3) |
| Chemistry | Record 10x version / Smart-seq2 / snRNA-seq; snRNA-seq triggers the §2.3 point 3 caveat |
| Malignant annotation | Malignant cells identified by CNV inference (e.g. inferCNV/CopyKAT) **plus** marker evidence, not markers alone |
| Cell counts | Pre-specify a minimum cell count per cell type per level; types below the floor are merged upward, and the merge is recorded |
| Doublets/ambient | Doublet removal and ambient-RNA correction applied and documented; ambient contamination directly inflates cross-type signal and mimics spillover |

中文 — **恶性细胞的处理是最关键的一步**。BayesPrism 期望恶性细胞作为一个带 `key` 的类型，其表达谱按样本更新；这正是 §7.4 能工作的机制。如果把恶性细胞折叠进"上皮"并固定其谱，就丧失了区分"每细胞表达"与"细胞数量"的能力。同时，**参考中恶性细胞的患者来源应尽量多样**，因为恶性谱的患者间异质性远大于免疫细胞。

EN — **Malignant-cell handling is the pivotal step.** BayesPrism expects malignant cells as a `key`-flagged type whose profile is updated per sample; that mechanism is exactly what makes §7.4 work. Folding malignant cells into a fixed "epithelial" profile destroys the ability to separate per-cell expression from cell abundance. Also, **draw malignant cells from as many patients as possible**, since between-patient heterogeneity of the malignant profile far exceeds that of immune cells.

### §3.3 参考 QC：伪 bulk 自检 / Reference QC via Pseudobulk Self-Test

中文 — 在把参考用于真实 bulk 之前，必须通过自检：用参考自身构建**异质性伪 bulk**（见 §8.1 关于"随机取细胞"为何不足），跑一遍全部方法，确认每个 L2 细胞类型的估计与已知真值达到预先设定的相关下限。未通过的参考不得进入主流程。

EN — Before a reference touches real bulk, it must pass a self-test: build **heterogeneous pseudobulk** from the reference itself (see §8.1 on why random cell sampling is insufficient), run all methods, and confirm each L2 cell type meets a pre-specified correlation floor against known truth. A reference failing the self-test does not enter the main pipeline.

---

## §4 各方法运行规程 / Per-Method Run Procedures

中文 — 以下为运行骨架。**所有代码块均为模板**：在执行前对照已安装的软件包版本核对 API，并把确认过的版本号写入 `env/versions.lock`。参数不得依赖默认值，全部在 `config/` 中显式声明。

EN — Run skeletons follow. **All code blocks are templates**: verify the API against your installed package versions before executing, and record confirmed versions in `env/versions.lock`. No parameter may rely on a default; declare all of them explicitly in `config/`.

### §4.1 InstaPrism / BayesPrism

```r
# methods/deconv_advanced/run/run_instaprism.R
suppressPackageStartupMessages({
  library(InstaPrism); library(Matrix); library(yaml)
})
cfg <- yaml::read_yaml("config/instaprism.yaml")
set.seed(cfg$seed)

# bulk: genes x samples, RAW COUNTS
# sc  : genes x cells, RAW COUNTS ; labels aligned to columns of sc
bulk <- readRDS(cfg$bulk_counts_rds)
sc   <- readRDS(cfg$sc_counts_rds)
ct   <- readRDS(cfg$cell_type_labels_rds)    # L1 / L2 / L3 depending on cfg$level
cs   <- readRDS(cfg$cell_state_labels_rds)   # finer labels nested within ct

genes <- intersect(rownames(bulk), rownames(sc))
stopifnot(length(genes) >= cfg$min_shared_genes)
bulk <- bulk[genes, , drop = FALSE]; sc <- sc[genes, , drop = FALSE]

res <- InstaPrism(
  input_type       = "raw",
  sc_Expr          = as.matrix(sc),
  bulk_Expr        = as.matrix(bulk),
  cell.type.labels = ct,
  cell.state.labels= cs,
  n.iter           = cfg$n_iter,
  update           = cfg$update_reference   # sensitivity axis, see 8.3
)

theta <- t(res@Post.ini.ct@theta)            # samples x cell types
saveRDS(theta, file.path(cfg$outdir, "theta_celltype.rds"))

# Compartment-resolved expression for 7.4.
# Full BayesPrism exposes get.exp(bp, state.or.type = "type", cell.name = <label>).
# InstaPrism stores Z compressed as a 2D scaling matrix and reconstructs it on demand;
# the accessor name differs between package versions, so bind it in config rather than
# hard-coding it, and assert the returned object is genes x samples before use.
get_Z <- match.fun(cfg$instaprism_Z_accessor)          # e.g. "get_Z_array" / "reconstruct_Z"
mal_expr <- get_Z(res, cell.type = cfg$malignant_label)
stopifnot(nrow(mal_expr) == length(genes), ncol(mal_expr) == ncol(bulk))
saveRDS(mal_expr, file.path(cfg$outdir, "Z_malignant.rds"))
```

```yaml
# config/instaprism.yaml
seed: 20260101
level: "L2"
n_iter: 100
update_reference: true
malignant_label: "Malignant"
instaprism_Z_accessor: "reconstruct_Z"   # pin to the accessor your installed version exposes
min_shared_genes: 5000
outlier_cut: 0.01          # applied upstream when filtering the bulk matrix
outlier_fraction: 0.10
exclude_gene_classes: ["MT-", "RPL", "RPS", "HB[ABDGEQZ]"]
outdir: "results/deconv/{cohort_id}/instaprism/L2"
```

中文 — 若同时运行完整 BayesPrism 作为等价性检查（推荐至少在一个队列上做），使用 `new.prism(...)` + `run.prism(...)`，`key = cfg$malignant_label`，并把两者的 $\theta$ 做逐细胞类型相关比较，结果记入 `results/qc/instaprism_vs_bayesprism.tsv`。

EN — If you also run full BayesPrism as an equivalence check (recommended on at least one cohort), use `new.prism(...)` + `run.prism(...)` with `key = cfg$malignant_label`, compare the two $\theta$ matrices per cell type, and record the comparison in `results/qc/instaprism_vs_bayesprism.tsv`.

### §4.2 CIBERSORTx

```bash
# methods/deconv_advanced/run/run_cibersortx.sh
set -euo pipefail
: "${CSX_USER:?set CSX_USER}" ; : "${CSX_TOKEN:?set CSX_TOKEN}"
COHORT="$1" ; LEVEL="$2"
IN="$PWD/work/csx/${COHORT}/${LEVEL}/in" ; OUT="$PWD/results/deconv/${COHORT}/cibersortx/${LEVEL}"
mkdir -p "$IN" "$OUT"

# Step 1 - fractions + signature matrix (S-mode: scRNA-seq-derived signature)
docker run --rm \
  -v "$IN":/src/data -v "$OUT":/src/outdir \
  cibersortx/fractions \
  --username "$CSX_USER" --token "$CSX_TOKEN" \
  --single_cell TRUE \
  --refsample refsample.txt \
  --mixture mixture_tpm.txt \
  --rmbatchSmode TRUE \
  --fraction 0 \
  --perm 100 \
  --QN FALSE \
  --verbose TRUE

# Step 2 - high-resolution cell-type-specific expression (input to 7.4 and EcoTyper)
docker run --rm \
  -v "$IN":/src/data -v "$OUT":/src/outdir \
  cibersortx/hires \
  --username "$CSX_USER" --token "$CSX_TOKEN" \
  --mixture mixture_tpm.txt \
  --sigmatrix CIBERSORTx_sigmatrix.txt \
  --cibresults CIBERSORTx_Results.txt \
  --threads 8 \
  --verbose TRUE
```

中文 — 配置决策必须显式记录：(a) `--rmbatchSmode` 对 scRNA-seq 来源签名，`--rmbatchBmode` 对芯片/外部签名，二选一；(b) `--QN FALSE` 用于 RNA-seq 混合物，`TRUE` 通常只在芯片场景考虑；(c) `--fraction 0` 表示 refsample 非流式纯化样本；(d) `--perm` 决定每样本 p 值，本项目中该 p 值**不用于任何筛选**，仅作 QC 记录。

EN — Configuration decisions must be recorded explicitly: (a) `--rmbatchSmode` for scRNA-seq-derived signatures versus `--rmbatchBmode` for array/external signatures — pick exactly one; (b) `--QN FALSE` for RNA-seq mixtures, with `TRUE` considered only in array settings; (c) `--fraction 0` declares that the refsample is not flow-purified; (d) `--perm` drives the per-sample p-value, which in this project is **never used for filtering** and is retained only as a QC record.

### §4.3 MuSiC / MuSiC2

```r
# methods/deconv_advanced/run/run_music.R
suppressPackageStartupMessages({
  library(MuSiC); library(MuSiC2); library(SingleCellExperiment); library(yaml)
})
cfg <- yaml::read_yaml("config/music.yaml"); set.seed(cfg$seed)

bulk <- as.matrix(readRDS(cfg$bulk_counts_rds))     # genes x samples, counts
sce  <- readRDS(cfg$sce_rds)                        # colData has $cellType, $sampleID
stopifnot(length(unique(sce$sampleID)) >= 5)        # 3.2 hard requirement

# --- MuSiC, tree-guided with an a-priori immunological grouping (1.3 point b)
prop <- music_prop(
  bulk.mtx  = bulk,
  sc.sce    = sce,
  clusters  = "cellType",
  samples   = "sampleID",
  select.ct = cfg$cell_types,
  verbose   = TRUE
)
saveRDS(prop$Est.prop.weighted, file.path(cfg$outdir, "theta_music.rds"))

# --- MuSiC2, when reference condition != bulk condition
# NOTE: the split encodes REFERENCE-LIKE vs QUERY condition, never responder status.
idx_ref_like <- readRDS(cfg$reference_like_sample_ids_rds)
prop2 <- music2_prop(
  bulk.control.mtx = bulk[, idx_ref_like,  drop = FALSE],
  bulk.case.mtx    = bulk[, setdiff(colnames(bulk), idx_ref_like), drop = FALSE],
  sc.sce           = sce,
  clusters         = "cellType",
  samples          = "sampleID",
  select.ct        = cfg$cell_types,
  n_resample       = cfg$n_resample,
  sample_prop      = cfg$sample_prop,
  cutoff_c         = cfg$cutoff_common,
  cutoff_D         = cfg$cutoff_rare
)
saveRDS(prop2$Est.prop, file.path(cfg$outdir, "theta_music2.rds"))
saveRDS(prop2$DE.genes, file.path(cfg$outdir, "music2_removed_DE_genes.rds"))
```

中文 — `music2_removed_DE_genes.rds` 必须保存并检查：如果 `TACSTD2` 出现在被剔除的细胞类型特异 DE 基因列表中，这是一个**重要的方法学事实**，说明参考与 bulk 在该基因上存在条件性差异，必须在 §7.4 的解释中显式声明（但不构成结果）。

EN — Always save and inspect `music2_removed_DE_genes.rds`: if `TACSTD2` appears among the removed cell-type-specific DE genes, that is a **material methodological fact** indicating condition-dependent divergence between reference and bulk for this gene, and it must be declared explicitly when interpreting §7.4 (while not itself constituting a result).

### §4.4 EcoTyper (Recovery)

```yaml
# config/ecotyper_recovery_bulk.yml
default:
  Input:
    Discovery dataset name: "Carcinoma"          # pre-built pan-carcinoma CS/CE model
    Expression matrix: "work/ecotyper/{cohort_id}/expr_log2_scaled.tsv"
    Annotation file: "work/ecotyper/{cohort_id}/annotation.tsv"
    Annotation file column to scale by: "TumorType"   # per-tumor-type unit-variance scaling
  Output:
    Output folder: "results/ecotyper/{cohort_id}"
  Pipeline settings:
    Number of threads: 8
    z-score cutoff: 1.65
```

```bash
# methods/deconv_advanced/run/run_ecotyper.sh
set -euo pipefail
Rscript EcoTyper_recovery_bulk.R -c config/ecotyper_recovery_bulk.yml
```

中文 — 前置步骤（必须在 EcoTyper 之前完成，且写成独立脚本以便审计）：把表达矩阵取 log2，然后**在每个癌种/数据集内**对每个基因做单位方差缩放。这一步若遗漏，recovery 仍会"成功"运行并输出数值，但数值无意义——这是该工具最危险的静默失败。恢复后，未通过 z-score 阈值的细胞状态整体剔除（§1.4 点 c），剔除清单写入 `results/ecotyper/{cohort_id}/dropped_states.tsv`。

EN — Prerequisite step (must run before EcoTyper, as a separate auditable script): log2-transform the expression matrix, then scale each gene to unit variance **within each tumor type / dataset**. Omitting this still lets recovery "succeed" and emit numbers, but the numbers are meaningless — this is the tool's most dangerous silent failure. After recovery, drop cell states failing the z-score threshold wholesale (§1.4 point c) and write the dropped list to `results/ecotyper/{cohort_id}/dropped_states.tsv`.

### §4.5 InstantDL (Image-Side, Role A)

中文 — 仅在 `has_matched_imaging == TRUE` 的队列上运行。产出两个 slide 级读数，作为 §8.1 的正交金标准。

EN — Run only on cohorts with `has_matched_imaging == TRUE`. Produces two slide-level readouts used as orthogonal ground truth in §8.1.

```yaml
# config/instantdl_cd8.yaml  (task 1: CD8+ cell density from mIF)
use_algorithm: "InstanceSegmentation"
path: "data/imaging/{cohort_id}/cd8"
pretrained_weights: "models/instantdl/mrcnn_nuclei.h5"
batchsize: 2
iterations_over_dataset: 100
loss_function: "binary_crossentropy"
Early_stopping_patience: 10
calculate_uncertainty: true      # per-slide confidence weight for 8.2
```

```yaml
# config/instantdl_tls.yaml  (task 2: TLS region detection from H&E)
use_algorithm: "SemanticSegmentation"
path: "data/imaging/{cohort_id}/tls"
pretrained_weights: "models/instantdl/unet_tissue.h5"
batchsize: 4
iterations_over_dataset: 150
loss_function: "dice_loss"
calculate_uncertainty: true
```

中文 — 下游聚合规程：(a) CD8 密度 = 分割得到的 CD8⁺ 细胞数 / 组织面积（mm²），并另算"CD8⁺ 占有核细胞比例"以便与**比例**型去卷积输出同尺度比较——这一点很重要，因为去卷积给的是比例而非密度，直接把密度与比例相关会引入组织细胞密度这个混杂；(b) TLS 读数 = 每 mm² 的 TLS 数量与 TLS 面积占比，并按成熟度（有无生发中心）分层，需要病理学家对模型输出做抽样复核；(c) 抽样复核比例、复核者、以及 Cohen's κ 的计算方式在 `methods/deconv_advanced/imaging_review_sop.md` 中另行规定。

EN — Downstream aggregation procedure: (a) CD8 density = segmented CD8⁺ cell count / tissue area (mm²), **plus** a "CD8⁺ share of nucleated cells" so that it is on the same scale as the **fraction**-valued deconvolution output — this matters because correlating a density against a fraction imports tissue cellularity as a confounder; (b) TLS readouts = TLS count per mm² and TLS area fraction, stratified by maturity (germinal center present or absent), with pathologist spot-review of model output; (c) the review sampling rate, reviewer identity, and Cohen's κ computation are specified separately in `methods/deconv_advanced/imaging_review_sop.md`.

---

## §5 输出统一化 / Output Harmonization

### §5.1 公共细胞类型本体 / Common Cell-Type Ontology

中文 — 各方法输出的标签命名不同（如 `T cells CD8`、`CD8_Tcell`、`CD8T`）。建立 `config/ontology_map.tsv`，把每个方法在每个粒度级别的原始标签映射到统一术语，**映射为多对一时求和**，且求和只允许在同一谱系内进行。映射表是版本化产物，任何修改必须走 git 提交并说明理由。

EN — Methods emit different labels (`T cells CD8`, `CD8_Tcell`, `CD8T`). Maintain `config/ontology_map.tsv` mapping each method's raw label at each granularity level onto canonical terms; **sum when the mapping is many-to-one**, and only within a lineage. The map is a versioned artifact: any edit goes through a git commit with a stated rationale.

```tsv
method       level  raw_label            canonical      lineage
instaprism   L2     CD8_T                CD8_T          T_NK
cibersortx   L2     T cells CD8          CD8_T          T_NK
music        L2     CD8Tcell             CD8_T          T_NK
cibersortx   L2     T cells follicular helper  CD4_T     T_NK
ecotyper     L2     CD8.T.cells          CD8_T          T_NK
```

### §5.2 相对与绝对 / Relative vs Absolute

中文 — 三类输出不可混用：

1. **相对比例**（和为 1，含或不含恶性细胞）— 大多数方法的默认。
2. **绝对/半绝对得分** — CIBERSORTx absolute mode、xCell 类富集得分。
3. **密度**（细胞数/面积）— 仅影像侧（§4.5）。

规则：主分析统一使用**含恶性细胞在内、和为 1 的相对比例**，因为这是 BayesPrism 与 MuSiC 的原生输出尺度，且是 §7.5 成分处理的前提。若某方法只给不含恶性细胞的免疫内比例，必须记录并在 §8 中单列，不得强行重标定。

EN — Three output kinds must not be mixed:

1. **Relative fractions** (sum to one, with or without the malignant compartment) — most methods' default.
2. **Absolute / semi-absolute scores** — CIBERSORTx absolute mode, xCell-style enrichment scores.
3. **Densities** (cells per area) — image side only (§4.5).

Rule: the primary analysis uses **malignant-inclusive relative fractions summing to one**, since that is the native scale for BayesPrism and MuSiC and the precondition for the compositional handling in §7.5. If a method only yields immune-internal fractions excluding malignant cells, record that and report it separately in §8; do not force-rescale it.

### §5.3 成分变换 / Compositional Transform

中文 — 这是 R0-1 的实现。比例向量 $\theta_n \in \Delta^{k-1}$ 位于单纯形上，欧氏空间的方法（Pearson、OLS、PCA）在其上不成立。规程：

1. **零值替换**：结构性零与取样零都会阻断对数。使用贝叶斯乘性替换（`zCompositions::cmultRepl`）而非加常数，因为加常数会依赖于常数大小改变所有比值。
2. **CLR 变换**：$\mathrm{clr}(\theta)_j = \log \theta_j - \frac{1}{k}\sum_{i=1}^{k}\log \theta_i$。CLR 后的坐标可进入线性模型，但注意 CLR 坐标之和恒为 0（奇异协方差），不可把全部 $k$ 个坐标同时放进一个回归。
3. **替代方案 — 参照比值（ALR / 比值终点）**：当有明确的生物学参照时，直接用 $\log(\theta_{\text{CD8}} / \theta_{\text{ref}})$ 更易解释。本项目预设两个参照：`log(CD8_T / CD4_T_conv)`（T 细胞内部组成，免疫浸润总量被约掉）与 `log(CD8_T / Malignant)`（每肿瘤细胞的 CD8 负荷）。后者与 §7.4 的区室视角天然对齐。
4. **替代方案 — beta 回归**：把单个比例作为 (0,1) 上的结局，用 `betareg`，配合 logit 链接。适合只关心一个细胞类型且希望保持原尺度可解释性的场景。

EN — This implements R0-1. A fraction vector $\theta_n \in \Delta^{k-1}$ lives on the simplex, where Euclidean methods (Pearson, OLS, PCA) are not valid. Procedure:

1. **Zero replacement**: both structural and sampling zeros block the logarithm. Use Bayesian multiplicative replacement (`zCompositions::cmultRepl`) rather than adding a constant, since an additive pseudo-count changes every ratio as a function of the constant chosen.
2. **CLR transform**: $\mathrm{clr}(\theta)_j = \log \theta_j - \frac{1}{k}\sum_{i=1}^{k}\log \theta_i$. CLR coordinates may enter linear models, but note they sum to zero (singular covariance), so never place all $k$ coordinates in one regression.
3. **Alternative — reference ratios (ALR / ratio endpoints)**: with a defensible biological denominator, $\log(\theta_{\text{CD8}} / \theta_{\text{ref}})$ is more interpretable. Two denominators are pre-specified here: `log(CD8_T / CD4_T_conv)` (composition within T cells, cancelling total immune infiltration) and `log(CD8_T / Malignant)` (CD8 burden per tumor cell). The latter aligns naturally with the compartment view of §7.4.
4. **Alternative — beta regression**: treat a single fraction as a (0,1)-valued outcome with `betareg` and a logit link. Suited to a single cell type of interest where the original scale should stay interpretable.

```r
# methods/deconv_advanced/lib/compositional.R
library(zCompositions); library(compositions)

clr_transform <- function(theta, label = "theta") {
  stopifnot(all(theta >= 0), all(abs(rowSums(theta) - 1) < 1e-6))
  if (any(theta == 0)) {
    theta <- zCompositions::cmultRepl(theta, label = 0, method = "CZM",
                                      output = "prop", suppress.print = TRUE)
  }
  as.matrix(compositions::clr(compositions::acomp(theta)))
}

ratio_endpoint <- function(theta, num, den) {
  eps <- min(theta[theta > 0]) / 2          # documented, reported, varied in 8.3
  log((theta[, num] + eps) / (theta[, den] + eps))
}
```

### §5.4 不确定性传播 / Uncertainty Propagation

中文 — 点估计的比例被下游当作已知常数使用，会低估标准误。规程：对每个队列 × 方法生成 $B$ 次 bootstrap 的 $\theta$（重抽样参考细胞，而非重抽样 bulk 样本 —— 前者反映参考不确定性，后者已由回归本身处理），把 §7 的关联估计在每个 bootstrap 副本上重跑一遍，最终区间取 bootstrap 分布的分位数。BayesPrism 另可直接使用后验抽样。InstaPrism 的速度使 $B \geq 200$ 可行；若使用完整 BayesPrism，$B$ 可降低并在报告中声明。

EN — Treating point-estimate fractions as known constants downstream understates standard errors. Procedure: for each cohort × method, generate $B$ bootstrap replicates of $\theta$ by resampling **reference cells** (not bulk samples — the former captures reference uncertainty, while the latter is already handled by the regression), rerun the §7 association on every replicate, and take interval endpoints from the bootstrap distribution. BayesPrism can alternatively use posterior draws directly. InstaPrism's speed makes $B \geq 200$ feasible; with full BayesPrism, reduce $B$ and declare it in the report.

---

## §6 TLS 量化模块 / TLS Quantification Module

中文 — bulk 数据无法直接观察三级淋巴结构（TLS）——TLS 是空间组织现象，而 bulk 是空间平均。所有 bulk TLS 读数都是**代理**。因此本模块的规程是：并行计算多个已发表签名，把它们当作同一潜变量的多个有噪测量，并在主分析中使用**预先指定的单一主签名 + 其余作为敏感性**。

EN — Bulk data cannot observe tertiary lymphoid structures directly — TLS is a spatial organizational phenomenon and bulk is a spatial average. Every bulk TLS readout is a **proxy**. The module therefore computes several published signatures in parallel, treats them as multiple noisy measurements of one latent variable, and uses **one pre-specified primary signature with the rest as sensitivity analyses**.

### §6.1 签名集合 / Signature Set

| ID | Composition | Notes |
|---|---|---|
| `TLS_12CK` | CCL2, CCL3, CCL4, CCL5, CCL8, CCL18, CCL19, CCL21, CXCL9, CXCL10, CXCL11, CXCL13 | The 12-chemokine score; the most widely applied across solid tumors |
| `TLS_CABRITA` | CD79B, CD1D, CCR6, LAT, SKAP1, CETP, EIF1AY, RBP5, PTGDS | Melanoma-derived; **contains the Y-linked `EIF1AY`** — see §6.2 |
| `TLS_HALLMARK` | CCL19, CCL21, CXCL13, CCR7, CXCR5, SELL, LAMP3 | Compact hallmark set |
| `TLS_TFH` | CXCL13, CD200, FBLN7, ICOS, SGPP2, SH2D1A, TIGIT, PDCD1 | Tfh-centric; overlaps T-cell activation, so it is not independent of CD8 readouts |
| `TLS_MEYLAN` | as published (RCC-derived TLS imprint) | Pin the exact gene list and its source table in `config/tls_signatures.tsv` |
| `TLS_BCELL` | BANK1, CD19, CD22, CD79A, CR1, CR2, FCRL2, MS4A1, PAX5, FCER2, MZB1 | B/plasma-lineage composition; strongly collinear with the B-cell fraction from §4 |

中文 — 主签名的预先指定规则：若队列有配套影像（§4.5），主签名选为**在该癌种中与影像 TLS 计数相关性最高者，且该选择必须在锁定关联分析之前、仅使用影像与签名（不涉及 `TACSTD2`）完成**。若无配套影像，默认主签名为 `TLS_12CK`（跨癌种应用最广、最易与文献对齐）。选择过程写入 `results/tls/primary_signature_selection.md`。

EN — Rule for pre-specifying the primary signature: if a cohort has matched imaging (§4.5), select as primary the signature with the highest correlation to image-based TLS counts in that tumor type — and make that selection **before locking the association analysis, using only imaging and signatures, never touching `TACSTD2`**. Without matched imaging, default to `TLS_12CK` (broadest cross-tumor use, easiest alignment with the literature). Record the selection in `results/tls/primary_signature_selection.md`.

### §6.2 强制性混杂检查 / Mandatory Confounder Checks

中文 — 两项检查是强制的：

1. **性染色体污染**：`TLS_CABRITA` 含 `EIF1AY`，位于 Y 染色体。在性别混合队列中，该签名分数会部分反映患者性别而非 TLS。**必须**在计算前对该签名做以下之一：剔除 `EIF1AY` 并记为 `TLS_CABRITA_noY`（推荐）；或在所有使用该签名的模型中强制包含性别协变量。同时保留原版用于与文献对齐，但原版不得作为主签名。
2. **与 CD8 读数的循环性**：`TLS_12CK` 含 `CXCL9/10/11`（IFN-γ 诱导，直接来自 T 细胞/髓系激活），`TLS_TFH` 含 `PDCD1`、`TIGIT`、`ICOS`。把这些签名与 CD8 比例相关，部分是在测量同一组基因。因此 §7 中 CD8 与 TLS 是**两个独立的主终点**，绝不构建 "CD8+TLS 复合评分"，且必须报告 TLS 签名与 CD8 比例之间的相关作为解释背景（这属于方法学诊断，不是结果）。

EN — Two checks are mandatory:

1. **Sex-chromosome contamination.** `TLS_CABRITA` includes `EIF1AY`, a Y-linked gene. In sex-mixed cohorts the score partly reflects patient sex rather than TLS. **Required** before scoring: either drop `EIF1AY` and record the variant as `TLS_CABRITA_noY` (preferred), or force a sex covariate into every model using the signature. Keep the original for literature alignment, but the original may not serve as the primary signature.
2. **Circularity with the CD8 readout.** `TLS_12CK` contains `CXCL9/10/11` (IFN-γ-induced, arising directly from T-cell/myeloid activation) and `TLS_TFH` contains `PDCD1`, `TIGIT`, `ICOS`. Correlating these signatures with a CD8 fraction partly re-measures the same genes. Consequently CD8 and TLS are **two separate primary endpoints** in §7, no "CD8+TLS composite score" is ever constructed, and the correlation between the TLS signature and the CD8 fraction must be reported as interpretive background (a methodological diagnostic, not a result).

### §6.3 打分方法 / Scoring Methods

```r
# methods/deconv_advanced/lib/tls_score.R
suppressPackageStartupMessages({ library(GSVA); library(singscore); library(matrixStats) })

# expr_log  : genes (HGNC) x samples, log2(TPM + 1)
# sig_full  : named list of complete signature gene sets, loaded from config/tls_signatures.tsv
# cov_floor : minimum retained fraction of signature genes; below it the signature is dropped
score_tls <- function(expr_log, sig_full, cov_floor, method = c("ssgsea", "singscore", "zmean")) {
  method <- match.arg(method)
  sig_list <- lapply(sig_full, function(g) intersect(g, rownames(expr_log)))
  coverage <- vapply(names(sig_full),
                     function(n) length(sig_list[[n]]) / length(sig_full[[n]]), numeric(1))
  dropped <- names(coverage)[coverage < cov_floor]
  if (length(dropped)) message("dropped for low coverage: ", paste(dropped, collapse = ", "))
  sig_list <- sig_list[setdiff(names(sig_list), dropped)]   # drop whole, never score partially
  attr(sig_list, "coverage") <- coverage                    # persisted to the run manifest
  switch(method,
    ssgsea   = GSVA::gsva(GSVA::ssgseaParam(expr_log, sig_list, normalize = TRUE)),
    singscore = {
      rk <- singscore::rankGenes(expr_log)
      t(vapply(sig_list, function(g)
        singscore::simpleScore(rk, upSet = g)$TotalScore, numeric(ncol(expr_log))))
    },
    zmean = {
      z <- (expr_log - rowMeans(expr_log)) / matrixStats::rowSds(expr_log)
      t(vapply(sig_list, function(g) colMeans(z[g, , drop = FALSE]), numeric(ncol(expr_log))))
    })
}
```

中文 — 打分方法的选择本身是一个敏感性轴（§8.3）。关键差异：`zmean` 与 `ssgsea` 的分数依赖于队列内其他样本（队列相对），因此**不可跨队列直接比较绝对值**，只能在队列内做关联然后元分析（这正是 R0-5 的另一个理由）；`singscore` 基于样本内基因秩，是样本自足的，可跨队列比较。若目标包含跨队列绝对可比性，主打分方法应为 `singscore`。覆盖度下限（例如签名基因在表达矩阵中的留存比例）必须预先设定；低于下限的签名整体剔除，不做部分打分。

EN — The scoring method is itself a sensitivity axis (§8.3). The key difference: `zmean` and `ssgsea` scores depend on the other samples in the cohort (cohort-relative), so **absolute values are not comparable across cohorts** — associate within cohort and then meta-analyze, which is another reason for R0-5. `singscore` uses within-sample gene ranks and is self-contained, hence cross-cohort comparable. If cross-cohort absolute comparability is required, make `singscore` the primary scorer. Set a coverage floor (fraction of signature genes retained in the expression matrix) a priori; signatures below the floor are dropped entirely rather than partially scored.

### §6.4 细胞组成侧的 TLS 代理 / Composition-Side TLS Proxies

中文 — 除签名外，从 §4 的输出构造两个组成侧代理，作为签名的独立佐证：(a) B/浆细胞比例（L2）；(b) 若 L3 可用，`Tfh` 比例与 `log(B / Malignant)`。EcoTyper 侧则记录 B 细胞与 CD4 的相关状态及其所属 ecotype 的丰度。这三条证据线（签名、比例、细胞状态）在 §8.4 中共同参与结论稳定性评估。

EN — Beyond signatures, build two composition-side proxies from §4 output as independent corroboration: (a) the B/plasma fraction at L2; (b) where L3 is available, the `Tfh` fraction and `log(B / Malignant)`. On the EcoTyper side, record the relevant B-cell and CD4 states and the abundance of their parent ecotype. These three evidence lines (signature, fraction, cell state) jointly enter the conclusion-stability assessment in §8.4.

---

## §7 TACSTD2 × CD8/TLS 关联协议 / The TACSTD2 Correlation Protocol

### §7.1 为什么朴素相关是错的 / Why the Naive Correlation Fails

中文 — 最直接的做法是：取 bulk 中 `TACSTD2` 的 log2(TPM+1)，取去卷积得到的 CD8 比例，算 Spearman。这个数几乎一定是"某个东西"，但它至少混合了四个来源：

1. **纯度混杂**：`TACSTD2` 主要由上皮/恶性细胞表达。样本肿瘤细胞越多，bulk `TACSTD2` 越高，同时**所有**免疫比例机械性越低。这一条足以单独产生一个可观的负相关，且与 TROP2 的任何免疫生物学无关。
2. **成分闭合**：即使没有纯度混杂，比例之和为 1 也使各细胞类型的估计相互约束。
3. **区室歧义**：bulk `TACSTD2` = (每恶性细胞表达) × (恶性细胞比例) + (其他区室贡献)。前两项无法从 bulk 分离。
4. **技术共变**：FFPE 降解、文库深度、3′ 偏倚会同时影响单基因表达与去卷积输出。

因此本协议的核心不是"算相关"，而是**先把估计量定义清楚**。

EN — The obvious move is: take bulk log2(TPM+1) of `TACSTD2`, take the deconvolved CD8 fraction, compute Spearman. That number is certainly *something*, but it mixes at least four sources:

1. **Purity confounding.** `TACSTD2` is expressed mainly by epithelial/malignant cells. More tumor cells means higher bulk `TACSTD2` and mechanically lower **every** immune fraction. This alone can manufacture a sizable negative correlation with no TROP2 immunobiology involved.
2. **Compositional closure.** Even absent purity confounding, the sum-to-one constraint couples the cell-type estimates.
3. **Compartment ambiguity.** Bulk `TACSTD2` = (expression per malignant cell) × (malignant fraction) + (contributions from other compartments). Bulk alone cannot separate the first two terms.
4. **Technical covariation.** FFPE degradation, library depth, and 3′ bias affect single-gene expression and deconvolution output simultaneously.

The core of this protocol is therefore not "computing a correlation" but **defining the estimand first**.

### §7.2 三个不同的估计量 / Three Distinct Estimands

中文 — 下面三个问题在文献中经常被混为一谈，但需要不同的分析：

EN — These three questions are routinely conflated in the literature and require different analyses:

| ID | Question | Exposure variable | Interpretation |
|---|---|---|---|
| **E1** | Do tumors with more bulk `TACSTD2` signal have less CD8 / less TLS? | bulk log2(TPM+1) of `TACSTD2` | Descriptive, **purity-confounded by construction**. Useful only as a biomarker-as-measured statement (e.g. "what an assay reading a bulk sample would see"), never as a mechanistic claim |
| **E2** | Holding tumor content constant, does `TACSTD2` track immune composition? | bulk `TACSTD2`, with purity as a covariate | Partial association; removes source 1 of §7.1 but not compartment ambiguity |
| **E3** | Do malignant cells expressing more `TACSTD2` **per cell** sit in tumors with less CD8 / less TLS? | $E^{\text{mal}}_{\texttt{TACSTD2}}$ from BayesPrism $Z$ or CIBERSORTx HiRes | The mechanistically interpretable estimand; the one that maps onto TROP2 cell-biology hypotheses |

中文 — **主估计量为 E3**，E2 为共同主要，E1 为强制报告的描述性对照（R0-2）。理由：现有文献在不同癌种中报告了方向相反的关联，而方向差异的一个直接的方法学候选解释就是 E1 与 E3 被当成了同一个量。把三者分开报告，使得"方向相反"这件事本身可被诊断，而不是被当作生物学结论。本手册不预设任何方向。

EN — **E3 is the primary estimand**, E2 is co-primary, and E1 is a mandatory descriptive control (R0-2). Rationale: the existing literature reports associations in opposite directions across tumor types, and one immediate methodological candidate explanation is that E1 and E3 have been treated as the same quantity. Reporting all three separately makes directional disagreement itself diagnosable rather than a biological conclusion. This playbook pre-specifies no direction.

### §7.3 纯度的获取与校正 / Obtaining and Adjusting for Purity

中文 — 纯度估计的来源按优先级排序，**优先级即循环性风险的倒序**：

1. **DNA 来源（首选）**：ABSOLUTE / FACETS / Sequenza / PureCN，基于 WES/WGS/SNP array。与表达数据正交，因此用它做协变量不会引入循环性。
2. **甲基化来源**：若有 450K/EPIC 数据，基于甲基化的纯度同样与表达正交。
3. **表达来源（次选）**：ESTIMATE 的 tumor purity。与暴露变量共享表达平台，存在部分循环性，但仍显著优于不校正。
4. **去卷积来源（谨慎）**：直接用同一次去卷积得到的恶性细胞比例。**这是循环的**：$\theta_{\text{mal}}$ 与 $\theta_{\text{CD8}}$ 来自同一次求解、共享闭合约束，把前者作为协变量去解释后者会诱导偏倚。仅当 1–3 全部不可得时使用，且必须在结果中显著标注，并配合 §7.5 的成分方法（成分方法本身已部分吸收闭合效应）。

规程：在 `data/purity/{cohort_id}.tsv` 中记录所有可得的纯度来源；主分析使用可得的最高优先级来源；把"纯度来源"作为 §8.3 的一条敏感性轴。

EN — Purity sources in priority order, where **priority is the inverse of circularity risk**:

1. **DNA-based (preferred)**: ABSOLUTE / FACETS / Sequenza / PureCN from WES/WGS/SNP array. Orthogonal to expression, so using it as a covariate introduces no circularity.
2. **Methylation-based**: with 450K/EPIC data available, methylation purity is likewise orthogonal to expression.
3. **Expression-based (second choice)**: ESTIMATE tumor purity. Shares the expression platform with the exposure, so partly circular, but still far better than no adjustment.
4. **Deconvolution-derived (use with care)**: the malignant fraction from the same deconvolution run. **This is circular**: $\theta_{\text{mal}}$ and $\theta_{\text{CD8}}$ come from one solve and share the closure constraint, so conditioning the latter on the former induces bias. Use only when 1–3 are all unavailable, flag it prominently in results, and pair it with the compositional approach in §7.5 (which already absorbs part of the closure effect).

Procedure: record every available purity source in `data/purity/{cohort_id}.tsv`; the primary analysis uses the highest-priority source available; carry "purity source" as a sensitivity axis in §8.3.

### §7.4 区室解析的 TACSTD2 / Compartment-Resolved TACSTD2 (E3)

中文 — 这是本协议最重要的一步，也是"高级去卷积"在这个问题上真正的增量价值：**不只要比例，还要区室特异表达**。

**路线 A — BayesPrism/InstaPrism 后验（首选）。** BayesPrism 的生成模型同时给出 $\theta$ 与细胞类型特异表达 $Z$。取恶性区室的 `TACSTD2`：

$$E^{\text{mal}}_{\texttt{TACSTD2},n} = \frac{Z_{\texttt{TACSTD2},\,\text{mal},\,n}}{\sum_{g} Z_{g,\,\text{mal},\,n}} \times 10^6$$

即在恶性区室内部重新归一化为 CPM，这一步不可省略——不归一化的 $Z$ 同时编码了区室大小，等于把纯度混杂又请回来了。

**路线 B — CIBERSORTx HiRes（次选/验证）。** HiRes 直接给出每样本每细胞类型的基因表达估计。取恶性/上皮类型的 `TACSTD2` 行。注意 HiRes 输出存在不可估计条目，必须记录每个样本的可估计状态；缺失不可用 0 填充（0 是一个有意义的表达值），应作为缺失进入模型或触发该样本在 E3 分析中的排除。

**一致性要求。** 路线 A 与路线 B 的 $E^{\text{mal}}_{\texttt{TACSTD2}}$ 必须先做相关比较，结果记入 `results/qc/e3_route_concordance.tsv`。若两路线不一致，E3 降级为探索性，且这一降级必须在报告中说明。

**可识别性告警。** 区室特异表达的可识别性依赖于该基因在区室间有足够的表达对比度、以及该区室在样本中有足够丰度。因此规程要求：(a) 预先设定恶性比例下限，低于该下限的样本在 E3 分析中排除（在极低纯度样本中，"每恶性细胞的表达"几乎不可估计）；(b) 检查 `TACSTD2` 在参考中的区室特异性——如果参考里免疫细胞也有非平凡表达，$Z$ 的分解会更不稳定，必须记录。

EN — This is the protocol's most important step and the genuine added value of "advanced" deconvolution for this question: **not just fractions, but compartment-specific expression**.

**Route A — BayesPrism/InstaPrism posterior (preferred).** BayesPrism's generative model yields $\theta$ and cell-type-specific expression $Z$ jointly. Take `TACSTD2` in the malignant compartment:

$$E^{\text{mal}}_{\texttt{TACSTD2},n} = \frac{Z_{\texttt{TACSTD2},\,\text{mal},\,n}}{\sum_{g} Z_{g,\,\text{mal},\,n}} \times 10^6$$

that is, renormalize to CPM *within* the malignant compartment. This step is not optional — unnormalized $Z$ still encodes compartment size, which reintroduces the purity confounder you just removed.

**Route B — CIBERSORTx HiRes (secondary / verification).** HiRes directly imputes per-sample, per-cell-type gene expression; take the `TACSTD2` row for the malignant/epithelial type. HiRes emits non-estimable entries, so record estimability per sample; **do not impute missing values with 0** (zero is a meaningful expression value) — carry them as missing or exclude the sample from the E3 analysis.

**Concordance requirement.** Correlate $E^{\text{mal}}_{\texttt{TACSTD2}}$ from Route A against Route B first and record it in `results/qc/e3_route_concordance.tsv`. If the routes disagree, E3 is demoted to exploratory and the demotion is stated in the report.

**Identifiability warning.** Compartment-specific expression is identifiable only when the gene has sufficient contrast between compartments and the compartment is sufficiently abundant in the sample. The procedure therefore requires: (a) a pre-specified floor on the malignant fraction, below which samples are excluded from E3 ("expression per malignant cell" is barely estimable at very low purity); and (b) a check of `TACSTD2` compartment specificity in the reference — if immune cells carry non-trivial expression there, the $Z$ decomposition is less stable and this must be recorded.

```r
# methods/deconv_advanced/lib/e3_compartment.R
# Route A: malignant-compartment CPM for a target gene
malignant_cpm <- function(Z_mal, gene = "TACSTD2") {   # Z_mal: genes x samples
  stopifnot(gene %in% rownames(Z_mal))
  cpm <- sweep(Z_mal, 2, colSums(Z_mal), "/") * 1e6
  log2(cpm[gene, ] + 1)
}

# Identifiability gate
e3_eligible <- function(theta, malignant_label, floor_frac) {
  theta[, malignant_label] >= floor_frac
}
```

### §7.5 成分侧的结局变量 / Compositional Outcome Variables

中文 — 结局侧按 §5.3 处理。本协议预先指定的结局如下表，**主终点两个，其余为次要或探索性**。预先指定的目的是防止在多种变换与多个细胞类型之间挑选。

EN — The outcome side follows §5.3. Pre-specified outcomes are in the table below: **two primary endpoints**, the rest secondary or exploratory. Pre-specification exists to prevent selecting among transforms and cell types after seeing results.

| Tier | Outcome | Definition |
|---|---|---|
| **Primary 1** | CD8 compartment | $\mathrm{clr}(\theta)_{\text{CD8\_T}}$ at L2 |
| **Primary 2** | TLS | primary TLS signature score from §6.1, scored by the primary scorer from §6.3 |
| Secondary | CD8 ratio endpoints | $\log(\theta_{\text{CD8\_T}} / \theta_{\text{CD4\_T}})$ and $\log(\theta_{\text{CD8\_T}} / \theta_{\text{mal}})$ |
| Secondary | B/plasma | $\mathrm{clr}(\theta)_{\text{B}}$, $\mathrm{clr}(\theta)_{\text{plasma}}$ |
| Secondary | TLS alternates | remaining signatures from §6.1, all scorers from §6.3 |
| Exploratory | L3 CD8 states | exhausted / effector-memory / tissue-resident CLR coordinates |
| Exploratory | EcoTyper | CD8-lineage cell-state abundances and parent ecotype abundance |
| Exploratory | L1 | coarse-level CLR coordinates, as the stability floor |

### §7.6 模型规格 / Model Specification

中文 — 单队列内的主模型：

$$\mathrm{clr}(\theta)_{\text{CD8},n} = \beta_0 + \beta_1 \cdot X_n + \beta_2 \cdot \pi_n + \gamma^\top C_n + \varepsilon_n$$

其中 $X_n$ 依估计量而定（E1: bulk `TACSTD2`；E2: bulk `TACSTD2`，且 $\pi_n$ 在场；E3: $E^{\text{mal}}_{\texttt{TACSTD2}}$）。协变量 $C_n$ 预先固定为：年龄、性别（若使用含 Y 连锁基因的签名则强制）、治疗线数、样本时间点（pre/on）、固定方式（FFPE/FF）、文库深度、以及 FFPE 质量指标。**协变量集在看到结果前冻结，不做逐步筛选。**

注意 E3 中 $\pi_n$ 的角色变化：既然 $X_n$ 已是区室内部的量，$\pi_n$ 不再是必需的混杂校正，而是作为**效应修饰候选**（低纯度样本的 $E^{\text{mal}}$ 估计噪声更大）。规程：E3 主模型仍包含 $\pi_n$ 作为协变量，并额外报告以 $1/\mathrm{Var}$ 或恶性比例为权重的加权拟合。

跨队列合并按 R0-5，两种可接受形式：

- **两阶段（首选）**：每队列独立拟合 → 提取 $\hat{\beta_1}$ 与其标准误 → `metafor::rma` 随机效应合并 → 报告合并估计、95% CI、$I^2$ 与 $\tau^2$。这一形式对队列间平台差异最稳健，且异质性本身是可解释的输出。
- **一阶段**：混合模型，队列为随机截距，必要时对 $X$ 的斜率也设随机效应（`lme4::lmer` 或 `glmmTMB` 的 beta 族）。仅在队列数较少、两阶段不稳定时使用。

**禁止**：把所有队列的表达矩阵拼接后当作单一数据集分析。

EN — Primary within-cohort model:

$$\mathrm{clr}(\theta)_{\text{CD8},n} = \beta_0 + \beta_1 \cdot X_n + \beta_2 \cdot \pi_n + \gamma^\top C_n + \varepsilon_n$$

with $X_n$ set by the estimand (E1: bulk `TACSTD2`; E2: bulk `TACSTD2` with $\pi_n$ present; E3: $E^{\text{mal}}_{\texttt{TACSTD2}}$). Covariates $C_n$ are fixed a priori as: age, sex (forced whenever a Y-linked-gene signature is used), line of therapy, sample timepoint (pre/on), fixation (FFPE/FF), library depth, and FFPE quality metrics. **The covariate set freezes before results are seen; no stepwise selection.**

Note how $\pi_n$'s role changes under E3: since $X_n$ is already a within-compartment quantity, $\pi_n$ is no longer a required confounder adjustment but a **candidate effect modifier** (low-purity samples have noisier $E^{\text{mal}}$). Procedure: the E3 primary model still includes $\pi_n$ as a covariate, and additionally reports a weighted fit with weights from $1/\mathrm{Var}$ or the malignant fraction.

Cross-cohort aggregation follows R0-5 in one of two acceptable forms:

- **Two-stage (preferred)**: fit each cohort independently → extract $\hat{\beta_1}$ and its standard error → pool with `metafor::rma` random effects → report the pooled estimate, 95% CI, $I^2$, and $\tau^2$. This form is most robust to cross-cohort platform differences, and the heterogeneity is itself an interpretable output.
- **One-stage**: a mixed model with cohort random intercepts, adding a random slope on $X$ where warranted (`lme4::lmer`, or `glmmTMB` with a beta family). Use only when the number of cohorts is small and the two-stage fit is unstable.

**Prohibited**: concatenating all cohort expression matrices and analyzing them as one dataset.

```r
# methods/deconv_advanced/analysis/associate.R
suppressPackageStartupMessages({ library(ppcor); library(metafor); library(betareg) })

fit_cohort <- function(df, estimand = c("E1", "E2", "E3"), covars) {
  estimand <- match.arg(estimand)
  x <- switch(estimand, E1 = "tacstd2_bulk", E2 = "tacstd2_bulk", E3 = "tacstd2_malignant")
  rhs <- switch(estimand,
                E1 = c(x, covars),                       # purity deliberately excluded
                E2 = c(x, "purity", covars),
                E3 = c(x, "purity", covars))             # purity as modifier, see 7.6
  f <- reformulate(rhs, response = "clr_cd8")
  m <- lm(f, data = df)
  s <- summary(m)$coefficients[x, ]
  data.frame(cohort = df$cohort[1], estimand = estimand,
             beta = s["Estimate"], se = s["Std. Error"], n = nrow(df))
}

# Rank-based companion (reported alongside every parametric fit)
partial_spearman <- function(df, x, y, z) {
  ppcor::pcor.test(df[[x]], df[[y]], df[[z]], method = "spearman")
}

pool <- function(per_cohort) {
  metafor::rma(yi = per_cohort$beta, sei = per_cohort$se, method = "REML")
}
```

### §7.7 治疗语境分层 / Treatment-Context Stratification

中文 — ICI 队列的异质性不是噪声，而是结构。以下分层在分析计划中预先声明，每层单独估计后再决定是否合并（由 §7.6 的 $I^2$ 指导）：

1. **时间点**：pre-treatment 与 on-treatment 必须分开。on-treatment 样本的免疫组成已被治疗改变，与 pre-treatment 混合会同时稀释和污染估计。
2. **药物类别**：抗 PD-1/PD-L1 与抗 CTLA-4 分开；联合治疗单列。
3. **癌种**：不同癌种分别估计。上皮结构与免疫组织在癌种间差异极大，而 `TACSTD2` 恰是上皮基因。跨癌种合并只在异质性可接受时进行。
4. **既往治疗**：化疗/放疗暴露记录为协变量。
5. **配对样本**：若队列含同一患者的 pre/on 配对，配对分析（患者随机截距）作为预先声明的次要分析，而非把配对样本当独立观测。

**与结局的交互属于另一个问题。** 若要问"`TACSTD2` 与 CD8 的关联是否随 ICI 响应而不同"，这是一个交互项分析（$X \times \text{response}$），其功效远低于主效应，必须预先声明为探索性并按此报告。绝不能在主分析失败后转而报告某个亚组的交互。

EN — Heterogeneity across ICI cohorts is structure, not noise. Declare these strata in the analysis plan, estimate within each, and only then decide on pooling (guided by $I^2$ from §7.6):

1. **Timepoint.** Pre-treatment and on-treatment must be separated. On-treatment immune composition has already been altered by therapy; mixing it with pre-treatment both dilutes and contaminates the estimate.
2. **Agent class.** Anti-PD-1/PD-L1 separate from anti-CTLA-4; combinations reported separately.
3. **Tumor type.** Estimate per tumor type. Epithelial architecture and immune organization differ sharply across tumor types, and `TACSTD2` is precisely an epithelial gene. Pool across types only when heterogeneity is acceptable.
4. **Prior therapy.** Record chemotherapy/radiotherapy exposure as a covariate.
5. **Paired samples.** Where a cohort contains pre/on pairs from the same patient, a paired analysis (patient random intercept) is a pre-declared secondary analysis — never treat paired samples as independent observations.

**Interaction with outcome is a different question.** Asking "does the `TACSTD2`–CD8 association differ by ICI response" is an interaction analysis ($X \times \text{response}$), far less powered than the main effect, and must be pre-declared as exploratory and reported as such. Never pivot to a subgroup interaction after the primary analysis fails.

### §7.8 对照与多重比较 / Controls and Multiplicity

中文 — R0-4 的实现。每次关联分析必须同时产出：

**阳性对照（验证管线可工作）**
- `CD8A` 与 `CD8B` 的 bulk 表达 vs CD8 比例：若这个相关不强，去卷积或本体映射有问题，整条管线在该队列上不可用。
- `MS4A1` vs B 细胞比例；`PTPRC` vs 免疫总比例；`COL1A1` vs 成纤维细胞比例。
- 影像可得时（§4.5）：影像 CD8 密度/占比 vs 去卷积 CD8 比例。

**阴性对照（估计经验零分布）**
- **基因侧经验零**：构造一个匹配的基因集合——与 `TACSTD2` 表达水平相近、表达方差相近、且与纯度相关性相近的 $\geq 500$ 个基因——对每个基因跑相同的模型，得到 $\beta_1$ 的经验零分布。`TACSTD2` 的统计量相对于该零分布定位，比相对于理论零假设定位更有意义，因为它自动吸收了"任何上皮基因都会有的"结构性相关。这是本节最重要的一条。
- **管家基因对照**：若干稳定表达基因，预期无关联。
- **置换检验**：在队列内置换样本标签，重估。

**多重比较**
- 主终点仅两个（CD8、TLS），预先设定 Bonferroni 或按 $\alpha/2$ 分配。
- 次要与探索性终点：在各自族内 BH-FDR 控制，且**报告时必须标注所属层级**。
- 跨方法（§8）不做 FDR 校正：方法间比较的目的是稳定性评估而非发现，把它当作多重发现处理是概念错误。

EN — This implements R0-4. Every association analysis emits, in the same run:

**Positive controls (pipeline works at all)**
- Bulk `CD8A` and `CD8B` versus the CD8 fraction: if this correlation is not strong, deconvolution or ontology mapping is broken and the whole pipeline is unusable in that cohort.
- `MS4A1` versus B fraction; `PTPRC` versus total immune fraction; `COL1A1` versus fibroblast fraction.
- Where imaging exists (§4.5): image CD8 density/share versus deconvolved CD8 fraction.

**Negative controls (empirical null)**
- **Gene-side empirical null**: assemble a matched gene set — $\geq 500$ genes matched to `TACSTD2` on mean expression, expression variance, and correlation with purity — and run the identical model for each, yielding an empirical null distribution of $\beta_1$. Positioning the `TACSTD2` statistic against *that* null is far more informative than against the theoretical null, because it automatically absorbs the structural correlation that **any** epithelial gene exhibits. This is the single most important item in this section.
- **Housekeeping controls**: stably expressed genes with no expected association.
- **Permutation**: shuffle sample labels within cohort and re-estimate.

**Multiplicity**
- Only two primary endpoints (CD8, TLS); pre-specify Bonferroni or an $\alpha/2$ split.
- Secondary and exploratory endpoints: BH-FDR within their own family, and **every reported value must be labeled with its tier**.
- No FDR correction across methods (§8): cross-method comparison assesses stability, not discovery, and treating it as multiple discovery is a conceptual error.

### §7.9 分析计划摘要 / Analysis Plan Summary

中文 — 一页纸版本，冻结后存为 `methods/deconv_advanced/analysis_plan_frozen.md` 并记录 git commit hash：

EN — The one-page version; freeze it as `methods/deconv_advanced/analysis_plan_frozen.md` with a recorded git commit hash:

| Item | Pre-specified choice |
|---|---|
| Primary estimand | E3 (malignant-compartment `TACSTD2`), co-primary E2 |
| Mandatory descriptive control | E1 (unadjusted bulk), always reported (R0-2) |
| Primary deconvolution method | Frozen in §8.0 before any association is computed |
| Primary granularity | L2 |
| Primary endpoints | $\mathrm{clr}(\theta)_{\text{CD8\_T}}$ ; primary TLS score |
| Primary purity source | Highest-priority available per §7.3 |
| Compositional handling | CLR with `cmultRepl` zero replacement |
| Cross-cohort | Two-stage random-effects meta-analysis |
| Stratification | Timepoint, agent class, tumor type (§7.7) |
| Controls | Positive + gene-matched empirical null (§7.8) |
| Uncertainty | Reference-resampling bootstrap, $B \geq 200$ (§5.4) |
| Success criterion for a stable claim | Consistent direction across methods per §8.4 |
| Two-seed / two-method contrast | `bayesprism_ecotyper_vs_tacstd2_cldn4.md` (E3-T, E3-C, E4–E8; patterns P-agree / P-gene-split / P-method-split / P-null) |

---

## §8 方法比较协议 / Comparison Protocol

中文 — 比较协议的目的不是评"哪个方法最好"这种脱离语境的排名，而是回答两个具体问题：**(Q1) 在我的数据条件下，哪个方法的细胞比例估计最可信？(Q2) 我关于 `TACSTD2` 的结论是否依赖于方法选择？** Q2 才是决定性的：如果结论在方法间翻转，那么真正的发现是"该结论不可判定"，而不是任何一个方法给出的方向。

EN — The comparison protocol does not aim at a context-free ranking of "which method is best." It answers two concrete questions: **(Q1) under my data conditions, whose cell-fraction estimates are most credible? (Q2) does my `TACSTD2` conclusion depend on the method choice?** Q2 is decisive: if the conclusion flips across methods, the real finding is that the conclusion is undecidable, not whichever direction one method produced.

### §8.0 预注册 / Pre-Registration

中文 — 在运行任何关联分析之前完成，写入 `methods/deconv_advanced/prereg.md`：

1. 声明主去卷积方法与理由。选择依据只能来自 §8.1–§8.3 的基准结果（这些不涉及 `TACSTD2`）与文献，不得来自任何关联结果。
2. 声明比较中纳入的方法臂与粒度级别。
3. 声明 §8.2 的主指标与 §8.4 的稳定性判据。
4. 声明所有敏感性轴（§8.3）及其取值。
5. 提交 git 并记录 commit hash。此后对计划的任何修改必须以新 commit 追加，并标注为"事后修改"。

EN — Complete before any association analysis runs; write it to `methods/deconv_advanced/prereg.md`:

1. Declare the primary deconvolution method and the rationale. The rationale may draw only on §8.1–§8.3 benchmark results (which never touch `TACSTD2`) and on the literature — never on an association result.
2. Declare the method arms and granularity levels in the comparison.
3. Declare the primary metric from §8.2 and the stability criterion from §8.4.
4. Declare every sensitivity axis (§8.3) and its levels.
5. Commit to git and record the hash. Any later change is appended as a new commit labeled "post hoc amendment."

### §8.1 金标准来源 / Ground-Truth Sources

中文 — 按可信度与可得性排序，尽可能用多个：

| Tier | Source | Notes |
|---|---|---|
| **G1** | Matched flow cytometry / CyTOF on the same specimens | Strongest, rarely available for archival ICI cohorts |
| **G2** | Matched mIF / IHC quantified by image analysis (§4.5, InstantDL Role A) | The realistic strong option; produces **densities**, so convert to share-of-nucleated-cells before comparing to fractions |
| **G3** | Matched scRNA-seq on the same specimens, with fractions computed from cell counts | Beware capture-efficiency bias per cell type — this is a biased ground truth, not a neutral one |
| **G4** | **Heterogeneous pseudobulk simulation** | The default workhorse; see the critical note below |
| **G5** | Purified/sorted bulk profiles run through the pipeline | The only clean way to measure spillover (§8.2) |
| **G6** | Spike-in / defined cell mixtures | Clean but distant from tumor tissue |

> 中文 — **关于 G4 的关键点。** 传统伪 bulk 通过从单细胞池里**随机抽取细胞**再求和构造混合物。2024 年的方法学工作明确指出这样构造的伪 bulk **缺乏真实 bulk 应有的生物学方差**：随机抽样把供体间与细胞状态间的异质性平均掉了，于是几乎所有方法在这种简单基准上都表现优异，基准失去区分力，并系统性高估参考匹配良好时的性能。规程要求使用**异质性伪 bulk**：按供体分层抽样、在细胞类型内引入状态组成的样本间变异、并使模拟的方差结构匹配真实 bulk 的观测方差。`SimBu` 一类工具支持这类设计（含 mRNA-bias 缩放）。**用随机取细胞的简单伪 bulk 作为唯一基准，本身即构成方法学缺陷。**
>
> EN — **The critical point about G4.** Conventional pseudobulk sums **randomly sampled cells** from a single-cell pool. Methodological work published in 2024 showed such pseudobulk **lacks the biological variance real bulk exhibits**: random sampling averages away donor-to-donor and cell-state heterogeneity, so nearly every method looks excellent on this easy benchmark, the benchmark loses discriminative power, and performance is systematically overstated in the well-matched-reference regime. The procedure therefore requires **heterogeneous pseudobulk**: stratify sampling by donor, inject between-sample variation in within-cell-type state composition, and match the simulated variance structure to variance observed in real bulk. Tools such as `SimBu` support this design (including mRNA-bias scaling). **Using simple random-cell pseudobulk as the only benchmark is itself a methodological defect.**

中文 — 模拟设计矩阵（每个格子生成预先设定数量的样本）：肿瘤纯度网格（覆盖低到高，因为多项基准显示纯度升高时正常上皮易被误判为恶性上皮，而这直接影响 `TACSTD2` 的区室归属）× 免疫浸润水平 × 缺失细胞类型情形（从混合物中移除一个参考里有的类型，反之亦然）× 参考-混合物平台错配（10x 参考 vs Smart-seq2 参考；scRNA vs snRNA）。

EN — Simulation design matrix (each cell generates a pre-specified number of samples): a tumor-purity grid spanning low to high (multiple benchmarks report that normal epithelium is increasingly mis-assigned to malignant epithelium as purity rises, which directly affects `TACSTD2` compartment assignment) × infiltration level × missing-cell-type scenarios (remove a type from the mixture that is present in the reference, and vice versa) × reference–mixture platform mismatch (10x versus Smart-seq2 reference; scRNA versus snRNA).

### §8.2 指标 / Metrics

中文 — 指标必须分开报告两个方向，混在一起是常见错误：

- **按细胞类型（跨样本）**：对每个细胞类型，估计值与真值在样本间的相关。这回答"能否追踪某细胞类型在患者间的变化"，**这才是 §7 关联分析所依赖的性质**。
- **按样本（跨细胞类型）**：单个样本内各细胞类型的估计与真值的相关。这回答"能否刻画一个样本的组成谱"。一个方法可以在后者上很好而在前者上无用（例如系统性压缩了某类型的样本间方差）。

EN — Report both directions separately; conflating them is a common error:

- **Per cell type (across samples)**: for each cell type, correlate estimate against truth across samples. This answers "can it track a cell type's variation between patients," and **it is the property the §7 association depends on**.
- **Per sample (across cell types)**: within one sample, correlate the estimated composition profile against truth. This answers "can it describe one sample's composition." A method can excel at the latter and be useless for the former (e.g. by systematically compressing between-sample variance for a type).

| Metric | Definition | Why included |
|---|---|---|
| Pearson $r$ | per cell type, across samples | Comparability with published benchmarks |
| Spearman $\rho$ | per cell type, across samples | Rank tracking; robust to calibration error |
| Lin's CCC | concordance correlation | Penalizes bias and scale error that $r$ ignores |
| RMSE / MAE | absolute error | Interpretable on the fraction scale |
| **Calibration slope** | regress truth on estimate | Detects systematic compression/inflation of between-sample range — the failure that silently attenuates §7 |
| **Spillover matrix** | run pure/sorted profiles (G5) through the pipeline; record the fraction assigned to each non-true type | The decisive metric for CD8: the recurring failure modes are CD4↔Treg confusion and NK→CD8 misassignment |
| Missing-type robustness | error induced when a mixture type is absent from the reference | Real references are always incomplete |
| Rank stability | rank correlation of per-method rankings across simulation scenarios | Whether a "winner" is scenario-specific |
| Runtime / memory | wall-clock and peak RSS at cohort scale | Governs whether bootstrap (§5.4) is affordable |

中文 — 主指标（在 §8.0 预先声明）建议为 **L2 `CD8_T` 的跨样本 Spearman $\rho$ 与该类型的溢出总量**，因为主终点就定义在这个量上。不要用"所有细胞类型的平均性能"作主指标——平均值会被丰度高、易估的类型（如恶性细胞）主导，而它与 CD8 的可估性无关。

EN — The recommended primary metric (declared in §8.0) is **across-sample Spearman $\rho$ for L2 `CD8_T` together with that type's total spillover**, because the primary endpoint is defined on exactly that quantity. Do not use "mean performance across all cell types" as the primary metric — the mean is dominated by abundant, easy types such as malignant cells and says nothing about CD8 estimability.

### §8.3 敏感性轴 / Sensitivity Axes

中文 — 完全交叉运行成本很高；规程为：先在**一个代表性队列**上跑全交叉，识别出影响最大的 2–3 条轴，再在其余队列上只跑这几条。所有轴与取值在 §8.0 预先声明。

EN — A full cross of all axes is expensive; the procedure is: run the full cross on **one representative cohort**, identify the 2–3 highest-impact axes, then run only those on the remaining cohorts. All axes and levels are declared in §8.0.

| Axis | Levels |
|---|---|
| Method | InstaPrism, CIBERSORTx, MuSiC, MuSiC2, (Scaden as DL arm) |
| Granularity | L1, L2, L3 |
| Reference | primary atlas, alternate atlas, down-sampled reference, donor-leave-one-out |
| Reference chemistry | scRNA vs snRNA (triggers the §2.3 `TACSTD2` intronless caveat) |
| Marker/feature selection | method default vs a fixed shared gene panel |
| Batch mode (CIBERSORTx) | S-mode vs B-mode |
| Reference update (InstaPrism) | `update = TRUE` vs `FALSE` |
| Fraction scale | relative vs absolute |
| Purity source | DNA / methylation / ESTIMATE / deconvolution-derived (§7.3) |
| Compositional transform | CLR vs ratio endpoints vs beta regression |
| TLS signature | all of §6.1 |
| TLS scorer | ssGSEA vs singscore vs z-mean |
| E3 route | BayesPrism $Z$ vs CIBERSORTx HiRes |
| Sample filters | malignant-fraction floor on/off; FFPE-only vs all |

### §8.4 结论稳定性：决定性终点 / Conclusion Stability: The Decisive Endpoint

中文 — 这是整个比较协议的落点。对每条敏感性轴的每个取值，重跑 §7 的主分析，得到一组 $\hat{\beta_1}$。然后按下表分级，**分级规则在 §8.0 预先声明，不得事后调整**：

EN — This is where the comparison protocol lands. For every level of every sensitivity axis, rerun the §7 primary analysis to obtain a set of $\hat{\beta_1}$ values, then grade with the table below. **The grading rule is declared in §8.0 and may not be adjusted afterwards**:

| Grade | Criterion | Permitted claim |
|---|---|---|
| **S1 — Robust** | Same sign across **all** method arms and **all** pre-declared sensitivity levels; confidence intervals excluding the null in the majority of arms; consistent direction in the meta-analysis with acceptable $I^2$ | A directional association may be stated, with the estimand (E1/E2/E3) named explicitly |
| **S2 — Conditional** | Same sign in the primary arm and most arms, but flipping or nulling on an identifiable axis (e.g. only at L3, or only under S-mode) | The claim must be stated **with its condition attached**, and the responsible axis named |
| **S3 — Unstable** | Sign flips across method arms with no identifiable governing axis | The permitted claim is that **the association is not determinable with these data and methods**. This is a legitimate, reportable methodological outcome, not a failure |
| **S4 — Artifact-dominated** | The `TACSTD2` statistic sits inside the gene-matched empirical null (§7.8), or the positive controls fail | No association claim at all; report the pipeline diagnostic instead |

中文 — 报告要求：稳定性分级必须与效应量同时出现在任何摘要句中。写"`TACSTD2` 与 CD8 负相关"是不完整的；完整形式是"在估计量 E3、主方法 X、稳定性 S2（在 L3 粒度下翻转）的条件下，观察到负向关联"。这条要求的目的是让读者能够判断结论的适用边界，而这正是当前该领域文献中方向性矛盾难以调和的原因之一。

EN — Reporting requirement: the stability grade must appear alongside the effect size in any summary sentence. Writing "`TACSTD2` is negatively associated with CD8" is incomplete; the complete form is "under estimand E3, primary method X, stability S2 (sign flips at L3 granularity), a negative association was observed." The point of this requirement is to let a reader see the boundary of applicability — one of the reasons directional contradictions in this literature are currently hard to reconcile.

### §8.5 共识与集成 / Consensus and Ensembles

中文 — 允许构建共识估计，但有严格约束：

1. 共识**不能替代**预先声明的主方法。主分析仍以单一冻结方法为准；共识作为预先声明的次要分析。
2. 共识的构造方式必须在 §8.0 固定：可接受的形式为 CLR 空间内的逐细胞类型中位数，或以 §8.2 主指标为权重的加权平均。**禁止**在看到关联结果后调整权重。
3. 报告共识时必须同时报告成员方法间的离散度；离散度大的共识是 S3 的伪装，不是 S1。
4. 集成不能修复系统性偏倚：若所有成员方法都存在同向的 CD4→CD8 溢出，共识只会更自信地重复该偏倚。因此共识必须与 §8.2 的溢出矩阵一起解读。

EN — Consensus estimates are permitted under strict constraints:

1. Consensus **does not replace** the pre-declared primary method. The primary analysis stays with one frozen method; consensus is a pre-declared secondary analysis.
2. The consensus construction is fixed in §8.0: acceptable forms are a per-cell-type median in CLR space, or a weighted mean with weights from the §8.2 primary metric. **Reweighting after seeing association results is prohibited.**
3. When reporting consensus, always report dispersion across member methods; a high-dispersion consensus is S3 in disguise, not S1.
4. Ensembling cannot repair systematic bias: if every member method carries the same directional CD4→CD8 spillover, the consensus merely repeats that bias with more confidence. Consensus must therefore be read together with the §8.2 spillover matrix.

### §8.6 比较协议执行顺序 / Execution Order

```
1. Build and QC references (§3), including the pseudobulk self-test (§3.3)
2. Run all method arms on simulated ground truth (§8.1 G4/G5) across the design matrix
3. Compute metrics (§8.2); produce the spillover matrix
4. Where imaging exists, calibrate against G2 (§4.5, §8.1)
5. FREEZE: choose the primary method, write prereg.md, commit (§8.0)
   ---- no association analysis has run up to this point ----
6. Run deconvolution on real cohorts, all arms (§4)
7. Harmonize outputs (§5); bootstrap for uncertainty (§5.4)
8. Score TLS (§6)
9. Run §7 primary analysis with the frozen method; run controls (§7.8)
10. Run §7 across all sensitivity axes (§8.3)
11. Grade conclusion stability (§8.4)
12. Report per §9
```

中文 — 第 5 步的冻结点是整个协议的完整性所在：在这一点之前，没有任何 `TACSTD2` 与免疫读数的关联被计算过，因此主方法的选择在结构上不可能被结果污染。执行时应有机制保证这一点（例如关联脚本在 `prereg.md` 存在且含有效 commit hash 之前拒绝运行）。

EN — The freeze at step 5 is where the protocol's integrity lives: before that point, no `TACSTD2`–immune association has been computed, so the primary-method choice is structurally incapable of being contaminated by results. Enforce this mechanically (for example, have the association script refuse to run until `prereg.md` exists and contains a valid commit hash).

---

## §9 可复现性与报告 / Reproducibility and Reporting

### §9.1 环境 / Environment

| Artifact | Content |
|---|---|
| `env/versions.lock` | Exact versions of R, Python, and every deconvolution package; container digests |
| `env/Dockerfile` | Pinned base image; CIBERSORTx runs from its own official image, referenced by digest |
| `env/LICENSES.md` | License terms per tool, notably CIBERSORTx's restricted terms |
| `env/seeds.yaml` | Every random seed: deconvolution, bootstrap, simulation, permutation |
| `env/hardware.md` | CPU/RAM/GPU, since GPU nondeterminism affects any DL arm |

### §9.2 运行清单 / Run Manifest

中文 — 每次运行写一条 `results/manifests/{run_id}.json`，包含：git commit hash、配置文件内容哈希、输入文件哈希、方法与版本、参数全集、随机种子、运行时长、输出文件哈希。**没有 manifest 的输出视为不存在**，不得进入报告。

EN — Every run writes `results/manifests/{run_id}.json` containing: git commit hash, config content hash, input file hashes, method and version, the full parameter set, random seeds, wall-clock duration, and output file hashes. **Output without a manifest does not exist** and may not enter the report.

### §9.3 报告检查表 / Reporting Checklist

中文 — 提交前逐项核对；每项在报告中都要能定位到具体位置：

EN — Verify item by item before submission; each must be locatable in the report:

- [ ] Cohort registry with n, platform, fixation, timepoint, agent class (§2.1)
- [ ] Reference provenance: donors, chemistry, malignant-annotation method, cell counts per type (§3.2)
- [ ] Reference self-test result (§3.3)
- [ ] Granularity levels analyzed, with L2 identified as primary (§3.1)
- [ ] Full parameters for every method, including batch mode and units (§4)
- [ ] Ontology map version (§5.1)
- [ ] Fraction scale (relative/absolute) and compositional transform (§5.2, §5.3)
- [ ] Uncertainty procedure and $B$ (§5.4)
- [ ] TLS signature list, primary selection rationale, scorer, coverage (§6)
- [ ] `EIF1AY` handling declared (§6.2)
- [ ] All three estimands E1/E2/E3 reported (§7.2, R0-2)
- [ ] Purity source and its circularity tier (§7.3)
- [ ] E3 route concordance and identifiability gate (§7.4)
- [ ] Model formula, covariates, aggregation method, $I^2$/$\tau^2$ (§7.6)
- [ ] Stratification results (§7.7)
- [ ] Positive controls and gene-matched empirical null (§7.8)
- [ ] Benchmark design, metrics, spillover matrix (§8.1, §8.2)
- [ ] Sensitivity axes actually run (§8.3)
- [ ] **Conclusion stability grade S1–S4 attached to every claim (§8.4)**
- [ ] `prereg.md` commit hash and any post hoc amendments (§8.0)
- [ ] Manifests for every reported figure and table (§9.2)
- [ ] InstantDL scope note wherever InstantDL is mentioned (§1.5)

---

## §10 目录结构 / Directory Layout

```
methods/deconv_advanced/
├── playbook.md                     # this file (methods only)
├── bayesprism_ecotyper_vs_tacstd2_cldn4.md   # §13 companion (two methods × two seeds)
├── README.md
├── prereg.md                       # §8.0, written at the freeze point
├── analysis_plan_frozen.md         # §7.9
├── imaging_review_sop.md           # §4.5 pathologist review SOP
├── config/
│   └── barrier_genes.yaml          # locked BARRIER5 contract
├── run/
│   ├── run_instaprism.R            # §4.1
│   ├── run_cibersortx.sh           # §4.2
│   ├── run_music.R                 # §4.3
│   ├── run_ecotyper.sh             # §4.4
│   ├── scale_for_ecotyper.R        # §4.4 / §B5.1 auditable scaling
│   └── run_instantdl.sh            # §4.5 (imaging, Role A)
├── lib/
│   ├── compositional.R             # §5.3
│   ├── tls_score.R                 # §6.3
│   ├── e3_compartment.R            # §7.4
│   └── controls.R                  # §7.8
├── analysis/
│   ├── associate.R                 # §7.6
│   ├── e4_e5_seeds.R               # §B4.4 TACSTD2 vs CLDN4 residuals
│   ├── bp_vs_ecotyper_pattern.R    # §B6.2 pattern classifier
│   ├── sensitivity.R               # §8.3
│   └── stability_grade.R           # §8.4
└── bench/
    ├── simulate_heterogeneous.R    # §8.1 G4
    ├── spillover.R                 # §8.2
    └── metrics.R                   # §8.2

config/                             # every parameter, no defaults relied upon
data/                               # cohorts, references, purity, imaging
results/                            # ALL numeric output lives here, never in methods/
```

中文 — 硬性分离：`methods/` 只放方法与代码，`results/` 只放数值产出。本手册中不出现任何 `results/` 的内容，这正是文首"METHODS ONLY"声明的实际含义。

EN — Hard separation: `methods/` holds methods and code only; `results/` holds numeric output only. Nothing from `results/` appears in this playbook, which is what the "METHODS ONLY" banner at the top actually means in practice.

---

## §11 陷阱清单 / Pitfall Register

| # | Pitfall | 中文说明 | Detection / Prevention |
|---|---|---|---|
| P1 | Naive correlation of bulk `TACSTD2` with immune fractions | 纯度混杂足以单独造出可观的相关 | §7.2 三估计量；§7.3 纯度校正 |
| P2 | Treating fractions as Euclidean | 单纯形上的数据用 Pearson/OLS | §5.3 CLR / ratio / beta |
| P3 | Adjusting for deconvolution-derived purity | 与结局同源，诱导偏倚 | §7.3 优先 DNA/甲基化来源 |
| P4 | Simple random-cell pseudobulk as the only benchmark | 缺乏真实方差，所有方法都"很好" | §8.1 异质性模拟 |
| P5 | Wrong units into a method | TPM 喂给需要 counts 的方法 | §2.2 单位契约 + 断言 |
| P6 | CIBERSORTx B-mode/S-mode confusion | 选错批次模式 | §4.2；在 manifest 中记录 |
| P7 | Skipping EcoTyper per-dataset scaling | 静默失败，仍输出数值 | §4.4 前置步骤为独立可审计脚本 |
| P8 | `EIF1AY` in the Cabrita TLS signature | 性别当成 TLS | §6.2 强制处理 |
| P9 | CD8/TLS circularity via CXCL9-11, PDCD1 | 同一批基因被测两次 | §6.2；两个独立主终点，不做复合评分 |
| P10 | Cross-cohort ComBat before deconvolution | 破坏线性混合尺度 | §2.3 点 4 |
| P11 | Naive pooling across cohorts | 忽略平台与人群异质性 | R0-5；§7.6 两阶段元分析 |
| P12 | Method chosen after seeing associations | 结果驱动的方法选择 | §8.0 冻结点 + 脚本级强制 |
| P13 | CD4↔Treg / NK→CD8 spillover | CD8 主终点被污染 | §8.2 溢出矩阵；L2 而非 L3 作主终点 |
| P14 | HiRes missing values imputed as 0 | 0 是有效表达值 | §7.4 路线 B |
| P15 | Unnormalized $Z$ used for E3 | 区室大小混入"每细胞表达" | §7.4 区室内 CPM 归一化 |
| P16 | Point-estimate fractions treated as known | 标准误被低估 | §5.4 bootstrap |
| P17 | On-treatment mixed with pre-treatment | 治疗已改变免疫组成 | §7.7 点 1 |
| P18 | EcoTyper state labels reused across runs | 标签无跨运行语义 | §1.4 点 b；绑定模型版本 |
| P19 | InstantDL treated as a deconvolution arm | 范畴错误 | §1.5 范围更正 |
| P20 | Paired pre/on samples treated as independent | 有效样本量被高估 | §7.7 点 5 |
| P21 | MuSiC2 `control`/`case` set to responder status | 结局信息泄漏进比例估计 | §1.3 点 c |
| P22 | snRNA reference for an intronless gene | `TACSTD2` 参考谱可能系统性偏移 | §2.3 点 3；§8.3 化学轴 |
| P23 | Reporting an effect without a stability grade | 读者无法判断适用边界 | §8.4 报告要求 |
| P24 | Interaction with response reported after a null main effect | 事后转向 | §7.7 预先声明为探索性 |
| P25 | Collapsing `TACSTD2` + `CLDN4` before concordance | 用复合分数掩盖单基因分歧 | companion §B1, R-B1 |
| P26 | Treating `CLDN4` as the TROP2 binding partner | 文献互作对象是 `CLDN7` | companion §B1.1 |
| P27 | Scoring BARRIER5 on bulk and calling it E6 | 纯度被重新引入"区室"分数 | companion §B4.5 |
| P28 | EcoTyper CE argmax as the E7 outcome | 丢掉丰度，在单纯形上赢者通吃 | companion §B5.3 |
| P29 | Carcinoma EcoTyper recovery outside supported histologies | CS/CE 标签不可转移 | companion §B2 |
| P30 | BayesPrism $\theta_{\text{CD8}}$ treated as EcoTyper CE9 | 比例 ≠ 群落 | companion §B6.2 |
| P31 | One empirical null shared by both seeds | 纯度相关结构不同，零假设错位 | companion §B6.5 |
| P32 | Promoting the hottest epithelial CS from E8 | 事后替换主终点 | companion §B6.3 |

---

## §12 参考文献 / References

中文 — 引用用于定位方法出处，不用于支持任何结果性陈述。使用前请核对版本与最新修订。

EN — Citations locate method provenance and support no result-level statement. Verify versions and latest revisions before use.

**Deconvolution methods**

1. Chu T, Wang Z, Pe'er D, Danko CG. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. *Nature Cancer* 3:505–517 (2022).
2. Hu M, Chikina M. InstaPrism: an R package for fast implementation of BayesPrism. *Bioinformatics* 40:btae440 (2024). DOI: 10.1093/bioinformatics/btae440
3. Newman AM, Steen CB, Liu CL, et al. Determining cell type abundance and expression from bulk tissues with digital cytometry (CIBERSORTx). *Nature Biotechnology* 37:773–782 (2019).
4. Wang X, Park J, Susztak K, Zhang NR, Li M. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference (MuSiC). *Nature Communications* 10:380 (2019).
5. Fan J, Lyu Y, Zhang Q, Wang X, Xiao R, Li M. MuSiC2: cell-type deconvolution for multi-condition bulk RNA-seq data. *Briefings in Bioinformatics* 23:bbac430 (2022).
6. Luca BA, Steen CB, Matusiak M, et al. Atlas of clinically distinct cell states and ecosystems across human solid tumors (EcoTyper). *Cell* 184:5482–5496 (2021).
7. Steen CB, Luca BA, Esfahani MS, et al. The landscape of tumor cell states and ecosystems in diffuse large B cell lymphoma. *Cancer Cell* 39:1422–1437 (2021).
8. Menden K, Marouf M, Oller S, et al. Deep learning-based cell composition analysis from tissue expression profiles (Scaden). *Science Advances* 6:eaba2619 (2020).
9. Chen Y, Wang S, et al. Deep autoencoder for interpretable tissue-adaptive deconvolution and cell-type-specific gene analysis (TAPE). *Nature Communications* 13:6735 (2022).

**Imaging pipeline (scope-corrected, §1.5)**

10. Waibel DJE, Shetab Boushehri S, Marr C. InstantDL: an easy-to-use deep learning pipeline for image segmentation and classification. *BMC Bioinformatics* 22:103 (2021). DOI: 10.1186/s12859-021-04037-3

**Benchmarking and simulation**

11. Hu M, Chikina M. Heterogeneous pseudobulk simulation enables realistic benchmarking of cell-type deconvolution methods. *Genome Biology* 25:127 (2024). DOI: 10.1186/s13059-024-03292-w
12. Dietrich A, Merotto L, Pelz K, et al. omnideconv: a unifying framework for using and benchmarking single-cell-informed deconvolution of bulk RNA-seq data. (bioRxiv 2024.06.10.598226; *Genome Biology*, PMC12837286)
13. Dietrich A, Sturm G, Merotto L, et al. SimBu: bias-aware simulation of bulk RNA-seq data with variable cell-type composition. *Bioinformatics* 38:ii141–ii147 (2022).
14. Sutton GJ, Poppe D, Simmons RK, et al. Comprehensive evaluation of deconvolution methods for human brain gene expression. *Nature Communications* 13:1358 (2022).
15. Performance of tumour microenvironment deconvolution methods in breast cancer using single-cell simulated bulk mixtures. *Nature Communications* 14:5678 (2023). DOI: 10.1038/s41467-023-41385-5
16. Hippen AA, Falco MM, Weber LM, et al. Performance of computational algorithms to deconvolve heterogeneous bulk ovarian tumor tissue depends on experimental factors. *Genome Biology* 24:239 (2023).
17. A robust workflow to benchmark deconvolution of multi-omic data. *Genome Biology* (2025). DOI: 10.1186/s13059-025-03897-9
18. Evaluating deconvolution methods using real bulk RNA-expression data for robust prognostic insights across cancer types. *Genome Biology* (2026). DOI: 10.1186/s13059-026-03942-1
19. DALE-Eval: a comprehensive cell-type-specific expression deconvolution benchmark for transcriptomics data. bioRxiv 2025.07.31.667984.

**TLS signatures and quantification**

20. Coppola D, Nebozhyn M, Khalil F, et al. Unique ectopic lymph node-like structures present in human primary colorectal carcinoma are identified by immune gene array profiling (12-chemokine signature). *American Journal of Pathology* 179:37–45 (2011).
21. Cabrita R, Lauss M, Sanna A, et al. Tertiary lymphoid structures improve immunotherapy and survival in melanoma. *Nature* 577:561–565 (2020).
22. Meylan M, Petitprez F, Becht E, et al. Tertiary lymphoid structures generate and propagate anti-tumor antibody-producing plasma cells in renal cell cancer. *Immunity* 55:527–541 (2022).
23. Gu-Trantien C, Loi S, Garaud S, et al. CD4⁺ follicular helper T cell infiltration predicts breast cancer survival. *Journal of Clinical Investigation* 123:2873–2892 (2013).
24. The 12-CK score: global measurement of tertiary lymphoid structures. *Frontiers in Immunology* 12:694079 (2021).
25. Comprehensive evaluation of gene expression-based signatures for detecting TLS in metastatic renal cell carcinoma. *Journal of Clinical Oncology* 43(5_suppl):591 (2025).

**TACSTD2 / TROP2 context (background for estimand design, §7.2)**

26. Genomic, immunologic, and prognostic associations of TROP2 (TACSTD2) expression in solid tumors. *The Oncologist* (2024). DOI: 10.1093/oncolo/oyae168
27. TROP2/claudin program mediates immune exclusion to impede checkpoint blockade in breast cancer. *Journal for ImmunoTherapy of Cancer* 14:e012265 (2026). DOI: 10.1136/jitc-2025-012265
28. Pan-cancer multi-omic integration of Trop2 reveals biological determinants and translational implications for ADC therapy. *npj Precision Oncology* (2026). DOI: 10.1038/s41698-026-01523-w
28a. Multicellular immune ecotypes within solid tumors predict real-world therapeutic benefits with immune checkpoint inhibitors. *Nature Communications* (2025). DOI: 10.1038/s41467-025-65016-3
28b. Bessede A, et al. TROP2 and atezolizumab primary resistance (OAK/POPLAR). *Clinical Cancer Research* (2024). PMID 38048058

**Statistical methods**

29. Aitchison J. *The Statistical Analysis of Compositional Data*. Chapman & Hall (1986).
30. Palarea-Albaladejo J, Martín-Fernández JA. zCompositions: R package for multivariate imputation of left-censored data under a compositional approach. *Chemometrics and Intelligent Laboratory Systems* 143:85–96 (2015).
31. Viechtbauer W. Conducting meta-analyses in R with the metafor package. *Journal of Statistical Software* 36:1–48 (2010).
32. Cribari-Neto F, Zeileis A. Beta regression in R. *Journal of Statistical Software* 34:1–24 (2010).
33. Yoshihara K, Shahmoradgoli M, Martínez E, et al. Inferring tumour purity and stromal and immune cell admixture from expression data (ESTIMATE). *Nature Communications* 4:2612 (2013).
34. Foroutan M, Bhuva DD, Lyu R, et al. Single sample scoring of molecular phenotypes (singscore). *BMC Bioinformatics* 19:404 (2018).
35. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics* 14:7 (2013).

---

## §13 BayesPrism / EcoTyper × TACSTD2 / CLDN4

中文 — `playbook.md` §7 规定**单个**上皮基因（`TACSTD2`）与 CD8/TLS 的关联，并把多种去卷积方法当作可比较的比例估计器。当问题收窄为 **BayesPrism 的区室表达 $Z$ 对 EcoTyper 的细胞状态/群落**，且暴露同时包括 **`TACSTD2` 与 `CLDN4`** 时，使用配套协议：

**[`bayesprism_ecotyper_vs_tacstd2_cldn4.md`](bayesprism_ecotyper_vs_tacstd2_cldn4.md)**

该文件规定：两个种子不可互换（`CLDN4` ≠ TROP2 结合伴侣；结合伴侣是 `CLDN7`）；E1–E3 对每个种子对称重跑；追加 E4（种子一致性）、E5（残差独有信息）、E6（锁定的 `BARRIER5` 区室分数）、E7（CE 丰度结局）、E8（上皮 CS 赋值，探索性）；以及四种预先定义的方法间格局（P-agree / P-gene-split / P-method-split / P-null）。主张必须同时带上 §8.4 的 S 级与该文件的 P 类。

EN — `playbook.md` §7 specifies the association of a **single** epithelial gene (`TACSTD2`) with CD8/TLS and treats deconvolution methods as comparable fraction estimators. When the question narrows to **BayesPrism compartment expression $Z$ versus EcoTyper cell states/communities**, with **both `TACSTD2` and `CLDN4`** as exposures, use the companion protocol:

**[`bayesprism_ecotyper_vs_tacstd2_cldn4.md`](bayesprism_ecotyper_vs_tacstd2_cldn4.md)**

That file specifies: the two seeds are not interchangeable (`CLDN4` is not the TROP2 binding partner; that partner is `CLDN7`); E1–E3 are rerun symmetrically per seed; E4 (seed concordance), E5 (residual unique information), E6 (locked `BARRIER5` compartment score), E7 (CE-abundance outcome), and E8 (epithelial CS assignment, exploratory) are added; and four pre-defined between-method patterns (P-agree / P-gene-split / P-method-split / P-null) license what may be claimed. Every claim must carry both an S grade from §8.4 and a P class from the companion.

---

*END OF PLAYBOOK — 本文件到此结束。任何数值结果请写入 `results/` 并引用本文件章节编号。/ Write all numeric results into `results/`, citing the section numbers above.*
