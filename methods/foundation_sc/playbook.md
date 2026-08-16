# Foundation & Deep Models for Single-Cell Lung Atlas Integration — A Methods Playbook
# 单细胞肺图谱整合的基础/深度模型方法手册

> **Scope / 范围**: Methods only. This document covers *how* and *when* to use scVI/scANVI, scGPT,
> Geneformer, and UCE for lung atlas integration and query projection, as of 2025–2026.
> It contains **no benchmark numbers** — where a published finding is referenced, it is named so you
> can read the primary source and run your own evaluation. Do not cite this file as a source of metrics.
>
> **范围**：仅方法。本文说明在 2025–2026 年阶段，如何以及何时使用 scVI/scANVI、scGPT、Geneformer、UCE
> 进行肺图谱整合与查询投影。**不含任何跑分数字**——凡引用已发表结论处均标明出处，请自行查阅原文并在自己的数据上评测。
> 请勿把本文件当作性能指标来源引用。

---

## 0. TL;DR / 速读

**EN.** For production lung-atlas *integration + reference mapping* today, the conditional VAE family
(scVI → scANVI, mapped via scArches) remains the default workhorse: it is purpose-built for batch
correction and semi-supervised label transfer, it exposes calibrated-ish uncertainty, and the Human
Lung Cell Atlas (HLCA) reference is distributed in exactly this format. Transformer/embedding
foundation models (scGPT, Geneformer, UCE) are most valuable as *complements*: cross-dataset /
cross-species embedding, zero-shot exploration of unseen data, gene-network and in-silico perturbation
analysis. All four models will happily assign a label to every query cell — so the dangerous failure
mode for a marker-defined population like **TACSTD2-high epithelium** is confident mislabeling, not a
crash. Treat every projected label as a hypothesis to be validated against marker genes and neighborhood
composition.

**中文.** 就当前生产级肺图谱**整合 + 参考映射**而言，条件变分自编码器家族（scVI → scANVI，经 scArches 映射）
仍是默认主力：它天生为批次校正与半监督标签迁移设计，能给出（大致）可用的不确定性，且人类肺细胞图谱（HLCA）
参考正是以该格式发布。Transformer/嵌入类基础模型（scGPT、Geneformer、UCE）最大价值在于**互补**：跨数据集/
跨物种嵌入、对未见数据的零样本探索、基因网络与虚拟扰动分析。四类模型都会给每个查询细胞强行分配标签——因此对
**TACSTD2 高表达上皮**这类以 marker 定义的群体，危险失败模式是"自信地贴错标签"，而非报错。请把每个投影标签
都当作待验证的假设，用 marker 基因和邻域组成去核对。

---

## 1. The models at a glance / 模型概览

### 1.1 scVI / scANVI (scvi-tools)

**EN.**
- **What it is.** scVI is a conditional variational autoencoder that models raw UMI counts with a
  (zero-inflated) negative binomial likelihood and learns a low-dimensional latent space while
  conditioning out batch/technical covariates. scANVI extends it semi-supervised: known cell-type labels
  shape the latent space and enable probabilistic label transfer to unlabeled cells.
- **Reference mapping.** scArches performs "architectural surgery": freeze the trained reference model
  and learn a small set of query-specific weights. This is the mechanism used to map new data onto the
  HLCA. It is lightweight and does not require re-sharing the reference training data.
- **Inputs.** Raw counts, a chosen batch key, and (for scANVI) a labels key. Highly-variable-gene (HVG)
  selection and consistent gene ordering matter a lot.
- **Where it shines for lung.** The HLCA reference is published as an scANVI model; batch integration
  across the many donors/protocols in lung studies is exactly its design target; it is fast enough to
  retrain and the projection step is cheap.

**中文.**
- **是什么。** scVI 是条件变分自编码器，用（零膨胀）负二项分布对原始 UMI 计数建模，学习低维潜空间并在其中
  剔除批次/技术协变量。scANVI 是其半监督扩展：已知细胞类型标签会塑造潜空间，并支持向未标注细胞做概率式标签迁移。
- **参考映射。** scArches 做"架构手术"：冻结已训练的参考模型，只学习一小组查询特异权重。这正是把新数据映射到
  HLCA 的机制。它轻量，且无需重新分发参考训练数据。
- **输入。** 原始计数、指定的 batch key，以及（scANVI 需要的）labels key。高变基因（HVG）选择与基因顺序一致性
  极其重要。
- **肺场景优势。** HLCA 参考即以 scANVI 模型发布；跨众多供者/建库方案的批次整合正是它的设计目标；重训练足够快，
  投影步骤开销很低。

### 1.2 scGPT

**EN.**
- **What it is.** A generative transformer pretrained on tens of millions of human cells (CELLxGENE).
  Genes are tokens; expression is value-binned. It supports zero-shot cell embeddings and fine-tuning
  for integration, cell-type annotation, gene-regulatory-network inference, and perturbation prediction.
- **Practical stance (2025–2026).** Fine-tuned scGPT is a credible integration/annotation option;
  zero-shot embeddings are more variable and have been reported to sometimes trail simple baselines
  (see §6). Treat zero-shot as exploratory, fine-tune when you need reliable integration.
- **Inputs.** Its own gene vocabulary and tokenizer; normalized/binned expression; a pretrained
  checkpoint. Genes outside the vocabulary are dropped.

**中文.**
- **是什么。** 在数千万人类细胞（CELLxGENE）上预训练的生成式 transformer。基因作为 token，表达量做分箱。支持零样本
  细胞嵌入，也可微调用于整合、细胞类型注释、基因调控网络推断和扰动预测。
- **实践取向（2025–2026）。** 微调后的 scGPT 是可信的整合/注释方案；零样本嵌入波动更大，且有报道其有时不及简单基线
  （见 §6）。把零样本当作探索手段，需要可靠整合时请微调。
- **输入。** 其自带基因词表与 tokenizer；归一化/分箱后的表达；预训练权重。词表外基因会被丢弃。

### 1.3 Geneformer

**EN.**
- **What it is.** A transformer pretrained on a very large human corpus (Genecorpus-30M, expanded to
  ~95M cells in the 2024 refresh, with larger multi-layer variants). Its key design choice is
  **rank-value encoding**: each cell is represented by genes ordered by expression normalized against a
  corpus-wide median, rather than by explicit expression values. Strong for network biology and
  in-silico perturbation; annotation is done by fine-tuning a classifier head.
- **Inputs.** Requires tokenization with Geneformer's gene-median dictionary and Ensembl gene IDs; the
  tokenizer step is CPU-bound preprocessing.

**中文.**
- **是什么。** 在超大人类语料上预训练的 transformer（Genecorpus-30M；2024 年更新扩至约 9500 万细胞，并提供更大的
  多层变体）。核心设计是**秩-值编码**：每个细胞用"按（相对全语料中位数归一化后的）表达排序的基因序列"表示，而非显式表达值。
  在网络生物学与虚拟扰动上很强；注释需微调一个分类头。
- **输入。** 需用 Geneformer 的基因中位数字典和 Ensembl 基因 ID 做 tokenization；该步骤是 CPU 密集的预处理。

### 1.4 UCE (Universal Cell Embeddings)

**EN.**
- **What it is.** A large transformer that produces **zero-shot** embeddings across tissues and species.
  Its trick is representing each gene by a protein language-model embedding (ESM2) of its product, so it
  can embed genes/species never seen in training without retraining. Trained on tens of millions of cells
  across species. Typically used embedding-only (no fine-tuning), which makes it attractive for novel or
  cross-species query data.
- **Inputs.** Counts + species; it maps genes via protein embeddings internally. It is the heaviest of
  the four to run.

**中文.**
- **是什么。** 一个大型 transformer，产出跨组织、跨物种的**零样本**嵌入。其巧思是用基因产物的蛋白语言模型嵌入
  （ESM2）来表示基因，因而无需重训练即可嵌入训练中从未见过的基因/物种。在跨物种的数千万细胞上训练。通常仅用于嵌入
  （不微调），这使它在新颖或跨物种查询数据上很有吸引力。
- **输入。** 计数 + 物种；内部通过蛋白嵌入映射基因。它是四者中运行开销最大的。

---

## 2. When they *help* lung atlas integration / 何时对肺图谱整合有帮助

**EN.** Pick the model by the job, not by novelty:

| Task / 任务 | First choice | Why |
|---|---|---|
| Integrate many donors/protocols into one lung reference | **scVI/scANVI** | Built for batch correction + label-aware latent space; HLCA precedent |
| Map a new query onto HLCA | **scANVI + scArches** | Native, lightweight surgery; per-cell transfer + uncertainty |
| Explore a brand-new dataset with no reference yet | **UCE** (zero-shot), scGPT (zero-shot) | Embedding without training; quick lay-of-the-land |
| Cross-species lung comparison (e.g. human ↔ mouse) | **UCE** | Protein-embedding genes generalize across species |
| Gene-network / in-silico perturbation hypotheses | **Geneformer**, scGPT | Designed for network + perturbation tasks |
| High-accuracy annotation with labeled reference | **scANVI** or **fine-tuned** scGPT/Geneformer | Supervision beats zero-shot for label fidelity |

**EN — concrete "help" scenarios.**
1. **Donor/batch heterogeneity dominates.** Lung atlases pool many chemistries, sites, disease states.
   Conditional VAEs remove this structured nuisance while preserving biological variation — the core win.
2. **You already have a curated reference (HLCA).** scArches lets you project without retraining and
   returns a transfer probability per cell — use that probability.
3. **No reference / exotic data.** Foundation-model embeddings give an immediate coordinate system for
   clustering and visualization before you commit to a reference.
4. **Cross-species or unseen genes.** UCE's protein-embedding gene representation is the one design that
   natively tolerates gene sets it never trained on.

**中文.** 按任务而非新鲜度选择模型：

| 任务 | 首选 | 原因 |
|---|---|---|
| 把众多供者/方案整合成单一肺参考 | **scVI/scANVI** | 为批次校正 + 标签感知潜空间而生；有 HLCA 先例 |
| 把新查询映射到 HLCA | **scANVI + scArches** | 原生、轻量手术；逐细胞迁移 + 不确定性 |
| 尚无参考时探索全新数据集 | **UCE**（零样本）、scGPT（零样本） | 无需训练即得嵌入；快速摸底 |
| 跨物种肺对比（如 人 ↔ 鼠） | **UCE** | 蛋白嵌入式基因表示可跨物种泛化 |
| 基因网络/虚拟扰动假设 | **Geneformer**、scGPT | 面向网络与扰动任务设计 |
| 有标注参考、要高精度注释 | **scANVI** 或**微调后**的 scGPT/Geneformer | 标签保真度上，监督优于零样本 |

**中文 — 具体"有帮助"的场景。**
1. **供者/批次异质性主导。** 肺图谱汇集多种化学、多中心、多疾病状态。条件 VAE 在保留生物学变异的同时剔除这种结构化
   噪声——这是核心收益。
2. **已有精编参考（HLCA）。** scArches 让你无需重训练即可投影，并返回逐细胞迁移概率——请使用该概率。
3. **无参考/异常数据。** 基础模型嵌入在你确定参考之前，先给出可聚类、可可视化的坐标系。
4. **跨物种或未见基因。** UCE 的蛋白嵌入式基因表示，是唯一天生能容忍"从未训练过的基因集"的设计。

---

## 3. When they *hallucinate* labels / 何时会"臆造"标签

**EN.** "Hallucination" here means a model producing a confident, plausible-looking label or embedding
neighborhood that is not biologically supported. Root causes and where each model is exposed:

1. **Every cell gets a label (closed-world assumption).** scANVI, and any fine-tuned classifier head on
   scGPT/Geneformer, assigns the *nearest reference label* even to states absent from the reference —
   novel disease epithelium, doublets, ambient-RNA artifacts, low-quality cells. There is no built-in
   "none of the above."
2. **Rare / aberrant populations collapse into common neighbors.** A small marker-defined population
   (e.g. TACSTD2-high aberrant basaloid) can be absorbed into the nearest abundant type (basal, club,
   AT2) in the latent space, especially with aggressive HVG filtering or heavy batch correction.
3. **Over-correction erases real biology.** Strong batch conditioning can merge genuinely distinct states
   if disease/condition is confounded with batch — the "integration removed my signal" failure.
4. **Zero-shot embeddings can encode nuisance, not biology.** Reported findings (see §6) show zero-shot
   foundation-model embeddings sometimes cluster by depth/protocol or underperform simple HVG+PCA on
   some tasks. A clean-looking UMAP is not evidence of correct structure.
5. **Vocabulary / gene-ID drift.** Genes outside a model's vocabulary (scGPT) or unmapped Ensembl IDs
   (Geneformer) are silently dropped; if your discriminating markers are dropped, the model literally
   cannot see the population that defines your query.
6. **Label transfer inherits reference errors.** If the reference mislabels a state, projection
   faithfully reproduces the mistake with high confidence.

**EN — detection & mitigation.**
- **Use and threshold uncertainty.** scANVI/scArches give per-cell transfer probability / entropy —
  flag low-confidence cells as "unassigned" instead of forcing a label.
- **Marker-anchor everything.** Independently verify each projected label against canonical markers
  (see §5) *in the query's own expression*, not just in the embedding.
- **Inspect neighborhood composition.** For each query cell, check the reference cell types among its
  nearest neighbors; heterogeneous neighborhoods signal an out-of-distribution or ambiguous cell.
- **Run an OOD / novelty check.** Reconstruction error (VAE), distance-to-reference in embedding space,
  or a dedicated novelty detector; treat outliers as "candidate novel," not misclassified.
- **Keep discriminating genes in-vocabulary.** Confirm your markers survive tokenization before trusting
  any transformer output.
- **Never present a single model's labels as ground truth.** Agreement across an independent method
  (marker scoring, or a second model) is the minimum bar.

**中文.** 此处"臆造"指模型给出自信、看似合理、但无生物学支持的标签或嵌入邻域。成因及各模型的暴露点：

1. **每个细胞都被贴标签（封闭世界假设）。** scANVI，以及 scGPT/Geneformer 上任何微调分类头，都会把参考中不存在的
   状态（新型疾病上皮、双胞体、环境 RNA 伪迹、低质量细胞）也贴上*最近的参考标签*。没有内建的"以上都不是"。
2. **稀有/异常群体塌缩进常见邻居。** 以 marker 定义的小群体（如 TACSTD2 高的异常基底样细胞）可能在潜空间被吸收进
   最近的高丰度类型（基底、club、AT2），在激进 HVG 过滤或强批次校正下尤甚。
3. **过度校正抹除真实生物学。** 若疾病/状态与批次混淆，强批次条件化会把本应区分的状态合并——即"整合把我的信号
   洗掉了"。
4. **零样本嵌入可能编码噪声而非生物学。** 已有报道（见 §6）显示零样本基础模型嵌入有时按测序深度/方案聚类，或在某些
   任务上不及简单 HVG+PCA。UMAP 好看不等于结构正确。
5. **词表/基因 ID 漂移。** 模型词表外基因（scGPT）或未映射的 Ensembl ID（Geneformer）会被静默丢弃；若你的判别性
   marker 被丢弃，模型根本"看不到"定义你查询的那个群体。
6. **标签迁移继承参考错误。** 若参考本身标错某状态，投影会高置信度地忠实复现该错误。

**中文 — 检测与缓解。**
- **使用并设阈值于不确定性。** scANVI/scArches 提供逐细胞迁移概率/熵——把低置信细胞标为"未分配"，而非强行贴标签。
- **凡事以 marker 锚定。** 在*查询自身的表达*中（而非仅在嵌入里）用经典 marker（见 §5）独立核验每个投影标签。
- **检查邻域组成。** 对每个查询细胞，查看其最近邻中的参考细胞类型；邻域混杂说明该细胞离群或含糊。
- **做 OOD/新颖性检查。** 重构误差（VAE）、嵌入空间到参考的距离，或专用新颖性检测器；把离群者当作"疑似新颖"而非
  错分。
- **确保判别基因在词表内。** 信任任何 transformer 输出前，先确认你的 marker 在 tokenization 后存活。
- **切勿把单一模型的标签当真值。** 至少要与一个独立方法（marker 打分或第二个模型）取得一致。

---

## 4. Projecting TACSTD2-high epithelium / 投影 TACSTD2 高表达上皮

**EN — why this population is a stress test.** TACSTD2 (TROP2) marks airway/epithelial states —
enriched in basal, secretory/club, and hillock-like airway epithelium, and it separates airway from
alveolar (AT1/AT2) epithelium. In disease (e.g. fibrosis) TROP2 also appears in metaplastic and
aberrant basaloid (KRT5–/KRT17+) epithelium. These are often **rare, condition-specific, and close in
transcriptome to abundant neighbors** — exactly the setting where §3's hallucination modes bite. So
"project TACSTD2-high epithelium" is really "map a rare, marker-defined, possibly-novel epithelial state
without letting the model collapse it into basal/club/AT2."

**EN — a defensive projection recipe.**
1. **Define the query population explicitly first.** In the query data, identify TACSTD2-high epithelial
   cells by co-expression (TACSTD2 with epithelial identity) *before* projecting, so you have an
   independent ground-truth mask to check the projection against.
2. **Confirm markers survive preprocessing.** Ensure `TACSTD2`/`TROP2` and companion genes are in the
   HVG set (scVI/scANVI), in the model vocabulary (scGPT), and map to valid Ensembl IDs (Geneformer). If
   `TACSTD2` is filtered out, your projection is blind to the very axis you care about.
3. **Choose the mapping strategy.**
   - *Reference exists and contains the relevant airway epithelial states (HLCA):* use **scANVI +
     scArches**; keep the per-cell transfer probability.
   - *Reference may lack the state (disease/aberrant epithelium):* prefer an **embedding-first** approach
     (scVI latent, UCE, or scGPT zero-shot) and cluster de novo; do not force a reference label.
4. **Decide compartment scope.** Projecting within the epithelial compartment sharpens fine airway
   distinctions; projecting against the full atlas guards against mislabeling a non-epithelial contaminant
   as epithelium. Doing both and comparing is the safe move.
5. **Validate the projection against the mask.** Overlay your independent TACSTD2-high mask on the
   embedding/labels. Expected: the mask forms a coherent neighborhood, not scatter across basal/club/AT2.
   Scatter = collapse/hallucination.
6. **Check neighborhood composition + uncertainty for the masked cells.** High entropy or mixed
   neighborhoods on TACSTD2-high cells is the signal that the reference does not truly contain this state
   — report it as candidate-novel rather than assigning basal/club/AT2.
7. **Corroborate with companion markers, not TACSTD2 alone.** Airway vs alveolar and normal vs aberrant
   need a panel (see §5), because TROP2 alone is expressed across several epithelial states.

**中文 — 为何该群体是压力测试。** TACSTD2（TROP2）标记气道/上皮状态——在基底、分泌/club、以及气道 hillock 样上皮中
富集，并把气道上皮与肺泡（AT1/AT2）上皮区分开。在疾病（如纤维化）中，TROP2 也出现于化生上皮及异常基底样
（KRT5–/KRT17+）上皮。这些群体常常**稀有、依赖状态、且转录组上贴近高丰度邻居**——正是 §3 各类臆造模式最易发作之处。
因此"投影 TACSTD2 高表达上皮"实质是"映射一个稀有、marker 定义、可能新颖的上皮状态，同时不让模型把它塌缩进
基底/club/AT2"。

**中文 — 防御式投影配方。**
1. **先明确定义查询群体。** 在查询数据里，先通过共表达（TACSTD2 + 上皮身份）识别 TACSTD2 高上皮细胞*再*投影，从而拿到
   一个独立的真值掩码来核对投影。
2. **确认 marker 在预处理中存活。** 保证 `TACSTD2`/`TROP2` 及配套基因进入 HVG 集（scVI/scANVI）、在模型词表内
   （scGPT）、并映射到有效 Ensembl ID（Geneformer）。若 `TACSTD2` 被过滤掉，你的投影对最关心的那条轴是"瞎的"。
3. **选择映射策略。**
   - *参考存在且含相关气道上皮状态（HLCA）：* 用 **scANVI + scArches**；保留逐细胞迁移概率。
   - *参考可能缺该状态（疾病/异常上皮）：* 优先**嵌入优先**（scVI 潜空间、UCE 或 scGPT 零样本）并从头聚类；不要强套参考
     标签。
4. **确定分区范围。** 在上皮分区内投影能锐化气道内部的精细区分；对全图谱投影则可防止把非上皮污染物误标为上皮。两者都做
   并对比是稳妥之举。
5. **用掩码验证投影。** 把独立的 TACSTD2 高掩码叠加到嵌入/标签上。预期：掩码形成连贯邻域，而非散布于基底/club/AT2。
   散布 = 塌缩/臆造。
6. **对掩码细胞检查邻域组成 + 不确定性。** TACSTD2 高细胞上的高熵或混杂邻域，正说明参考中并不真正含该状态——应报告为
   疑似新颖，而非贴上基底/club/AT2。
7. **用配套 marker 佐证，勿只靠 TACSTD2。** 气道 vs 肺泡、正常 vs 异常都需要一组 marker（见 §5），因为 TROP2 单独一个
   在多个上皮状态中都有表达。

---

## 5. Marker panel for validating epithelial projections / 上皮投影验证 marker 面板

**EN.** Use these as an *independent* check on any projected label. This is a starting panel, not an
exhaustive taxonomy — confirm against your own reference's annotations.

| Compartment / state | Illustrative markers |
|---|---|
| Airway epithelial identity (incl. TROP2+) | `TACSTD2` (TROP2), `EPCAM`, `KRT8`, `KRT18` |
| Basal | `KRT5`, `KRT17`, `TP63` |
| Club / secretory | `SCGB1A1`, `SCGB3A2` |
| Goblet / mucous | `MUC5AC`, `MUC5B` |
| Ciliated | `FOXJ1`, `TUBB4B` |
| Alveolar type 2 (AT2) | `SFTPC`, `SFTPB`, `LAMP3` |
| Alveolar type 1 (AT1) | `AGER`, `PDPN` |
| Aberrant basaloid (disease) | `KRT17`, `KRT5`(low), `CDKN1A`, `KRT8` |

**中文.** 把以上作为对任何投影标签的*独立*核查。这是起步面板，非穷尽分类——请与你自己参考的注释对齐确认。

| 分区/状态 | 示例 marker |
|---|---|
| 气道上皮身份（含 TROP2+） | `TACSTD2`（TROP2）、`EPCAM`、`KRT8`、`KRT18` |
| 基底 | `KRT5`、`KRT17`、`TP63` |
| Club/分泌 | `SCGB1A1`、`SCGB3A2` |
| 杯状/黏液 | `MUC5AC`、`MUC5B` |
| 纤毛 | `FOXJ1`、`TUBB4B` |
| 肺泡 II 型（AT2） | `SFTPC`、`SFTPB`、`LAMP3` |
| 肺泡 I 型（AT1） | `AGER`、`PDPN` |
| 异常基底样（疾病） | `KRT17`、`KRT5`(低)、`CDKN1A`、`KRT8` |

---

## 6. Reading the evidence honestly / 客观看待证据

**EN.** This playbook deliberately gives you no numbers. When you need to decide, read the primary
sources and — critically — **run the comparison on your own lung data**, because performance is
task-, dataset-, and preprocessing-dependent. Landmarks worth reading:
- **HLCA** (Sikkema et al., *Nat Med* 2023) — the reference and the scANVI/scArches mapping workflow.
- **scVI / scANVI / scArches** (Lopez 2018; Xu 2021; Lotfollahi 2022) — the model and mapping mechanics.
- **scGPT** (Cui et al., *Nat Methods* 2024) — architecture and fine-tuning tasks.
- **Geneformer** (Theodoris et al., *Nature* 2023; 2024 corpus/model refresh) — rank-value encoding,
  perturbation.
- **UCE** (Rosen et al., 2023–2024) — cross-species zero-shot embeddings via protein embeddings.
- **Zero-shot limits** (Kedzierska et al., 2023, "Assessing the limits of zero-shot foundation models in
  single-cell biology") — why zero-shot embeddings can underperform simple baselines; motivates §3.
- **Integration benchmarking** (Luecke et al., scIB, 2022) — a framework/metrics to run your *own*
  batch-correction vs bio-conservation evaluation.

**EN — how to benchmark without faking it.** Hold out donors/batches (not random cells). Report both
batch-mixing *and* biology-conservation — a method can win one by sacrificing the other. Include a simple
baseline (HVG + PCA + Harmony/BBKNN). Validate rare populations (like TACSTD2-high) separately from bulk
accuracy, since global metrics hide rare-cell collapse. Never quote a metric you did not compute on the
data in question.

**中文.** 本手册刻意不给数字。需要决策时，请读原文，并且——最关键——**在你自己的肺数据上跑对比**，因为性能依赖于任务、
数据集与预处理。值得读的里程碑：
- **HLCA**（Sikkema 等，*Nat Med* 2023）——参考本体与 scANVI/scArches 映射流程。
- **scVI / scANVI / scArches**（Lopez 2018；Xu 2021；Lotfollahi 2022）——模型与映射机制。
- **scGPT**（Cui 等，*Nat Methods* 2024）——架构与微调任务。
- **Geneformer**（Theodoris 等，*Nature* 2023；2024 语料/模型更新）——秩-值编码、扰动。
- **UCE**（Rosen 等，2023–2024）——经蛋白嵌入的跨物种零样本嵌入。
- **零样本局限**（Kedzierska 等，2023，"评估单细胞生物学中零样本基础模型的极限"）——为何零样本嵌入会不及简单基线；
  §3 的动因。
- **整合基准**（Luecke 等，scIB，2022）——一套可用来跑*你自己*的批次校正 vs 生物学保留评测的框架/指标。

**中文 — 如何不作假地做基准。** 按供者/批次留出（而非随机留细胞）。同时报告批次混合*与*生物学保留——一个方法可能靠牺牲另一个
来赢。纳入简单基线（HVG + PCA + Harmony/BBKNN）。对稀有群体（如 TACSTD2 高）单独验证，别只看整体准确率，因为全局指标会掩盖
稀有细胞的塌缩。绝不引用你未在相关数据上亲自计算过的指标。

---

## 7. Practical GPU vs CPU notes / GPU 与 CPU 实务

**EN — general.**
- **VAEs (scVI/scANVI):** GPU gives a large training speedup, but they run acceptably on CPU for modest
  datasets, and the **scArches query mapping step is cheap** (small adapter weights) — often fine on CPU.
  Memory footprint is modest. Apple Silicon (MPS) support is partial; CUDA is the smooth path.
- **Transformers (scGPT, Geneformer, UCE):** GPU is strongly recommended even for inference. CPU-only
  zero-shot embedding is possible for small cell numbers but slow; batch/atlas-scale on CPU is impractical.
  Use mixed precision (fp16/bf16), tune batch size to VRAM, and enable efficient attention kernels if your
  stack supports them.
- **UCE is the heaviest** of the four to run (large model + protein-embedding gene mapping); budget the
  most VRAM/time for it.
- **Geneformer tokenization is CPU-bound preprocessing** (building rank-value tokens with the gene-median
  dictionary); parallelize it on CPU and cache the tokenized dataset to disk so you tokenize once.

**EN — sizing heuristics (order-of-magnitude, verify on your hardware).**
- Inference/zero-shot embedding of the transformer models is comfortable on a single modern data-center or
  high-end consumer GPU with ~16 GB+ VRAM; reduce batch size if you hit OOM.
- Fine-tuning transformers needs more memory than inference; use gradient accumulation / smaller batches /
  gradient checkpointing to fit.
- For scVI/scANVI, wall-clock scales with cell count and epochs; the projection (scArches) is far cheaper
  than the initial reference training.
- **Reproducibility:** pin library and checkpoint versions, set seeds, and record the exact gene
  vocabulary / HVG list — results shift across model releases (e.g. Geneformer corpus updates).

**中文 — 通则.**
- **VAE（scVI/scANVI）：** GPU 带来显著训练加速，但在中等规模数据上 CPU 也可接受，且 **scArches 查询映射步骤开销很小**
  （只学小适配器权重）——常常在 CPU 上就行。显存占用不高。Apple Silicon（MPS）支持不完整；CUDA 才是顺畅路径。
- **Transformer（scGPT、Geneformer、UCE）：** 即便只做推理也强烈建议用 GPU。仅 CPU 的零样本嵌入在小细胞数下可行但很慢；
  图谱级规模用 CPU 不现实。使用混合精度（fp16/bf16），按显存调 batch size，若技术栈支持则启用高效注意力核。
- **UCE 是四者中运行开销最大的**（大模型 + 蛋白嵌入基因映射）；为其预留最多显存/时间。
- **Geneformer 的 tokenization 是 CPU 密集预处理**（用基因中位数字典构建秩-值 token）；在 CPU 上并行并把 tokenized 数据
  缓存到磁盘，做到只 tokenize 一次。

**中文 — 规模经验（数量级，请在自己硬件上核实）。**
- Transformer 模型的推理/零样本嵌入，在单张约 16 GB+ 显存的现代数据中心或高端消费级 GPU 上通常从容；遇到 OOM 就调小
  batch size。
- 微调 transformer 比推理更吃显存；用梯度累积/更小 batch/梯度检查点来塞下。
- scVI/scANVI 的墙钟时间随细胞数与 epoch 增长；投影（scArches）远比初始参考训练便宜。
- **可复现性：** 固定库与权重版本、设随机种子、记录确切的基因词表/HVG 列表——不同模型版本间结果会变（如 Geneformer 语料
  更新）。

---

## 8. A defensible default workflow / 一套站得住脚的默认流程

**EN.**
1. QC and preprocess consistently; keep raw counts; confirm discriminating markers (incl. `TACSTD2`)
   survive filtering/vocabulary.
2. Integrate/map with **scANVI + scArches** onto HLCA; retain per-cell transfer probability.
3. In parallel, compute an independent embedding (scVI latent and/or one foundation model) for
   cross-checking structure.
4. Threshold uncertainty → mark low-confidence cells "unassigned."
5. Validate every label of interest against the §5 marker panel in the query's own expression.
6. For rare/marker-defined targets (TACSTD2-high), follow §4: define the mask first, check neighborhood
   composition, and report OOD candidates instead of forcing a reference label.
7. If you need numbers, benchmark on held-out donors with scIB metrics and a simple baseline (§6). Report
   rare-population accuracy separately.

**中文.**
1. 一致地做 QC 与预处理；保留原始计数；确认判别性 marker（含 `TACSTD2`）在过滤/词表后存活。
2. 用 **scANVI + scArches** 映射到 HLCA；保留逐细胞迁移概率。
3. 并行地算一个独立嵌入（scVI 潜空间和/或某个基础模型）以交叉核对结构。
4. 对不确定性设阈值 → 把低置信细胞标为"未分配"。
5. 在查询自身表达中用 §5 marker 面板验证每个关注标签。
6. 对稀有/marker 定义目标（TACSTD2 高），照 §4：先定义掩码、检查邻域组成、报告 OOD 候选而非强套参考标签。
7. 若需数字，用留出供者 + scIB 指标 + 简单基线做基准（§6）。单独报告稀有群体准确率。

---

*Methods only. No fabricated benchmarks. Verify every claim on your own data and cite primary sources.*
*仅方法。无编造跑分。请在自己数据上验证每条主张并引用原始文献。*
