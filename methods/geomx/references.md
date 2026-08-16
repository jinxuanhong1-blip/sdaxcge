# Annotated references / 参考文献（含注解）

Every peer-reviewed entry below was resolved against PubMed (PMID and DOI given) at the time of writing.
Entries without a PMID are software documentation, vendor manuals, preprints or conference abstracts and are
labelled as such — cite them accordingly.

**中文:** 以下所有同行评议文献均已在 PubMed 核对（给出 PMID 与 DOI）。没有 PMID 的条目为软件文档、厂商手册、预印本或
会议摘要，已明确标注，引用时请注意其证据等级。

Recency note: the brief asked for 2024–2026 DSP statistics work. Section A is the 2024–2026 core; Section B
holds the foundational (pre-2024) methods that the 2024–2026 papers build on and that you will still cite.

---

## A. 2024–2026 core / 2024–2026 核心文献

**A1 · standR — the reference GeoMx DSP workflow**
Liu N, Bhuva DD, Mohamed A, Bokelund M, Kulasinghe A, Tan CW, Davis MJ.
*standR: spatial transcriptomic analysis for GeoMx DSP data.*
**Nucleic Acids Research** 2024;52(1):e2. PMID 37953397. doi:10.1093/nar/gkad1026
Surveys published DSP analyses (2020–2022) and shows the community default — ProbeQC + LOQ filtering, Q3
normalization, t-test or LMM — underperforms. Proposes the SpatialExperiment-based workflow used throughout
this playbook: QC → TMM/other normalization → RUV4 batch correction → **limma-voom with
`duplicateCorrelation` blocking on patient**. This is the single most useful citation for justifying your
pipeline to a reviewer.
**中文:** 系统梳理 2020–2022 年 DSP 文献，指出 "ProbeQC+LOQ 过滤 → Q3 → t 检验/LMM" 的社区默认流程表现不佳，
提出基于 SpatialExperiment 的推荐流程（本手册主线）：质控 → TMM 等归一化 → RUV4 批次校正 →
**limma-voom + 以患者为 block 的 `duplicateCorrelation`**。是向审稿人论证流程合理性的首选引用。

**A2 · Library size is confounded with biology in spatial data**
Bhuva DD, Tan CW, Salim A, Marceaux C, Pickering MA, Chen J, et al.
*Library size confounds biology in spatial transcriptomics data.*
**Genome Biology** 2024;25(1):99. PMID 38637899. doi:10.1186/s13059-024-03241-7
Demonstrates that in spatial assays library size carries biological signal (cell number, cell type, tissue
density), so scaling normalization removes signal along with noise. For GeoMx this is directly the
tumor-AOI-vs-immune-AOI size problem. The reason §4.2 of the playbook tells you to decide *deliberately*
whether to normalize compartments jointly or separately.
**中文:** 证明空间数据中文库大小本身携带生物学信息（细胞数、细胞类型、组织密度），缩放归一化会连同信号一起消除。
对 GeoMx 而言即 "肿瘤 AOI 大、免疫 AOI 小" 的问题，是手册 §4.2 要求明确选择合并/分别归一化的依据。

**A3 · smiDE — contamination-aware, spatially correlated DE**
Vasconcelos AG, McGuire D, Simon N, et al.
*Differential expression analysis for spatially correlated data using smiDE.*
**Genome Biology** 2026;27(1):21. PMID 41519776. doi:10.1186/s13059-025-03867-1
Developed for CosMx single-cell spatial data, but two components transfer directly to GeoMx segment analysis:
(i) explicit modelling and filtering of genes prone to **segmentation/contamination bias** — the formal
version of the spillover controls in playbook §6.3; (ii) benchmarking of strategies for spatial correlation,
showing that ignoring it inflates false discoveries.
**中文:** 虽为 CosMx 单细胞空间数据开发，但两点可直接迁移到 GeoMx 分区分析：①对**分割/污染偏倚**易感基因的显式
建模与过滤（即手册 §6.3 渗漏对照的正式版本）；②系统比较空间相关性的处理策略，证明忽略它会抬高假阳性。

**A4 · The GSE271689-class study itself**
Aung TN, Monkman J, Warrell J, Vathiotis I, Bates KM, Gavrielatou N, et al., Kulasinghe A, Rimm DL.
*Spatial signatures for predicting immunotherapy outcomes using multi-omics in non-small cell lung cancer.*
**Nature Genetics** 2025;57(10):2482–2493. PMID 41073787. doi:10.1038/s41588-025-02351-7
The published analysis of the Yale/UQ/Athens NSCLC cohorts (GEO: GSE271689, GSE292098). Design template for
this playbook: FFPE TMA, **4 ROIs per tumor**, morphology markers **PanCK (tumor) / CD45 (leukocyte) /
CD68 (macrophage)** plus SYTO 13, WTA readout, 131 patients in the transcriptomic arm. Methodologically note:
quantile-normalized matrices, CIBERSORTx/LM22 deconvolution of the stromal compartment, LASSO-derived
signatures, **tertile cut point fixed in training and applied unchanged to validation**, two-sided tests in
discovery and one-sided in validation. Worth reading specifically for how they handled the
train/validate split.
**中文:** Yale/UQ/雅典三队列 NSCLC 研究的正式发表版本（GEO: GSE271689、GSE292098），也是本手册的设计模板：FFPE 组织
芯片、**每瘤 4 个 ROI**、形态学标记 **PanCK/CD45/CD68** + SYTO 13、WTA、转录组部分 131 例。方法学要点：分位数归一化、
对基质分区做 CIBERSORTx/LM22 解卷积、LASSO 构建 signature、**三分位切点在训练集固定后原封不动用于验证集**、发现集双侧
检验而验证集单侧检验。其训练/验证划分方式尤其值得参考。

**A5 · Applied template: LMM for markers, mixed-effects Cox for survival**
Mateiou C, Lokhande L, Diep LH, Knulst M, et al.
*Spatial tumor immune microenvironment phenotypes in ovarian cancer.*
**npj Precision Oncology** 2024;8(1):148. PMID 39026018. doi:10.1038/s41698-024-00640-8
A clean, citable precedent for the exact statistical stack in this playbook: immune and tumor data
**normalized separately** to avoid AOI-size scaling bias, `lmerTest` LMMs with **patient ID as random effect**
for group comparisons, and **`coxme` mixed-effects Cox models** for survival with repeated ROIs.
**中文:** 与本手册统计栈高度一致的可引用先例：免疫与肿瘤数据**分别归一化**以避免 AOI 面积造成的缩放偏倚；用
`lmerTest` 以**患者 ID 为随机效应**做组间比较；用 **`coxme` 混合效应 Cox 模型**处理多 ROI 的生存分析。

**A6 · Applied template: standR + voomLmFit in a TMA cohort**
Han YB, Lee S, Lee JO, Jeong SI, et al.
*Spatial transcriptomics reveal high T cell and monocyte status as predictive and prognostic markers in
pancreatic cancer.*
**Journal of Translational Medicine** 2025;23(1):576. PMID 40410886. doi:10.1186/s12967-025-06599-9
PanCK/CD45-segmented TMA analysed with standR: TMM normalization, `edgeR::filterByExpr` gene filtering,
**`edgeR::voomLmFit` with duplicate correlation**, BH control, then clusterProfiler GSEA. A concise methods
paragraph to model your own on.
**中文:** PanCK/CD45 分区的 TMA，使用 standR 流程：TMM 归一化、`edgeR::filterByExpr` 过滤、
**`edgeR::voomLmFit` + duplicate correlation**、BH 校正、clusterProfiler GSEA。其方法学段落可作为写作范本。

**A7 · Applied template: standR QC thresholds spelled out**
Sadeghirad H, Monkman J, Tan CW, et al.
*Spatial dynamics of tertiary lymphoid aggregates in head and neck cancer: insights into immunotherapy
response.*
**Journal of Translational Medicine** 2024;22(1):677. PMID 39049036. doi:10.1186/s12967-024-05409-y
Useful because it states its filters numerically (ROI detection count < 350,000; nuclei < 250; > 3% low-expressing
genes; RLE/PCA outlier removal) and uses `voomLmFit` with sample weights.
**中文:** 价值在于把过滤阈值写成了具体数字（ROI 检出计数 < 350,000、细胞核 < 250、低表达基因 > 3%、RLE/PCA 剔除离群），
并使用带样本权重的 `voomLmFit`。

**A8 · TROP2/TACSTD2 is epithelium-restricted — the biological prior**
Notini G, Galbardi B, Viale G, et al.
*Pan-cancer multi-omic integration of Trop2 reveals biological determinants and translational implications
for ADC therapy.*
**npj Precision Oncology** 2026;10(1). PMID 42297906. doi:10.1038/s41698-026-01523-w
Across single-cell datasets, > 70–95% of malignant epithelial cells express `TACSTD2` versus < 2% of stromal,
immune and endothelial cells. This is the prior against which your CD45/CD68 segment signal should be judged:
appreciable immune-AOI `TACSTD2` is more likely spillover than biology. Also reviews why **membrane-localized**
TROP2 (not total mRNA/protein) tracked response in TROPION-Lung01 exploratory analyses — a limit on what any
transcriptomic assay can claim about ADC benefit.
**中文:** 跨单细胞数据集，`TACSTD2` 在 >70–95% 的恶性上皮细胞中表达，而在基质/免疫/内皮细胞中 <2%。这是判断
CD45/CD68 分区信号的先验：免疫 AOI 中出现可观的 `TACSTD2` 更可能是渗漏而非生物学。该文亦综述了 TROPION-Lung01
探索性分析中**膜定位** TROP2（而非总 mRNA/蛋白）才与疗效相关——这划定了转录组检测在 ADC 获益预测上的能力边界。

**A9 · CLDN4 as a compartment-specific spatial target**
Wang J, Seo JW, Kare AJ, et al.
*Spatial transcriptomic analysis drives PET imaging of tight junction protein expression in pancreatic cancer
theranostics.*
**Nature Communications** 2024;15(1):10751. PMID 39737976. doi:10.1038/s41467-024-54761-6
Uses spatial transcriptomics to select an imaging target: cancer cell-surface markers are spatially
correlated with each other and give specific cancer localization, while their correlation with
immune/fibroblast markers is low; `CLDN4` is ~16-fold higher in cancer than normal pancreas. The
methodological point to borrow: **target selection is a spatial-specificity argument, not a fold-change
argument** — report co-localization with epithelial markers and anti-correlation with immune/stromal markers.
**中文:** 用空间转录组筛选显像靶点：癌细胞表面标记彼此空间相关且定位特异，与免疫/成纤维标记相关性低；`CLDN4` 在
癌组织较正常胰腺高约 16 倍。可借鉴的方法学要点是：**靶点选择是空间特异性论证，而非单纯的 fold change 论证**——应报告
与上皮标记的共定位和与免疫/基质标记的反相关。

**A10 · ADC-target spatial heterogeneity and sampling sufficiency**
Sabtan D, Eich ML, Loch F, et al.
*Spatial heterogeneity of antibody-drug conjugate targets in pancreatic ductal adenocarcinoma.*
**Journal of Pathology: Clinical Research** 2026;12(2):e70083. PMID 41837391. doi:10.1002/2056-4538.70083
62 patients × 6 cores = 1,116 cores scored for c-MET, NECTIN4 and TROP-2, with an explicit **sampling
sufficiency simulation** ("how many cores are needed to reach the maximum score?"). This is the design
question behind `heterogeneity_score` in playbook §7 and the ROI-count planning in §10 — and a template for
answering it with your own ROIs.
**中文:** 62 例 × 6 个 core = 1,116 个 core 评估 c-MET/NECTIN4/TROP-2，并做了明确的**取样充分性模拟**（"需要几个 core
才能达到最高评分"）。这正是手册 §7 中 `heterogeneity_score` 与 §10 中 ROI 数量规划背后的设计问题，可作为用自有 ROI
回答该问题的模板。

---

## B. Foundational methods still required / 仍需引用的基础方法学

**B1 · Merritt CR, Ong GT, Church SE, et al.** *Multiplex digital spatial profiling of proteins and RNA in
fixed tissue.* **Nature Biotechnology** 2020;38(5):586–599. PMID 32393914. doi:10.1038/s41587-020-0472-9
The platform paper; cite for the assay, UV-cleavable barcode chemistry and segmentation principle.
**中文:** 平台原始论文；引用于检测原理、紫外可切割条码化学与分割机制。

**B2 · Danaher P, Kim Y, Nelson B, et al.** *Advances in mixed cell deconvolution enable quantification of
cell types in spatial transcriptomic data.* **Nature Communications** 2022;13(1):385. PMID 35046414.
doi:10.1038/s41467-022-28020-5
SpatialDecon and the `safeTME` matrix, built from 10,377 TCGA samples to exclude cancer-expressed genes. The
`is_pure_tumor` mechanism — deriving tumor profiles from PanCK+ AOIs and appending them to the profile matrix
— is the principled way to stop tumor spillover being scored as immune abundance.
**中文:** SpatialDecon 与 `safeTME` 矩阵（基于 10,377 例 TCGA 样本构建以排除癌细胞表达基因）。其 `is_pure_tumor`
机制（从 PanCK+ AOI 推导肿瘤表达谱并并入矩阵）是防止肿瘤渗漏被误计为免疫细胞丰度的规范做法。

**B3 · van Hijfte L, Geurts M, Vallentgoed WR, et al.** *Alternative normalization and analysis pipeline to
address systematic bias in NanoString GeoMx Digital Spatial Profiling data.* **iScience** 2023;26(1):105760.
PMID 36590163. doi:10.1016/j.isci.2022.105760
Benchmarks Q3 against modified CPM, DESeq2 size factors, gamma-fit and quantile normalization on three
criteria (inter-sample distribution similarity, MA-plot deviation, SNR–logFC correlation). Q3 and other
size-factor methods fail; **quantile normalization scored best**. Cite whenever you deviate from Q3, or when
you keep Q3 and need to acknowledge its limits.
**中文:** 在三项标准（样本间分布相似性、MA 图偏离、信噪比与 logFC 的相关）下比较 Q3 与改良 CPM、DESeq2 size factor、
gamma 拟合和分位数归一化：Q3 及其他 size factor 方法均不达标，**分位数归一化最优**。偏离 Q3 时应引用；坚持用 Q3 时
也应引用以说明其局限。

**B4 · Hoffman GE, Roussos P.** *dream: powerful differential expression analysis for repeated measures
designs.* **Bioinformatics** 2021;37(2):192–201. PMID 32730587. doi:10.1093/bioinformatics/btaa687
Gene-level linear mixed models with precision weights and Kenward–Roger inference. Use when
`duplicateCorrelation`'s single consensus correlation is too crude: crossed random effects (patient **and**
slide), very unbalanced ROI counts, or genes whose patient-clustering differs from the panel average.
**中文:** 带精度权重与 Kenward–Roger 推断的基因级线性混合模型。当 `duplicateCorrelation` 的单一共识相关过于粗糙时使用：
交叉随机效应（患者**与**玻片）、ROI 数严重不平衡、或某些基因的患者聚集程度明显偏离全panel 平均。

**B5 · Zimmerman KD, Espeland MA, Langefeld CD.** *A practical solution to pseudoreplication bias in
single-cell studies.* **Nature Communications** 2021;12(1):738. PMID 33531494. doi:10.1038/s41467-021-21038-1
Written for single-cell data but the argument is identical for multi-ROI DSP: treating sub-samples of the same
subject as independent inflates type-I error dramatically; pseudobulk or mixed models fix it. The cleanest
citation for playbook §1.2.
**中文:** 虽针对单细胞数据，但论证与多 ROI 的 DSP 完全一致：把同一受试者的子样本当作独立观测会大幅抬高 I 类错误，
pseudobulk 或混合模型可解决。是手册 §1.2 最直接的引用。

**B6 · Bergholtz H, Carter JM, Cesano A, et al.** *Best Practices for Spatial Profiling for Breast Cancer
Research with the GeoMx Digital Spatial Profiler.* **Cancers** 2021;13(17):4456. PMID 34503266.
doi:10.3390/cancers13174456
Community best-practice paper; states plainly that cohort-level DSP studies sampling multiple ROIs per patient
require mixed-effect models and that t-tests/ANOVA should be limited to designs without within-slide
replication — including for TMAs spanning slides or sites.
**中文:** 社区最佳实践文章，明确指出：每位患者取多个 ROI 的队列型 DSP 研究必须使用混合效应模型；t 检验/ANOVA 仅适用于
玻片内无重复的设计；即便是跨玻片、跨中心的 TMA，也应使用混合效应模型以吸收玻片间变异。

**B7 · Altman DG, Lausen B, Sauerbrei W, Schumacher M.** *Dangers of using "optimal" cutpoints in the
evaluation of prognostic factors.* **JNCI** 1994;86(11):829–835. PMID 8182763.
The original demonstration that minimum-p-value cut points produce inflated significance and biased effect
sizes. Still the citation of record for playbook §8.3.
**中文:** 最早证明 "最小 p 值切点" 会造成显著性虚高与效应量偏倚的文献，至今仍是手册 §8.3 的标准引用。

**B8 · Polley MC, Dignam JJ.** *Statistical Considerations in the Evaluation of Continuous Biomarkers.*
**Journal of Nuclear Medicine** 2021;62(5):605–611. PMID 33579807. doi:10.2967/jnumed.120.251520
Modern, readable treatment: keep biomarkers continuous, use splines for non-linearity, avoid the
"select a cut point then test survival difference in the same data" loop, and understand how optimistic bias
scales with the number of variables considered.
**中文:** 现代、易读的综述：生物标志物应保持连续、非线性用样条处理、避免 "在同一数据里选切点再检验生存差异" 的闭环，
并说明乐观偏倚如何随考察变量数增加而放大。

**B9 · Foroutan M, Bhuva DD, Lyu R, Horan K, Cursons J, Davis MJ.** *Single sample scoring of molecular
phenotypes.* **BMC Bioinformatics** 2018;19(1):404. PMID 30400809. doi:10.1186/s12859-018-2435-4
`singscore`: rank-based, single-sample, stable when new samples are added. Preferred over cohort-relative
scoring for a biomarker you intend to lock and validate externally.
**中文:** `singscore`：秩基、单样本、加入新样本时已有分值不变。若打算锁定评分并做外部验证，优于依赖队列的相对评分方法。

**B10 · Hänzelmann S, Castelo R, Guinney J.** *GSVA: gene set variation analysis for microarray and RNA-seq
data.* **BMC Bioinformatics** 2013;14:7. PMID 23323831. doi:10.1186/1471-2105-14-7
Cohort-relative pathway scoring; fine for exploration, inappropriate as a locked biomarker because scores
change when the cohort changes.
**中文:** 依赖队列的通路评分；用于探索没问题，但因队列变化会改变分值，不适合作为锁定的生物标志物。

---

## C. Software, manuals and non-peer-reviewed sources / 软件、手册与非同行评议资源

Label these as such in a manuscript. **中文:** 在论文中应标明其性质。

| Resource / 资源 | Type / 类型 | Use / 用途 |
|---|---|---|
| `GeomxTools` / `GeoMxWorkflows` — *Analyzing GeoMx-NGS RNA Expression Data with GeomxTools* (Bioconductor workflow vignette) | Software vignette | Canonical QC thresholds, LOQ formula, `mixedModelDE` with `(1 + testRegion \| slide)` random-slope guidance |
| `standR` + `GeoMXAnalysisWorkflow` (Davis Laboratory, Bioconductor) | Software | Implementation of A1; `findNCGs`, `geomxBatchCorrection`, RUV4 |
| `GeoDiff` (Bioconductor) — Poisson background model, background score test, NB-threshold size factors; preprint bioRxiv 2022.05.26.493637 | Software + preprint | Background-aware QC/normalization when many targets sit near LOQ (small CD68+ AOIs) |
| `SpatialDecon` (Bioconductor), `safeTME`, `is_pure_tumor` | Software | Implementation of B2 |
| `variancePartition::dream`, `lme4`, `lmerTest`, `emmeans` | Software | Mixed models, Satterthwaite/Kenward–Roger inference, contrasts |
| `survival` (`coxph`, `cluster()`, `frailty()`, `cox.zph`), `coxme`, `rms`, `timeROC`, `maxstat` | Software | §8 survival stack |
| GeoMx DSP Data Analysis User Manual (MAN-10154-01, Bruker/NanoString) | Vendor manual | Source of the "TMA with 1 ROI/patient → t-test; multiple ROIs/patient → LMM" rule and the random-slope guidance |
| GeoMx WTA design/performance white paper (Bruker/NanoString) | Vendor white paper | Panel content: 19,505 targets, 99.5% of HUGO protein-coding genes, 23 high-expressors intentionally removed |
| GEO **GSE271689** / **GSE292098** | Data | The GeoMx DSP arms of A4. Note that the GEO protocol text for GSE271689 lists S100B as the tumor morphology marker while the Nature Genetics methods describe PanCK for the NSCLC cohorts — reconcile against the published methods before reusing the annotation verbatim / GSE271689 的 GEO 协议文本把肿瘤形态学标记写作 S100B，而 Nature Genetics 方法学部分对 NSCLC 队列描述的是 PanCK；直接复用注释前请以正式发表的方法学为准 |
| ASCO/JCO conference abstracts on TACSTD2 tumor-vs-TME DSP contrasts (e.g. anal cancer, *J Clin Oncol* 2026;44(16_suppl):3519) | Abstract | Useful magnitude anchor (reported log fold change ≈ 3.49 tumor vs TME segments, adjusted p < 0.0001; DSP-vs-IHC correlation r ≈ 0.43) but abstract-level evidence — do not cite as a primary methods reference / 可作量级参考，但属摘要级证据，不应作为主要方法学引用 |

---

## D. How to cite this playbook's choices in a methods section / 如何在方法学部分引用

A defensible minimal citation set for a GSE271689-class analysis:

- QC/normalization/DE pipeline → **A1** (+ **B3** or **A2** if you deviate from Q3)
- Patient random effect / pseudoreplication → **B5**, **B6**, GeoMx DSP User Manual (Section C)
- Mixed-model DE implementation → **A1** (voom + duplicateCorrelation) or **B4** (dream)
- Compartment spillover control → **B2** (+ **A3** for the formal treatment)
- Target biology priors → **A8** (TACSTD2), **A9** (CLDN4), **A10** (heterogeneity/sampling)
- Survival with repeated ROIs → **A5** (`coxme` precedent), **B7**/**B8** (cut points, continuous biomarkers)
- Design/validation strategy → **A4**

**中文:** GSE271689 类分析的最小可辩护引用集合：质控/归一化/DE 流程用 **A1**（偏离 Q3 时补 **B3** 或 **A2**）；
患者随机效应与伪重复用 **B5**、**B6** 及 GeoMx 用户手册；混合模型实现用 **A1** 或 **B4**；分区渗漏对照用 **B2**
（正式处理见 **A3**）；靶点生物学先验用 **A8**、**A9**、**A10**；多 ROI 生存分析用 **A5** 与 **B7**/**B8**；
设计与验证策略用 **A4**。
