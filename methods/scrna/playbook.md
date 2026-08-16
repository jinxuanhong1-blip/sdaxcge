# Lung-Cancer ICI scRNA-seq Playbook (2024–2026)
# 肺癌免疫检查点治疗单细胞转录组最佳实践手册（2024–2026）

**Focus / 主题:** Tumor epithelial **TACSTD2 / CLDN4 / junction** program vs. immune outcomes
(CD8 cytotoxicity, TLS, CXCL13) under immune checkpoint inhibition (ICI), joined to
**MPR / RECIST** clinical response.

> **Scope note / 范围说明.** This is a *methods* deliverable. Everything lives under
> `methods/scrna/`. It contains (1) this bilingual playbook, (2) runnable Python/R
> templates in `scripts/`, (3) environment notes in `env/`, and (4) a minimal demo in
> `demo/` run on real GEO data (**GSE207422**). No statistics anywhere in this repo are
> fabricated — every number in `demo/` is computed from the downloaded matrix. Where a
> step is illustrative (e.g. small-n response comparisons), it is labeled *descriptive*.

---

## 0. TL;DR decision tree / 快速决策树

| Question / 问题 | Recommended default (2024–2026) / 推荐默认 |
|---|---|
| Ambient RNA / 环境游离RNA | **CellBender** `remove-background` (raw `.h5`), fallback **SoupX** / SoupX |
| Doublets / 双细胞 | **scDblFinder** (R) or **DoubletFinder**; per-sample, before integration / 每样本、整合前 |
| Integration / 批次整合 | **scVI** (deep, counts) if compute allows; else **Harmony** (fast, PCA) / scVI 或 Harmony |
| Label transfer / 标签迁移 | **scANVI** / **CellTypist** (Immune All Low + custom) / CellTypist |
| **Primary DE** / 主要差异分析 | **Pseudobulk** (edgeR-QLF / DESeq2), sample = replicate / 伪bulk，**样本**为重复单位 |
| Cells-as-replicates Wilcoxon | **Only** for marker discovery, never for cross-condition inference / 仅用于marker，不用于组间推断 |
| Module scoring / 模块打分 | `AddModuleScore` / `sc.tl.score_genes` / **UCell** (rank-based) / UCell 更稳健 |
| Cell–cell comm / 细胞通讯 | **LIANA+** (consensus) or **CellChat v2**, **secondary/hypothesis-generating only** / 仅作次要假设生成 |
| Composition shift / 组成变化 | **scCODA** / **sccomp** / **milo** (differential abundance) / 差异丰度 |
| Response label join / 临床标签对接 | Join at **sample** level; keep MPR *and* RECIST separately / 在样本层对接，MPR 与 RECIST 分开 |

---

## 1. Study design & the response-label problem / 研究设计与响应标签

### 1.1 What the labels mean / 标签含义
- **RECIST 1.1** (radiologic, on measurable disease): CR / PR / SD / PD. Often collapsed to
  **responder = CR+PR** vs **non-responder = SD+PD**, or ORR. Applies to advanced/metastatic
  and neoadjuvant-with-imaging settings.
  RECIST 为影像学标准，常合并为“应答(CR+PR)/非应答(SD+PD)”。
- **Pathologic response** (neoadjuvant/resection): **pCR** (0% viable tumor), **MPR**
  (≤10% residual viable tumor), **NMPR** (>10%). MPR/pCR are the accepted surrogate
  endpoints for neoadjuvant ICI in NSCLC (e.g. CheckMate-816). 病理缓解：pCR / MPR(≤10%残存) / NMPR。
- **They are not interchangeable.** A tumor can be RECIST-PR but pathologically NMPR
  (the GSE207422 metadata shows exactly this discordance). Report both; pre-register which
  is primary. 二者并不等价，需分别报告并预先指定主要终点。

### 1.2 Where labels attach / 标签归属层级
Labels are **patient/sample-level**, not cell-level. The single most common statistical
error in this field is treating cells as independent replicates of a *patient-level* label.
Keep a tidy sample-metadata table keyed by `sample_id` and join once (see
`scripts/07_join_clinical.py` / `scripts/07_join_clinical.R`).
标签属于患者/样本层级；**切勿把细胞当作独立重复**。用 `sample_id` 维护整洁的样本表，仅在一处 join。

Minimum columns / 最少字段:
`sample_id, patient_id, timepoint(pre/post), tissue_site, histology, ici_agent,
chemo, recist, path_response(pCR/MPR/NMPR), residual_tumor_pct, batch, n_cells`.

### 1.3 Timepoint & site confounds / 时间点与取样部位混杂
- Pre- vs post-treatment and biopsy vs resection change composition dramatically. Do **not**
  pool timepoints in one DE contrast without a covariate. 治疗前/后、活检/手术差异巨大，勿直接合并。
- Metastatic-site vs primary samples (GSE205335 spans multiple sites) must be modeled or
  stratified. 多部位样本需分层或建模。

---

## 2. Ingestion & QC / 数据读取与质控

### 2.1 Ambient RNA / 环境RNA (soup)
- **CellBender** `remove-background` on the **raw** (unfiltered) 10x `.h5`. It also calls
  cells; run per capture lane. Save the cleaned `.h5` and use it downstream.
- **SoupX** alternative when only filtered + raw are available; estimates contamination
  fraction and subtracts. Always sanity-check that lineage markers (e.g. HBB, ALB, SFTPC)
  stop bleeding into other clusters after correction. 校正后确认谱系marker不再“外溢”。
- Pitfall / 陷阱: over-aggressive correction removes real lowly-expressed genes; keep the
  correction fraction plot and eyeball it. Ambient removal is **per-lane**, before merging.

### 2.2 Doublets / 双细胞
- **scDblFinder** (Bioconductor, fast, well-benchmarked) or **DoubletFinder**. Run
  **per sample/lane** *before* integration — never on the integrated object, and never
  after aggressive filtering that distorts the neighborhood. 每个样本单独跑、整合前。
- Combine with biological sanity: clusters co-expressing exclusive lineages (EPCAM+ &
  CD3E+) are suspect. 双谱系共表达簇需警惕。

### 2.3 Per-cell QC thresholds / 单细胞质控阈值
- Filter on `n_genes_by_counts`, `total_counts`, `pct_counts_mt`, `pct_counts_ribo`,
  and (for nuclei) `pct_counts_hb`. **Set thresholds per-dataset from the distributions,
  not by copy-pasted magic numbers.** Prefer **MAD-based** cutoffs (e.g. `scater::isOutlier`,
  median ± 3–5 MAD on log scale) computed **per sample**. 用每样本 MAD 自适应阈值，而非固定数字。
- Typical starting window (tumor tissue, 3' 10x): `200 ≤ n_genes ≤ 6000–8000`,
  `pct_mt < 15–20%` (tumor/necrotic tissue tolerates higher mito than blood). 仅为起点，需按分布调整。
- Doublet-suspected high-`total_counts` tail: handle with the doublet caller, not a hard cap.

### 2.4 Sex/quality covariates / 协变量
Compute cell-cycle scores (`sc.tl.score_genes_cell_cycle`) and, if needed, regress or
supply as scVI covariate. XIST/Y-genes can flag sample swaps. 计算细胞周期，必要时作为协变量。

---

## 3. Normalization, HVG, integration / 归一化、高变基因、整合

### 3.1 Normalization / 归一化
- Standard: `normalize_total(1e4)` + `log1p` for visualization/scoring. For model-based DE
  we go back to **raw counts** (pseudobulk). 可视化/打分用 log-CP10k；DE 用原始 counts。
- **Analytic Pearson residuals** / **sctransform v2** are good alternatives for HVG & PCA.

### 3.2 HVG / 高变基因
- 2000–3000 HVGs, **computed with `batch_key`** so batch-specific genes don't dominate.
  Exclude TCR/BCR variable genes (TRAV/TRBV/IGHV…) and, if studying tumor programs,
  optionally mito/ribo/heat-shock. 计算高变基因时加 batch_key，并剔除 TCR/BCR 可变区基因。

### 3.3 Integration / 整合
- **scVI** (scvi-tools): probabilistic, integrates on **raw counts**, handles batch +
  continuous covariates, yields a latent space robust for clustering and a denoised
  expression for visualization. Best when you have >20–30 samples or strong batch. Set a
  fixed seed; monitor ELBO. scVI：基于 counts 的深度整合，样本多/批次强时首选。
- **Harmony**: fast PCA-space correction; excellent default when compute is limited or the
  batch is mild. Harmony：PCA 空间快速校正，算力有限时的稳妥默认（本仓库 demo 即用 Harmony）。
- **scANVI** extends scVI with partial labels for semi-supervised integration + label
  transfer. scANVI 可半监督整合并迁移标签。
- Evaluate integration with **scIB** metrics (batch mixing vs bio conservation: kBET,
  iLISI/cLISI, ARI, silhouette). Don't over-integrate — you can erase the tumor-vs-immune
  or responder biology you care about. 用 scIB 指标平衡“批次混合”与“生物学保留”，勿过度整合。

Pitfall / 陷阱: integrating **malignant epithelial** cells across patients is dangerous —
CNV-driven patient-specific programs are real biology, not batch. Common practice: cluster
all cells for lineage, then analyze malignant cells **per patient** or with methods that
respect patient structure (see §5). 恶性上皮细胞跨患者整合需谨慎，CNV 程序是真实信号而非批次。

---

## 4. Annotation / 注释

### 4.1 Strategy / 策略
1. Cluster on the integrated latent space (**Leiden**, sweep resolution 0.4–1.5).
2. **Coarse lineage** by canonical markers (below), then **subcluster** each lineage.
3. **Automated priors** to speed/steady it: **CellTypist** (`Immune_All_Low`,
   `Immune_All_High`), **scANVI** label transfer from a reference (e.g. a NSCLC atlas),
   or Azimuth. Always reconcile automated calls with markers and DE. 自动注释需与marker互证。

### 4.2 Canonical lung-TME markers / 常用肺TME标记
| Lineage / 谱系 | Markers |
|---|---|
| T/NK | CD3D, CD3E, TRAC, CD8A, CD4, IL7R, FOXP3, NKG7, GNLY, KLRD1 |
| B / Plasma | MS4A1, CD79A, CD79B, BANK1 / MZB1, IGHG1, XBP1 |
| Myeloid | LYZ, CD68, CD14, FCGR3A, C1QA/B/C, SPP1, MRC1, FCN1 |
| DC | CLEC9A (cDC1), CD1C (cDC2), LAMP3 (mregDC), LILRA4 (pDC) |
| Mast | TPSAB1, TPSB2, CPA3, MS4A2 |
| Epithelial (normal+malignant) | EPCAM, KRT8/18/19, CDH1, SFTPC/SFTPB (AT2), SCGB1A1 (club), FOXJ1 (ciliated) |
| **Tumor program of interest** | **TACSTD2 (TROP2), CLDN4, CLDN3/7, CLDN18, ELF3, KRT17** |
| Endothelial | PECAM1, VWF, CLDN5, CLEC14A |
| Fibroblast / CAF | COL1A1/2, DCN, LUM, PDGFRB, ACTA2 (myCAF), FAP |
| TLS / follicular | CXCL13, CCL19, CCL21, CR2, CXCR5, LTB |

### 4.3 Malignant-cell identification / 恶性上皮判定
- Separate malignant from normal epithelial via **inferCNV** / **CopyKAT** / **numbat**
  (numbat adds allele info). Reference = immune/stromal cells. Malignant cells show broad
  CNV; this also explains why they don't integrate across patients. 用 inferCNV/CopyKAT/numbat 判定恶性并解释患者特异性。

---

## 5. Differential expression — pseudobulk first / 差异表达：伪bulk优先

> **This is the single most important methodological stance in this playbook.**
> **这是本手册最重要的方法学立场。**

### 5.1 Why not Wilcoxon-on-cells as primary / 为何不以“细胞级Wilcoxon”为主
- Cells from one patient are **not** independent replicates of a patient-level label
  (MPR/RECIST). Treating them as such inflates the effective n from ~15 patients to ~90,000
  cells, producing **anti-conservative** p-values and near-guaranteed "significance." This is
  documented pseudoreplication (Squair et al., *Nat Commun* 2021; Murphy & Skene 2022;
  Zimmerman et al. 2021; Crowell et al. *muscat* 2020).
  一个患者的细胞并非该患者标签的独立重复；细胞级检验把 n 从~15患者夸大到~9万细胞，导致 p 值反保守、几乎必然“显著”。
- Cell-level Wilcoxon is fine for **within-dataset marker discovery** (what defines a
  cluster) — that's an internal, descriptive use, not cross-condition inference. 细胞级检验仅适用于簇marker发现。

### 5.2 Pseudobulk recipe / 伪bulk流程
1. Choose the cell population (e.g. **malignant epithelial**, or **CD8 T**). 选定细胞群。
2. **Sum raw counts** per `sample_id` within that population → sample × gene matrix.
   `decoupler.get_pseudobulk` / `muscat::aggregateData` / manual `groupby.sum`. 按样本求和原始counts。
3. Filter low-count samples (e.g. `< 10 cells` in that population; `< 10` counts genes via
   `edgeR::filterByExpr`). 过滤细胞过少的样本与低表达基因。
4. Fit **edgeR-QLF** or **DESeq2** with the design
   `~ batch + covariates + response`. Use TMM/median-of-ratios normalization. 用 edgeR/DESeq2 建模。
5. Report log2FC, FDR (BH), and **show per-sample points** — not a volcano built from cells.
   报告 log2FC、FDR，并展示每样本点图。

### 5.3 Design & power / 设计与效能
- With ~15 patients split MPR/NMPR you have **low power**; be honest. Add covariates only if
  df allow. Consider histology stratification (Adeno vs Squamous). Consider
  **timepoint as blocking** if paired. 小样本效能低，需诚实报告；必要时按组织学分层、按配对设区组。
- For continuous outcomes (residual tumor %), use a continuous predictor rather than a
  dichotomy. 连续终点(残存肿瘤%)优先用连续预测变量。

### 5.4 Alternatives that respect hierarchy / 尊重层级的替代法
- **muscat** (pseudobulk + mixed models), **NEBULA** / **glmmTMB** mixed models when you
  truly need cell-level with random patient effects, **Milo/scCODA/sccomp** for composition.

`scripts/05_pseudobulk_de.py` (decoupler) and `scripts/05_pseudobulk_de.R`
(edgeR + DESeq2) implement this; `demo/pseudobulk_de.R` runs it on the demo matrix.

---

## 6. Epithelial gene-module scoring / 上皮基因模块打分

### 6.1 The TACSTD2 / CLDN4 / junction module / 目标模块
Rationale / 依据: **TACSTD2 (TROP2)** is the target of the ADC sacituzumab govitecan and a
marker of epithelial/tumor cells; **CLDN4** and other claudins (CLDN3/7/18) + tight/adherens
junction genes describe an epithelial-integrity / low-EMT program relevant to barrier
function and immune exclusion. Suggested gene set (edit to taste):
建议基因集（可调整）:

```
TACSTD2, CLDN4, CLDN3, CLDN7, CLDN18, CDH1, TJP1 (ZO-1), OCLN, F11R (JAM-A), EPCAM, ELF3
```

### 6.2 Scoring methods / 打分方法
- `sc.tl.score_genes` (Scanpy) / `AddModuleScore` (Seurat): control-gene-corrected mean.
  Fast, but sensitive to library size and the random control bins. 需注意文库大小影响。
- **UCell** (rank-based, per-cell Mann–Whitney U statistic): **robust to depth and
  normalization**, recommended for cross-sample comparison. 推荐 UCell（基于排序，更稳健）。
- **AUCell** for "is the program on" (binary-ish activity). 
- Always score on **malignant epithelial cells only** for a tumor program, and compare at
  the **sample level** (mean/median per sample), then relate to response. 仅在恶性上皮细胞打分，样本层比较。

### 6.3 Pitfalls / 陷阱
- Don't compare module scores across cell types (scores aren't calibrated across very
  different transcriptomes). 不同细胞类型间的分数不可直接比较。
- Depth confound: verify the score isn't just tracking `total_counts` (plot score vs depth).
  验证分数并非只反映测序深度。
- Regressing/using scVI-denoised expression can stabilize scores in sparse data.

`scripts/06_module_scores.py` and `scripts/06_module_scores.R` (UCell) provide templates.

---

## 7. Immune neighborhood scores: CD8 / TLS / CXCL13 / 免疫邻域评分

### 7.1 Signatures / 信号集 (starting points, cite & adapt / 起点，需引用与调整)
- **CD8 cytotoxicity / effector:** CD8A, CD8B, GZMB, GZMK, GZMH, PRF1, IFNG, NKG7, KLRG1.
  **Exhaustion / dysfunction:** PDCD1, HAVCR2, LAG3, TIGIT, CTLA4, TOX, ENTPD1 (CD39). CD39+
  (ENTPD1) CD8 marks tumor-reactive T cells. CD8毒性/耗竭信号。
- **Tissue-resident memory (T_RM):** ITGAE (CD103), ZNF683, CXCR6 — associated with ICI benefit.
- **TLS (tertiary lymphoid structures):** CXCL13, CCL19, CCL21, CR2, CXCR5, LTB, SELL, plus
  B-cell (MS4A1) + follicular DC signals. A widely used **12-chemokine TLS signature**
  exists; also the Cabrita/Meylan B-cell/TLS signatures. TLS 与 ICI 获益相关。
- **CXCL13:** a single, high-value gene — marks T follicular-helper-like / exhausted CD8 and
  predicts ICI response in several cohorts. Track both as a **module** and a **single gene**.
  CXCL13 既作模块也作单基因追踪。

### 7.2 How to use / 使用方式
- Score CD8 signatures **within CD8 T cells**; TLS signatures are best assessed as
  **spatial/neighborhood** features — in dissociated scRNA you approximate via co-occurrence
  of B, T-fh (CXCL13+), and CCL19/21 stroma at the **sample** level. TLS 本质是空间结构，解离scRNA只能近似。
- Relate **sample-level** mean scores to MPR/RECIST (descriptive with small n), and/or use
  as pseudobulk covariates. 样本层评分关联响应；小样本仅作描述。
- **Neighborhood analysis proper** needs spatial data (Visium/Xenium/CosMx) or, within
  scRNA, **niche inference** methods; state clearly when you're only approximating. 邻域分析严格需空间数据。

---

## 8. Cell–cell communication — SECONDARY only / 细胞通讯：仅作次要

> CCC results are **hypothesis-generating**, not confirmatory. Report them after DE/scoring,
> clearly flagged. CCC 结果仅用于生成假设，不作确证，须明确标注。

- **LIANA+** (Python/R): runs multiple methods (CellPhoneDB, NATMI, SingleCellSignalR,
  Connectome, logFC) and a **consensus rank**; integrates with pseudobulk and MOFA. Preferred
  for robustness. LIANA+ 提供多方法共识排名，稳健性更好。
- **CellChat v2**: pathway-level, nice visualization, includes spatial mode.
- **Good practice / 良好实践:**
  - Compare conditions with **differential** CCC (e.g. LIANA on MPR vs NMPR), not just one
    global network. 用差异CCC比较组别。
  - Aggregate to the population you trust; sparse/rare cell types give noisy L-R calls.
  - Validate priorities against DE and, if possible, spatial proximity. 与DE/空间邻近互证。
  - **Pitfall:** L-R "significance" from permutation over cells inherits the same
    pseudoreplication risk — treat p-values skeptically and confirm at sample level. 排列检验同样有伪重复风险。

`scripts/08_ccc_liana.py` and `scripts/08_ccc_cellchat.R` are templates (not run in demo).

---

## 9. Joining MPR / RECIST and testing associations / 对接临床标签并检验关联

1. Build `sample_metadata.csv` (see §1.2) once; validate every `sample_id` matches the
   AnnData/Seurat `obs`. 先构建样本元数据并校验ID匹配。
2. Join at sample level. For each readout (module score, pseudobulk gene, composition):
   - **Categorical** (MPR vs NMPR; responder vs non): Wilcoxon/​t at the **sample** level, or
     logistic regression `response ~ score + covariates`. 样本层Wilcoxon/逻辑回归。
   - **Continuous** (residual tumor %): Spearman / linear model. 连续终点用Spearman/线性模型。
   - **Ordinal RECIST**: ordered logistic or collapse to responder/non a priori. 有序RECIST需预先处理。
3. **Multiple testing:** correct across genes/signatures (BH). Report effect sizes + CIs, not
   just p. 多重检验校正，报告效应量与置信区间。
4. **Honesty with n:** with ~4 MPR vs ~8 NMPR, most tests are underpowered; present as
   descriptive/exploratory and, where possible, replicate in an independent cohort
   (GSE207422 ↔ GSE205335 / bulk TCGA-LUAD / OAK-POPLAR). 小样本需诚实，尽量外部验证。

---

## 10. Reproducibility & reporting / 可复现与报告

- Pin versions (`env/`), set seeds, log package versions in outputs. 固定版本与随机种子。
- Save the integrated object (`.h5ad`/`.rds`), the pseudobulk matrix, and the sample design
  table as first-class artifacts. 保存整合对象、伪bulk矩阵、样本设计表。
- Report: cells before/after QC per sample, doublet/ambient parameters, integration method +
  scIB metrics, DE design formula, exact gene sets with a citation, and every n. 报告全部关键参数与n。
- Follow community reporting norms; deposit code. 遵循社区报告规范并公开代码。

---

## 11. Common pitfalls checklist / 常见陷阱清单

- [ ] Cells-as-replicates for a patient-level contrast (pseudoreplication). 细胞当重复。
- [ ] Over-integration erasing tumor/responder biology. 过度整合抹去生物学。
- [ ] Integrating malignant cells across patients as if CNV programs were batch. 恶性细胞跨患者整合。
- [ ] Copy-pasted QC thresholds instead of per-sample MAD. 固定阈值。
- [ ] Module scores confounded by sequencing depth. 打分受深度混杂。
- [ ] Comparing module scores across cell types. 跨细胞类型比分数。
- [ ] Ambient correction per merged object instead of per lane. 合并后再去环境RNA。
- [ ] Treating CCC permutation p-values as confirmatory. CCC p 值当确证。
- [ ] Mixing RECIST and pathologic response as if identical. 混用RECIST与病理缓解。
- [ ] Pooling pre/post or multi-site samples without a covariate. 合并时点/部位不加协变量。
- [ ] No multiple-testing correction across signatures/genes. 不做多重检验校正。
- [ ] Overclaiming from n≈15 without external replication. 小样本过度声称。

---

## 12. Key tools & references / 关键工具与参考

**Tools (current) / 工具:** Scanpy, AnnData, scvi-tools (scVI/scANVI), Harmony/harmonypy,
scDblFinder, DoubletFinder, CellBender, SoupX, CellTypist, inferCNV/CopyKAT/numbat,
edgeR, DESeq2, muscat, decoupler, UCell, AUCell, LIANA+, CellChat v2, scCODA/sccomp, Milo,
scIB, Seurat v5.

**Method references / 方法学参考 (verify DOIs before citing / 引用前核对):**
- Squair et al. "Confronting false discoveries in single-cell differential expression."
  *Nat Commun* 12:5692 (2021). — pseudoreplication.
- Crowell et al. "muscat detects subpopulation-specific state transitions…" *Nat Commun*
  11:6077 (2020). — pseudobulk DS.
- Zimmerman, Espeland, Langefeld. "A practical solution to pseudoreplication bias in
  single-cell studies." *Nat Commun* 12:738 (2021).
- Lopez et al. "Deep generative modeling for single-cell transcriptomics (scVI)."
  *Nat Methods* 15:1053 (2018); Xu et al. scANVI, *Mol Syst Biol* 2021.
- Korsunsky et al. "Fast, sensitive and accurate integration with Harmony." *Nat Methods*
  16:1289 (2019).
- Germain et al. scDblFinder, *F1000Research* (2021/2022).
- Fleming et al. "CellBender remove-background." *Nat Methods* 20:1323 (2023).
- Young & Behjati. "SoupX." *GigaScience* 9:giaa151 (2020).
- Domínguez Conde et al. "CellTypist / Cross-tissue immune cell atlas." *Science* 376 (2022).
- Andreatta & Carmona. "UCell." *Comput Struct Biotechnol J* 19:3796 (2021).
- Dimitrov et al. "LIANA / comparison of CCC methods." *Nat Commun* 13:3224 (2022);
  LIANA+ *Nat Cell Biol* 2024.
- Jin et al. "CellChat." *Nat Commun* 12:1088 (2021); CellChat v2 (2024/2025).
- Büttner et al. "scIB benchmarking of integration." *Nat Methods* 19:41 (2022).
- Luecken & Theis. "Current best practices in single-cell RNA-seq analysis: a tutorial."
  *Mol Syst Biol* 15:e8746 (2019) + sc-best-practices.org (2023–2025, living book).

**Datasets used/considered / 数据集:**
- **GSE207422** — Hu et al. "Tumor microenvironment remodeling after neoadjuvant
  immunotherapy in NSCLC…" *Genome Medicine* 15:14 (2023). 15 pts, neoadjuvant PD-1+chemo,
  MPR/NMPR/pCR + RECIST. **Used for the demo** (`demo/`).
- **GSE205335** — "Single-cell transcriptome profiles… lung cancer receiving ICI"; Kim et al.
  *eLife* 12:RP98366 (2024). Processed matrix is a **499 MB `.rds`** (< 2 GB); raw in EGA.
  Response = responder(PR)/non-responder(SD/PD). Suitable as an external replication cohort.

> Citations are provided for orientation; **verify each DOI/version at use time**. 引用请在使用时核实。
