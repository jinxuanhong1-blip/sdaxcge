<!--
methods/cnv/playbook.md
METHODS ONLY. No results, no interpretation of a specific dataset.
This is a reusable TEMPLATE. Fill in every {{PLACEHOLDER}} before use.
Bilingual: English (EN) and 简体中文 (ZH) are provided section-by-section.
本文件仅描述方法（METHODS ONLY），不包含任何具体数据结果，是可复用的模板。
-->

# Playbook — Separating Malignant vs Normal Epithelium before Scoring TACSTD2 / CLDN4
# 操作手册 —— 在评分 TACSTD2 / CLDN4 之前区分恶性上皮与正常上皮

**Scope / 适用范围:** single-cell (and single-nucleus) RNA-seq. CNV inference with **inferCNV**, **CopyKAT**, **Numbat**.
**Type / 类型:** METHODS ONLY template. 仅方法模板。
**Status / 状态:** `{{DRAFT | REVIEWED | LOCKED}}`

---

## 0. Template metadata / 模板元数据

```yaml
# Fill this block per run. 每次运行填写。
project_id:        "{{PROJECT_ID}}"
analyst:           "{{NAME}}"
date:              "{{YYYY-MM-DD}}"
assay:             "{{scRNA-seq | snRNA-seq}}"
platform:          "{{10x 3' v3 | 10x 5' | Smart-seq2 | ...}}"
tissue:            "{{tumor site}}"
indication:        "{{cancer type}}"
treatment_context: "{{treatment-naive | ICI on-treatment | post-progression}}"
timepoint:         "{{baseline | on-treatment C{{n}}D{{n}} | EOT}}"
targets_scored:    ["TACSTD2", "CLDN4"]     # Trop-2 / Claudin-4
cnv_tools:         ["inferCNV", "CopyKAT", "Numbat"]
genome_build:      "{{GRCh38 | GRCh37}}"
software_versions:
  infercnv:        "{{x.y.z}}"
  copykat:         "{{x.y.z}}"
  numbat:          "{{x.y.z}}"
  cellsnp_lite:    "{{x.y.z}}"      # for Numbat allele counts
  eagle2:          "{{x.y.z}}"      # for Numbat phasing
```

---

## 1. Purpose & rationale / 目的与原理

**EN.** TACSTD2 (Trop-2) and CLDN4 (Claudin-4) are antibody–drug-conjugate (ADC) relevant surface targets that are **also expressed by normal/reactive epithelium**. Any expression score computed over an undifferentiated "epithelial" compartment therefore mixes malignant and normal signal and is not interpretable as a tumor-target metric. The purpose of this playbook is to define a reproducible procedure that (a) identifies bona-fide malignant epithelial cells using copy-number-variation (CNV) evidence from three orthogonal tools, (b) reconciles the calls into a consensus label, and (c) restricts TACSTD2/CLDN4 scoring to confirmed malignant epithelium, with an internal normal-epithelium comparison.

**ZH.** TACSTD2（Trop-2）与 CLDN4（Claudin-4）是抗体偶联药物（ADC）相关的表面靶点，但**正常/反应性上皮同样表达**。若在未区分的"上皮"整体上直接计算表达评分，会把恶性与正常信号混在一起，无法作为肿瘤靶点指标解读。本手册的目的是给出可复现流程：(a) 用三种正交工具的拷贝数变异（CNV）证据鉴定真正的恶性上皮细胞；(b) 将三者结果整合为共识标签；(c) 仅在确认的恶性上皮上评分 TACSTD2/CLDN4，并与内部正常上皮作对照。

**Non-goals / 非目标:** clinical decision-making, cutoff certification, or dataset-specific conclusions. 不涉及临床决策、阈值认证或针对某数据集的结论。

---

## 2. Inputs & prerequisites / 输入与前置条件

**EN.**
- Per-sample **raw UMI count matrix** (genes × cells). Keep raw counts; CNV tools must not receive batch-integrated or scaled values.
- Cell-level QC already applied and documented (see §3). Record thresholds; do **not** re-QC silently.
- A cell-type annotation with at least: epithelial, T, B/plasma, myeloid, fibroblast, endothelial. Epithelial identity from `EPCAM`, `KRT8/18/19`, tissue-specific keratins.
- Gene-position/ordering file matching `genome_build` (chromosome, start, end) for inferCNV/CopyKAT.
- For **Numbat**: aligned BAM(s), a phasing reference panel (e.g., 1000G) and genetic map, plus an aggregated normal expression reference.
- **Run per sample**, never on a merged/integrated object (see §9 pitfalls). Aggregate only after per-sample calls.

**ZH.**
- 每个样本的**原始 UMI 计数矩阵**（基因 × 细胞）。保留原始计数；CNV 工具不得使用批次整合或缩放后的值。
- 已完成并记录细胞级 QC（见 §3）。记录所有阈值，**不要**私自重新 QC。
- 细胞类型注释，至少包含：上皮、T、B/浆细胞、髓系、成纤维、内皮。上皮身份用 `EPCAM`、`KRT8/18/19` 及组织特异角蛋白判定。
- 与 `genome_build` 匹配的基因位置/排序文件（染色体、起、止），供 inferCNV/CopyKAT 使用。
- **Numbat** 额外需要：比对后的 BAM、分型参考面板（如 1000G）与遗传图谱，以及聚合的正常表达参考。
- **按样本分别运行**，切勿在合并/整合对象上运行（见 §9 陷阱）。仅在得到各样本结果后再聚合。

---

## 3. QC gate (record, do not improvise) / QC 关卡（记录，勿临时改）

| Parameter / 参数 | Template value / 模板值 | Note / 说明 |
|---|---|---|
| min genes/cell `nFeature` | `{{≥ 500}}` | tissue-dependent 视组织而定 |
| max mito % | `{{≤ 15–20%}}` | on-treatment biopsies run higher; do not blanket-relax (§9) |
| min counts/cell | `{{≥ 1000}}` | |
| doublet removal | `{{scDblFinder | DoubletFinder}}` | tumor–immune doublets bias CNV (§9) |
| ambient RNA correction | `{{CellBender | SoupX}}` | **required** before target scoring (§8) |

**EN.** Ambient-RNA correction is mandatory here because TACSTD2/CLDN4 are highly/secreted-membrane expressed and leak into droplets, inflating apparent positivity in non-epithelial cells. Correct ambient RNA for scoring, but consider running CNV tools on the pre-correction raw counts if the correction method distorts the count distribution — document the choice.

**ZH.** 此处必须做环境 RNA 校正：TACSTD2/CLDN4 高表达/膜分泌，易泄漏进液滴，导致非上皮细胞出现假阳性。评分用校正后计数；若校正会扭曲计数分布，则 CNV 工具可用校正前原始计数运行——务必记录选择。

---

## 4. Diploid reference selection / 二倍体参考细胞选择

**EN.** Choose confidently non-epithelial, non-malignant cells as the diploid baseline:
- **Use:** T cells, B/plasma cells, myeloid cells (immune) and/or fibroblasts + endothelial (stromal).
- **Do NOT use** any epithelial cells as reference.
- Ensure the reference is **present within the same sample** and reasonably abundant (`{{≥ 100–200}}` cells). If a sample lacks reference cells, flag it — CopyKAT's auto-baseline and inferCNV's relative signal both degrade.
- Keep reference composition **consistent across samples** in a study to avoid baseline-driven differences.

**ZH.** 选择确定为非上皮、非恶性的细胞作为二倍体基线：
- **使用：** T 细胞、B/浆细胞、髓系（免疫）和/或成纤维 + 内皮（间质）。
- **不得**使用任何上皮细胞作参考。
- 参考细胞需**存在于同一样本内**且数量充足（`{{≥ 100–200}}` 个）。若某样本缺乏参考细胞，需标记——CopyKAT 的自动基线与 inferCNV 的相对信号都会退化。
- 全研究中参考组成应**跨样本一致**，避免基线差异被误读为生物学差异。

---

## 5. Tool A — inferCNV / 工具 A —— inferCNV

**EN.** Expression-smoothing CNV inference; produces per-cell/subcluster CNV profiles relative to the reference.

```r
# TEMPLATE — adapt paths and params.
library(infercnv)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = "{{counts.tsv | matrix}}",
  annotations_file  = "{{cell_annotations.tsv}}",   # cell -> group label
  gene_order_file   = "{{gene_order.{{build}}.tsv}}",
  ref_group_names   = c("{{T}}", "{{Myeloid}}", "{{Fibroblast}}", "{{Endothelial}}")
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff              = {{0.1}},          # 0.1 for 10x, ~1 for Smart-seq
  out_dir             = "{{out/infercnv/SAMPLE}}",
  cluster_by_groups   = FALSE,            # cluster within observations
  analysis_mode       = "subclusters",    # resolve intratumor subclones
  tumor_subcluster_partition_method = "leiden",
  denoise             = TRUE,
  HMM                 = TRUE,
  HMM_type            = "i6",             # 6-state; use i3 if unstable
  num_threads         = {{N}}
)
```

- **Malignant call:** flag observation subclusters with clear arm-level CNV and elevated CNV burden vs reference. Optionally quantify per-cell CNV score (mean squared deviation of the smoothed signal) and threshold `{{describe method}}`.
- Save: CNV heatmap, HMM state matrix, per-cell CNV score → `infercnv_label ∈ {malignant, normal, ambiguous}`.

**ZH.** 基于表达平滑的 CNV 推断；输出相对参考的每细胞/亚簇 CNV 图谱。
- **恶性判定：** 标记出现明显臂级 CNV、且 CNV 负荷显著高于参考的观测亚簇。可计算每细胞 CNV 评分（平滑信号的均方偏差）并设阈值 `{{说明方法}}`。
- 保存：CNV 热图、HMM 状态矩阵、每细胞 CNV 评分 → `infercnv_label ∈ {malignant, normal, ambiguous}`。

---

## 6. Tool B — CopyKAT / 工具 B —— CopyKAT

**EN.** Bayesian segmentation with an integrative aneuploid/diploid classifier.

```r
# TEMPLATE
library(copykat)

ck <- copykat(
  rawmat        = {{raw_count_matrix}},     # genes x cells, raw UMIs
  id.type       = "S",                      # "S" symbol / "E" ensembl
  ngene.chr     = {{5}},
  win.size      = {{25}},
  KS.cut        = {{0.1}},                  # segmentation sensitivity
  distance      = "euclidean",
  norm.cell.names = c({{known_normal_barcodes}}),  # anchor the diploid baseline
  genome        = "{{hg20 | hg38}}",
  n.cores       = {{N}},
  sam.name      = "{{SAMPLE}}"
)
# ck$prediction$copykat.pred ∈ {aneuploid, diploid, not.defined}
```

- **Malignant call:** `copykat.pred == "aneuploid"` **and** epithelial identity.
- **Baseline safeguard:** always pass `norm.cell.names` from §4 so the diploid baseline is anchored on known immune/stromal cells rather than auto-detected — critical for low-purity ICI samples (§9).
- Handle `not.defined` as `ambiguous`.

**ZH.** 贝叶斯分段 + 整合式非整倍体/二倍体分类器。
- **恶性判定：** `copykat.pred == "aneuploid"` **且**具上皮身份。
- **基线保护：** 始终传入 §4 的 `norm.cell.names`，把二倍体基线锚定在已知免疫/间质细胞上，而非自动检测——对低纯度 ICI 样本尤为关键（§9）。
- `not.defined` 记为 `ambiguous`。

---

## 7. Tool C — Numbat / 工具 C —— Numbat

**EN.** Integrates **allele** (phased SNP) and expression evidence; strongest for low-purity samples and for detecting copy-neutral LOH.

```bash
# Step 1: allele counts per cell (cellsnp-lite) + phasing (Eagle2). TEMPLATE.
cellsnp-lite -s {{sample.bam}} -b {{barcodes.tsv}} -O {{out/cellsnp}} \
  -R {{1000G_snps.vcf.gz}} -p {{N}} --minMAF 0.1 --minCOUNT 20 --UMItag {{UB}} --cellTAG {{CB}}
# then phase with Eagle2 + reference panel + genetic map (numbat helper).
```

```r
# Step 2: run Numbat. TEMPLATE.
library(numbat)
out <- run_numbat(
  count_mat  = {{raw_count_matrix}},        # genes x cells
  lambdas_ref= {{aggregated_normal_ref}},   # expression reference profile
  df_allele  = {{phased_allele_df}},         # from step 1
  genome     = "{{hg38}}",
  t          = 1e-5,
  ncores     = {{N}},
  out_dir    = "{{out/numbat/SAMPLE}}"
)
# outputs: clone posterior, p_cnv (aneuploidy prob) per cell, phylogeny
```

- **Malignant call:** high `p_cnv` / assignment to an aneuploid clone (not the normal clone), **and** epithelial identity.
- Allelic evidence makes Numbat the tie-breaker when expression-only tools disagree, and it is less prone to the interferon/MHC artifact of §9.

**ZH.** 融合**等位基因**（分型 SNP）与表达证据；对低纯度样本、以及检测拷贝数中性 LOH 最有优势。
- **恶性判定：** `p_cnv` 高 / 归入非整倍体克隆（非正常克隆），**且**具上皮身份。
- 等位基因证据使 Numbat 成为表达类工具分歧时的仲裁者，且较少受 §9 干扰素/MHC 伪影影响。

---

## 8. Consensus labeling & target scoring / 共识标签与靶点评分

**EN. Consensus rule (template — adjust and record):**
1. Restrict to epithelial cells (marker + not in reference groups).
2. Collect three labels: `infercnv_label`, `copykat_pred`, `numbat_call`.
3. **Malignant** if aneuploid/malignant in **≥ 2 of 3** tools. **Normal epithelium** if diploid/normal in **≥ 2 of 3**. Otherwise **ambiguous** (excluded from primary scoring; reported separately).
4. Weight Numbat as tie-breaker for low-purity or IFN-high samples.
5. Emit a per-cell table: `barcode, epithelial, infercnv, copykat, numbat, consensus, purity_flag`.

**Target scoring (only after consensus):**
- Use ambient-corrected, log-normalized expression (§3).
- Report for **TACSTD2** and **CLDN4**, stratified by `{malignant, normal epithelium}`:
  - fraction of positive cells (`{{detection cutoff}}`),
  - mean/median normalized expression,
  - **pseudobulk** per sample (sum raw counts over malignant cells, then normalize) for cross-sample comparison.
- Always show the **normal-epithelium internal control** alongside malignant scores.

**ZH. 共识规则（模板——可调整并记录）：**
1. 限定为上皮细胞（标记阳性且不属于参考组）。
2. 汇总三个标签：`infercnv_label`、`copykat_pred`、`numbat_call`。
3. **恶性**：三工具中 **≥ 2** 判为非整倍体/恶性。**正常上皮**：**≥ 2** 判为二倍体/正常。其余为 **ambiguous**（不进入主评分，单独报告）。
4. 低纯度或 IFN 高的样本，以 Numbat 作仲裁并加权。
5. 输出每细胞表：`barcode, epithelial, infercnv, copykat, numbat, consensus, purity_flag`。

**靶点评分（仅在共识后）：**
- 使用环境校正后的对数归一化表达（§3）。
- 对 **TACSTD2** 与 **CLDN4**，按 `{恶性, 正常上皮}` 分层报告：
  - 阳性细胞比例（`{{检测阈值}}`）、
  - 平均/中位归一化表达、
  - 每样本 **pseudobulk**（对恶性细胞求原始计数之和后归一化），用于跨样本比较。
- 恶性评分旁**始终并列正常上皮内部对照**。

---

## 9. Pitfalls in ICI on-treatment biopsies / ICI 治疗中活检的陷阱

**EN.** On-treatment biopsies under immune-checkpoint inhibition are the hardest case. Watch for:

1. **Low tumor cellularity / high immune infiltrate.** Responding lesions may retain few malignant cells; CNV baselines become unstable and CopyKAT auto-diploid can mis-anchor. → Anchor baseline on §4 references; lean on Numbat allele evidence; flag samples below `{{min malignant cells}}`.
2. **Interferon / MHC artifact.** ICI induces strong ISG and MHC-I/II programs. Coordinated expression at the MHC locus (chr6p) and other ISG clusters can masquerade as focal CNV in expression-smoothing tools (inferCNV/CopyKAT). → Cross-check chr6p and ISG-dense regions against Numbat's allele signal; regress or mask known ISG programs when defining CNV score.
3. **Transcriptional plasticity / EMT.** Treatment can downregulate `EPCAM`/keratins, so malignant cells may fail epithelial gating and be lost or mislabeled stromal. → Use CNV/aneuploidy identity, not epithelial markers alone, to recover malignant cells; keep an "aneuploid but marker-low" bucket.
4. **Stressed / apoptotic / low-quality cells.** Post-treatment necrosis raises mito% and ambient RNA; naively relaxing QC lets debris distort CNV. → Keep QC thresholds documented; do not blanket-relax mito% (§3).
5. **Tumor–immune doublets.** Produce chimeric CNV signals and false "intermediate" cells. → Aggressive doublet removal before CNV.
6. **Antigen modulation of the targets themselves.** TACSTD2/CLDN4 expression may shift on treatment (down-modulation vs subclone selection). Longitudinal decreases must **not** be read as pure "target loss" without confirming the malignant denominator is stable. → Report target scores together with malignant-cell counts, purity, and clone structure per timepoint.
7. **Reference contamination by reactive/regenerating epithelium.** Treatment-injured normal epithelium carries stress signatures that can mimic weak CNV. → Never use epithelium as reference; verify the "normal epithelium" cluster co-embeds with reference and lacks arm-level events.
8. **Batch / timepoint confounding.** Merging baseline and on-treatment cells before CNV creates baseline artifacts. → Run CNV per sample/timepoint; integrate only for visualization, never for the CNV call.
9. **Ambient target leakage.** Lysed tumor cells release TACSTD2/CLDN4 transcripts that contaminate immune cells, inflating apparent expression outside epithelium. → Mandatory ambient correction (§3) and strict gating on consensus-malignant cells.

**ZH.** ICI 治疗中活检是最难的情形，需警惕：

1. **肿瘤细胞少 / 免疫浸润高。** 应答病灶恶性细胞可能很少，CNV 基线不稳，CopyKAT 自动二倍体易锚错。→ 用 §4 参考锚定基线；依赖 Numbat 等位基因证据；对低于 `{{最少恶性细胞数}}` 的样本作标记。
2. **干扰素 / MHC 伪影。** ICI 强烈诱导 ISG 与 MHC-I/II 程序。MHC 位点（chr6p）及其他 ISG 簇的协同表达，在表达平滑类工具（inferCNV/CopyKAT）中可伪装成局灶 CNV。→ 用 Numbat 等位基因信号核对 chr6p 与 ISG 密集区；定义 CNV 评分时回归或屏蔽已知 ISG 程序。
3. **转录可塑性 / EMT。** 治疗可下调 `EPCAM`/角蛋白，使恶性细胞漏出上皮门控或被误判为间质。→ 用 CNV/非整倍体身份而非单靠上皮标记来回收恶性细胞；保留"非整倍体但标记低"分组。
4. **应激 / 凋亡 / 低质量细胞。** 治疗后坏死抬高线粒体比例与环境 RNA；随意放宽 QC 会让碎片扭曲 CNV。→ 记录 QC 阈值；不要一刀切放宽线粒体阈值（§3）。
5. **肿瘤–免疫双细胞。** 产生嵌合 CNV 信号与假"中间态"细胞。→ CNV 前严格去双细胞。
6. **靶点本身的抗原调变。** 治疗中 TACSTD2/CLDN4 表达可能改变（下调 vs 亚克隆选择）。纵向下降**不可**在未确认恶性分母稳定时直接解读为"靶点丢失"。→ 报告靶点评分时并列每时间点的恶性细胞数、纯度与克隆结构。
7. **反应性/再生上皮污染参考。** 治疗损伤的正常上皮带应激特征，可能像弱 CNV。→ 绝不以上皮作参考；确认"正常上皮"簇与参考共嵌且无臂级事件。
8. **批次 / 时间点混杂。** 在 CNV 前合并基线与治疗中细胞会产生基线伪影。→ 按样本/时间点分别运行 CNV；整合仅用于可视化，绝不用于 CNV 判定。
9. **环境靶点泄漏。** 裂解的肿瘤细胞释放 TACSTD2/CLDN4 转录本污染免疫细胞，抬高上皮外表达。→ 强制环境校正（§3），并严格门控在共识恶性细胞上。

---

## 10. Validation / QC of the malignant call / 恶性判定的验证与质控

**EN.**
- Concordance matrix across the three tools (Cohen's κ per pair); record disagreement rate.
- Confirm the consensus **normal-epithelium** cluster co-embeds with the diploid reference and shows no arm-level CNV.
- Compare inferred arm-level events to known driver arms for the indication and, if available, to bulk WGS/WES/SNP-array.
- Report per-sample tumor purity (fraction consensus-malignant among epithelium) and the ambiguous fraction.
- Visualize: CNV heatmaps, UMAP colored by consensus label and by `p_cnv`/CNV score.

**ZH.**
- 三工具一致性矩阵（两两 Cohen's κ）；记录不一致率。
- 确认共识**正常上皮**簇与二倍体参考共嵌且无臂级 CNV。
- 将推断的臂级事件与该癌种已知驱动臂比较；若有，则与 bulk WGS/WES/SNP-array 比较。
- 报告每样本肿瘤纯度（上皮中共识恶性占比）与 ambiguous 比例。
- 可视化：CNV 热图、按共识标签与 `p_cnv`/CNV 评分着色的 UMAP。

---

## 11. Deliverables & reporting checklist / 交付物与报告清单

**EN.**
- [ ] Per-cell consensus table (§8.5).
- [ ] Per-sample tool parameters (§0 + §5–7) and QC thresholds (§3).
- [ ] CNV heatmaps + concordance matrix (§10).
- [ ] Tumor purity + ambiguous fraction per sample.
- [ ] TACSTD2/CLDN4 scores stratified by malignant vs normal epithelium, single-cell and pseudobulk (§8).
- [ ] ICI pitfalls addressed, per sample, with the mitigations from §9 documented.
- [ ] Decision log for every deviation from this template.

**ZH.**
- [ ] 每细胞共识表（§8.5）。
- [ ] 每样本工具参数（§0 + §5–7）与 QC 阈值（§3）。
- [ ] CNV 热图 + 一致性矩阵（§10）。
- [ ] 每样本肿瘤纯度 + ambiguous 比例。
- [ ] 按恶性 vs 正常上皮分层的 TACSTD2/CLDN4 评分，含单细胞与 pseudobulk（§8）。
- [ ] 逐样本处理 §9 的 ICI 陷阱并记录对应缓解措施。
- [ ] 每一处偏离本模板的决策记录。

---

## 12. References / 参考

- Tickle T, et al. **inferCNV** of the Trinity CTAT Project. Broad Institute. https://github.com/broadinstitute/inferCNV
- Gao R, et al. Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes. **CopyKAT**, *Nat Biotechnol* 2021. https://github.com/navinlabcode/copykat
- Gao T, et al. Haplotype-aware analysis of somatic copy number variations from single-cell transcriptomes. **Numbat**, *Nat Biotechnol* 2023. https://github.com/kharchenkolab/numbat
- Huang X, et al. **cellsnp-lite**: efficient genotyping of single cells. *Bioinformatics* 2021.

<!-- END OF TEMPLATE / 模板结束 -->
