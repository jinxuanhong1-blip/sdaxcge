<!--
  methods/perturbation/playbook.md
  Bilingual (English + 中文) methods playbook for knockdown / knockout RNA-seq.
  Every GEO accession named here was verified live via NCBI E-utilities.
  No accession in this document is invented.
-->

# Knockdown / Knockout RNA-seq Analysis Playbook<br>敲低 / 敲除 RNA-seq 分析手册

**Scope / 适用范围:** siRNA / shRNA / CRISPRi / CRISPR knockout perturbation
RNA-seq, differential expression with **DESeq2** and **edgeR**, **GSEA**
(Hallmark interferon, EMT, tight junction), directional "opposite-gene" checks
(e.g. *CLDN4* knockdown → *TACSTD2/TROP2*?), batch effects, the pitfalls of
**n = 2–3**, and **mouse vs human** symbol harmonization. Includes copy-paste
templates and a reproducible way to mine GEO for *CLDN4/TACSTD2* perturbations.

> **Golden rules / 基本原则**
> 1. **Never invent an accession.** Query GEO/SRA live; paste only IDs the query returned. / **绝不编造 accession**，只粘贴查询真实返回的 ID。
> 2. **DE tests need RAW integer counts**, not FPKM/TPM/RPKM. / 差异分析必须用**原始整数 counts**，不能用 FPKM/TPM/RPKM。
> 3. **Model batch in the design; don't "clean" counts before testing.** / 把批次放进设计矩阵建模，**不要**在检验前"清洗"counts。
> 4. **A sanity check first:** does the perturbed gene actually go down? / **先做 sanity check**：被敲的基因真的下调了吗？

---

## Table of contents / 目录
1. [Experimental design & controls / 实验设计与对照](#1-design)
2. [From reads to a count matrix / 从测序到 counts 矩阵](#2-quant)
3. [QC before DE / 差异分析前的质控](#3-qc)
4. [Differential expression: DESeq2 & edgeR / 差异表达](#4-de)
5. [GSEA: Hallmark IFN, EMT, tight junction / 基因集富集](#5-gsea)
6. [The "opposite gene" check (CLDN4→TACSTD2) / "反向基因"检验](#6-opposite)
7. [Batch effects / 批次效应](#7-batch)
8. [The n = 2–3 pitfalls / 小样本量陷阱](#8-nsmall)
9. [Mouse vs human symbols / 小鼠与人类基因符号](#9-symbols)
10. [Mining GEO for CLDN4/TACSTD2 perturbations / 在 GEO 中检索](#10-geo)
11. [Templates & worked example / 模板与实例](#11-templates)

---

<a name="1-design"></a>
## 1. Experimental design & controls / 实验设计与对照

**EN.** The single most important design decision is the **control**, because
"differential expression" is always *relative to it*:

| Perturbation | Correct control | Why |
|---|---|---|
| siRNA | **non-targeting / scrambled siRNA (NTC)** at the same total siRNA dose | separates gene-specific effects from transfection + RNAi machinery load |
| shRNA (lenti) | **shControl / shScramble / shLuc / shGFP**, matched MOI + selection | controls for viral integration, puromycin selection, hairpin load |
| CRISPRi (dCas9-KRAB) | **non-targeting sgRNA (NTG)** in the same dCas9 line | controls for dCas9-KRAB tethering |
| CRISPR-KO (Cas9) | **safe-harbor sgRNA** (e.g. targeting *AAVS1/ROSA26*) or NTG; ideally a **rescue** arm | a cut without a gene target controls for DNA-damage response |

Other design rules:
- **Replicates are biological**, not technical: independent transfections /
  infections / clones. Two lanes of one sample are *not* n = 2.
- **Independent guides/siRNAs** (≥2 per gene) guard against off-target effects —
  a hit reproducible across independent guides is far more credible.
- **Randomize/balance** condition across batches (see §7). Never put all
  controls in batch 1 and all knockdowns in batch 2.
- **Time & efficiency matter:** an siRNA at 48 h and a stable KO clone are
  different biology (acute vs adapted). Record and report knockdown efficiency
  (qPCR/WB), and expect the RNA readout of the target itself to correlate.

**中文.** 最关键的设计是**对照**，因为"差异表达"永远是**相对对照而言**的：siRNA 用
**非靶向/乱序 siRNA (NTC)** 且总剂量相同；shRNA 用 **shControl/shScramble/shLuc/shGFP**
且 MOI 与筛选一致；CRISPRi 用**非靶向 sgRNA**；CRISPR 敲除用**安全港位点 sgRNA**
（如 *AAVS1/ROSA26*）或非靶向 sgRNA，最好再加一个**回补(rescue)**组。其他规则：
**重复必须是生物学重复**（独立转染/感染/克隆，同一样本两条 lane 不算 n=2）；**每个基因
≥2 条独立 guide/siRNA** 以排除脱靶；**跨批次随机化/平衡**条件（见 §7）；记录并报告敲低
效率（qPCR/WB），并预期靶基因自身在 RNA 层面同向变化。

---

<a name="2-quant"></a>
## 2. From reads to a count matrix / 从测序到 counts 矩阵

**EN.** Two mainstream routes; either yields the **gene-level raw counts** that
DESeq2/edgeR require:

- **Alignment-free (recommended for speed):** `salmon` or `kallisto` →
  transcript quant → `tximport` (R) collapses to gene-level counts +
  `countsFromAbundance="lengthScaledTPM"` for DE.
- **Alignment-based:** `STAR` → BAM → `featureCounts` (Subread) or `htseq-count`
  with the matching GTF.

```bash
# salmon (per sample) -- index built from a transcriptome + genome decoy
salmon quant -i salmon_index -l A \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --validateMappings --gcBias -p 8 -o quant/sample
# then in R: tximport(files, type="salmon", tx2gene=tx2gene)  -> gene counts
```

Keep the **exact genome build + annotation version** (e.g. GENCODE vXX) in your
methods — gene symbols and IDs change between releases.

**中文.** 两条主流路线，都产出 DESeq2/edgeR 所需的**基因水平原始 counts**：
(1) **免比对**（推荐、快）：`salmon`/`kallisto` 定量转录本 → R 里用 `tximport` 汇总到基因，
DE 用 `countsFromAbundance="lengthScaledTPM"`；(2) **比对**：`STAR` → BAM →
`featureCounts` 或 `htseq-count`（配套 GTF）。务必在方法学里写清**基因组版本+注释版本**
（如 GENCODE vXX），因为不同版本符号/ID 会变。

---

<a name="3-qc"></a>
## 3. QC before DE / 差异分析前的质控

**EN.** Do these *before* trusting any p-value:
- **Read QC:** `fastqc` + `multiqc`; check adapter content, duplication, rRNA %.
- **Alignment/assignment rate:** low % assigned to genes → annotation mismatch.
- **Library-size normalization + variance stabilization** for plotting: DESeq2
  `vst()`/`rlog()` or edgeR log-CPM.
- **Sample-relationship plots:** PCA and a sample-distance heatmap on the VST
  matrix. You are looking for (a) replicates clustering, (b) the biggest axis of
  variation being condition — or, if it is **batch**, that tells you to model it
  (§7). An **outlier** replicate shows up here.
- **On-target check (do not skip):** confirm the perturbed gene's counts drop in
  the KD/KO arm. If *CLDN4* is not down in a *CLDN4* knockdown, stop and debug
  the sample sheet / labels before interpreting anything else.

**中文.** 在相信任何 p 值之前先做：**读段质控**（`fastqc`+`multiqc`：接头、重复率、rRNA%）；
**比对/计数分配率**（偏低多为注释不匹配）；用 DESeq2 `vst()/rlog()` 或 edgeR log-CPM 做
**方差稳定**后画 **PCA** 与**样本距离热图**（看重复是否聚类、最大变异轴是条件还是**批次**——
若是批次则需建模，见 §7；离群重复也在此暴露）；以及**靶基因验证（不可跳过）**：确认被敲基因
在 KD/KO 组 counts 下降。若 *CLDN4* 敲低实验里 *CLDN4* 没下调，先排查样本表/标签再谈其他。

---

<a name="4-de"></a>
## 4. Differential expression: DESeq2 & edgeR / 差异表达

**EN.** Both fit a **negative-binomial** model to raw counts; for the small
replicate numbers typical of perturbation work they behave similarly. Choices
that matter more than "DESeq2 vs edgeR":

- **Reference level = control**, explicitly (`relevel`/`ref_level`). Otherwise R
  picks alphabetically and your sign flips.
- **Pre-filter** low-count genes (DESeq2: ≥10 counts in ≥ smallest-group size;
  edgeR: `filterByExpr`). Fewer tests → more power after FDR.
- **Shrink log2 fold-changes** for ranking/plotting/GSEA (DESeq2 `lfcShrink`
  type `apeglm`; edgeR reports `logFC` from the QL fit). Shrinkage tames the
  wild LFCs of low-count genes that dominate naïve rankings.
- **edgeR:** use the **quasi-likelihood** F-test (`glmQLFit`+`glmQLFTest`), which
  controls type-I error better than the LRT at small n.
- **Report** log2FC, an FDR (BH-adjusted p, "padj"/"FDR"), and base mean/CPM.
  Threshold thoughtfully (e.g. padj < 0.05 **and** |log2FC| > 1), and remember
  padj = NA (independent filtering / all-zero) is not "significant".

Templates: [`templates/deseq2_template.R`](templates/deseq2_template.R),
[`templates/edger_template.R`](templates/edger_template.R), and a pure-Python
DESeq2 via [`templates/pydeseq2_template.py`](templates/pydeseq2_template.py).

```r
# DESeq2 essentials (full script in the template)
dds <- DESeqDataSetFromMatrix(cts, coldata, design = ~ batch + condition)
dds$condition <- relevel(dds$condition, ref = "control")
dds <- dds[rowSums(counts(dds) >= 10) >= min(table(dds$condition)), ]
dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "condition_knockdown_vs_control", type = "apeglm")
```

**中文.** 两者都对原始 counts 拟合**负二项**模型；在扰动实验常见的小重复量下表现相近。
比"DESeq2 还是 edgeR"更重要的选择：**显式把参照设为对照**（否则 R 按字母序选，符号会反）；
**预过滤**低表达基因（DESeq2：≥10 counts 且出现在≥最小组样本数；edgeR：`filterByExpr`）；
**收缩 log2FC** 用于排序/画图/GSEA（DESeq2 `lfcShrink(type="apeglm")`；edgeR 用 QL 的
`logFC`）；**edgeR 用准似然 F 检验**（`glmQLFit`+`glmQLFTest`），小 n 下 I 型错误控制更好；
**报告** log2FC、FDR（BH 校正的 padj/FDR）与 baseMean/CPM，阈值要合理（如 padj<0.05 **且**
|log2FC|>1），且 padj=NA 不等于"显著"。模板见上方链接。

---

<a name="5-gsea"></a>
## 5. GSEA: Hallmark IFN, EMT, tight junction / 基因集富集

**EN.** For perturbation data, **pre-ranked GSEA** is the workhorse: rank *all*
genes by a DE metric and ask whether a gene set concentrates at the top/bottom.
This avoids arbitrary cutoffs and works even when few genes pass FDR.

- **Ranking metric:** shrunken **log2FC** (clean, monotone) or **signed
  −log10(p)** (`sign(logFC) * -log10(PValue)`, rewards significance). Pick one
  and state it.
- **Gene sets to prioritize here** (all real MSigDB/KEGG names):
  - `HALLMARK_INTERFERON_ALPHA_RESPONSE`
  - `HALLMARK_INTERFERON_GAMMA_RESPONSE`
  - `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`
  - Tight junction: `KEGG_TIGHT_JUNCTION` (MSigDB C2:CP:KEGG) or GO:CC
    `GOCC_TIGHT_JUNCTION` / KEGG `hsa04530`.
  Why these: claudins (incl. *CLDN4*) are core **tight-junction** components;
  tight-junction loss is coupled to **EMT**; and interferon programs frequently
  move when epithelial adhesion/identity is perturbed.
- **Read the sign correctly:** positive NES = enriched among genes **up in the
  knockdown/knockout**. Always report NES + FDR (q), not just "enriched".
- **Tools:** R `fgsea`/`clusterProfiler::GSEA` + `msigdbr`; Python `gseapy.prerank`.

Templates: [`templates/gsea_fgsea_clusterprofiler.R`](templates/gsea_fgsea_clusterprofiler.R),
[`templates/gseapy_prerank_template.py`](templates/gseapy_prerank_template.py).
See the worked example (§11) for a real IFN-α/γ result on *CLDN4* KO.

**中文.** 扰动数据首选**预排序 GSEA**：用 DE 指标给**所有**基因排序，检验某基因集是否富集
在头/尾，避免人为阈值，即使很少基因过 FDR 也能用。**排序指标**：收缩后的 **log2FC**，或
**带符号 −log10(p)**（`sign(logFC)*-log10(p)`），二选一并写明。**本主题优先的基因集**（均为
真实 MSigDB/KEGG 名称）：`HALLMARK_INTERFERON_ALPHA_RESPONSE`、
`HALLMARK_INTERFERON_GAMMA_RESPONSE`、`HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`、
紧密连接 `KEGG_TIGHT_JUNCTION`（或 GO:CC `GOCC_TIGHT_JUNCTION` / KEGG `hsa04530`）。
原因：claudin（含 *CLDN4*）是紧密连接核心，其缺失与 **EMT** 耦联，干扰素程序也常随上皮
黏附/身份被扰动而变化。**符号判读**：NES 为正 = 富集于**在 KD/KO 中上调**的基因；务必报告
NES+FDR(q)。工具：R `fgsea`/`clusterProfiler::GSEA`+`msigdbr`，Python `gseapy.prerank`。

---

<a name="6-opposite"></a>
## 6. The "opposite gene" check (CLDN4 → TACSTD2/TROP2) / "反向基因"检验

**EN.** A recurring, legitimate question: *when I knock down gene A, does a
functionally linked gene B move in a predictable direction?* Here the concrete
case is **knock down `CLDN4` → does `TACSTD2` (TROP2) go up?** The biology behind
the hypothesis: both are apical epithelial surface molecules enriched at cell–cell
junctions, and TROP2 has documented cross-talk with claudins and tight-junction
integrity — so one might expect a **compensatory up-regulation** of *TACSTD2* when
*CLDN4* is lost. **This is a hypothesis, not a fact — you must test it.**

How to test it rigorously (turn "it went up" into a defensible statement):
1. Pull the DE row for the responder: **log2FC**, **padj**.
2. Place it in context: its **rank / percentile** in the genome-wide LFC
   distribution ("top 3% up-regulated" ≫ "log2FC = +0.4, padj = 0.3").
3. Check **direction + significance + reproducibility** (both independent guides,
   both cell lines if applicable).
4. State the verdict plainly: **supported / wrong direction / not significant.**
5. Also check the **paralog/related set** (`EPCAM`, other claudins `CLDN1/3/7`,
   `TJP1/OCLN`) so you don't cherry-pick one gene.

Reusable tool: [`templates/opposite_gene_check.py`](templates/opposite_gene_check.py).

```bash
python templates/opposite_gene_check.py de_results.tsv \
    --target CLDN4 --opposite TACSTD2 --expect up
```

**Reality check from the worked example (GSE207704, CLDN4 KO):** *TACSTD2* went
**down** together with *CLDN4* (mean log2FC −0.81), i.e. the compensatory-up
hypothesis was **not supported** in that dataset. That negative result is exactly
why the check exists — report it honestly rather than assuming the expected sign.

**中文.** 一个常见且合理的问题：**敲低基因 A 时，功能相关的基因 B 是否会朝可预期的方向变化？**
本例的具体问题是 **敲低 `CLDN4` → `TACSTD2`(TROP2) 会上调吗？** 假设的生物学依据：二者都是顶端
上皮表面分子、富集于细胞连接处，且 TROP2 与 claudin/紧密连接完整性存在已知 cross-talk，因此
可能预期 *CLDN4* 缺失时 *TACSTD2* 出现**代偿性上调**。**这只是假设，必须检验。** 严谨做法：
(1) 取该基因的 **log2FC** 与 **padj**；(2) 放进全基因组 LFC 分布看**排名/百分位**（"上调前 3%"
远强于"log2FC=+0.4, padj=0.3"）；(3) 检查**方向+显著性+可重复性**（两条独立 guide、多个细胞系）；
(4) 明确给出结论：**支持 / 方向相反 / 不显著**；(5) 同时看**旁系/相关集**（`EPCAM`、其他 claudin
`CLDN1/3/7`、`TJP1/OCLN`），避免只挑一个基因。工具见上。**实例现实结果（GSE207704，CLDN4 敲除）：**
*TACSTD2* 与 *CLDN4* **一起下调**（平均 log2FC −0.81），代偿性上调假设**不成立**——这正是要做
该检验的原因：如实报告，而非默认符合预期。

---

<a name="7-batch"></a>
## 7. Batch effects / 批次效应

**EN.**
- **Detect:** colour the PCA by batch (sequencing run, transfection day, clone,
  operator, kit lot). If a non-condition axis separates samples, you have a batch.
- **The correct fix is to MODEL it**, not to erase it: `~ batch + condition` in
  DESeq2/edgeR. The test then estimates the condition effect *within* batch.
- **`limma::removeBatchEffect`** produces a corrected matrix **for visualization
  only** (PCA/heatmap). **Never** feed it into DESeq2/edgeR — you'd be testing
  doctored data and understating variance.
- **`ComBat_seq` (sva)** adjusts raw counts and is appropriate mainly when you
  must **merge datasets** quantified separately and cannot otherwise model batch.
- **Unknown/latent batch:** estimate surrogate variables with `sva`, or use
  `RUVSeq` with control genes, and add them to the design.
- **The unfixable case — confounding:** if every control is batch 1 and every
  knockdown is batch 2, batch and condition are **aliased** and *no* method can
  separate them. The only fix is a balanced re-run. **Prevent this by design.**

Template: [`templates/batch_correction.R`](templates/batch_correction.R)
(includes a confounding check that stops early).

**中文.** **发现**：用批次（测序批、转染日、克隆、操作者、试剂批号）给 PCA 上色，若非条件轴
分开样本即存在批次。**正确做法是建模而非抹除**：DESeq2/edgeR 用 `~ batch + condition`，在
批次**内部**估计条件效应。**`limma::removeBatchEffect`** 的矫正矩阵**仅用于可视化**（PCA/热图），
**绝不能**送进 DESeq2/edgeR。**`ComBat_seq`(sva)** 调整原始 counts，主要用于**合并**分别定量的
数据集且无法在设计中建模批次时。**未知/潜在批次**：用 `sva` 估计代理变量，或用 `RUVSeq` 配合
对照基因，加入设计。**无法挽救的情形——混杂**：若对照全在批次1、敲低全在批次2，批次与条件
**完全共线**，任何方法都无法分离，唯一解是平衡地重做——**靠设计预防**。

---

<a name="8-nsmall"></a>
## 8. The n = 2–3 pitfalls / 小样本量陷阱

**EN.** Perturbation studies often run at **n = 2–3 per group**. What that means:
- **n = 2 is the floor** for DESeq2/edgeR (you cannot estimate within-group
  variance at n = 1). n = 3 is markedly more robust; prefer it if you can.
- **Dispersion estimates are shaky** at small n. Both tools **borrow strength
  across genes** (empirical-Bayes shrinkage), which is exactly why you should use
  DESeq2/edgeR and **not** a naïve per-gene t-test or "fold-change only".
- **A single outlier dominates.** DESeq2 flags/soft-replaces outliers via Cook's
  distance (`refit_cooks`), but with n = 2 there is no majority to appeal to —
  inspect PCA/counts for the driver gene and don't let one library set the story.
- **Ranking beats hard thresholds.** With few significant genes, lean on
  **pre-ranked GSEA** (§5) and effect-size ranking rather than counting how many
  genes cleared padj < 0.05.
- **Independent guides are worth more than more replicates of one guide** for
  ruling out off-target artefacts.
- **FPKM/TPM without replicates → no valid p-value.** If a processed matrix has
  one column per condition (common in GEO supplements), you can rank and run
  pre-ranked GSEA, but you **cannot** claim DESeq2/edgeR significance. Say so.
  (This is exactly the situation in the worked example — see §11.)
- **Power is limited:** report effect sizes and treat borderline hits as
  hypothesis-generating; validate the key ones (qPCR/independent cohort).

**中文.** 扰动实验常见 **每组 n=2–3**：**n=2 是 DESeq2/edgeR 的下限**（n=1 无法估计组内方差），
**n=3 明显更稳**，能则优先。小 n 下**离散度估计不稳**，两者都**跨基因借力**（经验贝叶斯收缩），
这正是要用 DESeq2/edgeR 而**非**朴素 t 检验或"只看 fold-change"的原因。**单个离群样本影响巨大**：
DESeq2 用 Cook 距离标记/替换离群（`refit_cooks`），但 n=2 时没有"多数"可依，需查 PCA/该基因
counts，别让单个文库主导结论。**排序优于硬阈值**：显著基因少时，多依赖**预排序 GSEA**（§5）与
效应量排序，而非数有多少基因过 padj<0.05。为排除脱靶，**多一条独立 guide 比同一 guide 多一个
重复更值钱**。**无重复的 FPKM/TPM → 无有效 p 值**：若处理后矩阵每条件仅一列（GEO 补充文件常见），
可排序并跑预排序 GSEA，但**不能**声称 DESeq2/edgeR 显著性——要写清楚（实例正是此情形，见 §11）。
**统计效力有限**：报告效应量，边缘结果作假设生成，关键结论另行验证（qPCR/独立队列）。

---

<a name="9-symbols"></a>
## 9. Mouse vs human symbols / 小鼠与人类基因符号

**EN.**
- **Casing convention:** human symbols are all-caps (`CLDN4`, `TACSTD2`, `EPCAM`);
  mouse are title-case (`Cldn4`, `Tacstd2`, `Epcam`). A quick `toupper()` is fine
  for eyeballing a handful of focus genes but **fails silently** for many genes
  (ORFs, riken clones, 1-to-many orthologs) — do not use it for tables that feed
  statistics.
- **Gene sets are human-centric:** MSigDB Hallmark and most KEGG sets use human
  symbols. To run GSEA on a **mouse** experiment, map mouse → human orthologs
  first (or use `msigdbr(species="Mus musculus")`, which returns mouse-symbol
  versions of the same sets).
- **Use a real ortholog table:** `babelgene` (offline, reproducible) or `biomaRt`
  `getLDS` (queries Ensembl, always current). Keep 1:many mappings explicit.
- **Harmonize aliases/deprecations before merging studies:** e.g. TROP2 → official
  `TACSTD2`; older claudin names → current `CLDN*`. Map to a fixed annotation
  release and record it.

Template: [`templates/ortholog_map_mouse_human.R`](templates/ortholog_map_mouse_human.R).

**中文.** **大小写规则**：人类符号全大写（`CLDN4`、`TACSTD2`、`EPCAM`），小鼠首字母大写
（`Cldn4`、`Tacstd2`、`Epcam`）。`toupper()` 只适合肉眼核对少数关注基因，对很多基因（ORF、
RIKEN 克隆、1对多同源）会**静默出错**，不能用于喂给统计的表格。**基因集以人类为主**：MSigDB
Hallmark 与多数 KEGG 用人类符号；跑**小鼠**实验的 GSEA 前先做小鼠→人类同源映射，或用
`msigdbr(species="Mus musculus")` 得到同一套集的小鼠符号版本。**使用真实同源表**：`babelgene`
（离线、可复现）或 `biomaRt::getLDS`（查 Ensembl、最新）；显式保留 1对多。**合并研究前先统一
别名/废弃名**：如 TROP2→官方 `TACSTD2`，旧 claudin 名→现行 `CLDN*`；固定到某注释版本并记录。

---

<a name="10-geo"></a>
## 10. Mining GEO for CLDN4/TACSTD2 perturbations / 在 GEO 中检索

**EN.** Use NCBI E-utilities (`db=gds`) or EDirect so you only ever cite **real**
accessions. Two provided tools do this and also list a series' supplementary
files with sizes, so you can pick a processed matrix under your download budget:

```bash
# stdlib-only Python (no installs); lists suppl files + sizes with --suppl
python templates/geo_query.py "CLDN4 knockout" "CLDN4 shRNA" "TACSTD2 knockdown" --suppl

# EDirect equivalent with fielded queries (RNA-seq series only)
bash templates/geo_query_edirect.sh
```

Refine on the GEO website with **fielded queries**, e.g.:

```
CLDN4[Description] AND "expression profiling by high throughput sequencing"[DataSet Type] AND "Homo sapiens"[Organism]
```

**Verified public accessions** (live E-utilities + GEO FTP; re-query before
citing). Full table: [`public_series_catalog.md`](public_series_catalog.md).

| Accession | Organism | Assay | Perturbation | Suppl size | Used |
|---|---|---|---|---:|---|
| **GSE207704** | human | RNA-seq | CLDN4 CRISPR KO vs WT (MCF7, T47D) | 1.0 MB FPKM | §11 |
| **GSE245459** | human | RNA-seq | shTACSTD2 vs shNC ± cisplatin (SKOV3, n=3) | 14 MB FPKM | §11 |
| **GSE334497** | mouse | RNA-seq | Trop2/Tacstd2 KO vs WT (4T1 tumors, n=5) | 1.2 MB norm. counts | §11 |
| **GSE22493** | human | microarray | CLDN4-silencing vs over-expression (SKOV-3) | 13 MB CEL | catalog only |

How to choose a *runnable* set: prefer a **Series with a processed count/FPKM
matrix in `suppl/`** that is small (MBs). If only FASTQs exist, budget for SRA
download + re-quantification (often ≫ 2 GB) — outside a quick <2 GB run.

> **No-fabrication rule:** if a search returns nothing for a gene/assay
> combination, **report "none found"** — do not paste a plausible-looking GSE.
> Dedicated *human cell-line TACSTD2 CRISPR-KO* RNA-seq was **not** found in
> this search (GSE334497 is mouse in-vivo KO; GSE245459 is shRNA). Additional
> human *CLDN4* KD/KO RNA-seq besides GSE207704 was **not** found (GSE22493 is
> array). Those absences are findings.

**中文.** 用 NCBI E-utilities（`db=gds`）或 EDirect，确保只引用**真实** accession。提供的两个
工具会检索并列出某 Series 的补充文件及大小，便于按下载预算挑选处理后矩阵（命令见上）。在 GEO
网站可用**带字段查询**精炼。**已核实的公开 accession**（完整表见
[`public_series_catalog.md`](public_series_catalog.md)）：**GSE207704**（人，CLDN4 CRISPR KO）、
**GSE245459**（人，SKOV3 shTACSTD2，n=3）、**GSE334497**（小鼠，4T1 Trop2 KO，n=5）、
**GSE22493**（人，芯片，CLDN4 沉默，仅编目）。**如何挑可运行的数据集**：优先选 `suppl/` 中带
处理后矩阵、体积小（MB 级）的 Series。
> **禁止编造原则**：查无结果就写"未找到"。本次检索**未找到**独立的人细胞系 TACSTD2 CRISPR-KO
> RNA-seq（GSE334497 是小鼠体内 KO；GSE245459 是 shRNA），也**未找到** GSE207704 以外的人
> CLDN4 KD/KO RNA-seq（GSE22493 是芯片）。缺失本身就是结论。

---

<a name="11-templates"></a>
## 11. Templates & worked example / 模板与实例

**Templates** (`templates/`):

| File | Purpose |
|---|---|
| [`deseq2_template.R`](templates/deseq2_template.R) | DESeq2 DE + LFC shrinkage + ranked file for GSEA |
| [`edger_template.R`](templates/edger_template.R) | edgeR quasi-likelihood DE + ranked file |
| [`pydeseq2_template.py`](templates/pydeseq2_template.py) | DESeq2 in pure Python (PyDESeq2), no R |
| [`gsea_fgsea_clusterprofiler.R`](templates/gsea_fgsea_clusterprofiler.R) | pre-ranked GSEA in R (fgsea + msigdbr), IFN/EMT/tight-junction |
| [`gseapy_prerank_template.py`](templates/gseapy_prerank_template.py) | pre-ranked GSEA in Python (gseapy), human & mouse |
| [`opposite_gene_check.py`](templates/opposite_gene_check.py) | CLDN4→TACSTD2-style directional hypothesis test |
| [`batch_correction.R`](templates/batch_correction.R) | model batch, ComBat-seq, removeBatchEffect + confounding check |
| [`ortholog_map_mouse_human.R`](templates/ortholog_map_mouse_human.R) | mouse↔human ortholog/symbol harmonization |
| [`geo_query.py`](templates/geo_query.py) / [`geo_query_edirect.sh`](templates/geo_query_edirect.sh) | mine GEO for perturbation series + suppl sizes |

**Worked examples** (all public, all < 2 GB, none invent accessions):

| Run | Series | Perturbation | On-target | Opposite gene | Hallmark IFN |
|---|---|---|---|---|---|
| [`example_run/`](example_run/README.md) | GSE207704 | CLDN4 KO (human lines) | CLDN4 −0.88 | TACSTD2 **−0.81** (not UP) | **down** (FDR 0.007 / 0.035) |
| [`example_gse245459/`](example_gse245459/README.md) | GSE245459 | shTACSTD2 (SKOV3, n=3) | TACSTD2 −2.63 | CLDN4 **−1.92** (not UP) | **down** (vehicle arm) |
| [`example_gse334497/`](example_gse334497/README.md) | GSE334497 | Tacstd2 KO (4T1, n=5) | Tacstd2 −3.82 (mean) | Cldn4 **−0.82** (not UP) | **up** (in-vivo mouse) |

Cross-series takeaway: the compensatory-up story (CLDN4↓→TACSTD2↑ or the reverse)
is **not supported** in any of the three public matrices. IFN direction is
**model-dependent** (down in two human cell-line KDs/KOs; up in mouse in-vivo
Trop2 KO) — do not pool. All three supplements are FPKM/normalized, so p-values
are exploratory; the count-model templates remain the path for raw-count DE.

```bash
pip install pandas numpy scipy gseapy
python example_run/run_gse207704.py
python example_gse245459/run_gse245459.py
python example_gse334497/run_gse334497.py
```

**中文.** **模板**见上表。**三个公开实例**（均 <2GB，无编造 accession）见上表：CLDN4 敲除、
shTACSTD2、小鼠 Tacstd2 敲除。跨数据集结论：代偿性上调（CLDN4↓→TACSTD2↑ 或其反向）在三个
公开矩阵中**均不成立**；干扰素方向**依赖模型**（两个人细胞系下调，小鼠体内 KO 上调）——不要合并。
三份补充文件都是 FPKM/标准化值，p 值仅作探索；原始 counts 的正规 DE 仍走模板。

---

### Environment / 环境

Python deps for the templates/example: `pip install pandas numpy scipy gseapy pydeseq2`.
R deps: `DESeq2`, `apeglm`, `edgeR`, `limma`, `sva`, `fgsea`, `clusterProfiler`,
`msigdbr`, `babelgene`/`biomaRt` (Bioconductor). All accessions in this document
were verified live via NCBI E-utilities on the date of writing; re-verify before
citing. / 模板与实例的依赖见上；本文所有 accession 均于撰写时经 NCBI E-utilities 实时核实，
引用前请再次验证。
