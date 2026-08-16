# Bulk RNA immune deconvolution and exclusion scores in NSCLC ICI
# 非小细胞肺癌免疫治疗队列的批量转录组免疫反卷积与排斥评分

**2024–2026 methods playbook** · outputs live only under `methods/bulk_immune/` ·
real statistics only (GSE126044, GSE135222, TCGA LUAD+LUSC, fetched 2026-08-16).

How to read / 读法: each section is English first, then 中文. Numbers in
§8 are copied from `results/demo/KEY_STATS.json`; if you re-run the demo and
they move, the JSON wins.

---

## 0. One-page verdict / 一页结论

**English.** TACSTD2 and CLDN4 are *not* interchangeable immune-exclusion
markers, and a 16- or 27-sample ICI cohort cannot decide the question.

- In **TCGA NSCLC** (n = 1017 primary tumours; 515 LUAD + 502 LUSC) the two
  genes correlate only modestly (Spearman r = 0.251, p = 4.0 × 10⁻¹⁶). Within
  histology the correlation is stronger (LUAD r = 0.460; LUSC r = 0.410) —
  the drop in the pooled number is a mean shift: LUSC has higher TACSTD2
  (median 13.11 vs 12.63 log2 RSEM) and lower CLDN4 (12.07 vs 13.15).
- **TACSTD2** is consistently anti-correlated with the CD8 / cytotoxic axis
  in both histologies (xCell CD8+ T cells: pooled r = −0.229, BH-FDR =
  2.7 × 10⁻¹²; LUAD r = −0.166, FDR = 1.8 × 10⁻³; LUSC r = −0.219, FDR =
  5.9 × 10⁻⁶). That is the most reproducible “exclusion-like” signal we
  measured.
- **CLDN4 vs TIDE Exclusion** is *not* that signal. The pooled r = −0.316
  (FDR = 7.9 × 10⁻²⁴) collapses once LUAD and LUSC are split: LUAD r =
  +0.040, p = 0.37; LUSC r = −0.189, FDR = 1.4 × 10⁻⁴. Pooling histologies
  here is a batch error, not a biological finding.
- In **GSE126044** (n = 16 pre-anti-PD-1 NSCLC; 5 responders / 11
  non-responders) no TACSTD2 or CLDN4 correlation with any exclusion or
  inflamed score survives family-wise BH-FDR. The nominal CLDN4–Exclusion
  r = +0.679 (p = 0.0038, FDR = 0.064) even has the *opposite sign* of the
  pooled TCGA number. Inflamed scores themselves *do* separate responders
  (CYT Cliff’s δ = 1.0, FDR = 0.0028); TACSTD2 does not (δ = −0.27, p =
  0.44).
- In **GSE135222** (n = 27, 21 PFS events) the same pattern: no
  family-wise FDR < 0.05 for TACSTD2/CLDN4 vs exclusion/inflamed; no
  survival association survives FDR.

**中文。** TACSTD2 与 CLDN4 **不能**当成同一个“免疫排斥”标记；16 或 27 例
的 ICI 队列也撑不起这个结论。

- TCGA NSCLC（n = 1017；LUAD 515 + LUSC 502）中两基因 Spearman r 仅 0.251。
  分组织学后升至 LUAD 0.460、LUSC 0.410——合并后变弱是均值平移：LUSC 的
  TACSTD2 更高、CLDN4 更低。
- TACSTD2 在两种组织学中都与 CD8/杀伤轴负相关（xCell CD8：合并 r = −0.229，
  FDR = 2.7×10⁻¹²）。这是本演示里最可重复的“偏排斥”信号。
- CLDN4 与 TIDE Exclusion 的合并 r = −0.316 **不是**同一信号：LUAD 内
  r = +0.040（p = 0.37），LUSC 内 r = −0.189。合并组织学在这里是批次错误。
- GSE126044（n = 16）中 TACSTD2/CLDN4 与任何排斥/炎症评分的相关都过不了
  家族内 BH-FDR。CLDN4–Exclusion 名义 r = +0.679，符号还与合并 TCGA 相反。
  炎症评分本身能分开缓解者（CYT Cliff’s δ = 1，FDR = 0.0028），TACSTD2 不能。
- GSE135222（n = 27）同样：目标基因相关与生存均无家族内 FDR < 0.05。

---

## 1. What is being estimated / 我们在估计什么

**English.** A bulk RNA-seq mixture is one number per gene. Every method
below compresses that mixture into a *score*. Three different things get
called “immune infiltration” and they are not interchangeable:

| Class | Examples | Units | Comparable across cell types? |
|---|---|---|---|
| Marker mean | MCP-counter | log2 expression | **No** |
| Rank enrichment | ssGSEA, ESTIMATE, xCell raw, TIP | arbitrary / cohort-relative | only after the same background |
| Constrained regression | CIBERSORT relative, quanTIseq-style | fraction of the *explained* mix | yes, inside one sample |
| Trained model | TIDE | correlation with a published residual | cohort-relative |

CIBERSORT *relative* mode always sums to 1. An immune-desert tumour still
gets a full set of “fractions”; they are fractions of the immune slice, not
of the tumour. quanTIseq-style NNLS with an “Other” compartment is the
estimator that is allowed to say “almost no immune cells”. TIDE Exclusion
is a correlation with a CAF/MDSC/M2 residual, not a cell count.

**中文。** 每一种方法都把混合表达压成一个分数，但分数的含义不同：MCP-counter
是 marker 的 log2 均值，跨细胞类型不可比；ssGSEA/ESTIMATE/xCell/TIP 是相对
队列背景的秩富集；CIBERSORT 相对模式强制归一化到 1，免疫荒漠样本仍会得到
“组分”；quanTIseq 风格的带 Other 的 NNLS 才允许“几乎没有免疫细胞”；TIDE
Exclusion 是与 CAF/MDSC/M2 残差的相关，不是细胞计数。

---

## 2. Method cards (2024–2026 practice) / 方法卡片

### 2.1 CIBERSORTx / LM22

Newman et al., *Nat Methods* 2015 (CIBERSORT); *Nat Biotechnol* 2019
(CIBERSORTx). LM22 is 547 genes × 22 leukocyte types, Affymetrix-derived.

**Do.** Download LM22 from https://cibersortx.stanford.edu after the
Stanford academic licence. Run *relative* mode on **linear** TPM (RNA-seq)
or RMA (array), `QN=FALSE` for RNA-seq and `QN=TRUE` for arrays. For a
paper that says “CIBERSORTx”, use the official Docker
(`scripts/R/run_cibersortx.sh`) so B-mode / S-mode / absolute mode are the
authors’ implementation.

**Do not.** Log-transform the mixture. ComBat the mixture and then
deconvolve. Report CIBERSORTx numbers from this repo’s nu-SVR. Treat
relative fractions as tumour-level abundance.

**This repo.** `bulkimmune.cibersort.cibersort` is classic nu-SVR
(ν ∈ {0.25, 0.5, 0.75}, clip negatives, rescale to 1, optional permutation
p). LM22 is **not** shipped. Without it the demo uses the GPL-licensed
quanTIseq TIL10 matrix (`constrained_ls_deconvolve` + mRNA scaling +
`Other`) and labels the block `quantiseq_style`, not `cibersort`.

**中文。** LM22 需从 CIBERSORTx 门户自行下载，本仓库不分发。RNA-seq 用线性
TPM、关闭分位数标准化；芯片用 RMA、打开 QN。论文里写 “CIBERSORTx” 必须跑
官方容器。本仓库的 nu-SVR 只复现 2015 年相对模式核心。演示在无 LM22 时改用
quanTIseq TIL10，并如实标注。

### 2.2 xCell

Aran, Hu & Butte, *Genome Biol* 2017. 489 signatures → 64 types, then
calibration (`fv`) and spillover compensation (`K`, α = 0.5).

**Do.** Restrict `cell_types_use` to lineages you expect in a lung tumour
(the demo uses 18). Use the RNA-seq spillover matrix on RNA-seq and the
array matrix on arrays. Report the gene-universe overlap (xCell refuses
below 5,000 of its 10,808 genes).

**Do not.** Compare xCell scores across cohorts that were scored
separately. Steps 4 and 6 subtract a *per-dataset minimum*, so the score
is defined only inside the matrix you handed it. Do not run all 64 types
on a biopsy and interpret osteoblasts.

**This repo.** `bulkimmune.xcell.xcell_analysis` against the official
`xCell.data` resources extracted by `scripts/00_fetch_resources.py`.

**中文。** xCell 分数在数据集内部相对，分开跑的两个队列不可比。必须用对应
平台的 spillover 矩阵，并限制细胞类型。本实现使用官方 489 个签名与 K/fv。

### 2.3 MCP-counter

Becht et al., *Genome Biol* 2016. Ten populations; score = mean log2 of
the official markers.

**Do.** Feed **log2** expression. Report `n_markers_found`. The official
CD8 signature is **one gene (CD8B)** — a 1-gene mean is still the
published estimator; do not silently drop it because `min_markers=2`.

**Do not.** Read the table as fractions. Compare fibroblast scores to CD8
scores. Use linear TPM (the highest-expressed marker then dominates).

**中文。** 必须用 log2；官方 CD8 签名只有 CD8B，不能因为“至少 2 个 marker”
就丢掉。跨细胞类型不可比。

### 2.4 ESTIMATE

Yoshihara et al., *Nat Commun* 2013. Stromal + Immune ssGSEA on a 10,412
gene background; ESTIMATEScore = sum; purity = cos(0.605 + 1.47×10⁻⁴ ×
score) **only on Affymetrix**.

**Do.** Filter to the common-gene background first (the enrichment walk
depends on it). On RNA-seq report Stromal / Immune / ESTIMATEScore and
treat purity as an uncalibrated ordinal proxy, or skip it.

**This repo.** Official SI gene sets from estimate 1.0.13 (R-Forge, GPL-2).
The demo RNA-seq runs do **not** emit TumorPurity.

**中文。** 先限制到 10,412 个 common genes。RNA-seq 不要报告校准过的纯度。

### 2.5 TIDE

Jiang et al., *Nat Med* 2018. Official package: `tidepy` (MIT), offline
weights. For NSCLC the signature SD is the mean of LUAD and LUSC.

**Do.** Pass linear expression; let tidepy log2(x+1) and **row-centre
across the samples in hand**. Use `cancer="NSCLC"`. Leave `pretreat=False`
on treatment-naive biopsies so CTL-high samples are scored on Dysfunction
rather than Exclusion.

**Do not.** Concatenate two cohorts, score once, then compare Exclusion
distributions. The row-centre makes TIDE cohort-relative. Do not treat
`Responder` (TIDE < 0) as a calibrated probability.

**中文。** 用官方 tidepy，`cancer=NSCLC`，初治样本 `pretreat=False`。行中心化
使分数依赖同批样本，不能跨队列直接比 Exclusion。

### 2.6 TIP

Xu et al., *Cancer Res* 2018. Seven-step cancer-immunity cycle. This repo
scores the published annotation table
(`signature annotation.txt` from biocc.hrbmu.edu.cn/TIP) as
ssGSEA(positive) − ssGSEA(negative), and expands step 4 by cell type.

**Do.** Report per-step overlap. Step 5 (infiltration) is the exclusion
family member; step 7 (killing) is the inflamed family member.

**Do not.** Call this “the TIP web-server score” unless you used the
server. Cite Xu 2018 and the annotation file version.

**中文。** 按发表的注释表做 ssGSEA(正)−ssGSEA(负)，不是网页服务器的私有权重。
Step5 归入排斥家族，Step7 归入炎症家族。

### 2.7 ssGSEA IFN / TLS / exclusion signatures

Curated lists (verbatim from the papers) live in
`resources/curated_signatures.json`. Hallmark sets are downloaded from
**MSigDB 2025.1.Hs** — do not paste 50-gene lists by hand.

| Signature | n | Source | Note |
|---|---|---|---|
| IFNG_Ayers6 | 6 | Ayers 2017 *JCI* | unweighted mean in the paper; we also emit ssGSEA |
| TcellInflamed_GEP18 / TIS | 18 | Ayers 2017; Danaher 2018 | clinical assay is a *weighted* NanoString sum; RNA-seq ssGSEA is a proxy, the −1.54 cut-off does **not** transfer |
| ExpandedImmune_Ayers18 | 18 | Ayers 2017 | **different 18 genes** from GEP18; the literature confuses them |
| TLS_12chemokine_Coppola | 12 | Coppola 2011 | overlaps IFN (CCL5, CXCL9/10) — a 12-CK/IFN correlation is partly definitional |
| TLS_Cabrita9 | 9 | Cabrita 2020 *Nature* | **EIF1AY is Y-linked**. GSE135222 is 22 M / 5 F. Use `TLS_Cabrita9_noY` in sex-imbalanced cohorts |
| CYT | 2 | Rooney 2015 | geometric mean of GZMA & PRF1 in **linear TPM**, not ssGSEA |
| Hallmark IFNG / EMT / TGF-β / hypoxia / angiogenesis | 50-ish | Liberzon 2015, MSigDB 2025.1.Hs | exclusion-family EMT/TGF-β/hypoxia; inflamed-family IFNG/IFNA |

**中文。** 小签名按原文基因表固化；Hallmark 从 MSigDB 2025.1.Hs 下载。GEP18
与 ExpandedImmune18 不是同一组 18 个基因。Cabrita 9 基因含 Y 染色体
EIF1AY，性别不平衡队列请用 `TLS_Cabrita9_noY`。CYT 是线性 TPM 上 GZMA/PRF1
的几何平均，不是 ssGSEA。12-CK 与 IFN 签名基因重叠，相关有一部分是定义造成的。

---

## 3. Input-scale contract / 输入尺度

| Method | Required matrix | Log? |
|---|---|---|
| CIBERSORT / quanTIseq-style | TPM / RMA, columns comparable | linear |
| MCP-counter | library-size normalised | **log2** |
| ESTIMATE, xCell, ssGSEA, TIP | any within-sample monotone scale | either (rank-based) |
| TIDE | linear in; tidepy logs and row-centres | log2 after |
| CYT | TPM | linear |

Rank-based methods are invariant to a *within-sample* monotone transform
and **not** invariant to between-sample normalisation. Library size still
matters. `bulkimmune.preprocess.detect_scale` records what the matrix
looked like; the CIBERSORT and MCP-counter entry points refuse the wrong
scale instead of returning a plausible-looking number.

Identifier failures are silent in every published tool: a score is still
emitted on the genes that happened to match. `genes.match_report` is
written next to every score block. GSE135222 identifiers are versioned
Ensembl (`ENSG….10`); they are stripped and mapped through HGNC. Symbol
aliases (e.g. ESTIMATE’s `FYB` → `FYB1`) are resolved the same way.

**中文。** 尺度错了分数仍会看起来合理。CIBERSORT/MCP-counter 在入口处拒绝错误
尺度。所有工具在基因 ID 对不上时仍会输出分数——必须同时报告签名覆盖率。
GSE135222 是带版本号的 Ensembl，经 HGNC 映射到符号。

---

## 4. How to correlate with TACSTD2 / CLDN4 / 如何与 TACSTD2、CLDN4 相关

**English.**

1. Score TACSTD2 and CLDN4 on the **same** log matrix as the rank-based
   immune scores (they are epithelial genes; do not z-score them against
   immune genes).
2. Primary test: Spearman, not Pearson. Fractions are zero-inflated;
   ssGSEA is skewed.
3. Secondary test: partial Spearman after residualising both ranks on
   (batch + ESTIMATEScore). This asks “does the association survive
   purity and batch”, not “is it causal”.
4. Do **not** add TACSTD2 and CLDN4 into a 2-gene “exclusion signature”
   until you have shown they agree. In these data they do not (see §8).
5. When an ICI endpoint exists, test the *immune scores* against the
   endpoint as a positive control (does this cohort even have an immune–
   response relationship?) before testing the epithelial genes.

**中文。** 主检验用 Spearman。偏相关只在秩上扣除批次与 ESTIMATEScore。在证明
TACSTD2 与 CLDN4 方向一致之前，不要把它们加成“排斥签名”——本数据里它们并不
一致。有 ICI 终点时，先用免疫评分对终点做阳性对照，再测上皮基因。

---

## 5. Batch, RNA-seq vs microarray / 批次与平台

**English.**

- **Correct the question, not the matrix.** Residualise the *scores* on
  batch + purity (`bulkimmune.batch.residualise`). ComBat on a CIBERSORT
  mixture breaks the linear mixing assumption. Rank-based scores computed
  on a ComBat-merged matrix are only valid if those samples will be
  analysed together afterwards.
- **RNA-seq and microarray are a method problem, not a ComBat problem.**
  LM22 is Affymetrix; xCell ships two spillover matrices; ESTIMATE purity
  is Affymetrix-only. Run the platform-appropriate configuration. Do not
  force the matrices onto one scale and compare scores.
- **Histology is a batch.** LUAD vs LUSC shifts TACSTD2 and CLDN4 in
  opposite directions. A pooled correlation can be large and FDR-significant
  and still be a mean shift (CLDN4 vs TIDE Exclusion, §8). Always stratify.
- **FFPE vs fresh is a batch.** GSE126044 is 11 fresh + 5 FFPE. Of 103
  scores, **zero** had BH-FDR < 0.05 against preservation, so we did not
  ComBat; we residualised. That is the correct order: test, then decide.
- **xCell and TIDE are cohort-relative.** Never concatenate, score, split.

**中文。** 优先在分数上对批次/纯度做残差，而不是 ComBat 表达矩阵再反卷积。
RNA-seq 与芯片要换方法配置，不是靠 ComBat 对齐。组织学（LUAD/LUSC）本身就是
批次：CLDN4–Exclusion 的合并显著相关在分层后消失。FFPE/新鲜：GSE126044 的
103 个评分对保存方式的 BH-FDR 全 > 0.05，因此不做 ComBat，只做残差。xCell
与 TIDE 是队列内相对分数，禁止合并后再拆开比。

---

## 6. Multiple testing / 多重检验

**English.** BH-FDR is applied **inside (target × family)**, not on the
whole table. Three families were specified *before* looking at p-values:

| Family | Members (short names) | Hypothesis |
|---|---|---|
| `exclusion` | TIDE, Exclusion, CAF, MDSC, TAM M2, fibroblasts, endothelium, StromalScore, TIP step 5, Hallmark EMT / TGF-β / hypoxia / angiogenesis | TACSTD2/CLDN4 mark an excluded epithelium |
| `inflamed` | IFNG/CD8/CTL/CYT, Ayers/GEP/TLS sets, ImmuneScore, MCP/xCell T/B/NK/CD8, TIP steps 4 & 7, Hallmark IFNG/IFNA/inflammatory/allograft | the competing “hot tumour” hypothesis |
| `other` | everything else | not a pre-specified claim |

A hit in `other` that does not survive its own FDR is a hypothesis for the
next cohort, not a finding. Comparing “Exclusion FDR = 0.01 vs B-cell FDR
= 0.04, therefore Exclusion wins” is not licensed by this procedure.

With n = 16 or 27, family-wise FDR < 0.05 on a Spearman correlation is
almost unavailable (r needs to be ≳ 0.75). Report the nominal r, the FDR,
and the TCGA direction, and stop. Do not “relax to p < 0.05”.

**中文。** BH-FDR 在（目标基因 × 家族）内做，不在整张表上做。`other` 家族里
未过 FDR 的结果只是下一队列的假说。n = 16/27 时 Spearman 几乎不可能过家族
FDR（大约需要 \|r\| ≳ 0.75）。应同时报告名义 r、FDR 和 TCGA 方向，而不是
放宽到 p < 0.05。

---

## 7. Runnable pipeline / 可运行流程

```bash
cd methods/bulk_immune
pip install -r requirements.txt          # numpy pandas scipy scikit-learn statsmodels tidepy
# Rscript is optional; only needed to extract xCell S4 signatures
bash scripts/05_run_demo.sh
```

| Script | Role |
|---|---|
| `scripts/00_fetch_resources.py` | MCP-counter, ESTIMATE 1.0.13, xCell.data, TIP annotation, TIL10, MSigDB 2025.1.Hs Hallmark, HGNC |
| `scripts/01_prepare_geo.py` | GSE126044 counts + response/FFPE; GSE135222 FPKM + PFS; ID mapping |
| `scripts/02_prepare_tcga.py` | Xena HiSeqV2 LUAD+LUSC, primary tumours only (`-01`) |
| `scripts/03_score_cohort.py` | every method, skip-and-record if a resource is missing |
| `scripts/04_correlate.py` | Spearman, partial, batch KW, two-group, survival; family-wise BH |
| `scripts/06_make_figures.py` | optional bar plots of the two pre-specified families |
| `scripts/07_concordance.py` | cross-method Spearman on CD8/IFN/exclusion/TLS axes; cross-cohort sign table |
| `scripts/R/run_cibersortx.sh` | official CIBERSORTx Docker (licence + token required) |

Python 3.12, no R required for the GEO/TCGA demo once xCell resources have
been extracted once. Re-run `00_fetch_resources.py` on a machine with
`Rscript` if `resources/xCell_signatures.gmt` is missing.

**中文。** 一条 `05_run_demo.sh` 跑完资源下载、GEO/TCGA 准备、评分与相关。
xCell 签名的首次解包需要 R；之后纯 Python。CIBERSORTx 官方镜像单独走
`run_cibersortx.sh`。

---

## 8. Demo results — real statistics only / 演示结果（只报实测）

Fetched 2026-08-16. Full tables: `results/demo/`. Figures:
`results/demo/*_highlight.png`.

### 8.1 Cohorts / 队列

| Cohort | Assay | n | Endpoint | Batch covariate | TACSTD2 & CLDN4 |
|---|---|---|---|---|---|
| GSE126044 (Cho 2020) | RNA-seq counts → log2 CPM | 16 | ICI response 5 / 11 | 11 fresh / 5 FFPE | both present |
| GSE135222 (Jung 2019) | RNA-seq (columns ≈ 1e6) | 27 | PFS, 21 events | sex 22 M / 5 F | both present |
| TCGA LUAD+LUSC (Xena HiSeqV2) | log2 RSEM | 1017 (515+502) | none (not ICI) | histology | both present |

TCGA is treatment-naive surgical NSCLC. It estimates correlation
*direction* with enough n to survive FDR. It is not ICI-response evidence.

**中文。** TCGA 是初治手术标本，只用来看相关方向，不能写成 ICI 疗效证据。

### 8.2 TACSTD2 vs CLDN4 / 两基因彼此

| Cohort | n | Spearman r | p |
|---|---|---|---|
| GSE126044 | 16 | 0.488 | 0.055 |
| GSE135222 | 27 | 0.350 | 0.073 |
| TCGA pooled | 1017 | 0.251 | 4.0 × 10⁻¹⁶ |
| TCGA LUAD | 515 | 0.460 | 2.7 × 10⁻²⁸ |
| TCGA LUSC | 502 | 0.410 | 8.9 × 10⁻²² |

They co-vary, weakly. They are not a single axis.

### 8.3 TCGA: pre-specified families / TCGA 预设家族

Pooled n = 1017, BH-FDR inside (target × family).

| Target | Score | Family | r | p | FDR |
|---|---|---|---|---|---|
| TACSTD2 | xCell CD8+ T cells | inflamed | −0.229 | 1.5 × 10⁻¹³ | 2.7 × 10⁻¹² |
| TACSTD2 | MCP-counter CD8 T cells | inflamed | −0.164 | 1.4 × 10⁻⁷ | 8.5 × 10⁻⁷ |
| CLDN4 | TIDE Exclusion | exclusion | −0.316 | 4.7 × 10⁻²⁵ | 7.9 × 10⁻²⁴ |
| CLDN4 | TIDE MDSC | exclusion | −0.314 | 1.2 × 10⁻²⁴ | 1.0 × 10⁻²³ |

**Histology split of the two headline numbers:**

| Contrast | LUAD r (FDR) | LUSC r (FDR) |
|---|---|---|
| TACSTD2 vs xCell CD8 | −0.166 (1.8 × 10⁻³) | −0.219 (5.9 × 10⁻⁶) |
| CLDN4 vs TIDE Exclusion | +0.040 (p = 0.37) | −0.189 (1.4 × 10⁻⁴) |

TACSTD2 vs CD8 **replicates in both histologies**. CLDN4 vs Exclusion
**does not**; the pooled star is a LUAD/LUSC mean shift (LUSC: higher
TACSTD2, lower CLDN4, and a different Exclusion distribution). After
residualising on cohort + ESTIMATEScore, TACSTD2 vs TIDE Exclusion is
r = −0.034, p = 0.28.

Partial Spearman (cohort + ESTIMATEScore) keeps a TACSTD2 association
with TIP step 5 (partial r = 0.262, FDR = 3.2 × 10⁻¹⁶) and a weak
negative with TIDE (partial r = −0.080, FDR = 0.026). Interpret step 5
cautiously: it is a small published gene set, not a cell count.

**中文。** TACSTD2–CD8 在 LUAD 与 LUSC 中方向一致；CLDN4–Exclusion 的合并
显著相关在 LUAD 内消失（r = 0.04），是组织学均值平移。对 cohort +
ESTIMATEScore 偏相关后，TACSTD2 与 TIDE Exclusion 不再相关（r = −0.034）。

### 8.4 GSE126044 — ICI response, n = 16 / 免疫治疗缓解

**Target vs scores.** 0 / (3 targets × 2 families) tests at FDR < 0.05.
Nominal (not findings):

| Target | Score | r | p | FDR |
|---|---|---|---|---|
| CLDN4 | TIDE Exclusion | +0.679 | 0.0038 | 0.064 |
| CLDN4 | MCP-counter NK | −0.688 | 0.0032 | 0.097 |
| TACSTD2 | TIDE Exclusion | +0.488 | 0.055 | 0.74 |

The CLDN4–Exclusion sign is the opposite of pooled TCGA and of LUSC. With
n = 16 that is what noise plus a 5-sample FFPE subset looks like.

**Endpoint positive control.** Inflamed scores separate the 5 responders
from the 11 non-responders; the epithelial genes do not.

| Score | Cliff’s δ (R vs NR) | p | FDR |
|---|---|---|---|
| TIDE CTL.flag | +1.00 | 1.4 × 10⁻⁴ | 0.0028 |
| CYT | +1.00 | 4.6 × 10⁻⁴ | 0.0028 |
| MCP-counter CD8 | +1.00 | 4.6 × 10⁻⁴ | 0.0028 |
| TcellInflamed_GEP18 | +0.96 | 9.2 × 10⁻⁴ | 0.0041 |
| TACSTD2 | −0.27 | 0.44 | 0.59 |
| CLDN4 | −0.49 | 0.11 | 0.26 |

Median log2 CPM: TACSTD2 6.17 (R) vs 6.29 (NR); CLDN4 2.54 vs 3.77.

**Batch.** 0 / 103 scores at BH-FDR < 0.05 for fresh vs FFPE. No ComBat.

**中文。** n = 16 时目标基因相关全部不过家族 FDR；CLDN4–Exclusion 名义正相关
还与 TCGA 相反。炎症评分能完美分开 5 例缓解者（CYT/CD8 Cliff’s δ = 1），
TACSTD2/CLDN4 不能。fresh/FFPE 对 103 个评分无 FDR 显著批次效应。

### 8.5 GSE135222 — PFS, n = 27 / 无进展生存

TACSTD2 vs CLDN4 r = 0.350, p = 0.073. Family-wise FDR < 0.05 for
TACSTD2/CLDN4 vs exclusion/inflamed: **none** (one `other`-family
correlation for TACSTD2 only). StromalScore vs TACSTD2+CLDN4 mean:
r = −0.547, p = 0.0032, exclusion-family FDR = 0.054 — report as
nominal. Median-split log-rank on 21 events: no score survives family
FDR (top nominal: TIP step-4 macrophage p = 0.0066, FDR = 0.24;
Hallmark TGF-β p = 0.021, FDR = 0.18). Cabrita TLS was computed with
and without EIF1AY because the cohort is 22 male / 5 female. The Coppola
12-CK set is **8/12 genes** here: CCL3, CCL4, CCL5 and CCL18 are absent
from the published FPKM matrix after Ensembl→symbol mapping. Do not quote
that TLS score as the published 12-gene signature.

**中文。** 目标基因与预设家族相关、以及生存，全部不过家族 FDR。Cabrita TLS
因 22 男 / 5 女同时报告去 Y 版本。本队列 12-CK 只覆盖 8/12（缺 CCL3/4/5/18）。

### 8.6 What the demo did *not* run / 演示未跑的部分

LM22 / official CIBERSORTx (no Stanford token in this environment).
quanTIseq-style TIL10 *was* run and is in `block_quantiseq_style.tsv`.
ESTIMATE TumorPurity was not emitted (RNA-seq). xCell used the 18
biopsy-plausible types, not all 64.

### 8.7 Cross-method concordance / 方法间一致性

Tables: `results/demo/method_concordance.tsv`, `CONCORDANCE.json`.
Median pairwise Spearman among the pre-specified axes:

| Axis | GSE126044 median r (pairs ≥0.5 / all) | GSE135222 | TCGA |
|---|---|---|---|
| CD8 / cytotoxic | 0.80 (42/45) | 0.76 (43/45) | 0.70 (36/45) |
| IFN / inflamed | 0.84 (15/15) | 0.81 (15/15) | 0.87 (15/15) |
| TLS / B cell | 0.76 (13/15) | 0.67 (13/15) | 0.79 (15/15) |
| Exclusion / stroma | 0.34 (10/28) | 0.55 (16/28) | 0.60 (17/28) |

**CD8.** MCP-counter CD8B, TIDE CD8 and TIDE CTL are almost the same
ranking in TCGA (MCP vs TIDE CD8 r = 0.960). xCell CD8 vs MCP CD8 is
0.682 — same direction, not interchangeable. quanTIseq-style TIL10
`T.cells.CD8` agrees with nobody (r = 0.07–0.20). That block is a
licence-free smoke test, not a CD8 readout; do not put it in a TACSTD2
figure next to MCP/xCell/TIDE.

**IFN.** Ayers 6-gene vs TIDE IFNG r = 0.982 in TCGA. GEP18 vs
ExpandedImmune18 r = 0.971 (different gene lists, same axis). Hallmark
IFNG vs ImmuneScore r = 0.869.

**Exclusion is not one axis.** TIDE Exclusion vs MDSC r = 0.805, vs CAF
r = 0.515, vs MCP fibroblasts r = 0.311, vs ESTIMATE StromalScore
r = **−0.170**. MDSC vs StromalScore r = −0.572. A “TACSTD2 vs
exclusion family” result that is carried by MDSC is not a fibroblast
result, and vice versa. Report the member, not the family name.

**Sign concordance of TACSTD2/CLDN4 vs headline scores**
(`results/demo/sign_concordance.tsv`) — three cohorts, sign only:

| Pair | 126044 r | 135222 r | TCGA r | signs agree? |
|---|---|---|---|---|
| TACSTD2 vs xCell CD8 | −0.26 | −0.25 | −0.229 | **yes, all −** |
| TACSTD2 vs ImmuneScore | −0.24 | −0.22 | −0.080 | **yes, all −** |
| TACSTD2 vs TIDE Exclusion | +0.49 | +0.07 | +0.099 | **yes, all +** |
| CLDN4 vs xCell CD8 | −0.31 | −0.30 | −0.052 | **yes, all −** |
| TACSTD2 vs MCP CD8 | −0.26 | +0.09 | −0.164 | no |
| CLDN4 vs TIDE Exclusion | +0.68 | −0.02 | −0.316 | no |

The only TACSTD2–immune statement that is sign-stable across the ICI
n=16, the ICI n=27 and TCGA n=1017 is: **higher TACSTD2, lower CD8 /
ImmuneScore, slightly higher TIDE Exclusion**. Magnitude in the ICI
cohorts is not distinguishable from noise after FDR. CLDN4 vs Exclusion
is the pair that flips.

**中文。** CD8/IFN/TLS 轴在方法间高度一致（TCGA 中位 r 0.70–0.87）；quanTIseq
风格的 TIL10 CD8 与谁都不一致（r 0.07–0.20），不能当 CD8 读数。Exclusion
家族不是一条轴：Exclusion–MDSC r = 0.81，Exclusion–StromalScore r = −0.17。
三队列符号一致的只有：TACSTD2↑ 伴随 xCell CD8↓ / ImmuneScore↓ / TIDE
Exclusion↑；CLDN4–Exclusion 符号不一致。

### 8.8 Signature coverage / 签名覆盖率

`results/demo/signature_overlap_all.tsv`. Curated IFN/GEP/TLS sets are
complete in GSE126044 and TCGA. Exceptions that change interpretation:

- GSE135222 Coppola 12-CK: **8/12** (CCL3, CCL4, CCL5, CCL18 missing from
  the published matrix).
- MCP-counter CD8 is officially 1 gene (CD8B) in all three cohorts.
- MCP myeloid DC misses WFDC21P (5/6) everywhere; T-cell set misses
  CHRM3-AS2 / MGC40069 (14/16).
- ESTIMATE common-gene background: 9821/10412 (GSE126044), 9781/10412
  (GSE135222), **10412/10412** (TCGA Xena symbols).

**中文。** GSE135222 的 12-CK 只有 8/12。MCP CD8 官方就是单基因 CD8B。
ESTIMATE 背景在 TCGA 上齐，在两个 GEO 上缺约 600 个 common genes。

---

## 9. What we do not claim / 不声称的东西

- That TACSTD2 or CLDN4 **predicts ICI response**. In the only ICI cohort
  with a response label, they do not, and the immune scores do.
- That the two genes are a joint exclusion signature. They correlate
  modestly and associate with different immune axes; CLDN4’s pooled
  Exclusion hit is a histology shift.
- That any ssGSEA / xCell / TIDE number is a cell fraction.
- That this Python CIBERSORT core is CIBERSORTx.
- That TCGA surgical expression is an ICI dataset.

**中文。** 不声称 TACSTD2/CLDN4 预测 ICI 疗效（有缓解标签的队列里它们不能，
免疫评分能）；不声称二者构成联合排斥签名；不把 ssGSEA/xCell/TIDE 当成细胞
比例；不把本仓库的 nu-SVR 写成 CIBERSORTx；不把 TCGA 手术标本当成 ICI 数据。

---

## 10. Citations / 文献

- Newman AM et al. *Nat Methods* 2015;12:453–457. CIBERSORT.
- Newman AM et al. *Nat Biotechnol* 2019;37:773–782. CIBERSORTx.
- Aran D, Hu Z, Butte AJ. *Genome Biol* 2017;18:220. xCell.
- Becht E et al. *Genome Biol* 2016;17:218. MCP-counter.
- Yoshihara K et al. *Nat Commun* 2013;4:2612. ESTIMATE.
- Jiang P et al. *Nat Med* 2018;24:1550–1558. TIDE.
- Xu L et al. *Cancer Res* 2018;78:6575–6584. TIP.
- Barbie DA et al. *Nature* 2009;462:108–112. ssGSEA.
- Hänzelmann S, Castelo R, Guinney J. *BMC Bioinformatics* 2013;14:7. GSVA.
- Ayers M et al. *J Clin Invest* 2017;127:2930–2940. IFN-γ / GEP.
- Danaher P et al. *J Immunother Cancer* 2018;6:63. TIS.
- Rooney MS et al. *Cell* 2015;160:48–61. CYT.
- Coppola D et al. *Am J Pathol* 2011;179:37–45. 12-CK TLS.
- Cabrita R et al. *Nature* 2020;577:561–565. 9-gene TLS.
- Liberzon A et al. *Cell Syst* 2015;1:417–425. MSigDB Hallmarks.
- Finotello F et al. *Genome Med* 2019;11:34. quanTIseq / TIL10.
- Johnson WE, Li C, Rabinovic A. *Biostatistics* 2007;8:118–127. ComBat.
- Cho JW et al. GEO GSE126044; related literature on anti-PD-1 NSCLC RNA-seq.
- Jung H et al. *Nat Commun* 2019;10:3847. GSE135222.
- Cancer Genome Atlas Research Network. LUAD 2014 / LUSC 2012; Xena HiSeqV2.

---

## 11. Re-run / 复现

```bash
cd methods/bulk_immune
PYTHONPATH=. python3 tests/test_ssgsea.py
bash scripts/05_run_demo.sh
python3 -c "import json; print(json.load(open('results/demo/KEY_STATS.json'))['TCGA_NSCLC']['highlighted'])"
```

If GEO or Xena is unreachable the committed `results/demo/*.tsv` still
hold the numbers in §8.
