# BayesPrism / EcoTyper × TACSTD2 / CLDN4

**Scope: METHODS ONLY — 仅方法学**

> 中文 — 本文件是 `playbook.md` 的**配套比较协议**，只规定 BayesPrism（含 InstaPrism）与 EcoTyper 如何对照 `TACSTD2` 与 `CLDN4`。此处**不包含任何数据、结果、效应量、p 值或结论**。文中数字均为协议参数。数值产出必须写入 `results/` 并回引本节编号。
>
> EN — Companion comparison protocol to `playbook.md`. It specifies how BayesPrism (including InstaPrism) and EcoTyper are to be read **against** `TACSTD2` and `CLDN4`. It contains **no data, no results, no effect sizes, no p-values, and no conclusions**. Every number is a protocol parameter. Numeric output belongs in `results/` and must cite the section IDs here.

| | |
|---|---|
| Path | `methods/deconv_advanced/bayesprism_ecotyper_vs_tacstd2_cldn4.md` |
| Parent | `methods/deconv_advanced/playbook.md` (§1.1, §1.4, §4.1, §4.4, §7, §8) |
| Related | `methods/wgcna/playbook.md` (data-driven junction module); `methods/exclusion_scores/playbook.md`; `methods/causal/playbook.md` |
| Languages | 中文 + English |

---

## 目录 / Contents

- [§B0 为什么这两个方法对这两个基因 / Why these two methods against these two genes](#b0-为什么这两个方法对这两个基因--why-these-two-methods-against-these-two-genes)
- [§B1 基因契约 / Gene contract](#b1-基因契约--gene-contract)
- [§B2 两个方法回答的不同问题 / What each method can estimate](#b2-两个方法回答的不同问题--what-each-method-can-estimate)
- [§B3 估计量 / Estimands](#b3-估计量--estimands)
- [§B4 BayesPrism 臂 / BayesPrism arm](#b4-bayesprism-臂--bayesprism-arm)
- [§B5 EcoTyper 臂 / EcoTyper arm](#b5-ecotyper-臂--ecotyper-arm)
- [§B6 对照协议 / Head-to-head protocol](#b6-对照协议--head-to-head-protocol)
- [§B7 稳定性分级 / Stability grading](#b7-稳定性分级--stability-grading)
- [§B8 预注册冻结清单 / Freeze checklist](#b8-预注册冻结清单--freeze-checklist)
- [§B9 陷阱 / Pitfalls](#b9-陷阱--pitfalls)
- [§B10 参考文献 / References](#b10-参考文献--references)

---

## §B0 为什么这两个方法对这两个基因 / Why these two methods against these two genes

中文 — `playbook.md` §7 把 `TACSTD2` 当作单一上皮暴露、把多种去卷积方法当作可互换的比例估计器。本文件收窄到一个更尖锐的对照：

1. **`TACSTD2` 与 `CLDN4` 不是同一个基因的两个别名。** `TACSTD2`（TROP2）是跨膜糖蛋白，文献中与 claudin-7 有物理互作并参与紧密连接稳定；`CLDN4` 是经典紧密连接 claudin，常与 `TACSTD2` 共表达，但**不是**已被证明的 TROP2 结合伴侣。把 `CLDN4` 当成 `CLDN7` 的替身，或把两者塌缩成一个"TROP2/claudin 分数"再去看结果，都是方法学错误。本协议强制**分开估计、再联合**。
2. **BayesPrism 与 EcoTyper 也不是同一类工具。** BayesPrism 给出连续的细胞比例 $\theta$ 与区室特异表达 $Z$。EcoTyper 给出细胞**状态**（CS）与多细胞**群落**（CE）的丰度。用 EcoTyper 的 CE 丰度去"验证" BayesPrism 的 CD8 比例，或反过来，是范畴错误。本协议把它们当作**两条正交证据线**，并预先规定何种一致/不一致可以写成什么级别的主张。
3. **两者的交点才是本文件的对象：** 恶性/上皮区室里的 `TACSTD2` 与 `CLDN4`（BayesPrism $Z$）是否落在特定的上皮 CS 与特定的 CE（EcoTyper recovery）上；以及这两条线对免疫读数（CD8、TLS、淋巴细胞缺乏型 CE）是否给出同方向的叙述。

EN — `playbook.md` §7 treats `TACSTD2` as a single epithelial exposure and the deconvolution methods as interchangeable fraction estimators. This file narrows to a sharper contrast:

1. **`TACSTD2` and `CLDN4` are not two names for one gene.** `TACSTD2` (TROP2) is a transmembrane glycoprotein with a published physical interaction with claudin-7 and a role in tight-junction stability. `CLDN4` is a canonical tight-junction claudin, frequently co-expressed with `TACSTD2`, but **not** the demonstrated TROP2 binding partner. Treating `CLDN4` as a stand-in for `CLDN7`, or collapsing both into one "TROP2/claudin score" *before* seeing concordance, is a methodological error. This protocol estimates them **separately, then jointly**.
2. **BayesPrism and EcoTyper are not the same class of tool.** BayesPrism yields continuous fractions $\theta$ and compartment-specific expression $Z$. EcoTyper yields cell-**state** (CS) and multicellular-**community** (CE) abundances. Using a CE abundance to "validate" a BayesPrism CD8 fraction, or the reverse, is a category error. This protocol treats them as **two orthogonal evidence lines** and pre-specifies which agreements and disagreements may support which grade of claim.
3. **The object of this file is their intersection:** whether malignant/epithelial-compartment `TACSTD2` and `CLDN4` (BayesPrism $Z$) land on particular epithelial CS and particular CE (EcoTyper recovery), and whether the two lines tell the same directional story about immune readouts (CD8, TLS, lymphocyte-deficient CE).

配套手册分工 / Division of labor with sibling playbooks:

| Question | Lives in |
|---|---|
| How to run BayesPrism / EcoTyper at all | `playbook.md` §4.1, §4.4 |
| How to correlate one epithelial gene with CD8/TLS | `playbook.md` §7 |
| How to *discover* a data-driven TACSTD2–CLDN4 co-expression module | `methods/wgcna/playbook.md` |
| How to score published exclusion signatures | `methods/exclusion_scores/playbook.md` |
| How BayesPrism $Z$ and EcoTyper CS/CE speak to **both** genes, and to each other | **this file** |

---

## §B1 基因契约 / Gene contract

### §B1.1 标识 / Identifiers

| Symbol | Ensembl (unversioned) | UniProt | Role in this protocol | Must not be treated as |
|---|---|---|---|---|
| `TACSTD2` | ENSG00000184292 | P09758 | Seed A; TROP2 | interchangeable with `CLDN4` |
| `CLDN4` | ENSG00000189143 | O14493 | Seed B; canonical TJ claudin | the TROP2 binding partner |
| `CLDN7` | ENSG00000113231 | O95471 | Pre-specified **secondary** (published physical partner of TROP2) | optional / post-hoc add-on |
| `CLDN1` | ENSG00000163347 | O95832 | Pre-specified secondary TJ | a seed |
| `OCLN` | ENSG00000197822 | Q16625 | Pre-specified secondary TJ | a seed |

中文 — 小鼠同源：`Tacstd2`、`Cldn4`、`Cldn7`、`Cldn1`、`Ocln`。跨物种分析见 `methods/mouse_human/playbook.md`，本文件只规定人源 bulk。

EN — Mouse orthologs: `Tacstd2`, `Cldn4`, `Cldn7`, `Cldn1`, `Ocln`. Cross-species mapping lives in `methods/mouse_human/playbook.md`; this file specifies human bulk only.

### §B1.2 锁定的屏障小集合 / Locked barrier mini-set

中文 — 本协议使用**预先锁定的 5 基因集合** `BARRIER5 = {TACSTD2, CLDN4, CLDN7, CLDN1, OCLN}`。这不是 WGCNA 模块，也不是从当前队列里筛出来的。扩大、缩小或替换该集合必须作为事后修改写入 `prereg.md`，并触发 §B7 的降级。数据驱动的模块发现（若要做）走 `methods/wgcna/`，其产出可以作为**预先声明的敏感性轴**投影到本协议，但不能回写 `BARRIER5`。

EN — This protocol uses a **pre-locked 5-gene set** `BARRIER5 = {TACSTD2, CLDN4, CLDN7, CLDN1, OCLN}`. It is not a WGCNA module and is not selected from the cohort under analysis. Expanding, shrinking, or replacing the set is a post-hoc amendment to `prereg.md` and triggers a demotion under §B7. Data-driven module discovery, if performed, lives in `methods/wgcna/`; its output may be projected here as a **pre-declared sensitivity axis** but must not rewrite `BARRIER5`.

```yaml
# methods/deconv_advanced/config/barrier_genes.yaml
set_id: BARRIER5
locked: true
seeds:
  - {symbol: TACSTD2, ensembl: ENSG00000184292, role: seed_A}
  - {symbol: CLDN4,   ensembl: ENSG00000189143, role: seed_B}
secondaries:
  - {symbol: CLDN7, ensembl: ENSG00000113231, role: published_TROP2_partner}
  - {symbol: CLDN1, ensembl: ENSG00000163347, role: tj}
  - {symbol: OCLN,  ensembl: ENSG00000197822, role: tj}
scoring:
  primary: singscore          # sample-self-contained; see playbook §6.3
  sensitivity: [zmean, ssgsea]
  coverage_floor: 0.80        # drop the whole set if retained fraction is below this
do_not:
  - treat_CLDN4_as_CLDN7
  - collapse_seeds_before_concordance
  - add_genes_after_seeing_associations
```

### §B1.3 两个种子的技术差异 / Technical differences between the two seeds

| Property | `TACSTD2` | `CLDN4` | Protocol consequence |
|---|---|---|---|
| Exon structure | Intronless (single coding exon) | Multi-exon | snRNA-seq references can shift `TACSTD2` more than `CLDN4` (`playbook.md` §2.3 point 3). Carry chemistry as a sensitivity axis **separately per gene**. |
| Normal-epithelium expression | Tissue-restricted | Broader in many epithelia | Normal-vs-malignant mis-assignment (a known BayesPrism/TME failure mode at high purity) contaminates `CLDN4` $Z$ more readily. The malignant-fraction floor in §B4.3 is therefore **stricter for `CLDN4` E3** than for `TACSTD2` E3 if a single floor must be chosen; default is to apply the same floor and report both. |
| DNA-contamination risk | Higher (intronless) | Lower | FFPE / high-gDNA libraries: inspect both genes against a DNA-contamination metric before E1. |
| Expected compartment | Epithelial / malignant | Epithelial / malignant | Both E3 analyses are restricted to the malignant (BayesPrism `key`) or Epithelial.cells (EcoTyper) compartment. Immune-compartment $Z$ for either gene is a QC diagnostic, not an exposure. |

中文 — 共表达是**待检验的假设**，不是本协议的前提。两个种子落入不同方向、或只有一个与免疫读数关联，是可报告的方法学结局，不是把它们平均掉的理由。

EN — Co-expression is a **hypothesis to test**, not a premise of this protocol. Opposite directions, or an association carried by only one seed, is a reportable methodological outcome — not a reason to average the two genes.

---

## §B2 两个方法回答的不同问题 / What each method can estimate

| Capability | BayesPrism / InstaPrism | EcoTyper recovery (Carcinoma model) |
|---|---|---|
| Cell-type fractions $\theta$ | yes (type and nested state) | no (CS abundance is not a type fraction) |
| Per-sample, per-type expression $Z$ | **yes** (native posterior) | no; consumes CIBERSORTx HiRes $Z$ only if discovery is run |
| Cell **states** within a type | only if state labels are supplied in the reference | **yes** (NMF states, version-bound) |
| Multicellular **communities** | no | **yes** (CE1–CE10 in the pre-built carcinoma model) |
| Handles `TACSTD2`/`CLDN4` as exposures | via $Z$ in the malignant compartment | via association of bulk or $Z$ with CS/CE abundance |
| Input units | raw counts | log2, then **unit-variance scaling per dataset / tumor type** |
| Minimum $n$ (practical) | no hard floor from the method; association floors follow `playbook.md` §7 | recovery: authors recommend $>25$ samples; below that, CE assignment is exploratory only |
| Supported histologies (recovery) | any, if a matching scRNA-seq reference exists | the 16 carcinoma types listed by the EcoTyper carcinoma portal (BLCA, BRCA, CESC, CHOL, COAD, ESCA, HNSC, LUAD, LUSC, OV, PAAD, PRAD, READ, STAD, THCA, UCEC). **Do not recover the carcinoma model in melanoma, GBM, SCLC, or hematologic disease** without a pre-declared, histology-matched model |
| What it must not be asked | "which ecotype is this sample?" | "how much `TACSTD2` does each malignant cell express?" |

中文 — 一句话：BayesPrism 回答**区室内部的基因剂量**；EcoTyper 回答**该剂量落在哪一种细胞状态与哪一种多细胞群落里**。本协议的主对照是这两句话是否指向同一方向，而不是谁的 CD8 数字更准（那是 `playbook.md` §8 的问题）。

EN — In one sentence: BayesPrism answers **gene dose inside a compartment**; EcoTyper answers **which cell state and which multicellular community that dose sits in**. The primary contrast in this file is whether those two sentences point the same way — not whose CD8 number is more accurate (that question lives in `playbook.md` §8).

### §B2.1 EcoTyper 预建模型中与本问题相关的标签 / Labels in the pre-built carcinoma model that this question may use

中文 — 仅使用 Luca et al. 2021 的 **Carcinoma** 预建模型做 recovery。下列标签是该模型的公开输出，**不是**本协议对新状态的命名。引用时必须绑定模型版本（git commit / 发布标签）。

EN — Recovery uses only the Luca et al. 2021 **Carcinoma** pre-built model. The labels below are that model's public outputs, **not** names this protocol assigns to new states. Bind every citation to a model version (git commit / release tag).

| Object | Public labels | How this protocol may use them |
|---|---|---|
| Cell types | `B.cells`, `PCs`, `CD8.T.cells`, `CD4.T.cells`, `NK.cells`, `Monocytes.and.Macrophages`, `Dendritic.cells`, `Mast.cells`, `PMNs`, `Fibroblasts`, `Endothelial.cells`, `Epithelial.cells` | Map `Epithelial.cells` CS to the barrier genes; map `CD8.T.cells` CS and B/PC CS to the immune side |
| Ecotypes | CE1–CE10 | Pre-declared immune-context axes: **CE1** (lymphocyte-deficient), **CE9** (IFN-γ / ICI-favorable in the discovery papers), **CE2** (least favorable in subsequent ICI recovery), **CE10** (proinflammatory). Other CEs are reported but are not primary endpoints |
| CS identifiers | integer states within each cell type | **no stable semantics across discovery runs**; recovery of the frozen carcinoma model is the only setting in which a CS id may be named in a claim |

中文 — 本协议**不预设**哪一个上皮 CS 高表达 `TACSTD2` 或 `CLDN4`。那是 E8 的检验对象，不是输入。

EN — This protocol **does not pre-assign** which epithelial CS is `TACSTD2`-high or `CLDN4`-high. That assignment is the object of E8, not an input.

---

## §B3 估计量 / Estimands

中文 — `playbook.md` §7.2 的 E1–E3 对**每一个种子基因各做一遍**。本文件追加四个联合估计量。主分析在看到任何关联之前冻结一张表：哪个是主、哪个是共同主、哪个是强制对照。默认如下，修改必须写入 `prereg.md`。

EN — E1–E3 from `playbook.md` §7.2 are run **once per seed gene**. This file adds four joint estimands. Freeze a table — primary, co-primary, mandatory control — before any association is computed. Defaults below; any change goes in `prereg.md`.

| ID | Question | Exposure | Outcome | Tier (default) |
|---|---|---|---|---|
| **E1-T / E1-C** | Does bulk gene dose track immune composition? | bulk log2(TPM+1) of `TACSTD2` or `CLDN4` | $\mathrm{clr}(\theta)_{\text{CD8}}$ ; primary TLS | Mandatory descriptive control (purity-confounded) |
| **E2-T / E2-C** | Same, holding purity fixed | bulk gene + $\pi$ | same | Co-primary |
| **E3-T / E3-C** | Does **per-malignant-cell** gene dose track immune composition? | $E^{\text{mal}}_{g}$ from BayesPrism $Z$ (CPM within malignant) | same | **Primary** (gene-wise) |
| **E4** | Do the two seeds agree inside the malignant compartment? | $E^{\text{mal}}_{\texttt{TACSTD2}}$ vs $E^{\text{mal}}_{\texttt{CLDN4}}$ | concordance (Spearman + Lin's CCC) | Mandatory diagnostic; **not** an immune claim |
| **E5** | Is one seed sufficient? | residual of $E^{\text{mal}}_{\texttt{TACSTD2}} \mid E^{\text{mal}}_{\texttt{CLDN4}}$ and the reverse | $\mathrm{clr}(\theta)_{\text{CD8}}$ ; TLS | Secondary (unique-information test) |
| **E6** | Does the locked barrier set track immune composition? | `BARRIER5` singscore computed on **malignant-compartment $Z$** (not on bulk) | $\mathrm{clr}(\theta)_{\text{CD8}}$ ; TLS | Secondary; only after E4 is reported |
| **E7** | Do the seeds track EcoTyper communities? | E3-T, E3-C, and E6 | CE1 / CE9 / CE2 / CE10 abundance (continuous) | **Co-primary community endpoint** |
| **E8** | Which recovered epithelial CS, if any, carries the seeds? | E3-T, E3-C | epithelial CS abundances from recovery | Exploratory assignment; pre-declared, not mined |

中文 — 规则：

- **R-B1** E4 必须在 E6 之前完成。种子不一致时，E6 降为探索性，不得用 `BARRIER5` 分数掩盖单基因分歧。
- **R-B2** E7 的结局是 **CE 丰度（连续）**，不是最高 CE 的硬分类。硬分类只作描述。
- **R-B3** E1/E2/E3 对两个种子**对称报告**。只报告"显著的那个基因"视为选择性报告。
- **R-B4** 免疫结局仍遵守 `playbook.md` R0-1（成分变换）与 R0-2（未校正 + 纯度校正成对出现）。E7 的 CE 丰度本身也是成分数据（CE1–CE10 之和为 1），进入线性模型前必须 CLR 或改用 Dirichlet / 单一 CE 的 beta 回归。
- **R-B5** 本文件不预设关联方向。

EN — Rules:

- **R-B1** E4 must finish before E6. If the seeds disagree, E6 is demoted to exploratory; a `BARRIER5` score must not paper over single-gene divergence.
- **R-B2** The E7 outcome is **CE abundance (continuous)**, not a hard assignment to the top CE. Hard assignment is descriptive only.
- **R-B3** E1/E2/E3 are reported **symmetrically** for both seeds. Reporting only "the significant gene" is selective reporting.
- **R-B4** Immune outcomes still obey `playbook.md` R0-1 (compositional transform) and R0-2 (unadjusted + purity-adjusted as a pair). CE abundances in E7 are themselves compositional (CE1–CE10 sum to one) and must be CLR-transformed or modeled with Dirichlet / single-CE beta regression before entering a linear model.
- **R-B5** This file pre-specifies no association direction.

---

## §B4 BayesPrism 臂 / BayesPrism arm

中文 — 运行骨架、单位、`key`、异常基因过滤见 `playbook.md` §4.1。本节只规定**针对双种子**的额外步骤。

EN — Run skeleton, units, `key`, and outlier-gene filters are in `playbook.md` §4.1. This section specifies only the **two-seed** extras.

### §B4.1 参考中的区室特异性检查 / Compartment-specificity check in the reference

中文 — 在跑 bulk 之前，对参考 scRNA-seq 计算每个 `BARRIER5` 基因在各 L1 细胞类型中的平均 counts。写入 `results/qc/barrier5_reference_specificity.tsv`（该文件是 QC，不是关联结果）。规则：

- 若 `TACSTD2` 或 `CLDN4` 在非上皮/非恶性类型中的均值达到上皮/恶性均值的预先设定比例（默认写在 config，建议 0.20），将该基因的 E3 标记为"低区室对比度"，E3 主张自动降为探索性。
- 免疫区室的非零表达保留为诊断，不删基因。

EN — Before touching bulk, compute each `BARRIER5` gene's mean counts per L1 cell type in the scRNA-seq reference. Write `results/qc/barrier5_reference_specificity.tsv` (QC, not an association result). Rule:

- If mean `TACSTD2` or `CLDN4` in a non-epithelial / non-malignant type reaches a pre-specified fraction of the epithelial / malignant mean (pin the fraction in config; suggested default 0.20), mark that gene's E3 as "low compartment contrast" and automatically demote E3 claims to exploratory.
- Non-zero immune-compartment expression is retained as a diagnostic; do not drop the gene.

### §B4.2 恶性区室 CPM / Malignant-compartment CPM

中文 — 对 `BARRIER5` 的每一个基因，用与 `playbook.md` §7.4 路线 A 相同的公式：

$$E^{\text{mal}}_{g,n} = \frac{Z_{g,\,\text{mal},\,n}}{\sum_{g'} Z_{g',\,\text{mal},\,n}} \times 10^6$$

然后取 $\log_2(E^{\text{mal}}_{g,n}+1)$。五个基因必须来自**同一次** InstaPrism/BayesPrism 运行、同一次参考更新设置。禁止把 `TACSTD2` 取自 `update=TRUE`、`CLDN4` 取自 `update=FALSE`。

EN — For every gene in `BARRIER5`, use the same formula as `playbook.md` §7.4 Route A, then $\log_2(E^{\text{mal}}_{g,n}+1)$. All five genes must come from **one** InstaPrism/BayesPrism run and one reference-update setting. Taking `TACSTD2` from `update=TRUE` and `CLDN4` from `update=FALSE` is prohibited.

### §B4.3 可识别性门 / Identifiability gates

| Gate | Default (pin in config) | Fail action |
|---|---|---|
| Malignant fraction floor | $\theta_{\text{mal}} \ge 0.20$ | exclude the sample from all E3/E4/E5/E6 |
| Shared-gene floor vs reference | $\ge 5000$ genes | abort the run |
| `BARRIER5` present in $Z$ | all 5 genes recovered | drop the missing gene from that sample's E6; never impute 0 |
| InstaPrism vs BayesPrism equivalence (optional, ≥1 cohort) | per-gene Spearman of $E^{\text{mal}}$ | record in `results/qc/`; do not pick the "better" one after seeing E3 |

### §B4.4 E4 与 E5 的计算顺序 / Order of E4 and E5

```r
# methods/deconv_advanced/analysis/e4_e5_seeds.R
# Template: confirm API against env/versions.lock before use.
# Inputs are already log2(malignant-CPM + 1). No association with CD8 is computed here.

e4_concordance <- function(e3_t, e3_c) {
  list(
    spearman = suppressWarnings(cor(e3_t, e3_c, method = "spearman")),
    ccc      = lin_ccc(e3_t, e3_c),          # implement or import; pin the definition
    n        = sum(is.finite(e3_t) & is.finite(e3_c))
  )
}

e5_residuals <- function(e3_t, e3_c) {
  # residuals in the malignant-compartment space; purity is not a covariate here
  # because both sides are already within-compartment (playbook §7.6 E3 note)
  list(
    tacstd2_given_cldn4 = resid(lm(e3_t ~ e3_c)),
    cldn4_given_tacstd2 = resid(lm(e3_c ~ e3_t))
  )
}
```

中文 — E4 的一致性本身**不构成免疫主张**。它只决定 E6 是否保持次要层级。E5 的残差随后进入与 §7.6 相同的免疫模型（CD8 CLR、TLS），并对称报告两个方向。

EN — E4 concordance is **not an immune claim**. It only decides whether E6 stays secondary. E5 residuals then enter the same immune models as §7.6 (CD8 CLR, TLS), reported symmetrically in both directions.

### §B4.5 恶性区室上的 BARRIER5 分数 / BARRIER5 score on the malignant compartment

中文 — 把 $Z_{\text{mal}}$（基因 × 样本，区室内 CPM，再 $\log_2(x+1)$）当作表达矩阵，用 `playbook.md` §6.3 的 `singscore` 对 `BARRIER5` 打分。关键：分数在**恶性区室内部**计算，因此不再编码区室大小。禁止用 bulk 表达打 `BARRIER5` 再把它称为 E6。

EN — Treat $Z_{\text{mal}}$ (genes × samples, within-compartment CPM, then $\log_2(x+1)$) as the expression matrix and score `BARRIER5` with `singscore` from `playbook.md` §6.3. The score is computed **inside the malignant compartment**, so it no longer encodes compartment size. Scoring `BARRIER5` on bulk and calling that E6 is prohibited.

---

## §B5 EcoTyper 臂 / EcoTyper arm

中文 — 运行、缩放、z-score 门槛见 `playbook.md` §4.4。本节只规定 recovery 产出如何对接双种子。

EN — Run, scaling, and z-score threshold are in `playbook.md` §4.4. This section specifies only how recovery output meets the two seeds.

### §B5.1 强制前置 / Mandatory prelude

1. 组织学必须在 EcoTyper carcinoma 支持列表中，否则停止 recovery，改用 BayesPrism-only 路径并在报告中声明。
2. 表达矩阵：HUGO symbol、唯一行名、非负；然后 **log2 + 按癌种/数据集单位方差缩放**。缩放脚本独立于 EcoTyper 可执行文件，写入 `methods/deconv_advanced/run/scale_for_ecotyper.R`。
3. $n < 25$：CE 丰度可计算，但 E7 自动标为探索性。
4. 未通过 permutation z-score（默认 1.65）的 CS 整体剔除，见 `playbook.md` §1.4。

### §B5.2 本协议使用的 EcoTyper 产出 / Outputs this protocol consumes

| File (recovery layout) | Use |
|---|---|
| `Carcinoma_Cell_States/Epithelial.cells/*_Cell_State_Abundance.txt` | E8 exposure-side mapping |
| `Carcinoma_Cell_States/CD8.T.cells/*_Cell_State_Abundance.txt` | secondary immune-state corroboration (not a substitute for BayesPrism $\theta_{\text{CD8}}$) |
| `Carcinoma_Cell_States/B.cells/*` and `PCs/*` | TLS-side corroboration (`playbook.md` §6.4) |
| `Carcinoma_Ecotypes/Ecotype_Abundance.txt` | **E7 primary community outcome** |
| `Carcinoma_Ecotypes/Ecotype_Assignment.txt` | descriptive only |
| `dropped_states.tsv` | audit |

### §B5.3 CE 作为成分结局 / CEs as a compositional outcome

中文 — `Ecotype_Abundance.txt` 的每一列（或每一行，依版本）是一个单纯形上的向量。规程：

1. 确认 10 个 CE 之和为 1（允许 $10^{-6}$ 容差）；否则中止。
2. 零值用 `zCompositions::cmultRepl` 替换（`playbook.md` §5.3）。
3. 主模型使用预先指定的四个 CE 的 **ALR 或单一 CE 的 CLR 坐标**：`CE1`、`CE9`、`CE2`、`CE10`。不得把全部 10 个 CLR 坐标同时放进一个回归。
4. 预先指定的社区主结局是 **CLR(CE1)** 与 **CLR(CE9)**（淋巴细胞缺乏 vs IFN-γ 群落）。CE2、CE10 为次要。其余 CE 只作描述。

EN — Each sample in `Ecotype_Abundance.txt` is a vector on the simplex. Procedure:

1. Confirm the 10 CEs sum to 1 (tolerance $10^{-6}$); otherwise abort.
2. Replace zeros with `zCompositions::cmultRepl` (`playbook.md` §5.3).
3. The primary model uses **ALR or single-CE CLR coordinates** for the four pre-specified CEs: `CE1`, `CE9`, `CE2`, `CE10`. Never place all 10 CLR coordinates in one regression.
4. Pre-specified community primary outcomes are **CLR(CE1)** and **CLR(CE9)** (lymphocyte-deficient vs IFN-γ community). CE2 and CE10 are secondary. Remaining CEs are descriptive.

### §B5.4 禁止用 EcoTyper 做的事 / What EcoTyper must not be asked to do here

- 从 EcoTyper 输出回推 `TACSTD2` 或 `CLDN4` 的每细胞表达。
- 把上皮 CS 丰度当作肿瘤纯度。
- 在本 ICI 队列上做 de novo discovery 然后用同一队列做 E7（泄漏）。Discovery 若要进行，必须在无 ICI 标签的大型队列（如 TCGA LUAD/LUSC，分开跑），锁定后再 recovery 到 ICI 队列——且即便如此，本协议的主路径仍是预建 Carcinoma 模型的 recovery。
- 在不支持的组织学上强行 recovery。

---

## §B6 对照协议 / Head-to-head protocol

### §B6.1 共同输入 / Shared inputs

中文 — 两条臂必须吃**同一份** bulk 矩阵、同一套样本 QC（`playbook.md` §2.4）、同一套纯度来源（§7.3）。样本交集在关联前冻结为 `data/processed/{cohort_id}/bp_ecotyper_shared_samples.txt`。只在一条臂上可用的样本不得进入对照，可进入该臂的单臂描述。

EN — Both arms must consume **the same** bulk matrix, the same sample QC (`playbook.md` §2.4), and the same purity source (§7.3). Freeze the sample intersection before association as `data/processed/{cohort_id}/bp_ecotyper_shared_samples.txt`. Samples available on only one arm must not enter the head-to-head; they may enter that arm's single-arm description.

### §B6.2 对照不是比赛 / This is not a bake-off

中文 — 不评"谁赢"。预先定义四类**一致性格局**，每类对应一种允许的写法。格局在看到数字之前用规则定义，用数字去归类，而不是反过来。

EN — There is no "winner." Four **concordance patterns** are defined a priori; each licenses a specific form of sentence. Patterns are defined by rules before numbers are seen; numbers classify a pattern, they do not invent one.

| Pattern | BayesPrism side (E3-T, E3-C, E5) | EcoTyper side (E7) | Permitted sentence |
|---|---|---|---|
| **P-agree** | Both seeds, same sign vs CD8/TLS (or E5 shows shared information) | Same-direction association with CLR(CE1) and/or opposite-direction with CLR(CE9) | The two methods agree on direction for the named estimands |
| **P-gene-split** | The two seeds disagree with each other (E4 weak or opposite) | CE associations follow one seed but not the other | Report the split; name the seed that carries the CE association; do not average |
| **P-method-split** | Seeds agree with each other and with CD8/TLS | CE associations null or opposite | State that compartment dose and community membership **diverge**; this is a first-class finding |
| **P-null** | Seeds vs CD8/TLS sit inside the gene-matched empirical null (`playbook.md` §7.8) | CE associations likewise inside their own null | No association claim; report the diagnostic |

中文 — `P-method-split` 不是失败。它的方法学含义是：恶性细胞上的屏障基因剂量与样本所属的多细胞群落可以脱钩——例如高 $E^{\text{mal}}_{\texttt{TACSTD2}}$ 同时落在 CE9。这种脱钩必须原样报告，不得通过改 CE 定义或改区室门去消除。

EN — `P-method-split` is not a failure. Its methodological meaning is that barrier-gene dose on malignant cells can decouple from the sample's multicellular community — for example high $E^{\text{mal}}_{\texttt{TACSTD2}}$ landing in CE9. Report the decoupling as-is; do not dissolve it by redefining CEs or moving the compartment floor.

### §B6.3 映射表（E8） / Mapping table (E8)

中文 — 对每个通过 z-score 门槛的上皮 CS，计算其丰度与 E3-T、E3-C 的 Spearman（队列内）。这是**赋值探索**，预先声明为探索性，使用 BH-FDR 在"上皮 CS × 两个种子"这一族内控制。禁止在看到热图后再挑选"最相关的 CS"作为主结局。

EN — For every epithelial CS that passes the z-score threshold, compute Spearman of its abundance against E3-T and E3-C (within cohort). This is an **assignment exploration**, pre-declared as exploratory, with BH-FDR inside the family "epithelial CS × two seeds." Do not, after seeing a heatmap, promote "the most correlated CS" to a primary endpoint.

### §B6.4 执行顺序 / Execution order

```
0. Freeze this file + playbook.md §8.0 prereg (no TACSTD2/CLDN4 associations computed yet)
1. Reference QC for BARRIER5 compartment specificity (§B4.1)
2. BayesPrism/InstaPrism on the shared sample set (playbook §4.1)
3. Extract E^mal for BARRIER5; apply identifiability gates (§B4.2–B4.3)
4. E4 concordance; then E5 residuals; then E6 only if E4 permits (§B4.4–B4.5)
5. EcoTyper recovery on the same samples (playbook §4.4, §B5)
6. CLR-transform CE abundances (§B5.3)
7. E3-T / E3-C / E2 / E1 vs CD8 and TLS (playbook §7.6), both seeds
8. E7: E3-T / E3-C / E6 vs CLR(CE1) and CLR(CE9)
9. E8 mapping table (exploratory)
10. Classify P-agree / P-gene-split / P-method-split / P-null (§B6.2)
11. Apply §B7 stability grade
12. Write results/ with section citations; nothing in this file changes
```

中文 — 第 0 步之前不得计算任何 `TACSTD2`/`CLDN4` 与免疫或 CE 的相关。关联脚本在 `prereg.md` 缺少有效 commit hash 时必须拒绝运行（与 `playbook.md` §8.6 同一机制）。

EN — No `TACSTD2`/`CLDN4` association with an immune or CE readout may be computed before step 0. The association script must refuse to run if `prereg.md` lacks a valid commit hash (same mechanism as `playbook.md` §8.6).

### §B6.5 对照指标 / Contrast metrics

| Metric | Applied to | Role |
|---|---|---|
| Spearman $\rho$ + Lin's CCC | E3-T vs E3-C (E4) | seed concordance |
| Partial Spearman | E3-T vs CD8 \| E3-C, and reverse (E5) | unique information |
| $\hat{\beta_1}$ from `playbook.md` §7.6 | E3-T, E3-C, E6 vs CLR(CD8) and TLS | gene-wise immune association |
| $\hat{\beta_1}$ from the same formula | E3-T, E3-C, E6 vs CLR(CE1), CLR(CE9) | community association (E7) |
| Sign concordance table | all of the above across cohorts | input to §B7 |
| Gene-matched empirical null | same null construction as `playbook.md` §7.8, matched on **both** mean and purity-correlation, built **separately** for `TACSTD2` and `CLDN4` | artifact gate |

中文 — 两个种子的经验零分布必须**分开构建**。`CLDN4` 与纯度的相关结构不必与 `TACSTD2` 相同；共用一个零分布会让其中一个基因的零假设错位。

EN — Build the gene-matched empirical null **separately** for each seed. `CLDN4`'s correlation-with-purity structure need not match `TACSTD2`; a shared null misplaces one of the two null hypotheses.

---

## §B7 稳定性分级 / Stability grading

中文 — 在 `playbook.md` §8.4 的 S1–S4 之上，本对照另加一维：**方法间格局**（§B6.2）。完整主张必须同时带上 S 级与 P 类。

EN — On top of S1–S4 in `playbook.md` §8.4, this contrast adds one dimension: the **between-method pattern** (§B6.2). A complete claim carries both an S grade and a P class.

| Grade | Additional requirement for a TACSTD2/CLDN4 claim in this file |
|---|---|
| **S1** | S1 as in §8.4, **and** pattern is P-agree for both seeds |
| **S2** | S2 as in §8.4, or P-gene-split with the carrying seed named, or P-method-split with the divergence named |
| **S3** | Sign flips across BayesPrism vs EcoTyper with no identifiable axis, or E4 and E7 contradict without fitting P-method-split (e.g. because E4 itself is unstable across cohorts) |
| **S4** | Either seed sits inside its own empirical null, or EcoTyper scaling/QC failed, or BayesPrism identifiability gates drop the cohort below the association floor |

中文 — 允许的摘要句模板（填空，不预填方向）：

> Under estimand {E3-T \| E3-C \| E5 \| E7}, primary method {InstaPrism \| EcoTyper-recovery}, stability {S1–S4}, pattern {P-agree \| P-gene-split \| P-method-split \| P-null}, a {positive \| negative \| null} association with {CLR(CD8) \| TLS \| CLR(CE1) \| CLR(CE9)} was observed in {cohort / meta-analysis}.

缺少 S 或 P 的句子视为不完整，不得进入摘要。

EN — Permitted summary-sentence template (fill in; do not pre-fill a direction):

> Under estimand {E3-T \| E3-C \| E5 \| E7}, primary method {InstaPrism \| EcoTyper-recovery}, stability {S1–S4}, pattern {P-agree \| P-gene-split \| P-method-split \| P-null}, a {positive \| negative \| null} association with {CLR(CD8) \| TLS \| CLR(CE1) \| CLR(CE9)} was observed in {cohort / meta-analysis}.

A sentence missing S or P is incomplete and may not enter an abstract.

---

## §B8 预注册冻结清单 / Freeze checklist

中文 — 写入 `methods/deconv_advanced/prereg.md` 的本对照专节，在第 0 步提交。

EN — Write this block into `methods/deconv_advanced/prereg.md` at step 0.

- [ ] Carcinoma EcoTyper model version (commit / release)
- [ ] BayesPrism / InstaPrism version and `update` setting
- [ ] Shared sample list path
- [ ] Purity source and circularity tier (`playbook.md` §7.3)
- [ ] Malignant-fraction floor
- [ ] `BARRIER5` locked; no gene added
- [ ] Primary estimands: E3-T, E3-C; co-primary: E2-T, E2-C, E7 (CLR(CE1), CLR(CE9))
- [ ] Mandatory controls: E1-T, E1-C, E4
- [ ] Secondary: E5, E6 (E6 conditional on E4), CE2, CE10
- [ ] Exploratory: E8
- [ ] Immune outcomes: L2 CLR(CD8), primary TLS signature (`playbook.md` §6)
- [ ] Cross-cohort: two-stage random-effects (`playbook.md` R0-5)
- [ ] Empirical-null construction: separate for each seed
- [ ] Histology restriction for EcoTyper declared
- [ ] Commit hash recorded

---

## §B9 陷阱 / Pitfalls

| ID | Pitfall | Why it matters here | Block |
|---|---|---|---|
| P-B1 | Collapsing `TACSTD2` and `CLDN4` into one score before E4 | Hides seed discordance; E6 then launders a split as a program | R-B1 |
| P-B2 | Treating `CLDN4` as the TROP2 binding partner | Published physical partner is `CLDN7`; `CLDN4` is a co-expression hypothesis | §B1.1 |
| P-B3 | Scoring `BARRIER5` on bulk and calling it compartment-resolved | Reimports purity into E6 | §B4.5 |
| P-B4 | Using EcoTyper CE assignment (argmax) as the E7 outcome | Discards abundance; winner-take-all on a simplex | §B5.3, R-B2 |
| P-B5 | Recovering the carcinoma model in melanoma / GBM / SCLC | Outside the supported histology list; CS/CE labels do not transfer | §B2 |
| P-B6 | De novo EcoTyper discovery on the ICI cohort used for E7 | Outcome-adjacent leakage; states become cohort-specific | §B5.4 |
| P-B7 | Comparing BayesPrism $\theta_{\text{CD8}}$ to EcoTyper CE9 as if they were the same quantity | Fraction ≠ community | §B2, §B6.2 |
| P-B8 | One empirical null shared by both seeds | Misplaced null for the gene whose purity correlation differs | §B6.5 |
| P-B9 | Different InstaPrism `update` flags for the two seeds | $Z$ is then not jointly generated | §B4.2 |
| P-B10 | Promoting the hottest epithelial CS from E8 to primary | Post-hoc endpoint substitution | §B6.3 |
| P-B11 | Adjusting E7 for BayesPrism $\theta_{\text{mal}}$ | Circularity: EcoTyper epithelial CS and BayesPrism malignant fraction both track tumor content | use DNA/methylation purity, `playbook.md` §7.3 |
| P-B12 | Reporting only the seed that "worked" | Selective reporting | R-B3 |
| P-B13 | Skipping EcoTyper per-dataset scaling | Silent invalid recovery (`playbook.md` P7) | §B5.1 |
| P-B14 | Writing a directional abstract without S and P | Reader cannot see the applicability boundary | §B7 |

---

## §B10 参考文献 / References

中文 — 仅定位方法与基因角色，不支持任何结果性陈述。

EN — Locate methods and gene roles only; support no result-level statement.

1. Chu T, et al. BayesPrism. *Nature Cancer* 3:505–517 (2022).
2. Hu M, Chikina M. InstaPrism. *Bioinformatics* 40:btae440 (2024).
3. Luca BA, et al. EcoTyper carcinoma cell states and ecosystems. *Cell* 184:5482–5496 (2021).
4. Steen CB, et al. EcoTyper in DLBCL. *Cancer Cell* 39:1422–1437 (2021).
5. Newman AM, et al. CIBERSORTx (HiRes input to EcoTyper discovery). *Nature Biotechnology* 37:773–782 (2019).
6. Multicellular immune ecotypes predict real-world ICI benefit (EcoTyper recovery in ORIEN). *Nature Communications* (2025). DOI: 10.1038/s41467-025-65016-3
7. TROP2/claudin program and immune exclusion (TROP2–claudin-7 physical interaction). *Journal for ImmunoTherapy of Cancer* 14:e012265 (2026). DOI: 10.1136/jitc-2025-012265
8. Pan-cancer Trop2 multi-omic integration (TACSTD2 co-expression with CLDN4/CLDN7). *npj Precision Oncology* (2026). DOI: 10.1038/s41698-026-01523-w
9. Bessede A, et al. TROP2 and atezolizumab primary resistance (OAK/POPLAR). *Clinical Cancer Research* (2024). PMID 38048058
10. Aitchison J. *The Statistical Analysis of Compositional Data* (1986). — applies to both $\theta$ and CE abundances.

---

*END — 数值结果写入 `results/`，引用 §B#。 / Write numeric results to `results/`, citing §B#.*
