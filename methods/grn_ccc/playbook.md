<!-- markdownlint-disable MD013 MD033 -->
# Playbook — GRN & cell–cell communication for the TACSTD2/CLDN4 → T-cell-recruitment question

> **Languages:** [English](#english) · [中文](#中文)
>
> **Scope of this document.** A practical, opinionated guide to deciding whether
> **TACSTD2/CLDN4-high lung epithelium shows a reduced T-cell recruitment
> signal** on 10x ICI scRNA-seq, using **LIANA/CellChat**, **NicheNet**, and
> **SCENIC/pySCENIC** — and, just as importantly, where each of those tools
> **overclaims**. Companion code lives in [`scripts/`](scripts/) and a runnable
> [`demo/`](demo/).

---

<a name="english"></a>

# English

## 0. TL;DR (the honest verdict up front)

- The four tool families answer **different sub-questions** and none of them, on
  its own, establishes "reduced T-cell recruitment". They all rest on
  **co-expression** or **generic priors**, not on observed signaling or cell
  movement.
- The claim you can actually defend comes from a **simple, patient-level test**:
  pseudobulk the recruitment chemokines (CXCL9/10/11, CCL5, CXCL16) in
  TACSTD2/CLDN4-high vs -low epithelium and compare **across patients**
  ([`scripts/05_pseudobulk_guardrail.py`](scripts/05_pseudobulk_guardrail.py)).
- Use **LIANA/CellChat** to say *which receptors on which T/NK subsets* and
  *which competing senders* are involved — as ranked hypotheses, not proof.
- Use **NicheNet** for a *different* question (does the high state reshape T-cell
  transcriptional **state**), because recruitment is chemotaxis, not a
  transcriptional target program NicheNet models.
- Use **SCENIC** to nominate the **transcription factors** of the TACSTD2/CLDN4
  state and to ask, as a hypothesis, whether chemokine loci are anticorrelated
  with those regulons.
- Then **triangulate** with orthogonal data (spatial, bulk/TCGA, IHC). Report
  effect sizes and patient counts, not cell-level p-values.
- Two public processed objects are **<2 GB** and usable: **GSE207422** (175 MB,
  BD Rhapsody, default real demo) and **GSE205335** (500 MB, 10x). Do not
  download raw FASTQ. Catalog: [`config/public_datasets.yaml`](config/public_datasets.yaml).

## 1. The question, stated precisely

"Reduced T-cell recruitment signal" is ambiguous; pin it to a measurable claim
before touching a tool. Three distinct readouts, in increasing order of how much
they overclaim:

1. **Ligand abundance (cleanest).** TACSTD2/CLDN4-high epithelium expresses
   *fewer* T-cell-recruiting chemokines than the low state. Directly measurable;
   no communication model needed.
2. **Inferred communication (LR tools).** The high state has *weaker inferred
   signaling* to T/NK cells through recruitment axes. Depends on a co-expression
   model.
3. **Functional recruitment (strongest, not answerable by scRNA alone).** Fewer
   T cells are actually recruited to high-state regions. Needs spatial / imaging
   / functional data.

This playbook centers readout (1) as the primary evidence and uses the tools to
enrich readouts (2) and, cautiously, hypotheses toward (3).

## 2. Biological priors you are encoding

- **TACSTD2 (TROP2)** and **CLDN4** mark a common malignant/reactive lung
  adenocarcinoma epithelial state (TROP2 is itself an ADC target). Also present
  on some *normal* reactive epithelium — hence you must separate malignant vs
  normal (§7).
- **T-cell recruitment chemokines** and their receptors (the readout genes,
  defined in [`config/gene_sets.yaml`](config/gene_sets.yaml)):

  | Axis | Ligands (sender = epithelium) | Receptor (receiver = T/NK) | Note |
  |------|------|------|------|
  | CXCR3 | CXCL9, CXCL10, CXCL11 | CXCR3 | Canonical effector-CD8/Th1 axis; strongest ICI-response link |
  | CCR5 | CCL5, CCL3, CCL4 | CCR5, CCR1 | Cytotoxic/effector-memory |
  | CXCR6 | CXCL16 | CXCR6 | Tissue-resident memory CD8 |
  | CCR7 | CCL19, CCL21 | CCR7 | Naive/Tcm homing, TLS |
  | CXCR5 | CXCL13 | CXCR5 | TLS / B–Tfh |

- **Immune-exclusion hypothesis:** epithelial states that silence the CXCL9/10/11
  →CXCR3 axis correlate with "cold"/excluded tumors and worse ICI benefit. That
  is the biology behind "reduced recruitment".
- **Counter-signals:** CCL17/CCL22→CCR4 (Treg/Th2) and CXCL1/2/5/8, CCL2 (myeloid
  /MDSC) recruit *suppressive* cells. A gain here is **not** "more anti-tumor
  recruitment" — track them so you don't misread direction.

## 3. Analysis plan (which tool for which sub-question)

```
                     ┌─────────────────────────────────────────────┐
 raw 10x lung ICI ──▶│ 00_preprocess: ambient+doublet removal, QC,  │
                     │ annotate, malignant call, define epi_state    │
                     └───────────────┬─────────────────────────────┘
                                     │
        ┌────────────────────────────┼───────────────────────────────┐
        ▼                            ▼                                 ▼
 PRIMARY EVIDENCE            HYPOTHESIS GENERATION              MECHANISM (TF)
 05_pseudobulk_guardrail     01_liana / 02_cellchat            04_pyscenic (+04b)
 patient-level chemokine     which receptors/senders;          which TFs define the
 high vs low                 03_nichenet: T-cell STATE effects  state; chemokines in/out
        │                            │                                 │
        └────────────┬───────────────┴─────────────────────────────────┘
                     ▼
        TRIANGULATE: spatial (Visium/Xenium proximity), bulk/TCGA
        (TACSTD2 vs CXCL9/10 & CD8 signature), IHC/mIF (CD8 infiltration)
                     ▲
        06_public_lung_loader  →  GSE207422 / GSE205335 into the same schema
        07_scenic_ccc_joint    →  conservative join of 05 + 01 + 04b
```

Decision rule: a "reduced recruitment" conclusion needs **agreement between the
pseudobulk test and at least one orthogonal modality**. LR/GRN tools sharpen the
mechanism; they do not license the conclusion by themselves. [`07_scenic_ccc_joint.py`](scripts/07_scenic_ccc_joint.py)
will only call a "consistent hypothesis" when the **primary** (05) supports and
at least one of LIANA/SCENIC agrees in direction — and it still refuses to call
that proof.

## 4. Cross-cutting overclaim traps (read before any tool)

These break *all* of LIANA/CellChat/NicheNet/SCENIC if ignored, and each can, on
its own, manufacture a fake "reduced recruitment" result.

1. **Co-expression ≠ communication.** Every LR tool infers interaction from mean
   ligand + receptor expression per cluster. It never observes secretion,
   diffusion, receptor engagement, or **spatial proximity**. Two clusters may be
   in different tissue compartments and still score highly.
2. **Ambient RNA (the #1 threat here).** Dissociated solid lung leaks
   highly-expressed epithelial transcripts — including keratins, surfactant, and
   TACSTD2/CLDN4 themselves — into every droplet. This both (a) fabricates
   epithelial-ligand signal in non-epithelial cells and (b) contaminates the
   epithelial compartment. Run **CellBender / SoupX / decontX** *before*
   clustering, and treat any "ligand expressed in epithelium" number as suspect
   until you have. A depth/contamination difference between the high and low
   states can look exactly like "reduced recruitment".
3. **Doublets.** Epithelial–T doublets create chimeric profiles and phantom LR
   pairs (a "T-cell" that expresses epithelial chemokines). Remove with
   **Scrublet / scDblFinder**.
4. **Pseudoreplication (the #1 statistical error).** Cells are **not**
   independent replicates. CellPhoneDB/LIANA/CellChat permutation p-values
   shuffle **cluster labels within one sample** — they test *specificity*, not
   reproducibility across patients. A "significant" result with n=1 patient is
   meaningless for a between-condition claim. Aggregate to **patient-level
   pseudobulk** and use mixed models / limma-voom / DESeq2 with patient as the
   unit. Require ≥3 patients per group.
5. **Compositional & abundance bias.** More cells → more "significant"
   interactions; magnitude (expression) and specificity (rank) are different
   axes. Do not compare raw interaction counts between conditions with different
   cell numbers.
6. **Batch effects & integration.** Run CCC/DE on **log-normalized,
   uncorrected** values; use Harmony/scVI **only** for the neighbor graph/embedding.
   Feeding imputed/corrected values into LR or DE inflates or erases signal.
7. **Dropout at low expression.** Chemokines are sparsely detected; a difference
   in detection can be dropout, not biology. Always report the **fraction of
   cells expressing** alongside the level.
8. **Directionality is assumed, not measured.** The database says ligand→receptor;
   scRNA cannot tell you the flux direction or that it is active.

## 5. Tool-by-tool

### 5.1 LIANA — ligand–receptor consensus ([`01_liana_lr.py`](scripts/01_liana_lr.py))

**What it computes.** A consensus rank over several LR methods (CellPhoneDB,
NATMI, Connectome, logFC, SingleCellSignalR, geometric mean). Each scores a
ligand in the sender cluster against a receptor in the receiver cluster,
producing `magnitude_rank` (how strong) and `specificity_rank` (how cluster-specific).

**How to map it to the question.** Split epithelium into `Epi_TACSTD2high` /
`Epi_TACSTD2low` as separate senders, run `rank_aggregate`, then compare the
same LR pair from high vs low into T/NK (see `01`, §3–4). For a defensible
between-patient statement, run `rank_aggregate.by_sample` and test the focus
pairs with a **patient-level** Mann–Whitney (`01`, §5).

**Where LIANA overclaims.**
- The consensus dampens method-specific quirks but inherits the shared
  **co-expression assumption** — it still cannot see space or signaling.
- Different member methods mix **magnitude vs specificity**; a top consensus rank
  can be driven by one axis. Say which.
- **Complexes** are reduced to the min of subunits; heteromeric receptors are
  approximate.
- Default analysis is **within-sample**; the naive high-vs-low comparison is
  descriptive. Only the `by_sample` + patient-level test supports a claim.
- It will pair a chemokine with **any** resource receptor that clears the
  expression threshold — including biologically implausible ones (the demo shows
  CXCL10 "interacting" with `ADRA2A`/`GRM7`). **Filter to the receptors you
  care about.**

### 5.2 CellChat — pathway-level communication ([`02_cellchat.R`](scripts/02_cellchat.R))

**What it computes.** A "communication probability" per LR pair via a
law-of-mass-action + Hill function over the curated **CellChatDB**, then
aggregates LR pairs into **signaling pathways** (e.g. the CXCL, CCL families).

**How to map it to the question.** Build one CellChat object per state
(high/low), extract significant LR from epithelium→T/NK on the recruitment
ligands, and compare pathway "information flow" with `rankNet` after
`mergeCellChat`.

**Where CellChat overclaims.**
- The **"probability" is a heuristic score, not a probability**, and moves with
  `type` (triMean is conservative, needs ~≥25% expressing; truncatedMean is
  looser), `population.size`, and `trim`. **Report these knobs**; results change
  with them.
- **Pathway aggregation** can create or mask a signal by summing LR pairs.
- Comparing two objects (`mergeCellChat`, `netVisual_diffInteraction`, `rankNet`)
  is **descriptive** — the permutation is per-object label shuffling, again **not**
  a cross-patient test.
- **CellChatDB is literature-curated** and incomplete; absence of a pathway is
  not evidence of absence.

### 5.3 NicheNet — ligand→target activity ([`03_nichenet.R`](scripts/03_nichenet.R))

**Read this first — scope.** NicheNet links a sender **ligand** to a
**transcriptional** response in the receiver via a generic prior
(ligand→signaling→target genes). **Recruitment is chemotaxis/migration, which is
largely not a transcriptional target program.** So NicheNet is the **wrong
primary tool for "recruitment"** and will rank chemokines for the wrong reasons.
Use it for the legitimate adjacent question: *does the TACSTD2/CLDN4-high state
send ligands that reshape the transcriptional STATE of infiltrating T cells*
(activation, dysfunction/exhaustion)?

**What it computes.** "Ligand activity" = how well a ligand's prior-predicted
target genes match a **receiver DE gene set** (AUPR/Pearson vs a background).

**Where NicheNet overclaims.**
- By default it does **not** require the ligand or its receptor to be **expressed**
  — a "top ligand" may be absent from your data. Always add expression filters
  (the template does).
- The prior model is **tissue/context-agnostic** and integrated from generic
  public networks. It predicts *potential* regulation, correlationally.
- **Everything hinges on `geneset_oi`** (the receiver DE set). A sloppy DE set →
  meaningless rankings. Define it with a proper contrast (ideally patient-level
  pseudobulk DE) and **re-run with a perturbed gene set to check stability**.
- It does not tell you the sender; you assign senders post hoc by expression.
- Never phrase output as "epithelium **causes** T-cell state X".

### 5.4 SCENIC / pySCENIC — TF regulons ([`04_pyscenic.sh`](scripts/04_pyscenic.sh), [`04b_scenic_downstream.py`](scripts/04b_scenic_downstream.py))

**What it computes.** GRNBoost2 (TF–target co-expression) → cisTarget (keep
targets bearing the TF motif) → AUCell (per-cell regulon activity). Output:
regulons (TF + pruned targets) and a cell×regulon activity matrix.

**How to map it to the question.** Find regulons whose activity separates the
high vs low epithelial state (patient-aware, `04b` §2), then ask whether the
recruitment chemokines are (a) **inside** any high-state regulon and (b)
**anticorrelated** with high-state regulon activity — a *hypothesis* that a TF
of the TACSTD2/CLDN4 program represses them.

**Where SCENIC overclaims.**
- **Correlation, not causation.** A regulon whose activity tracks the state does
  not prove the TF drives (or represses) the state or the chemokines.
- **GRNBoost2 is stochastic.** Run ≥3 seeds and keep a **consensus**; single-run
  regulons are unstable.
- **cisTarget databases are genome- and TF-specific** (hg38/hg19/mm10), cover
  only TFs with curated motifs, and model **activation** — **repression is not
  directly modeled**, yet "reduced recruitment via repression" is exactly a
  repression hypothesis. You cannot conclude repression from SCENIC alone; you
  need motif presence in the chemokine locus (cisTarget/ATAC) + anticorrelation
  + orthogonal support.
- **Dropout** degrades GRN inference; AUCell binarization thresholds are heuristic.
- SCENIC is **intracellular** — it says nothing about cell–cell communication by
  itself.

## 6. Practicality on 10x lung ICI specifically

- **Malignant vs normal epithelium.** TACSTD2/CLDN4 are high in malignant *and*
  some reactive normal epithelium. Run **inferCNV / CopyKAT / numbat** and define
  the state within the malignant compartment, or you will confound tumor biology
  with normal alveolar/airway epithelium.
- **Per-patient heterogeneity & CNV.** Malignant cells are aneuploid and
  patient-private; the TACSTD2/CLDN4-high state may look different per patient.
  Prefer within-patient high-vs-low contrasts, then meta-analyze across patients.
- **Ambient contamination is severe** in dissociated solid tumor — §4.2 is not
  optional here.
- **Compute.** Lung ICI atlases are large (10^5–10^6 cells). GRNBoost2 and
  cisTarget are heavy: **subsample epithelium** (e.g. ≤20k cells) for GRN
  inference, then AUCell-score everyone. cisTarget DBs are multi-GB downloads.
  LIANA/CellChat scale better but still benefit from restricting to relevant
  lineages.
- **Small N patients.** ICI cohorts are often tens of patients, sometimes with
  pre/post and responder/non-responder splits. Power for between-group tests is
  limited — report confidence intervals and treat single-cohort results as
  provisional.
- **Sparsity.** Chemokines are lowly expressed; pseudobulk is far more reliable
  than per-cell tests.

## 6b. Public lung ICI objects — what is actually downloadable

The playbook rule is **"demo only if a small scRNA <2 GB is available"**. Raw
10x FASTQ / Cell Ranger outs fail that rule. Authors' **processed** matrices
sometimes pass. Checked 2026-08-16; full table in
[`config/public_datasets.yaml`](config/public_datasets.yaml).

| Accession | Platform | ICI? | Epithelium? | Processed size | Demo? |
|-----------|----------|------|-------------|----------------|-------|
| **GSE207422** (Hu 2023) | BD Rhapsody | neoadj. PD-1+chemo, 15 pts, MPR/NMPR | yes (must marker-gate; GEO deposited **no** cell-type column) | **175 MB** UMI + 11 KB sample xlsx | **yes — default** |
| **GSE205335** | 10x 3' | ICI, RECIST on some | yes (CellIdentity table) | **500 MB** `.rds` | yes, needs R→h5ad |
| GSE179994 | 10x | pembro+chemo | **no** (T-cell-only object) | 421 MB | no (no epithelial sender) |
| GSE131907 | 10x | **no** (LUAD atlas) | yes | raw is multi-GB | no (not ICI; sibling PR #96 skipped raw) |
| Caushi 2021 | 10x | neoadj. anti-PD-1 | yes | EGA controlled | no (needs DAC) |

**GSE207422 caveats (do not hide these).**
- It is **BD Rhapsody, not 10x**. Gene space / UMI scale differ; do not drop it
  into a 10x embedding without integration. It is still the smallest public
  *ICI + epithelium + T/NK + response labels* object under 2 GB.
- Metadata xlsx is **sample-level** (15 patients). Lineages must be
  marker-gated ([`demo/run_gse207422_demo.py`](demo/run_gse207422_demo.py)).
- Authors already QC + Scrublet; **ambient RNA cannot be re-estimated**. Record
  that in `.uns['ambient_correction']`.
- n=15, MPR n=4: underpowered for between-response claims. Prefer
  **within-patient** high vs low, then describe.

**GSE205335** is the 10x-native replication (499.5 MB RDS + cell identity).
Convert with `zellkonverter` / `sceasy` before the Python stack.

Loader: [`scripts/06_public_lung_loader.py`](scripts/06_public_lung_loader.py).

**What the GSE207422 demo actually returned** (see
[`demo/gse207422_example/RESULTS.md`](demo/gse207422_example/RESULTS.md)):
13 paired patients; 8/10 recruitment ligands lower in the high state; **0/10
FDR < 0.05**. LIANA's only focus edge (CXCL16–CXCR6) is *stronger* from the
high sender. Per §6c: primary does not support → **do not claim** reduced
T-cell recruitment from this cohort.

## 6c. Joint SCENIC × communication readout

[`scripts/07_scenic_ccc_joint.py`](scripts/07_scenic_ccc_joint.py) joins the
three already-computed tables and issues a **conservative** call:

- primary 05 supports **and** LIANA or SCENIC agrees in direction →
  "consistent hypothesis" (still not proof of recruitment or of TF repression);
- primary 05 supports, others missing/disagree → report the pseudobulk only;
- primary 05 does **not** support → **do not claim** reduced recruitment from
  LIANA / CellChat / NicheNet / SCENIC alone.

That last line is the whole point of this playbook.

## 7. The guardrail analysis (do this regardless) — [`05_pseudobulk_guardrail.py`](scripts/05_pseudobulk_guardrail.py)

1. Within epithelium (malignant), sum **raw counts** per (patient × epi_state) →
   pseudobulk.
2. CPM/log or, better, **limma-voom / DESeq2** with `~ patient + state`.
3. Test the recruitment chemokines high vs low, **paired by patient**.
4. Report **fraction expressing** per state next to the level (dropout check).
5. Only if this shows a real, patient-level reduction *and* an orthogonal
   modality agrees do you write "TACSTD2/CLDN4-high epithelium shows reduced
   T-cell recruitment signal".

Reference DE snippet (R), the version to use for a paper:

```r
# pb: genes x pseudobulk-samples counts; meta: patient, state
library(edgeR); library(limma)
d <- DGEList(pb); d <- calcNormFactors(d)
design <- model.matrix(~ patient + state, meta)   # patient as blocking factor
v <- voom(d, design); fit <- eBayes(lmFit(v, design))
topTable(fit, coef = "stateTACSTD2_CLDN4_high", sort.by = "P")
```

## 8. Reporting checklist (claim / do-not-claim)

- [ ] Ambient RNA correction run and named (tool + version).
- [ ] Doublets removed; malignant epithelium separated from normal.
- [ ] `epi_state` definition (score, method, threshold) recorded in `.uns`.
- [ ] Every between-condition claim is **patient-level** (≥3/group), with effect
      size + CI — **no cell-level p-values as evidence**.
- [ ] LR/GRN results framed as **ranked hypotheses**, with tool knobs reported.
- [ ] Fraction-expressing shown for every "reduced" gene.
- [ ] At least one **orthogonal** validation (spatial / bulk-TCGA / IHC).
- [ ] Direction sanity-checked against suppressive chemokines (not just CXCR3 axis).
- [ ] **Do not** claim causation from NicheNet/SCENIC, or "recruitment" from
      NicheNet, or "communication" from co-expression alone.

## 9. Resources

- LIANA-py — <https://liana-py.readthedocs.io> · Dimitrov et al., *Nat Commun* 2022.
- CellChat v2 — <https://github.com/jinworks/CellChat> · Jin et al., *Nat Commun* 2021/2023.
- NicheNet — <https://github.com/saeyslab/nichenetr> · Browaeys et al., *Nat Methods* 2020; prior model Zenodo 7074291.
- pySCENIC — <https://pyscenic.readthedocs.io> · Aibar et al., *Nat Methods* 2017; Van de Sande et al., *Nat Protoc* 2020.
- Ambient: CellBender (Fleming et al.), SoupX (Young & Behjati), decontX (Yang et al.).
- Malignant calling: inferCNV, CopyKAT (Gao et al.), numbat (Gao et al.).
- Pseudobulk/pseudoreplication: Squair et al., *Nat Commun* 2021; Murphy & Skene 2022.

---
---

<a name="中文"></a>

# 中文

## 0. 一句话结论（先把话说明白）

- 这四类工具回答的是**不同的子问题**，任何单独一个都**不能**证明"T 细胞招募信号
  下降"。它们本质上都建立在**共表达**或**通用先验网络**之上，而**不是**观测到的
  真实信号传导或细胞迁移。
- 真正**站得住脚**的结论来自一个**简单的、以病人为单位的检验**：把招募型趋化因子
  （CXCL9/10/11、CCL5、CXCL16）在 TACSTD2/CLDN4-高 与 -低 上皮中做 pseudobulk，
  **跨病人**比较（见 [`scripts/05_pseudobulk_guardrail.py`](scripts/05_pseudobulk_guardrail.py)）。
- 用 **LIANA/CellChat** 说明*涉及哪些 T/NK 亚群上的受体、哪些竞争性发送细胞*——
  作为**排序后的假设**，而非证据。
- 用 **NicheNet** 回答*另一个*问题（高状态是否重塑 T 细胞的转录**状态**），因为
  招募是趋化/迁移，并不是 NicheNet 所建模的转录靶基因程序。
- 用 **SCENIC** 提名 TACSTD2/CLDN4 状态的**转录因子**，并作为假设去问：趋化因子
  位点是否与这些 regulon 反相关。
- 最后用**正交数据**（空间转录组、bulk/TCGA、IHC）**三角验证**。报告效应量与病人
  数，而不是细胞级 p 值。
- 两份公开处理后对象 **<2 GB** 可用：**GSE207422**（175 MB，BD Rhapsody，默认真实
  demo）与 **GSE205335**（500 MB，10x）。不要下原始 FASTQ。目录见
  [`config/public_datasets.yaml`](config/public_datasets.yaml)。

## 1. 把问题说精确

"T 细胞招募信号下降"含义模糊，动手前先落到可测量的命题。三种读出，按"越往下越容易
过度解读"排序：

1. **配体丰度（最干净）**：TACSTD2/CLDN4-高上皮表达的招募型趋化因子比低状态**更少**。
   可直接测量，不需要任何通讯模型。
2. **推断的通讯（LR 工具）**：高状态经招募轴向 T/NK 的**推断信号更弱**。依赖共表达模型。
3. **功能性招募（最强，scRNA 单独无法回答）**：高状态区域实际招募到的 T 细胞更少。
   需要空间/成像/功能数据。

本手册以读出 (1) 作为**主要证据**，用工具去丰富读出 (2)，并谨慎地为读出 (3) 提出假设。

## 2. 你正在编码的生物学先验

- **TACSTD2（TROP2）**与 **CLDN4** 标记一种常见的肺腺癌恶性/反应性上皮状态（TROP2
  本身也是 ADC 靶点）；部分**正常反应性上皮**也表达，因此必须区分恶性与正常（见 §7）。
- **T 细胞招募趋化因子**及其受体（读出基因，定义于
  [`config/gene_sets.yaml`](config/gene_sets.yaml)）：

  | 轴 | 配体（发送=上皮） | 受体（接收=T/NK） | 说明 |
  |----|------|------|------|
  | CXCR3 | CXCL9、CXCL10、CXCL11 | CXCR3 | 经典效应 CD8/Th1 轴；与 ICI 应答关系最强 |
  | CCR5 | CCL5、CCL3、CCL4 | CCR5、CCR1 | 细胞毒/效应记忆 |
  | CXCR6 | CXCL16 | CXCR6 | 组织驻留记忆 CD8 |
  | CCR7 | CCL19、CCL21 | CCR7 | 初始/中央记忆归巢、TLS |
  | CXCR5 | CXCL13 | CXCR5 | TLS / B–Tfh |

- **免疫排斥假说**：沉默 CXCL9/10/11→CXCR3 轴的上皮状态与"冷"/排斥型肿瘤及较差的
  ICI 获益相关——这正是"招募下降"背后的生物学。
- **反向信号**：CCL17/CCL22→CCR4（Treg/Th2）与 CXCL1/2/5/8、CCL2（髓系/MDSC）招募
  的是**抑制性**细胞。这些升高**不等于**"抗肿瘤招募增加"，务必一并追踪以免读反方向。

## 3. 分析方案（哪个工具回答哪个子问题）

```
                     ┌─────────────────────────────────────────────┐
 原始 10x 肺 ICI ───▶│ 00_preprocess：去环境RNA+去双胞、QC、注释、    │
                     │ 恶性判定、定义 epi_state                       │
                     └───────────────┬─────────────────────────────┘
                                     │
        ┌────────────────────────────┼───────────────────────────────┐
        ▼                            ▼                                 ▼
   主要证据                    假设生成                          机制（TF）
 05_pseudobulk_guardrail    01_liana / 02_cellchat            04_pyscenic (+04b)
 病人级 趋化因子             哪些受体/发送方；                  哪些 TF 定义该状态；
 高 vs 低                    03_nichenet：T 细胞“状态”效应       趋化因子在/不在 regulon
        │                            │                                 │
        └────────────┬───────────────┴─────────────────────────────────┘
                     ▼
        三角验证：空间（Visium/Xenium 邻近）、bulk/TCGA
        （TACSTD2 vs CXCL9/10 及 CD8 signature）、IHC/mIF（CD8 浸润）
                     ▲
        06_public_lung_loader  →  把 GSE207422 / GSE205335 装进同一套 schema
        07_scenic_ccc_joint    →  对 05 + 01 + 04b 做保守汇总
```

判定规则：得出"招募下降"的结论，需要 **pseudobulk 检验与至少一种正交模态相互印证**。
LR/GRN 工具用来**锐化机制**，本身**不足以**支撑结论。[`07_scenic_ccc_joint.py`](scripts/07_scenic_ccc_joint.py)
只有在**主检验**（05）支持、且 LIANA/SCENIC 至少一个方向一致时，才会给出"一致假设"，
并且**拒绝**把这称为证明。

## 4. 贯穿全局的"过度解读"陷阱（用任何工具前必读）

以下问题若被忽视，会**同时**摧毁 LIANA/CellChat/NicheNet/SCENIC，且每一个都能**单独**
制造出假的"招募下降"。

1. **共表达 ≠ 通讯**：所有 LR 工具都用"每个簇的配体+受体平均表达"推断互作，从不观测
   分泌、扩散、受体结合或**空间邻近**。两个簇可能位于不同组织区室，却仍得到高分。
2. **环境 RNA（此处头号威胁）**：实体肺组织解离会把高表达的上皮转录本——包括角蛋白、
   表面活性蛋白，以及 TACSTD2/CLDN4 本身——泄漏进每个液滴，既 (a) 在非上皮细胞里
   伪造"上皮配体"信号，又 (b) 污染上皮区室。务必在聚类**之前**运行 **CellBender /
   SoupX / decontX**；在此之前，任何"上皮表达某配体"的数字都应视为可疑。高/低状态间
   的深度或污染差异，看起来会**和"招募下降"一模一样**。
3. **双胞（doublet）**：上皮–T 双胞会产生嵌合表达和幽灵 LR 对（一个"表达上皮趋化因子
   的 T 细胞"）。用 **Scrublet / scDblFinder** 去除。
4. **伪重复（头号统计错误）**：细胞**不是**独立重复。CellPhoneDB/LIANA/CellChat 的
   置换 p 值在**单个样本内打乱簇标签**，检验的是*特异性*，不是跨病人的可重复性。
   n=1 病人的"显著"结果对组间结论毫无意义。应聚合到**病人级 pseudobulk**，以病人为
   单位用混合模型 / limma-voom / DESeq2，且每组 ≥3 病人。
5. **构成与丰度偏倚**：细胞越多"显著"互作越多；幅度（表达）与特异性（排名）是不同的轴。
   不要在细胞数不同的两个条件间直接比较原始互作计数。
6. **批次效应与整合**：CCC/DE 要在**未校正的对数归一化**值上跑；Harmony/scVI **只**用于
   邻接图/嵌入。把插补/校正后的值喂给 LR 或 DE 会放大或抹掉信号。
7. **低表达 dropout**：趋化因子检出稀疏，检出率差异可能是 dropout 而非生物学。务必在
   报告表达水平的同时报告**表达细胞比例**。
8. **方向是假定的，不是测得的**：数据库规定配体→受体；scRNA 无法告诉你通量方向或其是否活跃。

## 5. 逐工具说明

### 5.1 LIANA —— 配体–受体共识 ([`01_liana_lr.py`](scripts/01_liana_lr.py))

**它算什么**：对多种 LR 方法（CellPhoneDB、NATMI、Connectome、logFC、
SingleCellSignalR、几何平均）做共识排名，得到 `magnitude_rank`（强度）与
`specificity_rank`（簇特异性）。

**如何对应本问题**：把上皮拆成 `Epi_TACSTD2high` / `Epi_TACSTD2low` 作为不同发送方，
运行 `rank_aggregate`，比较同一 LR 对从高/低发向 T/NK 的排名（见 `01` §3–4）。要得到
可辩护的跨病人结论，用 `rank_aggregate.by_sample` 并对重点 LR 对做**病人级** Mann–Whitney
（`01` §5）。

**LIANA 在哪里过度解读**：
- 共识能削弱单一方法的怪癖，但继承了共同的**共表达假设**——依然看不到空间与信号传导。
- 各成员方法混合了**幅度 vs 特异性**；靠前的共识排名可能由其中一个轴主导，需说明。
- **复合体**被简化为亚基最小值；异源受体只是近似。
- 默认是**样本内**分析；朴素的高–低比较只是描述性的，只有 `by_sample`+病人级检验才支撑结论。
- 它会把某趋化因子与**任何**越过表达阈值的资源受体配对——包括生物学上不合理者
  （demo 中 CXCL10 与 `ADRA2A`/`GRM7` "互作"）。**务必只筛你关心的受体。**

### 5.2 CellChat —— 通路级通讯 ([`02_cellchat.R`](scripts/02_cellchat.R))

**它算什么**：基于策展的 **CellChatDB**，用质量作用定律 + Hill 函数给每个 LR 对一个
"通讯概率"，再把 LR 对聚合成**信号通路**（如 CXCL、CCL 家族）。

**如何对应本问题**：为每个状态（高/低）建一个 CellChat 对象，提取上皮→T/NK 在招募配体上
的显著 LR，`mergeCellChat` 后用 `rankNet` 比较通路"信息流"。

**CellChat 在哪里过度解读**：
- **"概率"是启发式分数而非概率**，并随 `type`（triMean 保守，约需 ≥25% 细胞表达；
  truncatedMean 更宽松）、`population.size`、`trim` 变化。**必须报告这些参数**。
- **通路聚合**通过对 LR 对求和，可能凭空造出或掩盖信号。
- 比较两个对象（`mergeCellChat`、`netVisual_diffInteraction`、`rankNet`）是**描述性**的
  ——置换仍是对象内标签打乱，**不是**跨病人检验。
- **CellChatDB 是文献策展**且不完整；某通路"缺失"不代表其不存在。

### 5.3 NicheNet —— 配体→靶基因活性 ([`03_nichenet.R`](scripts/03_nichenet.R))

**先读这里——适用范围**：NicheNet 通过通用先验（配体→信号→靶基因）把发送方**配体**
与接收方的**转录**响应联系起来。**招募是趋化/迁移，基本不是转录靶基因程序**，所以
NicheNet 是**回答"招募"的错误主力工具**，会因错误原因给趋化因子排名。它适合的相邻问题是：
*TACSTD2/CLDN4-高状态是否分泌能重塑浸润 T 细胞转录**状态**（激活、功能失调/耗竭）的配体*。

**它算什么**："配体活性" = 该配体先验预测的靶基因与**接收方 DE 基因集**的匹配程度
（相对背景的 AUPR/Pearson）。

**NicheNet 在哪里过度解读**：
- 默认**不要求**配体或其受体**表达**——"头部配体"可能压根不在你的数据里。务必加表达过滤
  （模板已加）。
- 先验模型**与组织/语境无关**，整合自通用公共网络，预测的是*潜在*调控、相关性层面。
- **一切取决于 `geneset_oi`**（接收方 DE 集）。DE 集马虎 → 排名毫无意义。用正规对比
  （最好病人级 pseudobulk DE）定义它，并**用扰动后的基因集复跑以检验稳定性**。
- 它不告诉你发送方；发送方需按表达事后指派。
- 绝不能把结果表述为"上皮**导致** T 细胞状态 X"。

### 5.4 SCENIC / pySCENIC —— 转录因子 regulon ([`04_pyscenic.sh`](scripts/04_pyscenic.sh)、[`04b_scenic_downstream.py`](scripts/04b_scenic_downstream.py))

**它算什么**：GRNBoost2（TF–靶共表达）→ cisTarget（保留带该 TF motif 的靶）→ AUCell
（每个细胞的 regulon 活性）。输出 regulon（TF+剪枝后的靶）和 细胞×regulon 活性矩阵。

**如何对应本问题**：找出活性能区分高/低上皮状态的 regulon（病人感知，`04b` §2），再问
招募趋化因子是否 (a) **落在**某个高状态 regulon 内、(b) 与高状态 regulon 活性**反相关**
——这是"TACSTD2/CLDN4 程序的某 TF 抑制它们"的*假设*。

**SCENIC 在哪里过度解读**：
- **相关而非因果**：regulon 活性随状态变化，并不证明该 TF 驱动（或抑制）该状态或趋化因子。
- **GRNBoost2 是随机的**：跑 ≥3 个种子取**共识**；单次 regulon 不稳定。
- **cisTarget 数据库依赖基因组与 TF**（hg38/hg19/mm10），只覆盖有策展 motif 的 TF，且
  建模的是**激活**——**不直接建模抑制**，而"经抑制导致招募下降"恰恰是抑制假设。**仅凭
  SCENIC 不能下抑制结论**；需要趋化因子位点存在该 TF 的 motif（cisTarget/ATAC）+ 反相关
  + 正交证据。
- **dropout** 损害 GRN 推断；AUCell 二值化阈值是启发式的。
- SCENIC 是**胞内**方法——本身不谈细胞间通讯。

## 6. 在 10x 肺 ICI 上的实操性

- **恶性 vs 正常上皮**：TACSTD2/CLDN4 在恶性*和*部分反应性正常上皮中都高。运行
  **inferCNV / CopyKAT / numbat**，在恶性区室内定义该状态，否则会把肿瘤生物学与正常
  肺泡/气道上皮混淆。
- **病人间异质性与 CNV**：恶性细胞非整倍、且高度病人特异，TACSTD2/CLDN4-高状态在不同
  病人中可能形态各异。优先做**病人内**高–低对比，再跨病人 meta 分析。
- **环境污染严重**：解离实体瘤时 §4.2 不是可选项。
- **算力**：肺 ICI 图谱很大（10^5–10^6 细胞）。GRNBoost2 与 cisTarget 很重：GRN 推断时
  **对上皮下采样**（如 ≤2 万细胞），再对全部细胞做 AUCell 打分；cisTarget 数据库是数 GB
  下载。LIANA/CellChat 扩展性更好，但也建议只保留相关谱系。
- **病人数 N 小**：ICI 队列常仅数十人，还可能分 pre/post 与 应答/非应答。组间检验效力有限
  ——报告置信区间，单队列结果视为暂定。
- **稀疏性**：趋化因子表达低，pseudobulk 远比逐细胞检验可靠。

## 6b. 公开肺 ICI 对象——真正能下到的是什么

手册规则是**"仅当存在 <2 GB 的小型 scRNA 时才做 demo"**。原始 10x FASTQ / Cell Ranger
输出不满足。作者**处理后的**矩阵有时可以。2026-08-16 核查；完整表见
[`config/public_datasets.yaml`](config/public_datasets.yaml)。

| 登录号 | 平台 | ICI？ | 有上皮？ | 处理后大小 | 做 demo？ |
|--------|------|-------|----------|------------|-----------|
| **GSE207422**（Hu 2023） | BD Rhapsody | 新辅助 PD-1+化疗，15 人，MPR/NMPR | 有（须用 marker 门控；GEO **未**提供细胞类型列） | **175 MB** UMI + 11 KB 样本 xlsx | **是 — 默认** |
| **GSE205335** | 10x 3' | ICI，部分有 RECIST | 有（CellIdentity 表） | **500 MB** `.rds` | 是，需 R→h5ad |
| GSE179994 | 10x | pembro+化疗 | **无**（仅 T 细胞对象） | 421 MB | 否（没有上皮发送方） |
| GSE131907 | 10x | **否**（LUAD 图谱） | 有 | 原始数 GB | 否（非 ICI；姊妹 PR #96 已跳过原始） |
| Caushi 2021 | 10x | 新辅助 anti-PD-1 | 有 | EGA 受控 | 否（需 DAC） |

**GSE207422 必须写明的限制。**
- 这是 **BD Rhapsody，不是 10x**。基因空间 / UMI 尺度不同；未经整合不要丢进 10x 嵌入。
  它仍是 2 GB 以下最小的公开 *ICI + 上皮 + T/NK + 应答标签* 对象。
- metadata xlsx 是**样本级**（15 人）。谱系必须用 marker 门控
  （[`demo/run_gse207422_demo.py`](demo/run_gse207422_demo.py)）。
- 作者已做 QC + Scrublet；**无法再估计环境 RNA**。记入
  `.uns['ambient_correction']`。
- n=15，MPR n=4：组间应答比较效力不足。优先做**病人内**高–低，再描述。

**GSE205335** 是 10x 原生的重复队列（499.5 MB RDS + 细胞身份）。先用
`zellkonverter` / `sceasy` 转成 Python 栈。

装载脚本：[`scripts/06_public_lung_loader.py`](scripts/06_public_lung_loader.py)。

**GSE207422 demo 的实际结果**（见
[`demo/gse207422_example/RESULTS.md`](demo/gse207422_example/RESULTS.md)）：
13 对病人；10 个招募配体中 8 个在高状态更低；**0/10 达到 FDR < 0.05**。LIANA
唯一通过的焦点边（CXCL16–CXCR6）在高发送方反而*更强*。按 §6c：主检验不支持 →
**不得**据此队列声称 T 细胞招募下降。

## 6c. SCENIC × 通讯的联合读出

[`scripts/07_scenic_ccc_joint.py`](scripts/07_scenic_ccc_joint.py) 汇总三张已算好的表，
给出**保守**判定：

- 主检验 05 支持 **且** LIANA 或 SCENIC 方向一致 → "一致假设"（仍不是招募或 TF 抑制的证明）；
- 主检验 05 支持、其余缺失/不一致 → 只报告 pseudobulk；
- 主检验 05 **不**支持 → **不得**仅凭 LIANA / CellChat / NicheNet / SCENIC 声称招募下降。

最后这一条就是本手册的全部要点。

## 7. 护栏分析（无论如何都要做）—— [`05_pseudobulk_guardrail.py`](scripts/05_pseudobulk_guardrail.py)

1. 在（恶性）上皮内，按 (病人 × epi_state) 对**原始 counts** 求和 → pseudobulk。
2. CPM/log，或更好地用 **limma-voom / DESeq2**，模型 `~ patient + state`。
3. 对招募趋化因子做高–低检验，**以病人配对**。
4. 在表达水平旁**报告每个状态的表达细胞比例**（dropout 核查）。
5. 只有当此检验显示真实的病人级下降**且**一种正交模态相互印证时，才写下
   "TACSTD2/CLDN4-高上皮表现出 T 细胞招募信号下降"。

论文级 DE 参考片段（R）：

```r
# pb：基因 × pseudobulk 样本 的 counts；meta：patient、state
library(edgeR); library(limma)
d <- DGEList(pb); d <- calcNormFactors(d)
design <- model.matrix(~ patient + state, meta)   # 以病人作为区组因子
v <- voom(d, design); fit <- eBayes(lmFit(v, design))
topTable(fit, coef = "stateTACSTD2_CLDN4_high", sort.by = "P")
```

## 8. 汇报清单（能说 / 不能说）

- [ ] 已运行并注明环境 RNA 校正（工具+版本）。
- [ ] 已去双胞；已把恶性上皮与正常分开。
- [ ] `epi_state` 定义（打分、方法、阈值）已记录在 `.uns`。
- [ ] 每个组间结论都是**病人级**（每组 ≥3），并给出效应量+置信区间——**不得以细胞级 p 值
      作为证据**。
- [ ] LR/GRN 结果表述为**排序后的假设**，并报告工具参数。
- [ ] 每个"下降"的基因都给出表达细胞比例。
- [ ] 至少一种**正交**验证（空间 / bulk-TCGA / IHC）。
- [ ] 结合抑制性趋化因子核查方向（不要只看 CXCR3 轴）。
- [ ] **不要**从 NicheNet/SCENIC 下因果结论、不要用 NicheNet 谈"招募"、不要仅凭共表达谈"通讯"。

## 9. 资源

- LIANA-py — <https://liana-py.readthedocs.io>；Dimitrov 等, *Nat Commun* 2022。
- CellChat v2 — <https://github.com/jinworks/CellChat>；Jin 等, *Nat Commun* 2021/2023。
- NicheNet — <https://github.com/saeyslab/nichenetr>；Browaeys 等, *Nat Methods* 2020；先验模型 Zenodo 7074291。
- pySCENIC — <https://pyscenic.readthedocs.io>；Aibar 等, *Nat Methods* 2017；Van de Sande 等, *Nat Protoc* 2020。
- 环境 RNA：CellBender（Fleming 等）、SoupX（Young & Behjati）、decontX（Yang 等）。
- 恶性判定：inferCNV、CopyKAT（Gao 等）、numbat（Gao 等）。
- Pseudobulk/伪重复：Squair 等, *Nat Commun* 2021；Murphy & Skene 2022。
