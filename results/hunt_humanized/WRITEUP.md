# Public CD34 / PBMC-humanized NSCLC IO RNA/scRNA: TACSTD2/CLDN4 vs CD8 nest entry

Hunt date: 2026-08-16. Outputs: `results/hunt_humanized/`. Reanalysis: `scripts/hunt_humanized/analyze_public.py`.

## English

### Question

In **public** CD34-humanized or PBMC-humanized NSCLC immuno-oncology RNA or scRNA, is there any evidence that **TACSTD2 (TROP2)** or **CLDN4** tracks **CD8 nest entry**?

Nest entry means CD8 T cells inside epithelial tumor nests, not stromal CD8 and not bulk `CD8A`/`CD8B` CPM.

### Methods (what was searched)

NCBI GEO eSearch (humanized / CD34 / PBMC × NSCLC / lung / H460 / A549 / H1975 / H358 / H226 × RNA-seq / scRNA-seq / pembrolizumab / PD-1). ArrayExpress E-MTAB-9891. Literature: Meraz 2019 CIR and 2024 eLife; Yu 2023 Cell Commun Signal (Shenmai); Na 2024 *Life* (BJIKT, GSE260575); Kim 2025 MCT (denfivontinib, GSE276724); Kwok 2025/2026 CCL20 DTP (GSE293914); SITC 2024 abstract 807 (TROP2 ADC + sabestomig); AACR 2024 abstract 4186; Chiorazzi 2023 autologous MISTRG6. False-positive filters: patient blood PBMC scRNA; “humanized” antibody; hHGF-NSG cytokine mice; syngeneic mouse GEO that only mention HIS in the paper.

Opened matrices: GSE260575 count files; GSE276724 xlsx. GSE293914 file list inspected (barcodes / features / 462 MB MTX, no metadata). Yu 2023 and SITC 807 have no public count matrix.

### Hard limit

**Nest entry is not encoded in dissociated RNA.** Bulk and 10x suspensions lose nest vs stroma. Without spatial coordinates or a nest/stroma IHC score linked to the same tumors, TACSTD2/CLDN4 vs CD8 RNA is at best a **purity-confounded infiltration** surrogate.

### Public HIS-NSCLC RNA that actually exists

Three deposits are real HIS-tumor transcriptomes. None are nest assays.

**1. GSE260575 — PBMC-humanized H460, ICI bulk (usable, empty for this question)**

- MHC I/II double-knockout NSG + human PBMC + H460, QuantSeq 3', mapped to hg19, n=12 (Control / Pembro / BJIKT / Combo).
- Paper IHC: lower PD-1+ and LAG-3+ fraction among CD4/CD8. Not nest vs stroma. Not TACSTD2/CLDN4.
- Reanalysis of deposited counts:
  - **TACSTD2 = 0 in every sample.** H460 in this library does not support a TROP2 test.
  - **CD8A** is 0–16 raw counts (~0–2 CPM). **CD8B** is 0 in 11/12 samples. Human T-cell RNA is at the detection floor (PBMC model + tumor-dominated 3' library + human-only align).
  - **CLDN4** is expressed (group-mean CPM Control 40, Pembro 56, BJIKT 63, Combo 93). vs CD8A Spearman ρ = −0.22, approximate p = 0.48, n=12, 5 CD8A zeros.
- Conclusion: this is the only tidy public PBMC-HIS NSCLC ICI bulk series. It cannot test TACSTD2, cannot measure nest, and barely detects CD8.

**2. GSE276724 — CD34-NSG NSCLC PDX bulk (usable with caveats, still not nest)**

- Paper (Kim et al., *Mol Cancer Ther* 2025): JAX **CD34-NSG** mice, YHIM-2004 PDX, denfivontinib ± pembrolizumab. GEO series text says RNA from a humanized-mouse tumor. Supplementary xlsx is **5 PDX × 3** raw counts with **no treatment labels** (GSM titles are only PDX names).
- TACSTD2 and CLDN4 are high in most PDXs. YHIM-2004 has the highest mean CD8A and the lowest mean TACSTD2/CLDN4, but one YHIM-2004 library (`SJ1RNAMFF004644`) is an **epithelial dropout** (KRT8 CPM ≈ 2, EPCAM ≈ 1, PTPRC ≈ 558). That sample inflates any “high immune / low TROP2” story.
- Spearman on CPM (approximate two-sided t p; **not independent mice** — 5 PDX clusters):

  | contrast | n=15 | drop dropout n=14 |
  |---|---|---|
  | TACSTD2 vs CD8A | ρ=−0.26, p=0.35 | ρ=−0.09, p=0.77 |
  | TACSTD2 vs CD8B | ρ=−0.53, p=0.044 | ρ=−0.42, p=0.14 |
  | CLDN4 vs CD8A | ρ=−0.35, p=0.20 | ρ=−0.20, p=0.48 |
  | EPCAM vs CD8A | ρ=−0.35, p=0.20 | ρ=−0.20, p=0.48 |

- Residualizing log1p CPM on log1p EPCAM does **not** create a nest result. TACSTD2 vs CD8B stays nominally negative (p≈0.02) but CD8B is low, n is clustered, and this is still bulk composition.
- Conclusion: direction is compatible with “more epithelial / less CD8 RNA,” which is the purity half of a Bessede-style infiltration story. It is **not** nest entry and should not be over-read.

**3. GSE293914 — PBMC-humanized H1975 scRNA (public, not nest-ready)**

- Kwok et al.: NSG + H1975 + human PBMC, anti-PD-1 ± anti-CCL20. One GEO sample titled with all four arms. Deposit is a single 10x MTX (449 MB) **without** cell-type or treatment barcodes.
- Paper multiplex IF reports more M1-like macrophages, pDC, memory B cells, and spatial segregation of Tregs after combo. That is not CD8-in-nest vs TACSTD2/CLDN4.
- De novo clustering of an unlabeled pool would still not recover nest coordinates. Not run.

### Closed or off-target (do not pretend they are public answers)

- **Yu 2023** (HuNCG CD34 + H226 + pembro ± Shenmai, 10x, 38,010 cells): the only on-paper CD34 HIS NSCLC tumor scRNA with epithelial + T/NK clusters. Data **upon request**, not GEO.
- **SITC 2024 #807**: CD34 NSG-SGM3 + HLA-partial-matched H358 + TROP2 ADC ± sabestomig. Bulk RNA + CellDive IF. **Not deposited.** This is the closest TROP2 + HIS + IO experiment; it is closed.
- **Meraz eLife 98258**: CD34 NSG + A549 KRAS/STK11 mets, NanoString PanCancer Immune n=12. “All data … in the manuscript and source data files.” No GEO. Immune panel; TACSTD2/CLDN4 not a nest readout. They report more CD8b after NPRL2, by design of that therapy.
- **E-MTAB-9891**: PBMC + PC9, osimertinib then RIG-I agonist. Public-ish (ENA), wrong pharmacology for an ICI nest test. Not reanalyzed.
- **GSE217706/GSE217722**: CD34 MISTRG lung ILC/NK or monocytes. No tumor.
- **GSE296572**: “humanized HGF-NSG” = human HGF knock-in, not HIS.
- **GSE194166**: murine KP/KPL scRNA; HIS mentioned only in the paper.
- **GSE285888** etc.: patient PBMCs, not PBMC-humanized mice.

### Honest answer

There is **no** public HIS-NSCLC RNA/scRNA result for TACSTD2/CLDN4 vs CD8 nest entry, because the phenotype was never deposited as a spatial assay.

The two open bulk matrices fail the question in different ways: GSE260575 has no TACSTD2 and no CD8; GSE276724 has both genes but only unlabeled cross-PDX bulk counts, one bad library, and a purity-shaped CD8 signal.

If a later slice needs a HIS functional test, the missing pieces are (i) Yu 2023 counts if the authors release them, (ii) SITC 807 CellDive if it ever lands, or (iii) a new spatial HIS study. Patient ICI spatial cohorts (outside this hunt) are the current place to score nest CD8.

Figures: `figures/fig0_verdict.png`, `fig1_gse260575_cpm.png`, `fig2_gse276724_scatter.png`.

---

## 中文

### 问题

公开的 **CD34 人源化** 或 **PBMC 人源化** NSCLC 免疫治疗 RNA/scRNA 里，有没有 **TACSTD2/CLDN4 对 CD8 巢内进入（nest entry）** 的证据？

巢内进入 = CD8 进上皮肿瘤巢，不是间质 CD8，也不是 bulk `CD8A` 表达。

### 硬限制

解离 RNA **测不到** nest。没有空间坐标或 nest/间质 IHC 评分，就不能回答这个问题。

### 实际打开的公开数据

- **GSE260575**（PBMC + H460 + 帕博利珠单抗 ± BJIKT，n=12）：**TACSTD2 全为 0**；CD8A/CD8B 接近检测下限；CLDN4 有表达，对 CD8A ρ≈−0.22、p≈0.48。不能测 TROP2，不能测 nest。
- **GSE276724**（论文为 CD34-NSG + PDX + denfivontinib/帕博利珠；仓库是 5 个 PDX × 3、无治疗组标签）：TACSTD2/CLDN4 高；对 CD8A 的负相关在去掉一个上皮脱落文库后变成 ρ≈−0.09。有效样本更接近 5 个 PDX，不是 15 只鼠。这是纯度/组成，不是 nest。
- **GSE293914**（PBMC + H1975 的 10x）：公开但 **无细胞注释**；论文多色 IF 不是 TACSTD2/CLDN4–CD8 nest 评分。
- Yu 2023 HuNCG H226 scRNA、SITC 2024 #807 TROP2 ADC CellDive、Meraz NanoString：**没有可用的公开矩阵**来做本题。

### 结论

**没有。** 公开 HIS-NSCLC RNA/scRNA 里不存在 TACSTD2/CLDN4 对 CD8 nest entry 的可引用结果。不要把 GSE276724 的 bulk 负相关写成 nest 证据。
