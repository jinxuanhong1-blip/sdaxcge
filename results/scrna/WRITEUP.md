# Malignant TACSTD2 / CLDN4 vs T/NK and ICI response (scRNA)

**Verdict:** TACSTD2 is **epithelial / malignant-restricted** (real in GSE207422, GSE205335, and the GSE131907 atlas). It is **not** significantly higher in ICI non-responders, and it is **not** anti-correlated with T/NK or CD8 at the claimed ρ ≈ −0.4 to −0.5. GSE207422 (MPR) is underpowered after the <10 malignant-cell rule (MPR n=2). GSE205335 (RECIST) is a null response test; TACSTD2 vs T/NK is weakly **positive** (ρ = +0.28 to +0.46), opposite the claim. GSE131907 (no ICI labels) sample-level epi TACSTD2 vs T/NK is also null (tLung n=11 ρ=+0.09; tumor sites n=36 ρ=+0.17).

## 中文摘要

在两个公开人 NSCLC ICI 单细胞队列里，**TACSTD2 明确限制在恶性/上皮细胞**（相对 T/NK：GSE207422 n=9，Wilcoxon p=0.0039；GSE205335 n=22，p=4.8×10⁻⁷）。  
**GSE207422（新辅助 PD-1+化疗，MPR）**：术后样本中恶性样细胞 TACSTD2，MPR n=2 vs NMPR n=7，Mann–Whitney p=0.50；与 T/NK 分数 Spearman n=9，ρ=+0.17，p=0.67。不支持 ρ≈−0.4～−0.5。MPR 样本恶性细胞很少，按作者 <10 细胞规则会丢掉 2/4 个 MPR，检验力不足。  
**GSE205335（姑息 ICI，仅 RECIST，不是 MPR）**：恶性 TACSTD2，PR n=6 vs SD/PD n=10，p=0.96；与 T/NK n=22，ρ=+0.28，p=0.20。CLDN4 与 CD8 在 R+NR（n=16）有负向趋势（ρ=−0.43，p=0.097），未过 0.05。  
**GSE131907（LUAD 图谱，无 ICI 标签）**：上皮 vs T/NK TACSTD2 %pos 74.5 vs 1.6；样本水平 tLung n=11 ρ=+0.09，p=0.79。  
GSE271689 GEO 只有 RTS 探针 DCC，没有基因注释和 OS；GSE154826 无 ICI 疗效标签。数字全部从 GEO 处理后的矩阵当场计算。

## Question

In **malignant / epithelial cells only**, is TACSTD2 (Trop-2) or CLDN4 (i) higher in ICI non-responders and (ii) negatively correlated with T/NK, CD8, or TLS-proxy infiltration at Spearman ρ ≈ −0.4 to −0.5?

## Datasets (processed GEO only)

| Accession | Design | ICI label | Cells used | Malignant definition |
|---|---|---|---|---|
| **GSE207422** | Neoadjuvant PD-1 + chemo, stage IIIA NSCLC (Hu et al. *Genome Med* 2023, PMID 36869384) | **MPR / pCR / NMPR** (+ RECIST) | 92,330 (all UMI ≥ 200) | Marker proxy: epithelial-lineage **and** normal-lung score ≤ 75th percentile of epithelial cells. GEO has **no** CopyKAT labels. |
| **GSE205335** | Palliative lung ICI atlas | **RECIST** PR / SD / PD / NE. **No MPR field** | 96,505 | Authors’ published `lineage.sub == Malignant cells` (28,512 cells). Normal tissues excluded from tests. |
| **GSE131907** | LUAD atlas (Kim et al. *Nat Commun* 2020, PMID 32385277) | **None** (treatment-naive / mixed; not an ICI trial) | 208,506 | Author `Cell_type == Epithelial cells` (36,467) and `Cell_subtype == Malignant cells` (24,784; mets/PE). |

Skipped as instructed: T-sorted-only (GSE176022, GSE99254); files >2 GB (GSE131907 log2TPM txt 2.9 GB; LuCA h5ad). **GSE154826** LCAM annots have no ICI/RECIST/MPR field (CD45-enriched LUAD/LUSC). **GSE271689** GEO SOFT has `treatment=Immunotherapy` and CK/CD45/CD68 segments, but DCC rows are **RTS probe IDs** with no gene map and **no OS/PFS** — not analyzed.

## Methods (sample / patient is the unit)

- Expression: mean **log1p CP10K** of TACSTD2 or CLDN4 inside malignant / malignant-like cells.
- Immune: T/NK, CD8, B/Plasma fractions of all cells in that sample/patient. B/Plasma is the TLS proxy; CXCL13/CCL19/CCL21/CXCL9/CXCL10 mean is a second TLS score (GSE207422).
- GSE207422 primary: **12 post-treatment surgical samples**. pCR counted as MPR. Samples with **<10 malignant-like cells** dropped from gene means (authors dropped one NMPR with <10 malignant cells). Pre-treatment biopsies reported only as sensitivity.
- GSE205335: **patient** unit; ≥20 malignant cells; PR = R, SD/PD = NR; NE kept only in “all evaluable” Spearman.
- Tests: two-sided Mann–Whitney U and Spearman. Cell-level tests are exploratory (pseudoreplication).
- GSE207422 lineages are a marker argmax, **not** CopyKAT. That is a real limitation.
- GSE131907: author annotations; metric = mean **log1p(raw UMI)** and %pos (library-size CP10K not computed; panel-only extract). Not an ICI-response test.

## GSE207422 results (MPR)

Lineage (marker argmax): T/NK 38,514; myeloid 25,383; epithelial 13,043; B/Plasma 11,327; endothelial 1,956; fibroblast 1,103; mast 1,004. Malignant-like = 9,782.

Post-tx with ≥10 malignant-like cells: **n=9 (MPR=2, NMPR=7)**. Dropped MPR: P11 (1 cell), P14 (8 cells). That is expected after a good pathologic response and **systematically removes MPR from the TACSTD2 test**.

| Test | n | Stat | p |
|---|---|---|---|
| Malignant TACSTD2, NMPR vs MPR | 7+2 | U=4.0; med 1.10 vs 0.43 | **0.50** |
| Malignant CLDN4, NMPR vs MPR | 7+2 | U=5.0; med 0.73 vs 0.49 | **0.67** |
| TACSTD2 vs T/NK fraction | 9 | ρ = **+0.167** | 0.67 |
| TACSTD2 vs CD8-like fraction | 9 | ρ = −0.017 | 0.97 |
| TACSTD2 vs B/Plasma (TLS proxy) | 9 | ρ = +0.067 | 0.86 |
| TACSTD2 vs TLS chemokine score | 9 | ρ = +0.25 | 0.52 |
| TACSTD2 vs residual tumor fraction | 9 | ρ = +0.41 | 0.27 |
| T/NK fraction, NMPR vs MPR (all 12 post) | 8+4 | U=19 | 0.68 |
| CD8-like fraction, NMPR vs MPR | 8+4 | U=23 | 0.28 |
| Paired malignant vs T/NK TACSTD2 | 9 | W=0; med 0.78 vs 0.009 | **0.0039** |

Sensitivity (min 1 malignant-like cell, keeps remnant MPR): TACSTD2 MPR n=4 vs NMPR n=8, p=0.15; vs T/NK n=12, ρ=+0.20, p=0.54. All-epithelial (no normal-lung filter): response p=0.68; vs T/NK ρ=+0.13, p=0.70.

**Does not support** “NMPR > MPR” or “ρ ≈ −0.4 to −0.5 vs T/NK”. Direction of the two remaining MPR means is lower TACSTD2, but n=2.

## GSE205335 results (RECIST, not MPR)

Tumor tissues; patients with ≥20 malignant cells: **n=22** (R=6, NR=10, NE=6). Four patients had ~0 malignant cells captured and were excluded (asymmetric dropout: more PR than PD).

| Test | n | Stat | p |
|---|---|---|---|
| Malignant TACSTD2, NR vs R | 10+6 | U=29; med 0.91 vs 1.10 | **0.96** |
| Malignant CLDN4, NR vs R | 10+6 | U=36; med 0.95 vs 1.40 | **0.56** |
| TACSTD2 vs T/NK (all evaluable) | 22 | ρ = **+0.284** | 0.20 |
| TACSTD2 vs T/NK (R+NR) | 16 | ρ = **+0.456** | 0.076 |
| TACSTD2 vs CD8 (all evaluable) | 22 | ρ = **+0.409** | 0.059 |
| TACSTD2 vs B/Plasma | 22 | ρ = +0.030 | 0.89 |
| CLDN4 vs T/NK (R+NR) | 16 | ρ = −0.347 | 0.19 |
| CLDN4 vs CD8 (R+NR) | 16 | ρ = −0.429 | 0.097 |
| NSCLC ADC+SQ TACSTD2, NR vs R | 8+4 | U=22 | 0.37 |
| NSCLC ADC+SQ TACSTD2 vs T/NK | 12 | ρ = −0.021 | 0.95 |
| Paired malignant vs T/NK TACSTD2 | 22 | W=0; med 0.93 vs 0.021 | **4.8×10⁻⁷** |

RECIST is **not** substituted for MPR. Cohort includes SCLC and NUT; NSCLC-only sensitivity is also null.

## What is supported

1. **Compartment:** TACSTD2 (and CLDN4) live in malignant / epithelial cells, not T/NK. This is the only robust result.
2. **Response and infiltration claims are not supported** at sample/patient level with honest n. GSE207422 is too small after the malignant-cell gate. GSE205335 TACSTD2–T/NK/CD8 trends are **positive**, not negative.
3. CLDN4 vs CD8 in GSE205335 R+NR (ρ=−0.43, p=0.097) is the only infiltration trend in the claimed direction; it is not significant and was not pre-specified as the primary.
4. GSE131907 confirms compartment restriction in a large annotated atlas and again fails to show sample-level TACSTD2–T/NK anti-correlation.

## GSE131907 results (atlas; not ICI)

208,506 cells, author labels. Epithelial 36,467; malignant subtype 24,784 (almost all mets/PE, **none in tLung**); T/NK 91,227. EPCAM %pos 82.5 (epi) vs 1.5 (T/NK); PTPRC 3.2 vs 67.7.

| Test | n | Stat | p |
|---|---|---|---|
| Epithelial vs T/NK TACSTD2 (cells) | 36,467+91,227 | %pos 74.5 vs 1.6; med log1p UMI 1.10 vs 0 | **<1e-300** (exploratory) |
| tLung epi vs T/NK TACSTD2 (cells) | 7,270+19,591 | %pos 86.0 vs 3.2 | **<1e-300** (exploratory) |
| Paired tumor-site epi vs T/NK TACSTD2 | 36 samples | W=1; med 1.35 vs 0.02 | **2.7×10⁻⁷** |
| tLung epi TACSTD2 vs T/NK fraction | 11 | ρ = **+0.091** | 0.79 |
| Tumor-site epi TACSTD2 vs T/NK | 36 | ρ = **+0.168** | 0.33 |
| Tumor-site malignant TACSTD2 vs T/NK | 21 | ρ = +0.082 | 0.72 |
| Tumor-site malignant CLDN4 vs T/NK | 21 | ρ = −0.396 | 0.076 |

No ICI / MPR / RECIST labels. Does not support ρ ≈ −0.4 vs T/NK at the sample level.

## Figures

- `gse207422_mal_tacstd2_cldn4_vs_mpr.png`
- `gse207422_mal_tacstd2_vs_immune.png`
- `gse207422_mal_cldn4_vs_immune.png`
- `gse207422_tacstd2_compartment.png`
- `gse207422_lineage_counts.png`
- `gse205335_mal_tacstd2_cldn4_vs_recist.png`
- `gse205335_mal_tacstd2_vs_immune.png`
- `gse205335_mal_cldn4_vs_immune.png`
- `gse131907_compartment.png`
- `gse131907_epi_tacstd2_vs_tnk.png`
- `gse131907_tacstd2_by_celltype.png`

Tables: `gse207422_sample_table.tsv`, `gse205335_patient_table.tsv`, `gse131907_sample_table.tsv`, `stats.tsv`.

## Reproducibility

```bash
bash scripts/scrna/00_download.sh
python3 scripts/scrna/01_extract_gse207422.py
python3 scripts/scrna/02_extract_gse205335.py
python3 scripts/scrna/03_analyze.py
python3 scripts/scrna/04_gse131907_atlas.py
```

Raw matrices stay in `/tmp/scrna_data/` (not committed). Gene panel: `scripts/scrna/gene_panel.tsv`.
