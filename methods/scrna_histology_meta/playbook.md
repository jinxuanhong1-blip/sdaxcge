# Split-and-pool public lung tumor scRNA by histology
# 按组织学拆分并合并公开肺癌肿瘤 scRNA

> **Additive only · 仅加法。** New folder. Does not edit sibling scRNA or TLS analyses.
> **Public only · 仅公开。** GEO processed matrices. No EGA / dbGaP / DAC download.

## 0. Why this exists · 为何做这一步

Bulk TLS work in this repository (`notes/opus_tls`) found that the TACSTD2–TLS/B/CD8 package after purity **holds in TCGA-LUSC and does not hold as TLS in TCGA-LUAD**. That is a histology interaction at bulk RNA.

This folder asks the same split at **single-cell**, with a larger public n than any one ICI scRNA series: **malignant / epithelial TACSTD2 or CLDN4 vs T/NK fraction**, computed **inside LUAD** and **inside LUSC**, then meta **within histology**. LUAD and LUSC are never pooled.

中文：这是 TLS 组织学交互的单细胞、更大 n 版本。LUAD 与 LUSC **分开算、分开合并**，禁止混池。

## 1. Estimand · 估计目标

| Item | Choice |
| --- | --- |
| Unit | Sample or patient. Cells are not replicates of a patient-level correlation. |
| Exposure | Mean TACSTD2 or CLDN4 in **malignant** cells when the author labeled them; else **epithelial**. Prefer log1p(CP10k) or the sibling table’s published mean. |
| Outcome | T/NK fraction of all cells in that sample/patient |
| Stratification | Author histology only: LUAD / Adeno / ADC vs LUSC / Squamous / SQ |
| Held out | ASC, SCLC, NUT, “NSCLC” without a LUAD/LUSC field |
| Association | Two-sided Spearman ρ |
| Meta | Fisher z, DerSimonian–Laird random effects **and** fixed effect, one model per gene × histology |
| Primary inclusion | n ≥ 6 independent samples, finite ρ, no overlapping patient sets |

Do **not** mix TACSTD2 with CLDN4. Do **not** mix T/NK fraction with a TLS signature or a CD8 module. Those are different estimands.

中文：样本/患者为单元；恶性（或上皮）基因均值对 T/NK 比例；只按作者组织学分层；ASC/SCLC/未标注 NSCLC 不进 LUAD 或 LUSC 池。

## 2. Series gate · 系列门控

A series enters the **primary** LUAD or LUSC pool only if all of the following are true:

1. Public processed counts exist (GEO suppl / author table already extracted in this repo).
2. TACSTD2 and CLDN4 are present in epithelium (not a CD45-only object).
3. A **LUAD and/or LUSC** label is on GEO, in the author metadata table, or in a harvested sibling table. Paper-only supplements that were not parsed are catalogued, not guessed.
4. After the cell-count gate, that histology arm has n ≥ 6.

Named series (user list + others surveyed):

| Series | Histology label | Primary? |
| --- | --- | --- |
| GSE207422 | Pathology Adeno / Squamous | No (LUAD n=4, LUSC n=5) |
| GSE241934 | Histology LUAD / ASC | Yes, LUAD (IIT n=10, RWC n=28) |
| GSE131907 | LUAD atlas | Yes, tLung n=11 |
| GSE253013 | LUAD | Yes, n=9 |
| GSE205335 | ADC / SQ / SCLC / NUT | Yes, LUAD n=14; LUSC n=3 not pooled |
| GSE291670 | GEO: NSCLC only | No — no LUAD/LUSC field |
| GSE148071 | Paper supplement, not on SOFT | Catalogued |
| GSE154826 | LUAD/LUSC but CD45-enriched | Not usable for epithelial TACSTD2 |
| LuCA / HLCA integrated | Mixed | Files too large; skipped |

## 3. Cell and sample gates · 细胞与样本门控

- Malignant definition is **author** when present (GSE205335 `Malignant cells`, GSE131907 `Malignant cells` / `Epithelial cells`, GSE241934 `Epi`). GSE207422 and GSE253013 use the sibling marker malignant-like rule; that limitation is inherited, not hidden.
- Drop samples with too few epithelial/malignant cells (sibling defaults: GSE207422 ≥10, GSE205335 ≥20, GSE131907 ≥20, GSE253013 ≥10, GSE241934 ≥10 epi and ≥20 T/NK).
- GSE207422 primary timing is **post-treatment surgery**, matching the sibling ICI scRNA slice.
- GSE131907 primary site is **tLung** (primary tumor). Mets / PE malignant-subtype rows are sensitivity only and are not in the primary meta (overlapping patients).
- GSE241934 IIT (EGFR-mut NEOTIDE) and RWC (WT LUAD/ASC) are **separate cohorts** (independent patients, different designs). A combined LUAD row is sensitivity only.

## 4. Meta · 合并

For each gene × histology:

1. Convert ρ to Fisher z = artanh(ρ), SE = 1/√(n−3).
2. Fixed-effect inverse-variance pool.
3. DerSimonian–Laird τ²; random-effects pool.
4. Back-transform with tanh. Report k, n_total, ρ_RE, 95% CI, p, I², Q.

If k = 0, write that. Do not fill LUSC with LUAD. Do not borrow a bulk ρ.

## 5. What not to claim · 不可声称

- Do not claim a LUSC scRNA meta. There is **no** public series here with LUSC n ≥ 6 and epithelial TACSTD2/CLDN4 vs T/NK.
- Do not claim the bulk LUSC TLS package is “replicated at single-cell.” The scRNA LUAD pool is null-to-weak (TACSTD2 RE ρ = −0.21, p = 0.11). That is compatible with the bulk LUAD TLS-null, not a LUSC confirmation.
- Do not treat GSE253013 ρ = −0.72 (n=9) as the series-level truth. It is one small cohort; the other four LUAD cohorts are near zero.
- Do not call ASC “LUAD” or “LUSC.”
- Do not assign GSE291670 a histology from the paper title.

## 6. Reproduction · 复现

See [README.md](README.md). Harvested sibling tables are committed so the Spearman / meta step does not re-download GSE207422 / GSE205335 / GSE131907 / GSE253013. GSE241934 MTX is public GEO; the extracted per-sample table is committed.
