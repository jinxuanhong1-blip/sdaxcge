# Public lung scRNA TLS / B-cell vs malignant TACSTD2 / CLDN4: methods playbook
# 公开肺肿瘤单细胞：TLS / B 细胞 vs 恶性 TACSTD2 / CLDN4 方法手册

## 1. Scope and estimand / 范围与目标效应

**Additive only.** User A6 / B6 immune-cold (T-cell / spatial neighborhood) is taken as given and is not re-cut here. This slice asks a different question: whether the same malignant TACSTD2 / CLDN4 axis also tracks **B-cell and TLS-like** readouts in public lung-tumor scRNA.

**仅增量。** 用户 A6 / B6 免疫冷（T 细胞 / 空间邻域）视为既定，本切片不重做。本切片只问：恶性 TACSTD2 / CLDN4 是否同时对应 **B 细胞 / TLS 样** 读数。

**Primary estimand.** Within each independent cohort, the patient-level Spearman correlation between (i) malignant (or, if malignancy is not labeled, epithelial / author gate) TACSTD2 or CLDN4 and (ii) a B / TLS-like score. Cohort-level Fisher-z values are then combined with a DerSimonian–Laird random-effects meta-analysis. Matrices are **not** concatenated across GEO accessions.

**主目标效应。** 每个独立队列内，患者水平 Spearman：恶性（若无恶性标签则上皮 / 作者 gate）TACSTD2 或 CLDN4 vs B / TLS 样分数。队列 Fisher-z 再做 DerSimonian–Laird 随机效应荟萃。禁止把各 GEO 表达矩阵拼成一张表。

**Experimental unit = patient**, not cell and not 10x lane. Multi-sample patients are collapsed (tumor tissue only).

**实验单位 = 患者**，不是细胞、不是 lane。多样本患者在肿瘤组织内合并。

**Effect direction.** Negative ρ means higher malignant TACSTD2 / CLDN4 goes with lower B / TLS-like scores. Positive always means the opposite. Direction is pre-specified and is not flipped after seeing results.

**效应方向。** 负 ρ = 恶性 TACSTD2 / CLDN4 高伴随 B / TLS 低。方向预先指定，见结果后不翻转。

This slice does **not** test ICI response, MPR, ORR, or spatial follicles. Treatment-naive and neoadjuvant cohorts are both eligible; labels are preserved, not recoded onto one scheme.

本切片 **不** 检验 ICI 反应、MPR、ORR 或组织学滤泡。初治与新辅助队列都可纳入；标签保持原样。

## 2. Eligible public series / 合格公开系列

Named series (public processed files only; no EGA / dbGaP FASTQ):

- GSE131907 (Kim 2020 LUAD atlas; author cell types; B cells present)
- GSE253013 (Sze/Xiang 2024 treatment-naive LUAD; mixed TME)
- GSE207422 (Hu 2023 neoadjuvant PD-1 + chemo; mixed TME)
- GSE241934 (NEOTIDE / real-world neoadjuvant IO; IIT and RWC are non-overlapping patients and are two strata)
- GSE154826 (Leader 2021 CITE-seq) **only if** epithelium or an author gate fraction can be scored. CD45-bead libraries are sort-biased; that is recorded, not hidden.
- Other public lung-tumor scRNA with B cells and a processed matrix that can be streamed (example: GSE148071 Wu 2021 advanced NSCLC).

A series is skipped (honest row, no substitution) if: no public processed matrix, file > ~10 GB and cannot be streamed, epithelium-only with no B/TLS genes and no second immune library, or n patients < 5 after floors.

无法纳入时写诚实 skip，不替换登录号。

## 3. Lineage and malignant calls / 谱系与恶性

Prefer **author labels**. If barcodes are unlabeled, lineage = argmax of mean log1p(UMI) marker scores (epithelial / T / NK / B / plasma / myeloid / fibroblast / endothelial). Unassigned if best < 0.15 or (best − second) < 0.05.

**Malignant-like** (when authors do not call it): epithelial AND low normal-lung markers (SFTPA / AGER / SCGB / TPPP3 / FOXJ1). This is a marker proxy, not CopyKAT. GSE154826 uses the authors’ `epi_endo_fibro_doublet` gate and is labeled `gate`, not malignant.

B = author B / `B lymphocytes` / marker B. Plasma is separate. `frac_B` uses B only; `frac_B_plasma` is secondary.

## 4. Scores / 分数

Report all three malignant / epithelial TACSTD2 and CLDN4 metrics; do not pick the “best”:

| Metric | Formula |
|---|---|
| `mean_log1p_cp10k` | mean of `log1p(UMI / total_UMI × 10⁴)` |
| `pct_pos` | % cells with UMI > 0 (no 1% / 10% cutoff) |
| `pseudobulk_cpm` | `sum(gene UMI) / sum(total UMI) × 10⁶` |

**Primary x** = malignant `mean_log1p_cp10k` (epithelial or gate if malignant is unavailable).

**Primary y**

1. `frac_B` = n_B / n_cells in the same tumor object (sort bias recorded).
2. `tls12_z` = mean of per-gene within-cohort z-scores of patient-mean log1p(CP10K) for the 12 chemokines. Report k/12 present.
3. Secondary: CXCL13 mean / %pos (all cells); MS4A1 mean; `frac_B_plasma`; CXCL13+ T-cell fraction.

No both-high TACSTD2∩CLDN4 gate. Optional epithelial module = median(TACSTD2, CLDN4) is exploratory only.

## 5. Floors and tests / 门槛与检验

- Eligible patient: tumor tissue, ≥20 malignant / epithelial / gate cells.
- `frac_B` allows zero B cells (zero is informative). B-intrinsic scores require ≥10 B cells.
- Spearman needs n ≥ 5 finite pairs. Underpowered → “underpowered; inconclusive,” never “null = no association.”
- Always report **n / ρ / p**. BH-FDR only inside the pre-specified primary list (4 contrasts × eligible cohorts).
- Meta: Fisher z, SE = 1/√(n−3), DerSimonian–Laird τ², random-effects pooled ρ = tanh(z). Forest one row per independent patient set.
- Extra table: TLS at single-cell (CXCL13+ and MS4A1+ rates by lineage). Not a follicle call.

## 6. What this cannot test / 本切片不能检验的内容

- A6 / B6 immune-cold T-cell claim (taken as given).
- Histologic TLS / spatial B-cell aggregates (dissociated 10x).
- ICI benefit as a function of TLS.
- Protein TROP2 (RNA only; GSE154826 has no TACSTD2 ADT).
- Private / EGA FASTQ.

## 7. Outputs / 产出

Scripts write only to `results/scrna_tls_meta/`. Playbook stays methods-only (no numbers).
