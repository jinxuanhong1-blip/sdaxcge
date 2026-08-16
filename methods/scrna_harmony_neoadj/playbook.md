# Harmony integration of public neoadjuvant lung scRNA (methods only)

Playbook only. No results, n, ρ, or biology claims.

## Purpose

Integrate **compatible public processed** neoadjuvant NSCLC tumor scRNA matrices that contain **epithelium and immune** cells onto a shared gene space with Harmony (or equivalent). Score malignant-like `TACSTD2` / `CLDN4` against MPR/NMPR with **cohort as a covariate**, and score per-patient malignant programs against T/NK fraction in the **joint embedding**. Split LUAD vs LUSC when the deposited labels allow it.

User A3 (GSE207422 CopyKAT slide) is **taken as given** and is not re-audited here.

## Inclusion

A series is eligible only if all of the following hold:

1. Public processed UMI / author matrix (MTX, TSV UMI, or equivalent). **No FASTQ**, no EGA/GSA raw, no file >2 GB compressed, no cumulative pull that forces raw sequence.
2. Tumor (or tumor-bed) single-cell RNA, not blood / CD45-sorted-only / T-sorted-only.
3. Epithelium **and** immune genes or author labels are present (`TACSTD2`, `CLDN4`, `EPCAM`, `PTPRC`, T/NK markers).
4. Neoadjuvant pathologic response (MPR / NMPR / pCR) is public **or** the series is inventoried and left out with a written reason.

Leave a matrix out when genes are missing, the public object is immune-only, or the endpoint is not neoadjuvant MPR (do not recode RECIST as MPR).

## Shared gene space

1. Take the intersection of deposited gene **symbols**.
2. Require `TACSTD2` and `CLDN4`. If either is absent, drop that series.
3. Keep a lineage/target core panel plus highly variable genes computed on one deposited UMI matrix among the intersection (default cap 800 genes) so Harmony fits in modest RAM.
4. Library size = author `nCount_RNA` when deposited, otherwise the sum of all gene UMIs in that matrix. Scores use `log1p(UMI / nUMI × 10⁴)`.

## Harmony

1. Normalize to CP10k on the shared genes; `log1p`; z-score genes; PCA (default 30 PCs).
2. Run `harmonypy.run_harmony` on PCA with batch = **dataset** (not patient). Write contiguous `X_pca_harmony`.
3. Do **not** treat malignant CNV programs as batch. Harmony is for chemistry/study mixing.
4. Neighbors + Leiden + UMAP on Harmony PCs. If RAM requires it, build the graph on a **sample-stratified subsample** and transfer cluster labels to all cells by kNN in PCA space. Report both total deposited cells and cells in the embedding.
5. Assign Leiden clusters to lineages by mean marker scores (epithelial, T, NK, B, myeloid, endothelial, fibroblast, mast). T/NK = T + NK. Malignant-like = epithelial clusters after dropping a high normal-lung program (`SFTPA2`, `AGER`, `SCGB1A1`, `SCGB3A1`, `TPPP3`) when those genes exist.

Author `major.cell.type` (when present) is a **sensitivity**, not a replacement for the joint-embedding definition used in the T/NK correlation.

## Estimands (pre-specify)

Patient (or post-treatment sample) is the unit. Do not report cell-level p-values for a patient-level label.

| Contrast | Rule |
|---|---|
| MPR | pCR ⊂ MPR. NMPR = non-MPR. RECIST is not MPR. |
| Floor | Drop patients with <20 malignant-like cells from tumor-cell tests; T/NK pairing also requires ≥20 T/NK. List dropouts. |
| TACSTD2 / CLDN4 | Mean `log1p(CP10k)` in malignant-like cells (primary). Also store %UMI>0. |
| (1) MPR | Two-sided MWU (exact Wilcoxon by enumeration when C(n,k)≤20 000). Cohort-adjusted OLS: `score ~ C(MPR) + C(cohort)`. |
| (2) T/NK | Spearman of malignant score vs T/NK fraction. Partial Spearman after residualizing cohort. |
| (3) Histology | LUAD vs LUSC when labeled. Do not recode ASC as LUSC. |

## Pairwise / leave-one-out (preferred over one forced object)

Do not require a single Harmony of every eligible series. Screen pairwise and leave-one-out combinations. Pre-specify KEEP as: malignant-like TACSTD2 NMPR median > MPR median (n≥2/arm) **or** TACSTD2 vs T/NK ρ<0 (n≥4). Extra figures only for KEEP. Report every combo with n/ρ/p, including those that fail. GSE205335 (RECIST) may enter T/NK tests only.

## Honesty

- Report n with every statistic.
- Null + small n = underpowered / inconclusive, not “no association.”
- Do not retune lineage or Harmony to recover a User A3 number.
- Do not concatenate accessions and treat accession as a biological replicate without the cohort covariate.

## Commands

```bash
python methods/scrna_harmony_neoadj/scripts/00_download.py
python methods/scrna_harmony_neoadj/scripts/run_all.py
```
