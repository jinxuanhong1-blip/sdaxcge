# OncoSG LUAD — CLDN4 vs CD8A / ImmuneScore / CD274

**Additive only.** A9 extra-cohort 3/3 catalog (OncoSG LUAD, CPTAC LUAD RNA, GSE31210 LUAD) is **taken as given** and is not re-audited. TACSTD2 vs CD8 / GEP / immune on this same public matrix is already in [PR 139](https://github.com/jinxuanhong1-blip/sdaxcge/pull/139) and is treated as given. This folder only measures **CLDN4**.

Public East-Asian surgical LUAD (Chen et al., *Nat Genet* 2020; cBioPortal [`luad_oncosg_2020`](https://www.cbioportal.org/study/summary?id=luad_oncosg_2020)). Unit is the **public z-score column**. No ICI arm. No slide was re-scored.

## Honest n

| item | n | rule |
|---|---:|---|
| cBioPortal RNA sample list `luad_oncosg_2020_rna_seq_v2_mrna` | 181 | portal description |
| public z-score matrix columns | **169** | `data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt` |
| RNA-list IDs absent from the public matrix | 12 | A008, A114, A122, A136, A139, A147, A184, A302, A435, A484, A489, A507 |
| clinical samples | 305 | `data_clinical_sample.txt` |
| clinical PURITY non-NA | 302 | includes many RNA-absent tumors |
| matrix samples with PURITY | **169** | all 169 columns have published PURITY |
| matrix samples with IMSIG T cells | **169** | published clinical column |
| CLDN4 / CD8A / CD274 finite | **169** | 0 missing z-scores |
| ONCOTREE LUAD on matrix | **169** | all columns are Lung Adenocarcinoma / Tumor |
| ESTIMATE ImmuneScore | **0** | skipped — deposit is z-scores, no raw RSEM |
| ICI / PD-1 labels | **0** | not an ICI-response cohort |

Primary tests use **n=169**. Do not write n=181. The 12 portal-listed RNA IDs are not in the open expression table and cannot enter a CLDN4 correlation.

## One-row table

| dataset | n | CLDN4–CD8A ρ (p) | CLDN4–ImmuneScore ρ (p) | CLDN4–CD274 ρ (p) | CLDN4–CD8A partial \| PURITY (p) | CLDN4–ImmuneScore partial \| PURITY (p) | CLDN4–CD274 partial \| PURITY (p) |
|---|---:|---|---|---|---|---|---|
| OncoSG LUAD public RNA | 169 | -0.416 (1.85e-08) | -0.432 (4.59e-09) | -0.363 (1.28e-06) | -0.285 (0.000181) | -0.308 (4.90e-05) | -0.242 (0.00156) |

ImmuneScore = A1 8-gene T-cell effector mean z (`CD8A GZMA GZMB IFNG EOMES CXCL9 CXCL10 TBX21`; 8/8 present). This is **not** Yoshihara ESTIMATE ImmuneScore. Full numeric rows: `tables/one_row.tsv`, `tables/stats.tsv`.

## ImmuneScore (honest skip of ESTIMATE)

cBioPortal exposes **only z-score** mRNA for this study (`luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores`). Datahub has no `data_mrna_seq_v2_rsem.txt`.

| score | status | why |
|---|---|---|
| ESTIMATE ImmuneScore / StromalScore / TumorPurity | **skipped** | official ssGSEA ranks genes *within a sample* on raw/log expression; ranking across-sample z-scores is a different transform; Affymetrix cosine purity is invalid on z-scores |
| xCell / MCP-counter | **skipped** | need counts or TPM |
| ImmuneScore used here | A1 8-gene mean z | same gene list as PR 139; Spearman is rank-based |
| published immune column | IMSIG T cells | deposited on `data_clinical_sample.txt`; not recomputed from genes |

Do not quote an ESTIMATE ImmuneScore for OncoSG from these public files.

**Covariate:** published sample-level `PURITY` (OncoSG clinical attribute). This is **not** TCGA ABSOLUTE and **not** ESTIMATE cosine purity.

**Statistic:** Spearman on deposited z-scores. Partial = Pearson of average-rank residuals after published PURITY; two-sided t, df = n−3; Fisher-z 95% CI with SE = 1/√(n−4). Same estimator as PR 139.

## Main Spearman (n=169)

Partial residualises on published PURITY.

| pair | n | ρ | p | partial ρ \| PURITY | partial p |
|---|---:|---:|---:|---:|---:|
| CLDN4 vs CD8A | 169 | -0.416 | 1.85e-08 | -0.285 | 0.000181 |
| CLDN4 vs ImmuneScore (A1 8-gene) | 169 | -0.432 | 4.59e-09 | -0.308 | 4.90e-05 |
| CLDN4 vs CD274 | 169 | -0.363 | 1.28e-06 | -0.242 | 0.00156 |
| CLDN4 vs IMSIG T cells (published) | 169 | -0.526 | 2.04e-13 | -0.414 | 2.34e-08 |
| CLDN4 vs TACSTD2 | 169 | 0.505 | 2.65e-12 | 0.462 | 2.79e-10 |
| CLDN4 vs PURITY | 169 | 0.359 | 1.64e-06 | NA | NA |
| CLDN4 vs TJ (no CLDN4) | 169 | 0.547 | 1.45e-14 | 0.484 | 2.87e-11 |
| CD8A vs ImmuneScore (positive control) | 169 | 0.910 | 8.53e-66 | NA | NA |
| CD274 vs ImmuneScore | 169 | 0.563 | 1.62e-15 | 0.434 | 4.33e-09 |
| TACSTD2 vs CD8A (companion; PR 139) | 169 | -0.380 | 3.38e-07 | -0.309 | 4.69e-05 |

CLDN4 vs PURITY Spearman ρ = 0.359 (p = 1.64e-06). Partialling purity shrinks the immune |ρ| but does not flip the CD8A / ImmuneScore / IMSIG T-cell signs.

CD8A vs ImmuneScore ρ = 0.910 is a positive control and is partly by construction (CD8A is one of the eight ImmuneScore genes).

CD274 (PD-L1 transcript) is **not** a CD8 / ImmuneScore substitute. On this matrix CLDN4–CD274 is weaker than CLDN4–CD8A.

## Extra scatter

Unadjusted scatters for the three requested pairs, plus extra residual and companion panels:

- `figures/fig1_cldn4_vs_cd8a.png`
- `figures/fig2_cldn4_vs_immunescore.png`
- `figures/fig3_cldn4_vs_cd274.png`
- `figures/fig4_spearman_forest.png`
- `figures/fig5_correlation_heatmap.png`
- `figures/fig6_extra_residual_scatter.png` — rank residuals after PURITY (CD8A / ImmuneScore / CD274)
- `figures/fig7_extra_imsig_tacstd2_scatter.png` — CLDN4 vs IMSIG T cells and vs TACSTD2
- `figures/fig8_extra_purity_control_scatter.png` — CLDN4 vs PURITY; CD8A vs ImmuneScore control

## What is not done

- No re-audit of the A9 3/3 catalog or of PR 139 TACSTD2 numbers.
- No ESTIMATE / xCell / MCP on z-scores.
- No ICI ORR / PFS model (labels are not deposited).
- No imputation of the 12 portal-listed RNA IDs that lack a public z-score column.
- This is East-Asian surgical LUAD, not an ICI-response cohort. A negative bulk correlation is not evidence that CLDN4-high tumors fail checkpoint blockade.

## Files

- `analyze.py` — Datahub download, complete-case n, Spearman, extra scatters
- `tables/one_row.tsv`, `stats.tsv`, `per_sample.tsv`, `coverage.tsv`, `gene_coverage.tsv`, `summary.json`
- `figures/fig1_cldn4_vs_cd8a.png` … `fig8_extra_purity_control_scatter.png`

```bash
python3 methods/oncosg_cldn4_immune/analyze.py
```
