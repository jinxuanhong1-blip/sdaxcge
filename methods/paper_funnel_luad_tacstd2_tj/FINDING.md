# PAPER FUNNEL — Bulk LUAD TACSTD2-high vs low (OncoSG + GSE31210)

**Rule:** never fabricate. Every number below was computed in `analyze.py`.
PPT claims A1 / A8 are corroborated or not from these numbers; empty/null is reported as such.

## Cohorts

| cohort | n | Q4 vs Q1 | matrix | ImmuneScore | purity covariate |
|---|---:|---|---|---|---|
| OncoSG_LUAD | **169** | 43 vs 43 | z-score RSEM (all-sample ref) | A1_8gene_mean_z | published PURITY |
| GSE31210_LUAD | **226** | 57 vs 57 | Affymetrix GPL570 MAS5 series matrix | ESTIMATE_Immune141_ssGSEA | ESTIMATEScore (rank residual; cosine TumorPurity often outside [0,1]) |

OncoSG portal RNA list is 181; public z-score matrix has **169** columns. Do not write n=181. ESTIMATE is **skipped** on OncoSG (z-scores only).

## PPT corroboration

| cohort | claim | call | detail |
|---|---|---|---|
| OncoSG_LUAD | A1 | **PASS** | TACSTD2–CD8A unadj ρ=-0.380 p=3.38e-07; partial ρ=-0.309 p=4.69e-05. TACSTD2–ImmuneScore unadj ρ=-0.387 p=1.98e-07; partial ρ=-0.318 p=2.70e-05. |
| OncoSG_LUAD | A8 | **PARTIAL** | KEGG_TJ NES=-1.00 FDR_prim=0.318; HALLMARK_APICAL_JUNCTION NES=-1.48 FDR_prim=0.0045; GOBP_KERATINIZATION NES=1.92 FDR_prim=0.002; HALLMARK_EMT NES=-2.33 FDR_prim=0.002. |
| OncoSG_LUAD | CLDN4_immune_companion | **PASS** | CLDN4–CD8A unadj ρ=-0.416 p=1.85e-08; partial ρ=-0.285 p=0.000181. CLDN4–ImmuneScore unadj ρ=-0.432 p=4.59e-09; partial ρ=-0.308 p=4.90e-05. |
| GSE31210_LUAD | A1 | **PARTIAL** | TACSTD2–CD8A unadj ρ=-0.289 p=9.76e-06; partial ρ=-0.107 p=0.108. TACSTD2–ImmuneScore unadj ρ=-0.306 p=2.75e-06; partial ρ=-0.023 p=0.735. |
| GSE31210_LUAD | A8 | **PASS** | KEGG_TJ NES=2.00 FDR_prim=0.00183; HALLMARK_APICAL_JUNCTION NES=-1.16 FDR_prim=0.0989; GOBP_KERATINIZATION NES=1.91 FDR_prim=0.00183; HALLMARK_EMT NES=-2.15 FDR_prim=0.00183. |
| GSE31210_LUAD | CLDN4_immune_companion | **PARTIAL** | CLDN4–CD8A unadj ρ=-0.341 p=1.45e-07; partial ρ=-0.127 p=0.0574. CLDN4–ImmuneScore unadj ρ=-0.366 p=1.50e-08; partial ρ=-0.031 p=0.639. |

### Call rules (locked)

- **A1 PASS:** TACSTD2 vs CD8A **and** vs ImmuneScore both have partial ρ < 0 and partial p < 0.05.
- **A1 PARTIAL:** unadjusted Spearmans are negative and significant for at least one endpoint, but purity-adjusted tests do not both pass.
- **A8 PASS:** KEGG TJ or Hallmark apical junction UP (NES>0, primary FDR<0.05) **and** GOBP keratinization UP **and** Hallmark EMT DOWN.
- **A8 PARTIAL:** keratin and/or TJ/apical UP at primary FDR<0.05, but the full triad (TJ-or-apical UP + keratin UP + Hallmark EMT DOWN) does not hold (includes cases where Hallmark apical junction is DOWN, or EMT fails).
- **CLDN4 companion PASS:** CLDN4 vs CD8A and vs ImmuneScore both negative after purity; PARTIAL = unadjusted only.

## TACSTD2 / CLDN4 vs ImmuneScore / CD8A

| cohort | TACSTD2–CD8A ρ (p) | TACSTD2–Immune ρ (p) | TACSTD2–CD8A partial (p) | TACSTD2–Immune partial (p) | CLDN4–CD8A ρ (p) | CLDN4–Immune ρ (p) |
|---|---|---|---|---|---|---|
| OncoSG_LUAD | -0.380 (3.38e-07) | -0.387 (1.98e-07) | -0.309 (4.69e-05) | -0.318 (2.70e-05) | -0.416 (1.85e-08) | -0.432 (4.59e-09) |
| GSE31210_LUAD | -0.289 (9.76e-06) | -0.306 (2.75e-06) | -0.107 (0.108) | -0.023 (0.735) | -0.341 (1.45e-07) | -0.366 (1.50e-08) |

ImmuneScore definitions differ by cohort (A1 8-gene mean z on OncoSG; ESTIMATE Immune141 on GSE31210). Do not pool the two ImmuneScore ρ values.

## GSEA primary sets (TACSTD2 Q4 vs Q1, Welch prerank)

| set | want | OncoSG | GSE31210 |
|---|---|---|---|
| KEGG_TIGHT_JUNCTION | UP | NES=-1.00 FDR=0.318 | NES=2.00 FDR=0.00183 |
| HALLMARK_APICAL_JUNCTION | UP | NES=-1.48 FDR=0.0045 | NES=-1.16 FDR=0.0989 |
| GOBP_KERATINIZATION | UP | NES=1.92 FDR=0.002 | NES=1.91 FDR=0.00183 |
| GOBP_TIGHT_JUNCTION_ORGANIZATION | UP | NES=1.07 FDR=0.17 | NES=1.79 FDR=0.00244 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | DOWN | NES=-2.33 FDR=0.002 | NES=-2.15 FDR=0.00183 |
| GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION | DOWN | NES=-1.91 FDR=0.002 | NES=-1.54 FDR=0.0132 |

BH-FDR is within the 12 primary sets per cohort. Full table: `tables/gsea_primary.tsv`.

## Focal TJ genes (Welch delta, Q4 − Q1)

Delta units differ: OncoSG = z-score mean difference; GSE31210 = MAS5 intensity mean difference (not log2FC). Do not compare absolute delta magnitudes across cohorts.

| gene | OncoSG delta / FDR | GSE31210 delta / FDR |
|---|---|---|
| CLDN1 | +0.977 / 0.000169 | +2277.132 / 0.000111 |
| CLDN4 | +1.275 / 2.13e-06 | +704.175 / 1.78e-09 |
| CLDN7 | +1.190 / 9.33e-06 | +196.745 / 0.0201 |
| F11R | +0.351 / 0.228 | +2911.496 / 1.04e-05 |
| PARD3 | +0.344 / 0.198 | +207.689 / 0.0832 |
| OCLN | +0.334 / 0.255 | +2746.699 / 6.22e-05 |
| TJP1 | -0.429 / 0.105 | +407.066 / 0.522 |
| CGN | +0.701 / 0.0148 | +824.422 / 4.06e-06 |
| CRB3 | +1.180 / 1.82e-05 | +117.139 / 0.00141 |

## What this is / is not

- East-Asian surgical LUAD (OncoSG + GSE31210). **Not** an ICI-response cohort.
- Does **not** re-audit TCGA-LUAD A8 (see PR 117): LUAD Hallmark EMT can be UP.
- Does **not** invent ESTIMATE on OncoSG z-scores.
- Does **not** merge private 8KL scRNA.
- PPT numbers are corroborated or contradicted by **these** public matrices only.

## Reproduce

```bash
pip install -r methods/paper_funnel_luad_tacstd2_tj/requirements.txt
python3 methods/paper_funnel_luad_tacstd2_tj/analyze.py
```

## Files

- `tables/immune_correlations.tsv`, `gsea_primary.tsv`, `gsea_all.tsv`
- `tables/deg_focal.tsv`, `ppt_verdicts.tsv`, `per_sample.tsv`, `summary.json`
- `figures/fig1_gsea_primary_nes.png` … `fig4_ppt_verdicts.png`
