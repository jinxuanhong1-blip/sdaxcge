# Concordant-4: TACSTD2 Q4 vs Q1 malignant DEG → TJ → T/NK

PAPER FUNNEL public evidence. **TACSTD2 (TROP2) only** as the split.
Locked cohorts only: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Not GSE148071 / GSE127465 / GSE154826 / GSE200563 / E-MTAB-13526.
Does not re-fit the locked CLDN4 %pos vs T/NK ρ = −0.531.

Honest unit = patient / donor / sample. Do not quote cell counts as n.
TACSTD2 split uses malignant **pseudobulk expression** (log2 TMM-CPM+1),
not cell-level %pos (no per-cell matrix in this page).

## Honest n

| cohort | unit | vector n | Q1 / Q4 |
|---|---|---:|---|
| GSE123902 | donor | 13 | 4 / 3 |
| GSE131907 | sample | 21 | 6 / 5 |
| GSE205335 | patient | 21 | 6 / 5 |
| GSE189357 | patient | 9 | 3 / 2 |

DE contrast = stacked within-cohort Q4 vs Q1: **19 vs 15** (do not quote n=65 as the DE n). Expression matrix has **64** units.

## 1. DEG (OLS on log2 TMM-CPM+1, cohort covariates)

Positive logFC = higher in TACSTD2 Q4. Genes tested: 24083.
TACSTD2 QC: logFC = 3.904, p = 7.71e-07.

### Family scores

| family | n_genes | logFC | p | FDR |
|---|---:|---:|---:|---:|
| TJ | 206 | +0.223 | 0.001163 | 0.004187 |
| CLAUDIN_PANEL | 19 | +0.426 | 0.004114 | 0.006171 |
| EPITHELIAL_ADHESION | 23 | +0.734 | 0.001396 | 0.004187 |
| KERATIN | 52 | +0.929 | 0.002725 | 0.00545 |
| IFN | 222 | +0.348 | 0.1123 | 0.1348 |
| MHC-I/APM | 21 | +0.189 | 0.567 | 0.567 |

### Focal genes

| gene | logFC | t | p | FDR |
|---|---:|---:|---:|---:|
| TACSTD2 | +3.904 | +6.27 | 7.709e-07 | 0.01857 |
| KRT19 | +2.230 | +2.98 | 0.005731 | 0.6137 |
| CLDN4 | +1.632 | +2.66 | 0.01256 | 0.6712 |
| OCLN | +0.989 | +2.66 | 0.0126 | 0.6712 |
| CDH1 | +1.032 | +2.36 | 0.02546 | 0.7284 |
| TJP1 | +0.593 | +2.18 | 0.03733 | 0.7499 |
| CLDN1 | +1.603 | +2.16 | 0.0388 | 0.7499 |
| KRT8 | +1.143 | +1.95 | 0.0615 | 0.7536 |
| ELF3 | +1.192 | +1.80 | 0.08232 | 0.7698 |
| F11R | +0.447 | +1.77 | 0.08706 | 0.7752 |
| CLDN7 | +0.713 | +1.76 | 0.08905 | 0.7752 |
| KRT18 | +0.865 | +1.25 | 0.2223 | 0.8045 |
| EPCAM | +0.552 | +1.05 | 0.3004 | 0.8241 |
| CLDN3 | +0.652 | +0.98 | 0.336 | 0.8351 |
| HLA-B | +0.561 | +0.90 | 0.3781 | 0.8504 |
| CXCL10 | +0.658 | +0.83 | 0.4122 | 0.8626 |
| B2M | +0.239 | +0.62 | 0.5432 | 0.9003 |
| STAT1 | +0.105 | +0.36 | 0.7213 | 0.9437 |
| IRF1 | -0.105 | -0.31 | 0.7561 | 0.9502 |
| HLA-A | +0.126 | +0.23 | 0.8221 | 0.9652 |

## 2. ORA / GSEA — is TJ / claudin / adhesion top or near-top?

ORA query rule: **p<0.01 & |logFC|>0.25 (FDR arm too thin)**. n_up = 274, n_down = 79.
Background = all genes in the DE table. TACSTD2 held out of the UP query.

### ORA top 5 (all tested A8 + custom sets)

| rank | term | overlap | enrichment | p | FDR |
|---:|---|---:|---:|---:|---:|
| 1 | HALLMARK_APICAL_JUNCTION | 13 | 5.89 | 3.79e-07 | 1.28e-05 |
| 2 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 13 | 5.74 | 5.07e-07 | 1.28e-05 |
| 3 | GOBP_KERATINIZATION | 7 | 13.98 | 6.11e-07 | 1.28e-05 |
| 4 | HALLMARK_INFLAMMATORY_RESPONSE | 12 | 5.44 | 2.48e-06 | 3.66e-05 |
| 5 | HALLMARK_KRAS_SIGNALING_UP | 12 | 5.35 | 2.91e-06 | 3.66e-05 |

Highest-enrichment barrier ORA term: **GOBP_KERATINIZATION** (rank 3 of all ORA terms, enrichment 13.98, FDR 1.28e-05).
Best-ranked barrier ORA term: **HALLMARK_APICAL_JUNCTION** (rank 1, enrichment 5.89, FDR 1.28e-05).
**ORA verdict: tight-junction / apical-junction / keratin-barrier is top or near-top.**

### Prerank GSEA (OLS *t*, TACSTD2 dropped from rank; 1000 gene-set perms, seed=42)

Positive NES = enriched at the TACSTD2-Q4 end.

| rank | term | NES | nom p | FDR |
|---:|---|---:|---:|---:|
| 1 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | +3.200 | 0.000999 | 0.00146 |
| 2 | GOBP_KERATINIZATION | +2.825 | 0.000999 | 0.00146 |
| 3 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | +2.818 | 0.000999 | 0.00146 |
| 4 | HALLMARK_INFLAMMATORY_RESPONSE | +2.803 | 0.000999 | 0.00146 |
| 5 | GOBP_KERATINOCYTE_DIFFERENTIATION | +2.717 | 0.000999 | 0.00146 |

Best barrier GSEA set: **GOBP_KERATINIZATION** (NES +2.825, FDR 0.00146, rank 2 among positive-NES sets).
**Verdict: TJ/barrier is near-top (rank ≤ 5 among positive NES).**

### Enrichr (GO BP 2023 / KEGG 2021 / Hallmark 2020)

Ran on the UP list. 1847 terms returned; 29 matched a TJ/claudin/adhesion/keratin regex. See `tables/enrichr_up.tsv` and `tables/enrichr_barrier_hits.tsv`.

## 3. TJ score vs T/NK

Broad TJ score = mean log2(TMM-CPM+1) of 206 TJ/adhesion genes (TACSTD2 held out). Core panel = CLDN1, CLDN3, CLDN4, CLDN7, OCLN, F11R, TJP1. Unit = patient/donor/sample.

| score | method | N | effect | p | I² |
|---|---|---:|---|---:|---:|
| broad TJ | DL meta Spearman | 64 | ρ = +0.137 (95% CI -0.365 to +0.577) | 0.605 | 70.5% |
| broad TJ | stacked Q4 vs Q1 MWU | 19/15 | r_rb = +0.116, Δmedian = +0.0310 | 0.579 | — |
| TJ core | DL meta Spearman | 64 | ρ = -0.191 (95% CI -0.571 to +0.257) | 0.406 | 61.5% |
| TJ core | stacked Q4 vs Q1 MWU | 19/15 | r_rb = -0.291, Δmedian = -0.1417 | 0.155 | — |
| TACSTD2 | DL meta Spearman | 64 | ρ = +0.155 | 0.261 | 0.0% |
| CLDN4 %pos (recomputed) | DL meta Spearman | 64 | ρ = -0.523 | 2.86e-05 | 0.0% |

Broad TJ member rhos: `GSE123902:-0.330,GSE131907:-0.208,GSE205335:0.583,GSE189357:0.433`. Core member rhos: `GSE123902:-0.275,GSE131907:-0.587,GSE205335:0.234,GSE189357:0.000`.

**Honest TJ→T/NK call: null on this page.** Broad TJ and the 7-gene core do not reproduce the locked CLDN4 %pos vs T/NK anti-correlation. GSE205335 alone is positive for broad TJ (see member rhos); the other three cohorts are negative or flat. Do not write “TJ-high tumors exclude T/NK” from this TACSTD2 funnel alone.

CLDN4 %pos vs T/NK remains the locked primary (ρ = −0.531, n = 65) from PR #503/#539. This page’s recomputed CLDN4 %pos meta on the expression-overlapping units is ρ = -0.523 (n = 64).
This page adds the TJ-score arm after the TACSTD2-high DEG → barrier enrichment step.

## Funnel verdict (computed, not wished)

1. **TACSTD2 Q4 vs Q1 DEG:** barrier / keratin / epithelial-adhesion family scores are up.
2. **ORA:** barrier term near-top (GOBP_KERATINIZATION, rank 3). Apical-junction / keratinization lead the ORA table.
3. **GSEA:** keratin/barrier near-top (GOBP_KERATINIZATION rank 2 among positive NES). Classical KEGG TJ is FDR-significant but not NES-top.
4. **TJ score vs T/NK:** null. The immune-exclusion arm stays on **CLDN4 %pos**, not on TACSTD2 expression or the broad TJ score.

## Not claimed

- Cell-level TACSTD2 %pos (not in this matrix).
- Private 8KL / KD coculture.
- Visium spatial exclusion.
- Merging non-concordant accessions.
- “TACSTD2-high = IFN/MHC-low” (family scores here are NS / slightly up).
- “Broad TJ-high = T/NK-low” on concordant-4 (null on this page).

## Reproduce

```bash
python3 methods/concordant4_tacstd2_malignant_deg/analyze.py
```

