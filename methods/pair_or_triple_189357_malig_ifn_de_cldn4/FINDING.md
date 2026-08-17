# Triple GSE131907+GSE205335+GSE189357: CLDN4-only malignant IFN/MHC/TJ DE

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. **Not infiltrate.**
Tumor-cell-intrinsic program only. No GSE148071.

The pair GSE123902+GSE189357 %pos n=22 ρ=−0.638 (PR #459) already has a
malignant IFN/MHC/TJ DE (`methods/pair_123902_189357_malig_ifn_de_cldn4`,
PR #469). This extra is the remaining 189357-axis triple that differs.

Given triple cut (PR #459, **not re-audited**): GSE131907 + GSE205335 +
GSE189357 malignant CLDN4 %pos vs same-unit T/NK, n=52, Spearman ρ=−0.497
(p=0.00035). This extra does **not** re-audit that T/NK ρ.

**Thesis (already correct):** CLDN4-high malignant cells should be IFN/MHC-I
down, TJ up on the malignant cell itself.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
GSE131907 is sample-level. GSE205335 / GSE189357 are patient-level.
p-values are descriptive.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE131907 + GSE205335 + GSE189357 | Pair 123902+189357 already filled; **no GSE148071** |
| Units | PR #459 %pos n=52 (21 samples + 22 patients + 9 patients) | Mixed units (sample vs patient); do not quote a cell-level n |
| Split | Within-cohort marker/author-malignant CLDN4 **%pos Q4 vs Q1** | Quartile tails; mid quartiles unused in binary DE |
| Continuous | CLDN4 %pos z, units in the malignant matrix, cohort covariates | Linear; not a causal model |
| Malignant | GSE131907/GSE205335 author-malig UMI-sum; GSE189357 marker-malig UMI-sum | Marker gate ≠ author annotation |
| Dropped from DE | P4001 (27 malignant cells) | In the locked n=52 T/NK vector; out of malignant matrix (n_mal≥30) |
| T/NK | **not a DE compartment** | T/NK ρ taken as given from PR #459; not re-audited |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (single) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I/APM (custom), TJ (KEGG/GO + junction focal; **CLDN4 held out**) | MHC-II and chemokine panels are not the claim |

## Honest n

PR #459 locked labels (malignant CLDN4 %pos; T/NK ρ not re-scored):

- GSE131907: n=21 tumor-origin samples (author malignant, n_mal≥20). Q1=6 Q4=5.
- GSE205335: n=22 patients (author malignant). Q1=6 Q4=6.
- GSE189357: n=9 patients (marker-malignant). Q1=3 Q4=2.

Locked Q1: EBUS_13, EBUS_15, EBUS_49, NS_02, NS_06, NS_16, P1015, P1062, P1063, P1090, P1119, P4001, TD2, TD4, TD7 (n=15).
Locked Q4: EBUS_19, EBUS_28, NS_03, NS_04, NS_07, P1016, P1025, P1037, P1084, P1089, P1115, TD6, TD9 (n=13).

Out of malignant DE: P4001 (GSE205335, 27 malignant cells, Q1).

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 combined | 14/13 | 16174 | cohort covariates; P4001 out |
| Q4 vs Q1 GSE131907 | 6/5 | 15738 | sample-level |
| Q4 vs Q1 GSE205335 | 5/6 | 18035 | P4001 out of matrix |
| Q4 vs Q1 GSE189357 | 3/2 | 0 | thin tail n_Q1=3 n_Q4=2 (need ≥3 each) |
| continuous combined | n=51 | 16980 | CLDN4 %pos z; units in the matrix (51 if P4001 out) |

Do not quote n=52 for the malignant DE. n=52 is the PR #459 T/NK cut.
Malignant DE n is the units actually in the UMI-sum matrices.

## Family DE (combined Q4 vs Q1) — n / logFC / p

Family score = mean log2(TMM-CPM+1) of family genes present. Positive logFC =
higher in CLDN4-high (Q4) than CLDN4-low (Q1). Machine table:
`tables/family_de.tsv`.

Expected under the thesis: IFN down, MHC-I/APM down, TJ up.

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 27 | 14 / 13 | 221 | -0.519 | 0.02262 | 0.06785 |
| MHC-I/APM | 27 | 14 / 13 | 21 | -0.740 | 0.04661 | 0.06991 |
| TJ | 27 | 14 / 13 | 195 | +0.061 | 0.3736 | 0.3736 |

Continuous family-score DE (units in the malignant matrix, CLDN4 %pos z):

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 51 | — | 221 | -0.117 | 0.2143 | 0.2143 |
| MHC-I/APM | 51 | — | 21 | -0.236 | 0.08356 | 0.1253 |
| TJ | 51 | — | 195 | +0.045 | 0.08091 | 0.1253 |

Per-cohort family scores (GSE189357 Q4 vs Q1 is thin and may be skipped):

| cohort | family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---|---:|---|---:|---:|---:|---:|
| GSE131907 | IFN | 11 | 6 / 5 | 221 | +0.079 | 0.7915 | 0.9122 |
| GSE131907 | MHC-I/APM | 11 | 6 / 5 | 21 | -0.063 | 0.9122 | 0.9122 |
| GSE131907 | TJ | 11 | 6 / 5 | 195 | +0.275 | 0.004363 | 0.01309 |
| GSE205335 | IFN | 11 | 5 / 6 | 221 | -1.243 | 0.002251 | 0.006754 |
| GSE205335 | MHC-I/APM | 11 | 5 / 6 | 21 | -1.605 | 0.01769 | 0.02653 |
| GSE205335 | TJ | 11 | 5 / 6 | 195 | -0.115 | 0.3446 | 0.3446 |

CLDN4 itself +1.77 (p=0.02965, FDR=0.2964) — direction check on the split gene (held out of the TJ family). n_Q1=14 n_Q4=13.

## Gene-level family members (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high.

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---|
| IFN | 218 | 54 (3/51) | 0 | -0.462 | CCL5 (-2.627, 3.72e-05, 0.05466) |
| MHC-I/APM | 21 | 5 (0/5) | 0 | -0.711 | TAP2 (-1.623, 0.0021, 0.1226) |
| TJ | 189 | 26 (20/6) | 0 | +0.063 | TBCD (+0.857, 0.0005054, 0.0772) |

Headline family genes (combined Q4 vs Q1, lowest p):

| family | gene | n_Q1 | n_Q4 | logFC | p | FDR |
|---|---|---:|---:|---:|---:|---:|
| IFN | CCL5 | 14 | 13 | -2.627 | 3.72e-05 | 0.05466 |
| IFN | GZMA | 14 | 13 | -2.348 | 5.43e-05 | 0.05607 |
| IFN | GBP4 | 14 | 13 | -2.518 | 9.85e-05 | 0.05607 |
| TJ | TBCD | 14 | 13 | +0.857 | 0.0005054 | 0.0772 |
| IFN | CSF2RB | 14 | 13 | -0.607 | 0.000577 | 0.0772 |
| IFN | GBP2 | 14 | 13 | -2.298 | 0.0009399 | 0.09336 |
| TJ | DLG3 | 14 | 13 | +1.050 | 0.001064 | 0.096 |
| IFN | CD86 | 14 | 13 | -1.113 | 0.001533 | 0.1112 |
| TJ | CLDN3 | 14 | 13 | +1.874 | 0.001559 | 0.1112 |
| IFN | LCP2 | 14 | 13 | -1.757 | 0.001574 | 0.1112 |
| IFN | CD69 | 14 | 13 | -2.144 | 0.001615 | 0.1125 |
| IFN | AUTS2 | 14 | 13 | +2.544 | 0.001761 | 0.1149 |

### Continuous CLDN4 %pos (cohort covariates)

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---|
| IFN | 221 | 40 (7/33) | 8 | -0.099 | VCAM1 (-0.785, 5.11e-05, 0.0101) |
| MHC-I/APM | 21 | 6 (0/6) | 2 | -0.195 | TAP1 (-0.667, 0.0002464, 0.02313) |
| TJ | 195 | 44 (33/11) | 11 | +0.018 | DLG3 (+0.492, 1.01e-05, 0.004616) |

## What this is not

- Not infiltrate / T/NK fraction DE and not a re-audit of PR #459 ρ.
- Not CellChat / LIANA / NicheNet.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not GSE148071 and not the already-filled GSE123902+GSE189357 pair DE.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ change.
- Genome-wide FDR on these n is expected to be thin; family scores and
  family median logFC / sign counts are the claim.

## Files

- `tables/family_de.tsv` — **headline family DE** (n / logFC / p)
- `tables/de_q4q1_combined_families.tsv` — gene-level family members
- `tables/family_summary.tsv` — IFN / MHC-I/APM / TJ counts
- `tables/de_families.tsv` — family rows, all contrasts
- `tables/de_all.tsv` — all genes, all contrasts
- `tables/n_honest.tsv`
- `tables/sample_inventory.tsv`
- `figures/` — volcano, family heatmap, key-gene boxes, family-score boxes, forests, family-median bar, CLDN4 strip, n bars

Reproduce:

```bash
python3 methods/pair_or_triple_189357_malig_ifn_de_cldn4/analyze.py
```
