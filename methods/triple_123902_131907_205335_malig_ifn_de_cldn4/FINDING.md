# Triple GSE123902+GSE131907+GSE205335 — malignant-cell-intrinsic CLDN4 DE

ADDITIVE. **CLDN4-only.** Tumor-cell-intrinsic program DE.
Not CellChat. Not T/NK infiltrate. No dual-high TACSTD2×CLDN4.
**Not the 7-pool. Not +GSE148071.**

The triple that already differs is taken from PR #459 (malignant CLDN4 **%pos**,
n=56, ρ=−0.522 vs T/NK, Q4 vs Q1 r=−0.735). That T/NK rho is **not re-audited**.
This PR asks a different question: in the same units, do CLDN4-high malignant
cells themselves have IFN / MHC-I/APM down and TJ up.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
Within-cohort CLDN4 %pos Q4 vs Q1, then stacked with cohort covariates.
GSE123902 is donor-level. GSE131907 is sample-level. GSE205335 is patient-level.
p-values are descriptive.

Thesis (not re-derived): CLDN4 KD / low in the malignant cell raises that
cell's own IFN / MHC-I / APM program and lowers TJ. Observationally,
CLDN4-high malignant cells should have IFN/MHC down, TJ up.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE123902 + GSE131907 + GSE205335 only | Not the 7-pool. Not +GSE148071 / 189357 / 127465 / 207422 |
| Units | PR #459 %pos n=56 (13 donors + 21 samples + 22 patients) | Mixed unit (donor / sample / patient) |
| Split | **Within-cohort** malignant CLDN4 **%pos Q4 vs Q1**, then stacked | Quartile tails; mid quartiles unused in binary DE |
| Continuous | CLDN4 %pos z-scored within the stacked matrix, cohort covariates | Linear; not a causal model |
| Malignant | GSE123902 marker-malig UMI-sum; GSE131907 / GSE205335 author-malig UMI-sum | Marker gate ≠ author annotation |
| T/NK | **not run** | Do not re-audit PR #459 ρ=−0.522 |
| Model | ~ cohort + CLDN4_Q4 (stacked); ~ CLDN4_Q4 (single) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I/APM (custom), TJ (KEGG/GO + focal, **CLDN4 held out**), keratin (KRT_EPITHELIAL), chemokine panel | Chemokine panel is compact, not MSigDB C2 |

## Honest n

PR #459 locked labels (malignant CLDN4 %pos):

- GSE123902: n=13 donors (tumor/met; marker-malignant; normals dropped). Q1=4 Q4=3.
- GSE131907: n=21 samples (tumor origins, n_mal≥20, author malignant). Q1=6 Q4=5.
- GSE205335: n=22 patients (author malignant). Q1=6 Q4=6. Dropped from malignant DE (not in UMI-sum): P4001.

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 stacked | 15/14 | 15615 | cohort covariates; malignant only |
| continuous stacked | n=55 | 16327 | CLDN4 %pos z |

Per-cohort Q4 vs Q1 n is in `tables/n_honest.tsv`. Do not quote n=56 as the DE n:
the stacked Q4 vs Q1 uses only the tails that are in the count matrices, and
P4001 (27 malignant cells) is a PR #459 label but not a malignant-DE unit.

## Family DE (stacked Q4 vs Q1) — n / logFC / p

Family score = mean log2(TMM-CPM+1) of family genes present. Positive logFC =
higher in CLDN4-high (Q4) than CLDN4-low (Q1). Machine table:
`tables/family_de.tsv`.

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 29 | 15 / 14 | 219 | -0.640 | 0.004109 | 0.01027 |
| MHC-I/APM | 29 | 15 / 14 | 21 | -0.851 | 0.015 | 0.02499 |
| TJ | 29 | 15 / 14 | 194 | +0.059 | 0.4287 | 0.4287 |
| keratin | 29 | 15 / 14 | 17 | -0.485 | 0.3353 | 0.4192 |
| chemokine | 29 | 15 / 14 | 24 | -1.038 | 0.0001427 | 0.0007137 |

Continuous family-score DE (units in the count matrices, CLDN4 %pos z):

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 55 | — | 219 | -0.166 | 0.06527 | 0.1026 |
| MHC-I/APM | 55 | — | 21 | -0.290 | 0.02772 | 0.0693 |
| TJ | 55 | — | 194 | +0.047 | 0.08209 | 0.1026 |
| keratin | 55 | — | 17 | -0.067 | 0.7191 | 0.7191 |
| chemokine | 55 | — | 24 | -0.327 | 0.004347 | 0.02173 |

Per-cohort family scores (within-cohort Q4 vs Q1; GSE123902 tails are thin):

| cohort | family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---|---:|---|---:|---:|---:|---:|
| GSE123902 | IFN | 7 | 4 / 3 | 219 | -0.759 | 0.07079 | 0.1735 |
| GSE123902 | MHC-I/APM | 7 | 4 / 3 | 21 | -0.844 | 0.08333 | 0.1735 |
| GSE123902 | TJ | 7 | 4 / 3 | 194 | +0.048 | 0.837 | 0.837 |
| GSE123902 | keratin | 7 | 4 / 3 | 17 | -0.788 | 0.3884 | 0.4854 |
| GSE123902 | chemokine | 7 | 4 / 3 | 24 | -1.419 | 0.1041 | 0.1735 |
| GSE131907 | IFN | 11 | 6 / 5 | 219 | +0.033 | 0.9149 | 0.9149 |
| GSE131907 | MHC-I/APM | 11 | 6 / 5 | 21 | -0.117 | 0.8364 | 0.9149 |
| GSE131907 | TJ | 11 | 6 / 5 | 194 | +0.228 | 0.006941 | 0.03471 |
| GSE131907 | keratin | 11 | 6 / 5 | 17 | +0.568 | 0.5512 | 0.9149 |
| GSE131907 | chemokine | 11 | 6 / 5 | 24 | -0.424 | 0.08339 | 0.2085 |
| GSE205335 | IFN | 11 | 5 / 6 | 219 | -1.239 | 0.00182 | 0.004549 |
| GSE205335 | MHC-I/APM | 11 | 5 / 6 | 21 | -1.589 | 0.01703 | 0.02838 |
| GSE205335 | TJ | 11 | 5 / 6 | 194 | -0.104 | 0.3543 | 0.3543 |
| GSE205335 | keratin | 11 | 5 / 6 | 17 | -1.347 | 0.08467 | 0.1058 |
| GSE205335 | chemokine | 11 | 5 / 6 | 24 | -1.413 | 0.001466 | 0.004549 |

CLDN4 itself (held out of TJ) combined Q4 vs Q1: logFC=+1.724, p=0.02474, n_Q1=15, n_Q4=14 — direction check on the split gene.

Stacked Q4 vs Q1 **directions match the thesis on IFN / MHC-I/APM / chemokine
(down) and TJ (up)**. Family-score p is not thin for IFN, MHC-I/APM, and
chemokine. TJ is the right sign but p=0.43 on the stacked split; continuous
TJ is +0.047 (p=0.082). GSE205335 carries the IFN/MHC down; GSE131907 carries
the TJ up; GSE123902 is the same signs on thin tails (4 vs 3). Keratin is
null / mixed. Do not quote the stacked TJ p as a hit.

## Gene-level family members (stacked Q4 vs Q1)

Positive logFC = higher in CLDN4-high.

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---|
| IFN | 219 | 74 (2/72) | 5 | -0.602 | CCL5 (-2.670, 4.95e-06, 0.008566) |
| MHC-I/APM | 21 | 8 (0/8) | 0 | -0.744 | TAP2 (-1.743, 0.0005317, 0.05125) |
| TJ | 187 | 41 (28/13) | 2 | +0.069 | TBCD (+0.956, 9.06e-05, 0.02283) |
| keratin | 16 | 2 (0/2) | 0 | -0.260 | KRT6A (-3.346, 0.006377, 0.142) |
| chemokine | 23 | 13 (0/13) | 4 | -1.017 | CCL5 (-2.670, 4.95e-06, 0.008566) |

Headline genes (stacked Q4 vs Q1, family members only, lowest p):

| family | gene | n_Q1 | n_Q4 | logFC | p | FDR |
|---|---|---:|---:|---:|---:|---:|
| IFN|chemokine | CCL5 | 15 | 14 | -2.670 | 4.95e-06 | 0.008566 |
| IFN | GZMA | 15 | 14 | -2.393 | 1.12e-05 | 0.009731 |
| IFN | CD69 | 15 | 14 | -2.795 | 4.00e-05 | 0.0166 |
| IFN | GBP2 | 15 | 14 | -2.687 | 7.48e-05 | 0.02058 |
| chemokine | CCL3 | 15 | 14 | -2.346 | 7.51e-05 | 0.02058 |
| TJ | TBCD | 15 | 14 | +0.956 | 9.06e-05 | 0.02283 |
| chemokine | XCL2 | 15 | 14 | -1.672 | 0.0001552 | 0.03046 |
| TJ | CLDN3 | 15 | 14 | +2.003 | 0.0003004 | 0.04114 |
| chemokine | CCL4 | 15 | 14 | -2.245 | 0.000305 | 0.04114 |
| IFN | GBP4 | 15 | 14 | -2.313 | 0.0004116 | 0.04613 |
| MHC-I/APM | TAP2 | 15 | 14 | -1.743 | 0.0005317 | 0.05125 |
| IFN | CSF2RB | 15 | 14 | -1.147 | 0.0006036 | 0.05343 |

## What this is not

- Not CellChat / LIANA / NicheNet and not a T/NK infiltrate test.
- Not a re-audit of PR #459 ρ=−0.522.
- Not the 7-cohort pool and not +GSE148071.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ/keratin/chemokine change.
- Genome-wide FDR on these n is expected to be thin; family scores are the claim.

## Files

- `tables/family_de.tsv` — **headline family DE** (n / logFC / p)
- `tables/de_all.tsv` — all genes, all contrasts
- `tables/de_families.tsv` — IFN / MHC-I/APM / TJ / keratin / chemokine rows
- `tables/de_q4q1_combined_families.tsv` — gene-level headline
- `tables/family_summary.tsv` — gene-level family counts
- `tables/n_honest.tsv`
- `tables/sample_inventory.tsv`
- `figures/` — volcano, family heatmap, key-gene boxes, family-score boxes, forests (combined + stacked per-cohort), CLDN4 strip, n bars

Reproduce:

```bash
python3 methods/triple_123902_131907_205335_malig_ifn_de_cldn4/analyze.py
```
