# Pair GSE189357+GSE205335 — malignant-cell-intrinsic CLDN4 DE

ADDITIVE. **CLDN4-only.** Tumor-cell-intrinsic program DE.
Not CellChat. Not T/NK infiltrate. No dual-high TACSTD2×CLDN4. No GSE148071.

The pair that already differs is taken from PR #459 (malignant CLDN4 **%pos**,
n=31, ρ=−0.478 vs T/NK; Q4 vs Q1 r=−0.750). That T/NK rho is **not
re-audited**. This PR asks a different question: in the same 31 patients,
do CLDN4-high malignant cells themselves have IFN / MHC-I/APM down and TJ up.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
muscat `pbDS` is the same collapse (patient × cell-type UMI-sum → bulk DE).
Both cohorts are **patient-level**. p-values are descriptive.

**Thesis (already correct):** CLDN4-high malignant = IFN/MHC down, TJ up
(KD-like low = IFN up).

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE189357 + GSE205335 only | Not GSE148071 / 123902 / 131907 / a bigger merge |
| Units | PR #459 %pos n=31 (9 + 22 patients) | Malignant definitions differ (marker vs author) |
| Split | Within-cohort malignant CLDN4 **%pos Q4 vs Q1** | Quartile tails; mid quartiles unused in binary DE |
| Continuous | CLDN4 %pos z-scored, all 31 patients, cohort covariate | Linear; not a causal model |
| Malignant | GSE189357 marker-malig UMI-sum; GSE205335 author-malig UMI-sum | Marker gate ≠ author annotation |
| T/NK | **not run** | Do not re-audit PR #459 ρ=−0.478 / Q4 r=−0.750 |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (single) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I/APM (custom), TJ (KEGG/GO + focal, **CLDN4 held out**); keratin / chemokine are extra panels | Chemokine panel is compact, not MSigDB C2 |

## Honest n

PR #459 locked labels (malignant CLDN4 %pos):

- GSE189357: n=9 patients (marker-malignant; one 10x sample / patient). Q1=3 Q4=2.
- GSE205335: n=22 patients (author `Malignant cells`; 4 extra GEO patients with 0 malignant already out). Q1=6 Q4=6. GSE205335 Q4 histology: ADC=2, SCLC=3, SQ=1.

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 combined | 9/8 | 17885 | cohort covariate; malignant only |
| continuous combined | n=31 | 19544 | CLDN4 %pos z |

Per-cohort Q4 vs Q1 n is in `tables/n_honest.tsv`. Cells are not n.
GSE189357 Q4 is a thin tail (need ≥3 each for single-cohort binary DE).

## Family DE (combined Q4 vs Q1) — n / logFC / p

Family score = mean log2(TMM-CPM+1) of family genes present. Positive logFC =
higher in CLDN4-high (Q4) than CLDN4-low (Q1). Machine table:
`tables/family_de.tsv`.

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 17 | 9 / 8 | 222 | -0.887 | 0.001772 | 0.00443 |
| MHC-I/APM | 17 | 9 / 8 | 21 | -1.134 | 0.01389 | 0.02315 |
| TJ | 17 | 9 / 8 | 194 | -0.012 | 0.8947 | 0.8947 |
| keratin | 17 | 9 / 8 | 16 | -0.779 | 0.183 | 0.2288 |
| chemokine | 17 | 9 / 8 | 25 | -1.204 | 0.000551 | 0.002755 |

On this pair the claim families **split**: IFN and MHC-I/APM are down in
CLDN4-high malignant cells (combined Q4 vs Q1 and continuous). **TJ is
flat** (logFC ≈ 0). GSE205335 Q4 is 3/6 SCLC, which is a plausible reason
the TJ arm does not rise. Chemokine (extra) is also down.

Continuous family-score DE (same 31 patients, CLDN4 %pos z):

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 31 | — | 222 | -0.341 | 0.00408 | 0.009957 |
| MHC-I/APM | 31 | — | 21 | -0.480 | 0.005974 | 0.009957 |
| TJ | 31 | — | 194 | -0.006 | 0.8609 | 0.8609 |
| keratin | 31 | — | 16 | -0.468 | 0.0668 | 0.0835 |
| chemokine | 31 | — | 25 | -0.424 | 0.00302 | 0.009957 |

Per-cohort family scores (GSE189357 Q4 vs Q1 is thin, n_Q4=2, skipped):

| cohort | family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---|---:|---|---:|---:|---:|---:|
| GSE205335 | IFN | 12 | 6 / 6 | 222 | -1.147 | 0.001723 | 0.004307 |
| GSE205335 | MHC-I/APM | 12 | 6 / 6 | 21 | -1.461 | 0.01801 | 0.03002 |
| GSE205335 | TJ | 12 | 6 / 6 | 194 | -0.004 | 0.973 | 0.973 |
| GSE205335 | keratin | 12 | 6 / 6 | 16 | -1.258 | 0.112 | 0.14 |
| GSE205335 | chemokine | 12 | 6 / 6 | 25 | -1.411 | 0.0005357 | 0.002679 |

GSE189357 continuous (n=9; the only single-cohort test that meets n):

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 9 | — | 222 | -0.059 | 0.5713 | 0.7455 |
| MHC-I/APM | 9 | — | 21 | -0.038 | 0.8156 | 0.8156 |
| TJ | 9 | — | 194 | -0.025 | 0.2908 | 0.727 |
| keratin | 9 | — | 16 | +0.174 | 0.09279 | 0.464 |
| chemokine | 9 | — | 25 | -0.137 | 0.5964 | 0.7455 |

CLDN4 itself (held out of TJ) combined Q4 vs Q1: logFC=+0.885, p=0.02644, n_Q1=9, n_Q4=8 — direction check on the split gene.

## Gene-level family members (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high.

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---|
| IFN | 219 | 76 (0/76) | 0 | -0.912 | CCL5 (-3.633, 0.0003119, 0.2584) |
| MHC-I/APM | 21 | 13 (0/13) | 0 | -1.288 | PSMB9 (-2.219, 0.003639, 0.2709) |
| TJ | 189 | 26 (11/15) | 0 | -0.017 | ACTR3 (-0.942, 0.00268, 0.2663) |
| keratin | 14 | 1 (0/1) | 0 | -0.746 | KRT6A (-3.654, 0.005154, 0.2959) |
| chemokine | 23 | 9 (0/9) | 0 | -0.974 | CCL5 (-3.633, 0.0003119, 0.2584) |

Headline genes (combined Q4 vs Q1, family members only, lowest p):

| family | gene | n_Q1 | n_Q4 | logFC | p | FDR |
|---|---|---:|---:|---:|---:|---:|
| IFN|chemokine | CCL5 | 9 | 8 | -3.633 | 0.0003119 | 0.2584 |
| IFN | OAS2 | 9 | 8 | -3.072 | 0.0003341 | 0.2584 |
| IFN | GBP4 | 9 | 8 | -3.156 | 0.0003584 | 0.2584 |
| IFN | GZMA | 9 | 8 | -3.322 | 0.0003919 | 0.2584 |
| IFN | SLAMF7 | 9 | 8 | -1.696 | 0.001046 | 0.2584 |
| IFN | SAMD9L | 9 | 8 | -2.403 | 0.00135 | 0.2584 |
| IFN | TMEM140 | 9 | 8 | -1.872 | 0.001545 | 0.2584 |
| IFN | MYD88 | 9 | 8 | -1.164 | 0.001581 | 0.2584 |
| IFN | GBP2 | 9 | 8 | -3.006 | 0.001621 | 0.2584 |
| IFN | EPSTI1 | 9 | 8 | -2.632 | 0.001864 | 0.2646 |
| chemokine | XCL2 | 9 | 8 | -2.353 | 0.002004 | 0.2648 |
| IFN | NFKB1 | 9 | 8 | -1.561 | 0.002211 | 0.2648 |

## What this is not

- Not CellChat / LIANA / NicheNet and not a T/NK infiltrate test.
- Not a re-audit of PR #459 ρ=−0.478 or Q4 r=−0.750.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not GSE148071 and not a merge beyond GSE189357+GSE205335.
- Not evidence that CLDN4 *causes* IFN/MHC/TJ/keratin/chemokine change.
- GSE205335 Q4 can be SCLC-heavy; histology is in `tables/sample_inventory.tsv`.
- Genome-wide FDR on these n is expected to be thin; family scores are the claim.

## Files

- `tables/family_de.tsv` — **headline family DE** (n / logFC / p)
- `tables/de_all.tsv` — all genes, all contrasts
- `tables/de_families.tsv` — IFN / MHC-I/APM / TJ / keratin / chemokine rows
- `tables/de_q4q1_combined_families.tsv` — gene-level headline
- `tables/family_summary.tsv` — gene-level family counts
- `tables/n_honest.tsv`
- `tables/sample_inventory.tsv`
- `figures/` — volcano, family heatmap, key-gene boxes, family-score boxes, forests, CLDN4 strip, n bars, family-score vs %pos

Reproduce:

```bash
python3 methods/pair_189357_205335_malig_ifn_de_cldn4/download.py
python3 methods/pair_189357_205335_malig_ifn_de_cldn4/build_malignant_pseudobulk.py
python3 methods/pair_189357_205335_malig_ifn_de_cldn4/analyze.py
```
