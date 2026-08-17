# Pair GSE123902+GSE131907 — malignant-cell-intrinsic CLDN4 DE

ADDITIVE. **CLDN4-only.** Tumor-cell-intrinsic program DE.
Not CellChat. Not T/NK infiltrate. No dual-high TACSTD2×CLDN4. No GSE148071.

The pair that already differs is taken from PR #459 (malignant CLDN4 **%pos**,
n=34, ρ=−0.575 vs T/NK). That T/NK rho is **not re-audited**. This PR asks a
different question: in the same 34 units, do CLDN4-high malignant cells
themselves have IFN / MHC-I/APM down and TJ up.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
muscat `pbDS` is the same collapse (patient × cell-type UMI-sum → bulk DE).
GSE123902 is donor-level. GSE131907 is sample-level. p-values are descriptive.

Thesis (not re-derived): CLDN4 KD / low in the malignant cell raises that
cell's own IFN / MHC-I / APM program and lowers TJ. Observationally,
CLDN4-high malignant cells should have IFN/MHC down, TJ up.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE123902 + GSE131907 only | Not GSE148071 / 205335 / 189357 / a bigger merge |
| Units | PR #459 %pos n=34 (13 donors + 21 samples) | GSE131907 is **sample-level**; GSE123902 is **donor-level** |
| Split | Within-cohort malignant CLDN4 **%pos Q4 vs Q1** | Quartile tails; mid quartiles unused in binary DE |
| Continuous | CLDN4 %pos z-scored, all 34 units, cohort covariate | Linear; not a causal model |
| Malignant | GSE123902 marker-malig UMI-sum; GSE131907 author-malig UMI-sum | Marker gate ≠ author annotation |
| T/NK | **not run** | Do not re-audit PR #459 ρ=−0.575 |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (single) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I/APM (custom), TJ (KEGG/GO + focal, **CLDN4 held out**), keratin (KRT_EPITHELIAL), chemokine panel | Chemokine panel is compact, not MSigDB C2 |

## Honest n

PR #459 locked labels (malignant CLDN4 %pos):

- GSE123902: n=13 donors (tumor/met; marker-malignant; normals dropped). Q1=4 Q4=3.
- GSE131907: n=21 samples (tumor origins, n_mal≥20, author malignant). Q1=6 Q4=5.

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 combined | 10/8 | 15311 | cohort covariate; malignant only |
| continuous combined | n=34 | 16490 | CLDN4 %pos z |

Per-cohort Q4 vs Q1 n is in `tables/n_honest.tsv`. Do not quote a pooled n that
ignores the donor-vs-sample unit difference.

## Family DE (combined Q4 vs Q1) — n / logFC / p

Family score = mean log2(TMM-CPM+1) of family genes present. Positive logFC =
higher in CLDN4-high (Q4) than CLDN4-low (Q1). Machine table:
`tables/family_de.tsv`.

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 18 | 10 / 8 | 216 | -0.255 | 0.3048 | 0.4082 |
| MHC-I/APM | 18 | 10 / 8 | 21 | -0.381 | 0.3265 | 0.4082 |
| TJ | 18 | 10 / 8 | 186 | +0.176 | 0.09342 | 0.2336 |
| keratin | 18 | 10 / 8 | 16 | +0.050 | 0.9414 | 0.9414 |
| chemokine | 18 | 10 / 8 | 24 | -0.783 | 0.02416 | 0.1208 |

Continuous family-score DE (same 34 units, CLDN4 %pos z):

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 34 | — | 216 | -0.028 | 0.7861 | 0.7861 |
| MHC-I/APM | 34 | — | 21 | -0.123 | 0.4207 | 0.5258 |
| TJ | 34 | — | 186 | +0.092 | 0.01721 | 0.08607 |
| keratin | 34 | — | 16 | +0.254 | 0.2445 | 0.4075 |
| chemokine | 34 | — | 24 | -0.240 | 0.1079 | 0.2696 |

Per-cohort family scores (thin tails; GSE123902 Q4 vs Q1 is 3 vs 4):

| cohort | family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---|---:|---|---:|---:|---:|---:|
| GSE123902 | IFN | 7 | 4 / 3 | 216 | -0.737 | 0.08065 | 0.1858 |
| GSE123902 | MHC-I/APM | 7 | 4 / 3 | 21 | -0.822 | 0.0947 | 0.1858 |
| GSE123902 | TJ | 7 | 4 / 3 | 186 | +0.065 | 0.7963 | 0.7963 |
| GSE123902 | keratin | 7 | 4 / 3 | 16 | -0.789 | 0.3997 | 0.4997 |
| GSE123902 | chemokine | 7 | 4 / 3 | 24 | -1.382 | 0.1115 | 0.1858 |
| GSE131907 | IFN | 11 | 6 / 5 | 216 | +0.047 | 0.8789 | 0.8789 |
| GSE131907 | MHC-I/APM | 11 | 6 / 5 | 21 | -0.103 | 0.8555 | 0.8789 |
| GSE131907 | TJ | 11 | 6 / 5 | 186 | +0.246 | 0.006958 | 0.03479 |
| GSE131907 | keratin | 11 | 6 / 5 | 16 | +0.578 | 0.5565 | 0.8789 |
| GSE131907 | chemokine | 11 | 6 / 5 | 24 | -0.406 | 0.08687 | 0.2172 |

CLDN4 itself (held out of TJ) combined Q4 vs Q1: logFC=+2.252, p=0.06049, n_Q1=10, n_Q4=8 — direction check on the split gene.

## Gene-level family members (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high.

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---|
| IFN | 213 | 27 (7/20) | 0 | -0.290 | PSMA2 (-0.932, 0.0002609, 0.08904) |
| MHC-I/APM | 21 | 2 (0/2) | 0 | -0.346 | TAP2 (-1.209, 0.01641, 0.2881) |
| TJ | 181 | 34 (23/11) | 0 | +0.122 | TBCD (+1.164, 0.0001795, 0.08201) |
| keratin | 14 | 0 (0/0) | 0 | +0.282 | KRT10 (+0.780, 0.05808, 0.4257) |
| chemokine | 20 | 6 (0/6) | 0 | -1.105 | CXCR4 (-2.110, 0.00136, 0.1504) |

Headline genes (combined Q4 vs Q1, family members only, lowest p):

| family | gene | n_Q1 | n_Q4 | logFC | p | FDR |
|---|---|---:|---:|---:|---:|---:|
| TJ | TBCD | 10 | 8 | +1.164 | 0.0001795 | 0.08201 |
| IFN | PSMA2 | 10 | 8 | -0.932 | 0.0002609 | 0.08904 |
| TJ | ZEB1 | 10 | 8 | -2.069 | 0.0003666 | 0.09942 |
| IFN | CD69 | 10 | 8 | -2.557 | 0.0008701 | 0.1359 |
| TJ | PRKAA2 | 10 | 8 | +2.225 | 0.001016 | 0.1465 |
| chemokine | CXCR4 | 10 | 8 | -2.110 | 0.00136 | 0.1504 |
| TJ | CLDN3 | 10 | 8 | +1.989 | 0.001388 | 0.1518 |
| IFN|chemokine | CCL5 | 10 | 8 | -1.837 | 0.001488 | 0.1557 |
| TJ | F11R | 10 | 8 | +0.960 | 0.001832 | 0.1564 |
| IFN | GZMA | 10 | 8 | -1.696 | 0.002497 | 0.1655 |
| TJ | MPP7 | 10 | 8 | +2.013 | 0.003021 | 0.1771 |
| chemokine | CCL3 | 10 | 8 | -2.233 | 0.003039 | 0.1771 |

## What this is not

- Not CellChat / LIANA / NicheNet and not a T/NK infiltrate test.
- Not a re-audit of PR #459 ρ=−0.575.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not GSE148071 and not a merge beyond GSE123902+GSE131907.
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
- `figures/` — volcano, family heatmap, key-gene boxes, family-score boxes, forests, CLDN4 strip, n bars

Reproduce:

```bash
python3 methods/pair_123902_131907_malig_ifn_de_cldn4/analyze.py
```
