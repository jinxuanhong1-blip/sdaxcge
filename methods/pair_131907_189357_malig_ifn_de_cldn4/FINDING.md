# Pair GSE131907+GSE189357 — malignant-cell-intrinsic CLDN4 DE

ADDITIVE. **CLDN4-only.** Tumor-cell-intrinsic program DE.
Not CellChat. Not T/NK infiltrate. No dual-high TACSTD2×CLDN4. No GSE148071.

The pair that already differs is taken from PR #459 (malignant CLDN4 **%pos**,
n=30, ρ=−0.542 vs T/NK; Q4 vs Q1 r=−0.619). That T/NK rho is **not
re-audited**. This PR asks a different question: in the same 30 units, do
CLDN4-high malignant cells themselves have IFN / MHC-I/APM down and TJ up.

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
muscat `pbDS` is the same collapse (patient × cell-type UMI-sum → bulk DE).
GSE131907 is sample-level. GSE189357 is patient-level. p-values are descriptive.

Thesis (not re-derived): CLDN4 KD / low in the malignant cell raises that
cell's own IFN / MHC-I / APM program and lowers TJ. Observationally,
CLDN4-high malignant cells should have IFN/MHC down, TJ up.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE131907 + GSE189357 only | Not GSE148071 / 123902 / 205335 / a bigger merge |
| Units | PR #459 %pos n=30 (21 samples + 9 patients) | GSE131907 is **sample-level**; GSE189357 is **patient-level** |
| Split | Within-cohort malignant CLDN4 **%pos Q4 vs Q1** | Quartile tails; mid quartiles unused in binary DE |
| Continuous | CLDN4 %pos z-scored, all 30 units, cohort covariate | Linear; not a causal model |
| Malignant | GSE131907 author-malig UMI-sum; GSE189357 marker-malig UMI-sum | Marker gate ≠ author annotation |
| T/NK | **not run** | Do not re-audit PR #459 ρ=−0.542 / Q4 r=−0.619 |
| Model | ~ cohort + CLDN4_Q4 (combined); ~ CLDN4_Q4 (single) | Small n; no muscat mixed model; no voom weights |
| Families | IFN (Hallmark IFNα/γ), MHC-I/APM (custom), TJ (KEGG/GO + focal, **CLDN4 held out**), keratin (KRT_EPITHELIAL), chemokine panel | Chemokine panel is compact, not MSigDB C2 |

## Honest n

PR #459 locked labels (malignant CLDN4 %pos):

- GSE131907: n=21 samples (tumor origins, n_mal≥20, author malignant). Q1=6 Q4=5.
- GSE189357: n=9 patients (marker-malignant; all 9 eligible). Q1=3 Q4=2.

| contrast | n_low / n_high | n_genes | note |
|---|---:|---:|---|
| Q4 vs Q1 combined | 9/7 | 15004 | cohort covariate; malignant only |
| continuous combined | n=30 | 15846 | CLDN4 %pos z |

GSE189357 Q4 vs Q1 is thin (Q4 n=2) and is **skipped** as a single-cohort
binary DE. Per-cohort n is in `tables/n_honest.tsv`. Do not quote a pooled n
that ignores the sample-vs-patient unit difference.

## Family DE (combined Q4 vs Q1) — n / logFC / p

Family score = mean log2(TMM-CPM+1) of family genes present. Positive logFC =
higher in CLDN4-high (Q4) than CLDN4-low (Q1). Machine table:
`tables/family_de.tsv`.

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 16 | 9 / 7 | 220 | -0.086 | 0.703 | 0.703 |
| MHC-I/APM | 16 | 9 / 7 | 21 | -0.217 | 0.5934 | 0.703 |
| TJ | 16 | 9 / 7 | 188 | +0.117 | 0.0803 | 0.2008 |
| keratin | 16 | 9 / 7 | 16 | +0.502 | 0.4684 | 0.703 |
| chemokine | 16 | 9 / 7 | 25 | -0.508 | 0.04311 | 0.2008 |

Continuous family-score DE (same 30 units, CLDN4 %pos z):

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 30 | — | 220 | +0.049 | 0.638 | 0.8322 |
| MHC-I/APM | 30 | — | 21 | -0.033 | 0.8322 | 0.8322 |
| TJ | 30 | — | 188 | +0.068 | 0.004934 | 0.02467 |
| keratin | 30 | — | 16 | +0.448 | 0.03883 | 0.09708 |
| chemokine | 30 | — | 25 | -0.034 | 0.7873 | 0.8322 |

Per-cohort family scores (GSE189357 Q4 vs Q1 skipped when n_Q4 < 3):

| cohort | family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---|---:|---|---:|---:|---:|---:|
| GSE131907 | IFN | 11 | 6 / 5 | 220 | +0.004 | 0.9896 | 0.9896 |
| GSE131907 | MHC-I/APM | 11 | 6 / 5 | 21 | -0.150 | 0.7933 | 0.9896 |
| GSE131907 | TJ | 11 | 6 / 5 | 188 | +0.205 | 0.01521 | 0.07604 |
| GSE131907 | keratin | 11 | 6 / 5 | 16 | +0.557 | 0.578 | 0.9634 |
| GSE131907 | chemokine | 11 | 6 / 5 | 25 | -0.416 | 0.07948 | 0.1987 |

Per-cohort continuous family scores (GSE189357 n=9 is thin):

| cohort | family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---|---:|---|---:|---:|---:|---:|
| GSE131907 | IFN | 21 | — | 220 | +0.077 | 0.5531 | 0.9218 |
| GSE131907 | MHC-I/APM | 21 | — | 21 | -0.025 | 0.9 | 0.9806 |
| GSE131907 | TJ | 21 | — | 188 | +0.090 | 0.001509 | 0.007546 |
| GSE131907 | keratin | 21 | — | 16 | +0.487 | 0.08319 | 0.208 |
| GSE131907 | chemokine | 21 | — | 25 | -0.003 | 0.9806 | 0.9806 |
| GSE189357 | IFN | 9 | — | 220 | -0.081 | 0.4571 | 0.6943 |
| GSE189357 | MHC-I/APM | 9 | — | 21 | -0.060 | 0.7181 | 0.7181 |
| GSE189357 | TJ | 9 | — | 188 | -0.045 | 0.05522 | 0.2761 |
| GSE189357 | keratin | 9 | — | 16 | +0.157 | 0.1267 | 0.3167 |
| GSE189357 | chemokine | 9 | — | 25 | -0.152 | 0.5554 | 0.6943 |

CLDN4 itself (held out of TJ) combined Q4 vs Q1: logFC=+2.351, p=0.07913, n_Q1=9, n_Q4=7 — direction check on the split gene.

Directions on the combined Q4 vs Q1 split match the thesis (IFN / MHC-I/APM /
chemokine down, TJ / keratin up). Family-score p is thin except chemokine.
TJ is a trend on the binary split and is the family that is up on the
continuous n=30 score (GSE131907 carries that TJ up). GSE189357 continuous
TJ is slightly down — an honest thin-n (n=9) limitation, not a second
independent TJ-up call. IFN is null on GSE131907 alone and flips slightly
positive on the continuous combined score. Keratin is up but not a claim.

## Gene-level family members (combined Q4 vs Q1)

Positive logFC = higher in CLDN4-high.

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---|
| IFN | 215 | 24 (11/13) | 0 | -0.021 | CSF2RB (-0.609, 0.003179, 0.311) |
| MHC-I/APM | 21 | 0 (0/0) | 0 | -0.190 | TAP2 (-0.986, 0.05545, 0.5186) |
| TJ | 183 | 23 (19/4) | 0 | +0.106 | ZEB1 (-1.505, 0.0003417, 0.2564) |
| keratin | 15 | 2 (2/0) | 0 | +0.619 | KRT10 (+0.931, 0.03898, 0.4762) |
| chemokine | 22 | 7 (0/7) | 0 | -0.638 | CXCR4 (-2.115, 0.005107, 0.3283) |

Headline genes (combined Q4 vs Q1, family members only, lowest p):

| family | gene | n_Q1 | n_Q4 | logFC | p | FDR |
|---|---|---:|---:|---:|---:|---:|
| TJ | ZEB1 | 9 | 7 | -1.505 | 0.0003417 | 0.2564 |
| TJ | ARL2 | 9 | 7 | -0.771 | 0.001796 | 0.311 |
| TJ | TBCD | 9 | 7 | +0.972 | 0.003062 | 0.311 |
| IFN | CSF2RB | 9 | 7 | -0.609 | 0.003179 | 0.311 |
| IFN | IRF8 | 9 | 7 | -0.936 | 0.00398 | 0.311 |
| chemokine | CXCR4 | 9 | 7 | -2.115 | 0.005107 | 0.3283 |
| TJ | DLG3 | 9 | 7 | +1.149 | 0.007236 | 0.3603 |
| TJ | CDC42 | 9 | 7 | +0.477 | 0.007798 | 0.3603 |
| IFN | PLA2G4A | 9 | 7 | +2.206 | 0.007904 | 0.3603 |
| TJ | CCND1 | 9 | 7 | +1.469 | 0.007908 | 0.3603 |
| IFN|chemokine | CCL5 | 9 | 7 | -1.701 | 0.01035 | 0.3949 |
| IFN | GZMA | 9 | 7 | -1.569 | 0.01142 | 0.4045 |

## What this is not

- Not CellChat / LIANA / NicheNet and not a T/NK infiltrate test.
- Not a re-audit of PR #459 ρ=−0.542 / Q4 r=−0.619.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not a dual-high TACSTD2×CLDN4 score.
- Not GSE148071 and not a merge beyond GSE131907+GSE189357.
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
- `figures/` — volcano, family heatmap, key-gene boxes, family-score boxes, forests, median-logFC bars, CLDN4 strip, n bars

Reproduce:

```bash
python3 methods/pair_131907_189357_malig_ifn_de_cldn4/analyze.py
```
