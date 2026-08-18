# Concordant-4 CLDN4-only — T/NK + malignant IFN/MHC

ADDITIVE. **CLDN4-only.** No TACSTD2∩CLDN4 dual-high. Not a mega-merge.
The four sets that already point the same way: **GSE123902 + GSE131907 +
GSE205335 + GSE189357**. Not GSE148071, GSE127465, GSE207422, GSE154826,
or CD45+/T-only extracts.

Thesis (already correct; not re-derived): CLDN4-high malignant cells have
lower own IFN/MHC-I, and patients have lower T/NK. CLDN4-low/KD opens IFN/MHC.
Matching extras here are IFN/MHC **DOWN** in CLDN4-high. Do not sell a
DoRothEA IFN-up-in-high as the KD direction.

The **new number is the four-set pooled n**. Pairwise rhos from PR #459 / #320
are comparison rows only and were not re-audited.

Patient / donor / sample is the unit. Do not quote cell counts as n.
GSE123902 = donor. GSE131907 = sample (tumor-bearing). GSE205335 = patient
(RECIST not required). GSE189357 = patient. p-values are descriptive.

## 1. T/NK — four-set pooled n

Primary score = malignant CLDN4 **%pos**. Mean is the matching extra.
Pooling = DerSimonian–Laird on Fisher-z of the four cohort Spearmans.
Q4 vs Q1 primary = **within-cohort quartiles stacked**, then Mann–Whitney
on T/NK (rank-biserial r). That uses the same Q labels as the DE.

| score | k | N | ρ (p, I², 95% CI) | stacked Q4 vs Q1 r (n_Q1/n_Q4, p) |
|---|---:|---:|---|---|
| %pos | 4 | 65 | -0.531 (1.65e-05, I²=0.0%, -0.697 to -0.312) | -0.724 (19/16, 0.0002879) |
| mean | 4 | 65 | -0.403 (0.00186, I²=0.0%, -0.602 to -0.157) | -0.592 (19/16, 0.00304) |

Honest N = 13 donors + 21 samples + 22 patients + 9 patients = **65 units**.

Within-cohort rank-then-qcut Q4 (combo-enum method) on %pos: r=-0.724, n_Q1/n_Q4=19/16, p=0.0002879.

### Leave-one-dataset-out (%pos)

| dropped | remaining | N | ρ (p, I²) | stacked Q4 vs Q1 r (n_Q1/n_Q4, p) |
|---|---|---:|---|---|
| GSE123902 | GSE131907+GSE205335+GSE189357 | 52 | -0.497 (0.0003493, I²=0.0%) | -0.682 (15/13, 0.002363) |
| GSE131907 | GSE123902+GSE205335+GSE189357 | 44 | -0.536 (0.0004013, I²=0.0%) | -0.804 (13/11, 0.0009587) |
| GSE205335 | GSE123902+GSE131907+GSE189357 | 43 | -0.580 (0.0001139, I²=0.0%) | -0.723 (13/10, 0.003929) |
| GSE189357 | GSE123902+GSE131907+GSE205335 | 56 | -0.522 (7.24e-05, I²=0.0%) | -0.723 (16/14, 0.0008186) |

LOO that drops GSE189357 recovers the PR #459 triple membership (n=56). That
is a robustness row of this four-set, not a re-audit of the triple rho.

### Singles (context for the pool; not the new number)

| cohort | unit | n | ρ | p | Q4 r (n_Q1/n_Q4, p) |
|---|---|---:|---:|---:|---|
| GSE123902 | donor | 13 | -0.659 | 0.01423 | -1.000 (4/3, 0.05714) |
| GSE131907 | sample | 21 | -0.522 | 0.0152 | -0.600 (6/5, 0.1255) |
| GSE205335 | patient | 22 | -0.435 | 0.04286 | -0.778 (6/6, 0.02597) |
| GSE189357 | patient | 9 | -0.600 | 0.08762 | — (nan/nan, —) |

### Comparison rows (PR #459 / #320; not re-audited)

| source | combo | analysis | N | effect | p |
|---|---|---|---:|---|---|
| PR320_given | GSE131907+GSE205335 | q4q1 pct | 23 | r=-0.705 | 0.000301 |
| PR459_given | GSE123902+GSE189357 | spearman pct | 22 | ρ=-0.638 | 0.003 |
| PR459_given | GSE123902+GSE131907 | spearman pct | 34 | ρ=-0.575 | 0.001 |
| PR459_given | GSE131907+GSE189357 | spearman pct | 30 | ρ=-0.542 | 0.003 |
| PR459_given | GSE123902+GSE205335 | spearman pct | 35 | ρ=-0.522 | 0.002 |
| PR459_given | GSE131907+GSE205335 | spearman pct | 43 | ρ=-0.479 | 0.002 |
| PR459_given | GSE189357+GSE205335 | spearman pct | 31 | ρ=-0.478 | 0.009 |

### GSE131907 nLung sensitivity (not primary)

Primary T/NK uses tumor-bearing sites only (PE, mBrain, mLN, tL/B, tLung;
n_mal≥20). All 11 nLung samples have **n_malignant=0**, so they
cannot join a malignant-CLDN4 vs T/NK test.

Epithelial CLDN4 %pos vs T/NK (sensitivity only):

| cut | n | ρ | p |
|---|---:|---:|---:|
| tumor-bearing epi (n_epi≥20) | 21 | -0.542 | 0.01123 |
| tumor-bearing epi + nLung | 32 | -0.415 | 0.01831 |
| nLung epi only | 11 | 0.209 | 0.5372 |

## 2. Tumor-cell-intrinsic — malignant Q4 vs Q1

**Method: patient-pseudobulk OLS on log2(TMM-CPM+1).** Not muscat.
muscat `pbDS` is the same collapse (patient × cell-type UMI-sum → bulk DE).
Within-cohort malignant CLDN4 %pos Q4 vs Q1, then stacked with cohort
covariates. Positive logFC = higher in CLDN4-high.

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | the four concordant sets | not the 7-pool; not +GSE148071 |
| Units | %pos N=65 | mixed donor / sample / patient |
| Split | within-cohort %pos Q4 vs Q1 | mid quartiles unused in binary DE |
| Malignant | marker-malig UMI-sum (123902, 189357); author-malig (131907, 205335) | marker gate ≠ author CNV |
| Model | ~ cohort + CLDN4_Q4 | small n; no muscat mixed model |
| Families | IFN Hallmark α∪γ; MHC-I/APM custom; chemokine panel; TJ KEGG/GO (CLDN4 held out) | chemokine panel is compact |

Honest DE n (units in the count matrices):

- GSE123902: n=13 donors. Q1=4 Q4=3.
- GSE131907: n=21 samples. Q1=6 Q4=5.
- GSE205335: n=22 patients. Q1=6 Q4=6. Out of malignant DE: P4001.
- GSE189357: n=9 patients. Q1=3 Q4=2.

| contrast | n_low / n_high | n_genes | note |
|---|---|---:|---|
| Q4 vs Q1 stacked | 18 / 16 | 15779 | cohort covariates; malignant only |
| continuous stacked | n=64 | 16419 | CLDN4 %pos z |

Do not quote N=65 as the DE n: stacked Q4 vs Q1 uses only the tails
that are in the count matrices (n=34).

### Family scores (stacked Q4 vs Q1)

Family score = mean log2(TMM-CPM+1) of family genes present.
Expected under the thesis: IFN down, MHC-I/APM down, chemokine down, TJ up.

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 34 | 18 / 16 | 221 | -0.584 | 0.002434 | 0.004869 |
| MHC-I/APM | 34 | 18 / 16 | 21 | -0.779 | 0.009883 | 0.01318 |
| TJ | 34 | 18 / 16 | 194 | +0.043 | 0.5386 | 0.5386 |
| chemokine | 34 | 18 / 16 | 25 | -0.964 | 7.73e-05 | 0.0003093 |

Continuous family-score DE (units in the count matrices):

| family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---:|---|---:|---:|---:|---:|
| IFN | 64 | — | 221 | -0.161 | 0.04215 | 0.0562 |
| MHC-I/APM | 64 | — | 21 | -0.273 | 0.01926 | 0.03851 |
| TJ | 64 | — | 194 | +0.037 | 0.168 | 0.168 |
| chemokine | 64 | — | 25 | -0.307 | 0.003142 | 0.01257 |

Per-cohort family scores (within-cohort Q4 vs Q1; GSE123902 tails are thin).
GSE189357 Q4 n=2: single-cohort binary family DE skipped (need ≥3 each);
TD6 and TD9 stay in the stacked n=34. P4001 is Q1 on the T/NK vector
(stacked T/NK n_Q1=19) but is not in the malignant UMI-sum (DE n_Q1=18).

| cohort | family | n | n_Q1 / n_Q4 | n_genes | logFC | p | FDR |
|---|---|---:|---|---:|---:|---:|---:|
| GSE123902 | IFN | 7 | 4 / 3 | 221 | -0.679 | 0.1022 | 0.1443 |
| GSE123902 | MHC-I/APM | 7 | 4 / 3 | 21 | -0.751 | 0.1082 | 0.1443 |
| GSE123902 | TJ | 7 | 4 / 3 | 194 | +0.129 | 0.62 | 0.62 |
| GSE123902 | chemokine | 7 | 4 / 3 | 25 | -1.355 | 0.1081 | 0.1443 |
| GSE131907 | IFN | 11 | 6 / 5 | 221 | -0.011 | 0.9708 | 0.9708 |
| GSE131907 | MHC-I/APM | 11 | 6 / 5 | 21 | -0.168 | 0.7643 | 0.9708 |
| GSE131907 | TJ | 11 | 6 / 5 | 194 | +0.185 | 0.02377 | 0.09509 |
| GSE131907 | chemokine | 11 | 6 / 5 | 25 | -0.424 | 0.07505 | 0.1501 |
| GSE205335 | IFN | 11 | 5 / 6 | 221 | -1.243 | 0.001759 | 0.003519 |
| GSE205335 | MHC-I/APM | 11 | 5 / 6 | 21 | -1.605 | 0.01607 | 0.02143 |
| GSE205335 | TJ | 11 | 5 / 6 | 194 | -0.115 | 0.3682 | 0.3682 |
| GSE205335 | chemokine | 11 | 5 / 6 | 25 | -1.377 | 0.001736 | 0.003519 |

CLDN4 itself (held out of TJ) stacked Q4 vs Q1: logFC=+1.678, p=0.0123, n_Q1=18, n_Q4=16 — direction check on the split gene.

### Gene-level family members (stacked Q4 vs Q1)

| family | n_tested | n p<0.05 (up/down) | n FDR<0.05 | median logFC | mean logFC | top gene (logFC, p, FDR) |
|---|---:|---|---:|---:|---:|---|
| IFN | 221 | 77 (3/74) | 12 | -0.513 | -0.542 | CCL5 (-2.457, 4.65e-06, 0.003859) |
| MHC-I/APM | 21 | 9 (0/9) | 2 | -0.686 | -0.729 | TAP2 (-1.506, 0.0005756, 0.03423) |
| TJ | 190 | 40 (30/10) | 6 | +0.109 | +0.084 | TBCD (+0.963, 9.17e-06, 0.005789) |
| chemokine | 24 | 13 (0/13) | 5 | -0.923 | -0.965 | CCL5 (-2.457, 4.65e-06, 0.003859) |

Headline genes (family members only, lowest p):

| family | gene | n_Q1 | n_Q4 | logFC | p | FDR |
|---|---|---:|---:|---:|---:|---:|
| IFN|chemokine | CCL5 | 18 | 16 | -2.457 | 4.65e-06 | 0.003859 |
| IFN | GZMA | 18 | 16 | -2.241 | 6.38e-06 | 0.004577 |
| TJ | TBCD | 18 | 16 | +0.963 | 9.17e-06 | 0.005789 |
| chemokine | XCL2 | 18 | 16 | -1.687 | 6.55e-05 | 0.01326 |
| IFN | GBP2 | 18 | 16 | -2.347 | 7.56e-05 | 0.01326 |
| TJ | CLDN3 | 18 | 16 | +1.885 | 0.0001132 | 0.01713 |
| IFN | CSF2RB | 18 | 16 | -1.112 | 0.0001212 | 0.01787 |
| IFN | GBP4 | 18 | 16 | -2.155 | 0.0001509 | 0.01905 |
| IFN | CD69 | 18 | 16 | -2.286 | 0.0001723 | 0.02039 |
| chemokine | CCL3 | 18 | 16 | -1.984 | 0.000311 | 0.02589 |
| TJ | PECAM1 | 18 | 16 | -1.826 | 0.0003635 | 0.02757 |
| chemokine | CCL4 | 18 | 16 | -1.968 | 0.0005091 | 0.03347 |

## What this is not

- Not a mega-merge of every lung scRNA set.
- Not GSE148071 / GSE127465 / GSE207422 / GSE154826 / CD45+ or T-only.
- Not a dual-high TACSTD2×CLDN4 score.
- Not a re-audit of PR #459 / #320 pairwise rhos.
- Not muscat mixed-model DE and not a cell-level Wilcoxon.
- Not DoRothEA / VIPER IFN activity sold as the KD direction.
- Not evidence that CLDN4 *causes* the T/NK or IFN/MHC change.
- Genome-wide FDR on these n is expected to be thin; family scores and
  gene-count direction (up/down) are the claim.

## Files

- `tables/tnk_pooled.tsv` — **headline T/NK table** (honest n)
- `tables/tnk_loo.tsv` — leave-one-dataset-out
- `tables/tnk_units.tsv` / `tables/tnk_singles.tsv`
- `tables/tnk_nlung_sensitivity.tsv`
- `tables/tnk_comparison_pr459_320.tsv` — given pairwise rows
- `tables/family_de.tsv` / `tables/family_summary.tsv` — **headline IFN/MHC table**
- `tables/de_q4q1_combined_families.tsv` / `tables/de_families.tsv` / `tables/de_all.tsv`
- `tables/n_honest.tsv` / `tables/sample_inventory.tsv`
- `figures/` — T/NK scatter / forests / Q4 box; malignant volcano / heatmap / family forests

Reproduce:

```bash
python3 methods/concordant4_123902_131907_205335_189357_cldn4/analyze.py
```
