# REWORK A1 — is the immune-cold TACSTD2 signal LUSC-driven?

**Self-contained. Public data only. Written to be read without the rest of the repo.**

**Verdict: PARTLY — feature-dependent, not a blanket LUSC effect**

GEP18 and ESTIMATE ImmuneScore are LUSC-only (LUAD is null / not negative). CD8 is negative in both histologies but stronger in LUSC. CYT is negative in both and the two rhos are not distinguishable. Pooled NSCLC TACSTD2–immune rho is therefore not a LUAD finding for GEP18/ESTIMATE, and not LUSC-only for CD8/CYT. Do not quote one pooled number.

## Why this rework exists

Claim A1 (as previously replicated on pooled TCGA NSCLC ± OncoSG) said TACSTD2
(TROP2) is negatively correlated with immune / cytotoxic / exhaustion signatures
and that the association stays negative after purity adjustment. That analysis
**pooled LUAD and LUSC**. Those are different diseases: different cells of
origin, different typical TACSTD2 levels, different typical immune infiltration.
A pooled negative rho can be:

1. a within-LUAD fact,
2. a within-LUSC fact,
3. both, or
4. a between-histology artifact (Simpson): LUSC higher TACSTD2 and colder
   immune scores, so the scatter leans negative even if each cloud is flat.

This file tests (1)–(4) on four pre-specified immune readouts, separately in
TCGA-LUAD and TCGA-LUSC, with ABSOLUTE partial Spearman.

## Analysis set

| Cohort | Primary tumors with RNA | With ABSOLUTE purity | With ESTIMATE ImmuneScore |
|---|---:|---:|---:|
| TCGA-LUAD | 515 | 502 | 515 |
| TCGA-LUSC | 502 | 493 | 501 |

Partial-correlation n (TACSTD2 + feature + ABSOLUTE) is ~502 LUAD and ~493 LUSC;
exact n is in every table cell because ESTIMATE is missing for a few samples.
Primary tumors only (`-01`). One row per 15-character barcode.

## Pre-specified features

| Name | Definition | Honest limitation |
|---|---|---|
| **CD8** | `CD8A` log2(RSEM+1) on Xena HiSeqV2 | Single gene, not a deconvolution fraction. |
| **CYT** | mean(log2 GZMA, log2 PRF1) = Rooney 2015 geometric mean on log2 data | Two genes; tracks cytotoxic mRNA, not protein or killing. |
| **GEP18** | unweighted mean of within-cohort z-scores of the 18 Ayers 2017 T-cell-inflamed GEP genes | Merck's NanoString weights are **not public**. This is the standard open surrogate, not the clinical assay. |
| **ESTIMATE** | official **ImmuneScore** from the MD Anderson RNAseqV2 tables (Yoshihara *Nat Commun* 2013) | ImmuneScore is built to track leukocyte/stromal content and is strongly (negatively) correlated with purity. That is why ABSOLUTE adjustment is not optional here. ESTIMATEScore (Immune+Stromal) is a secondary row, not the primary endpoint. |

Covariate: PanCanAtlas **ABSOLUTE** purity. Method: first-order partial Spearman
(algebraic formula, same as the original A1 script). Residual-rank partial
Spearman is a sensitivity column; the two methods agreed in sign and magnitude
here. 95% CIs are Fisher-z with variance `1/(n-4)` for partial correlations.
LUAD vs LUSC difference: Fisher z-test of two independent partial correlations.
BH-FDR is across the 8 primary partial tests (4 features × 2 histologies).

## Direct answer

**No single yes/no.** Pooled TCGA NSCLC looks immune-cold for all four
features. Split by histology, that is not one fact:

- **GEP18:** LUSC-only. LUAD ρ ≈ 0.
- **ESTIMATE ImmuneScore:** LUSC-only. LUAD is weakly *positive* and not significant.
- **CD8 (CD8A):** negative in both; significantly stronger in LUSC.
- **CYT:** negative in both; LUAD and LUSC rhos are not distinguishable.

So: the *T-cell-inflamed / ESTIMATE* half of the A1 story is LUSC-driven in
TCGA. The *cytotoxic mRNA* half is not. Quoting a pooled NSCLC rho hides that.

## Primary result — ABSOLUTE partial Spearman

| Feature | LUAD partial ρ | LUSC partial ρ | Pooled LUAD+LUSC | Fisher z p (LUAD vs LUSC) | Call |
|---|---|---|---|---|---|
| CD8 | -0.104 (p=0.020, n=502) | -0.244 (p=4.30e-08, n=493) | -0.191 (p=1.31e-09, n=995) | 0.023 | both_negative_LUSC_stronger |
| CYT | -0.116 (p=0.009, n=502) | -0.172 (p=1.27e-04, n=493) | -0.151 (p=1.78e-06, n=995) | 0.370 | both_negative_similar |
| GEP18 | -0.013 (p=0.775, n=502) | -0.226 (p=3.95e-07, n=493) | -0.152 (p=1.58e-06, n=995) | 6.38e-04 | LUSC_only_negative |
| ESTIMATE | 0.072 (p=0.109, n=502) | -0.131 (p=0.003, n=493) | -0.107 (p=7.38e-04, n=995) | 0.001 | LUSC_only_negative |

BH-FDR across the 8 primary partial tests is in `correlations.tsv`
(`partial_fdr_8tests`). LUAD CD8 and CYT remain FDR < 0.05; LUAD GEP18 and
ESTIMATE do not. All four LUSC primary tests remain FDR < 0.05.

Unadjusted Spearman (no purity) is in `correlations.tsv`. No primary sign flips.
One nuance: LUSC ESTIMATE ImmuneScore is **not** significant unadjusted
(ρ = -0.055, p = 0.216) and becomes significant only after ABSOLUTE (ρ = -0.131, p = 0.003).
That is expected — ImmuneScore is built as a purity/leukocyte composite —
and is why the partial, not the raw, number is the ESTIMATE claim.
CD8/CYT/GEP18 in LUSC are already negative before adjustment; partialling
makes them slightly more negative.

### How to read the calls

- `LUSC_only_negative`: significant negative in LUSC; LUAD not significant.
- `both_negative_LUSC_stronger`: both significant negative, and the two rhos differ (Fisher p < 0.05) with |LUSC| > |LUAD|.
- `both_negative_similar`: both significant negative, rhos not distinguishable.
- `neither_significant`: neither histology is significant after ABSOLUTE.
- `discordant_sign`: opposite significant signs — do not average these.

## Is the pooled signal a between-histology artifact?

If LUSC is TACSTD2-higher **and** immune-colder than LUAD, a pooled scatter
can look immune-cold even when each histology is internally weak. That is
checked two ways: (i) Mann–Whitney location shifts, (ii) whether |pooled ρ|
exceeds both within-histology |ρ|.

| Feature | median LUAD | median LUSC | LUSC − LUAD | Mann–Whitney p | pooled \|ρ\| > both within? |
|---|---:|---:|---:|---|---|
| TACSTD2 | 12.6 | 13.1 | 0.479 | 7.70e-11 | — |
| CD8 | 7.95 | 7.69 | -0.252 | 0.005 | no |
| CYT | 7.51 | 7.39 | -0.127 | 0.053 | no |
| GEP18 | 0.0331 | 0.0014 | -0.0317 | 0.919 | no |
| ESTIMATE | 980 | 418 | -563 | 7.81e-19 | no |
| purity | 0.45 | 0.5 | 0.05 | 1.87e-05 | — |

LUSC has **higher** median TACSTD2 than LUAD in this freeze.
LUSC is also ESTIMATE-colder (median ImmuneScore 418 vs 980, p = 7.81e-19) and slightly CD8-lower. That between-histology geometry **biases a pooled scatter toward a negative TACSTD2–immune rho**. It does not inflate |pooled ρ| past both within-histology |ρ| (Simpson flag is false for all four features): the pooled number is a blend, not a fake correlation from two flat clouds. The within-histology partial rhos are still the numbers that answer the question.

## TACSTD2 vs ABSOLUTE purity (context, not a feature)

If TACSTD2 is just an epithelial/purity gene, partialling purity should kill the
immune associations. That is not a substitute for the partial table, but it is
useful context.

| Histology | TACSTD2 vs ABSOLUTE ρ | p | n |
|---|---|---|---|
| LUAD | 0.007 | 0.869 | 502 |
| LUSC | -0.074 | 0.102 | 493 |

## Sensitivity (not used for the verdict)

| Feature | LUAD partial ρ | LUSC partial ρ | Why it exists |
|---|---|---|---|
| CD8_AB | -0.115 (p=0.010, n=502) | -0.265 (p=2.22e-09, n=493) | mean(CD8A, CD8B) instead of CD8A alone |
| GEP18_common | -0.014 (p=0.762, n=502) | -0.227 (p=3.43e-07, n=493) | GEP18 z-scored on LUAD+LUSC together, not within histology |
| ESTIMATE_Score | 0.041 (p=0.364, n=502) | -0.220 (p=8.68e-07, n=493) | official ESTIMATEScore = ImmuneScore + StromalScore |
| ESTIMATE_StromalScore | 0.007 (p=0.869, n=502) | -0.272 (p=8.67e-10, n=493) | stromal, not immune; included so ImmuneScore is not silently swapped |

Residual-rank partial Spearman (Pearson of rank residuals) is in
`correlations.tsv` as `partial_rho_residual_method`. It is a method check,
not a second discovery pass.

## Honest interpretation

1. **Headline.** GEP18 and ESTIMATE ImmuneScore are LUSC-only (LUAD is null / not negative). CD8 is negative in both histologies but stronger in LUSC. CYT is negative in both and the two rhos are not distinguishable. Pooled NSCLC TACSTD2–immune rho is therefore not a LUAD finding for GEP18/ESTIMATE, and not LUSC-only for CD8/CYT. Do not quote one pooled number.
2. **Effect sizes.** Even where p is small, |ρ| in bulk RNA is modest. A
   significant LUSC rho of −0.2 is not “immune desert because of TROP2”.
   It is a weak-to-moderate rank association in mixed tissue.
3. **Purity.** ESTIMATE ImmuneScore is almost a purity inverse. The
   ABSOLUTE-partial number is the one that is allowed to be called
   “immune” rather than “not tumor”. CD8/CYT/GEP18 are also compositionally
   entangled with purity; partialling helps, it does not prove tumor-intrinsic biology.
4. **GEP18 is not the Merck assay.** Unweighted z-mean of the 18 genes is
   the public approximation. Do not write “T-cell-inflamed GEP (NanoString)”
   as if the clinical weights were used.
5. **This is not OncoSG, not protein, not ICI response.** Original A1 also
   used OncoSG LUAD (stronger negative rhos in that cohort). This rework is
   TCGA histology only. TROP2 protein (the ADC target) was not measured.
6. **No causality.** Bulk correlation, even purity-adjusted and histology-split,
   does not say TACSTD2 excludes T cells.

## Reproduce

```
pip install -r requirements.txt
python scripts/rework_A1_histology.py
```

Downloads ~65 MB of public tables into `data/` (gitignored) on first run.

## Files

- `correlations.tsv` — unadjusted + ABSOLUTE-partial Spearman, both methods, CIs
- `luad_vs_lusc_fisher.tsv` — Fisher z tests + per-feature class + Simpson flag
- `histology_location_shifts.tsv` — Mann–Whitney LUAD vs LUSC on each score
- `sample_table.tsv` — per-sample TACSTD2, features, ABSOLUTE, ESTIMATE
- `summary.json` — machine-readable verdict
- `provenance.json` — URLs, md5, gene lists
- `figures/forest_partial_rho.png`
- `figures/scatter_by_histology.png`
- `figures/box_by_histology.png`

## Data

- Expression: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` and `TCGA.LUSC.sampleMap/HiSeqV2`.
- Purity: GDC `4f277128-f793-4354-a13d-30cc7fe9f6b5` (PanCanAtlas ABSOLUTE).
- ESTIMATE: MD Anderson official RNAseqV2 tables
  (`lung_adenocarcinoma_RNAseqV2.txt`, `lung_squamous_cell_carcinoma_RNAseqV2.txt`).

