# Claim A1 — ESTIMATE / ABSOLUTE recompute

**Claim (user):** TACSTD2 is inverse to immune / cytotoxic / exhaustion in
TCGA LUAD + LUSC + OncoSG, **still negative after purity**, **n ≈ 1031**.

**This update** recomputes the purity-partial Spearman under published
**ESTIMATE** (Aran et al. Nat Commun 2015), **ABSOLUTE** (GDC PanCanAtlas),
**CPE** (Aran consensus of ESTIMATE+ABSOLUTE+LUMP+IHC), and
**ESTIMATE-or-ABSOLUTE** (ESTIMATE preferred, ABSOLUTE fallback).
OncoSG has **no public ESTIMATE table**; its clinical `PURITY` attribute is
used wherever a non-ESTIMATE method includes OncoSG. That limitation is not
hidden.

Public data only. Code: `code/run_claim_A1.py`.

---

## Verdict

| Piece of the claim | Result |
|---|---|
| Still **negative** after purity | **MATCH** — every signature, every purity method, every cohort with data |
| Significant (p < 0.05) after purity | **MATCH** in all pooled analyses |
| **n ≈ 1031** for LUAD+LUSC+**OncoSG** | **MISMATCH** — three-cohort ABSOLUTE / ESTIMATE-or-ABSOLUTE n = **1145 / 1162** |
| **n ≈ 1031** for LUAD+LUSC ESTIMATE (no OncoSG) | **MATCH (closest)** — Firehose ESTIMATE **n = 1014** (Δ −17); Firehose ESTIMATE-or-ABSOLUTE **n = 1017** (Δ −14) |

The inverse association **replicates**. The reported n ≈ 1031 is the
**TCGA LUAD+LUSC ESTIMATE/CPE-complete set**, not the three-cohort pool.
Adding OncoSG (169 tumors with clinical purity) moves n to ~1145–1186.

---

## Primary request: LUAD+LUSC+OncoSG, ESTIMATE/ABSOLUTE partial Spearman

PanCanAtlas LUAD + LUSC + OncoSG. Purity = ESTIMATE if present, else GDC
ABSOLUTE; OncoSG = clinical PURITY.

| Signature | n | ρ (partial) | p | vs claim |
|---|---|---|---|---|
| Immune (T-cell effector) | 1162 | **−0.187** | 1.3e-10 | sign MATCH; n 1162 ≠ 1031 |
| Cytotoxic | 1162 | **−0.161** | 3.7e-8 | sign MATCH; n 1162 ≠ 1031 |
| Exhaustion / checkpoint | 1162 | **−0.165** | 1.6e-8 | sign MATCH; n 1162 ≠ 1031 |

Same pool, **ABSOLUTE only** (GDC + OncoSG clinical):

| Signature | n | ρ (partial) | p |
|---|---|---|---|
| Immune | 1145 | **−0.213** | 3.8e-13 |
| Cytotoxic | 1145 | **−0.188** | 1.6e-10 |
| Exhaustion | 1145 | **−0.194** | 3.9e-11 |

Unadjusted Spearman on the same 1163 expression samples is already negative
(−0.200 / −0.177 / −0.175). Purity adjustment does not create the effect.

---

## Closest n to 1031: Firehose LUAD+LUSC (OncoSG not in the n)

No public ESTIMATE exists for OncoSG, so an ESTIMATE-only analysis is
TCGA-only. Firehose RNA-seq + Aran ESTIMATE is the published ESTIMATE
complete set for LUAD+LUSC:

| Method | n | Immune ρ (p) | Cytotoxic ρ (p) | Exhaustion ρ (p) | Δ vs 1031 |
|---|---|---|---|---|---|
| ESTIMATE | **1014** | −0.167 (8.5e-8) | −0.126 (6.1e-5) | −0.142 (5.6e-6) | −17 |
| ESTIMATE or ABSOLUTE | **1017** | −0.160 (2.9e-7) | −0.118 (1.6e-4) | −0.131 (3.0e-5) | −14 |
| CPE | **1016** | −0.192 (7.1e-10) | −0.152 (1.1e-6) | −0.168 (6.8e-8) | −15 |
| ABSOLUTE (GDC) | **997** | −0.200 (1.8e-10) | −0.165 (1.6e-7) | −0.177 (2.0e-8) | −34 |

PanCanAtlas LUAD+LUSC is the same story at slightly smaller n (ESTIMATE 990,
CPE 994, ABSOLUTE 976).

Aran Supp Data 1, before any RNA-seq join, has **1014** LUAD+LUSC ESTIMATE
values and **1034** CPE values. 1031 sits on that CPE/ESTIMATE universe;
it is not a three-cohort count.

---

## Per-cohort (PanCan + OncoSG), ESTIMATE vs ABSOLUTE

| Cohort | Method | Immune ρ (p, n) | Cytotoxic | Exhaustion |
|---|---|---|---|---|
| TCGA LUAD | ESTIMATE | −0.151 (6.8e-4, 506) | −0.163 (2.3e-4, 506) | −0.129 (3.8e-3, 506) |
| TCGA LUAD | ABSOLUTE | −0.120 (7.4e-3, 497) | −0.127 (4.7e-3, 497) | −0.062 (0.17, 497) ns |
| TCGA LUSC | ESTIMATE | −0.208 (4.0e-6, 484) | −0.101 (0.027, 484) | −0.188 (3.1e-5, 484) |
| TCGA LUSC | ABSOLUTE | −0.271 (1.7e-9, 479) | −0.201 (9.3e-6, 479) | −0.266 (3.5e-9, 479) |
| OncoSG LUAD | clinical PURITY | −0.318 (2.7e-5, 169) | −0.363 (1.3e-6, 169) | −0.389 (1.9e-7, 169) |
| OncoSG LUAD | ESTIMATE | — | — | n = 0 (no public ESTIMATE) |

Every point estimate with data is negative. The only non-significant cell is
TCGA-LUAD exhaustion under ABSOLUTE (ρ −0.062, p = 0.17); under ESTIMATE the
same cell is significant (ρ −0.129, p = 3.8e-3).

---

## Match / mismatch (honest)

1. **Direction after purity — MATCH.** Not one pooled or per-cohort partial
   correlation with n ≥ 6 is positive.
2. **n ≈ 1031 for “TCGA+OncoSG” — MISMATCH.** Three-cohort n is 1145
   (ABSOLUTE) or 1162 (ESTIMATE-or-ABSOLUTE). That is ~110–130 samples above
   1031, which is the OncoSG RNA+purity set (169) minus TCGA samples that
   lack the chosen purity call.
3. **n ≈ 1031 as TCGA LUAD+LUSC ESTIMATE/CPE — MATCH.** Firehose ESTIMATE
   n = 1014; CPE n = 1016; ESTIMATE-or-ABSOLUTE n = 1017. Aran CPE universe
   (no RNA filter) = 1034.
4. **OncoSG ESTIMATE — cannot compute from public tables.** cBioPortal
   OncoSG has only z-score profiles and a clinical PURITY column. Official
   ESTIMATE needs the Yoshihara stromal/immune transcriptome calibration;
   we do not invent it.

See `match_mismatch.csv` for every pool × method × signature.

---

## Reproduce

```
pip install -r requirements.txt
python code/run_claim_A1.py
python code/make_figures.py
```

Purity tables (already in `data/`):

- `data/abs_purity.txt` — GDC PanCanAtlas ABSOLUTE
  (`https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5`)
- `data/aran_purity.tsv` — Aran et al. 2015 Supp Data 1 (ESTIMATE, ABSOLUTE, CPE)

## Files

- `pooled_correlations.csv` — all pools × purity methods
- `per_cohort_correlations.csv` — per cohort × method
- `match_mismatch.csv` — n / ρ / p / verdict vs n ≈ 1031
- `meta_analysis.csv` — Fisher-z meta of PanCan LUAD, LUSC, OncoSG
- `individual_gene_checks.csv` — CD8A, GZMA, PDCD1, CTLA4, HAVCR2, CD274
- `summary.json` — machine-readable primary numbers
- `sample_data_*.csv` — per-sample TACSTD2, signatures, all purity columns
- `forest_partial_spearman.png`, `scatter_pooled_rank.png`

## Caveats

- Bulk-tumor correlation is an association, not a causal tumor-intrinsic
  claim, even after a purity covariate.
- ESTIMATE is itself an RNA-derived stromal/immune score; partialling
  ESTIMATE out of an immune-signature correlation is a conservative test
  (it can remove true immune signal). The association remaining negative
  after ESTIMATE is therefore stronger evidence, not weaker.
- Signature scores are mean per-gene z-scores; Spearman is rank-based and
  invariant to that monotonic transform.
