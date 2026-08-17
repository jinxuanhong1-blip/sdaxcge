# FINDING — CLDN4-only ADC+ICI quadrant (no dual-high)

**Additive only. CLDN4 is the only gate. TACSTD2 is never a filter.**

This slice counts, on four public ICI-adjacent lung bulks, how many **CLDN4 Q4** tumors also sit on **IFN-γ Q4 or CD274 Q4** (combo-shaped) versus how many sit on **IFN-low (IFN Q1)**. It does not test an ADC arm. It does not invent dual-high.

**Verdict.** On the two cohorts with usable Q4 n, a real slice of CLDN4-high tumors is already IFN-γ-high or CD274-high, and a smaller slice is IFN-low. **GSE285029** (n=234): CLDN4 Q4 = **59**; combo **32**; IFN Q1 **9**. Combo is IFN-driven (IFN Q4 overlap OR=3.05, p=8.5×10⁻⁴); CD274 Q4 overlap is NS (OR=1.60, p=0.17). **GSE218989** (n=355): CLDN4 Q4 = **89**; combo **44**; IFN Q1 **22** (5 of those 22 are CD274 Q4, so they sit in both bins). Combo is CD274-driven (CD274 Q4 overlap OR=3.37, p=8.0×10⁻⁶); IFN Q4 overlap is chance (OR=1.23, p=0.48). **GSE126044** (n=16, Q4=4) and **GSE166449** (n=22, Q4=6) are honest small-n counts only — 1 vs 2 and 1 vs 1. Do not pool the four rows.

---

## Quadrant table (primary)

Equal-count quartiles: `rank(method='first')` then `qcut(4)`. Combo = CLDN4 Q4 ∩ (IFN-compact Q4 ∪ CD274 Q4). IFN-low = CLDN4 Q4 ∩ IFN-compact Q1. The two bins are **not** a partition of CLDN4 Q4 (mid-IFN Q4 tumors are neither; IFN Q1 ∩ CD274 Q4 is both).

| cohort | n | CLDN4 Q4 | ∩ IFN-γ Q4 | ∩ CD274 Q4 | ∩ (IFN-γ **or** CD274 Q4) | ∩ IFN Q1 | ∩ IFN ≤ median | overlap (combo ∩ IFN Q1) | Q4 neither |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GSE285029 | **234** | **59** | 25 | 19 | **32** | **9** | 22 | 0 | 18 |
| GSE218989 | **355** | **89** | 25 | 39 | **44** | **22** | 41 | **5** | 28 |
| GSE126044 | **16** | **4** | 1 | 0 | **1** | **2** | 3 | 0 | 1 |
| GSE166449 | **22** | **6** | 1 | 1 | **1** | **1** | 3 | 0 | 4 |

Machine table: `tables/quadrant_table.tsv`. Per-sample calls: `tables/*_sample_quadrants.tsv`.

---

## Honest n

| cohort | matrix n | CLDN4 complete | Q1 / Q2 / Q3 / Q4 | Q4 used | TACSTD2 gate | response on GEO |
|---|---:|---:|---|---:|---|---|
| GSE285029 | 234 | 234 | 59 / 58 / 58 / 59 | **59** | **none** | **no** |
| GSE218989 | 355 | 355 | 89 / 89 / 88 / 89 | **89** | **none** | 168 R / 187 NR |
| GSE126044 | 16 | 16 | 4 / 4 / 4 / 4 | **4** | **none** | 5 R / 11 NR |
| GSE166449 | 22 | 22 | 6 / 5 / 5 / 6 | **6** | **none** | 7 R / 15 NR |

Do not write n=234 or n=355 as the quadrant n. The claim unit is **CLDN4 Q4**. GSE126044 Q4 n=4 and GSE166449 Q4 n=6 cannot support a rate.

---

## Methods (fixed before looking at the counts)

- **Anchor:** CLDN4 only. No TACSTD2, no dual-high, no two-gene mean.
- **IFN-γ primary:** unweighted mean of per-gene z-scores of the A11 compact set `IFNG STAT1 IRF1 CXCL9 CXCL10 CXCL11 IDO1 GBP1` (8/8 present on every matrix used here).
- **CD274:** single gene, same transform as CLDN4 on that matrix.
- **IFN-low:** IFN-compact **Q1** (symmetric with Q4). IFN ≤ median is extra only.
- **Quartiles:** complete-case per cohort; `rank(method='first')` + `qcut(4)` so ties still fill four bins.
- **Enrichment (supporting, not the count):** Fisher exact for CLDN4 Q4 ∩ feature Q4 vs independence (expected ~25% of Q4). Spearman on the full matrix is supporting only.
- **Transforms (already used in this repo; not re-tuned):**
  - GSE285029: harvested A11 scores (`log2(pmax(x,0)+1)` of the author WTS matrix).
  - GSE218989: `log2(TPM+1)` from the public protein-coding TPM matrix. Harvested CLDN4 / IFNG / Ayers6; CD274 and IFN-compact computed from an 11-gene extract of the same file.
  - GSE126044: `log2(CPM+1)` from the public count matrix. Harvested CLDN4 matches the recompute.
  - GSE166449: deposited matrix used **as-is** (filename says TPM; values already match the harvested CLDN4 column — do not re-log).
- **Not done:** TACSTD2 gate, dual-high, cutpoint search, pooling, meta-OR, purity residual, ADC-response test.

Sources: `harvested/PROVENANCE.md`.

---

## Per-cohort reading

### GSE285029 — Koh 2025 pre-ICI NSCLC WTS (n=234)

CLDN4 Q4 = 59. Combo **32 / 59** (54%). IFN Q1 **9 / 59** (15%). Below-median IFN **22 / 59**.

| supporting | n | result |
|---|---:|---|
| CLDN4 vs IFN-compact Spearman | 234 | ρ=+0.204, p=0.0018 |
| CLDN4 vs CD274 Spearman | 234 | ρ=+0.167, p=0.011 |
| CLDN4 Q4 ∩ IFN Q4 | 25 / 59 | OR=3.05, Fisher p=8.5×10⁻⁴ |
| CLDN4 Q4 ∩ CD274 Q4 | 19 / 59 | OR=1.60, Fisher p=0.17 |

This is the same direction as the A11 Q4-overlap file on this matrix (IFN yes, CD274 Q4 no). GEO has no response / histology / purity. Extra combination-rationale only: about half of CLDN4-high tumors already sit on a high IFN-γ program; CD274-Q4 co-occupancy is not enriched.

### GSE218989 — Kang 2024 SMC–KAIST ICI TPM (n=355)

CLDN4 Q4 = 89. Combo **44 / 89** (49%). IFN Q1 **22 / 89** (25%). **Five** of the 22 IFN-Q1 tumors are CD274 Q4 — they are counted in both bins and listed in `tables/GSE218989_sample_quadrants.tsv`.

| supporting | n | result |
|---|---:|---|
| CLDN4 vs IFN-compact Spearman | 355 | ρ=+0.055, p=0.30 |
| CLDN4 vs CD274 Spearman | 355 | ρ=+0.148, p=0.0053 |
| CLDN4 Q4 ∩ IFN Q4 | 25 / 89 | OR=1.23, Fisher p=0.48 |
| CLDN4 Q4 ∩ CD274 Q4 | 39 / 89 | OR=3.37, Fisher p=8.0×10⁻⁶ |

IFN-γ Q4 among CLDN4 Q4 is what independence predicts. The combo count here is **CD274 Q4**, not IFN-γ Q4. That matches the already-published near-zero CLDN4–Ayers6 Spearman on this table. Histology is not deposited. Response exists (168 / 187) but is **not** the quadrant claim.

### GSE126044 — Cho pre-anti-PD-1 counts (n=16)

CLDN4 Q4 = **4**. Combo **1**. IFN Q1 **2**. Spearman vs IFN ρ=−0.49 (p=0.057); vs CD274 ρ=−0.64 (p=0.008). Opposite sign of GSE285029, n=16. All four Q4 tumors are GEO non-responders; 5/5 FFPE are NR so response is confounded with sample type. **Count only.**

### GSE166449 — Lee 2021 LUAD pembrolizumab TPM (n=22)

CLDN4 Q4 = **6**. Combo **1** (a non-responder). IFN Q1 **1** (a non-responder). Spearman vs IFN ρ=−0.14 (p=0.53); vs CD274 ρ=−0.02 (p=0.93). **Count only.**

---

## Extra scatter

Quartile lines on CLDN4 vs IFN-compact and CLDN4 vs CD274. Orange = CLDN4 Q4 ∩ (IFN Q4 ∪ CD274 Q4). Dark blue = CLDN4 Q4 ∩ IFN Q1 (and not combo). Light blue = other CLDN4 Q4. Grey = not CLDN4 Q4.

| file | n |
|---|---:|
| `figures/GSE285029_extra_scatter.png` | 234 |
| `figures/GSE218989_extra_scatter.png` | 355 |
| `figures/GSE126044_extra_scatter.png` | 16 |
| `figures/GSE166449_extra_scatter.png` | 22 |
| `figures/quadrant_counts.png` | four bars |

---

## What this is not

- Not dual-high. TACSTD2 is not used.
- Not an ADC trial and not an ADC+ICI efficacy test.
- Not a claim that CLDN4-high tumors are generally IFN-hot: GSE218989 IFN overlap is chance; GSE285029 IFN overlap is enriched; the two small sets go the other way or are null.
- Not a pooled rate. The four Q4 denominators are 59, 89, 4, and 6.

Reproduce: `python3 methods/cldn4_adc_ici_quadrant/analyze.py`
