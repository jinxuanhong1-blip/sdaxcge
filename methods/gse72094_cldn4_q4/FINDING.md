# Finding — GSE72094 CLDN4 Q4 vs Q1 vs CD8A / CD274

**Additive extra cut only.** Schabath / Moffitt **GSE72094** (PMID 26477306): resected **lung adenocarcinoma**, **GPL15048** HuRSTA, author IRON-normalized series matrix. All matrix samples are `source_name = lung adenocarcinoma`. Patient = array: **442 unique `patient_id` / 442 arrays** (multi-sample patients = 0).

**Already known (PR 302 / `methods/gse72094_cldn4`, not re-claimed):** continuous CLDN4 vs CD8A, n=442, Spearman ρ=**-0.219**, p=**3.35e-06**. After ESTIMATEScore residual the same pair is **not significant** (partial ρ=**-0.051**, p=**0.282**).

Two extra questions:

1. Continuous **CLDN4 vs CD274** (CD274 was not a primary row on the residual page).
2. **Q4 vs Q1** — CLDN4 quartiles on **CD8A** and **CD274**. Honest n is the quartile arms, not 442.

Reproduce: `python3 methods/gse72094_cldn4_q4/analyze.py` (GEO files cached under `/tmp/gse72094_cldn4_q4/`).

---

## Verdict

| Cut | n | Metric | Effect | p | Holds? |
|---|---:|---|---:|---:|---|
| CLDN4 vs CD8A after ESTIMATEScore (known residual) | 442 | partial Spearman ρ | -0.051 | 0.282 | **no (already NS)** |
| CLDN4 vs CD8A unadjusted (known, not extra) | 442 | Spearman ρ | -0.219 | 3.35e-06 | already known |
| **CLDN4 vs CD274** (continuous extra) | 442 | Spearman ρ | -0.034 | 0.477 | **no** |
| **CLDN4 Q4 vs Q1, endpoint CD8A** | **111 vs 111** | rank-biserial (Q4>Q1) | -0.328 | 2.47e-05 | **yes** |
| **CLDN4 Q4 vs Q1, endpoint CD274** | **111 vs 111** | rank-biserial (Q4>Q1) | -0.015 | 0.851 | **no** |

**What holds.** CLDN4-high (Q4) tumors have lower **CD8A** than CLDN4-low (Q1), n=**111 vs 111**, Δ median -0.683, rank-biserial -0.328, p=2.47e-05. That is the same unadjusted immune-low direction as the known continuous CD8A row, not a residual-independent axis.

**What does not hold.** The already-reported residual **CLDN4 vs CD8A after ESTIMATEScore remains not significant** (n=442, partial ρ=-0.051, p=0.282). This PR does not upgrade that residual. Continuous CLDN4 vs **CD274** is null at n=442 (ρ=-0.034, p=0.477). CLDN4 Q4 vs Q1 vs **CD274** is not significant at n=111 vs 111 (r=-0.015, p=0.851). Do not write n=442 for the Q4 vs Q1 tests. Those tests drop Q2+Q3 and use **111 vs 111**.

---

## Honest n

| Item | n |
|---|---:|
| GEO series / matrix columns | **442** |
| `source_name` = lung adenocarcinoma | **442** |
| Unique `patient_id` | **442** |
| Patients with >1 array | **0** |
| CLDN4 / CD8A / CD274 / TACSTD2 non-NA | **442 / 442 / 442 / 442** |
| CLDN4 Q1 / Q2 / Q3 / Q4 | **111 / 110 / 110 / 111** |
| Q4 vs Q1 used | **111 vs 111** |

Quartiles are `pd.qcut(rank(method="first"), 4)` among the 442 tumors with non-NA CLDN4. Do not write n=442 for the quartile tests.

---

## Extra continuous — CLDN4 vs CD274

Spearman on all 442 tumors. Author IRON log2.

| Endpoint | n | ρ | p | Kind |
|---|---:|---:|---:|---|
| CD8A | 442 | -0.219 | 3.35e-06 | known residual-page unadjusted |
| **CD274** | 442 | **-0.034** | **0.477** | extra |
| TACSTD2 (companion) | 442 | +0.549 | 3.33e-36 | companion only |

Supporting residual (not the Q4 extra): after harvested ESTIMATEScore, CLDN4 vs CD274 flips to partial ρ=+0.156 p=0.000987 (n=442). Unadjusted and Q4 vs Q1 CD274 stay null, so this page does not upgrade the residual sign-flip into a PD-L1 claim.

---

## Extra Q4 vs Q1

Two-sided Mann–Whitney U. Rank-biserial r = 2U/(n4 n1) − 1 (positive = Q4 higher).

### CLDN4 quartiles

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| **CD8A** | **111 / 111** | 9.476 | 10.159 | -0.683 | 4142 | **-0.328** | **2.47e-05** |
| **CD274** | **111 / 111** | 8.147 | 8.270 | -0.123 | 6070 | **-0.015** | **0.851** |
| TACSTD2 (companion) | 111 / 111 | 12.663 | 12.036 | +0.628 | 11232 | +0.823 | 3.02e-26 |

### CD8A quartiles (endpoint = CLDN4; supporting)

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| CLDN4 | **111 / 111** | 10.189 | 10.752 | -0.563 | 4004 | -0.350 | 6.61e-06 |

The reverse cut is the same sign as the known continuous CD8A row. Do not treat it as a second independent discovery.

---

## Methods (this slice)

- **Matrix:** GEO `GSE72094_series_matrix.txt.gz` (author-processed IRON + RNA-quality batch correction). Probe × sample values used as published.
- **Annotation:** GPL15048 `GeneSymbol` (NCBI platform table). First symbol if `///`. **Max-mean** probe collapse to HUGO, same rule as `scripts/gse72094_cldn4.py`.
- **CD8** = `CD8A`. **CD274** = PD-L1 RNA (not protein). TACSTD2 is a companion only.
- **Quartiles:** `pd.qcut(rank(method="first"), 4)` on the 442 non-NA CLDN4 values.
- **Continuous test:** two-sided Spearman.
- **Quartile test:** two-sided Mann–Whitney U, rank-biserial as above.
- **Residual reprint:** ESTIMATEScore from `harvested/estimate_scores.tsv` (PR 302 ssGSEA). This page does not re-fit ESTIMATE.
- Surgical LUAD. No ICI-response labels on this series.

---

## What this does not claim

- It does not re-open the PR 302 residual (already NS) as a new finding.
- It does not claim CLDN4 vs PD-L1 **protein**. CD274 is microarray RNA.
- It does not upgrade the supporting residual sign-flip (CLDN4 vs CD274 partial ρ=+0.156) into a PD-L1 claim. Unadjusted and Q4 vs Q1 CD274 are null.
- It does not treat n=442 as the Q4 vs Q1 n.
- Q4 vs Q1 CD8A is a coarsened version of the known unadjusted continuous test, not a new cohort and not a purity-residual claim.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `tables/q4_vs_q1.tsv` — MWU extra cuts
- `tables/continuous.tsv` — Spearman extra + known CD8A reprint
- `tables/n_table.tsv` — honest n
- `tables/verdict.tsv`
- `tables/samples.tsv` — per-array genes + quartiles
- `harvested/estimate_scores.tsv` — PR 302 ESTIMATEScore (residual reprint only)
- `figures/cldn4_q4q1_cd8a.png`, `cldn4_q4q1_cd274.png`
