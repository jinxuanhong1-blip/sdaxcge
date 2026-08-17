# Additive GSE182328 — CLDN4 vs CD274/HLA and Q4 vs Q1 CD8

**Extra cuts only.** Derosa / Routy **GSE182328**: advanced NSCLC tumor RNA-seq, ICI-treated, author gene-count matrix, **n=44** (one tumor per patient). GEO has Akkermansia detectable vs not. **No RECIST / DCB / PFS / MPR field on GEO.**

**Already known (leftover PR, not re-claimed here):** continuous CLDN4 vs CD8A, n=44, Spearman ρ=**−0.403**, p=**0.0067**. Same matrix, same log2(CPM+1). This page reprints that row so the extra cuts sit next to it.

Two extra questions:

1. Continuous **CLDN4 vs CD274** and **CLDN4 vs HLA** (class I genes and HLA-A/B/C mean-z).
2. **Q4 vs Q1 CD8** — CLDN4 quartiles on CD8A, and the reverse CD8A quartiles on CLDN4.

Reproduce: `python3 methods/gse182328_cldn4_extra/analyze.py` (GEO files cached under `/tmp/gse182328_cldn4_extra/`).

---

## Verdict

| Cut | n | Metric | Effect | p | Holds? |
|---|---:|---|---:|---:|---|
| CLDN4 vs CD8A (known, not extra) | 44 | Spearman ρ | **−0.403** | **0.0067** | already known |
| **CLDN4 vs CD274** | 44 | Spearman ρ | −0.169 | 0.273 | **no** |
| **CLDN4 vs HLA-A / HLA-B / HLA-C** | 44 | Spearman ρ | −0.040 / −0.121 / +0.124 | 0.80 / 0.43 / 0.42 | **no** |
| **CLDN4 vs HLA-I (HLA-A/B/C mean-z)** | 44 | Spearman ρ | +0.005 | 0.972 | **no** |
| **CLDN4 Q4 vs Q1, endpoint CD8A** | **11 vs 11** | rank-biserial (Q4>Q1) | **−0.57** | **0.026** | **yes, modest n** |
| CD8A Q4 vs Q1, endpoint CLDN4 | **11 vs 11** | rank-biserial (Q4>Q1) | −0.45 | 0.076 | same sign, not significant |

**What holds.** The leftover continuous CD8A anti-correlation survives the quartile extra cut: CLDN4-high (Q4) tumors have lower CD8A than CLDN4-low (Q1), n=**11 vs 11**, Δ median −0.73 log2(CPM+1), p=0.026. That is the same immune-low direction, not a new axis.

**What does not hold.** Extra continuous CLDN4 vs **CD274** and vs **HLA-A/B/C** (single genes or mean-z) are null at n=44. Do not upgrade the leftover MHC-I signature row (ρ=−0.21, p=0.16; HLA-A/B/C + B2M + TAP1 + TAP2) into an HLA-gene claim. HLA-A/B/C mean-z here is ρ=+0.005.

Do not write n=44 for the Q4 vs Q1 tests. Those tests drop Q2+Q3 and use **11 vs 11**.

---

## Honest n

| Item | n |
|---|---|
| GEO series / author count-matrix columns | **44** |
| Joined patient table (title = sample ID, 1:1) | **44** |
| Unique patients | **44** |
| CLDN4, CD8A, CD274, HLA-A, HLA-B, HLA-C non-NA | **44 / 44** (no dropout) |
| CLDN4 Q1 / Q2 / Q3 / Q4 | **11 / 11 / 11 / 11** |
| CD8A Q1 / Q2 / Q3 / Q4 | **11 / 11 / 11 / 11** |
| Q4 vs Q1 used | **11 vs 11** |
| RECIST / DCB / PFS / MPR on GEO | **0** |
| HLA-I genes in the mean-z | HLA-A, HLA-B, HLA-C (3/3) |

Quartiles are `pd.qcut(rank(method="first"), 4)` among the 44 tumors. Ties still fill four bins of 11.

---

## Extra continuous — CLDN4 vs CD274 / HLA

Spearman on all 44 tumors. log2(CPM+1).

| Endpoint | n | ρ | p | Kind |
|---|---:|---:|---:|---|
| CD8A | 44 | −0.403 | 0.0067 | known leftover |
| **CD274** | 44 | **−0.169** | **0.273** | extra |
| **HLA-A** | 44 | **−0.040** | **0.798** | extra |
| **HLA-B** | 44 | **−0.121** | **0.435** | extra |
| **HLA-C** | 44 | **+0.124** | **0.421** | extra |
| **HLA-I (HLA-A/B/C mean-z)** | 44 | **+0.005** | **0.972** | extra |
| HLA-E | 44 | −0.196 | 0.203 | supporting |
| HLA-F | 44 | −0.037 | 0.810 | supporting |
| HLA-G | 44 | −0.048 | 0.759 | supporting |
| HLA-DRA | 44 | −0.139 | 0.369 | supporting |
| B2M | 44 | −0.385 | 0.0098 | supporting; not an HLA gene |
| TACSTD2 (companion) | 44 | +0.613 | 9.7e-06 | companion only |

CD274 and the classical HLA class I genes do not track CLDN4. B2M is the only antigen-presentation gene in this extra list with a nominal p<0.05; leftover MHC-I mixed B2M/TAP into HLA-A/B/C, which is why that older signature looked weakly negative. This page keeps HLA-A/B/C separate.

Same-axis CD8 extras (not requested; not a CD274/HLA claim): CD8B ρ=−0.420, p=0.0045; CD8A+CD8B mean-z ρ=−0.457, p=0.0018. Same direction as the known CD8A row.

---

## Extra Q4 vs Q1 CD8

Two-sided Mann–Whitney U. Rank-biserial r = 2U/(n4 n1) − 1 (positive = Q4 higher).

### CLDN4 quartiles (endpoint = CD8 / CD274 / HLA)

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| **CD8A** | **11 / 11** | 3.92 | 4.65 | −0.73 | 26 | **−0.57** | **0.026** |
| CD8sig (CD8A+CD8B mean-z) | 11 / 11 | −0.28 | 0.46 | −0.74 | 23 | −0.62 | 0.015 |
| CD274 | 11 / 11 | 3.96 | 4.30 | −0.34 | 46 | −0.24 | 0.358 |
| HLA-I mean-z | 11 / 11 | −0.008 | −0.017 | +0.009 | 58 | −0.041 | 0.895 |
| HLA-A | 11 / 11 | 5.35 | 7.00 | −1.65 | 55 | −0.091 | 0.741 |
| HLA-B | 11 / 11 | 0.72 | 3.81 | −3.09 | 39 | −0.36 | 0.165 |
| HLA-C | 11 / 11 | 2.78 | 1.96 | +0.82 | 73 | +0.21 | 0.419 |

### CD8A quartiles (endpoint = CLDN4)

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| CLDN4 | **11 / 11** | 6.90 | 7.69 | −0.79 | 33 | −0.45 | 0.076 |

The reverse cut (CD8-high vs CD8-low) is the same sign as the known continuous row and does not reach p<0.05 at 11 vs 11. Do not treat that as a second independent discovery.

---

## Methods (this slice)

- Source: GEO `GSE182328_Gene_counts_matrix.txt.gz` (author HUGO × 44 counts) and series matrix for GSM / title / Akkermansia characteristics.
- Unit: patient = GEO title = count-matrix column. 44/44 join.
- Expression: `log2(CPM+1)` from library-size CPM. Same transform as the leftover GSE182328 CD8A row.
- Predictor: CLDN4. TACSTD2 is a companion only.
- HLA-I score: unweighted mean of per-gene z-scores for HLA-A, HLA-B, HLA-C (3/3 present). Not the leftover MHC-I list (that one added B2M/TAP1/TAP2).
- Quartiles: `pd.qcut(rank(method="first"), 4)` on the 44 non-NA values.
- Continuous test: two-sided Spearman.
- Quartile test: two-sided Mann–Whitney U, rank-biserial as above.
- No ICI-response test. GEO does not deposit RECIST/DCB/PFS/MPR.

---

## What this does not claim

- It does not re-audit leftover CLDN4 vs CD8A / IFN / ImmuneScore / Akkermansia.
- It does not claim CLDN4 vs PD-L1 protein or HLA protein. These are RNA rows.
- It does not claim CLDN4 vs HLA class I after the leftover MHC-I signature already failed (p=0.16). The gene-level extra cut is also null.
- It does not treat n=44 as the Q4 vs Q1 n.
- 11 vs 11 is modest. The CD8 quartile hit is a coarsened version of the known continuous test, not a new cohort.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `tables/extra_cuts.tsv` — every extra (and the known CD8A reprint)
- `tables/n_table.tsv` — honest n
- `tables/samples.tsv` — per-patient log2CPM + quartiles
- `tables/summary.json`
- `figures/cldn4_vs_cd274.png`, `cldn4_vs_hla_i.png`, `cldn4_q4q1_cd8a.png`, `cd8a_q4q1_cldn4.png`
