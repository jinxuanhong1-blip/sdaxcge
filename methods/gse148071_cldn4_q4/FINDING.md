# Finding — GSE148071 malignant CLDN4 Q4 vs Q1 vs same-patient T/NK and CXCL13+

ADDITIVE. **CLDN4 only.** No TACSTD2 gate. Patient is the unit (42 advanced
NSCLC tumors, Wu *et al.* *Nat Commun* 2021; one sample each). Continuous
Spearman ρ is **taken as given** from the existing TISCH T/NK extract and
the PR #274 / #320 malignant TLS extract. **Q4 vs Q1 is extra**
(Mann–Whitney, rank-biserial *r*). p-values are descriptive.

Existing processed scores only. The 10x matrices were not re-downloaded.
No ICI response / RECIST / MPR labels are in these extracts — do not read
T/NK or CXCL13+ as immunotherapy outcome.

## Verdict

| Contrast | Given continuous ρ | Extra Q4 vs Q1 r | Honest n |
|---|---|---|---|
| TISCH T/NK (locked eligible) | n=25 ρ=+0.135 p=0.519 | n_Q1=7 n_Q4=6 r=+0.048 p=0.945 | 25 eligible / 42; 3 epithelial-like inside the 25 |
| TISCH T/NK (malignant only) | n=22 ρ=+0.080 p=0.725 | n_Q1=6 n_Q4=6 r=-0.056 p=0.937 | 22; drops P5/P35/P39 epithelial-like |
| TLS CXCL13+ among T (locked malignant) | n=31 ρ=-0.409 p=0.0225 | n_Q1=8 n_Q4=8 r=-0.609 p=0.0444 | 31 malignant with a T-cell CXCL13 call |
| Same-patient T/NK (shared TLS CLDN4 cut) | n=21 ρ=-0.048 p=0.836 | n_Q1=6 n_Q4=5 r=-0.133 p=0.792 | 21 patients in both extracts |
| Same-patient CXCL13+ (shared TLS CLDN4 cut) | n=21 ρ=-0.447 p=0.0422 | n_Q1=6 n_Q4=5 r=-0.667 p=0.0821 | same 21 patients as the T/NK intersection |

**What holds:** On the locked TLS malignant extract, higher malignant CLDN4
tracks **lower CXCL13+ among T** (given ρ=-0.409, p=0.0225,
n=31). The extra Q4 vs Q1 cut stays negative
(r=-0.609, p=0.0444, 8 vs 8).

**What does not hold:** Same-patient T/NK is **not** a CLDN4-low-TME story.
The locked TISCH continuous ρ is weakly **positive** (n=25, ρ=+0.135,
p=0.519). Extra Q4 vs Q1 on that set is also not a negative
hit (r=+0.048, p=0.945). Joining TLS malignant CLDN4
to TISCH T/NK on the 21 overlapping patients does not flip T/NK
to a significant negative quartile effect. Do not write n=42. Do not pool
T/NK and CXCL13+ into one immune score.

## Honest n

| item | n | note |
|---|---:|---|
| GEO patients (Wu 2021) | 42 | one sample each; not the test n |
| TISCH units in extract | 42 | all 42 samples present |
| TISCH eligible (≥20 scored epi + ≥20 T/NK) | 25 | locked given-ρ set |
| TISCH eligible but epithelial-like fallback | 3 | P35,P39,P5 |
| TISCH eligible and Malignant | 22 | honest malignant T/NK recut |
| TISCH ineligible | 17 | usually <20 T/NK |
| TLS eligible (author extract) | 41 | 41/42 |
| TLS malignant compartment | 35 | epithelial / insufficient dropped |
| TLS malignant + finite CXCL13+ among T | 31 | locked given-ρ set; 4 malignant lack a T-cell CXCL13 call |
| TLS malignant + CXCL13 mean | 35 | whole-sample CXCL13; n=35 |
| Same-patient TLS CXCL13+ ∩ TISCH eligible | 21 | shared CLDN4 Q4/Q1 |
| Same-patient ∩ TISCH Malignant only | 21 | epi-like TISCH units are not in the TLS malignant+CXCL13+ set |

Quartiles are assigned **inside** each analysis set. n_compared = n_Q1 + n_Q4
(the tails), not the GEO n=42 and not the continuous n. Same-patient shared
cut on n=21: Q1=6, Q2=5, Q3=5, Q4=5.

## Methods (this slice)

- Predictor: malignant CLDN4 only. TACSTD2 is not a gate.
- T/NK (TISCH): `n_TNK / n_cells`. T/NK labels = CD8T, CD8Tex, CD4Tconv,
  Treg, Tprolif, TMKI67, NK, Tcell, NKT, ILC. Eligible = ≥20 scored
  epithelial/malignant cells and ≥20 T/NK (locked given set).
- CXCL13+ (TLS PR #274): `frac_CXCL13pos_T` among T in the same patient;
  CLDN4 = mean log1p(CP10k) in the **malignant** compartment. Epithelial
  and insufficient rows are dropped.
- Given continuous: Spearman on the locked complete sets above. Recomputed
  here only to confirm the file still matches; the ρ is not a new audit.
- Extra cut: `pd.qcut(rank(method='average'), 4)` on pairwise-complete
  CLDN4. Two-sided Mann–Whitney U, Q4 vs Q1. Rank-biserial
  `r = 2U/(n4 n1) − 1` (positive = Q4 higher).
- Same-patient extra: inner join on `P1…P42`. One TLS malignant CLDN4
  quartile cut; both endpoints tested on those tails.
- Thin flag: n<8 or either tail <3. p-values are descriptive.

## Extra Q4 vs Q1 (primary table)

| Contrast | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | r | p |
|---|---|---:|---:|---:|---:|---:|
| TISCH T/NK, locked eligible | 6 vs 7 | 0.1509 | 0.07454 | +0.07638 | +0.048 | 0.945 |
| TISCH T/NK, malignant only | 6 vs 6 | 0.1509 | 0.2033 | -0.05241 | -0.056 | 0.937 |
| TLS CXCL13+ among T, locked malignant | 8 vs 8 | 0.1076 | 0.2738 | -0.1662 | -0.609 | 0.0444 |
| Same-patient T/NK, shared TLS CLDN4 cut | 5 vs 6 | 0.1198 | 0.2033 | -0.08351 | -0.133 | 0.792 |
| Same-patient CXCL13+, shared TLS CLDN4 cut | 5 vs 6 | 0.16 | 0.2738 | -0.1138 | -0.667 | 0.0821 |
| TLS CXCL13 mean (secondary) | 9 vs 9 | 0.02897 | 0.04274 | -0.01377 | -0.407 | 0.153 |

## Given continuous Spearman (supporting; not re-audited)

| Contrast | n | ρ | p | source |
|---|---:|---:|---:|---|
| TISCH CLDN4 mean vs T/NK | 25 | +0.135 | 0.519 | TISCH pool, eligible units |
| TISCH CLDN4 %pos vs T/NK | 25 | +0.044 | 0.835 | same n=25 |
| TISCH CLDN4 mean vs T/NK, malignant only | 22 | +0.080 | 0.725 | drop 3 epi-like |
| TLS CLDN4 mean vs CXCL13+ among T | 31 | -0.409 | 0.0225 | PR #274/#320 malignant |
| TLS CLDN4 %pos vs CXCL13+ among T | 31 | -0.222 | 0.23 | same n=31 |
| TLS CLDN4 mean vs CXCL13 mean | 35 | -0.281 | 0.102 | malignant n=35 |
| Same-patient TLS CLDN4 vs TISCH T/NK | 21 | -0.048 | 0.836 | join; not the locked TISCH ρ |
| Same-patient TLS CLDN4 vs CXCL13+ | 21 | -0.447 | 0.0422 | join subset of the n=31 |

Q4 vs Q1 throws away Q2+Q3. That is why a continuous CXCL13+ hit can weaken
on the intersection tails. Do not upgrade the n=31 Spearman to an n=42 claim,
and do not upgrade the T/NK null to a negative Q4 story.

## What was not done

- No dual-high TACSTD2×CLDN4 score.
- No new 10x download or re-annotation of GSE148071.
- No ICI response test (none in these extracts).
- TISCH and TLS annotations are different; cell counts do not match 1:1.
  Same-patient uses the patient ID join, not a cell-level merge.
- T/NK was not invented on the TLS table (no T/NK column there).

Tables: `tables/q4q1_table.tsv`, `tables/honest_n.tsv`,
`tables/same_patient_quartiles.tsv`.

Figures: `figures/q4q1_tisch_tnk.png`, `figures/q4q1_tls_cxcl13.png`,
`figures/q4q1_same_tnk.png`, `figures/q4q1_same_cxcl13.png`.

Reproduce: `python3 methods/gse148071_cldn4_q4/analyze.py`
