# Finding — CPTAC CLDN4 protein Q4 vs Q1 vs ImmuneScore / GEP18 / CD8A

**CLDN4 protein only.** Public CPTAC TMT freeze v1.2, **LUAD** and **LSCC** separately. **No TACSTD2 gate** at any step. Quartiles are assigned among tumors with quantified CLDN4 protein. Endpoints are RNA-layer: ESTIMATE ImmuneScore, Ayers GEP18, CD8A. Honest pairwise-complete n.

These are treatment-naive surgical cohorts (Gillette *Cell* 2020; Satpathy *Cell* 2021). **No ICI response labels.** Do not read OS/PFS or these scores as immunotherapy outcomes.

## Verdict

| Endpoint | LUAD Q4 vs Q1 (20 vs 20) | LSCC Q4 vs Q1 (20 vs 20) |
|---|---|---|
| ImmuneScore | Δ median **−606**; r=−0.16; p=**0.39** | Δ median **−2189**; r=−0.62; p=**8.4×10⁻⁴** |
| GEP18 | Δ median **−0.290**; r=−0.29; p=**0.13** | Δ median **−0.785**; r=−0.67; p=**3.0×10⁻⁴** |
| CD8A RNA | Δ median **−0.514**; r=−0.21; p=**0.27** | Δ median **−1.15**; r=−0.58; p=**1.8×10⁻³** |

**What holds:** In **LSCC**, CLDN4-protein Q4 is lower than Q1 on all three endpoints (n=20 vs 20). Continuous Spearman on the full CLDN4-complete set (n=78) is the same direction and stronger-powered (ImmuneScore ρ=−0.432, p=7.9×10⁻⁵; GEP18 ρ=−0.461; CD8A ρ=−0.437).

**What does not hold:** In **LUAD**, Q4 vs Q1 is the same sign but **not significant** at n=20 vs 20 (all p≥0.13). Continuous Spearman n=79 is also mostly null (ImmuneScore ρ=−0.168, p=0.14; CD8A ρ=−0.175, p=0.12). GEP18 is the only LUAD nominal hit (ρ=−0.237, p=0.036) and it does **not** survive the requested quartile contrast.

Do not pool LUAD+LSCC TMT values. Histology is not interchangeable here.

## Honest n

| | LUAD | LSCC |
|---|---:|---:|
| Tumors with protein ∩ RNA ∩ phenotype | 110 | 108 |
| CLDN4 protein quantified | **79** | **78** |
| CLDN4 protein NA (TMT dropout) | **31 (28%)** | **30 (28%)** |
| Q1 / Q2 / Q3 / Q4 | 20 / 20 / 19 / 20 | 20 / 19 / 19 / 20 |
| Q4 vs Q1 used (all three endpoints) | **20 vs 20** | **20 vs 20** |
| ImmuneScore / GEP18 / CD8A non-NA among CLDN4-complete | 79 / 79 / 79 | 78 / 78 / 78 |
| GEP18 genes used | 18 / 18 | 18 / 18 |
| TACSTD2 filter | **none** | **none** |

CLDN4 protein is `ENSG00000189143.9`. Missingness is real TMT dropout, not a join error. Quartiles are **not** computed on the 110 / 108 matrix n. Do not write n=110 or n=108 for the Q4 vs Q1 tests.

## Methods (this slice)

- Source: open S3 `cptac-pancancer-data / data_freeze_v1.2_reorganized/{LUAD,LSCC}/`. Filenames from the LinkedOmics CPTAC-pancan index; HEAD HTTP 200.
- Predictor: CLDN4 protein only. TACSTD2 is not a gate, not a covariate, and not a co-filter.
- Quartiles: `pd.qcut(rank(method="first"), 4)` on non-NA CLDN4 protein so ties still fill four bins.
- ImmuneScore = freeze `ESTIMATE_ImmuneScore` (RNA-derived).
- GEP18 = mean of per-gene z-scores on tumor RNA (Ayers 18-gene T-cell inflamed; z on the full tumor RNA matrix).
- CD8A = tumor RNA `ENSG00000153563.15` (log2 RSEM coding UQ 1500).
- Primary test: two-sided Mann–Whitney U, Q4 vs Q1. Rank-biserial `r = 2U/(n4 n1) − 1` (positive = Q4 higher).
- Supporting: Spearman on all CLDN4-protein-complete tumors (not Q1/Q4 only).
- Cohorts kept separate. No meta-analysis.

## LUAD (Gillette 2020; freeze v1.2)

110 tumors. CLDN4 protein n=**79** (31 NA). Q4 vs Q1 = **20 vs 20**. Protein cuts among the 79: Q1≤22.76, Q4≥23.29 (log2 reference-intensity).

### Q4 vs Q1 (primary)

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| ImmuneScore | 20 / 20 | 8179 | 8786 | −606 | 168 | −0.16 | 0.39 |
| GEP18 | 20 / 20 | −0.021 | 0.269 | −0.290 | 143 | −0.29 | 0.13 |
| CD8A RNA | 20 / 20 | 8.11 | 8.63 | −0.514 | 159 | −0.21 | 0.27 |

### Spearman (supporting, n=79)

| Endpoint | n | ρ | p |
|---|---:|---:|---:|
| ImmuneScore | 79 | −0.168 | 0.14 |
| GEP18 | 79 | −0.237 | **0.036** |
| CD8A RNA | 79 | −0.175 | 0.12 |

Q4 vs Q1 throws away Q2+Q3 (39 tumors). That is why GEP18 can be nominally negative on the continuous test and still fail the quartile contrast. Do not upgrade the LUAD GEP18 Spearman to a Q4-vs-Q1 claim.

## LSCC (Satpathy 2021; freeze folder LSCC / TCGA synonym LUSC)

108 tumors. CLDN4 protein n=**78** (30 NA). Q4 vs Q1 = **20 vs 20**. Protein cuts among the 78: Q1≤22.99, Q4≥23.78.

### Q4 vs Q1 (primary)

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| ImmuneScore | 20 / 20 | 7004 | 9193 | −2189 | 76 | −0.62 | **8.4×10⁻⁴** |
| GEP18 | 20 / 20 | −0.489 | 0.296 | −0.785 | 66 | −0.67 | **3.0×10⁻⁴** |
| CD8A RNA | 20 / 20 | 7.60 | 8.75 | −1.15 | 84 | −0.58 | **1.8×10⁻³** |

### Spearman (supporting, n=78)

| Endpoint | n | ρ | p |
|---|---:|---:|---:|
| ImmuneScore | 78 | −0.432 | **7.9×10⁻⁵** |
| GEP18 | 78 | −0.461 | **2.1×10⁻⁵** |
| CD8A RNA | 78 | −0.437 | **6.4×10⁻⁵** |

LSCC is the cohort where high CLDN4 protein lines up with lower ImmuneScore, lower GEP18, and lower CD8A. That is a protein-vs-RNA-score association in untreated resected squamous tumors, not an ICI-resistance result.

## What this does not claim

- It does not use a TACSTD2-high, TACSTD2-low, or TACSTD2-complete gate.
- It does not test CLDN4 RNA as the predictor.
- It does not test ICI response, PFS, or OS.
- It does not pool LUAD and LSCC.
- It does not treat n=110 / n=108 as the CLDN4 protein n.
- ImmuneScore, GEP18, and CD8A are RNA-derived; this is cross-layer, not CLDN4 protein vs CD8 protein.
- Q4 vs Q1 n=20 vs 20 is modest. The LUAD null is underpowered relative to a median split, not proof of a zero effect.

## Outputs

- `results/n_table.tsv` — honest n
- `results/q4_vs_q1.tsv` — MWU Q4 vs Q1
- `results/spearman.tsv` — continuous Spearman
- `results/sample_scores.tsv` — per-tumor CLDN4 protein, quartile, endpoints
- `results/summary.json`
- `results/fig_q4_vs_q1.png`

```bash
python3 methods/cptac_cldn4_q4/download.py --outdir data/cptac_cldn4_q4
python3 methods/cptac_cldn4_q4/analyze.py --data data/cptac_cldn4_q4 --outdir methods/cptac_cldn4_q4
```
