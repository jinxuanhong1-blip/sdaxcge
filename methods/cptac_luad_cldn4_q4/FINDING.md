# Finding — CPTAC LUAD extra: CLDN4 protein Q4 vs ImmuneScore / CD274

**LUAD only. Additive.** Public CPTAC TMT freeze v1.2 (Gillette *Cell* 2020). Predictor is **CLDN4 protein only**. **No TACSTD2 gate.** Quartiles among tumors with quantified CLDN4 protein.

**LSCC Q4 is already reported** in `methods/cptac_cldn4_q4` (ImmuneScore / GEP18 / CD8A; Q4 vs Q1 = 20 vs 20; ImmuneScore p=8.4×10⁻⁴). This slice does **not** re-run LSCC.

Primary endpoints here: ESTIMATE ImmuneScore and **CD274 RNA** (PD-L1). Supporting: CD274 protein on the same TMT matrix. Honest pairwise-complete n.

Treatment-naive surgical LUAD. **No ICI labels.** Do not read these scores as immunotherapy outcomes.

## Verdict

| Endpoint | Layer | LUAD Q4 vs Q1 | Δ median (Q4−Q1) | r | p |
|---|---|---|---:|---:|---:|
| ImmuneScore | RNA (ESTIMATE) | **20 vs 20** | −606 | −0.16 | **0.39** |
| CD274 RNA | RNA | **20 vs 20** | −0.663 | −0.33 | **0.072** |
| CD274 protein | TMT (supporting) | **20 vs 20** | +0.052 | +0.13 | **0.49** |

**What holds:** The LUAD ImmuneScore Q4 vs Q1 number matches the already-shipped joint slice (`methods/cptac_cldn4_q4`): same 20 vs 20, same Δ=−606, same p=0.39. CD274 RNA is the same sign (Q4 lower) and closer to a threshold, but it is **not** significant at n=20 vs 20 (p=0.072). Continuous Spearman n=79 is the same story (ImmuneScore ρ=−0.168, p=0.14; CD274 RNA ρ=−0.206, p=0.069).

**What does not hold:** A LUAD CLDN4-protein Q4 “cold / PD-L1-low” claim. CD274 **protein** is null and slightly opposite (Δ=+0.052, p=0.49; Spearman ρ=+0.056, p=0.63). Do not substitute CD274 RNA for CD274 protein. Do not import the LSCC Q4 hit into LUAD.

## Honest n

| | LUAD |
|---|---:|
| Tumors with protein ∩ RNA ∩ phenotype | 110 |
| CLDN4 protein quantified | **79** |
| CLDN4 protein NA (TMT dropout) | **31 (28%)** |
| Q1 / Q2 / Q3 / Q4 | 20 / 20 / 19 / 20 |
| ImmuneScore non-NA among CLDN4-complete | **79** |
| CD274 RNA non-NA among CLDN4-complete | **79** |
| CD274 protein non-NA among CLDN4-complete | **79** |
| Q4 vs Q1 used (ImmuneScore) | **20 vs 20** |
| Q4 vs Q1 used (CD274 RNA) | **20 vs 20** |
| Q4 vs Q1 used (CD274 protein) | **20 vs 20** |
| TACSTD2 filter | **none** |
| LSCC in this slice | **no** |

CLDN4 protein is `ENSG00000189143.9`. CD274 is `ENSG00000120217.14` on both RNA and TMT. Quartiles are **not** computed on the 110-tumor matrix n. Do not write n=110 for the Q4 vs Q1 tests. CD274 protein has **no** TMT dropout in this freeze (110/110 tumors); the n=79 / 20 vs 20 constraint is CLDN4 missingness, not CD274.

## Methods (this slice)

- Source: open S3 `cptac-pancancer-data / data_freeze_v1.2_reorganized/LUAD/`. Filenames from the LinkedOmics CPTAC-pancan index; HEAD HTTP 200.
- Predictor: CLDN4 protein only. TACSTD2 is not a gate, not a covariate, and not a co-filter.
- Quartiles: `pd.qcut(rank(method="first"), 4)` on non-NA CLDN4 protein so ties still fill four bins. Same rule as `methods/cptac_cldn4_q4`.
- ImmuneScore = freeze `ESTIMATE_ImmuneScore` (RNA-derived).
- CD274 RNA = tumor RNA `ENSG00000120217.14` (log2 RSEM coding UQ 1500).
- CD274 protein = same-gene TMT abundance (log2 reference-intensity). Supporting only.
- Primary test: two-sided Mann–Whitney U, Q4 vs Q1. Rank-biserial `r = 2U/(n4 n1) − 1` (positive = Q4 higher).
- Supporting: Spearman on all CLDN4-protein-complete tumors (not Q1/Q4 only).
- LSCC not downloaded and not tested.

## LUAD (Gillette 2020; freeze v1.2)

110 tumors. CLDN4 protein n=**79** (31 NA). Q4 vs Q1 = **20 vs 20**. Protein cuts among the 79: Q1≤22.76, Q4≥23.29 (log2 reference-intensity). Same cuts as the joint LUAD+LSCC Q4 slice.

### Q4 vs Q1 (primary)

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| ImmuneScore | 20 / 20 | 8179 | 8786 | −606 | 168 | −0.16 | 0.39 |
| CD274 RNA | 20 / 20 | 8.57 | 9.24 | −0.663 | 133 | −0.33 | 0.072 |
| CD274 protein | 20 / 20 | 22.10 | 22.05 | +0.052 | 226 | +0.13 | 0.49 |

### Spearman (supporting, n=79)

| Endpoint | n | ρ | p |
|---|---:|---:|---:|
| ImmuneScore | 79 | −0.168 | 0.14 |
| CD274 RNA | 79 | −0.206 | 0.069 |
| CD274 protein | 79 | +0.056 | 0.63 |

Q4 vs Q1 throws away Q2+Q3 (39 tumors). CD274 RNA is the strongest LUAD extra signal and it still fails both the quartile contrast and a conventional p<0.05 on the continuous test. Do not upgrade p=0.072 to a Q4 claim.

CD274 RNA vs CD274 protein among the 79 CLDN4-complete tumors: Spearman ρ=0.652 (p=7.6×10⁻¹¹). The layers agree with each other; they do **not** agree on a CLDN4-protein association. That is a real discordance, not a join error.

## Already reported (LSCC; not this PR)

From `methods/cptac_cldn4_q4` (same freeze, same quartile rule, same ImmuneScore definition):

| Endpoint | LSCC Q4 vs Q1 (20 vs 20) |
|---|---|
| ImmuneScore | Δ=−2189; r=−0.62; p=**8.4×10⁻⁴** |
| GEP18 | Δ=−0.785; r=−0.67; p=**3.0×10⁻⁴** |
| CD8A RNA | Δ=−1.15; r=−0.58; p=**1.8×10⁻³** |

LSCC remains the histology where high CLDN4 protein lines up with a colder RNA immune score. This LUAD extra does not add a CD274 exception that rescues LUAD.

## What this does not claim

- It does not re-analyze LSCC.
- It does not use a TACSTD2-high, TACSTD2-low, or TACSTD2-complete gate.
- It does not test CLDN4 RNA as the predictor.
- It does not test ICI response, PFS, or OS.
- It does not treat n=110 as the CLDN4 protein n.
- ImmuneScore and CD274 RNA are RNA-derived; CD274 protein is the only same-layer endpoint and it is null.
- Q4 vs Q1 n=20 vs 20 is modest. The LUAD nulls are underpowered relative to a median split, not proof of a zero effect.

## Outputs

- `results/n_table.tsv` — honest n
- `results/q4_vs_q1.tsv` — MWU Q4 vs Q1
- `results/spearman.tsv` — continuous Spearman
- `results/sample_scores.tsv` — per-tumor CLDN4 protein, quartile, endpoints
- `results/summary.json`
- `results/fig_q4_vs_q1.png`

```bash
python3 methods/cptac_luad_cldn4_q4/download.py --outdir data/cptac_luad_cldn4_q4
python3 methods/cptac_luad_cldn4_q4/analyze.py --data data/cptac_luad_cldn4_q4 --outdir methods/cptac_luad_cldn4_q4
```
