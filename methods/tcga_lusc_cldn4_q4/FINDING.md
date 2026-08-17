# Finding — TCGA-LUSC CLDN4 Q4 vs Q1 vs ImmuneScore / CD8 / CD274

**CLDN4 RNA only.** Public UCSC Xena TCGA-LUSC HiSeqV2. **LUSC only** (no LUAD, no pooled NSCLC). **No TACSTD2 gate** at any step. Quartiles are assigned among primary tumors (`-01`) with non-NA CLDN4. Endpoints: official ESTIMATE ImmuneScore, CD8A (CD8), CD274 (PD-L1 RNA). Honest pairwise-complete n.

Treatment-naive surgical TCGA. **No ICI response labels.** Do not read these scores as immunotherapy outcomes. CD274 here is **RNA**, not IHC CPS / TPS.

## Verdict

| Endpoint | Q4 vs Q1 | Δ median | r | p | verdict |
|---|---|---:|---:|---:|---|
| ImmuneScore | **126 vs 125** | −25.5 | −0.048 | **0.51** | **NO_EVIDENCE** |
| CD8A (CD8) | **126 vs 126** | −0.256 | −0.082 | **0.26** | **NO_EVIDENCE** |
| CD274 (PD-L1) | **126 vs 126** | −0.0016 | +0.011 | **0.88** | **NO_EVIDENCE** |

**What holds:** Nothing on the requested Q4 vs Q1 contrast. Continuous Spearman on the full primary set is the same: ImmuneScore ρ=−0.004 (n=501, p=0.92); CD8A ρ=−0.017 (n=502, p=0.70); CD274 ρ=+0.024 (n=502, p=0.58).

**What does not hold:** A claim that CLDN4-high LUSC is ImmuneScore-low, CD8-low, or CD274-high/low. All three tests are null. This matches the prior multi-cohort Q4 recut for LUSC ImmuneScore / CD8 (PR 311: p=0.50 / 0.26) and adds **CD274**, which was not on that table.

Do not borrow the CPTAC-LSCC **protein** Q4 immune-low result (n=20 vs 20) onto this TCGA **RNA** slice.

## Honest n

| Item | n |
|---|---:|
| Xena HiSeqV2 primary tumors (`-01`, 15-char) | **502** |
| Official ESTIMATE RNAseqV2 primaries | **501** |
| CLDN4 non-NA | **502** (0 NA) |
| CD8A non-NA | **502** (0 NA) |
| CD274 non-NA | **502** (0 NA) |
| ImmuneScore join | **501** |
| ImmuneScore NA among CLDN4-complete | **1** (`TCGA-33-4579-01`, Q1) |
| Q1 / Q2 / Q3 / Q4 | **126 / 125 / 125 / 126** |
| ImmuneScore Q4 vs Q1 used | **126 vs 125** |
| CD8A Q4 vs Q1 used | **126 vs 126** |
| CD274 Q4 vs Q1 used | **126 vs 126** |
| TACSTD2 filter | **none** |

Do **not** write n=502 for the ImmuneScore test. ESTIMATE is missing one primary; that case is in Q1, so ImmuneScore is 126 vs 125. CD8A and CD274 are complete on all 502, so those two tests are 126 vs 126.

CLDN4 RNA cuts among the 502: Q1 ≤ 11.31, median 12.07, Q4 ≥ 12.75 (log2(RSEM+1)).

## Methods (this slice)

- Source: UCSC Xena `TCGA.LUSC.sampleMap/HiSeqV2` (log2(RSEM normalized_count + 1)); official MD Anderson ESTIMATE RNAseqV2 `Immune_score`.
- Samples: TCGA barcode sample type `01` (primary solid tumor). Duplicate 15-char barcodes averaged (none after collapse).
- Predictor: CLDN4 RNA only. TACSTD2 is not a gate, not a covariate, and not a co-filter.
- Quartiles: `pd.qcut(rank(method="first"), 4)` on non-NA CLDN4 so ties still fill four bins.
- ImmuneScore = official ESTIMATE `Immune_score` (RNA-derived ssGSEA). Not recomputed.
- CD8 = `CD8A` on the same Xena matrix (single gene, not CIBERSORT / xCell).
- CD274 = `CD274` on the same Xena matrix (PD-L1 RNA, not protein / IHC).
- Primary test: two-sided Mann–Whitney U, Q4 vs Q1. Rank-biserial `r = 2U/(n4 n1) − 1` (positive = Q4 higher). 95% CI = percentile bootstrap, 2000 resamples, seed 20260817.
- Supporting: Spearman on all CLDN4-complete primaries with that endpoint (not Q1/Q4 only).
- LUSC only. No LUAD, no meta-analysis.

## Q4 vs Q1 (primary)

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| ImmuneScore | 126 / 125 | 363.3 | 388.8 | −25.5 | 7498 | −0.048 | 0.51 |
| CD8A (CD8) | 126 / 126 | 7.473 | 7.728 | −0.256 | 7289 | −0.082 | 0.26 |
| CD274 (PD-L1) | 126 / 126 | 6.399 | 6.401 | −0.0016 | 8029 | +0.011 | 0.88 |

## Spearman (supporting)

| Endpoint | n | ρ | p |
|---|---:|---:|---:|
| ImmuneScore | **501** | −0.004 | 0.92 |
| CD8A (CD8) | **502** | −0.017 | 0.70 |
| CD274 (PD-L1) | **502** | +0.024 | 0.58 |

Q4 vs Q1 throws away Q2+Q3 (250 tumors). The continuous test is better powered and is still null. Do not upgrade a median-split or a LUAD result onto this LUSC table.

## Figures

| File | Content |
|---|---|
| `figures/fig1_q4_vs_q1.png` | Q1 vs Q4 violins — ImmuneScore, CD8A, CD274 (n labeled on axes) |
| `figures/fig2_forest_extra.png` | **Extra figure:** rank-biserial forest with 95% bootstrap CI and honest n per arm |

## What this does not claim

- It does not use a TACSTD2-high, TACSTD2-low, or TACSTD2-complete gate.
- It does not test CLDN4 protein (see CPTAC LSCC protein Q4, different layer).
- It does not test ICI response, PFS, or OS.
- It does not pool LUAD and LUSC.
- It does not treat n=502 as the ImmuneScore n.
- ImmuneScore is RNA-derived and circular with ESTIMATE purity; it is reported as requested, not as an orthogonal DNA immune fraction.
- CD274 RNA is not PD-L1 IHC.
- A null at n=126 vs 126 is not proof of a zero effect; it is no evidence on this public matrix.

## Outputs

- `results/n_table.tsv` — honest n
- `results/q4_vs_q1.tsv` — MWU Q4 vs Q1
- `results/spearman.tsv` — continuous Spearman
- `results/sample_scores.tsv` — per-tumor CLDN4, quartile, endpoints
- `results/summary.json`
- `figures/fig1_q4_vs_q1.png` (+ pdf)
- `figures/fig2_forest_extra.png` (+ pdf)

```bash
pip install -r methods/tcga_lusc_cldn4_q4/requirements.txt
python3 methods/tcga_lusc_cldn4_q4/download.py --outdir data/tcga_lusc_cldn4_q4
python3 methods/tcga_lusc_cldn4_q4/analyze.py --data data/tcga_lusc_cldn4_q4 --outdir methods/tcga_lusc_cldn4_q4
```
