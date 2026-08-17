# Finding — TCGA-LUAD CLDN4 Q4 vs Q1 vs ImmuneScore / CD8 / CD274

**Single-gene extra.** Predictor is **CLDN4 RNA only**. TJ7 vs CD8 is already known (PR #69 / #90: TJ-claim 7-gene z-mean vs CD8 ρ = −0.31, n = 515) and is **not recomputed**. No TACSTD2 gate.

Primary matrix: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2`, log2(RSEM norm_count+1), primary tumors (`-01`). ImmuneScore = MD Anderson ESTIMATE `Immune_score` on the matching RNAseqV2 freeze. **Honest n = 515** (HiSeqV2 ∩ ESTIMATE). Do not write n = 502: that drop is ABSOLUTE-only and is not required here.

TCGA-LUAD is treatment-naive surgical RNA. **No ICI response labels.** Do not read ImmuneScore / CD8 / CD274 as immunotherapy outcomes.

## Verdict

| Endpoint | Q4 vs Q1 (129 vs 129) | Spearman (n=515) |
|---|---|---|
| ImmuneScore | Δ median **−244**; r=−0.12; p=**0.084** | ρ=−0.065; p=0.14 |
| CD8A | Δ median **−0.185**; r=−0.20; p=**5.9×10⁻³** | ρ=−0.118; p=**0.0072** |
| CD8 score (CD8A/CD8B mean) | Δ median **−0.364**; r=−0.21; p=**3.4×10⁻³** | ρ=−0.115; p=**0.0088** |
| CD274 (PD-L1) | Δ median **+0.044**; r=+0.03; p=**0.66** | ρ=+0.045; p=0.31 |

**What holds:** CLDN4-high (Q4) is **CD8-low** vs Q1. The effect is small (|r| ≈ 0.20; |ρ| ≈ 0.12). Same sign as the known TJ7–CD8 inverse, weaker than the 7-gene module.

**What does not hold:** ImmuneScore is the same sign and **not significant** at 129 vs 129 (p=0.084) or on the continuous test (p=0.14). CD274 is **null**.

Do not upgrade this to “CLDN4-high LUAD is immune-cold.” CD8 is a small negative; total ImmuneScore and PD-L1 RNA are not.

## Honest n

| | n |
|---|---:|
| HiSeqV2 primary tumors (`-01`) | **515** |
| MDACC ESTIMATE primaries (`-01`) | **515** |
| HiSeqV2 ∩ ESTIMATE (analysis set) | **515** |
| CLDN4 / CD8A / CD8B / CD274 / ImmuneScore non-NA | **515 / 515 / 515 / 515 / 515** |
| Q1 / Q2 / Q3 / Q4 | **129 / 129 / 128 / 129** |
| Q4 vs Q1 used (all four endpoints) | **129 vs 129** |
| ABSOLUTE-restricted set used elsewhere | 502 (13 dropped; **not used**) |
| TJ7 genes in this test | **0** |
| TACSTD2 filter | **none** |

Quartiles are assigned on all 515 CLDN4 values. Do not write n=515 for the MWU tests — those are **129 vs 129**. Do not write n=502 for this slice.

## Q4 vs Q1 (primary)

CLDN4 cuts among the 515: Q1 ≤ 12.52, Q4 ≥ 13.66 (HiSeqV2 log2(norm_count+1)). Two-sided Mann–Whitney U. Rank-biserial r = 2U/(n4 n1) − 1 (positive = Q4 higher).

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| ImmuneScore | 129 / 129 | 894 | 1138 | −244 | 7283 | −0.12 | 0.084 |
| CD8A | 129 / 129 | 7.91 | 8.09 | −0.185 | 6671 | −0.20 | **5.9×10⁻³** |
| CD8 score | 129 / 129 | 6.82 | 7.18 | −0.364 | 6567 | −0.21 | **3.4×10⁻³** |
| CD274 | 129 / 129 | 6.05 | 6.00 | +0.044 | 8588 | +0.03 | 0.66 |

## Spearman (supporting, n=515)

| Endpoint | n | ρ | p |
|---|---:|---:|---:|
| ImmuneScore | 515 | −0.065 | 0.14 |
| CD8A | 515 | −0.118 | **0.0072** |
| CD8 score | 515 | −0.115 | **0.0088** |
| CD274 | 515 | +0.045 | 0.31 |

Q4 vs Q1 and Spearman agree: CD8 yes (small), ImmuneScore no, CD274 no.

## Methods (this slice)

- Expression: UCSC Xena HiSeqV2 `https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz`. Primary `-01` only; one column per barcode already.
- ImmuneScore: MD Anderson ESTIMATE RNAseqV2 `https://bioinformatics.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt`.
- Predictor: CLDN4 only. TJ7 (`CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN`) is not scored. TACSTD2 is not a gate or covariate.
- Quartiles: `pd.qcut(rank(method="first"), 4)` so ties still fill four bins.
- CD8A = HiSeqV2 `CD8A`. CD8 score = mean(`CD8A`, `CD8B`), the B3-style definition, reported as a sensitivity column.
- CD274 = HiSeqV2 `CD274` (PD-L1 RNA).
- Primary test: two-sided Mann–Whitney U, Q4 vs Q1. Supporting: Spearman on all 515.

## Already known (not this extra)

| Claim | Result | Source |
|---|---|---|
| TJ7 vs CD8 in TCGA-LUAD | ρ = −0.31, p = 3.3×10⁻¹³, n=515 | PR #90 |
| TJ15 vs CD8, ESTIMATE-partial | ρ = −0.26, p = 1.7×10⁻⁹, n=511 | PR #69 / #90 |

This folder asks whether **CLDN4 alone**, as a Q4 vs Q1 cut, also moves ImmuneScore, CD8, and CD274. CD8 yes, small. The other two no.

## Sensitivity (not primary)

A prior recut on the **STAR TPM ∩ ABSOLUTE n=502** table (PR #311) reported CLDN4 Q4 vs Q1 ImmuneScore p=0.007 (126 vs 126) and CD8A p=0.036. That table drops 13 HiSeqV2+ESTIMATE primaries that lack ABSOLUTE and uses a different quantification. STAR CD274 on that same table is still not a clear hit (p=0.059, opposite-sign trend). **Primary numbers above use HiSeqV2 n=515**, which is the freeze ESTIMATE was computed on. Do not pool the two p-values.

ImmuneScore is a term in ESTIMATE_score. This test does not residualize ImmuneScore on ESTIMATE purity (circular). No purity adjustment is claimed.

## What this does not claim

- It does not re-test TJ7 vs CD8.
- It does not use a TACSTD2-high / TACSTD2-low gate.
- It does not test ICI response, PFS, or OS.
- It does not treat n=502 as the analysis n.
- It does not treat n=515 as the Q4 vs Q1 n (that is 129 vs 129).
- It does not call CLDN4-high LUAD immune-excluded or PD-L1-high.

## Outputs

- `results/n_table.tsv` — honest n
- `results/q4_vs_q1.tsv` — MWU Q4 vs Q1
- `results/spearman.tsv` — continuous Spearman
- `results/sample_scores.tsv` — per-tumor CLDN4, quartile, endpoints
- `results/summary.json`
- `figures/fig_q4_vs_q1.png`

```bash
# optional refresh of the 29 MB HiSeqV2 matrix
curl -fsSL -o data/tcga_luad_cldn4_q4/TCGA.LUAD.HiSeqV2.gz \
  https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz
curl -fsSL -o data/tcga_luad_cldn4_q4/MDACC_estimate_LUAD_RNAseqV2.txt \
  https://bioinformatics.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt

pip install -r methods/tcga_luad_cldn4_q4/requirements.txt
python3 methods/tcga_luad_cldn4_q4/analyze.py
```

The slim gene matrix `harvested/hiseqv2_cldn4_cd8_cd274.tsv` is enough to reproduce without re-downloading HiSeqV2.
