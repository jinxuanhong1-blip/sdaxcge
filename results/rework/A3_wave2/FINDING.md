# A3 wave-2 — GSE207422 malignant TACSTD2 (CopyKAT-like CNV)

**Claim is not recovered.** Recomputing aneuploid-like epithelium on the public UMI does not turn A3 into a significant NMPR>MPR result or a per-patient vs T/NK ρ in [−0.50, −0.40] on the pre-specified tests.

Author CopyKAT barcode IDs are still not public. Nothing here is that object.

## What was run

Public GEO UMI only: **92,330** cells × 24,292 genes (user PPT ~90,652; deposited matrix is 92,330). Marker lineages: epithelial **10,669**, stromal **1,303**, T **34,719**, NK **3,702**.

CopyKAT-like score (not the R package): 11,572 positioned genes, win=25, stromal reference n=1,303 (no immune fallback needed), score = sum |smoothed|. Primary malignant = epithelial AND score > stromal p95 (**6,080** cells). Sensitivities: p90 (7,065), mean+2 SD (4,260), 2-means (7,975), p95 minus normal-lung program (5,033), DRMref malignant (2,051), all epithelial, CNV-low epithelial.

Unit of every test is the **12 post-treatment patients** (MPR n=4 including pCR P06; NMPR n=8). Cell-level p-values are not reported.

## Primary (pre-specified)

| Test | n | Result | Match |
|---|---|---|---|
| CNV-p95 malignant mean log1p(CP10k) TACSTD2, NMPR vs MPR | 8 vs 4 | NMPR 1.86 vs MPR 1.45; Δ=+0.41; U=25; **p=0.15** | **PARTIAL_direction_NS** |
| same vs lineage T/NK fraction | 12 | **ρ=+0.077, p=0.81** | **MISMATCH** (wrong sign) |
| same vs CD3E/CD8A/NKG7 UMI≥1 fraction | 12 | **ρ=+0.12, p=0.71** | **MISMATCH** |
| same vs CD3E+ or (NKG7+ CD3E−) | 12 | ρ=+0.13, p=0.68 | **MISMATCH** |
| OLS residual on epithelial fraction vs lineage T/NK | 12 | ρ=+0.25, p=0.43 | **MISMATCH** |

Closest prior numbers: DRMref malignant NMPR vs MPR p=0.15; vs DRMref T/NK ρ=−0.49 p=0.21 (Student-t). This run’s CNV-p95 NMPR vs MPR is the same p and U as that DRMref mean test. The T/NK half is not: marker/UMI T/NK correlations are **positive**.

## %positive TACSTD2 (UMI ≥1 / ≥2 / ≥3) vs MPR

CNV-p95 malignant, n=8 vs 4, exact Wilcoxon:

| Metric | NMPR | MPR | Δ | p | Match |
|---|---|---|---|---|---|
| % UMI≥1 | 0.975 | 0.952 | +0.023 | 0.36 | PARTIAL_direction_NS |
| % UMI≥2 | 0.923 | 0.870 | +0.053 | 0.72 | PARTIAL_direction_NS |
| % UMI≥3 | 0.888 | 0.774 | +0.114 | 0.30 | PARTIAL_direction_NS |

Positivity is high in both groups. Cutoffs do not create a significant NMPR>MPR split.

## Other malignant definitions vs MPR (mean log1p CP10k)

All n=8 vs 4 except p95_nnl (several MPR samples have 0 cells after the normal-lung filter).

| Malignant def | NMPR | MPR | p | Match |
|---|---|---|---|---|
| CNV p90 | 1.85 | 1.40 | 0.15 | PARTIAL_direction_NS |
| CNV mean+2 SD | 1.83 | 1.44 | 0.28 | PARTIAL_direction_NS |
| CNV 2-means | 1.82 | 1.37 | 0.073 | PARTIAL_direction_NS |
| CNV p95, not normal-lung | 2.07 | 0.91 | 0.25 | PARTIAL_direction_NS (n_MPR with data = 1–4) |
| DRMref malignant | 1.57 | 1.12 | 0.15 | PARTIAL_direction_NS (wave-1 repeat) |
| all epithelial | 1.60 | 1.26 | 0.11 | PARTIAL_direction_NS |
| CNV-low epithelial | 1.42 | 1.14 | 0.37 | PARTIAL_direction_NS |

No definition reaches p<0.05. Direction is always NMPR higher. That is not a claim recovery.

## T/NK definitions and residuals

Pre-specified MATCH window: ρ ∈ [−0.50, −0.40]. NEAR: [−0.55, −0.35].

| TACSTD2 | T/NK def | Family | n | ρ | p | Match |
|---|---|---|---|---|---|---|
| CNV-p95 mean | lineage T/NK | raw | 12 | +0.077 | 0.81 | MISMATCH |
| CNV-p95 mean | CD3E/CD8A/NKG7 | raw | 12 | +0.12 | 0.71 | MISMATCH |
| CNV-p95 mean | DRMref T/NK | raw | 12 | −0.35 | 0.27 | MISMATCH |
| CNV-p95 %≥1 | DRMref T/NK | raw | 12 | −0.39 | 0.21 | NEAR |
| CNV-p95 mean | DRMref T/NK | partial epi frac | 12 | **−0.41** | 0.19 | MATCH window, **p NS** |
| CNV-p95 %≥1 | DRMref T/NK | partial epi frac | 12 | **−0.45** | 0.14 | MATCH window, **p NS** |
| DRMref mean | DRMref T/NK | raw | 12 | **−0.49** | 0.11 | MATCH window, **p NS** |
| DRMref mean | DRMref T/NK | residual epi frac | 12 | −0.34 | 0.29 | MISMATCH |
| all epithelial mean | lineage T/NK | raw | 12 | −0.098 | 0.76 | MISMATCH |

Author cell types do not exist on GEO. The only ρ values that land in the user window use **DRMref T/NK** (third-party labels, 30,877 cells) and are still non-significant. Wave-1 reported the same DRMref ρ=−0.4895 with a Student-t p=0.21; `scipy.stats.spearmanr` here gives p=0.11. Same coefficient, different p convention.

OLS residual of CNV-p95 TACSTD2 on epithelial fraction vs lineage/UMI T/NK is **positive** (ρ +0.25 to +0.31). Residualizing on epithelial fraction does not produce the claimed anti-correlation under public T/NK definitions.

## Honest limits

1. This is a window-smoothed expression-CNV score. It is not Hu et al. CopyKAT v0.1.0 and not inferCNV HMM.
2. n=12 (4 vs 8) is underpowered for ρ≈−0.45. A true ρ=−0.45 has two-sided Spearman p≈0.14 at n=12.
3. MPR residual tumor is 0–0.09. CNV-high cells in P06/P11/P14 can be noisy or residual epithelium; they were not dropped from the primary (all four MPR samples have ≥19 p95 cells).
4. User ~90,652 cells is not the public matrix (92,330). No public subset of 90,652 was found.
5. Do not cite the DRMref ρ=−0.49 as a new CNV recovery. It is the wave-1 third-party label, still NS.

## Files

- `stats.tsv` — every test (n, ρ/p, match)
- `verdict.tsv` — compact primary + requested sensitivities
- `per_sample.tsv` — per-patient TACSTD2 and T/NK fractions
- `cell_calls.tsv.gz` — 92,330 barcodes with lineage + CNV + flags
- `cnv_info.json` / `summary.json`
- `fig_cnv_score.png` / `fig_tacstd2_nmpr_tnk.png`
- `METHODS.md`
