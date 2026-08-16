# Open BLCA ICI cohorts with CLDN4 — OR / n / p

Median split of CLDN4 (response-blind). Endpoint = CR/PR vs SD/PD except GSE111636
(depositor binary responder/progressor). Odds ratio is *response* for CLDN4-high
vs CLDN4-low, so OR < 1 is the direction claim B5 asserts (claimed OR = 0.42).

## Per-cohort

| Cohort | Open source | n | responders | OR (95% CI) | p |
|---|---|---|---|---|---|
| IMvigor210 | Nature 2018 / IMvigor210CoreBiologies | **298** | 68 | **1.47 (0.82–2.64)** | 0.214 |
| BACI | GEO GSE176307 | **87** | 16 | **1.03 (0.30–3.53)** | 1.000 |
| Snyder 2017 | Zenodo 7058399 (PredictIO) | **21** | 7 | **0.76 (0.08–6.52)** | 1.000 |
| UC-GENOME | cBioPortal `blca_bcan_hcrn_2022` | **89** | 34 | **0.71 (0.27–1.82)** | 0.515 |
| GSE111636 | GEO GSE111636 (array; coarse endpoint) | **11** | 6 | **6.44 (0.33–490)** | 0.242 |

## Pooled

| Pool | k | n | OR (95% CI) | p | vs claimed 0.42 |
|---|---|---|---|---|---|
| **Prespecified primary** (IMvigor210 + BACI + Snyder) | 3 | **406** | **1.31 (0.82–2.10)** | 0.258 | rejected, z=4.74, p=2.1×10⁻⁶ |
| Post-hoc + UC-GENOME (RECIST-like only) | 4 | **495** | **1.14 (0.75–1.72)** | 0.540 | rejected, z=4.73, p=2.3×10⁻⁶ |

I² = 0% in both pools. No cohort-level test reaches p < 0.05. Adding the one
open cohort that points in the claimed direction (UC-GENOME, OR 0.71) moves the
pooled estimate from 1.31 toward 1, not toward 0.42.

UC-GENOME and GSE111636 were added **after** the 3-cohort primary was locked
(`analysis_plan.md`, commit `b1e1159`). They are reported because they are open
and eligible; they were not used to chase 0.42.

Machine-readable copy: `or_n_p.csv`.
