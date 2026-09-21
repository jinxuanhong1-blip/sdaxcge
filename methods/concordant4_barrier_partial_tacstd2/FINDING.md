# FINDING — TACSTD2 barrier outgoing after holding CLDN4 fixed

ADDITIVE. Does not replace the locked CLDN4 communication table.
Concordant four only: GSE123902, GSE131907, GSE205335, GSE189357.
GSE148071, GSE127465, GSE154826, GSE200563, and E-MTAB-13526 are not added.
Senders are malignant cells. Receiver is T/NK.
Barrier ligands are F11R, NECTIN2, CDH1, and LGALS9.
The patient is the unit. This is an expression ligand–receptor contrast, not a spatial test.

The shrink rule was fixed before scoring: a positive crude mean falls by at least 20%, and the paired Wilcoxon p on (crude − adjusted) is < 0.05. Co-primary scores are the expression-proportion gap (percentage points) and the CellChat Hill probability with no population-size weight.

## Headline

TACSTD2 Q4 vs Q1 barrier outgoing **shrinks after CLDN4 is held fixed, and a positive partial remains.**

On the same 63 patients, the expression-proportion family mean is **+25.19** percentage points crude and **+10.44** after CLDN4-stratum residualization (58.6% smaller). Paired Wilcoxon p on the drop = 1.16×10⁻¹¹. The residual mean is still above zero (p = 3.59×10⁻¹¹, 95% of patients positive, 4/4 cohorts positive). Sign-flip one-sided p = 1.00×10⁻⁴ (0 of 10,000 flips).

CellChat Hill moves the same way: **+0.0663** crude to **+0.0270** after the same residual (59.3% smaller; paired p = 2.84×10⁻¹¹; residual p = 2.52×10⁻¹⁰; 4/4 cohorts).

Linear residualization on CLDN4 log1p agrees: expression proportion **+12.05** (52.2% smaller), Hill **+0.0302** (54.4% smaller).

The reverse partial is the same size. CLDN4 crude expression proportion is **+25.60** (n=64; this matches the locked Q4 vs Q1 result). After TACSTD2-stratum residualization it is **+12.02** (53.0% smaller; residual p = 6.21×10⁻¹²; 4/4 cohorts). Hill goes from **+0.0651** to **+0.0280** (57.0% smaller). Each gene’s outgoing shrinks when the other is held fixed. Each remainder stays positive.

## Why a remainder is identifiable

Inside malignant cells, TACSTD2 and CLDN4 log1p have mean Spearman **0.432** (median 0.439, n=63). CLDN4 strata explain mean R² **0.223** of TACSTD2 (median 0.195). Most TACSTD2 variation sits inside a CLDN4 stratum.

In the crude TACSTD2 arms the CLDN4 log1p gap is **+0.959** and the TACSTD2 gap is **+2.022**. After 1:1 matching inside CLDN4 strata the CLDN4 gap is **+0.030** and the TACSTD2 gap is **+1.992**. Mean match rate is 0.55 (pairs / n_Q4). Five patients had fewer than 10 pairs and are out of the matched row (n=58). On the matched cells the expression-proportion gap is **+14.32** (42.7% below the crude mean of those 58 patients; p = 5.90×10⁻¹¹; 4/4 cohorts) and Hill is **+0.0410**.

A within-stratum median split of TACSTD2, averaged across CLDN4 strata, is a local contrast rather than Q4 vs Q1. Its expression-proportion mean is **+7.88** (n=60, 4/4 cohorts, p = 2.32×10⁻¹¹). Hill is **+0.0193**.

## Cohort means, Q4 vs Q1

Expression proportion, family mean of the four ligands. Order: GSE123902 / GSE131907 / GSE205335 / GSE189357.

| contrast | n | mean | cohorts | p Wilcoxon | I² |
|---|---:|---:|---|---|---:|
| TACSTD2 crude | 63 (13+20+21+9) | +25.19 | +18.90 / +21.98 / +32.18 / +25.12 | 7.23×10⁻¹² | 45% |
| TACSTD2 \| CLDN4 strata | 63 | +10.44 | +6.84 / +11.62 / +12.03 / +9.29 | 3.59×10⁻¹¹ | 35% |
| TACSTD2 \| CLDN4 linear | 63 | +12.05 | +9.77 / +12.23 / +13.55 / +11.46 | 2.05×10⁻¹¹ | 0% |
| TACSTD2 matched on CLDN4 | 58 (11+19+19+9) | +14.32 | +9.89 / +13.81 / +17.65 / +13.76 | 5.90×10⁻¹¹ | 14% |
| CLDN4 crude | 64 (13+21+21+9) | +25.60 | +21.33 / +21.77 / +32.44 / +24.79 | 4.91×10⁻¹² | 34% |
| CLDN4 \| TACSTD2 strata | 64 | +12.02 | +10.74 / +12.96 / +13.26 / +8.80 | 6.21×10⁻¹² | 22% |

CellChat Hill, same patients:

| contrast | mean | cohorts | p Wilcoxon |
|---|---:|---|---|
| TACSTD2 crude | +0.0663 | +0.0410 / +0.0483 / +0.0986 / +0.0675 | 1.55×10⁻¹¹ |
| TACSTD2 \| CLDN4 strata | +0.0270 | +0.0174 / +0.0240 / +0.0353 / +0.0283 | 2.52×10⁻¹⁰ |
| CLDN4 crude | +0.0651 | +0.0436 / +0.0495 / +0.0946 / +0.0640 | 6.21×10⁻¹² |
| CLDN4 \| TACSTD2 strata | +0.0280 | +0.0194 / +0.0282 / +0.0367 / +0.0198 | 4.28×10⁻¹¹ |

CellPhoneDB `lr_means` and LIANA log2FC are the same ligand-mean contrast on two scales. TACSTD2 crude → strata residual: lr_means +0.107 → +0.0487, log2FC +0.308 → +0.141 (both 54.4% smaller, 4/4 cohorts).

## Ligands, expression proportion, TACSTD2 Q4 vs Q1

All four ligands stay positive in all four cohorts after CLDN4-stratum residualization.

| ligand | crude | strata residual | matched |
|---|---:|---:|---:|
| F11R | +25.85 | +9.96 | +14.25 |
| NECTIN2 | +28.46 | +10.96 | +15.28 |
| CDH1 | +32.26 | +13.58 | +18.69 |
| LGALS9 | +14.20 | +7.24 | +9.05 |

## Decile sensitivity

Outer 10%, same residual. TACSTD2 expression proportion +26.14 → +11.34 (n=59, 4/4). Hill +0.0743 → +0.0328. CLDN4 expression proportion +27.76 → +13.20 (n=60). The locked decile CLDN4 crude mean is +27.8; this run is +27.76.

## Who is in the n

Inventory units: 65. P4001 has 27 malignant cells and stays out, as in the locked run. CLDN4 Q4 vs Q1 keeps n=64. EBUS_13 has no TACSTD2 counts in malignant cells, so the TACSTD2 gate is n=63. Spearman is undefined there and is computed on the other 63.

Label QC matched n_mal and n_tnk on all 65 units. The CLDN4 crude expression-proportion mean is +25.60 (n=64), the locked figure.

## What this does not say

- The percentage-point gap is not a CellChat probability. Hill is reported on its own scale. Population-size weighting is not reapplied.
- Matching uses the cells that share a CLDN4 stratum (about half of TACSTD2 Q4). The strata residual uses every Q4 and Q1 cell.
- The within-stratum median split is a local TACSTD2 contrast. Its mean is not a fraction of the Q4 vs Q1 gap.
- No dual-high sender gate. No added cohort. No spatial distance.
