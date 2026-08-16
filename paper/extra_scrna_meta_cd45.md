# Extra table — immune-compartment TACSTD2/CLDN4 leak vs MPR

Additive. User A3 (malignant epithelium) is taken as given.

Public **CD45+ / immune-only** lung ICI or neoadjuvant scRNA was scored for
patient-level TACSTD2 and CLDN4 **detection fraction**. This is residual /
ambient / epithelial leak in immune libraries, **not malignant RNA**.

| Gene | Series | n (MPR vs non-MPR) | Median frac | MWU p | Hedges’ g (MPR−non) |
| --- | --- | ---: | --- | ---: | --- |
| TACSTD2 | GSE243013 CD45+ atlas | 130 vs 112 | 0.00618 vs 0.01085 | 9.76×10⁻⁴ | −0.44 (−0.70, −0.19) |
| TACSTD2 | GSE229353 CD45+ beads | 2 vs 4 | 0.0130 vs 0.0245 | 0.53 | −0.55 (−2.28, 1.18) |
| TACSTD2 | RE meta | 132 vs 116 | — | 5.85×10⁻⁴ | −0.44 (−0.70, −0.19); I²=0 |
| CLDN4 | GSE243013 CD45+ atlas | 130 vs 112 | 0.00204 vs 0.00523 | 4.79×10⁻⁵ | −0.52 (−0.78, −0.26) |
| CLDN4 | GSE229353 CD45+ beads | 2 vs 4 | 0.00956 vs 0.0210 | 0.80 | −0.56 (−2.28, 1.17) |
| CLDN4 | RE meta | 132 vs 116 | — | 6.2×10⁻⁵ | −0.52 (−0.77, −0.26); I²=0 |

GSE154826 has no ICI labels on the scRNA patients. GSE229353 P03 MTX is
truncated on GEO (one MPR library dropped). GSE243013 RECIST CR/PR vs SD/PD
is not significant (TACSTD2 p=0.54; CLDN4 p=0.11).

Write-up: `results/extra/scrna_meta_cd45/`.
