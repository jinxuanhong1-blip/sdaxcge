# Milo/scVI concordant-4: TACSTD2 T/NK neighbourhoods conditional on CLDN4

ADDITIVE. Same 65 locked units (GSE123902, GSE131907, GSE205335, GSE189357). The locked malignant CLDN4 % positive vs T/NK Spearman is unchanged. This folder asks whether T/NK neighbourhoods that are down in TACSTD2-high units stay down after CLDN4 is in the same model. The graph is a transcriptional kNN on the scVI latent, not a tissue distance and not spatial exclusion.

## Design

scVI 1.3.3, n_latent=20, n_layers=2, negative binomial, patient batch, seed 1, 200 epochs. Cells in the graph: 28643. Each QC-pass cell was kept with probability 0.10 so embedded T/NK counts stay proportional to the unit. The 220/140/40 cap used for the earlier scVI integration was not used: that cap binds for almost every unit and would erase T/NK abundance. Embedded T/NK count vs full-unit T/NK count Spearman ρ=0.993 (n=65 units present in the graph).

Milo k=30 on all 20 latent dimensions, prop=0.1, seed=1. Tested neighbourhoods (cells from ≥5 units): 2158. edgeR 4.0.16, TMM, quasi-likelihood F, k-distance SpatialFDR. Design rank of the joint model: 6/6 on 65 samples.

CLDN4 high-rich is the locked Q3+Q4 (31 vs 34). TACSTD2 high-rich is the within-cohort quartile Q3+Q4 of malignant TACSTD2 % positive (31 vs 34). Discordant units: TACSTD2-only 9, CLDN4-only 9, both high 22, neither 25. VIF of TACSTD2-high on CLDN4-high plus dataset = 1.25.

## Patient-level companion (not Milo)

- Malignant CLDN4 %pos vs T/NK: ρ=-0.531 (p=1.65e-05, I²=0.0%, -0.697 to -0.312, N=65).
- Malignant TACSTD2 %pos vs T/NK: ρ=-0.112 (p=0.46, I²=15.4%, -0.388 to 0.183, N=65).
- TACSTD2 %pos vs CLDN4 %pos, within-cohort DL: ρ=0.535 (p=0.00583, I²=56.3%, 0.171 to 0.770, N=65).
- Partial Spearman, TACSTD2 vs T/NK given CLDN4, within-cohort DL: ρ=0.195 (p=0.15, I²=0.0%, -0.072 to 0.436, N=65).
- Partial Spearman, CLDN4 vs T/NK given TACSTD2, within-cohort DL: ρ=-0.543 (p=9.5e-06, I²=0.0%, -0.705 to -0.327, N=65).

## Neighbourhood DA

SpatialFDR < 0.05. log2FC is the high arm versus the low arm.

| model | class | tested | down | up | median log2FC |
|---|---|---:|---:|---:|---:|
| TACSTD2 high | T/NK | 977 | 1 | 2 | 0.261 |
| TACSTD2 high | malignant | 494 | 4 | 22 | 0.532 |
| TACSTD2 high | mixed | 29 | 0 | 1 | -0.125 |
| TACSTD2 high | other | 658 | 0 | 0 | 0.035 |
| CLDN4 high | T/NK | 977 | 56 | 1 | -0.470 |
| CLDN4 high | malignant | 494 | 2 | 269 | 1.568 |
| CLDN4 high | mixed | 29 | 5 | 3 | -0.529 |
| CLDN4 high | other | 658 | 67 | 7 | -0.500 |
| TACSTD2 | CLDN4 | T/NK | 977 | 0 | 20 | 0.632 |
| TACSTD2 | CLDN4 | malignant | 494 | 6 | 11 | -0.004 |
| TACSTD2 | CLDN4 | mixed | 29 | 0 | 1 | 0.061 |
| TACSTD2 | CLDN4 | other | 658 | 0 | 6 | 0.352 |
| CLDN4 | TACSTD2 | T/NK | 977 | 97 | 0 | -0.800 |
| CLDN4 | TACSTD2 | malignant | 494 | 2 | 238 | 1.560 |
| CLDN4 | TACSTD2 | mixed | 29 | 3 | 2 | -0.355 |
| CLDN4 | TACSTD2 | other | 658 | 66 | 3 | -0.690 |

## Does the TACSTD2 T/NK-down set shrink conditional on CLDN4?

TACSTD2-only T/NK-down neighbourhoods: n=1, median log2FC -1.415. The same coefficient with CLDN4 in the model: median log2FC -0.935 (median paired change 0.480). 100.0% have a smaller absolute coefficient. 0 remain down at SpatialFDR < 0.05. Fewest units in that set: 5.

All tested T/NK neighbourhoods, not only the hits: n=977, median TACSTD2 log2FC 0.261 alone and 0.632 conditional on CLDN4 (median paired change 0.420, Wilcoxon p=1.37e-141). That p is the paired shift across every tested T/NK neighbourhood. The median change is the quantity.

Symmetric set, CLDN4-only T/NK-down neighbourhoods: n=56, median log2FC -1.585 alone and -1.716 conditional on TACSTD2 (36 of those 56 still down; Wilcoxon p=0.458). The CLDN4 coefficient in the joint model calls 97 T/NK neighbourhoods down, against 56 in the CLDN4-only model.

Continuous sensitivity (log2FC per 10 percentage points of TACSTD2 %pos). T/NK-down set n=5, median -0.325 alone and -0.126 with CLDN4 %pos in the model (1 still down).

TACSTD2-high calls 1 T/NK neighbourhood down at SpatialFDR < 0.05. A shrinkage fraction is not a result when the unadjusted set has fewer than 5 neighbourhoods. The CLDN4 T/NK-down median log2FC does not move toward 0 when TACSTD2 is added.

## What this does not say

- A neighbourhood on the scVI latent is not a radius around a tumour cell.
- The 10% sample is not the full tissue. Library sizes are about a tenth of the unit.
- GSE131907 is sampled from the annotation and then QC-filtered. The other three cohorts are QC-filtered and then sampled. Both use the same 10% hash.
- No fifth cohort. No private 8KL matrix. No Visium co-localisation written as exclusion.

