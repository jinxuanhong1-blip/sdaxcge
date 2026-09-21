# Concordant-4: maximum |ρ| for ELF3 and a 4-gene module

Malignant cells, locked concordant-4 only (GSE123902 + GSE131907 + GSE205335 + GSE189357). The unit is the patient / donor / sample (**n = 65**). Not GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526. `frac_tnk = n_tnk / n_cells` is the locked denominator. Do not quote cell counts as n.

The association is a DerSimonian–Laird Spearman. A score is called concordant when the cohort coefficient is negative in all four cohorts. The grid is 1,030 specifications after dropping module rows that lose a member inside a cohort. p-values on the maximum are descriptive.

## Pipeline check

Malignant CLDN4 % of cells with UMI > 0 versus T/NK fraction: ρ = −0.531 (p = 1.65×10⁻⁵, I² = 0%, N = 65). The patient-level percentages match the locked table. This is the published concordant-4 result, not a new claim.

## ELF3

Across 28 marginal ELF3 scores, the largest |ρ| is the same score as before: the percent of malignant cells with ELF3 UMI ≥ 1.

| | ρ | p | I² | 95% CI |
| --- | ---: | ---: | ---: | --- |
| ELF3 %pos, UMI ≥ 1 | −0.459 | 0.00526 | 36% | −0.688 to −0.147 |

Cohorts: GSE123902 −0.736, GSE131907 −0.483, GSE205335 −0.119, GSE189357 −0.533. Stacked Q4 versus Q1 rank-biserial r = −0.487 (19/16, p = 0.015).

Raising the UMI or CP10k threshold does not increase |ρ|. The CP10k ≥ 1 call is −0.459 as well (p = 0.0023, I² = 25%). Mean log1p(UMI) is −0.345 and is positive in GSE205335. The unconstrained maximum is the same row as the concordant maximum.

## Pre-specified four-gene module

Membership fixed as ELF3 + TACSTD2 + CLDN4 + CLDN7. The previous summary, the mean of within-cohort z-scores of % positive, is ρ = −0.449 (p = 0.032, I² = 60%). GSE205335 is +0.012, so that summary is not concordant.

The concordant maximum inside this membership is the mean of within-cohort ranks of the same % positive values.

| | ρ | p | I² | 95% CI |
| --- | ---: | ---: | ---: | --- |
| Rank mean of %pos | −0.448 | 0.014 | 47% | −0.699 to −0.099 |

Cohorts: −0.700, −0.546, −0.031, −0.519. Q4 versus Q1 r = −0.530 (21/16, p = 0.0067). |ρ| does not move past the previous z-score summary. It is the only CLDN7-module specification in the grid that stays negative in every cohort. Of 66 CLDN7-module rows, it ranks 8th among the 114 concordant modules of any fourth gene.

Cohort-median “high” calls looked stronger (about −0.51) and were discarded. In GSE123902, GSE205335, and parts of GSE189357 the median CP10k of TACSTD2 or CLDN7 is zero, so “percent of cells at or above the median” is 100% in every unit. Z-scoring then drops that gene. Those rows are not four-gene scores.

## Searched fourth gene

Sixteen genes measured in all 65 units were swapped in as the fourth member. The concordant maximum is the rank mean of % positive for **ELF3 + TACSTD2 + CLDN4 + KRT19**.

| | ρ | p | I² | 95% CI |
| --- | ---: | ---: | ---: | --- |
| ELF3+TACSTD2+CLDN4+KRT19 rank mean of %pos | −0.483 | 0.0042 | 40% | −0.710 to −0.165 |

Cohorts: −0.786, −0.439, −0.186, −0.519. Q4 versus Q1 r = −0.549 (19/14, p = 0.0083). The unconstrained maximum is the same row.

KRT19 alone is ρ = −0.323 (p = 0.015, I² = 1%). The module is larger than KRT19 and larger than ELF3 (−0.459). It remains smaller than CLDN4 %pos (−0.531, I² = 0).

The next concordant rows are PC1 of %pos for KRT18 (ρ = −0.462, I² = 41%) and for CLDN3 (ρ = −0.460, p = 3.7×10⁻⁴, I² = 3%, CI −0.647 to −0.220). CLDN3 is the low-heterogeneity neighbor. It is not the |ρ| maximum.

## What this does not say

- The searched p-values are not a confirmatory test. The grid is in `results/tables/sweep_grid.tsv`.
- KRT19 is an epithelial keratin. The gain over the CLDN7 rank mean is 0.035 in |ρ|.
- No cohort was dropped. GSE205335 stays in the locked n.
- No private 8KL matrix, no Visium, no dual-high gate.

## Reproduce

```bash
bash methods/concordant4_elf3_maxrho/scripts/download.sh /tmp/geo_c4
python3 methods/concordant4_elf3_maxrho/scripts/extract_panel.py --geo /tmp/geo_c4
python3 methods/concordant4_elf3_maxrho/scripts/sweep.py --panel /tmp/geo_c4/panel
```
