# Concordant-4 ELF3–TACSTD2–CLDN4–CLDN7

ADDITIVE. Malignant cells in the locked concordant-4 only:
GSE123902 + GSE131907 + GSE205335 + GSE189357.
Not GSE148071, GSE127465, GSE207422, or GSE154826.
TACSTD2 is not a gate. The unit is the patient / donor / sample.
Do not quote cell counts as n.

DoRothEA A+B+C is a curated regulon (weighted mean on malignant
log1p CP10k). It is not lung ChIP and not a SCENIC/cisTarget run.
ELF3 → GRHL2 and GRHL2 → CLDN4 are confidence-C edges in that
network. TACSTD2 and CLDN7 are not ELF3 or GRHL2 targets there.

## Pipeline check

Malignant CLDN4 %pos versus T/NK fraction, same 65 units:
ρ = -0.531 (p = 1.65e-05, I² = 0.0%, N = 65).
Stacked within-cohort Q4 versus Q1 rank-biserial r = -0.724
(19/16, p = 0.0002879).
This matches the locked concordant-4 result and is not a new claim.

| cohort | n | ρ | p |
|---|---:|---:|---:|
| GSE123902 | 13 | -0.659 | 0.01423 |
| GSE131907 | 21 | -0.522 | 0.0152 |
| GSE205335 | 22 | -0.435 | 0.04286 |
| GSE189357 | 9 | -0.600 | 0.08762 |

## Honest n

- **n_units = 65** (13 donors + 21 samples + 22 patients + 9 patients).
- Malignant cells in those units: 72281. Not the test n.
- Within-cell Spearman uses units with at least 30 malignant cells (P4001 is out of that layer only).

## 1. Coexpression in malignant cells

Within-cell ρ is one Spearman per unit, on UMI counts. The Wilcoxon p
treats units as replicates. Pseudobulk ρ is DerSimonian–Laird across
the four cohort Spearmans of malignant mean log1p(UMI).

| pair | within-cell median ρ | units | fraction ρ>0 | Wilcoxon p | pseudobulk ρ | p | I² | 95% CI |
|---|---|---|---|---|---|---|---|---|
| ELF3–TACSTD2 | 0.687 | 63 | 1.00 | 5.17e-12 | 0.607 | 0.008971 | 71.7% | 0.174 to 0.843 |
| ELF3–CLDN4 | 0.771 | 63 | 1.00 | 5.17e-12 | 0.818 | 1.92e-06 | 64.9% | 0.590 to 0.925 |
| ELF3–CLDN7 | 0.653 | 63 | 1.00 | 5.17e-12 | 0.725 | 2.26e-11 | 0.0% | 0.571 to 0.830 |
| TACSTD2–CLDN4 | 0.746 | 63 | 1.00 | 5.17e-12 | 0.545 | 0.0004804 | 34.4% | 0.262 to 0.742 |
| TACSTD2–CLDN7 | 0.668 | 63 | 1.00 | 5.17e-12 | 0.427 | 0.002044 | 11.4% | 0.165 to 0.632 |
| CLDN4–CLDN7 | 0.695 | 64 | 1.00 | 3.53e-12 | 0.583 | 1.17e-06 | 0.0% | 0.379 to 0.734 |
| ELF3–GRHL2 | 0.311 | 60 | 1.00 | 1.63e-11 | 0.527 | 1.98e-05 | 0.0% | 0.307 to 0.694 |
| GRHL2–CLDN4 | 0.310 | 60 | 1.00 | 1.63e-11 | 0.518 | 3.00e-05 | 0.0% | 0.295 to 0.687 |
| GRHL2–TACSTD2 | 0.290 | 60 | 0.95 | 3.13e-11 | 0.681 | 1.50e-09 | 0.0% | 0.509 to 0.800 |
| GRHL2–CLDN7 | 0.296 | 60 | 0.98 | 1.71e-11 | 0.496 | 0.01694 | 60.4% | 0.097 to 0.757 |

Secondary GRHL1 / GRHL3 context:

| pair | within-cell median ρ | units | pseudobulk ρ | p |
|---|---|---|---|---|
| GRHL1–ELF3 | 0.352 | 63 | 0.408 | 0.03135 |
| GRHL1–CLDN4 | 0.349 | 64 | 0.358 | 0.008056 |
| GRHL3–ELF3 | 0.115 | 53 | 0.197 | 0.1454 |
| GRHL3–CLDN4 | 0.125 | 53 | 0.245 | 0.06811 |

## 2. Four-gene co-detection

Ratio is P(ELF3, TACSTD2, CLDN4, and CLDN7 all > 0) divided by the
product of the four detection rates. The all-unit median ratio is
3.05 (median joint detection 0.315, 65 units).

| cohort | units | median P(all four >0) | median ratio vs independence | ELF3 | TACSTD2 | CLDN4 | CLDN7 |
|---|---|---|---|---|---|---|---|
| GSE123902 | 13 | 0.206 | 4.09 | 0.58 | 0.37 | 0.43 | 0.35 |
| GSE131907 | 21 | 0.520 | 1.63 | 0.79 | 0.71 | 0.81 | 0.70 |
| GSE205335 | 22 | 0.352 | 3.32 | 0.63 | 0.58 | 0.61 | 0.50 |
| GSE189357 | 9 | 0.221 | 7.63 | 0.47 | 0.43 | 0.45 | 0.36 |
| ALL | 65 | 0.315 | 3.05 | 0.63 | 0.57 | 0.62 | 0.51 |

## 3. Regulon

| TF | gene | in DoRothEA ABC | confidence | targets in ABC |
|---|---|---|---|---|
| ELF3 | TACSTD2 | no | — | 53 |
| ELF3 | CLDN4 | no | — | 53 |
| ELF3 | CLDN7 | no | — | 53 |
| ELF3 | CDH1 | no | — | 53 |
| ELF3 | GRHL2 | yes | C | 53 |
| ELF3 | KRT8 | yes | C | 53 |
| GRHL2 | TACSTD2 | no | — | 79 |
| GRHL2 | CLDN4 | yes | C | 79 |
| GRHL2 | CLDN7 | no | — | 79 |
| GRHL2 | CDH1 | yes | C | 79 |
| GRHL2 | KRT8 | no | — | 79 |

Coherence is the within-cohort Spearman of each ABC target versus its
TF, on malignant mean log1p. A median near zero means the curated
target list does not move together in these malignant cells.

| TF | cohort | targets with a correlation | median ρ vs TF | fraction ρ>0 |
|---|---|---|---|---|
| ELF3 | GSE123902 | 52 | 0.505 | 0.90 |
| ELF3 | GSE131907 | 52 | 0.340 | 0.92 |
| ELF3 | GSE205335 | 53 | 0.362 | 0.89 |
| ELF3 | GSE189357 | 52 | 0.417 | 0.92 |
| GRHL2 | GSE123902 | 75 | 0.382 | 0.88 |
| GRHL2 | GSE131907 | 74 | 0.307 | 0.91 |
| GRHL2 | GSE205335 | 76 | 0.476 | 0.95 |
| GRHL2 | GSE189357 | 75 | 0.783 | 0.96 |

ELF3 wmean uses 52 targets.
GRHL2 wmean uses 75 targets,
74 after CLDN4 is removed.

## 4. Where the junction genes sit in the malignant correlation list

Percentile is within cohort, among genes with a finite Spearman against
the TF (100 = the top of the positive tail). The table is the median
of the four cohorts.

| TF | gene | median ρ across cohorts | median percentile |
|---|---|---|---|
| ELF3 | GRHL2 | 0.521 | 87.3 |
| ELF3 | TACSTD2 | 0.427 | 86.3 |
| ELF3 | CLDN4 | 0.801 | 99.7 |
| ELF3 | CLDN7 | 0.740 | 98.2 |
| ELF3 | CDH1 | 0.788 | 98.7 |
| ELF3 | EPCAM | 0.745 | 99.1 |
| ELF3 | KRT8 | 0.679 | 97.3 |
| ELF3 | KRT5 | -0.091 | 15.1 |
| ELF3 | PTPRC | -0.092 | 17.0 |
| GRHL2 | ELF3 | 0.521 | 93.4 |
| GRHL2 | TACSTD2 | 0.666 | 97.1 |
| GRHL2 | CLDN4 | 0.497 | 89.8 |
| GRHL2 | CLDN7 | 0.634 | 78.1 |
| GRHL2 | CDH1 | 0.693 | 97.6 |
| GRHL2 | EPCAM | 0.509 | 86.6 |
| GRHL2 | KRT8 | 0.287 | 66.6 |
| GRHL2 | KRT5 | 0.032 | 32.1 |
| GRHL2 | PTPRC | 0.276 | 67.1 |

## 5. Patient-level association with T/NK

`frac_tnk = n_tnk / n_cells` on the locked denominator.
The four-gene module is the mean of within-cohort z-scores.
The three-gene module drops CLDN4 so the T/NK number is not only the
pipeline check. Q4 versus Q1 is the stacked within-cohort quartile
contrast (rank-biserial r).

| score | kind | N | ρ | p | I² | 95% CI | Q4 vs Q1 r (nQ1/nQ4, p) |
|---|---|---|---|---|---|---|---|
| ELF3 %pos | primary | 65 | -0.459 | 0.00526 | 36.4% | -0.688 to -0.147 | -0.487 (19/16, 0.01494) |
| TACSTD2 %pos | primary | 65 | -0.112 | 0.4603 | 15.4% | -0.388 to 0.183 | -0.257 (19/16, 0.2024) |
| CLDN4 %pos (pipeline check) | check | 65 | -0.531 | 1.65e-05 | 0.0% | -0.697 to -0.312 | -0.724 (19/16, 0.0002879) |
| CLDN7 %pos | primary | 65 | -0.409 | 0.08864 | 68.4% | -0.732 to 0.066 | -0.382 (19/16, 0.05691) |
| GRHL2 %pos | primary | 65 | -0.021 | 0.8954 | 25.7% | -0.329 to 0.290 | -0.007 (19/16, 0.9868) |
| four-gene %pos module | primary | 65 | -0.449 | 0.03173 | 59.6% | -0.728 to -0.042 | -0.526 (19/16, 0.008476) |
| ELF3+TACSTD2+CLDN7 %pos (CLDN4 out) | primary | 65 | -0.412 | 0.063 | 63.1% | -0.717 to 0.024 | -0.480 (19/16, 0.01636) |
| ELF3 DoRothEA wmean | primary | 65 | -0.219 | 0.5084 | 82.0% | -0.708 to 0.411 | -0.079 (19/16, 0.7033) |
| GRHL2 DoRothEA wmean | primary | 65 | -0.268 | 0.2675 | 66.6% | -0.641 to 0.208 | -0.138 (19/16, 0.4973) |
| GRHL2 wmean, CLDN4 out | primary | 65 | -0.178 | 0.4789 | 68.1% | -0.589 to 0.307 | -0.112 (19/16, 0.5848) |
| four-gene mean-log1p module | companion | 65 | -0.337 | 0.08995 | 52.3% | -0.639 to 0.055 | -0.395 (19/16, 0.04881) |
| three-gene mean-log1p module | companion | 65 | -0.309 | 0.191 | 65.7% | -0.663 to 0.158 | -0.316 (19/16, 0.1157) |

Within malignant cells, the four genes are coexpressed: every primary pair
has a positive Spearman in essentially every unit, with median ρ about
0.65–0.77 (63–64 units; P4001 is already out). The same pairs stay
positive on the patient pseudobulk (DL ρ 0.43–0.82). GRHL2 RNA sits
with them, more strongly across patients than inside a single cell
(within-cell median ρ about 0.3). GRHL1 is weaker. GRHL3 is near the
background.

Against T/NK, malignant ELF3 %pos is negative (ρ = −0.459, p = 0.0053),
in the same direction as the locked CLDN4 check. The four-gene %pos
module is also negative (ρ = −0.449, p = 0.032). Dropping CLDN4, the
continuous pool is −0.412 (p = 0.063, CI crosses zero) and the stacked
Q4 versus Q1 contrast is r = −0.480 (p = 0.016). TACSTD2 %pos, GRHL2
%pos, and both DoRothEA wmeans have intervals that include zero.
GSE205335 changes sign for several of those non-CLDN4 scores, which is
why their I² is high. ELF3's own detection tracks low T/NK more clearly
than the average of its DoRothEA targets (52 or 53 symbols, depending on the matrix).

CLDN4 is at the top of the malignant ELF3 correlation list (median
percentile about 99.7), with EPCAM and KRT8. TACSTD2 is in the upper
tail (about the 86th percentile) and is the stronger GRHL2 neighbor
(about the 97th). KRT5 and PTPRC sit in the lower tail of the ELF3 list.

## What this does not say

- A public prior edge is not lung ELF3 or GRHL binding.
- Pseudobulk coexpression is shared malignant-program variation across
  patients. Within-cell ρ is the same-cell measurement.
- GSE205335 mixes histologies. That cohort stays in the locked n.
- No dual-high TACSTD2∩CLDN4 split, no private 8KL matrix, no Visium.

## Reproduce

```bash
bash methods/concordant4_elf3_regulon/scripts/download.sh /tmp/geo_c4
python3 methods/concordant4_elf3_regulon/scripts/analyze.py --geo /tmp/geo_c4
```
