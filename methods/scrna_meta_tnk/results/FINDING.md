# Patient-level meta: TACSTD2 / CLDN4 vs T/NK (public lung scRNA)

GSE207422 A3 is taken as given (n=12, TACSTD2 ρ=−0.490, p=0.106).

The goal is a combinatorial search over cohort subsets and definitions
(see `combinatorial/FINDING.md`). The 8-cohort merge below is **one row**.

Eight cohorts, **N=145 patients**.

| Gene | RE ρ [95% CI] | p | I² | Stouffer z (p) |
|---|---|---:|---:|---|
| TACSTD2 | −0.109 [−0.339, +0.133] | 0.378 | 39% | −0.61 (0.544) |
| CLDN4 | −0.137 [−0.306, +0.040] | 0.129 | 0% | −1.53 (0.125) |

GSE146100 skipped (n=1). Blood-only and EGA/dbGaP raw skipped. Leftover 2023–2026: GSE325414 included; other hits lacked a public malignant+T/NK table or were blood/bulk.

Forest: `forest_TACSTD2.png`, `forest_CLDN4.png`.
