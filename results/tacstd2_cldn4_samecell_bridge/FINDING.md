# TACSTD2–CLDN4 same-cell coexpression BRIDGE (not mediation)

Public-only. Unit for concordant-4 is patient/donor/sample. Cell counts are not n.
Does **not** recompute locked CosMx CLDN4 cytotoxic ratios 0.36/0.52 or
concordant-4 CLDN4 %pos vs T/NK ρ = −0.531. Does **not** claim mediation
(PR #718: TACSTD2→CLDN4→immune matrix is null).

## Slide verdict (PPT)

| Layer | Same-cell coexpression? | Overlap vs independence | Cold niche overlap |
|---|---|---|---|
| Concordant-4 malignant | **YES** — median within-cell ρ = **0.746** (100% units ρ>0; n_units≥30 = 63) | median P(both)/[P·P] = **1.31**; median OR = **13.65** | Dual-pos % vs T/NK Spearman = **-0.254** (p=4.42e-02) — descriptive bridge, not a new lock |
| CosMx He2022 tumor | **YES** — median ρ = **0.180** (range 0.093–0.257); **8/8** ρ>0 | median detection ratio = **1.17**; 8/8 sections >1 | Both-high immune frac @10 µm **0.033** vs both-low **0.091**; HH<LL **8/8** (PR #726 quadrants) |
| DepMap Gygi protein | **YES** — Spearman **0.693** (n=45, p=1.31e-07) | continuous protein coexpression (rounds to **0.69**) | N/A (cell lines; no immune niche) |

## Framing

- **BRIDGE:** the two genes mark overlapping malignant programs / cold neighborhoods.
- **NOT mediation:** holding CLDN4 does not account for TACSTD2 short-range exclusion on CosMx (PR #726), and the observational mediation matrix on CosMx + concordant-4 is null (PR #718).

## Concordant-4 cohort table

| cohort | units ≥30 | median ρ | median ratio vs indep | median OR |
|---|---:|---:|---:|---:|
| GSE123902 | 13 | 0.678 | 1.52 | 13.20 |
| GSE131907 | 20 | 0.674 | 1.12 | 12.48 |
| GSE189357 | 9 | 0.733 | 1.64 | 17.47 |
| GSE205335 | 21 | 0.823 | 1.31 | 15.62 |
| ALL | 63 | 0.746 | 1.31 | 13.65 |

## CosMx detection co-occurrence

Tumor cells = 295877 across 8 sections / 5 patients. Detection = raw count > 0.

| section | patient | n | ρ | P(both) | P(tac)·P(cld) | ratio | OR |
|---|---|---:|---:|---:|---:|---:|---:|
| LUAD-5 R1 | Lung5 | 17837 | 0.210 | 0.521 | 0.468 | 1.112 | 2.99 |
| LUAD-5 R2 | Lung5 | 17914 | 0.175 | 0.486 | 0.438 | 1.111 | 2.61 |
| LUAD-5 R3 | Lung5 | 15893 | 0.206 | 0.430 | 0.371 | 1.158 | 2.84 |
| LUSC-6 | Lung6 | 66196 | 0.171 | 0.130 | 0.092 | 1.417 | 2.36 |
| LUAD-9 R1 | Lung9 | 38864 | 0.167 | 0.229 | 0.192 | 1.191 | 2.90 |
| LUAD-9 R2 | Lung9 | 94876 | 0.184 | 0.209 | 0.164 | 1.279 | 2.94 |
| LUAD-12 | Lung12 | 18267 | 0.257 | 0.206 | 0.129 | 1.599 | 4.41 |
| LUAD-13 | Lung13 | 26030 | 0.093 | 0.291 | 0.256 | 1.137 | 2.39 |

## DepMap protein checksum

Gygi `_LUNG` complete cases: n=45, Spearman ρ=0.693
(95% bootstrap CI 0.48–0.82), p=1.31e-07.
Rounds to the slide value 0.69. Slide n=118 is not a complete-case row (PR #575).

## Do not claim

- Mediation / “CLDN4 explains TACSTD2 immune association”
- Private 8KL / KD co-culture
- Visium same-spot correlation as spatial exclusion
- CosMx TACSTD2 at 50/100 µm as 8/8 (it is not; locked CLDN4 0.36/0.52 stays)

## Reproduce

```bash
export GEO_C4=/tmp/geo_c4
# downloads: scripts/download.sh
python3 scripts/analyze_bridge.py
```
