# Methods — TACSTD2–CLDN4 same-cell BRIDGE

## Question

Do TACSTD2 and CLDN4 co-occur in the **same** malignant/tumor cells above
independence, and do dual-high cells sit in a colder immune niche?
Frame as **bridge only**, not mediation.

## Data

| Layer | Source | Cell definition |
|---|---|---|
| Concordant-4 | GSE123902, GSE131907, GSE205335, GSE189357 | Marker-malignant (EPCAM/KRT & PTPRC−) for 123902/189357; author Malignant for 131907/205335. Locked 65 units. |
| CosMx | figshare 25976224 He2022 clustered h5ad | Author-matched tumor label per section (tumor 5/6/9/12/13), same map as PR #726. n_tumor = 295,877. |
| DepMap protein | Gygi CCLE MS `_LUNG` complete cases | Pairwise TACSTD2–CLDN4 (n=45). |

## Statistics

- Detection: raw UMI/count > 0.
- Independence ratio: P(TACSTD2+ ∩ CLDN4+) / [P(TACSTD2+) · P(CLDN4+)].
- Odds ratio + Fisher exact on the 2×2 detection table.
- Spearman: concordant-4 on raw UMI; CosMx on log1p CP10k (to match PR #726 inventory).
- Unit for concordant-4 summaries: patient/donor/sample with ≥30 malignant cells.
- Cold niche: CosMx median-quadrant immune-neighbor fraction at 10 µm from PR #726 (not re-fit).

## Self-checks

- Concordant-4 median within-cell TACSTD2–CLDN4 ρ matches PR #638 (0.746).
- CosMx section Spearmans match PR #726 inventory (max |Δ| ≤ 1e−6).
- DepMap protein ρ matches PR #575 (0.693).

## Out of scope

- Mediation / partial attenuation (see PR #718 null; PR #726 CosMx dependence).
- Locked CLDN4 CosMx 0.36/0.52 and concordant-4 ρ=−0.531 (quoted only).
