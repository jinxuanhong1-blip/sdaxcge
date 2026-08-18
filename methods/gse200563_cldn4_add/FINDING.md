# Finding — GSE200563 CLDN4-only add-if-same-sign

ADDITIVE. **CLDN4 only.** No dual-high TACSTD2×CLDN4. No GSE148071.
Public matrix only: `GSE200563_processed_data.txt.gz`.

This accession is **not** Wu 2021 scRNA (that is GSE148071, excluded).
GSE200563 is Zhang/Abdo 2022 GeoMx DSP of NSCLC primary + paired brain
metastasis (PMID 36216799). PR #459 skipped it as “spatial, not scRNA”.
This folder tests the **compartment analog**: tumor-core (L/LB) CLDN4
vs same-patient TIME (TIME-L / TIME-B) T/NK. Honest n = patients.

## Decision: GSE200563 **does not join** the concordant pool

T/NK Spearman sign **flips** relative to the concordant sets (ρ=+0.129, n=16, p=0.633) and tumor-cell IFN/MHC is **not** down in CLDN4-high (IFN Δ=0.044, MHC Δ=0.649). Stop after the solo table. GSE200563 does not join.

## Solo table (this accession)

| item | n / effect | note |
|---|---|---|
| catalog patients (exclude BC) | **35** | numeric IDs on ROI names |
| patients with L or LB tumor core | **34** | epithelial/malignant analog |
| patients with tumor core **and** TIME ROI | **16** | Spearman unit |
| CLDN4 present | yes | Q3 range 32.1–1387.5 |
| TACSTD2 present | yes | audit only; never a gate |
| T/NK genes used | CD3D, CD3E, CD3G, CD2, CD8A, CD8B, NKG7, GNLY, KLRD1 | mean log2(Q3+1) |
| tumor CLDN4 vs TIME T/NK Spearman | n=16 ρ=0.129 p=0.633 | primary solo |
| same, Q4 vs Q1 T/NK | n_Q1/Q4=4/4 r=0.250 p=0.686 | same 16 patients; tails thin |

Patient is the unit. Tumor-core = L (primary) and/or LB (brain met).
TIME = CD45/immune AOIs labeled TIME-L / TIME-B. mLN and TBME are
sensitivity only. BC controls are dropped. Values are depositor Q3.

## Per-patient TIME-paired units

| patient | L | LB | TIME-L | TIME-B | CLDN4 tumor | T/NK TIME | IFN tumor | MHC tumor | TJ tumor |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 05 | 1 | 1 | 1 | 1 | 9.574 | 6.102 | 6.566 | 8.227 | 6.466 |
| 12 | 1 | 1 | 1 | 1 | 8.401 | 6.719 | 6.697 | 8.471 | 6.427 |
| 14 | 1 | 1 | 0 | 1 | 9.102 | 6.144 | 6.927 | 8.543 | 6.496 |
| 15 | 1 | 1 | 1 | 1 | 9.451 | 6.374 | 6.328 | 7.004 | 6.456 |
| 18 | 1 | 1 | 1 | 0 | 8.222 | 5.727 | 6.379 | 7.727 | 6.422 |
| 19 | 1 | 1 | 1 | 1 | 8.569 | 6.396 | 6.309 | 7.383 | 6.484 |
| 20 | 1 | 1 | 0 | 1 | 8.171 | 5.874 | 6.315 | 7.487 | 6.440 |
| 24 | 1 | 1 | 1 | 0 | 7.516 | 6.525 | 6.665 | 8.015 | 6.501 |
| 26 | 1 | 0 | 1 | 0 | 10.073 | 6.739 | 6.601 | 7.807 | 6.460 |
| 29 | 1 | 0 | 1 | 0 | 7.798 | 6.274 | 6.863 | 8.543 | 6.523 |
| 30 | 1 | 0 | 1 | 0 | 7.615 | 6.303 | 6.844 | 8.915 | 6.396 |
| 31 | 0 | 1 | 0 | 1 | 8.529 | 5.669 | 6.308 | 8.518 | 6.500 |
| 32 | 1 | 0 | 1 | 0 | 10.439 | 6.772 | 6.926 | 9.020 | 6.500 |
| 35 | 1 | 1 | 1 | 1 | 8.507 | 6.457 | 6.323 | 7.297 | 6.278 |
| 40 | 1 | 1 | 1 | 0 | 7.987 | 6.622 | 6.437 | 7.088 | 6.402 |
| 43 | 1 | 0 | 1 | 0 | 7.942 | 6.547 | 6.592 | 7.461 | 6.507 |

Honest Spearman n = **16** (not 35 cases, not 120 ROIs, not 109 paper ROIs).

## Sensitivity T/NK (not the join key unless primary n is unusable)

| contrast | n | ρ | p |
|---|---:|---:|---:|
| L vs TIME-L (lung pair) | 13 | 0.368 | 0.216 |
| LB vs TIME-B (brain-met pair) | 8 | 0.405 | 0.320 |
| tumor-core CLDN4 vs tumor-core T/NK (same ROI) | 34 | -0.119 | 0.504 |
| tumor+mLN CLDN4 vs TIME T/NK | 16 | 0.206 | 0.444 |

## Tumor-core IFN / MHC / TJ  (Q4 vs Q1 if n_units ≥ 8)

Units = patients with L or LB (n=34). Family score = mean log2(Q3+1)
of A8 genes present on the matrix. TJ holds CLDN4 out. Positive delta =
higher in CLDN4-high.

| family | n | n_Q1 / n_Q4 | Δ median (Q4−Q1) | rank-biserial r | p | Spearman vs CLDN4 ρ (p) |
|---|---:|---|---:|---:|---:|---|
| IFN | 34 | 9/9 | 0.044 | 0.111 | 0.724 | 0.141 (0.426) |
| MHC-I/APM | 34 | 9/9 | 0.649 | 0.259 | 0.377 | 0.147 (0.406) |
| TJ (CLDN4 held out) | 34 | 9/9 | 0.015 | 0.086 | 0.791 | 0.160 (0.367) |

IFN/MHC **down** in CLDN4-high means Q4−Q1 delta < 0 or continuous ρ < 0.

## Concordant-4 locked pool (PR #459; not re-audited as singles)

Members: GSE123902 + GSE131907 + GSE205335 + GSE189357.
Score = malignant CLDN4 **mean** vs same-unit T/NK fraction.
Pooled ρ = DerSimonian–Laird Fisher-z (same as PR #459).

| cohort | n | ρ | p |
|---|---:|---:|---:|
| GSE123902 | 13 | -0.654 | 0.015 |
| GSE131907 | 21 | -0.396 | 0.075 |
| GSE205335 | 22 | -0.200 | 0.371 |
| GSE189357 | 9 | -0.517 | 0.154 |
| **concordant-4** | **65** | **-0.403** | 0.002 (I²=0%) |

## Stacked pool

Not computed. Join rule failed (sign flip / IFN-MHC not down / n unusable).

## Methods (locked)

- CLDN4 only. TACSTD2 is an audit gene, never a gate.
- Matrix: GEO `GSE200563_processed_data.txt.gz` only (depositor Q3).
- Tumor-core ROIs = columns `L##` and `LB##`. TIME ROIs = `TIME-L##`, `TIME-B##`.
- Patient ID = the numeric suffix. Replicates (a/b) are averaged.
- CLDN4 score = mean log2(Q3+1) across a patient's tumor-core ROIs.
- T/NK score = mean log2(Q3+1) of CD3D/CD3E/CD3G/CD2/CD8A/CD8B/NKG7/GNLY/KLRD1
  on TIME ROIs (same genes as the PR #459 marker gate, scored as a mean).
- Spearman requires n≥4 finite pairs. Q4 vs Q1 requires n_units≥8 and two tails.
- Family scores use A8 Hallmark IFN-α/γ, custom MHC-I/APM, KEGG+GOBP TJ (CLDN4 held out).
- Join if solo T/NK ρ < 0 **or** tumor IFN or MHC Q4−Q1 delta < 0 (n≥8).
- Stop if n<4 for T/NK **and** IFN/MHC Q4 is unusable, or if T/NK ρ>0 and IFN/MHC not down.
- Stacked ρ = Fisher-z DerSimonian–Laird of cohort Spearmans. Honest n = sum of units.
- Concordant-4 vectors are the PR #459 locked tables. GSE148071 is not used.
- Reproduce: `python3 methods/gse200563_cldn4_add/analyze.py`

## How to read this

- This is a **GeoMx patient analog**, not a single-cell T/NK fraction.
- Do not quote n=35 cases or n=120 ROIs as the Spearman n.
- Do not call this Wu scRNA. Do not add GSE148071.
- No dual-high. No CellChat.

