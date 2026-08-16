# Finding · Milo neighborhood DA on GSE207422

**No neighborhood is differentially abundant at SpatialFDR < 0.1 or < 0.2.** That is the result, not a skipped test.

The user thesis is taken as given. This folder asked whether T/NK-depleted neighborhoods next to TACSTD2-high epithelium are *more abundant* in TACSTD2-high (or MPR) patients. Compositional neighborhoods of that type exist. They do not pass SpatialFDR at this n.

## n (the unit is the sample)

| Contrast | Arm | n samples | Why this n |
|---|---|---|---|
| TACSTD2 high vs low | high / low | **5 vs 5** | Median split of malignant-like mean log1p(CP10k) TACSTD2. 5 samples had <10 malignant-like cells and were unlabeled (P06=1, P11=0, P13=3, P14=0, P15=4). |
| MPR vs NMPR | MPR / NMPR | **4 vs 8** | Post-treatment only. pCR P06 counted as MPR. Three TN biopsies excluded. |

Cells in the graph: 87,841 epithelium+immune / 92,330. Neighborhoods: 6,285 (PCA), 5,971 (Harmony). Median nhood size 66–74.

## SpatialFDR (all testable nhoods)

| Embedding | Contrast | testable | P<0.05 | min P | min SpatialFDR | SpatialFDR<0.1 | SpatialFDR<0.2 | BH<0.1 |
|---|---|---|---|---|---|---|---|---|
| PCA | TACSTD2 high vs low | 3,782 | 384 | 7.7e-4 | 0.378 | **0** | **0** | **0** |
| PCA | MPR vs NMPR | 4,408 | 476 | 2.1e-4 | 0.235 | **0** | **0** | **0** |
| Harmony | TACSTD2 high vs low | 5,246 | 612 | 1.3e-3 | 0.452 | **0** | **0** | **0** |
| Harmony | MPR vs NMPR | 5,609 | 370 | 1.7e-4 | 0.744 | **0** | **0** | **0** |

Nominal P<0.05 is not empty, so the GLM is not stuck at p=1. Multiple-testing + n=5 vs 5 / 4 vs 8 is what removes every hit. Ordinary BH agrees (also 0 at 0.1).

## Kept neighborhoods (composition, not DA)

24 interface nhoods per embedding have TACSTD2-high epithelium and T/NK fraction ≤ the interface 25th percentile (mean T/NK fraction ≈ 0.05). “Next to” = KNN co-membership, not histology.

- Most of those 24 are **not testable** on PCA (only 5–6 have enough samples). Harmony tests 16–21 of them.
- None of the testable kept nhoods have SpatialFDR < 0.1. Best kept SpatialFDR is 0.38 (PCA TACSTD2) / 0.45 (Harmony TACSTD2).
- 4 index cells are kept on **both** embeddings for the TACSTD2 contrast. That is a composition overlap, not a DA replicate.

## What this does not say

- It does not overturn the sample-level A3 null (here ρ=−0.22, p=0.47, n=13 malignant-like TACSTD2 vs T/NK fraction).
- It does not say the thesis is false. It says this public matrix, at neighborhood resolution, with honest SpatialFDR and n=5 vs 5 / 4 vs 8, does not support a DA claim.
- Extreme logFC (±20) appear where one arm is all zeros. Those nhoods still fail SpatialFDR.
- Malignant-like is a marker proxy. Author CopyKAT IDs are not public.
- GSE241934 was not run.

See `da_counts.tsv`, `kept_neighborhoods.tsv`, `summary.json`, and `methods/scrna_milo/README.md`.
