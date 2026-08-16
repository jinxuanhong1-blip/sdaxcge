# Hunt: public paired pre/post ICB sets where TACSTD2/Tacstd2 rises

Question: are there public **paired** (same subject, pre vs on/post ICB) lung or pan-cancer
mouse/human expression or IHC sets in which TROP2 (`TACSTD2` / `Tacstd2`) goes up after
checkpoint blockade?

User-cited starting points (not treated as ground truth):

- TISMO: “49/64”
- Zhejiang IHC: “94 → 121”

This folder reports **what we could actually recompute from public files**, with
denominators and failure reasons. Scripts live in `scripts/`.

## TISMO (mouse syngeneic, not within-animal pairs)

TISMO in vivo ICB samples are treated vs control **arms** of the same cell line in the
same study. Every ICB-treated sample has `Baseline=0`. This is a model-level
after-vs-without comparison, not a longitudinal biopsy.

Honest counts from `tismo/tismo_summary.json` (gene present; 1,491/1,518 samples mapped
to the expression matrix):

| definition | n | up | frac up | sign-test p |
|---|---:|---:|---:|---:|
| study × cell line × treatment × timepoint, matched control arm | 61 | 43 | 0.705 | 0.0019 |
| same, but Tacstd2 above floor (max arm mean ≥ 0.5 log2) | 38 | 24 | 0.632 | 0.14 |
| study × cell line × treatment (timepoints pooled) | 48 | 33 | 0.688 | 0.013 |
| study × cell line × ICB class | 33 | 24 | 0.727 | 0.014 |
| study × cell line (all ICB pooled) | 31 | 24 | 0.774 | 0.0033 |

We **do not reproduce 49/64**. Closest raw numbers: 64 contrasts were *attempted*, 3 had
no matched control, and 43 of the remaining 61 went up. Median Δ is +0.07 log2; only
7/61 contrasts are nominally p<0.05 up. Lung carcinoma contributes **one** contrast
(GSE155972 LLC anti-PD1+anti-CTLA4, Δ +0.46, MWU p=0.037).

## GEO paired ICB (in progress)

`geo/geo_candidate_series.json` is the E-utilities union (718 GSE: 436 human, 260 mouse).
`scripts/geo_extract.py` then tries a curated seed of published ICB series plus any
extra accessions passed on the command line. A series is counted as confirmatory only
if TACSTD2 is quantified **and** ≥3 subjects have both a pre and an on/post sample.

## Zhejiang IHC (“94 → 121”)

No public patient-level IHC table was found that we can recompute. Closest published
numbers in the literature are different claims (e.g. 94/110 stage III/IV TROP2-positive
in Inomata et al., *Thorac Cancer* 2025; other papers report TROP2 as largely stable
after mixed anti-cancer therapy). Until a deposit or supplement with paired H-scores
is located, this remains a **literature citation, not a recomputed result**.
