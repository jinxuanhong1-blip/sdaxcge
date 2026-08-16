# Patient-level meta: malignant/epithelial TACSTD2 and CLDN4 vs T/NK

Public human lung scRNA only. Processed GEO / ArrayExpress. EGA and dbGaP raw are not downloaded.

**GSE207422 User A3 (malignant TACSTD2 vs T/NK) is taken as given.** This folder does not re-annotate or re-audit that slide.

## What is pooled

For each cohort, one patient-level Spearman:

- **x** = TACSTD2 or CLDN4 in malignant (or residual / tumor-site epithelial) cells — mean log1p or mean log1p(CP10k)
- **y** = T/NK fraction of all cells in that patient (CD8 fraction is a sensitivity)

Primary combination: DerSimonian–Laird **random-effects** on Fisher-z(ρ), back-transformed to ρ, with I². Companions: signed Stouffer (weight √(n−3)) and Fisher combined p.

Minimum n per cohort is 4 (Fisher-z variance 1/(n−3)). GSE146100 (n=1) is documented and not pooled.

## Pooled numbers (primary mean-expression metric)

| Gene | Cohorts | N patients | RE ρ [95% CI] | RE p | I² | Stouffer z | Stouffer p |
|---|---:|---:|---|---:|---:|---:|---:|
| TACSTD2 | 8 | 145 | −0.109 [−0.339, +0.133] | 0.378 | 39% | −0.61 | 0.544 |
| CLDN4 | 8 | 145 | −0.137 [−0.306, +0.040] | 0.129 | 0% | −1.53 | 0.125 |

## Per-cohort primary Spearmans

| Cohort | Compartment | n | TACSTD2 ρ | TACSTD2 p | CLDN4 ρ | CLDN4 p |
|---|---|---:|---:|---:|---:|---:|
| GSE207422 (A3 given) | DRMref malignant | 12 | −0.490 | 0.106 | −0.091 | 0.779 |
| GSE205335 | author malignant | 22 | +0.284 | 0.200 | −0.200 | 0.371 |
| GSE241934 IIT | author residual Epi | 11 | −0.073 | 0.832 | −0.064 | 0.853 |
| GSE241934 Real | author residual Epi | 24 | −0.167 | 0.436 | +0.178 | 0.405 |
| GSE291670 | marker malignant | 6 | −0.543 | 0.266 | −0.829 | 0.042 |
| GSE253013 | malig-like tumor | 9 | −0.717 | 0.030 | −0.333 | 0.381 |
| GSE131907 tumor sites | author epithelium | 36 | +0.168 | 0.326 | −0.210 | 0.218 |
| GSE325414 (2026 leftover) | author malignant | 25 | −0.073 | 0.728 | −0.119 | 0.570 |

GSE207422 CLDN4 is **not** on the A3 slide; it is the Spearman of `malig_CLDN4_mean` vs `frac_tnk` on the same given 12-patient table.

## Must-try disposition

| Accession | Disposition |
|---|---|
| GSE207422 | Taken as given (A3). Not re-audited. |
| GSE205335 | Included. Author malignant labels on processed GEO. |
| GSE241934 | Included. Fresh download of processed MTX + author meta; IIT and Real kept separate. |
| GSE291670 | Included. n=6. |
| GSE253013 | Included. n=9 tumor patients. |
| GSE131907 | Included. Tumor sites only (nLung/nLN out). |
| GSE146100 | Skipped. n=1 patient. |
| 2023–2026 leftovers | GSE325414 included (author labels, n=25). Others: blood-only, n=1, bulk FPKM, or 10x MTX without a public cell-type table — see `results/leftover_decisions.tsv`. |

## Run

```bash
python3 methods/scrna_meta_tnk/download.py          # GSE241934 processed GEO
python3 methods/scrna_meta_tnk/analyze_gse241934.py
python3 methods/scrna_meta_tnk/ingest_cohorts.py
python3 methods/scrna_meta_tnk/meta_combine.py
# optional:
python3 methods/scrna_meta_tnk/leftover_search.py
```

## Files

| Path | Role |
|---|---|
| `given_a3.json` | Locked A3 numbers |
| `data/existing/` | Per-patient tables used for non-A3 cohorts |
| `analyze_gse241934.py` | Residual Epi TACSTD2/CLDN4 vs T/NK |
| `ingest_cohorts.py` | Cohort registry of Spearmans |
| `lib_stats.py` / `meta_combine.py` | RE / Stouffer / Fisher + forest |
| `results/cohort_effects.tsv` | All rows (primary + sensitivity + skips) |
| `results/meta_pooled.tsv` | Pooled ρ / p / I² |
| `results/forest_TACSTD2.png` / `forest_CLDN4.png` | Forest plots |
| `results/leftover_decisions.tsv` | Leftover / skip log |
