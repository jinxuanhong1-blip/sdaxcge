# GSE248249 — CLDN4 / TACSTD2 vs pre/post acquired IO resistance

Public **array**, not NanoString. Memon et al., *Cancer Cell* 2024 ([PMID 38215748](https://pubmed.ncbi.nlm.nih.gov/38215748/)). GEO [GSE248249](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248249), platform [GPL23126](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL23126) Affymetrix Clariom D Human (transcript/gene version). Values = RMA from the deposited series matrix (Affymetrix Expression Console).

**CLDN4 is on the panel.** One official-symbol transcript cluster: `TC0700007993.hg.1` (`NM_001305 // CLDN4`). TACSTD2 is also on the panel: `TC0100014340.hg.1` (`NM_002353 // TACSTD2`). Each gene has exactly one official-symbol probe on GPL23126. Matrix feature count = **138,745**.

## What can be tested

GEO annotates **timepoint** (Pre-treatment / Post-treatment), patient, tumor site, sex. It does **not** annotate RECIST, responder, or a sensitive vs resistant label.

The series design (GEO `!Series_overall_design`) is: 29 NSCLC patients treated with PD-1 blockade; **13** FFPE tumors before therapy; **29** FFPE tumors at acquired resistance. Every patient in this molecular set has acquired resistance. The only public contrast is **pre vs post AR**, not sensitive vs resistant.

## n (counted from the series matrix)

| | Count |
|---|---:|
| Samples | 42 |
| Patients | 29 |
| Pre-treatment | 13 |
| Post-treatment (acquired resistance) | 29 |
| Patients with both timepoints | 13 |
| Same anatomic site in the pair | 4 (Pts 01 lung, 12 adrenal, 19 adrenal, 23 lymph node) |
| Sensitive vs resistant in GEO | **0** (no such field) |

Primary unit = sample for the unpaired test; patient for the paired test. Sixteen post-only patients have no pre sample.

## Scores (RMA; no invented numbers)

Two-sided Mann–Whitney U (unpaired) and Wilcoxon signed-rank (paired). Cliff’s δ = P(post>pre) − P(post<pre).

| Gene | Contrast | n | Pre median | Post median | Δ median (post−pre) | Stat | p |
|---|---|---|---:|---:|---:|---|---:|
| **CLDN4** | unpaired post vs pre | 13 vs 29 | 5.66 | 4.83 | −0.83 | U=111; δ=−0.41 | **0.036** |
| **CLDN4** | paired post−pre | 13 pairs | 5.66 | 4.83 | −0.89 | W=12; 11/13 post lower | **0.017** |
| **CLDN4** | paired, same site only | 4 pairs | 5.20 | 4.74 | −0.74 | W=0; 4/4 post lower | **0.125** |
| **TACSTD2** | unpaired post vs pre | 13 vs 29 | 4.84 | 5.07 | +0.24 | U=205; δ=+0.09 | **0.66** |
| **TACSTD2** | paired post−pre | 13 pairs | 4.84 | 4.75 | −0.16 | W=44; 6 up / 7 down | **0.95** |
| **TACSTD2** | paired, same site only | 4 pairs | 4.43 | 4.08 | −0.02 | W=5; 2 up / 2 down | **1.00** |

Means (same n): CLDN4 unpaired pre 5.79 vs post 5.04 (Δ −0.75); paired Δ mean −0.74. TACSTD2 unpaired pre 5.29 vs post 5.65 (Δ +0.36); paired Δ mean +0.09.

Same-site n=4 cannot give a two-sided Wilcoxon p below 0.125. Direction for CLDN4 matches the n=13 paired test; it is not an independent confirmation.

## CLDN4 vs TACSTD2 (sample-level Spearman)

| Stratum | n | ρ | p |
|---|---:|---:|---:|
| All samples | 42 | +0.55 | 1.4×10⁻⁴ |
| Pre | 13 | +0.66 | 0.013 |
| Post | 29 | +0.55 | 0.0019 |

## Read this as

- **CLDN4 is present and scoreable.** It is **lower at acquired-resistance (post)** than at pre-treatment in this public matrix (unpaired p=0.036; paired 11/13 down, p=0.017). That is a two-gene, pre-specified contrast on n=13 pairs, not a genome-wide claim.
- **TACSTD2 does not move** on the same samples (unpaired p=0.66; paired p=0.95).
- **Sensitive vs resistant cannot be scored.** GEO has no such labels. Do not treat post as “resistant vs sensitive.”
- **Site is a real confounder.** Only 4/13 pairs are the same organ. Several pairs jump lung↔liver/bone/spine/adrenal. Bulk FFPE RMA is also purity-mixed; this analysis does not adjust purity.
- Immune genes the paper highlights (GZMA, CXCL9, B2M) trend higher post on the same matrix (unpaired medians 4.49→5.96, 9.33→11.78, 9.58→11.33) but unpaired p=0.073 / 0.25 / 0.16. They are a processing check, not a target result.

## Files

- `analyze.py` — download the series matrix to `/tmp/gse248249/GSE248249_series_matrix.txt.gz`, then run
- `results/pheno.tsv` / `sample_level.tsv` / `tests.tsv` / `paired_sites.tsv` / `summary.json`
