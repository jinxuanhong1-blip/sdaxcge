# Hunt: ELF3 / GRHL1 / KLF4 / TFAP2A with TACSTD2 and CLDN4, NKX2-1 opposite

Honest public-data search. Co-expression is not evidence of direct transcriptional regulation.

## Claim being tested

In public **human and mouse lung RNA**, the user TF set **ELF3, GRHL1, KLF4, TFAP2A** co-correlates with **TACSTD2** and **CLDN4**, and **NKX2-1** anti-correlates.

## What would count as support (fixed before opening expression matrices)

Pair (raw Spearman on log2(CPM+1), n ≥ 12 after filters):

- positive pair: ρ ≥ 0.30 and two-sided p < 0.05
- NKX2-1 pair: ρ ≤ −0.20 and two-sided p < 0.05

Cohort call:

- **SUPPORTED** — ≥ 6 of 8 positive pairs **and** both NKX2-1 pairs anti
- **PARTIAL** — ≥ 4 of 8 positive pairs (NKX2-1 not required)
- **NOT_SUPPORTED** — otherwise
- **UNINFORMATIVE** — targets/TFs mostly undetected, or n < 12

The same rules are applied after rank-partialling **EPCAM** (epithelial content) and, separately, **SFTPC** (AT2 / alveolar content). A raw SUPPORTED call that becomes NOT_SUPPORTED after EPCAM residualization is **composition-confounded**, not a validated network.

Score (continuous, used for ranking only):  
`mean ρ(4 TFs × 2 targets) − mean ρ(NKX2-1 × 2 targets)`.

A 300-shuffle permutation of TF sample labels gives an empirical p for that score. It is a within-cohort sanity check, not a study-wide FDR.

## Search space

Source: [recount3](https://rna.recount.bio/) open data (human Gencode v26 / G026, mouse Gencode vM23 / M023). Same pipeline for both species.

1. All SRA studies with ≥ 6 runs in the recount3 project index (6,831 human, 8,195 mouse).
2. Per-study SRA metadata text is searched for lung / airway / pulmonary disease terms.
3. Runs are labelled tissue / culture-organoid / cell line / sorted-or-fluid. Single-cell studies are flagged and not downloaded.
4. Download list (pre-specified, not chosen by looking at correlations):
   - **Anchors (always):** GTEx lung, TCGA LUAD, TCGA LUSC
   - **SRA tissue:** ≥ 12 tissue lung RNA-seq runs, single-cell fraction < 0.25, median spots ≥ 5×10^5, and either ≥ 40% of runs are lung-tagged or the lung tag is run-level
   - **SRA culture / cell line:** same n and depth filters; reported separately as composition-controlled tests
   - Sorted cells, blood, BAL, platelets, and isolated immune studies are excluded

This is a screen of public lung RNA that recount3 already processed. It is not every GEO series, and it is not single-cell ATAC or spatial data.

## Why bulk lung can fake this network

TACSTD2 and CLDN4 are high in airway / regenerating epithelium. NKX2-1 and SFTPC are high in AT2 / alveolar lineage. ELF3, GRHL1, KLF4 and TFAP2A are epithelial transcription factors with airway-leaning published roles. In mixed lung tissue the same sample-to-sample swing in **airway vs alveolar fraction** will raise ELF3/GRHL1/KLF4/TFAP2A/TACSTD2/CLDN4 together and lower NKX2-1, without any of those TFs writing TACSTD2 or CLDN4.

That is why EPCAM- and SFTPC-partial correlations, and the culture/cell-line slice, are required before calling the network supported.

## Results

Filled after `scripts/03_fetch_matrices.py` and `scripts/04_analyze.py` finish. See `tables/cohort_scores.tsv` and `tables/pair_correlations.tsv`.

## How to rerun

```bash
# metadata cache (once): scripts/00_fetch_metadata.sh
python3 scripts/01_select_lung_studies.py
python3 scripts/02_choose_cohorts.py
python3 scripts/03_fetch_matrices.py
python3 scripts/04_analyze.py
python3 scripts/05_figures.py
```

Gene IDs are in `scripts/genes.json`.
