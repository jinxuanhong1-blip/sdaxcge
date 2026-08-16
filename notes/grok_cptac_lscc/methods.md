# Methods — CPTAC LSCC TACSTD2 / CLDN4 slice

## Scope
CPTAC lung squamous cell carcinoma only (freeze folder **LSCC**; TCGA synonym **LUSC**). Protein + RNA + immune deconvolution. No GEO ICI cohorts, no LUAD, no spatial, no mouse.

## Identifiers
- TACSTD2 / TROP2: `ENSG00000184292` (matrix row `ENSG00000184292.7`)
- CLDN4: `ENSG00000189143` (matrix row `ENSG00000189143.9`)
Match by Ensembl prefix (version stripped). Same IDs as the LUAD sibling slice.

## Data
Open CPTAC pan-cancer **data_freeze_v1.2_reorganized** on AWS S3 (`cptac-pancancer-data`, us-west-2). Filenames taken from the LinkedOmics CPTAC-pancan-LSCC download table and confirmed HTTP 200. Files >2 GB would be skipped; none of the used matrices exceed that.

## Treatment status
Satpathy et al., Cell 2021 (PMID 34358469): newly diagnosed resected LSCC, **no prior chemotherapy or radiotherapy**. The freeze has **no ICI response labels**. OS/PFS are treatment-naive follow-up, not immunotherapy outcomes.

## Statistics
- Spearman ρ, pairwise-complete; n < 4 → NA
- Mann–Whitney U two-sided + rank-biserial `r = 1 − 2U/(n_a n_b)`; n < 2/group → NA
- Paired Wilcoxon signed-rank for tumor vs overlapping NAT; n < 6 → NA
- BH-FDR within each result table
- RNA/protein signature score = mean of per-gene z-scores; coverage `used/total` reported
- Survival: Kaplan–Meier + log-rank (median split) and Cox PH on the continuous value; underpowered OS/PFS event counts are reported honestly

## How to rerun
```bash
python3 scripts/grok_cptac_lscc/00_download_freeze.py
python3 scripts/grok_cptac_lscc/01_analyze_lscc.py
```
Matrices land in `data/grok_cptac_lscc/` (not committed). Outputs: `results/grok_cptac_lscc/`.
