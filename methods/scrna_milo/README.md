# methods/scrna_milo

Additive high-end scRNA: Milo-style kNN neighbourhood differential abundance on **public processed UMI**.

User trend taken as given: TROP2-high neighbourhoods should be immune-poor. This folder tests that on GSE207422 (primary) and GSE241934 (if the public MTX is in `data/`).

**miloR is not used.** See `playbook.md` for the documented kNN neighbourhood t-test / Spearman fallback and the SpatialFDR formula.

## Results

- `results/GSE207422/` — post-treatment Hu et al. 2023 UMI
- `results/GSE241934/IIT/` — NEOTIDE EGFR-mutant trial tumours
- `results/GSE241934/RWC/` — real-world neoadjuvant IO tumours (optional)
- `FINDING.md` — honest SpatialFDR / *n* and whether TACSTD2-high neighbourhoods are T/NK-depleted

## Scripts

| File | Role |
|---|---|
| `scripts/knn_nhood.py` | Graph, refined indices, DA, SpatialFDR |
| `scripts/run_gse207422.py` | Stream GEO UMI, marker lineages, run DA |
| `scripts/run_gse241934.py` | Author-labelled MTX cohorts |
| `scripts/download.py` | Public GEO supplementary files only |
