# Methods — Harmony n=73 CLDN4-only high-end LR

## Given (not re-run)

Harmony naive joint object (`methods/scrna_harmony_naive`): GSE131907 + GSE253013 + GSE148071 + GSE127465.
Malignant-like *CLDN4* vs T/NK **n=73, ρ=−0.27**. Donor table copied to `data/given/per_donor_metrics.tsv`.

## Why per-cohort LR + meta

The joint Harmony embedding is not stored. Rebuilding it needs the GSE253013 Garnett RDS (~9 GB), which is not downloaded. Ligand–receptor scores are therefore computed on each contributing cohort that has a public processed matrix, then **patient deltas** are meta-analyzed (one FINDING).

## Lineage / malignant-like

Same rule as the given Harmony table: marker-argmax on the lineage panel; malignant-like = epithelial and near-zero normal-lung score (`SFTPA2/SFTPC/AGER/SCGB1A1/SCGB3A1/TPPP3/FOXJ1`). TACSTD2 is not a gate.

## Patient split

Within each given donor, malignant-like cells are split at that donor’s median *CLDN4* (log1p CP10k, or log1p of deposited normalized counts for GSE127465). A donor enters the paired test if it has ≥8 cells in **both** bins and ≥20 T/NK cells.

## Scores

- **LIANA-style (primary):** minimum subunit mean; pair score = mean of the two partner means (CellPhoneDB / LIANA documentation).
- **CellChat-style:** 10% truncated means, geometric mean of subunits, Hill P = LR / (0.5 + LR). CellChat R is not run.

Outgoing direction is malignant → that patient’s own T/NK.

## Meta

Wilcoxon signed-rank on patient Δ (high−low). BH-FDR within method. Optional DerSimonian–Laird on cohort mean Δ when ≥2 cohorts have ≥4 patients.

## Not run

- GSE253013 cell-level LR (9 GB RDS).
- Milo-style neighborhoods (no stored joint embedding).
- Dual-high TACSTD2+CLDN4.
