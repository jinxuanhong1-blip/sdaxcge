# Concordant-4 CD8/NK state vs malignant CLDN4

Patient-level effector and exhaustion scores in CD8 cells and in NK cells, high vs low malignant CLDN4, on the locked concordant-4 units only.

Exposure (malignant CLDN4 % positive, within-cohort quartile) is the table from the concordant-4 patient analysis. This folder does not re-estimate the T/NK fraction.

```bash
bash methods/concordant4_cd8nk_state/scripts/download.sh /tmp/concordant4_raw
gzip -dc /tmp/concordant4_raw/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz \
  | gzip -dc > /tmp/concordant4_raw/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds
Rscript methods/concordant4_cd8nk_state/scripts/export_gse205335.R
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE123902
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE189357
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE131907
python3 methods/concordant4_cd8nk_state/scripts/extract_scores.py --dataset GSE205335 \
  --cells205 /tmp/concordant4_state/gse205335_cells.tsv.gz
python3 methods/concordant4_cd8nk_state/scripts/analyze_state.py
```

Write-up: `FINDING.md`. Methods: `METHODS.md`.
