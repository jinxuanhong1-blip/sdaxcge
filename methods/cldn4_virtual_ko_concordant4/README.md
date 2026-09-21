# CLDN4 virtual knockout, concordant-4 malignant cells

Public-data only. Cohorts are GSE123902, GSE131907, GSE205335, and GSE189357.
The private KL matrices are not used.

Stock CellOracle knocks out transcription factors through a motif base GRN.
CLDN4 is not a TF, so that knockout has no outgoing edges. The script that
runs is the same one-step ridge shift CellOracle uses for a regulator, with
CLDN4 entered beside expressed TFs. The fit is inside each patient or sample
that has at least 80 QC-pass malignant cells and CLDN4 detected in at least
5% of those cells. A detection-matched housekeeping gene is knocked out in
the same model as a control.

## Reproduce

```bash
bash methods/cldn4_virtual_ko_concordant4/scripts/download.sh /tmp/concordant4_raw
# GSE205335 is gzip-wrapped twice:
gzip -dc /tmp/concordant4_raw/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz | gzip -dc > /tmp/gse205335_matrix.rds
python3 methods/cldn4_virtual_ko_concordant4/scripts/extract_malignant.py
python3 methods/cldn4_virtual_ko_concordant4/scripts/virtual_ko.py
```

The finished numbers are in `RESULTS.md` and `results/tables/`.
