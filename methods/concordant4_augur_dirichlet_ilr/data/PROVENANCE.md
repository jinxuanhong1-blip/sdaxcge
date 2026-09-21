# Inputs

`locked_patient_units.tsv` is the patient table from the concordant-4 Seurat run (65 units). Malignant CLDN4 percent-positive, quartiles, and the locked T/NK counts are taken from that table and not recomputed.

`GSE205335_gsm_map.tsv` maps `orig.ident` to patient and tissue so Normal* samples can be dropped the same way as that run.

`gene_panel.txt` is the shared activity panel for the Augur-style classifier. CLDN4 is listed so extractors can confirm it was measured, then every classifier drops it.

Count matrices are downloaded to `/tmp/geo_c4` by `scripts/download.sh` and are not committed.
