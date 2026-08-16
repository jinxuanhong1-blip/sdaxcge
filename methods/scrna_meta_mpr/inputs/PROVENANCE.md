# Input provenance

All rows are from public GEO processed objects. No unpublished barcode tables.

| File | Source |
|---|---|
| GSE207422_per_sample.tsv | Public GSE207422 UMI matrix + GEO sample sheet (Hu 2023). Epithelial and malignant-like means. |
| GSE207422_sample_metadata.tsv | Same GEO sample sheet (MPR / NMPR / TN). |
| GSE241934_scrna_per_sample.csv | Public IIT/Real MTX + author Meta (`major.cell.type`, Pathological Response). TACSTD2. |
| GSE241934_cldn4_per_sample.csv | Same MTX streamed for CLDN4 in author Epi. |
| GSE291670_per_sample.tsv | Public 10x MTX in GSE291670_RAW.tar; GEO titles MPR vs Non-MPR. |
| GSE205335_per_patient_tacstd2.csv | Public RDS + author CellIdentity (`Malignant cells`) + GEO RECIST. |
| GSE205335_per_patient_cldn4.csv | Same RDS, CLDN4 row, same patient/RECIST map. |
| GSE205335_sample_map.csv | orig.ident → patient / RECIST from GEO characteristics. |
| prior_hunt_inventory.tsv | Earlier public-matrix inventory (leftover bookkeeping). |
