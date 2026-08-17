# Committed inputs

| File | Source |
| --- | --- |
| `drmref_cell_annotation.tsv.gz` | DRMref public `meta.data` for GSE207422 Tor + Sin (Liu et al., *NAR* 2024). Post-tx barcodes. |
| `geo_scRNAseq_sample_metadata.tsv` | Hu et al. GEO sample sheet (pathologic response, timing). |
| `lr_network.tsv` | NicheNet-v2 human LR network (Zenodo 7074291), converted to TSV. |

Runtime downloads (GEO UMI matrices, ligand–target RDS) stay in `/tmp/multinichenet_cldn4/` and are not committed.
