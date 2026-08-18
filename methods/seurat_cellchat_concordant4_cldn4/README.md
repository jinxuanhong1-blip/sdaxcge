# Seurat + CellChat — concordant-four CLDN4-high vs low senders

ADDITIVE. Public processed GSE123902 + GSE131907 + GSE205335 + GSE189357.
CLDN4 only. No dual-high. No GSE148071. Patient / locked sample is the unit.
Primary engine is **R + Seurat + CellChat**.

```bash
Rscript methods/seurat_cellchat_concordant4_cldn4/scripts/install_packages.R
bash methods/seurat_cellchat_concordant4_cldn4/scripts/download.sh /tmp/concordant4_raw
Rscript methods/seurat_cellchat_concordant4_cldn4/scripts/run_cellchat.R --raw=/tmp/concordant4_raw
```

See `FINDING.md` for the patient-level ligand call and
`results/ligand_table.tsv` for the table.
