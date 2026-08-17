# scripts

| Script | Role |
| --- | --- |
| `extract_gse253013.py` | Stream the 9.3 GB GSE253013 RDS; keep the LR gene panel |
| `assemble_gse253013.py` | CSC row extract + patient/tissue metadata |
| `convert_prior.py` | NicheNet-v2 ligand–target RDS → parquet |
| `analyze.py` | CLDN4-only CellChat-style + NicheNet-style on the n=15 combo |
| `gene_sets.py` | Lineage, IFN, cytotoxicity, ligand classes |
