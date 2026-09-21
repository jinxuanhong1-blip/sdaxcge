# Concordant-4 RNA velocity gate (scVelo)

CLDN4-only. The four public cohorts are GSE123902 (13), GSE131907 (21), GSE205335 (22), and GSE189357 (9). Catalog unit count is 65. This folder asks whether epithelial cells have an RNA-velocity flow toward a CLDN4-high barrier state.

scVelo is run only when a public spliced layer and a public unspliced layer exist. They do not. The script records that gate. It does not estimate velocity and it does not replace velocity with pseudotime.

```bash
python3 methods/concordant4_scrna_velocity/inventory.py
```

Stdlib only. The script reads GEO directory listings, short file prefixes, and SRA runinfo. It does not download the count matrices.
