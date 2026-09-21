# Yamamoto OV mutation burden and PARPi interferon NES

Public TCGA-OV test of Yamamoto et al., Mol Cancer Ther 2022, Fig. 4F, and preranked GSEA of acute PARP-inhibitor RNA-seq in ovarian models.

Raw matrices are downloaded into `OV_MAX_DATA` (default `/tmp/ov_max`) and are not committed.

```bash
pip install -r methods/yamamoto_ov_parpi_max/requirements.txt
python3 methods/yamamoto_ov_parpi_max/scripts/max_rho_mutation_cldn4.py
python3 methods/yamamoto_ov_parpi_max/scripts/max_isg_nes_parpi.py
```

The rho script expects the cBioPortal datahub OV files, the GDC ABSOLUTE purity file, the PanCan clinical sample table, and `TCGA-OV.star_tpm.tsv.gz` under that cache. The NES script expects the GEO count matrices listed in `FINDING.md`. Paths are the ones used to produce `results/`.

Numbers and the sample filters are in `FINDING.md`.
