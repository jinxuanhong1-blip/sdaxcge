# TROP2 (TACSTD2) / CLDN4 protein and spatial vs immune contexture

Reproducible downloads and stats for CPTAC LUAD/LSCC protein, PXD042091, PXD059688, GSE271689 GeoMx WTA, and E-MTAB-13530 Visium WTA.

```bash
python3 -m pip install pandas numpy scipy matplotlib seaborn statsmodels lifelines h5py
bash run_all.sh
```

Local caches go to `/workspace/data/protein_spatial` (not committed). Raw MS/FASTQ and files >2 GB are skipped. CPTAC cohorts are treatment-naive (no ICI labels).
