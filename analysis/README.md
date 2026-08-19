# CLDN4 tumor-geography pipeline

Official CosMx NSCLC (NanoString/Bruker S3) and GSE307534 Visium LUAD (invasive titles only).

```bash
python3 analysis/download_official_data.py --out /tmp/spatial_data
python3 analysis/cldn4_tumor_geography.py --data /tmp/spatial_data --out results --n-perm 199
```

CLDN4 is the only fence gene. Outputs: `RESULTS.md`, infiltration-depth curves, barrier index, rotate/shift permutation.
