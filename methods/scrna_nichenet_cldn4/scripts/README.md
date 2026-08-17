# Scripts

Run from the repo root after `pip install -r methods/scrna_nichenet_cldn4/env/requirements.txt`.

```bash
python3 methods/scrna_nichenet_cldn4/scripts/00_download.py
python3 methods/scrna_nichenet_cldn4/scripts/01_convert_prior.py
python3 methods/scrna_nichenet_cldn4/scripts/02_extract_panel.py
python3 methods/scrna_nichenet_cldn4/scripts/03_analyze.py
```

Large files land in `/tmp/scrna_nichenet_cldn4/` (GEO UMI matrix, NicheNet RDS/parquet).
Committed inputs are DRMref barcodes + the GEO sample sheet + LR TSV under `data/`.
