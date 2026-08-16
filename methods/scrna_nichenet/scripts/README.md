# Scripts

Run from the repo root after `pip install -r methods/scrna_nichenet/env/requirements.txt`.

```bash
python3 methods/scrna_nichenet/scripts/00_download.py
python3 methods/scrna_nichenet/scripts/01_convert_prior.py
python3 methods/scrna_nichenet/scripts/02_extract_panel.py
python3 methods/scrna_nichenet/scripts/03_analyze.py
```

Large files land in `/tmp/scrna_nichenet/` (GEO UMI matrix, NicheNet RDS/parquet).
Committed inputs are DRMref barcodes + the GEO sample sheet under `data/`.
