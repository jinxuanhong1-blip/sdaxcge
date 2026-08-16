# w200 HPA lung TACSTD2 / CLDN4

Re-download official Human Protein Atlas tables and write honest public slices to `results/w200/HPA_lung/`.

```bash
python3 scripts/w200_hpa_lung/download_and_extract.py
```

Requires network access to `proteinatlas.org` and `v23.proteinatlas.org`. Large full TSVs stay in `results/w200/HPA_lung/.cache/` (gitignored). Only gene-filtered tables are committed.
