# opus_tls — TACSTD2 / CLDN4 vs TLS in public lung RNA

```
export OPUS_TLS_DATA=/tmp/opus_tls_data
bash 00_download.sh
python3 01_prep_tcga.py
python3 02_prep_cohorts.py
python3 03_core_association.py   # ~few minutes; writes the main tables
python3 04_finish_and_outcomes.py
python3 05_figures.py
```

Raw matrices stay in `$OPUS_TLS_DATA` and are not committed.
Outputs: `results/opus_tls/{tables,figures}/`.
Write-up: `notes/opus_tls/WRITEUP.md`.
