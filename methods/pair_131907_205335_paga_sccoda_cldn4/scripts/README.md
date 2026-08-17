# Scripts

1. `run_sccoda.py` — table-only composition (no GEO download).
2. `download.py` — GSE131907 raw UMI + GSE205335 UMI RDS. Skips TPM / GSE148071.
3. `extract.py` — per-unit subsample of malignant + T/NK, write h5ad.
4. `run_paga.py` — Harmony/PAGA/DPT on the subsample.
5. `write_finding.py` — fill `FINDING.md` from tables.
6. `run_all.py` — 1 then 2–5.
