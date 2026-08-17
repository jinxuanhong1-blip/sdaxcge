# Rerun — LinkedOmics CLDN4 OS/PFS

```bash
python3 -m pip install -r methods/linkedomics_cldn4_surv/scripts/requirements.txt
python3 methods/linkedomics_cldn4_surv/scripts/00_download.py
python3 methods/linkedomics_cldn4_surv/scripts/01_analyze.py
```

`00_download.py` HEADs each URL and refuses files that are not HTTP 200 or are >2 GB.
Local matrices stay in `methods/linkedomics_cldn4_surv/data/` (gitignored).
Do not download phenotype / ImmuneScore files for this slice.

Expected outputs:

- `methods/linkedomics_cldn4_surv/FINDING.md`
- `methods/linkedomics_cldn4_surv/results/survival.tsv`
- `methods/linkedomics_cldn4_surv/results/figures/km_cldn4_os_pfs.png`
