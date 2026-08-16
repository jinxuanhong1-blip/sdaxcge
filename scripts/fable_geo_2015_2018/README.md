# scripts/fable_geo_2015_2018

Reproducible pipeline for the GEO 2015–2018 human-lung ICI / PD-1 / PD-L1 / CTLA-4
mining + TACSTD2/CLDN4 analysis. Pure Python; uses only public NCBI endpoints
(no API key). Run in order from the repository root:

```bash
python3 scripts/fable_geo_2015_2018/01_search.py            # GEO gds search -> 47 candidate GSE
python3 scripts/fable_geo_2015_2018/02_verify.py            # verify EVERY accession (SOFT + FTP listings)
python3 scripts/fable_geo_2015_2018/03_download.py          # download open processed files < 2 GB
python3 scripts/fable_geo_2015_2018/04_build_clinical.py    # per-sample clinical + label flags
python3 scripts/fable_geo_2015_2018/05_target_gene_analysis.py  # TACSTD2/CLDN4 vs response/survival
python3 scripts/fable_geo_2015_2018/06_curate.py            # master_summary.tsv
```

`eutils.py` holds the shared, rate-limited E-utilities helpers.

Dependencies: `pandas`, `numpy`, `scipy`, `matplotlib` (see
`requirements.txt`). All outputs land under `results/fable_geo_2015_2018/`;
the full narrative (EN + 中文) is `notes/fable_geo_2015_2018/WRITEUP.md`.
