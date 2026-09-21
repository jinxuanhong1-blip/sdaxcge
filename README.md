# Public TROP2 / TACSTD2 omics

Reanalysis of public TROP2-ADC and TACSTD2 knockdown/knockout matrices for CLDN4 and interferon / effector programs.

- Catalog and numbers: [RESULTS.md](RESULTS.md)
- Script: [analysis/reanalyze_public_trop2_cldn4_immune.py](analysis/reanalyze_public_trop2_cldn4_immune.py)
- Tables: [results/](results/)

```bash
pip install -r requirements.txt
python3 analysis/reanalyze_public_trop2_cldn4_immune.py
```

The script downloads GEO supplementary matrices to `/tmp/trop2_data` and rewrites `results/`.
