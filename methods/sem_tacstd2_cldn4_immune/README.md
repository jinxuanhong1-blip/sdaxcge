# TACSTD2, CLDN4, and immune: three linear graphs

Compares three Gaussian SEMs on concordant-4 (n = 65) and on TCGA primary tumors:

- TACSTD2 → CLDN4 → immune
- CLDN4 → TACSTD2 → immune
- independent (no TACSTD2–CLDN4 edge; both point at immune)

The writeup is `FINDING.md`. The definitions are `METHODS.md`.

```bash
pip install -r methods/sem_tacstd2_cldn4_immune/requirements.txt
python3 methods/sem_tacstd2_cldn4_immune/test_semcore.py
python3 methods/sem_tacstd2_cldn4_immune/scripts/analyze.py
```

`analyze.py` reads the two tables in `data/` and rewrites `results/` and `FINDING.md`. Rebuilding the TCGA extract from the Xena GDC hub:

```bash
python3 methods/sem_tacstd2_cldn4_immune/scripts/extract_tcga.py
```
