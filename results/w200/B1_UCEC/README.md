# B1 analog — TCGA-UCEC: rank CLDN4 among surface genes coexpressed with TACSTD2

Honest ranking. No invented statistics. Primary metric is Spearman.

## Answer

**No.** In TCGA-UCEC primary tumours, `CLDN4` is **not** the top surfaceome partner of `TACSTD2`.

| setting | n | universe | CLDN4 Spearman rank | rho | actual #1 |
| --- | --- | --- | --- | --- | --- |
| Primary (all -01 samples; symbol + Ensembl fallback) | 549 | 2751 | **#22** (99.237th pct) | 0.314 | NECTIN4 (0.465) |
| Symbol-only (2018 UniProt names; BRCA-identical rule) | 549 | 2678 | #21 | 0.314 | SLC52A3 (0.412) |
| One sample per patient | 545 | 2751 | #25 | 0.312 | NECTIN4 (0.462) |

Pearson (secondary) puts CLDN4 at **#4** (r = 0.442). That does not change the Spearman-primary verdict.

Independent `scipy.stats` check of the TACSTD2–CLDN4 pair: Spearman 0.313860 (p = 5.142e-14); Pearson 0.442489 (p = 1.003e-27). Matches the ranked table.

## Files

| file | what |
| --- | --- |
| `summary.md` | human-readable answer + top 25 |
| `summary.json` | machine-readable ranks, mapping, sensitivities |
| `coexpression_TACSTD2_surfaceome.csv` | full ranked surfaceome |
| `top200.csv` | reporting window (`w200`) |
| `scatter_TACSTD2_vs_CLDN4.png` | pair scatter |
| `top_partners_TACSTD2.png` | top-20 bar chart (CLDN4 highlighted if in window) |

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/w200/B1_UCEC/download_data.py
python3 scripts/w200/B1_UCEC/coexpression.py
```

Inputs (git-ignored under `data/`): Xena GDC `TCGA-UCEC.star_fpkm-uq.tsv.gz`, GENCODE v36 probeMap, Bausch-Fluck 2018 table S3 (steveneschrich/surfaceome copy).
