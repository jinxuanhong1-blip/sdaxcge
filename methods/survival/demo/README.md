# Demo: TACSTD2 / CLDN4 vs PFS on two public GEO ICI cohorts

Open data only. No fabricated statistics.

| Cohort | Treatment | Endpoint deposited | n | Events (computed) | TACSTD2 / CLDN4 |
|---|---|---|---|---|---|
| [GSE135222](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222) | anti–PD-1/PD-L1 NSCLC | PFS (days) | 27 | 21 | RNA-seq TPM |
| [GSE190265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190265) | anti–PD-1 NSCLC | PFS (months) | 43 | 35 | RNA-seq TPM |

Neither series deposits OS. Continuous Cox for TACSTD2 and CLDN4 is null on
both clocks (every 95% CI includes 1). See `../playbook.md` §9 for why that
cannot establish a biomarker, §10 for why OAK/POPLAR are the right next
dataset and are not used here, and §13 for the computed numbers.

```bash
python3 01_fetch_geo.py          # downloads GEO files into data/raw/ (gitignored)
python3 02_analysis.py           # writes results/demo_results.md and figures
```

Committed files under `data/` are the small harmonised tables so the analysis can be re-run offline after the first fetch.
