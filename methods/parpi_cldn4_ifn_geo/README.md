# PARPi RNA-seq: CLDN4, tight junction, IFN

Scores public PARP-inhibitor RNA-seq for a pre-specified CLDN4 / tight-junction / IFN gene set. The CLDN4-to-NHEJ link (Yamamoto 2022) has no deposited RNA-seq. The DNA-repair-to-IFN link is the scored PARPi matrices plus Ding 2018, Pantelidou 2019, and Sen 2019.

```bash
python3 scripts/score_parpi_cldn4_ifn.py
python3 scripts/plot_contrasts.py
```

Requires pandas, scipy, openpyxl, xlrd, matplotlib. Raw downloads go to `cache/` and are gitignored. Read `FINDING.md` for the numbers and the limits.
