# GSE137244 KL vs KP gene-set scores

Public cell-line RNA-seq from Deng et al., *Nat Cancer* 2021 (PMID 34142094).

```bash
python3 analyze.py
```

Downloads `GSE137244_counts.fpkm.csv.gz` from GEO into `data/` (gitignored) and writes `tables/` and `figures/`.

The contrast is five KL libraries minus five KP libraries. Normal lung is held out. The score is the mean of log2(FPKM+1). The directional test is a one-sided Welch t for KL > KP, with an Epcam gate as a second pass. See `FINDING.md`.
