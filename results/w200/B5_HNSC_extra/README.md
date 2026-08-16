# B5: extra public HNSCC ICI cohort

This is a deliberately conservative reanalysis of the HNSCC subset of Prat et
al., *Cancer Research* (2017), PMID
[28487385](https://pubmed.ncbi.nlm.nih.gov/28487385/), GEO
[GSE93157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93157).
The repository had no coverage manifest, so prior coverage could not be checked.

## Reproduce

Python 3 is the only requirement:

```sh
python3 results/w200/B5_HNSC_extra/analyze.py
```

The script extracts the five `HEADNECK` cases from the GEO series matrix.
Objective response is defined prospectively as CR/PR; SD/PD are non-OR. For
each complete gene, it reports the responder-minus-nonresponder mean difference
on the deposited log-scale values and a complete-label permutation p-value
(all 10 allocations of two responders among five patients). PFS association is
an exact Spearman permutation test (all 120 allocations). BH correction is
performed separately for the two screens. The prespecified cytolytic score is
the mean of GZMA and PRF1.

## Files

- `RESULTS.md`: short interpretation
- `patients.csv`: auditable HNSCC clinical subset
- `gene_statistics.csv`: all 725 complete-gene exploratory results
- `cytolytic_score.csv`: per-patient GZMA/PRF1 score
- `summary.json`: machine-readable headline results
- `data/`: unmodified compressed GEO downloads

Source file SHA-256:

```text
acc63e3818d2c92182339db8fac9c69d6c3228000b17e4bf93882004439fab34  GSE93157_series_matrix.txt.gz
220a635a55ea32b771ee0c94238d17906de903f7d2c6ff133d75ba0f5eb67e93  GSE93157_family.soft.gz
```

This dataset is useful as a worked micro-cohort, not as an independent
validation cohort. No model was fitted and no threshold was optimized.
