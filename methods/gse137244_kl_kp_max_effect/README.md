# GSE137244 KL vs KP: Epcam and leave-one filters (n = 5 vs 5)

Public cell-line RNA-seq from Deng et al. (GEO [GSE137244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137244)). Five KrasG12D;Lkb1 lines versus five KrasG12D;Trp53 lines. `normal-lung-RNA` is not a tumor cell line and is excluded.

The script downloads `GSE137244_counts.fpkm.csv.gz` from the NCBI FTP, then writes `tables/`, `figures/`, and `FINDING.md`. Primary numbers are computed there. They are not typed in by hand.

```bash
python3 methods/gse137244_kl_kp_max_effect/analyze.py
```

Scale is log2(FPKM+1), the same scale as the locked Tacstd2 and Cldn4 contrast. No library is removed from that 5-vs-5 contrast. Leave-one-library rows are reported with their reduced n and are not the result.
