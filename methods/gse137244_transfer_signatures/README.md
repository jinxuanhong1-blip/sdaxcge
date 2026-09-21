# GSE137244 transferable KL vs KP signatures

Public GEO RNA-seq of nodule-derived cell lines (Deng et al., Nature Cancer 2021, PMID 34142094). Five KL libraries versus five KP libraries. Normal lung is excluded.

```bash
python3 analyze.py
```

The script downloads the three GEO supplementary count files into `$GSE137244_DIR` (default `/tmp/gse137244`) and writes `tables/` and `figures/`.

Score the GMT lists on another mouse scRNA matrix as the mean of log-normalized expression of the genes present in both. Do not re-cut the lists on the query. The KL-versus-KP direction of each list on this series is in `FINDING.md`.
