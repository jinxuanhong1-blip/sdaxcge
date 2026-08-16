# GSE131907 no-skip analysis

Lung adenocarcinoma scRNA atlas (Kim et al. 2020, GSE131907).
Question: TACSTD2 in author epithelial cells versus immune cells.

Size is not a skip reason. Both processed GEO text matrices are used:

- `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`
- `GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz` (2.9 GB; previously skipped)

RDS copies are the same matrices. This pipeline reads the text supplements
because it does not depend on R.

```bash
bash analysis/01_download.sh
python3 analysis/02_extract.py
python3 analysis/03_analyze.py
```

Outputs: `results/noskip/GSE131907/`.
