# FINDING — merged GSE131907 + GSE205335 LIANA/LR from CLDN4-high malignant to T/NK

**Status:** pipeline is on the branch. `results/lr_table.tsv` is written by
`scripts/01_analyze.py` after the processed UMI matrices are downloaded.

PR #320 Q4 vs Q1 author %pos vs T/NK on this slice (**n=23, r=−0.705**) is
given and is not re-audited. This folder tests outgoing CLDN4-high vs
CLDN4-low malignant → same-patient T/NK (patient-level paired Wilcoxon).
No dual-high. No GSE207422.

Reproduce:

```bash
python3 methods/merge_131907_205335_liana_cldn4/scripts/00_download.py
python3 methods/merge_131907_205335_liana_cldn4/scripts/01_analyze.py
```
