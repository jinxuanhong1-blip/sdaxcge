# Winning-pair DoRothEA TF activity (CLDN4-only)

ADDITIVE **CLDN4 only** on **GSE131907 + GSE205335** (PR #320 winning
pair). No dual-high. GSE207422 and GSE148071 are not merged in.

Patient-level malignant CLDN4-high vs low (Q4 vs Q1) for IFN / MHC /
TJ / keratin TF activities. decoupleR is not installed; scoring is
documented DoRothEA wmean (`METHODS.md`).

```bash
python3 methods/winpair_131907_205335_dorothea_cldn4/scripts/download.py
python3 methods/winpair_131907_205335_dorothea_cldn4/scripts/analyze.py
```

Primary writeup: `FINDING.md`. TF table: `results/tf_table.tsv`.
