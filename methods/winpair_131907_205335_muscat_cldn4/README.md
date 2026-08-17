# Winning-pair CLDN4-only patient-pseudobulk DE

GSE131907 + GSE205335 only. CLDN4-only (no dual-high). Patient-level
malignant and T/NK UMI-sum, then TMM + OLS on log2(CPM+1). Not muscat.

```bash
python3 methods/winpair_131907_205335_muscat_cldn4/build_tnk_pseudobulk.py
python3 methods/winpair_131907_205335_muscat_cldn4/analyze.py
```

Writeup: `FINDING.md`. Headline DE table: `tables/de_q4q1_combined_families.tsv`.
