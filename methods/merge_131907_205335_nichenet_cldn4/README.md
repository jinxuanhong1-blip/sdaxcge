# Merged GSE131907 + GSE205335 NicheNet-style CLDN4 ligand activity

ADDITIVE. CLDN4-only. No dual-high. Does **not** redo GSE207422 NicheNet.

```bash
python3 methods/merge_131907_205335_nichenet_cldn4/scripts/00_download.py
python3 methods/merge_131907_205335_nichenet_cldn4/scripts/01_convert_prior.py
python3 methods/merge_131907_205335_nichenet_cldn4/scripts/02_extract.py
python3 methods/merge_131907_205335_nichenet_cldn4/scripts/03_analyze.py
```

Primary table: `results/ligand_activity_table.tsv`.
Writeup: `FINDING.md`.
