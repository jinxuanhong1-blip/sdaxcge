# FINDING — MultiNicheNet / multi-sample NicheNet on winning pair GSE131907+GSE205335

**Status:** pipeline is on the branch. `results/ligand_activity_table.tsv` and this file are written by `scripts/03_analyze.py` after the public processed matrices are scored.

GSE207422 is not used (PR #334 NS). CLDN4 only. Patient is the unit.

```bash
python3 methods/winpair_131907_205335_multinichenet_cldn4/scripts/00_download.py
python3 methods/winpair_131907_205335_multinichenet_cldn4/scripts/01_convert_prior.py
python3 methods/winpair_131907_205335_multinichenet_cldn4/scripts/02_extract.py
python3 methods/winpair_131907_205335_multinichenet_cldn4/scripts/03_analyze.py
```
