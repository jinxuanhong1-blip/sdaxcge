# Winning-pair MultiNicheNet (CLDN4-only)

Additive multi-sample NicheNet on **GSE131907 + GSE205335** only.

- Sender = CLDN4-high vs low malignant (no dual-high).
- Receiver = same-patient T/NK.
- Patient is the unit.
- GSE207422 is not used (PR #334 NS).
- If `nichenetr` / `multinichenetr` are missing, documented NicheNet-v2 ligand–target prior scoring is used.

See `FINDING.md` and `results/ligand_activity_table.tsv`.

```bash
pip install -r methods/winpair_131907_205335_multinichenet_cldn4/env/requirements.txt
python3 methods/winpair_131907_205335_multinichenet_cldn4/scripts/00_download.py
python3 methods/winpair_131907_205335_multinichenet_cldn4/scripts/01_convert_prior.py
python3 methods/winpair_131907_205335_multinichenet_cldn4/scripts/02_extract.py
python3 methods/winpair_131907_205335_multinichenet_cldn4/scripts/03_analyze.py
```
