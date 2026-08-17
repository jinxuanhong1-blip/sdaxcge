# Locked patient-level extracts

Existing public tables (not re-audited):

- `GSE131907_samples.tsv` — PR #279 / #320 author-malignant T/NK extract (sample-level).
- `GSE205335_patients.tsv` — PR #279 / #320 author-malignant patients.
- `GSE207422_drmref_patients.tsv` — locked A3 DRMref 12-patient table (CLDN4 mean only).
- `GSE148071_tisch_units.tsv` — TISCH NSCLC_GSE148071 units (PR #281 / #335).
- `GSE127465_tisch_units.tsv` — TISCH NSCLC_GSE127465 patient units (PR #281).

Scored in this PR from GEO processed files <2 GB:

- `GSE123902_marker_units.tsv` — Laughney 2020 dense CSV (`GSE123902_RAW.tar`, 90 MB).
- `GSE189357_marker_units.tsv` — Zhu/Wang AIS–IAC 10x MTX (`GSE189357_RAW.tar`, 624 MB).

Re-score the GEO pair:

```bash
# tars in /tmp/geo_dl/
python3 methods/scrna_cldn4_combo_enum/scripts/score_geo_markers.py
```
