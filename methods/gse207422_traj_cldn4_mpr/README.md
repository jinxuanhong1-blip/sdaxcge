# GSE207422 — CLDN4-only PAGA / DPT vs MPR

Additive public slice. Patient is the unit. No dual-high. T/NK is not re-audited.

```bash
pip install -r methods/gse207422_traj_cldn4_mpr/requirements.txt
python3 methods/gse207422_traj_cldn4_mpr/scripts/download.py
python3 methods/gse207422_traj_cldn4_mpr/scripts/extract_epithelium.py
python3 methods/gse207422_traj_cldn4_mpr/scripts/analyze_paga_dpt.py
```

Primary tables: `tables/patient_cldn4_vs_dpt.tsv`, `tables/patient_dpt_vs_mpr.tsv`.
Write-up: `FINDING.md`.
