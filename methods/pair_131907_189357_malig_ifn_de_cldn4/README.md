# Inputs

Locked **CLDN4-only** units from PR #459 (GSE131907+GSE189357 %pos n=30,
ρ=−0.542, Q4 r=−0.619). That T/NK rho is **not re-audited**.
No dual-high. No GSE148071. No CellChat.

- `GSE131907_samples.tsv` — Kim 2020 sample-level author-malignant CLDN4 %pos (PR #459).
- `GSE189357_marker_units.tsv` — Zhu/Wang AIS–IAC patient-level marker-malignant CLDN4 %pos (PR #459).
- `GSE131907_malignant_counts.tsv.gz` — author malignant UMI-sum, 21 tumor samples (same matrix as PR #456 / #472).
- `GSE189357_malignant_counts.tsv.gz` — marker-malignant UMI-sum, 9 patients (same matrix as PR #469).
- `a8_sets.json` — Hallmark IFN, custom MHC-I/APM, KEGG/GO tight junction, KRT_EPITHELIAL.

Patient is the unit (GSE131907: sample; GSE189357: patient).

```bash
python3 methods/pair_131907_189357_malig_ifn_de_cldn4/analyze.py
```
