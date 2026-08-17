# Locked inputs (PR #459 pair GSE131907 + GSE189357)

Existing public tables (not re-audited):

- `GSE131907_samples.tsv` — PR #459 / #279 / #320 author-malignant sample-level CLDN4 %pos.
- `GSE189357_marker_units.tsv` — PR #459 marker-malignant patient-level CLDN4 %pos
  (gate: `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`).
- `GSE131907_malignant_counts.tsv.gz` — author-malignant UMI-sum, 21 tumor samples
  (same matrix as PR #456 / #472; `author_Cell_subtype==Malignant cells`; drop nLung/nLN).
- `GSE189357_malignant_counts.tsv.gz` — marker-malignant UMI-sum, 9 patients
  (same matrix as PR #469; rebuild with `build_gse189357_pseudobulk.py`).
- `GSE131907_malignant_meta.tsv` / `GSE189357_malignant_meta.tsv` — cells summed per unit.
- `a8_sets.json` — Hallmark IFN-α/γ, custom MHC-I/APM, KEGG/GO tight junction, KRT_EPITHELIAL.

Rebuild GSE189357 only (public GEO processed MTX, `GSE189357_RAW.tar` ~624 MB, <2 GB):

```bash
python3 methods/pair_131907_189357_malig_ifn_de_cldn4/download.py
python3 methods/pair_131907_189357_malig_ifn_de_cldn4/build_gse189357_pseudobulk.py
```

GSE131907 counts are the committed author-malignant UMI-sum from the public
processed Kim 2020 extract used in PR #456 / #472. They are not re-derived here.
