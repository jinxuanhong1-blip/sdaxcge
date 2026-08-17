# Leftover TISCH NSCLC merge — CLDN4-only combo + high-end

ADDITIVE slice. **CLDN4-only.** No dual-high with TACSTD2.

Hunt leftover public TISCH2 NSCLC objects that are **not** in the
GSE131907 / GSE148071 / GSE205335 / GSE207422 set:

`GSE117570`, `GSE139555`, `GSE146100`, `GSE149655`, `GSE143423`, `GSE127471`.

Inputs: public TISCH2 `expression.h5` + `CellMetainfo_table.tsv` only.
Skip files >2 GB. Skip a series when both the CLDN4 compartment
(epithelial/malignant) and T/NK are absent.

1. Patient-level malignant (else epithelial-like) CLDN4 vs T/NK per series.
2. The same on leftover merges. Spearman / Q4 vs Q1 when n allows (n≥8 for a merge claim).
3. CellChat-style focused LR only if ρ<0 with p<0.05 or a clear Q4 vs Q1 drop.

```bash
python3 -m pip install -r methods/tisch_leftover_merge_cldn4/requirements.txt
python3 methods/tisch_leftover_merge_cldn4/download.py
python3 methods/tisch_leftover_merge_cldn4/analyze.py
```

Writeup: [`FINDING.md`](FINDING.md). Combo table:
[`tables/highlighted_combos.tsv`](tables/highlighted_combos.tsv)
(empty hunt is allowed). Extra figures: `figures/`.
