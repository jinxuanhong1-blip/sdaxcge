# methods/quad_207422_scenic_cldn4

Additive **CLDN4-only** SCENIC/AUCell on the QUAD malignant merge:

`GSE207422 + GSE131907 + GSE148071 + GSE205335`

GSE207422 is **included**. A10 ELF3–CLDN4 is **taken as given** and is not re-tested as a success criterion. There is **no TACSTD2 dual-high** split.

Question: at the **patient** level, do CLDN4-high vs CLDN4-low malignant cells differ in IFN / MHC / TJ regulon AUCell?

AUCell uses public TF–target priors (TRRUST + DoRothEA + CollecTRI) plus hardcoded ISG / MHC-I APM / TJ program sets. This is **not** cisTarget and **not** ChIP.

## Run

```bash
python3 methods/quad_207422_scenic_cldn4/scripts/download_priors.py
python3 methods/quad_207422_scenic_cldn4/scripts/download_data.py --outdir /tmp/quad_scenic
python3 methods/quad_207422_scenic_cldn4/scripts/analyze.py
```

Matrices stay under `/tmp/quad_scenic/` and are not committed.

## Outputs

- [`FINDING.md`](FINDING.md) — verdict, honest n, primary table
- [`tables/regulons.tsv`](tables/regulons.tsv) — regulon gene lists
- [`tables/patient_scores.tsv`](tables/patient_scores.tsv)
- [`tables/tests.tsv`](tables/tests.tsv)
