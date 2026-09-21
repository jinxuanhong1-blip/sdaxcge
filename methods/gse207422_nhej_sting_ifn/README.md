# GSE207422 NHEJ / STING / IFN split

Public Hu et al. NSCLC series (GEO GSE207422). The run is gated on TACSTD2 and CLDN4 both being present in the matrix. They are. NHEJ, proximal cGAS–STING, and IFN effectors are scored as three non-overlapping modules.

Primary unit is the post-treatment patient. Within A3-malignant cells, CLDN4 Q4 vs Q1 (and, separately, TACSTD2 Q4 vs Q1). See `FINDING.md`.

A 64-contrast sweep hunted the opposite pattern (CLDN4-low with higher IFN, higher STING, and lower NHEJ). None of the contrasts had that joint sign. Status for this series: **FINAL discordant**.

Raw matrices stay in `data/GSE207422/` and are gitignored.

```bash
python3 methods/gse207422_nhej_sting_ifn/scripts/download.py
python3 methods/gse207422_nhej_sting_ifn/scripts/extract.py
python3 methods/gse207422_nhej_sting_ifn/scripts/analyze.py
python3 methods/gse207422_nhej_sting_ifn/scripts/analyze_bulk.py
python3 methods/gse207422_nhej_sting_ifn/scripts/sweep_thesis.py
```
