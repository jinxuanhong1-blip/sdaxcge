# GSE207422 Palantir destinies — CLDN4-only, MPR labels

Additive Palantir + PAGA on public GSE207422 **A3-malignant-like** epithelium. CLDN4 is the readout. No dual-high gate. Root is not CLDN4-high. Primary products are **Palantir destinies** (fate probabilities) and destiny-vs-MPR tables — not DPT alone.

```bash
pip install -r methods/gse207422_palantir_cldn4_mpr/requirements.txt
python3 methods/gse207422_palantir_cldn4_mpr/scripts/download.py
python3 methods/gse207422_palantir_cldn4_mpr/scripts/extract.py
python3 methods/gse207422_palantir_cldn4_mpr/scripts/analyze.py
```

Done when `results/tables/destiny_vs_mpr.tsv` and `results/tables/patient_destiny.tsv` exist, plus `FINDING.md`.
