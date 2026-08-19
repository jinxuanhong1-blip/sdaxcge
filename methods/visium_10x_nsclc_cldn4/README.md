# Public Visium CLDN4-only (10x NSCLC demos + open GEO LUAD)

Additive public Visium. CLDN4 only. No private 8-KL. No claim language.

## Data

1. Official 10x Genomics Visium CytAssist FFPE lung cancer demos (gene expression + coordinates):
   - [Human Lung Cancer (FFPE)](https://www.10xgenomics.com/datasets/human-lung-cancer-ffpe-2-standard) — LUSC, 3,858 spots
   - [Human Lung Cancer, 11 mm Capture Area (FFPE)](https://www.10xgenomics.com/datasets/human-lung-cancer-11-mm-capture-area-ffpe-2-standard) — neuroendocrine carcinoma, 6,195 spots
2. Open GEO Visium LUAD finished after that:
   - GSE189487 (6 sections, AIS/MIA/IAC)
   - GSE273378 and GSE300676 if the processed matrices + coordinates extracted successfully

Skipped: JGAS/HUM0394/EGA; Visium HD 8.9 GB bins; GSE307534 / GSE277206 (dedicated hunts); mouse LLC GSE303162.

## Run

```bash
python3 methods/visium_10x_nsclc_cldn4/scripts/analyze_cldn4_visium.py
```

See `RESULTS.md` (also copied to repo root).
