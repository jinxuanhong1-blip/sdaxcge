# GSE316655 — CLDN / LILRB2 myeloid IO (public MTX/TSV)

Additive public re-score of Liu et al. *Sci Immunol* 2026 (PMID 41931598).

**Species / model:** human CD45+ FACS scRNA from **SK-MEL-5** tumors (and matched PB) grown in **NSG-SGM3** mice humanized with human cord-blood CD34+ cells. Cell Ranger `GRCh38_and_mm10-2020-A`. Not patient tumors. Not CLDN18.2 gastric xenografts (those are described in the paper; this GEO deposit is SK-MEL-5).

**n:** 4 public 10x libraries, 1 labeled donor (`ND`), 1 library per tissue × treatment cell. No sample-level inference.

```
python3 download.py --outdir data/raw
python3 analyze.py --raw data/raw --outdir .
```

See `FINDING.md`.
