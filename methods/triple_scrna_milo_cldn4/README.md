# Triple-merge Milo / neighbourhood DA vs malignant CLDN4

Additive **CLDN4-only** neighbourhood differential abundance on the public
processed files from **GSE131907 + GSE148071 + GSE205335**.

- No dual-high (TACSTD2 is not a gate).
- No GSE207422 (a different agent).
- Independent unit = **patient / sample**, not the cell.
- Graph is built **per dataset** (per-sample PCA centering). Joint Harmony
  of the three public matrices is not feasible at ~15 GB RAM.
- miloR / edgeR are not used. DA = sample-level Spearman + SpatialFDR
  (Dann et al. 2022 k-distance formula). See `playbook.md`.

```bash
pip install -r methods/triple_scrna_milo_cldn4/requirements.txt
python3 methods/triple_scrna_milo_cldn4/scripts/download.py
python3 methods/triple_scrna_milo_cldn4/scripts/analyze.py
```

Deliverable: `FINDING.md` and `tables/nhoods.tsv`.
