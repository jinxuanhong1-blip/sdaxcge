# Merged GSE131907 + GSE205335 Milo vs malignant CLDN4

Additive **CLDN4-only** neighbourhood DA on the public UMIs of Kim et al.
2020 (GSE131907) and Ahn / Lee eLife (GSE205335).

PR #320 T/NK ρ is taken as given and is not re-audited. No dual-high
TACSTD2×CLDN4 score. GSE207422 is not run (that Milo was SpatialFDR-null
at n=7).

**Not miloR.** Graph + k-distance SpatialFDR follow Dann et al. 2022; DA
is a sample/patient Spearman on neighbourhood proportions. Harmony is an
optional joint graph if both HVG matrices fit; otherwise graphs are
per-dataset.

```bash
pip install -r methods/merge_131907_205335_milo_cldn4/requirements.txt
python3 methods/merge_131907_205335_milo_cldn4/scripts/download.py
python3 methods/merge_131907_205335_milo_cldn4/scripts/run_merge.py
```

`FINDING.md` and `tables/nhoods.tsv` are written by the runner.
