# methods/scrna_milo_cldn4

Additive Milo-style kNN neighbourhood DA vs **malignant CLDN4** on the public
GSE207422 processed UMI (Hu et al., *Genome Med* 2023). Prior TACSTD2 Milo
(`methods/scrna_milo`) is taken as given and is not re-run.

**miloR is not used.** Same documented fallback as the TACSTD2 folder: sample-level
Spearman / Welch on neighbourhood proportions, SpatialFDR with k-distance weights.

## Reproduce

```bash
pip install -r methods/scrna_milo_cldn4/requirements.txt
python3 methods/scrna_milo_cldn4/scripts/download.py --dataset GSE207422
python3 methods/scrna_milo_cldn4/scripts/run_gse207422.py
```

`FINDING.md` is written from `results/GSE207422/summary.json` at the end of the run.
The independent unit is the post-treatment sample. Do not cite cell count as *n*.
