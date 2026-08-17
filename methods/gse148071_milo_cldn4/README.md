# methods/gse148071_milo_cldn4

Additive Milo-style kNN neighbourhood DA vs **malignant CLDN4** on the public
GSE148071 per-sample raw-count matrices (Wu et al., *Nat Commun* 2021;
PMID 33953163). This is a different cohort from the GSE207422 TACSTD2 / CLDN4
Milo folders. Those runs are taken as given and are not re-run.

**miloR is not used.** Same documented fallback: sample-level Spearman / Welch
on neighbourhood proportions, SpatialFDR with k-distance weights.

The independent unit is the **patient biopsy** (one sample per patient in the
GEO deposit). Do not cite cell count as *n*.

## Reproduce

```bash
pip install -r methods/gse148071_milo_cldn4/requirements.txt
python3 methods/gse148071_milo_cldn4/scripts/download.py
python3 methods/gse148071_milo_cldn4/scripts/run_gse148071.py
```

`FINDING.md` is written from `results/GSE148071/summary.json` at the end of the run.
