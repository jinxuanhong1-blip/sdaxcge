# Maximum |ρ| for the locked CLDN4-high 221-gene signature

The gene list is the PR 590 signature: concordant-4 malignant CLDN4-high versus CLDN4-low, then genes that still track CLDN4 in TCGA-LUAD after KRT8, KRT18, KRT19, and ABSOLUTE purity. This folder does not rebuild that list.

The sweep scores OncoSG LUAD and the same GEO LUAD bulks (GSE273377, GSE282774, GSE233774) by z-mean or ssGSEA, by prefix length, and by purity adjustment. The reported maximum is the largest absolute cross-study Spearman versus CD8A inside that grid. ImmuneScore is a second endpoint. GSE10072, GSE11969, and GSE248378 are not opened.

```bash
python3 methods/cldn4_sig_max_rho/max_effect.py
```

Matrices download to `/tmp/cldn4sig` and are not committed. The write-up is `FINDING.md`.
