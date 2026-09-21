# Mixscape-like CLDN4 signature vs Hallmark IFN

ADDITIVE. CLDN4-only. Concordant-4 malignant cells only
(GSE123902, GSE131907, GSE205335, GSE189357). Not CRISPR.

CLDN4-low malignant cells are treated as a KD-like query. CLDN4-high
cells are the NT-like neighbor pool. The local residual follows
Seurat `CalcPerturbSig` (k nearest NT-like neighbors subtracted in
PCA space, within each patient / donor / sample). That residual is
scored on Hallmark IFNα ∪ IFNγ, with a detection-matched control and
Hallmark spermatogenesis as a negative hallmark.

`RunMixscape` is not used. There are no gRNAs.

Numbers: `FINDING.md`. Unit table: `results/tables/unit_summary.tsv`.

```
bash methods/mixscape_like_concordant4_cldn4/scripts/download.sh /tmp/geo_mixscape
python3 methods/mixscape_like_concordant4_cldn4/scripts/mixscape_like.py
python3 methods/mixscape_like_concordant4_cldn4/scripts/summarize.py
```

The unit is the patient / donor / sample. Do not quote cell counts as n.
