# FINDING — GSE205335 neighbourhood DA vs malignant CLDN4

Placeholder. Filled after `scripts/analyze.py` finishes. Until then there are no SpatialFDR numbers or honest-*n* counts to cite.

Primary question: are transcriptional neighbourhoods differentially abundant with patient-level malignant CLDN4, and are CLDN4-high neighbourhoods T/NK-poor?

Rules locked before the run:

- Unit is the **patient**, not the cell and not the neighbourhood.
- Author `lineage.sub == Malignant cells` and `lineage.total == T/NK cells`.
- MPR/NMPR is not labelled on GEO; RECIST is not substituted for MPR.
- miloR / edgeR are not used. DA is patient-level Spearman / Welch + SpatialFDR (k-distance).
