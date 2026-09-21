# Concordant-4 malignant hdWGCNA-like modules (CLDN4 strata)

Additive analysis on the locked concordant-4 cohort only:

- GSE123902 (donor)
- GSE131907 (sample)
- GSE205335 (patient)
- GSE189357 (patient)

Not GSE148071, GSE127465, GSE207422, GSE154826, or E-MTAB-13526.

Malignant cells only. Metacells are built inside patient × CLDN4 status (UMI > 0 vs UMI = 0) and are never mixed across patients. Two signed co-expression networks are built, one per stratum. TJ/barrier vs IFN is assigned by enrichment against the locked gene sets, before any T/NK test. The immune test is a patient-level Spearman of the malignant module score vs the locked T/NK fraction, meta-analyzed across the four cohorts.

This does not replace the locked CLDN4 % positive vs T/NK result.

```bash
bash methods/hdwgcna_concordant4_cldn4/scripts/download.sh /tmp/geo_hdwgcna
Rscript methods/hdwgcna_concordant4_cldn4/scripts/export_gse205335.R
python3 methods/hdwgcna_concordant4_cldn4/scripts/extract_counts.py
python3 methods/hdwgcna_concordant4_cldn4/scripts/run_hdwgcna.py
```

Numbers are written to `FINDING.md` by the last script.
