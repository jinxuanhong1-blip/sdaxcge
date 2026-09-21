# Concordant-4 malignant NHEJ / STING / MHC-I

Observational patient-level check on the locked concordant-4 cohorts only:
GSE123902, GSE131907, GSE205335, GSE189357.

Question: among malignant cells, how CLDN4 tracks an NHEJ module, a STING/IFN
module, and MHC-I, and whether NHEJ-low still sits with IFN-high inside CLDN4
strata. Not a causal mediation. Not a re-fit of the locked T/NK result.

```bash
bash methods/nhej_sting_concordant4/scripts/download.sh /tmp/geo_nhej
python3 methods/nhej_sting_concordant4/scripts/extract_sums.py
Rscript methods/nhej_sting_concordant4/scripts/export_gse205335.R /tmp/geo_nhej /tmp/geo_nhej/gse205335_export
python3 methods/nhej_sting_concordant4/scripts/analyze_associations.py
```

Python: numpy, scipy, pandas, matplotlib. R: base + Matrix, only to read the
GSE205335 RDS. Numbers are in `FINDING.md`.
