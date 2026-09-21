# Concordant-4 CNV epithelial re-test

CopyKAT and InferCNV malignant calls on the locked concordant-4
(GSE123902, GSE131907, GSE205335, GSE189357). CLDN4 % positive is then
correlated with the T/NK fraction only inside CNV-positive epithelial
cells, and compared with the annotation-only malignant call.

```bash
bash methods/concordant4_cnv_epithelial_cldn4/scripts/download.sh /tmp/concordant4_raw
python3 methods/concordant4_cnv_epithelial_cldn4/scripts/analyze.py --raw /tmp/concordant4_raw
```

GSE205335 needs R (`Matrix`) to read the double-gzipped GEO RDS.
Details: `METHODS.md`. Numbers: `FINDING.md`.
