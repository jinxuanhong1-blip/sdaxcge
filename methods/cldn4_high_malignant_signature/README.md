# CLDN4-high malignant signature

Additive public analysis. The signature is built from the locked concordant-4 malignant Q4 vs Q1 contrast (PR 503) and kept only if those genes still track CLDN4 in TCGA-LUAD after KRT8/KRT18/KRT19 and ABSOLUTE purity. Inverse immune association is tested in OncoSG LUAD and in open GEO LUAD bulks that are not the locked single-gene sets OncoSG / GSE10072 / GSE11969 / GSE248378.

```bash
python3 methods/cldn4_high_malignant_signature/analyze.py
```

Raw matrices are downloaded to `/tmp/cldn4sig` and are not committed. The write-up is `FINDING.md`. `weak_geo.py` is the size, ssGSEA, ESTIMATE, and histology grid for GSE282774 and GSE233774; `analyze.py` calls it.
