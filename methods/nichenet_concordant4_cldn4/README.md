# NicheNet ligand activity — concordant-4 CLDN4

Patient-aware NicheNet activity from CLDN4-high vs CLDN4-low malignant senders to T/NK receivers.

Cohorts: GSE123902, GSE131907, GSE205335, GSE189357. No discordant accession is added.

Read `RESULTS.md` for the numbers. `METHODS.md` for the definitions.

```bash
bash methods/nichenet_concordant4_cldn4/scripts/download.sh /tmp/concordant4_raw /tmp/nichenet_prior
Rscript methods/nichenet_concordant4_cldn4/scripts/extract_gse205335.R \
  --raw=/tmp/concordant4_raw --here=methods/nichenet_concordant4_cldn4 \
  --out=/tmp/nichenet_work/extract
EXTRACT_WHICH=123902,189357,131907 python3 methods/nichenet_concordant4_cldn4/scripts/extract_rest.py
python3 methods/nichenet_concordant4_cldn4/scripts/analyze.py
```
