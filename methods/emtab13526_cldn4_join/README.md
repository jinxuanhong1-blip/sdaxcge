# E-MTAB-13526 CLDN4 join test

ADDITIVE. **CLDN4-only.** Public De Zuani / Cvejic NSCLC atlas
([E-MTAB-13526](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13526);
PMID 38782901). Treatment-naive resections. Not ICI.

The 45.51 GB / 58.74 GB author-annotated h5ads are **not downloaded**.
Usable public processed files are the CD235a− tumor Cell Ranger matrices
(15 lanes, ~2.38 GB total, each < 350 MB). File decisions:
`tables/public_file_catalog.tsv`.

```bash
python3 methods/emtab13526_cldn4_join/download.py
python3 methods/emtab13526_cldn4_join/extract.py
python3 methods/emtab13526_cldn4_join/analyze.py
```

Join rule vs concordant pool (GSE123902+GSE131907+GSE205335+GSE189357):
sign must match (T/NK inverse and/or malignant IFN/MHC down in CLDN4-high).
No dual-high. Honest n. Decision: `FINDING.md`.
