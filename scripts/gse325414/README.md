# GSE325414 CLDN4 vs T/NK

Public BD Rhapsody WTA from [GSE325414](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325414) (Jimenez et al., pulsed-electric-field treat-and-resect NSCLC).

Download into `data/GSE325414/` (gitignored):

- `GSE325414_matrix.mtx.gz`
- `GSE325414_features.tsv.gz`
- `GSE325414_barcodes.tsv.gz`
- `GSE325414_metadata_individual_cells.txt.gz`

from `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325414/suppl/`.

```bash
python3 scripts/gse325414/analyze_cldn4_tnk.py
```

Donor clinical fields are in `donor_clinical.tsv` (taken from GEO sample characteristics; constant within donor). The patient is the unit. Primary compartment is RESRU (resection, untreated tumor).
