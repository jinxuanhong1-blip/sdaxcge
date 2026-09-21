# Concordant-4 surface genes vs patient T/NK

Head-to-head of malignant **CLDN1, CLDN4, CLDN7, and EPCAM** against the same-patient
T/NK fraction on the locked concordant-4 unit set (n = 65). Keratin
(KRT8 / KRT18 / KRT19) is partialled on both sides. MUC1 is extracted as
context because the locked bulk surface ranking named it; it is not part of
the pin rule.

Cohorts, and only these: GSE123902, GSE131907, GSE205335, GSE189357.
The unit, the malignant gate, and `frac_tnk` match PR #539. This script does
not refit Harmony.

```bash
python3 methods/concordant4_surface_tnk_keratin/analyze.py
```

GEO matrices are not committed. Set `GEO_DIR` (default `/tmp/geo_c4`) to a
folder that contains:

- `GSE123902_RAW.tar`
- `GSE189357_RAW.tar`
- `GSE131907_Lung_Cancer_cell_annotation.txt.gz`
- `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`
- `GSE205335_Lung_IO_CellIdentity.txt.gz`
- `GSE205335_family.soft.gz`
- `GSE205335_Lung_IO_UMI_matrix.rds.gz`

`extract_gse205335.R` double-gunzips the RDS and writes the eight gene rows.
Python does the patient scores, the DerSimonian–Laird partial Spearmans, and
the within-cohort permutation test of CLDN4 versus CLDN7.

The run refuses to write `FINDING.md` if recomputed CLDN4 % positive does not
match `data/pr539_cldn4_units.tsv`.
