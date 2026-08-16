# Templates

| Script | Role |
| --- | --- |
| `01_scirpy_clonotypes.py` | 10x contig CSVs → `clone_id`, within-sample expansion |
| `02_immunarch_clonality.R` | same folder → clonality / diversity tables |
| `03_patient_endpoints.py` | GSE243013 author TCR+meta (or scirpy `cells.tsv`) → E1 table; E2 skip unless a tumor-score file is passed |
| `config.template.yml` | freeze estimands before EDA |
| `metadata.template.tsv` | sample-level join key |

These scripts write tables only. They do not download controlled FASTQ and they
do not run the E1/E2 hypothesis tests.
