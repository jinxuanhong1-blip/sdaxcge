# Example inputs (synthetic)

These files exist so the template can be smoke-tested without any real
`results/*.csv`. **Every effect size is invented.** They are not estimates
from GSE126044, GSE135222, or any other real cohort.

| File | What it is |
| --- | --- |
| `results/synthetic_open_ici_effects.csv` | 8 open + 1 restricted (must be dropped) + 1 overlapping + 1 wrong-endpoint rows, both genes |
| `cohort_registry.csv` | Same 8 included, plus 3 synthetic EGA/dbGaP missing trials and one gene-absent panel |
| `known_open_ici_cohorts.csv` | Catalog of **real** accessions already seen in this repo, **without effect sizes** |

When sibling analyses write real per-cohort CSVs under `results/`, point
`--results-dir` at that directory. Do not copy numbers from this folder into a
results table.
