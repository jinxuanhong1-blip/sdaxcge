# Inputs for the observational CLDN4-quartile proxy

These files are the concordant-4 malignant pseudobulks and unit tables already
used for the stacked Q4 vs Q1 family scores. They are not a new GEO download
and not a new malignant-cell call.

Copied from `methods/concordant4_123902_131907_205335_189357_cldn4/data`:

- `GSE123902_marker_units.tsv`, `GSE123902_malignant_counts.tsv.gz`
- `GSE131907_samples.tsv`, `GSE131907_malignant_counts.tsv.gz`
- `GSE205335_patients.tsv`, `GSE205335_malignant_counts.tsv.gz`
- `GSE189357_marker_units.tsv`, `GSE189357_malignant_counts.tsv.gz`
- `gene_sets.json` — Hallmark IFN-α, Hallmark IFN-γ, custom MHC-I/APM, KEGG tight junction, GO tight-junction organization. Subset of that folder's `a8_sets.json`.

Quartiles are recomputed here with the same rule: within-cohort malignant CLDN4 %pos, rank then `qcut` into Q1–Q4. GSE123902 and GSE189357 store %pos as a fraction and are rescaled to percent when the cohort maximum is ≤ 1.5.

P4001 is absent from the GSE205335 malignant UMI-sum. That absence is inherited from the matrix, not re-decided here.
