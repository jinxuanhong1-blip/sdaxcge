# Data provenance

Malignant patient-pseudobulk UMI sums and locked unit tables were copied from the
concordant-4 CLDN4 page (PR #503 / `methods/concordant4_123902_131907_205335_189357_cldn4/data`).

| file | role |
|---|---|
| `GSE*_malignant_counts.tsv.gz` | malignant UMI-sum per unit |
| `GSE*_malignant_meta.tsv` | matrix column metadata |
| `GSE123902_marker_units.tsv` | donor-level T/NK + CLDN4 |
| `GSE131907_samples.tsv` | sample-level T/NK + CLDN4 |
| `GSE205335_patients.tsv` | patient-level T/NK + CLDN4 |
| `GSE189357_marker_units.tsv` | patient-level T/NK + CLDN4 |
| `a8_sets.json` | Hallmark / KEGG / GO gene sets |

These files are the same inputs as PR #741. TACSTD2 expression is computed from the count matrices. The only pre-stored %pos column is CLDN4.
