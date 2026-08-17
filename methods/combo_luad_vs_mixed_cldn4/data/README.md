# Locked public processed extracts (no new scRNA download)

CLDN4-only histology cut of three already-processed series. Matrices were
not re-downloaded.

| File | Source | Unit | Histology in file |
|---|---|---|---|
| `GSE131907_samples.tsv` | PR #279 T/NK extract (Kim *et al.* LUAD atlas) | sample | series is LUAD; no LUSC/SCLC |
| `GSE205335_patients.tsv` | PR #279 T/NK extract | patient | `cancer_subtype`: ADC / SQ / SCLC / NUT |
| `tisch_GSE148071_units.tsv` | TISCH NSCLC pool / PR GSE148071 Q4 | patient (1 sample each) | none (GEO SOFT is age/sex; paper is mixed NSCLC) |
| `tisch_GSE148071_given_spearman.tsv` | same TISCH pool | — | locked continuous ρ (taken as given) |

GSE148071 is **not** entered as LUAD. The locked extract has no LUAD/LUSC
column. Supplementary Table 1 of Wu *et al.* is a PDF figure, not a
machine table in these extracts. LUAD-only therefore uses GSE131907 plus
GSE205335 ADC only.

Do not write n=42 for GSE148071 or n=58 for GSE131907. Eligible n is
computed in `analyze.py`.
