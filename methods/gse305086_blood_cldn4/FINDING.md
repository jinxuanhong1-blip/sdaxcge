# GSE305086 blood CLDN4 / TACSTD2 — detection

Epithelial genes in blood are expected to be low. CLDN4 and TACSTD2 are epithelial / tumour-surface transcripts; in PAXgene whole blood they read as shed-tumour or ambient signal, not as a circulating immune program. This slice measures that on the public array and stops there.

## Dataset

| Item | Value |
|---|---|
| GEO | [GSE305086](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305086) |
| Title | Blood gene expression under immunotherapy as potential biomarker for immune checkpoint blockade in non-small-cell lung cancer |
| Platform | GPL570 Affymetrix Human Genome U133 Plus 2.0 |
| Compartment | PAXgene **whole blood** (not PBMC-sorted) |
| Arrays | **n = 173** (26,452 probes in the processed RMA matrix) |
| Timepoints | baseline 76 + cycle-4 follow-up 76 + age-matched controls 21 |
| Public regimen | 1LIO 58 arrays (29 patients × 2) / 2LIO 50 / CHTIO 44 / control 21 |
| Clinical cohort in the series text | 165 aNSCLC patients (arrays on 76/165) |

The series text describes PFS (RECIST v1.1) and long-term response in the parent study. Those labels are **not** in the GEO sample characteristics. Public fields are `sample type` (baseline / follow-up / control) and `group` (1LIO / 2LIO / CHTIO / control). See `label_inventory.tsv`.

## Platform confirmation

Both named probes are on GPL570 and map to the expected Entrez symbols (official GEO annotation, 2016-08-09):

| Probe | Expected | GPL570 symbol | Entrez | Title | In processed matrix |
|---|---|---|---|---|---|
| `201428_at` | CLDN4 | CLDN4 | 1364 | claudin 4 | yes |
| `202286_s_at` | TACSTD2 | TACSTD2 | 4070 | tumor-associated calcium signal transducer 2 | yes |

Same-gene extras on the platform: CLDN4 `1569421_at` (in matrix); TACSTD2 `202285_s_at`, `202287_s_at` (in matrix), `227128_s_at` (on GPL570, **absent** from this filtered matrix). Full table: `platform_probe_confirm.tsv`.

## Detection (n = 173)

RMA has no MAS5 present/absent calls. Detection rate = fraction of arrays in which the probe ranks above that array’s median intensity (26,452 probes). Immune lineage probes are the positive control; EPCAM is the epithelial background floor.

| Probe | Gene | Mean log2 | Array-wide percentile of mean | Median per-sample percentile | Detection rate (> sample median) | Detection rate (> sample p75) |
|---|---|---:|---:|---:|---:|---:|
| `201428_at` | **CLDN4** | 3.42 | **26.4** | 26.4 | **0 / 173 (0%)** | **0 / 173 (0%)** |
| `202286_s_at` | **TACSTD2** | 4.25 | **42.0** | 33.8 | **57 / 173 (33%)** | **3 / 173 (1.7%)** |
| `201839_s_at` | EPCAM | 3.05 | 18.1 | 15.9 | 2 / 173 (1.2%) | 0 / 173 (0%) |
| `212587_s_at` | PTPRC | 11.97 | 99.3 | 99.2 | 173 / 173 (100%) | 173 / 173 (100%) |
| `213539_at` | CD3D | 8.82 | 92.7 | 92.5 | 173 / 173 (100%) | 173 / 173 (100%) |
| `205758_at` | CD8A | 8.29 | 89.8 | 89.9 | 173 / 173 (100%) | 171 / 173 (99%) |

Named-probe CLDN4 sits with EPCAM in the left tail. Named-probe TACSTD2 is below the array median, with a thin right tail (3 arrays above p75). That is the expected shed-tumour / ambient pattern, not a blood-expressed target.

The other CLDN4 probe `1569421_at` ranks at the 71st percentile (detection > median 173/173; > p75 15/173). The two CLDN4 probes do not agree. This slice uses the named probe `201428_at`.

Full numbers: `detection_table.tsv`.

## Labels that would support an IO-benefit test

| Label | Public on GEO? | Honest n if used |
|---|---|---|
| IO benefit / LTR / RECIST | no | — |
| PFS time / event | no | — |
| Cycle-4 timepoint | yes (paired follow-up) | 76 paired patients; **1LIO 29 paired** |

Because the named probes are at / near array background, this slice does not build an IO-benefit or PFS model from regimen or from the cycle-4 pair. Cycle-4 deltas for the named probes are in `cycle4_named_probes.tsv` (1LIO n=29: CLDN4 median Δ −0.03 log2, Wilcoxon p=0.36; TACSTD2 median Δ +0.10 log2, p=0.39). Those shifts sit inside the background range above.

## What this adds

A platform-confirmed, n-honest detection table for blood CLDN4 / TACSTD2 on the largest public NSCLC 1L-IO Affy whole-blood series. The measurement is the low blood signal.

## Reproduce

```bash
# downloads (not committed)
curl -fL -o /tmp/geo/GSE305086_Expression_matrix_final.csv.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305086/suppl/GSE305086_Expression_matrix_final.csv.gz
curl -fL -o /tmp/geo/GSE305086_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305086/matrix/GSE305086_series_matrix.txt.gz
curl -fL -o /tmp/geo/GPL570.annot.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz
python3 methods/gse305086_blood_cldn4/analyze.py
```

## Outputs

- `detection_table.tsv` — named probes + lineage controls
- `platform_probe_confirm.tsv` — GPL570 symbol check
- `label_inventory.tsv` — public vs missing outcome fields
- `sample_annotation.tsv` — 173 arrays
- `cycle4_named_probes.tsv` — paired n only; not an IO-benefit model
- `summary.json`
