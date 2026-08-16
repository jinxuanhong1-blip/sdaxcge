# Verification log — GEO 2020 leftover TACSTD2 / CLDN4

Date of search: 2026-08-16. NCBI E-utilities + GEO FTP HTTPS.

## 1. Prior-wave inventory

`results/w200/GEO_2020/prior_gse_coverage.tsv` was built from a git-grep of
every remote branch of `jinxuanhong1-blip/sdaxcge` (595 unique GSE mentions;
72 with a dedicated path). Prior GEO / ICI waves that already touched 2020
lung series: `fable-geo-2019-2021`, `fable-ici-bulk`, `ici-lung-response`,
`geo-remaining-sweep`.

Already-analyzed 2020 (or 2020-adjacent) accessions **not** re-claimed here:
GSE126044, GSE136961. Later-year analyses (GSE135222, GSE166449, GSE182328,
GSE207422) are out of this date window.

## 2. Search

Script: `scripts/w200/GEO_2020/01_search_geo.py`.

Primary query: lung terms AND ICI terms AND `"Homo sapiens"[Organism]` AND
`"gse"[Entry Type]` AND `("2020/01/01"[PDAT] : "2020/12/31"[PDAT])`.

Result: **count=28, 28 UIDs fetched**. Per-drug add-ons (nivolumab,
pembrolizumab, atezolizumab, durvalumab, ipilimumab, avelumab, cemiplimab)
added **0** new UIDs.

UIDs: `results/w200/GEO_2020/search_uids.json`.
Metadata: `results/w200/GEO_2020/candidates_metadata.json` (esummary).

The 28 accessions match the PDAT-2020 rows of the 2019–2021
`candidates_metadata.json`. No extra 2020 series appeared.

## 3. Rule triage (28)

| category | n |
|---|---|
| ALREADY_ANALYZED | 2 (GSE126044, GSE136961) |
| LUNG_ICI_LEFTOVER_CANDIDATE | 5 |
| LUNG_ICI_single_cell | 4 |
| LUNG_ICI_methylation | 2 (GSE126043, GSE126045 — methylation arms of GSE126044) |
| LUNG_ICI_cellline_invitro | 3 |
| EXCLUDED_not_ICI | 1 |
| EXCLUDED_not_lung_cancer | 4 |
| EXCLUDED_not_lung | 7 |

Leftover candidates: GSE124885, GSE141479, GSE150972, GSE154286, GSE99995.

Force-probed as well (2019–2021 `LUNG_ICI_other` / controls): GSE139327,
GSE148944, GSE159785, GSE159787, GSE136961, GSE126044.

## 4. Record-level probe

FTP listings + series-matrix headers:
`results/w200/GEO_2020/tables/leftover_probe.csv` and per-accession
`downloads/<GSE>/probe.json`.

Characteristics tables:
`results/w200/GEO_2020/clinical/*_geo_characteristics.csv`.

### GSE141479 (PMID 31855576)

- Platform GPL23126 (Clariom D Human). 135,750 probes, 74 samples.
- `tissue=Blood`, `cell_type=CD8+`, `source_name=Peripheral blood, CD8+ sorted`
  for **all 74** samples. 41 individuals; 41 pre, 33 post-nivolumab.
- Characteristic scan for respon / recist / pfs / benefit / dcb: **0 hits**.
- Gene assignment field (not `Gene symbol`) maps
  `TC0100014340.hg.1` → TACSTD2 and `TC0700007993.hg.1` → CLDN4.
- Linear intensities: TACSTD2 mean 17.74 (23.8th percentile of probe-means);
  CLDN4 mean 10.53 (3.7th percentile). Array median-of-means = 31.25.
- Pre vs post Mann–Whitney p = 0.40 / 0.62; paired Wilcoxon n=33 p = 0.49 / 0.76.
- **Verdict:** leftover with genes, wrong compartment, no ICI endpoint in GEO.

### GSE154286 (PMID 32767058)

- GPL28856 custom NanoString, 201 IDs. `treatment_cohort` = P / PB;
  `timepoint` = Diagnosis / Progression. Protocol text: chemotherapy, not ICI.
- Panel IDs: `tables/GSE154286_panel_ids.txt`. TACSTD2 / TROP2 / CLDN4: absent.
- **Verdict:** not ICI; genes not measured.

### GSE150972

- n = 1. `!series_matrix_table` present but 0 expression rows.
- Suppl: RAW.tar only. **Verdict:** unusable.

### GSE99995 (PMID 29516506)

- 12 IIIA LUAD, IFN-γ high/low × PD-L1 high/low. No ICI drug, no ICI outcome.
- **Verdict:** not an ICI leftover. Not analyzed as a surrogate.

### GSE124885

- n = 4, pediatric autoimmune / post-infection. **Verdict:** not lung-cancer ICI.

### Already done (reconfirmed, not re-run)

- GSE126044: prior tumor R vs NR analysis stands.
- GSE136961: Oncomine 395-gene panel, targets absent.

## 5. What was not done (on purpose)

- No FASTQ / CEL / RAW.tar download.
- No invented response labels from the GSE141479 paper.
- No re-analysis of GSE126044 presented as a new leftover finding.
- No IFN-γ (GSE99995) dressed up as an ICI endpoint.
- GPL23126 full SOFT (~1 GB) and the 56 MB GSE141479 series matrix were used
  locally and are gitignored. The two-gene extract is committed.

## 6. Bottom line

Calendar-2020 GEO leftover work is **negative**. The parked 2020
`LUNG_ICI_other` list is closed. The only 2020 tumor ICI series that measures
TACSTD2 / CLDN4 is GSE126044, already analyzed, n = 16, non-significant.
