# HPA lung TACSTD2 / CLDN4 IHC and protein

Official Human Protein Atlas public tables only. Retrieved 2026-08-16 from [proteinatlas.org](https://www.proteinatlas.org) v25.1 (IHC knowledge tables, gene TSV/JSON/XML, MS, CPTAC slice) and [v23.proteinatlas.org](https://v23.proteinatlas.org) (legacy `normal_tissue.tsv` / `pathology.tsv` for a same/different check).

Nothing here is an H-score, a paired TACSTD2+CLDN4 core table, or an ICI outcome table. Those columns do not exist in HPA. Gaps are listed in `tables/13_not_available.tsv`.

## Start here

| File | What it is |
| --- | --- |
| `tables/02_normal_ihc_knowledge_respiratory.tsv` | Official cell-type IHC for lung / bronchus / nasopharynx |
| `tables/03_cancer_ihc_knowledge_lung.tsv` | Official lung-cancer TMA ordinal counts |
| `tables/07_antibody_cancer_ihc_patients.tsv` | Official per-patient ordinal scores + SNOMED + image URLs |
| `tables/14_cancer_ihc_counts_by_antibody.tsv` | Recount of those patients vs the knowledge table |
| `tables/13_not_available.tsv` | What HPA does not publish |

Re-run: `python3 scripts/w200_hpa_lung/download_and_extract.py`

License: HPA copyrightable data are [CC BY 4.0](https://www.proteinatlas.org/about/licence). Cite Uhlén et al., *Science* 2015, and Human Protein Atlas version 25.1.
