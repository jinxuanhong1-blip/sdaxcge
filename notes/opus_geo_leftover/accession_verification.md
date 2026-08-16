# Accession verification (Entrez `gds`, live lookup)

Source: `scripts/opus_geo_leftover/05_download.py` → `results/opus_geo_leftover/download_manifest.json`.

Query for each row: `GSE########[Accession]` on `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds`. The returned document’s `accession` field was required to equal the query. All rows below passed.

| Accession | verified | n_samples | gdstype | taxon | platform GPL | BioProject | GEO record |
|---|---|---|---|---|---|---|---|
| GSE253564 | yes | 32 | Expression profiling by high throughput sequencing | Homo sapiens | 24676 | PRJNA1066291 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564 |
| GSE248378 | yes | 29 | Expression profiling by high throughput sequencing | Homo sapiens | 24676 | PRJNA1043632 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378 |
| GSE161537 | yes | 82 | Expression profiling by high throughput sequencing | Homo sapiens | 18573 | PRJNA678660 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161537 |
| GSE162520 | yes | 92 | Expression profiling by high throughput sequencing | Homo sapiens | 18573 | PRJNA682139 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162520 |
| GSE309652 | yes | 72 | Expression profiling by array | Homo sapiens | 31904 | PRJNA1336689 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309652 |
| GSE182328 | yes | 44 | Expression profiling by high throughput sequencing | Homo sapiens | 18573 | PRJNA755863 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182328 |
| GSE248249 | yes | 42 | Expression profiling by array | Homo sapiens | 23126 | PRJNA1043107 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248249 |
| GSE216297 | yes | 286 | Expression profiling by high throughput sequencing | Homo sapiens | 20301 | PRJNA892949 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216297 |
| GSE111414 | yes | 20 | Expression profiling by high throughput sequencing | Homo sapiens | 21290 | PRJNA436957 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111414 |

Named exclusions (not downloaded in this slice): GSE126044, GSE135222, GSE136961, GSE166449, GSE93157, GSE207422, GSE205335. GSE126045 is the SuperSeries of GSE126044.

Files > 2 GB that were refused rather than fetched: `GPL23126_family.soft.gz` (4,849,895,138 bytes). No FASTQ/SRA payload was downloaded.
