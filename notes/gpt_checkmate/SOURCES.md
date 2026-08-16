# Source and search log

Cut-off: 2026-08-16 UTC.

## Trial identifiers searched

| Trial | Sponsor ID / registry ID |
|---|---|
| CheckMate-017 | CA209-017 / NCT01642004 |
| CheckMate-057 | CA209-057 / NCT01673867 |
| CheckMate-227 | CA209-227 / NCT02477826 |
| CheckMate-9LA | CA209-9LA / NCT03215706 |
| CheckMate-816 | CA209-816 / NCT02998528 |
| CheckMate-153 | CA209-153 / NCT02066636 |

## Authoritative data sources

- EGA public metadata API documentation: <https://ega-archive.org/discovery/metadata/public-metadata-api/>
- EGA CheckMate-153 WES study: <https://ega-archive.org/studies/EGAS00001007508>
- EGA CheckMate-153 RNA-seq study: <https://ega-archive.org/studies/EGAS00001007509>
- EGA shared controlled dataset: <https://ega-archive.org/datasets/EGAD00001011302>
- EGA machine-readable dataset metadata: <https://metadata.ega-archive.org/datasets/EGAD00001011302>
- Nature Medicine article and supplements: <https://doi.org/10.1038/s41591-024-03240-y>
- Open PMC mirror: <https://pmc.ncbi.nlm.nih.gov/articles/PMC12066197/>
- NCBI GEO DataSets E-utilities: <https://eutils.ncbi.nlm.nih.gov/entrez/eutils/>

## Search method and limits

`catalog_ega.py` downloads the complete public study and dataset collections from
the EGA metadata API, then case-insensitively matches explicit CheckMate and
CA209 aliases in titles, descriptions, and abstracts. This is stronger than a
web-search-only check but cannot see unreleased/private records or records that
omit every searched trial alias.

`search_geo.py` records exact trial-name searches in NCBI GEO DataSets. Exact
searches returned zero for all six trials at the cut-off. This does not exclude
indirectly described, misannotated, or future deposits.

Web and publication searches also used trial names, CA209 identifiers, registry
identifiers, and terms such as `RNA-seq`, `WES`, `omics`, `processed data`,
`supplementary data`, `GEO`, and `EGA`. The only verified open, individual-level
processed molecular tables were the two CheckMate-153 article supplements.

Clinical efficacy tables, protocols, PDFs, figures, and data available only
through sponsor request platforms were excluded because they are not open
processed omics. EGA FASTQ files were catalogued but not downloaded because
EGA marks the dataset controlled; file size alone does not make controlled data
eligible.
