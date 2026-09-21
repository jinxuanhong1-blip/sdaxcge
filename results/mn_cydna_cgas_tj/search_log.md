# Search log

NCBI E-utilities, 21 Sep 2026. GDS identifiers convert to GEO series as GSE plus the digits after the `200` prefix.

## Queries that define the set

| Database | Query | Count |
| --- | --- | --- |
| gds | `(micronuclei OR micronucleus OR cGAS OR CGAS OR "cytosolic DNA" OR "cytoplasmic DNA" OR MB21D1) AND (claudin OR CLDN OR "tight junction" OR TJP1 OR occludin) AND gse[ETYP]` | 2 (GSE295033 false positive; GSE279172 real) |
| gds | `("tight junction" OR claudin OR CLDN2 OR CLDN4) AND (cGAS OR STING OR "mitochondrial DNA" OR micronuclei) AND gse[ETYP]` | 1 (GSE279172) |
| pubmed | `claudin-4[Title] AND (micronuclei OR genome OR cGAS OR autophagy)` | 25 titles, of which the micronucleus papers have no GEO link |
| pubmed | paper-title elink to gds | Mackenzie 28738408 → GSE100771; Dou 28976970 → GSE99028; Glück 28759028 → GSE100102; Bakhoum 29342134 → GSE98183; 53BP1 39107288 → GSE237615; cyDNA CRC 41722050 → GSE313807; EAC CIN 41811963 → GSE316062, GSE316127, GSE315942; IFNL tight junction 40890487 → GSE279172; IFNL1 follow-up 41525418 → GSE296527; CLDN4 micronucleus papers 38867360, 39625235, 38293054 → no gds link |

Broad gds counts, used only to find the landmark series, not as the analysis set: micronuclei 115, cGAS 609, "cytosolic DNA" 113, chromosomal instability AND (cGAS OR STING OR micronuclei) 38.

## What was downloaded

Count matrices from the GEO supplementary HTTPS site into `/tmp/geo_mn` (not committed). Sample order for GSE279172 and GSE313807 was taken from the series-matrix `!Sample_description` lines. GSE296527 high-density WT, IFNL2/3 KO, and IFNLR KO columns are the same integers as GSE279172 samples 1–3, 7–9, and 13–15.
