# C3_PDL1

Honest verdict: **C3 is not supported.**

PD-L1/CD274 does not have a public literature or omics case for “rises after CLDN4 loss.” TROP2 ADC has one 2-day CX-1 FPKM bump (GSE312098, FPKM still <1) and a flat 29-day CRC PDX (GSE311016). That is not a general rise.

Full write-up: `WRITEUP.md`.

| Path | What |
|---|---|
| `catalog/records.tsv` | 2,235 Europe PMC records (union of 25 queries) |
| `catalog/screened.tsv` | Title/abstract buckets |
| `catalog/human_reviewed.tsv` | Papers actually read for the claim |
| `omics/inventory.tsv` | Every accession considered |
| `omics/cd274_contrasts/cd274_only.tsv` | CD274 numbers |
| `omics/cd274_contrasts/key_stats.json` | Machine-readable verdict |
| `omics/tcga_cbioportal/spearman.tsv` | TCGA association (not C3) |
| `omics/search_log.md` | Live GEO / AE / OmicsDI / SRA log |
