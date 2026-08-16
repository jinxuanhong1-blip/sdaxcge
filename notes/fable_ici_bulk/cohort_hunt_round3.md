# Extra-cohort hunt (round 3) — open lung ICI bulk with TACSTD2 **and** CLDN4

Date: 2026-08-16. NCBI eutils `esearch` on GEO DataSets (`expression profiling by high throughput sequencing`, Homo sapiens).

Queries: NSCLC/lung + anti-PD-1 / pembrolizumab / nivolumab / atezolizumab / immune checkpoint / neoadjuvant PD-1.

**39 unique GDS IDs.** After dropping already-used series, titles were inspected. None of the new hits is an open **tumor bulk RNA-seq** lung ICI cohort that measures both TACSTD2 and CLDN4 **and** deposits a response label:

| Accession | Why not used |
|-----------|----------------|
| GSE141383, GSE179994, GSE190905, GSE210997 | scRNA / glioma / other tissue |
| GSE143791 | prostate bone mets |
| GSE150972, GSE152590, GSE176021 | single-patient / CD8 PBL / T-cell programs — not tumor bulk |
| GSE179934 | n=9, not ICI-response labeled bulk |
| GSE212521 | FAP-4-1BB agonist trial, not lung ICI bulk |
| GSE217451, GSE224216 | TLS / hMENA, n=9 |
| GSE225620 | **blood**-based neoadjuvant PD-1 (not tumor) |
| GSE235048 | spectral flow, PBMC |
| GSE237818–GSE244946, GSE273409, GSE318080 | SCLC / cell-line / mouse-human mix |
| GSE260575 | PBMC-injected xenograft + pembrolizumab |
| GSE283829 | PD1–PD-L1 interaction assay, not bulk matrix with both genes |
| GSE285029 | MDSC / PD-L1 signaling, not a lung ICI response cohort |
| GSE295969 | **blood** immunomap, metastatic NSCLC anti-PD-1 |
| GSE309452/3 | EBV epithelial, not lung ICI |
| remaining | n≤18 mechanistic / not ICI-response bulk lung |

**Durvalumab GSE253564 / GSE248378:** both genes present. GEO `series_matrix` has tissue / histology / arm only. Cell Reports Medicine (2024) Table S1 has MPR for a 16-sample arm-2 subset, but no public per-sample ID↔MPR table was retrievable (PMC blocked; no GEO supplement). **Not used for response tests** — fabricating labels is disallowed.

**GSE110390** (durvalumab NSCLC/UC): 21-gene IFN-γ panel — TACSTD2/CLDN4 absent — skipped.

**Core usable for response-after-purity (both genes + label):**
GSE126044, GSE166449, GSE207422 (binary); GSE135222 (PFS; secondary DCB = PFS≥180d).
