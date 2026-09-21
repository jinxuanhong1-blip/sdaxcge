# Search log — NHEJ-STING logic wave

Date of the live queries: 2026-09-21.

Question: any public CLDN4/Cldn4 knockdown or overexpression RNA-seq, microarray, or RPPA that is not already in the catalog (GSE50927, GSE207704, GSE22493, and the PR #116 reject list), with priority on Bitler-lab HGSOC / PARPi experiments. GSE252340 is scored separately as a STING positive-control gene set.

## GEO DataSets (E-utilities, db=gds)

| Query | Count |
|---|---|
| `CLDN4[All Fields] AND gse[Entry Type]` | 21 |
| `Cldn4[All Fields] AND gse[Entry Type]` | 21 (same series) |
| `"claudin-4"[All Fields] AND gse[Entry Type]` | 10 |
| `"claudin 4"[All Fields] AND gse[Entry Type]` | 10 |
| CLDN4/Cldn4 plus knockdown/knockout/siRNA/shRNA/CRISPR/silencing/overexpression | 4 series: GSE22493, GSE26055 (CLDN7), GSE48443 (CLDN18), GSE50927 |
| CLDN4 plus overexpression / forced / ectopic expression | 1 series: GSE22493 |
| CLDN4 or Cldn4 plus RPPA / protein array / proteome / mass spectrometry | 0 |
| CLDN4 plus RPPA or "reverse phase" | 0 |

The GSE-level union is the same set already screened in PR #116 on 2026-08-16. No new CLDN4 series has been indexed since that hunt. GSE252340 is not in that union; it was fetched by accession because the series summary names TREX1 and STING.

## SRA

`CLDN4 AND (knockdown OR knockout OR shRNA OR siRNA OR CRISPR OR overexpression)` returned 3 runs (SRR11836422, SRR11836427, SRR11836428). All three are 4C-seq at the CLDN4 locus in MSH2-knockout gastric lines (study SRP263109). They are not CLDN4 perturbation transcriptomes.

## PubMed

Bitler-author query `Bitler[Author] AND (CLDN4 OR claudin-4)`: 9 PMIDs (41214101, 39625235, 39282307, 38867360, 38293054, 36237976, 35373300, 30606772, 28455726). None deposits a new expression or RPPA accession.

`CLDN4[Title/Abstract] AND (knockdown OR knockout OR siRNA OR shRNA OR CRISPR) AND (RNA-seq OR transcriptome OR microarray OR RPPA)`: 21 PMIDs. The only RNA-seq claims without an accession remain PMID 41016339 (H1688 CRISPR KO) and PMID 40892111 (pancreatitis shRNA), both already in the PR #116 unavailable list. PMID 37059993 is the GSE207704 paper.

Overexpression plus transcriptome/RPPA: 13 PMIDs. The only perturbation matrix is GSE22493. PMID 38477100 overexpresses Cldn4 by lentivirus and reads it out by qPCR, not by RNA-seq.

## Paper data-availability statements read in full

- PMID 35373300 (PMC8988515): "Data are available upon request from the corresponding author."
- PMID 41214101: "All raw data is available upon request of the corresponding author."
- PMID 39625235 (PMC11705808): TCGA (dbGaP PHS000178); "All other data will made available by the corresponding author upon request."
- PMID 38867360 (PMC11218812): "Data generated in this study are included in this manuscript and in its Supplementary Material. Data mining and flow cytometry raw data are available upon request."

No GEO, SRA, ArrayExpress, PRIDE, or figshare expression/RPPA matrix is named in those statements.
