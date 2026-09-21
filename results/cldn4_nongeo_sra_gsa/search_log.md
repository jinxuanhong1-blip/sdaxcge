# Search log — CLDN4 shRNA/CRISPR RNA-seq not in GEO

Date: 2026-09-21. Live NCBI E-utilities, Europe PMC, GSA, GSA-Human, and OMIX. No accession in the table was invented.

## Eligibility

A hit is a CLDN4/Cldn4 loss-of-function transcriptome (shRNA, siRNA, CRISPR, or knockout) deposited in SRA/BioProject or CNCB GSA / GSA-Human / OMIX, with no GEO series.

Counted as not downloadable when the paper describes the experiment but names no accession, or when the archive says author-on-request.

Not counted:

- The same experiment already in GEO (GSE207704, GSE50927, GSE22493).
- Observational RNA-seq in which CLDN4 is a differentially expressed gene rather than the perturbed gene.
- CRISPR screens, 4C-seq, ATAC-seq, or knockdown of a different gene in a CLDN4-expressing line.
- TACSTD2/TROP2 perturbations.

## NCBI counts

| Database | Query | Count |
|---|---|---|
| bioproject | CLDN4 / Cldn4 / "claudin 4" / "claudin-4" | 31 |
| sra | same tokens | 45 |
| sra | CLDN4/Cldn4 plus shRNA/knockdown/CRISPR/knockout/siRNA/silencing | 3 (all 4C-seq titles that contain both CLDN4 and MSH2KO) |
| sra | "CLDN4 shRNA", "CLDN4 knockout", "CLDN4 knockdown", shCLDN4, sgCLDN4 | 0 each |
| bioproject | "CLDN4 knockout" | 1 (PRJNA856719) |
| bioproject | "CLDN4 KO" or "Cldn4 KO" | 1 (PRJNA219385) |
| gds | CLDN4/Cldn4/"claudin 4" and gse[Entry Type] | 25 |
| pubmed | CLDN4/Cldn4 knockdown/knockout/CRISPR/shRNA/siRNA phrases, unioned with the RNA-seq abstract query | 87 |

`CLDN4-/-` was not used as a query. In Entrez the hyphen is a NOT operator, and that search returned unrelated records.

SRA runs returned by the CLDN4 token, grouped by BioProject, were all public (`cluster_name=public`):

| BioProject | Runs indexed | Public bytes | GEO |
|---|---|---|---|
| PRJNA219385 | 5 | 34,419,008,894 | GSE50927 |
| PRJNA856719 | 4 (CLDN4-/- titles only) | 1,812,472,134 | GSE207704 |
| PRJNA634686 | 6 | 4,771,025,110 | none (4C-seq) |
| PRJNA273003 | 8 | 25,032,941,413 | GSE65107, lymphocytic colitis, not a knockout |
| PRJNA428887 | 16 | 944,700,907 | other assay, GSM titles, not a CLDN4 knockout |
| PRJDB17476 | 6 | 31,256,234,691 | vasopressin-receptor intestine RNA-seq; CLDN4 is a text hit, not the perturbation |

Separate check: `PRJNA1424384[BioProject]` in SRA = 6 public sgRNA-seq runs, 229,032,623 bytes, GEO count 0.

## Papers

Europe PMC core records were pulled for 38 PMIDs that are direct CLDN4 perturbations or that came up in the RNA-seq abstract query. Full text XML was fetched when a PMCID existed, and accessions matching `GSE`, `PRJNA`, `PRJEB`, `PRJCA`, `PRJDB`, `CRA`, `HRA`, `OMIX`, `E-MTAB` were extracted. The only new accession in that set was PRJNA1424384 (BFT CRISPR screen). GSE207704 and GSE50927 reappeared as the known GEO series.

Publisher page for PMID 40892111 (Springer) states data are available from the corresponding author on reasonable request. The open ESM docx contains primer sequences only.

## CNCB

- `https://ngdc.cncb.ac.cn/gsa/search?searchTerm=CLDN4` — Total Items 11, all SRX/INSDC. Same 11 for `Cldn4`. No CRA.
- `https://ngdc.cncb.ac.cn/gsa-human/json/finished.json` — 7,312 records. Zero strict CLDN4/claudin-4/Cldn4 matches.
- `https://ngdc.cncb.ac.cn/omix/releaseList` — complete HTML, 10,741 distinct `OMIX` accessions, zero CLDN4/claudin-4 strings.
- BioStudies API `CLDN4 AND (knockout OR knockdown OR shRNA OR CRISPR OR siRNA)` — 84 hits, literature plus E-GEOD-22493 and E-GEOD-50927. No new expression accession.

## Download decision

No record met both conditions (CLDN4/Cldn4 shRNA or CRISPR RNA-seq, and open outside GEO). Public FASTQ of the GEO series and of the 4C-seq / sgRNA-seq near-misses was left in the archive.
