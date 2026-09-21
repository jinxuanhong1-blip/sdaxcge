# CLDN4 knockdown RNA-seq outside GEO

Search date: 2026-09-21.

Question: is there a public CLDN4/Cldn4 knockdown or knockout RNA-seq study in BioProject, SRA, ENA, or CNCB/GSA-Human that is not already a GEO GSE?

**No.** Nothing new was downloaded.

## What would have been downloaded

A hit had to be all of the following:

- CLDN4 or Cldn4 itself knocked down or knocked out (siRNA, shRNA, CRISPR, or germline KO)
- RNA-seq (a gene-expression transcriptome), not a microarray, 4C, ChIP, or methylation array
- deposited in BioProject/SRA/ENA or CNCB GSA / GSA-Human / OMIX
- not already mirrored as a GEO GSE
- open, so the reads or a gene-level matrix could be fetched

## True CLDN4-loss RNA-seq that is already a GSE

Both public CLDN4-loss RNA-seq studies are open and already mirrored. FASTQs were left in place.

| Study | Runs | What it is | GEO mirror | Openness |
|---|---|---|---|---|
| PRJNA856719 / SRP385386 | SRR20029118–SRR20029125 (8 RNA-seq runs) | T47D and MCF7 CLDN4 knockout vs wild type | GSE207704 | Open SRA FASTQ (~0.48–0.53 GB each) and GEO processed TXT |
| PRJNA219385 / SRP029974 | SRR988118–SRR988127 | Mouse lung Cldn4 knockout ± ventilator injury | GSE50927 | Open SRA FASTQ (about 1.8–4.8 GB per file) and GEO processed CSV |

The ovarian CLDN4 siRNA series PRJNA128653 is GSE22493. It is an expression microarray, not RNA-seq.

## Not deposited

Two 2025 papers still have no BioProject, SRA, ENA, DDBJ/GEA, or GSA accession for a CLDN4-loss transcriptome.

- PMID 41016339 (Kashiwagi et al., BBRC). CRISPR CLDN4 knockout RNA-seq in NCI-H1688 small-cell lung cancer cells, with SAA1 reported as an upregulated effector. Europe PMC data links for this paper are only RRID:Addgene_98293 and RRID:Addgene_52961. No sequence study was found under the author, Dokkyo, H1688 plus CLDN4, or SAA1 plus SCLC.
- PMID 40892111 (Zheng et al.). The RNA-seq in the abstract is cerulein acute pancreatitis versus control, used to nominate CLDN4. Knockdown is the later functional assay. Europe PMC lists zero data links. GSA-Human and OMIX searches for CLDN4 returned no study.

## Open records that are not CLDN4 knockdown RNA-seq

These are public and not a GSE, so they were checked and not downloaded.

| Accession | Openness | Why it is not the requested dataset |
|---|---|---|
| PRJDB17476 / DRR528490–DRR528495 | Open paired FASTQ on ENA/DDBJ (about 2.7–3.9 GB per file). No GSE link. No GEA matrix. | Small-intestine epithelial RNA-seq of wild-type versus vasopressin-receptor V1a/V1b double-deficient mice. Claudin-4 was the FACS marker used to collect the cells, not the gene that was knocked out. |
| PRJNA634686 4C runs at the CLDN4 locus (SRR11836422, SRR11836423, SRR11836427–SRR11836430) | Open FASTQ. The BioProject has no GSE link. | 4C-seq at the CLDN4 locus in MSH2-perturbed gastric lines. The RNA-seq in the same project is MSH2 siRNA or knockout, not CLDN4. |
| HRA003416 / PRJCA013084 | Controlled. Request goes to HDAC001946 (Peking University First Hospital). | GSA-Human RNA-seq of a claudin-low clear-cell renal carcinoma subtype (20 tumor/normal pairs). Not a CLDN4 knockdown. |

GSA text search for CLDN4 returned 11 INSDC experiments, all from GSE50927, GSE207704, or the 4C-seq above. GSA-Human and OMIX returned zero hits for CLDN4, Cldn4, claudin-4, and CLDN4 knockout. DDBJ GEA has no CLDN4 knockdown submission.

## Queries

Live counts and the accept/reject table are in `accessions.tsv`, `near_misses.tsv`, and `query_counts.tsv`.
