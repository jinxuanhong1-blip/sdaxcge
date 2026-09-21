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

- PMID 41016339 (Kashiwagi et al., BBRC). CRISPR CLDN4 knockout RNA-seq in NCI-H1688 small-cell lung cancer cells, with SAA1 reported as an upregulated effector. Europe PMC data links for this paper are only RRID:Addgene_98293 and RRID:Addgene_52961. A second pass on 2026-09-21 (author site, Dryad, lab GitHub, WeChat/CN supplements) also found no matrix. See below.
- PMID 40892111 (Zheng et al.). The RNA-seq in the abstract is cerulein acute pancreatitis versus control, used to nominate CLDN4. Knockdown is the later functional assay. Europe PMC lists zero data links. GSA-Human and OMIX searches for CLDN4 returned no study.

## PMID 41016339 second pass: final negative

Checked 2026-09-21. Corresponding author Korehito Kashiwagi, Dokkyo Medical University pathology (k-kore@dokkyomed.ac.jp), with Takuya Yazawa and Hideki Chiba (Fukushima). Nothing open was downloaded.

| Place | Result |
|---|---|
| Lab site https://dept.dokkyomed.ac.jp/dep-m/pathology/ | Citation only. No data, GitHub, or download link on the research, publications, or members pages. |
| researchmap https://researchmap.jp/ayakore0303/published_papers/51267947 | DOI only. The 資料公開 section lists no files. |
| Dokkyo repository https://dmu.repo.nii.ac.jp/records/2000795 | One-page conference abstract (Dokkyo Journal of Medical Sciences 52(2), K-4, 2025-12-25). It restates the H1688 knockout RNA-seq and SAA1 result. No accession and no count table. |
| KAKEN JP18K15958 final report (2022-01-27) and JP21K08164 | Claudin-4 in SCLC is described. No RNA-seq accession and no repository URL. |
| Dryad API | DOI search: 0. "H1688": 0. "Yazawa CLDN4": 0. The "CLDN4" hit is doi:10.5061/dryad.47d7wm3pc (unrelated lung epithelium). The "Kashiwagi" hit is Xenopus strain data, doi:10.5061/dryad.m6f93. |
| Zenodo DOI search and Figshare "Kashiwagi CLDN4" | 0 records. |
| GitHub | No dokkyomed, fmu-chiba, or hidchiba user. `poojascis/CLDN4_Data` is SNP tables (nsSNV and UTR), not this RNA-seq. |
| Elsevier supplements | Crossref `relation` is empty. Unpaywall: closed. `mmc1`–`mmc5` probes for PII S0006291X25014263 returned no file. |
| WeChat (Sogou) | "CLDN4 H1688 SAA1 RNA-seq 柏木": no articles. "CLDN4 小细胞肺癌 SAA1": two unrelated posts (cell-therapy news; RET-fusion investment note). |
| CNCB GSA | Kashiwagi, "H1688 CLDN4", and "SAA1 H1688": 0. |
| DDBJ | Kashiwagi: 0 studies. H1688 studies are KLF9 knockdown, FOXM1 knockdown, and drug treatment, not CLDN4 knockout. |

Final inventory for this paper: the H1688 CLDN4 knockout RNA-seq is described in the article and in a one-page Dokkyo abstract, and it is not in a public archive, author site, Dryad, lab GitHub, or Chinese/WeChat supplement.

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
