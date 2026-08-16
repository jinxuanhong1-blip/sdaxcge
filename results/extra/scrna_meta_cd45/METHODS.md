# Methods: CD45+ / immune-library TACSTD2/CLDN4 detection vs MPR

Additive extra. Not a re-analysis of User A3 (GSE207422 malignant epithelium).

## Question

In public **CD45+ or immune-only** lung ICI / neoadjuvant scRNA, is the
patient-level detection fraction of TACSTD2 or CLDN4 lower in MPR (or RECIST
response) than in non-responders?

This is an **immune-library leak** score. It is not cancer-cell TACSTD2/CLDN4.

## Inclusion

Must-try series were opened on GEO (public processed files only; no FASTQ):

| Series | Library | ICI / neoadjuvant | MPR or RECIST on a public object | Used in meta |
| --- | --- | --- | --- | --- |
| GSE243013 | Deposited CD45+ immune MTX (1,254,749 cells) | Neoadjuvant anti-PD-1 ± chemo | Yes (GEO `pathological_response`, `radiological_response`) | Yes |
| GSE229353 | CD45-bead 10x tumor libraries (7 patients) | NAC or NAPC | MPR in public paper Table S1, not on GEO | Yes, after dropping truncated P03 MTX |
| GSE154826 | CD45+ / CITE-seq early NSCLC | No ICI on these lesions | No | No |
| 2024–2026 leftovers (GSE280232, GSE303680) | T-cell or cytokine studies | Not CD45-all leak libraries | Not used | No |
| GSE207422 / GSE241934 / GSE291670 | Unsorted TME | Yes | Yes | No (not immune-only; A3 taken as given) |

Hui et al. 2022 (*Cell Death Dis*) is a 12-patient CD45+ neoadjuvant set from the
same group as GSE229353. The article has no GEO accession.

## GSE243013

- Metadata: `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz` (independent GEO
  download; 1,254,749 cells, 243 `sampleID`s).
- Genes: `GSE243013_genes.csv.gz`. TACSTD2 and CLDN4 each appear once.
- Expression: nonzero MTX entries for those two genes from the public
  `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz` (7.12 GB; 2,010,550,708 nnz).
  The two-gene extract is the public-matrix scan (22,966 nonzero rows).
  Patient-level fractions were recomputed here by joining that extract to the
  independently downloaded metadata.
- One sample (`P433`, label `unknowm`) was excluded. Analyzed **n=242**
  (pCR 85 + MPR 45 = MPR-any 130; non-MPR 112). This matches the paper’s
  MPR definition (RVT ≤10%; pCR ⊂ MPR). The publication reports 234 patients;
  we use deposited labels, not a forced n=234.
- Secondary RECIST: deposited `radiological_response` CR+PR vs SD+PD.
  Unevaluable / non-English labels were not recoded into CR/PR/SD/PD.

## GSE229353

- Per-patient 10x MTX from `GSE229353_RAW.tar`.
- Detection = fraction of barcodes with count > 0 for TACSTD2 or CLDN4.
- MPR from Hui et al., *npj Precis Oncol* 2023, Supplementary Table 1
  (public PDF). GEO has treatment (Chemo vs anti-PD1+Chemo) only.
- GSM7159185 / P03 MTX is **truncated on GEO** (gzip incomplete;
  332,846 / 4,500,643 declared nnz). P03 is an MPR library and was dropped.
  Remaining complete libraries: **2 MPR vs 4 non-MPR**.

## Statistics

Unit = patient. Two-sided Mann–Whitney U. Cliff’s δ and Hedges’ g are
MPR minus non-MPR (negative = lower leak in MPR). Meta-analysis is
inverse-variance DerSimonian–Laird random effects on Hedges’ g.
No multiple-testing theater. Cell-level p-values are not used.

## What this is not

- Not malignant TACSTD2/CLDN4.
- Not evidence that immune cells “express TROP2/claudin-4 as a resistance
  program.” Residual viable tumor is higher in non-MPR by definition;
  ambient RNA and misassigned epithelial barcodes are live alternatives
  that a CD45+ matrix cannot rule out.
