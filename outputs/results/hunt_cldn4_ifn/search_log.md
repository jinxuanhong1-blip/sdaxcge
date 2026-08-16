# Search and analysis log

## Question and eligibility

Question: among public transcriptome-scale experiments with direct CLDN4/Cldn4
KD, KO, shRNA, or siRNA, where do IFN, MHC-I, or antigen-processing machinery
(APM) genes rise, prioritizing IFI27, OAS2, IFIT1, MX1, ISG15, and HLA-A?

Included:

- Direct genetic/RNAi loss of CLDN4/Cldn4.
- Transcriptome-scale RNA-seq or expression array.
- Public accession plus downloadable gene-level processed values sufficient to
  establish direction.

Excluded:

- CLDN4 correlation, overexpression alone, ligand treatment, or perturbation of
  another claudin.
- Targeted qPCR/protein-only studies.
- A paper saying RNA-seq was done without public accession/matrix.
- Methylation, ChIP-seq, ATAC-seq, and non-expression assays.

## Repository search (2026-08-16)

NCBI GEO DataSets Entrez searches, restricted to GSE entry type:

```text
CLDN4[All Fields] AND gse[Entry Type]          21 hits
Cldn4[All Fields] AND gse[Entry Type]          21 hits
"claudin 4"[All Fields] AND gse[Entry Type]    10 hits
"claudin-4"[All Fields] AND gse[Entry Type]    10 hits
```

After unioning, there were 25 unique GSE records. Titles, summaries, overall
designs, sample/channel annotations, and supplementary files were screened.
The direct eligible records were GSE50927, GSE22493, and GSE207704.

Cross-checks:

- GEO FTP and SRA links for processed/raw availability.
- ArrayExpress/BioStudies and OmicsDI aliases.
- Web/PubMed combinations of `CLDN4`, `Cldn4`, `claudin 4`, `claudin-4` with
  `knockdown`, `knockout`, `KD`, `KO`, `siRNA`, `shRNA`, `RNA-seq`,
  `transcriptome`, `microarray`, `accession`, and `dataset`.
- LINCS GSE92742 and GSE70138 `sig_info` metadata were searched
  case-insensitively for CLDN4: zero records. The 24,187-signature L1000KD2
  curated CSV was also downloaded and had zero `pert_iname == CLDN4` rows.

This is a systematic indexed-resource search, not a guarantee that an
unindexed lab website, newly released record, controlled-access archive, or
unreported matrix does not exist.

## Direction and panels

All effects are oriented as CLDN4 loss / control. This was checked against
strongly negative CLDN4/Cldn4 values in each eligible dataset.

Priority panel:

```text
IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A
```

Human MHC-I/APM panel:

```text
HLA-A HLA-B HLA-C B2M NLRC5 TAP1 TAP2 TAPBP
PSMB8 PSMB9 PSMB10 ERAP1 ERAP2 CALR CANX PDIA3
```

Mouse proxies used:

```text
IFI27 -> Ifi27l2a (one-to-many family proxy; alternatives retained in report)
HLA-A -> H2-K1
HLA-B -> H2-D1
HLA-C -> H2-Q7
```

Other same-symbol mouse genes were matched case-insensitively.

## Dataset-specific processing

### GSE50927

Used submitter-deposited edgeR tables:

- `GSE50927_Cldn4lungWTvsKOgenes.csv.gz`
- `GSE50927_VILIwtkoloGenes.csv.gz`
- `GSE50927_VILIwtkohiGenes.csv.gz`

`logFC`, `PValue`, and BH `FDR` are reported exactly as deposited. No
reanalysis was possible because replicate-level counts are not deposited.

### GSE207704

Used `GSE207704_CLDN4_RNAseq.txt.gz`. The file has one FPKM column for each
cell-line/condition combination despite two replicate GSMs per group. For each
available gene, the displayed effect is:

```text
log2((KO FPKM + 0.1) / (WT FPKM + 0.1))
```

The 0.1 pseudocount only stabilizes low abundance. Duplicate gene records are
collapsed by median. No p-value is calculated from collapsed values.

### GSE22493

Used the GEO series matrix plus GPL10555 annotations in the family SOFT file.
The deposited values are normalized paired sample/control log-ratios from
three two-color arrays; positive is CLDN4-siRNA/control. Each gene summary is
the median of all available probe-replicate log-ratios. Probe-level values are
retained separately because HLA-A and several APM genes have discordant probes.
No p-value is invented.

## Interpretation rules

- Direction counts are descriptive, not enrichment tests.
- `FDR < 0.05` is only counted where the submitter deposited an FDR
  (GSE50927).
- A contrast is not called a broad MHC-I response merely because several
  chaperone/immunoproteasome genes rise; HLA-I heavy chains and TAP genes are
  assessed separately.
- Outcome-selected VILI-high lungs are explicitly considered confounded.
- Missing genes are `NA`, never treated as unchanged.

## Reproducibility

`analyze_public_sets.py` downloads immutable-named GEO processed files into
`_cache/` and rebuilds all quantitative TSVs using Python's standard library.
The cache is not evidence and may be deleted safely.
