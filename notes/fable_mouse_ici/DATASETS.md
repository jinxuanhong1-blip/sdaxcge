# Dataset inventory & verification — mouse lung ICI RNA slice

All ten accessions were queried live against NCBI GEO (`acc.cgi`,
`series_matrix`) and EBI BioStudies/ArrayExpress on 2026-08-16. Every accession
exists, is *Mus musculus*, is lung-tumour RNA, and involves an
immune-checkpoint-inhibitor (ICI) context. Raw JSON evidence:
`accession_verification.json`, `file_inventory.json`, `sample_metadata.json`.

Task rule "skip files > 2 GB": **no supplementary file exceeded 2 GB**, so
nothing was skipped for size. Only processed data were used (never re-aligned
FASTQ).

| Accession | Organism | Model / tissue | ICI context | Data type used | Processed file | Size | Groups (n) | Stats? |
|---|---|---|---|---|---|---|---|---|
| GSE239485 | Mmu | LLC lung carcinoma, tumour | anti-PD-1 (+Poly I:C ± anti-C5aR1) | bulk RNA-seq, log2-normalised matrix | `GSE239485_Processed_data.xlsx` | 7.1 MB | Control(8), PolyIC+aPD1(8), PolyIC+aPD1+aC5aR1(8) | yes (n=8) |
| GSE297630 | Mmu | LLC, s.c. tumour | anti-PD-1 tolerant vs control | Clariom S array, per-sample log2 (RMA) | series family SOFT + `..._processed_data.xlsx` | 5.3 MB | Control(3), anti-PD-1(3) | yes (n=3) |
| E-MTAB-13704 | Mmu | Lung GEMM, in-situ tumour | anti-PD-L1 combinations | bulk RNA-seq raw counts → CPM | `GEMMS_raw_counts.csv` (+ SDRF) | 13.3 MB | vehicle(5), aPD-L1(5), ATRi(4), ATRi/aPD-L1(5), Cis/aPD-L1/aCTLA4(4), VEGFRi/aPD-L1(4) | yes |
| GSE241978 | Mmu | CMT167 lung carcinoma | AhR-KO (PD-L1/IDO axis) vs Cas9 control | bulk RNA-seq, normalised-log block | `..._Sherr_analysis_CMT_KO_vs_Cas9Ctrl.xlsx` | 14.4 MB | Control(3), AhR-KO(3) | yes (n=3) |
| GSE330658 | Mmu | Egfr-mutant lung tumour | anti-PD-L1/VEGF programme (RNA arm: PTX/anti-VEGF) | bulk RNA-seq per-sample TPM | 8× per-sample `.xlsx` in `_RAW.tar` | 40.6 MB | Control(2), PTX(2), anti-VEGF(2), PTX-anti-VEGF(2) | yes (n=2, low power) |
| GSE197260 | Mmu | Egfr-TKI lung tumour | gefitinib + anti-PD-1 (4H2) ± anti-VEGFR2 (DC101) | bulk RNA-seq TPM | `GSE197260_RNAseqTPM_MM_EGFR-TKI-CD8.txt.gz` | 0.3 MB | 7 arms, **1 sample each** | descriptive only |
| GSE133604 | Mmu | KrasG12D;p53−/− (KP) lung tumour | anti-PD-1 ± Asf1a-KO | scRNA-seq 10x → pseudobulk | `GSE133604_RAW.tar` (+ genes.tsv) | 77.7 MB | Ctrl, Ctrl+PD1, KO, KO+PD1 — **1 each** | descriptive only |
| GSE129297 | Mmu | SCLC lung tumour | anti-PD-1 + CDK7i (YKL) | scRNA-seq 10x (unfiltered) → pseudobulk | `GSE129297_RAW.tar` (+ features.tsv) | 155 MB | Ctrl, PD1, YKL, combo — **1 each** | descriptive only |
| GSE297632 | Mmu | LLC, s.c. tumour | anti-PD-1 tolerant vs control | scRNA-seq 10x → pseudobulk | `GSE297632_RAW.tar` | 231 MB | Control, anti-PD-1 — **1 each** | descriptive only |
| GSE222158 | Mmu | Lung cancer, sorted immune cells | anti-PD-1 + "21DC" ± combo | scRNA-seq 10x (CD45/CD3 sorts) → pseudobulk | `GSE222158_RAW.tar` | 162 MB | CD45: Ctl/PD1/21DC/Combo; CD3: Ctl/Combo — **1 each** | descriptive only |

## Per-dataset notes / caveats

- **GSE297630**: the processed `.xlsx` is a C-vs-P summary table keyed by Clariom
  transcript-cluster IDs with **no gene symbols**. Probes were mapped to Cldn4
  (`TC0500003384.mm.2`) and Tacstd2 (`TC0600002427.mm.2`, RefSeq NM_020047) via
  the platform `mrna_assignment` field inside the series family SOFT file; the
  same file provided the 6 per-sample log2 (RMA) values used for our own t-test.
  Our result agrees with the authors' reported Cldn4 fold change (−1.36,
  p≈9e-4).
- **GSE241978**: comparison is a genotype knockout (AhR-KO vs Cas9 control) in
  the PD-L1/IDO signalling axis, not an antibody-ICI arm. Tacstd2 sits at the
  dataset detection floor in all six samples (not expressed), so no Tacstd2 test
  is possible; Cldn4 is down in AhR-KO but not significant.
- **GSE330658**: the RNA arm deposited here contains Control/PTX/anti-VEGF/
  PTX+anti-VEGF (2 each); the anti-PD-L1 antibody arm of the study is not in this
  RNA subset. n=2 gives very low power — treat as descriptive-leaning.
- **GSE197260**: 7 conditions with a single library each (a treatment/time
  course), so no replication → descriptive only.
- **Single-cell datasets** (GSE133604, GSE129297, GSE297632, GSE222158): counts
  were summed across all cells per sample (pseudobulk) and CPM-normalised. Each
  condition is a single 10x library, so no inferential statistics are reported.
  GSE129297 matrices are **unfiltered** (737,280 barcodes incl. empty droplets),
  making its pseudobulk noisier. In GSE222158 the CD45+/CD3+ immune-sorted
  samples largely exclude epithelial tumour cells, so the epithelial markers
  Tacstd2/Cldn4 are expectedly low — a compartment caveat, not a treatment
  effect.

## Target genes
- **Tacstd2** (Trop2; Ensembl `ENSMUSG00000051397`; RefSeq NM_020047)
- **Cldn4** (claudin-4; Ensembl `ENSMUSG00000041378`; RefSeq NM_009903)
Matched case-insensitively on symbol/alias (Tacstd2/Trop2, Cldn4) for symbol-
keyed tables and by Ensembl ID for E-MTAB-13704.
