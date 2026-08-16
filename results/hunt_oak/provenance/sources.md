# Sources & provenance

## Papers
- **Bessede A, Marabelle A, Guégan JP, et al.** TROP2 Is Associated with Primary
  Resistance to Immune Checkpoint Inhibition in Patients with Advanced NSCLC.
  *Clin Cancer Res* 2024;30(4):779-785. doi:10.1158/1078-0432.CCR-23-2566.
  (Analyzed 891 POPLAR+OAK tumors; TACSTD2/TROP2-high -> worse PFS/OS on atezo but
  not docetaxel; HR≈1.31 PFS, ≈1.26 OS; TROP2-high = reduced T-cell infiltration.)
- **Patil NS, Nabet BY, Müller S, et al.** Intratumoral plasma cells predict
  outcomes to PD-L1 blockade in NSCLC. *Cancer Cell* 2022;40(3):289-300.
  (Deposited the OAK/POPLAR RNA-seq at EGA: EGAS00001005013.)

## Controlled (NOT usable openly) — the actual OAK/POPLAR data
- EGA study **EGAS00001005013** and all sub-datasets — see `ega_access_map.md`.
  DAC: EGAC00001002120 (Genentech). Also gated: vivli.org (Roche).

## Open leftovers of OAK/POPLAR that exist but cannot test TACSTD2
- **Bessede 2024 supplements** (AACR figshare collection doi:10.1158/1078-0432.c.7077754):
  S1–S4 + Fig S1–S2. Downloaded 2026-08-16. Extracted text:
  `data/derived/bessede_supplements_extracted.md`. Aggregated characteristics
  and multivariate HRs only — no per-sample TACSTD2, no per-sample outcomes.
- **Gandara 2018 Nat Med ESM** (doi:10.1038/s41591-018-0134-3)
  `41591_2018_134_MOESM3_ESM.xlsx`: blood-TMB variants + OAK/POPLAR clinical
  (arm, PFS, OS). Schema: `data/derived/gandara2018_schema.json`.
  0 TACSTD2 variant rows; no RNA-seq. Used by github.com/Askir/pbmf-reproduction.

## Open data actually used in this hunt (proxies, not OAK)
### IMvigor210 — atezolizumab (anti-PD-L1), metastatic urothelial carcinoma, n=348
- Package: **IMvigor210CoreBiologies** (Mariathasan S, et al. *Nature* 2018;554:544-548,
  doi:10.1038/nature25501). Original host research-pub.gene.com/IMvigor210CoreBiologies
  is **offline (HTTP 404 as of 2026-08-16)**.
- Authoritative object used: `data/cds.RData` (a DESeq `CountDataSet` with raw counts,
  Entrez/symbol feature annotation, DESeq sizeFactors, and clinical pData incl.
  `binaryResponse`, `os`, `censOS`) from community mirror
  https://github.com/SiYangming/IMvigor210CoreBiologies (file `data/cds.RData`).
  sha256 in `sha256.txt`.
- NOTE: a CSV mirror (github.com/mimifp/tfm_mUC) was tried first and **rejected** —
  its expression matrix rows are censored to `gene_N` labels whose order does not
  match its re-sorted `fData`, so gene identity is lost (housekeeping ACTB summed to
  ~36 counts across 348 samples). Gene-level extraction from it is invalid.

### GSE135222 — NSCLC, anti-PD-1/PD-L1, n=27
- GEO **GSE135222** (Jung H, et al. *Nat Commun* 2019, DNA methylation / immune
  evasion RNA-seq subset). TPM matrix `GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz`;
  PFS event + pfs.time parsed from `GSE135222_series_matrix.txt`.
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222

## Reproduce
```
bash    results/hunt_oak/scripts/fetch_data.sh
Rscript results/hunt_oak/scripts/00_extract_imvigor210_cds.R
python3 results/hunt_oak/scripts/01_imvigor210.py
python3 results/hunt_oak/scripts/02_gse135222.py
python3 results/hunt_oak/scripts/03_figures.py
```
Gene identifiers: TACSTD2 = Entrez 4070 = ENSG00000184292; CD274/PD-L1 = 29126 /
ENSG00000120217; CD8A = 925 / ENSG00000153563.
