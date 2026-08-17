# Methods — CLDN4-only AUCell / TF-target proxy (GSE148071)

Additive folder. Prior GSE207422 SCENIC / CLDN4-TJ-IFN write-ups are taken as given and are **not** re-run. **This folder uses GSE148071 only** and splits putative malignant cells on **CLDN4 only**.

## Dataset

GEO `GSE148071` (Wu et al., *Nature Communications* 2021, PMID 33953163). Forty-two diagnostic biopsies from stage III/IV NSCLC (Singleron GEXSCOPE). Public processed files: `GSE148071_RAW.tar` (42 `GSM*_P*_exp.txt.gz` count matrices) + series matrix (age/sex). Author Seurat / inferCNV / CopyKAT objects are **not** on GEO.

GSE207422, GSE131907, and raw FASTQ were not used.

## What was run (and what was not)

| Layer | What it is | What it is not |
| --- | --- | --- |
| Public prior | Union of TRRUST v2 + DoRothEA (OmniPath) + CollecTRI (OmniPath) targets per TF | Lung ChIP, motif scan, or a binding map |
| AUCell | Aibar 2017 linear recovery of set members in the top 5% ranks of `log1p(CP10k)` | R `AUCell` binary / pySCENIC CLI |
| Gene-set AUCell | Same recovery curve on curated ISG / MHC-I APM / TJ lists | A TF regulon |
| Mean module | Mean `log1p(CP10k)` of the same lists | SCENIC |

**pySCENIC was not run.** cisTarget feather rankings are multi-GB. Full `grn → ctx → aucell` would be used to *imply* motif support; that is too close to a binding claim we cannot make on this public archive.

**CLDN4 is held out** of every regulon and of the TJ gene set (circular otherwise). Public priors list CLDN4 as a GRHL2 and STAT2 target; those edges are dropped before scoring. TACSTD2 is recorded and is **never** a gate.

## TFs (pre-specified)

- **IFN:** STAT1, STAT2, IRF1, IRF3, IRF7, IRF9
- **MHC:** NLRC5, CIITA, RFX5, RFXANK, RFXAP (IRF1 also scored in an MHC+IRF1 union)
- **TJ:** GRHL1, GRHL2, ELF3, KLF4, OVOL1, OVOL2
- **Control TF:** HIF1A (hypoxia; not a TJ/IFN claim)
- **Control gene set:** compact OXPHOS list

A TF-specific prior regulon is scored only if ≥5 targets are present in that sample’s ranking universe. Combined **IFN / MHC / TJ unions** are the primary GRN-proxy tests.

## Lineage and malignant proxy

No published barcode labels are in the GEO archive. Lineage = argmax of marker-module `log1p(CP10k)` (same panel as the GSE148071 CellChat folder). **Putative malignant** = epithelial lineage. A stricter **mal-like** slice drops the top quartile of a normal-lung score (surfactant / club / ciliated genes) within each sample’s epithelium. That is a proxy, not CopyKAT.

QC: cells with total UMI < 200 are dropped.

## CLDN4 split and inference

Within each patient, epithelial (or mal-like) cells are split on CLDN4 `log1p(CP10k)` **Q4 vs Q1**. Primary test = paired Wilcoxon of mean AUCell (high vs low) across patients. **Unit of inference is the patient.**

Eligibility floor (locked): ≥20 cells in the compartment and ≥8 cells in each tail **and** a real CLDN4 split (`Q3 > Q1`). Patients with all-zero CLDN4 quartiles meet the cell-count floor but have no high/low contrast; they are reported and excluded. Cell-level Spearman is exploratory (pseudoreplication).

## Honest n

Header n is 42 deposited biopsies. The claim n is the number of patients meeting the occupancy floor, not 42. One biopsy can dominate cell-pooled statistics; those are not the primary test.

## Software

Python (numpy / pandas / scipy). Public processed counts only.
