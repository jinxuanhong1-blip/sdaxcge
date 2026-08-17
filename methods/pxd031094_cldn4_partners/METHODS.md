# Methods — PXD031094 CLDN4 CoIP partner ranking

## Source

- PRIDE [PXD031094](https://www.ebi.ac.uk/pride/archive/projects/PXD031094) (Suarez-Artiles et al., *Cell Rep.* 2022).
- Public **SEARCH** table only: `proteinGroups_Cldn4.txt` from  
  `https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/10/PXD031094/proteinGroups_Cldn4.txt`
- SHA-1 `4dbde6c0d827105412ef1f34cbb969b8dd8acb36` (PRIDE file checksum).
- RAW spectra (610 files) were not downloaded. PRISMA `proteinGroups_PRISMA.txt` was not used (peptide-matrix C-termini, not full-length CLDN4 bait CoIP).

`scripts/00_download.py` fetches the SEARCH file and checks the SHA-1.

## Filters

From the MaxQuant proteinGroups table:

- drop `Reverse`, `Potential contaminant`, `Only identified by site`
- keep `Unique peptides >= 2`
- keep LFQ > 0 in at least 2 of 4 `Cldn4_*` columns

Gene symbol = MaxQuant `Gene names` if present (human bait CLDN4), else first UniProt `GN=` in `Fasta headers`. DLA class I without `GN=` on the lead accession is labeled **DLA88** from a secondary header `GN=DLA88`. Empty `GN=` is left empty.

## Quantification and tests

- LFQ columns: `Cldn4_1..4` vs `GFPctrl_1..4`
- log2 transform; zeros imputed at global minimum positive LFQ / 2
- Welch two-sample t-test on the four vs four log2 values
- Simple variance-moderated t: residual variance shrunk toward the median with `df_prior=10` (not limma)
- Benjamini–Hochberg FDR over all groups that passed the peptide/contaminant filter

**Primary partner** = (FDR_mod < 0.05 and log2FC > 1) or exclusive (LFQ in ≥2 Cldn4 and 0 GFP).  
Rank: bait, then stringent, then exclusive, then log2FC.

This is a reanalysis of the deposited table. It is not the authors’ exact moderated-t / GFP-second-cutoff call (they reported 24 CLDN4 interactions).

## Category mapping

Annotation only. A protein must already be in the MaxQuant table.

| Label | Sets |
| --- | --- |
| IFN | Hallmark IFNα / IFNγ; Reactome interferon signaling |
| MHC | Reactome MHC class I peptide loading (29 genes) + curated DLA88 / HLA / B2M / TAP* |
| TJ | Reactome tight junction; KEGG tight junction; GO TJ assembly/organization; Hallmark apical junction; curated CLDN/TJP/OCLN/PARD/MPP/PATJ |
| trafficking | Reactome vesicle-mediated transport, membrane trafficking, RAB regulation, clathrin endocytosis, ER–Golgi, lysosome vesicle; plus RAB* symbols |
| immune | Hallmark inflammatory / allograft / complement / IL6-JAK-STAT3, plus IFN and MHC |

Gene-set JSON: `data/category_sets.json` (MSigDB v2023.2 Hs subsets + a8 TJ sets already in this repo).

## Software

Python 3: pandas, numpy, scipy, statsmodels.  
`scripts/01_rank_partners.py` writes the result tables.

## Out of scope

ICI response, tumor abundance, RAW re-search, invented gene names.
