# GSE316655 is anti-LILRB2 scRNA, not a CLDN4 perturbation

Liu et al., *Science Immunology* (2026) 11:eadt7832 (PMID 41931598; PMC13282697). Data availability names one RNA-seq accession: **GSE316655**.

**No CLDN4 knockdown, knockout, or overexpression is in the GEO files or in the paper’s functional experiments.** The claudin that is genetically perturbed is **CLDN18.2 only**. Recombinant CLDN4 is one coated protein in a family-wide LILRB reporter screen. The screen’s stored CLDN4 values, and a CLDN4-only neighbor score on the one open cohort that was fully downloaded, are below. There is no CLDN4 EC50.

## Deposited files

GEO file list (`file_inventory.tsv`), public 2026-01-16. One supplementary archive plus 12 Cell Ranger files. No series matrix, no tumor-cell library, no CLDN genotype column.

| File | Bytes |
|---|---:|
| GSE316655_RAW.tar | 135,751,680 |
| GSM9457798 ND_anti_B2 barcodes / features / matrix | 56,954 / 593,990 / 45,372,697 |
| GSM9457799 ND_CTR barcodes / features / matrix | 46,506 / 685,596 / 39,712,106 |
| GSM9457800 PB_anti_B2 barcodes / features / matrix | 14,382 / 685,596 / 18,572,837 |
| GSM9457801 PB_CTR barcodes / features / matrix | 25,748 / 685,596 / 29,286,251 |

SRA study SRP662842 / BioProject PRJNA1403544. Four experiments, each with four paired-end NextSeq 500 runs (SRR36853505–SRR36853520).

| GSM | Library | GEO title | Tissue | Treatment | Protocol line | SRX |
|---|---|---|---|---|---|---|
| GSM9457798 | ND_anti_B2 | ND tumor LILRB2 treatment | tumor | anti-LILRB2 | SK-MEL-5 | SRX31826394 |
| GSM9457799 | ND_CTR | ND tumor Isotype | tumor | isotype | SK-MEL-5 | SRX31826395 |
| GSM9457800 | PB_anti_B2 | PB LILRB2 treatment | peripheral blood | anti-LILRB2 | SK-MEL-5 | SRX31826396 |
| GSM9457801 | PB_CTR | PB Isotype | peripheral blood | isotype | SK-MEL-5 | SRX31826397 |

Design, from the series and from the scRNA-seq methods paragraph: NSG-SGM3 mice reconstituted with human cord-blood CD34+ cells, then subcutaneous **SK-MEL-5**. After tumors reached 50–60 mm³, LALAPG-mutated anti-LILRB2 or isotype at 10 mg/kg twice weekly. Live single human CD45+ cells from tumor and peripheral blood were FACS-sorted and run on Chromium. The series-level overall design also names MIA PaCa-2 as an alternative inoculum; the sample protocols and the scRNA methods name SK-MEL-5 only. “ND” is the library prefix in the GEO title. It is not expanded in the SOFT file.

SOFT processing line: Cell Ranger **5.0.0**, assembly **GRCh38_and_mm10-2020-A**. The paper methods say Cell Ranger **8.0.0**. The deposited features match the 2020-A dual genome (36,601 `GRCh38_` genes + 32,285 `mm10___` genes = 68,886). Both version strings are recorded. They are not reconciled here.

Unfiltered Cell Ranger matrices, barcodes equal matrix columns:

| Library | Cells | Nonzero entries |
|---|---:|---:|
| ND_anti_B2 | 13,779 | 14,133,560 |
| ND_CTR | 9,086 | 10,240,782 |
| PB_anti_B2 | 2,596 | 4,892,215 |
| PB_CTR | 4,835 | 7,673,498 |
| Total | 30,296 | |

Human CLDN4 is on the reference (`GRCh38_CLDN4`, ENSG00000189143). Counts in these CD45+ libraries:

| Library | CLDN4 cells (UMIs) | CLDN18 cells (UMIs) | mm10 Cldn4 UMIs | EPCAM UMIs |
|---|---|---|---|---|
| ND_anti_B2 | 3 (3) | 0 (0) | 1 | 2 |
| ND_CTR | 2 (2) | 1 (1) | 0 | 2 |
| PB_anti_B2 | 0 (0) | 0 (0) | 0 | 4 |
| PB_CTR | 1 (1) | 0 (0) | 0 | 0 |

Human CLDN4 total is **6 UMIs / 30,296 cells**. That is reference presence in an immune sort, not a CLDN4 perturbation. Full gene table: `gene_detection.tsv`. These counts are not a cluster-level test of anti-LILRB2 versus isotype, and they are not an exclusion statistic. One library per arm.

## CLDN4 in the paper versus CLDN18.2

Main text and figure legends do not name CLDN4. Supplement table S4 does.

**Binding screen, not a perturbation.** Fig. 1A coats recombinant claudins and scores LILRB chimeric reporters. The text says all tested CLDN family members, including CLDN18.2 and its ECL1/ECL2, activated LILRB2 and LILRB5 and not the other LILRBs or mouse PirB/gp49B1. Table S4 is a reagent catalog, not an EC50 table. **CLDN4 is on that list** (Abnova H00001364-P01), next to CLDN1, 2, 3, 5–11, 14–17, 18.1, 18.2, 19, 20, and 22. There is no CLDN4 EC50, Kd, alanine scan, or coculture. EC50 and Kd numbers are for **CLDN18.2 ECL1/ECL2** only. The authors then write that CLDN18.2 was used as the representative claudin for subsequent studies. The GLW-C-C motif in ECL1 is described as family-conserved, and the discussion says the LILRB2 interface may be shared. That is their hypothesis. It is not a CLDN4 experiment.

**Genetic and peptide perturbations are CLDN18.2.** Plasmids section: human CLDN18.2 (NM_001002026) only. Constructs used in cells and mice:

- CLDN18.2 full length
- CLDN18.2-BD (binding-defective; deletion of residues 48–57 on the ΔECL2 background; G48A / L49A / W50A / C53A)
- CLDN18.2-ΔICD (binds LILRB2, lacks the C-terminal signaling tail)

Hosts: MIA PaCa-2, MC38, B16F10, LLC. Peptides: CLDN18.2 ECL1 and ECL2 only. In vivo antibody perturbation in the deposited scRNA is anti-LILRB2 versus isotype, on SK-MEL-5 tumors that are not described as CLDN18.2-transduced. The CLDN18.2-ΔICD humanized-mouse growth curves are MIA PaCa-2, and those tumors are not the GSE316655 libraries.

Table S2 lists qPCR primers for many claudins, including CLDN4. Table S3 lists an APC anti-CLDN4 antibody (R&D 382321). Figure S10A is the qPCR that uses those primers: CLDN mRNA in unmodified SK-MEL-5, GAPDH-normalized, three technical replicates. That is expression, not a knockdown.

## Quantitative CLDN4 rows

Table S4 has no EC50 and no reporter percent. The numbers are in the xlsx filed as `NIHMS2177808-supplement-Data_File_S1.xlsx` (sheets named by figure). The PDF filed as Data file S2 is western-blot images; extracted text has no CLDN4. The supplement captions swap those descriptions (S1 called immunoblots, S2 called tabulated data). Content, not the caption, is what was used.

Fig. 1A in that sheet stores one number per ligand × receptor. The block does not print a unit. Fig. 1B in the same sheet labels the dose curves “% GFP+ cells”, and the methods define this assay as percent GFP+ reporter cells. Replicates are not in the Fig. 1A block (one cell per pair). Source spellings CLND6, CLND9, CLND10, CLND14, CLND16, CLND17, CLND18.1, CLND19, CLND20, and CLND22 are kept in `fig1a_reporter_source.tsv`.

CLDN4 row, as stored: LILRB1 0, LILRB2 26, LILRB3 1, LILRB4 1, LILRB5 21, PirB 0, Gp49B 1. Negative control on LILRB2/LILRB5 is 3/1. Positive control is 53/57. CLDN18.2 is 32/15. ECL1 and ECL2 (the CLDN18.2 loops in this block) are 45/42 and 37/28.

Among the 20 full-length claudin rows, sorting the single stored LILRB2 numbers puts CLDN4 9th (26, behind CLDN1 42, CLDN15 38, CLDN3 37, CLDN14 32, CLDN18.2 32, CLDN6 30, CLDN10 28, CLDN5 27). On LILRB5, CLDN4 is tied for 5th (21, with CLDN5; CLDN3 is 49; CLDN18.2 is 15). The sheet has no test of CLDN4 versus CLDN18.2. Those ranks are a sort of one number per protein.

EC50 rows in the workbook (`ec50_source.tsv`): CLDN18.2 ECL1–LILRB2 293.2 nM, ECL2–LILRB2 1039 nM, ECL1–LILRB5 1725 nM, ECL2–LILRB5 7156 nM. Fig. 1F Kd is 175 nM for CLDN18.2 ECL1 versus LILRB2. No CLDN4 EC50 or Kd row exists.

Fig. S10A qPCR, SK-MEL-5, GAPDH-normalized, three technical replicates (`skmel5_cldn_qpcr_s10a.tsv`): CLDN4 is 1.483038, 1.211964, 0.904919 (mean of those three values 1.200). CLDN1, CLDN16, and CLDN19 are stored as N.D. CLDN6 is the highest of the panel (2.497311, 3.197886, 3.360633). This is baseline cell-line mRNA, not a perturbation and not a spatial measurement.

## CLDN4-only subset of the pooled proximity cohorts

Fig. 3H pools all CLDNs. The source sheet stores one enrichment score per sample and a log-rank p of 0.0164 for that pooled score versus good/poor outcome. SpaMAP’s Zenodo deposit (10.5281/zenodo.18252447) is the mapping script only. It has no per-gene output and no CLDN4 flag, so their z-score was not recomputed.

| Fig. 3H cohort | Samples in the source sheet | Gene-level spatial matrix | CLDN4-only score in this note |
|---|---|---|---|
| Ovarian (Denisenko; GSE211956) | GSM6506110, 6111, 6112, 6114, 6116, 6117 | Open Visium matrices and `tissue_positions_list.csv` | Scored |
| Breast (Wu; Zenodo 4739739) | CID4535, CID4290, CID44971, CID4465, 1142243F, 1160920F, plus D1, E1, M1–M8, M10 | The six CID/114/116 sections are the open Wu Visium deposit. D1, E1, and M* are not in that six-section record. Download of the Zenodo files returned HTTP 403 here | Not scored |
| Colorectal (Wang; GSE225857) | C1–C4, L1, L2 | Open GEO supplementary archive, 607 MB, not downloaded | Not scored |
| Kidney (Li) | PD45814, PD45816, PD47512, PD43948 cores | Visium raw data are EGA EGAD00001008781 (controlled). Not downloaded | Not scored |
| Pancreatic (Cui Zhou) | HTA12_22_4, HTA12_31_3, HTA12_24_6, HTA12_26_6, HTA12_29_3, HTA12_28_5 | No open GEO matrix located in this sweep. Not downloaded | Not scored |

Ovarian score, pre-specified in `scripts/ovarian_cldn4_lilrb2_neighbors.py` and not SpaMAP: in-tissue spots; index spots are EPCAM > 0; CLDN4-high versus low is a median split of CLDN4 on those spots (CLDN4 > 0 versus 0 when the median is 0); a neighbor is an in-tissue hex ring-1 spot with LILRB2 UMI > 0. GEO has no macrophage label, so LILRB2+ spots are not called macrophages. Self is excluded. Same-spot correlation is not the metric.

EPCAM > 0 is a loose gate (585–2,327 spots per section). Section deltas, CLDN4-high minus CLDN4-low mean LILRB2+ neighbor count:

| GSM | Fig. 3H outcome | High | Low | Δ |
|---|---|---:|---:|---:|
| GSM6506110 | Poor | 0.087 | 0.059 | +0.028 |
| GSM6506111 | Good | 0.258 | 0.247 | +0.011 |
| GSM6506112 | Good | 0.111 | 0.083 | +0.029 |
| GSM6506114 | Good | 0.596 | 0.448 | +0.148 |
| GSM6506116 | Poor | 0.115 | 0.058 | +0.058 |
| GSM6506117 | Poor | 0.092 | 0.079 | +0.013 |

Six of six sections are positive. Exact two-sided sign-flip p of the six deltas is 2/64 = 0.03125. The same sign appears in Good and Poor sections, so this proxy does not separate the Fig. 3H outcome labels. A pooled-CLDN sum on the same spots (CLDN1–12, 14–20, 22) is also positive in 6/6, with larger deltas (about +0.019 to +0.168). Direction is more LILRB2+ neighbors around CLDN4-high EPCAM+ spots, which is a proximity proxy, not cytotoxic exclusion, and not their SpaMAP enrichment score. Full table: `ovarian_cldn4_lilrb2_neighbors.tsv`.

An erratum was issued (Sci Immunol 2026 Jun 26, eaej3848, PMID 42361201). PubMed has no correction text. Nothing here assumes what it changes.

## Relevance to CLDN4 immune exclusion

This accession does not move the locked CLDN4 exclusion numbers (CosMx neighbor deficit without effector-transcript shutdown; concordant-4 malignant CLDN4+ versus T/NK; GSE137244 KL>KP Cldn4; TCGA keratin residual; TISMO Tacstd2). It is a different molecule, a different perturbation, and a different geometry.

What the paper actually shows, and how it sits next to exclusion:

1. **Contact suppression, for CLDN18.2.** CLDN18.2 ECL1 on MDSCs lowers T-cell proliferation and IFN-γ/TNF, raises STAT1/STAT3, PD-L1, and IDO1, and lowers NF-κB, through LILRB2 ITIMs and SHP-1. That is inhibition of immune cells that are in the coculture. The locked CosMx sentence is the other pattern: CLDN4-high tumor has fewer cytotoxic neighbors, and the effectors that are close are not lower for GZMB/PRF1/NKG7/IFNG. Do not use GSE316655 to rewrite that sentence.

2. **In vivo composition shift is also CLDN18.2.** In LILRB2-transgenic mice, MC38-CLDN18.2-ΔICD tumors grow faster with more MDSCs and fewer T, NK, and NKT cells; anti-LILRB2 reverses that. Fig. S10 describes the humanized scRNA (this GEO series) as LILRB2 blockade shifting myeloid cells and increasing memory T, CD8, and NK. Those are antibody-arm descriptions in a SK-MEL-5 CD45 sort. They are not a CLDN4-high versus CLDN4-low neighbor count. n = 1 library per arm in the deposit, so this note does not add a new paired test.

3. **Their spatial result is proximity of LILRB2-high macrophages to CLDN-positive cancer cells**, SpaMAP on five cohorts, all CLDNs pooled. The ovarian CLDN4-only spot-neighbor proxy above is the same direction (more LILRB2+ neighbors, 6/6) and is still not cytotoxic exclusion and not their statistic. The other four cohorts were not rescored.

4. **The matrices cannot host a tumor-side CLDN4 exclusion test.** EPCAM is 8 UMIs across 30,296 cells. There is no epithelial compartment and no CLDN4 genotype.

CLDN4 remains pinned by the private KD coculture, which is not in this repository. GSE316655 does not substitute for it. The public fact to keep is narrower: a recombinant CLDN4 protein was included in a family screen that activated LILRB2/LILRB5 reporters, and every functional tumor experiment after that screen used CLDN18.2.

```bash
python3 scripts/gse316655_inventory.py --mtx-dir /path/to/mtx --filelist /tmp/gse316655/filelist.txt
python3 scripts/ovarian_cldn4_lilrb2_neighbors.py --mtx-dir /path/to/GSE211956
```

Matrices are downloaded from the GEO supplementary URLs in `file_inventory.tsv` and are not committed.
