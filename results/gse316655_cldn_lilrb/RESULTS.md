# GSE316655 is anti-LILRB2 scRNA, not a CLDN4 perturbation

Liu et al., *Science Immunology* (2026) 11:eadt7832 (PMID 41931598; PMC13282697). Data availability names one RNA-seq accession: **GSE316655**.

**No CLDN4 knockdown, knockout, or overexpression is in the GEO files or in the paper’s functional experiments.** The claudin that is genetically perturbed is **CLDN18.2 only**. Recombinant CLDN4 is one coated protein in a family-wide LILRB reporter screen. That screen is not a CLDN4 match to immune exclusion.

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

**Binding screen, not a perturbation.** Fig. 1A coats recombinant claudins and scores LILRB chimeric reporters. The text says all tested CLDN family members, including CLDN18.2 and its ECL1/ECL2, activated LILRB2 and LILRB5 and not the other LILRBs or mouse PirB/gp49B1. Table S4 lists the proteins. **CLDN4 is on that list** (Abnova H00001364-P01), next to CLDN1, 2, 3, 5–11, 14–17, 18.1, 18.2, 19, 20, and 22. The paper does not report a CLDN4 EC50, Kd, alanine scan, or coculture. EC50 and Kd numbers in the text are for **CLDN18.2 ECL1/ECL2** (ECL1–LILRB2 EC50 293.2 nM; BLI Kd 175 nM). The authors then write that CLDN18.2 was used as the representative claudin for subsequent studies. The GLW-C-C motif in ECL1 is described as family-conserved, and the discussion says the LILRB2 interface may be shared. That is their hypothesis. It is not a CLDN4 experiment.

**Genetic and peptide perturbations are CLDN18.2.** Plasmids section: human CLDN18.2 (NM_001002026) only. Constructs used in cells and mice:

- CLDN18.2 full length
- CLDN18.2-BD (binding-defective; deletion of residues 48–57 on the ΔECL2 background; G48A / L49A / W50A / C53A)
- CLDN18.2-ΔICD (binds LILRB2, lacks the C-terminal signaling tail)

Hosts: MIA PaCa-2, MC38, B16F10, LLC. Peptides: CLDN18.2 ECL1 and ECL2 only. In vivo antibody perturbation in the deposited scRNA is anti-LILRB2 versus isotype, on SK-MEL-5 tumors that are not described as CLDN18.2-transduced. The CLDN18.2-ΔICD humanized-mouse growth curves are MIA PaCa-2, and those tumors are not the GSE316655 libraries.

Table S2 lists qPCR primers for many claudins, including CLDN4. Table S3 lists an APC anti-CLDN4 antibody (R&D 382321). No figure legend reports a CLDN4 qPCR or CLDN4 flow result.

An erratum was issued (Sci Immunol 2026 Jun 26, eaej3848, PMID 42361201). PubMed has no correction text. Nothing here assumes what it changes.

## Relevance to CLDN4 immune exclusion

This accession does not move the locked CLDN4 exclusion numbers (CosMx neighbor deficit without effector-transcript shutdown; concordant-4 malignant CLDN4+ versus T/NK; GSE137244 KL>KP Cldn4; TCGA keratin residual; TISMO Tacstd2). It is a different molecule, a different perturbation, and a different geometry.

What the paper actually shows, and how it sits next to exclusion:

1. **Contact suppression, for CLDN18.2.** CLDN18.2 ECL1 on MDSCs lowers T-cell proliferation and IFN-γ/TNF, raises STAT1/STAT3, PD-L1, and IDO1, and lowers NF-κB, through LILRB2 ITIMs and SHP-1. That is inhibition of immune cells that are in the coculture. The locked CosMx sentence is the other pattern: CLDN4-high tumor has fewer cytotoxic neighbors, and the effectors that are close are not lower for GZMB/PRF1/NKG7/IFNG. Do not use GSE316655 to rewrite that sentence.

2. **In vivo composition shift is also CLDN18.2.** In LILRB2-transgenic mice, MC38-CLDN18.2-ΔICD tumors grow faster with more MDSCs and fewer T, NK, and NKT cells; anti-LILRB2 reverses that. Fig. S10 describes the humanized scRNA (this GEO series) as LILRB2 blockade shifting myeloid cells and increasing memory T, CD8, and NK. Those are antibody-arm descriptions in a SK-MEL-5 CD45 sort. They are not a CLDN4-high versus CLDN4-low neighbor count. n = 1 library per arm in the deposit, so this note does not add a new paired test.

3. **Their spatial result is proximity of LILRB2-high macrophages to CLDN-positive cancer cells**, SpaMAP on five published cohorts (pancreatic, kidney, colorectal, ovarian, breast), with all CLDNs pooled, associated with worse outcome. That is shorter distance, not exclusion, and it is not CLDN4-specific. It is not entered here as a Visium same-spot exclusion result.

4. **The matrices cannot host a tumor-side CLDN4 exclusion test.** EPCAM is 8 UMIs across 30,296 cells. There is no epithelial compartment and no CLDN4 genotype.

CLDN4 remains pinned by the private KD coculture, which is not in this repository. GSE316655 does not substitute for it. The public fact to keep is narrower: a recombinant CLDN4 protein was included in a family screen that activated LILRB2/LILRB5 reporters, and every functional tumor experiment after that screen used CLDN18.2.

```bash
python3 scripts/gse316655_inventory.py --mtx-dir /path/to/mtx --filelist /tmp/gse316655/filelist.txt
```

Matrices are downloaded from the GEO supplementary URLs in `file_inventory.tsv` and are not committed.
