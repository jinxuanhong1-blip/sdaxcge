# Cell painting and DNA-damage morphology: CLDN4 inventory

Exploratory public catalog. No same-well CLDN4 RNA or protein measurement sits next to a DNA-damage morphology assay in the sets checked here. Private 8-KL matrices were not used.

Cell Painting’s DNA channel (Hoechst or DAPI) is nuclear morphology. It is not a γH2AX or comet readout. Feature counts from that channel are not a DNA-damage call.

## Where CLDN4 is a named perturbation

Three public Cell Painting resources name CLDN4 as a genetic reagent. None of them measure CLDN4 expression in the imaged wells.

**JUMP (cpg0016), U2OS.** `orf.csv.gz` (15,132 reagents) contains one CLDN4 ORF: `JCP2022_905637`, transcript `NM_001305.4`, pLX_304, insert 627 nt, protein match 99.5. `crispr.csv.gz` (7,977 genes) contains one CLDN4 knockout id: `JCP2022_801379`. In `well.csv.gz` (1,151,808 wells) each id is plated five times:

- ORF, source_4, plates `BR00123512`–`BR00123516`, well B05
- CRISPR, source_13, plates `CP-CC9-R1-06`–`CP-CC9-R5-06`, well G18

The JUMP-Target-1 ORF plate (175 rows), which is the compound-target subset, does not include CLDN4. The Morphmap paper attaches DepMap 23Q4 U2OS expression as a gene-level annotation. That matrix was not re-downloaded, and it is not a per-well transcriptome.

**RxRx3-core, HUVEC.** `metadata_rxrx3_core.csv` (222,601 wells) labels CLDN4 on 89 CRISPR wells, all `Query guides`, experiments `gene-003` (45) and `gene-176` (44). Five guides: `CLDN4_guide_1` through `_4` have 18 wells each; `CLDN4_guide_5` has 17. There is no expression column. Embeddings were not scored.

**PERISCOPE (cpg0021), A549 and HeLa.** CLDN4 is in the screen. `A549_NGS_Counts.csv` has four CLDN4 barcodes (counts 19, 19, 14, and 2). The authors’ six published 1% FDR hit lists do not contain CLDN4 (A549 whole-cell 399 genes, A549 compartment 690, HeLa DMEM whole-cell 1,039, HeLa DMEM compartment 891, HeLa HPLM whole-cell 597, HeLa HPLM compartment 956). CLDN6 is on the HeLa DMEM whole-cell list, so a claudin can be called; CLDN4 was not.

HeLa per-feature tables place CLDN4 among expressed genes, not the 0 TPM table. Counts of features the authors called significant:

| screen | Mito | ConA | DAPI | WGA | Phalloidin | Sum | median Sum | fraction of genes with Sum at or below CLDN4 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HeLa DMEM | 19 | 9 | 14 | 10 | 33 | 88 | 99 | 0.393 |
| HeLa HPLM | 29 | 17 | 28 | 57 | 31 | 157 | 142 | 0.605 |

Those sums sit near the middle of 16,920 expressed genes. They are not a 1% FDR gene hit, and the DAPI column is a DNA-stain feature count.

The only CLDN4 expression annotation in PERISCOPE is baseline DepMap Public 21Q4 `CCLE_expression.csv`, subsetted by the authors to A549 (`ACH-000681`) and HeLa (`ACH-001086`). Stored values, which that file reports as log2(TPM+1): A549 2.118 (TPM 3.34), HeLa 4.145 (TPM 16.69). This is not RNA from the imaging wells.

## Expression paired with Cell Painting, without a measured CLDN4 landmark

LINCS 2020 `geneinfo_beta.txt` (12,328 genes) lists CLDN4 (`gene_id` 1364, ENSG00000189143) as **best inferred**, not as one of the 978 landmark genes. Haghighi et al. 2022 (cpg0003) matched Cell Painting to L1000 for CDRP, CDRP-bio, the LUAD allele set, TA-ORF, and the LINCS pilot, and describe the released expression profiles as landmark genes. The 8.5 GB cpg0003 matrices were not opened, so this pass does not claim those files were searched row by row. A measured CLDN4 column is not expected in a landmark-only matrix.

Rohban TA-ORF (cpg0017 / BBBC037): the 302-row sequence-match table has no CLDN4, and the eLife 24060 text does not mention it. Caicedo LUAD alleles (cpg0031): `gene_funcs.csv` has 590 allele rows and 135 gene prefixes, and CLDN4 is not one of them. Any CLDN4 value in the matched L1000 matrix would be an inferred response to other alleles, not a CLDN4 perturbation.

## DNA-damage morphology

**idr0080 Cell Health** (Way et al. 2021) is the public set that actually measures γH2AX together with Cell Painting. CRISPR in A549, ES2, and HCC44; the cell-cycle panel counts γH2AX spots. Supplementary table 1 has 127 rows and no CLDN4. There is no expression matrix. γH2AX values predicted for other Cell Painting compound sets (LINCS cpg0004; later models applied to cpg0012) are model outputs, not measurements, and they are not CLDN4 expression.

Other morphology sets checked only as catalogs, without a CLDN4 expression matrix:

- cpg0012, Bray/Wawer, ~30,000 compounds in U2OS (GigaDB 100351). DNA channel is not γH2AX.
- cpg0004, 1,571 compounds at 6 doses in A549.
- idr0133, Dahlin U2OS injury reference (218 concentration-response compounds and 283 single-dose compounds). Cell Painting injury profiles, not γH2AX.
- cpg0022 cmQTL, 297 iPSC donors. Morphology plus controlled-access WGS (`phs002032.v1.p1`). No public matched RNA-seq.
- Zenodo 7673199, γH2AX / Hoechst valinomycin pilot. Record does not list CLDN4; the archive was not unpacked.
- Zenodo 13683162, nanomaterial HTS with DAPI, γH2AX, 8-oxoG, and caspase. Record does not list CLDN4; the archive was not unpacked.

## Not found

No dataset in this pass has CLDN4 RNA or protein measured in the same wells as γH2AX, a comet assay, or 8-oxoG. Image archives of tens of terabytes were not downloaded, and CLDN4 profiles in JUMP or RxRx3 were not scored. Hoechst texture was not treated as DNA damage.

Tables: `results/tables/inventory_cellpainting.tsv`, `results/tables/cldn4_checks.tsv`. Rebuild with `python3 scripts/inventory_cellpainting_cldn4.py`.
