# GSE243013 results: TACSTD2 / CLDN4 vs MPR

Computed from the public processed immune count matrix. Nothing here is a tumor-epithelial TACSTD2/CLDN4 measurement.

## Data that were used

| Object | Role |
| --- | --- |
| `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz` | Processed counts (cells × genes, 1,254,749 × 31,831, 2,010,550,708 nonzero entries) |
| `GSE243013_genes.csv.gz` | Gene order. TACSTD2 = column 1057. CLDN4 = column 11797. Each symbol once. |
| `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz` | Cell barcodes, `sampleID`, `pathological_response`, library size |
| `GSE243013_series_matrix.txt.gz` | Empty expression table (`data_row_count` = 0). Not used for expression. |
| `GSE243013_RAW.tar` | Per-sample TCR archives only. Not an expression matrix. |

Full MTX scan: declared nnz = scanned nnz = 2,010,550,708. Kept 22,966 nonzero entries for the two genes. Extract: `extracted_tacstd2_cldn4_cells.tsv`.

GEO metadata has **243** `sampleID`s. The paper reports **234** patients. Tests use deposited labels, not a forced n=234. One sample (`P433`, label `unknowm`) was excluded. Analyzed n = **242** (pCR 85, MPR 45, non-MPR 112). Primary contrast treats GEO `pCR` + `MPR` as MPR-any (n=130), matching the paper’s MPR definition (RVT ≤10%, pCR ⊂ MPR).

## Detection (immune compartment)

| Gene | Positive cells | Rate | Patients with ≥1 positive cell (of 242) | Cells positive for both genes |
| --- | ---: | ---: | ---: | ---: |
| TACSTD2 | 15,084 | 1.202% | 239 | 461 |
| CLDN4 | 7,882 | 0.628% | 228 | 461 |

Among nonzero entries, the median count is 1 for both genes (TACSTD2 max 133; CLDN4 max 22). This is sparse immune-compartment signal, not a bulk tumor score.

TACSTD2-positive cells are enriched in the myeloid annotation (8,336 / 191,099 = 4.36%) versus T/NK (0.64%) and B (0.62%). The most frequent TACSTD2+ subtypes are `cDC2_CD1C` and `cDC1_CLEC9A`. CLDN4 is more evenly spread (myeloid 1.02%, B 0.58%, T/NK 0.55%). These are deposited annotations, not a re-clustering.

Patient-level TACSTD2 and CLDN4 detection fractions are correlated (Spearman ρ = 0.791, n=242). That is compatible with a shared sample-level source (including residual tumor / ambient RNA) and does not prove independent immune programs.

## Primary test: MPR-any vs non-MPR

Two-sided Mann–Whitney U. Positive Cliff’s delta means higher in MPR-any.

| Metric | Median MPR-any (n=130) | Median non-MPR (n=112) | p | Cliff’s δ | AUC |
| --- | ---: | ---: | ---: | ---: | ---: |
| TACSTD2 fraction positive | 0.00618 | 0.01085 | 9.76×10⁻⁴ | −0.246 | 0.377 |
| TACSTD2 mean CPM | 2.45 | 3.48 | 5.98×10⁻³ | −0.205 | 0.397 |
| CLDN4 fraction positive | 0.00204 | 0.00523 | 4.79×10⁻⁵ | −0.303 | 0.348 |
| CLDN4 mean CPM | 0.638 | 1.64 | 1.21×10⁻⁴ | −0.287 | 0.357 |

Both genes are **lower** in MPR-any than in non-MPR. Mean raw count gives the same direction (TACSTD2 p=0.00385; CLDN4 p=2.82×10⁻⁵). pCR vs MPR-excluding-pCR is not significant for any metric (all p>0.42).

Patients with detection fraction ≥0.05: TACSTD2 n=14 (11 non-MPR); CLDN4 n=9 (all non-MPR). The shift is not only those tails (medians already differ), but a few high non-MPR samples contribute to the mean–median gap.

## Histology (secondary)

GEO `cancer_type` in the 242: LUSC 179 (MPR-any 109, non-MPR 70); LUAD 63 (MPR-any 21, non-MPR 42). Direction is the same in both histologies. LUAD n is small.

| Histology | Metric | Median MPR-any | Median non-MPR | p | Cliff’s δ |
| --- | --- | ---: | ---: | ---: | ---: |
| LUSC | TACSTD2 frac | 0.00692 | 0.01071 | 0.0220 | −0.203 |
| LUSC | CLDN4 frac | 0.00211 | 0.00525 | 2.00×10⁻⁴ | −0.330 |
| LUAD | TACSTD2 frac | 0.00452 | 0.01219 | 0.00437 | −0.444 |
| LUAD | CLDN4 frac | 0.00113 | 0.00493 | 0.0204 | −0.362 |

## What this is not

- Not cancer-cell TACSTD2 or CLDN4. The matrix is CD45+ immune scRNA after neoadjuvant chemo-immunotherapy.
- Not a skip. A processed matrix and MPR labels are public.
- Not evidence that immune cells biologically “express TROP2/claudin-4 as a resistance program.” Residual viable tumor is higher in non-MPR by definition; ambient/misassigned epithelial transcripts are a live alternative explanation that this matrix cannot rule out.
- Not the paper’s exact 234-patient analysis set. Sample IDs were not subsetted to match the publication.

Tables: `patient_level_tacstd2_cldn4.tsv`, `stats_tacstd2_cldn4_vs_mpr.tsv`, `stats_histology_stratified.tsv`, `celltype_detection.tsv`, `detection_summary.tsv`.
