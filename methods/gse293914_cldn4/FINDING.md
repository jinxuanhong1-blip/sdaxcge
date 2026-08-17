# ADDITIVE — GSE293914 H1975 EGFR xenograft (anti-PD-1 ± anti-CCL20)

**Label: XENOGRAFT** (humanized NSG + human PBMC + NCI-H1975, EGFR L858R/T790M). Not a patient tumor.

Public **series matrix + MTX only**. Kwok et al., PMID [42549034](https://pubmed.ncbi.nlm.nih.gov/42549034). GEO [GSE293914](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293914) / [GSM8893263](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8893263). One Cell Ranger 7.1 Flex/fixed-RNA **raw** matrix (GRCh38). Series matrix has **0 expression rows**.

GEO sample title lists isotype, anti-PD-1, anti-CCL20, and combination. Characteristics field is `treatment: isotype`. Protocol is 10x Flex multiplex (target 20,000 nuclei/sample). Four dominant Flex probe barcodes (BC001–BC004) are present. **GEO does not map probe barcode → treatment.** SRA was not used.

## Honest n

| | Count |
|---|---:|
| GEO samples / libraries | **1** |
| Biological n for any contrast | **1** (do not treat 4 probe groups as n=4) |
| Paper growth assay (not this MTX) | n=6 mice/group |
| Raw barcodes | 2,032,489 |
| QC cells (genes 200–7999, MT% <10, author filter) | 42,911 |
| Tumor (epithelial argmax) | 38,110 |
| Lineage T/NK | **5** (T=3, NK=2) |
| UMI proxy CD3E/CD8A/NKG7 ≥1 | 79 |
| PTPRC (CD45) UMI ≥1 | **2** |

T/NK are **not present** as a usable compartment. The five lineage calls include EGFR-high barcodes and are consistent with misassigned tumor / empty-ish nuclei. PTPRC is effectively zero. “Plasma / fibroblast / myeloid” calls still carry TACSTD2/CLDN4/EGFR and are not trusted as human immune or mouse stroma (matrix is human-only).

## Human CLDN4 / TACSTD2

Unit is the **single xenograft library**. Cell-level p-values are not reported (pseudoreplication). No treatment test.

| Compartment | n cells | TACSTD2 mean log1p(CP10k) | TACSTD2 %UMI≥1 | CLDN4 mean log1p(CP10k) | CLDN4 %UMI≥1 |
|---|---:|---:|---:|---:|---:|
| Tumor | 38,110 | **1.23** | **50.7** | **1.19** | **49.9** |
| T/NK | 5 | 0.65 | 20 | **0** | **0** |

Tumor EGFR mean log1p(CP10k) = 0.85 (sanity that the epithelial compartment is H1975-like).

### Unlabeled Flex arms (descriptive only; not treatments)

| Probe | QC cells | Tumor n | T/NK n | Tumor TACSTD2 mean | Tumor CLDN4 mean |
|---|---:|---:|---:|---:|---:|
| BC001 | 4,154 | 3,651 | 0 | 0.86 | 1.04 |
| BC002 | 8,503 | 7,015 | 1 | 1.52 | 1.39 |
| BC003 | 13,037 | 11,659 | 4 | 1.16 | 1.18 |
| BC004 | 17,217 | 15,785 | 0 | 1.23 | 1.13 |

BC002 is highest for both genes. That is **not** a treatment effect: the arm is unlabeled, and n remains 1 library.

## Read this as additive, not as an ICI contrast

- **XENOGRAFT.** H1975 cell-line tumor in a humanized mouse. Do not pool with patient scRNA.
- **n=1.** One GSM. Four probe barcodes ≠ four biological replicates. Probe→isotype/aPD-1/aCCL20/combo map is absent from the public series matrix and MTX.
- **Tumor CLDN4 and TACSTD2 are on** in this H1975 xenograft (~half of epithelial nuclei UMI≥1; mean log1p(CP10k) ≈ 1.2).
- **T/NK are not present** (5 cells; PTPRC≈0). Tumor vs T/NK is not a real contrast here.
- No anti-PD-1 ± anti-CCL20 score is possible from the public files.

## Files

- `fig_xenograft_cldn4_tacstd2.png` — tumor vs T/NK (n=5) and unlabeled probe arms
- `stats.tsv` / `per_probe.tsv` / `lineage_counts.tsv` / `summary.json`
- `cell_calls.tsv.gz` — QC barcodes
- `analyze.py` / `METHODS.md`
- `metadata/GSE293914_series_matrix.txt.gz`
