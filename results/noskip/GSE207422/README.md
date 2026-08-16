# GSE207422 malignant TACSTD2 (no-skip)

Hu et al., *Genome Medicine* 2023 (PMID 36869384; GEO GSE207422).
Neoadjuvant PD-1 + chemotherapy NSCLC scRNA-seq. Full author-processed UMI matrix analyzed (nothing skipped for size).

## Verdict

- **Claimed per-patient malignant TACSTD2 vs T/NK ρ = −0.40 to −0.50 is not supported.** Post-treatment Spearman ρ = **−0.021** (p = 0.95, n = 12) for malignant mean log1p(CP10K) vs T/NK fraction. All other framings are between ρ = −0.02 and −0.30 and none are significant.
- **NMPR > MPR malignant TACSTD2 is directional only.** Post-treatment malignant mean log1p(CP10K): NMPR median 1.72 (n=8) vs MPR median 1.33 (n=4, pCR P06 counted as MPR); Mann-Whitney p = 0.21. Percent TACSTD2+ is essentially identical (88% vs 86%, p = 0.81). Underpowered.
- **TACSTD2 is malignant-restricted vs T/NK.** In all 15 patients, malignant ≫ T/NK (paired Wilcoxon p = 6.1×10⁻⁵; median 86.5% vs 0.69% positive).

## Data used (nothing skipped for size)

| File | Role | Size |
|---|---|---|
| `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` | Author UMI counts, 24,292 genes × **92,330 cells** | 175.5 MB |
| `GSE207422_NSCLC_scRNAseq_metadata.xlsx` | Sample-level clinical labels (15 patients) | 11 KB |
| Paper Additional files 1–4 | Clinical tables + module gene lists | no barcodes |

- Matrix cell count and median genes/cell (1,253) match the paper (92,330 cells; median 1,256 genes).
- Raw FASTQ lives in GSA-Human HRA001033. Not needed; the processed matrix is complete.
- **No barcode-level author annotation file exists** on GEO or in the paper supplements. Labels were reconstructed from the authors' published scheme (Fig. 1B / Methods): major lineages by canonical markers; malignant = epithelial cells that are not a clear alveolar/club/ciliated program, using stromal cells (fibroblast + endothelial) as the CopyKAT-like CNV reference, as the authors did.

## Labels and sanity

| Lineage | Cells |
|---|---|
| T | 34,398 |
| Myeloid | 18,367 |
| Epithelial | 10,669 |
| Neutrophil | 10,053 |
| B | 7,604 |
| Plasma | 4,744 |
| NK | 4,023 |
| Fibroblast | 1,004 |
| Mast | 916 |
| Endothelial | 299 |
| Unassigned | 253 |
| **Malignant (of epithelial)** | **9,867** |
| **T/NK** | **38,421** |

Sanity: EPCAM+ 92.1% malignant vs 0.8% T/NK; PTPRC+ 5.2% malignant vs 93.9% T/NK. All 15 patients have ≥10 malignant cells (paper dropped one NMPR with <10; our stromal-referenced call is slightly more inclusive).

## 1. Malignant TACSTD2: NMPR vs MPR (post-treatment)

pCR (P06) grouped with MPR, matching the paper. Unit = patient.

| Metric | NMPR n=8 median | MPR n=4 median | Mann-Whitney p |
|---|---|---|---|
| mean log1p(CP10K) | 1.72 | 1.33 | 0.214 |
| % TACSTD2+ | 88.0 | 86.1 | 0.808 |
| pseudobulk CPM | 853 | 423 | 0.154 |

All-epithelial sensitivity (no CNV filter) is the same direction (p = 0.11 for mean log1p). Including the 3 pre-treatment biopsies does not change the conclusion (malignant mean p = 0.14).

## 2. Per-patient malignant TACSTD2 vs T/NK (claimed ρ −0.40 to −0.50)

| Contrast | n | Spearman ρ | p |
|---|---|---|---|
| **post malignant log1p vs T/NK fraction (primary)** | **12** | **−0.021** | **0.95** |
| post malignant %pos vs T/NK fraction | 12 | −0.238 | 0.46 |
| all malignant log1p vs T/NK fraction | 15 | −0.114 | 0.69 |
| post epithelial log1p vs T/NK fraction | 12 | −0.098 | 0.76 |
| post epithelial %pos vs T/NK fraction | 12 | −0.301 | 0.34 |
| all epithelial log1p vs T/NK fraction | 15 | −0.211 | 0.45 |
| post malignant log1p vs T/NK count | 12 | −0.161 | 0.62 |
| post residual tumor vs malignant log1p | 12 | +0.221 | 0.49 |

No contrast falls in −0.40 to −0.50. n=12 is too small to stabilize a ρ of that size; the point estimate is near zero.

## 3. TACSTD2 malignant vs T/NK (paired)

| Metric | Malignant median | T/NK median | Wilcoxon p |
|---|---|---|---|
| mean log1p(CP10K) | 1.63 | 0.011 | 6.1×10⁻⁵ |
| % positive | 86.5 | 0.69 | 6.1×10⁻⁵ |

P07 (NMPR, 53% malignant cells) is the only T/NK outlier (~16% TACSTD2+), consistent with ambient RNA from a tumor-rich sample, not true T/NK expression.

## Caveats

- Barcode labels are reconstructed from the paper's marker scheme, not the unpublished Seurat/CopyKAT object. Lineage sanity checks match expected markers; malignant counts are in the same ballpark as the paper's epithelial/CNV split.
- 4 vs 8 post-treatment patients. A non-significant NMPR>MPR shift should not be over-interpreted.
- P07 dominates malignant-cell counts (4,841 / 9,867). Patient-level tests treat each patient equally.

## Reproduce

```bash
python3 scripts/noskip/GSE207422/download.py --outdir data/GSE207422
python3 scripts/noskip/GSE207422/analyze.py \
  --matrix data/GSE207422/GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz \
  --metadata data/GSE207422/GSE207422_NSCLC_scRNAseq_metadata.xlsx \
  --gene-chr scripts/noskip/GSE207422/gene_chr.tsv \
  --outdir results/noskip/GSE207422
```
