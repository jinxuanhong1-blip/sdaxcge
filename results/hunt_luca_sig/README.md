# LuCA/Salcher/Leader TROP2 signature hunt

## Bottom line

For an `immune-low / TROP2-high` case, public processed Leader, Salcher, and
LuCA data support **TACSTD2 as an epithelial/malignant RNA signal**. They do
**not** support calling that case TRN-high or LCAM-high from TACSTD2 alone.

What can be stated from measured tables:

1. **TACSTD2 is not a member of the published immune signatures.** It is
   absent from Salcher TAN, NAN, TRN, and the generic-neutrophil signature,
   and from both sides of the Leader LCAM gene set.
2. **Leader's own processed tables treat TACSTD2 as epithelial, not LCAM.**
   Immune-versus-epithelial log2 fold-change is `-5.53` (immune
   `1.13e-5`, epithelial `5.71e-4`), in the same direction as `EPCAM`
   (`-7.83`) and `KRT8` (`-7.08`), opposite `CD3D` (`+4.29`) and `CXCL13`
   (`+2.93`). Leader LCAM construction drops genes with this ratio `<= -1`.
   In Leader LCAM-hi versus LCAM-lo pseudobulk DE, TACSTD2 `log2FC = 0.327`,
   `p = 0.251`, `adj. p = 0.359`.
3. **LuCA single-cell TACSTD2 is compartmentalized.** Detected in `55.7%` of
   malignant cells versus `4.14%` of neutrophils. Mean Census-normalized
   expression was `15.6`-fold higher in malignant cells. In BD Rhapsody
   alone, detection was `87.1%` malignant versus `4.85%` neutrophil.
4. **Bulk TCGA does not make TACSTD2 a proxy for LCAM.** On `992` primary
   LUAD/LUSC tumors, Spearman `TACSTD2` vs LCAM was `rho = -0.041`,
   `p = 0.196`.
5. **Bulk TCGA shows only a modest TACSTD2-TRN association, and it collapses
   in the immune-low TROP2-high slice.** All primary tumors:
   `rho = 0.221`, `p = 2.01e-12`, `n = 992`. The
   `immune-low and TROP2-high` slice (`n = 63`) was
   `rho = -0.155`, `p = 0.225` versus TRN and `rho = 0.212`, `p = 0.095`
   versus LCAM. That slice is too small and too mixed to assign this
   patient a TRN or LCAM state.

This patient's TRN and LCAM remain unmeasured. No patient matrix was
supplied.

## Leader processed result

Source files from
[`effiken/Leader_et_al`](https://github.com/effiken/Leader_et_al)
revision
[`4a88416`](https://github.com/effiken/Leader_et_al/commit/4a884161d50ed768963603c8a0aea38ea4c9299b):

- `input_tables/immune_vs_ep_de.csv`
- `input_tables/DE_LCAMhi_vs_LCAMlo_pseudobulk.rd`
- `scripts/get_LCAM_scores.R`
- `scripts/figure_6ac_s6abcdefjklmnopqr_7abcd_s7ab.R`

`immune_vs_ep_de.csv` columns are `l2fc`, `fg_exprs`, `bg_exprs`. The last
two are expression values, not p-values. Script polarity is immune over
epithelial: genes with `l2fc > -1` are kept when LCAM genes are selected.
TACSTD2 would be excluded by that filter.

TACSTD2 appears in the Figure 2A DC heatmap gene list next to `CLEC9A`.
That is a plotted gene, not an LCAM signature gene.

## Salcher / LuCA single-cell result

Final Salcher signatures:

- TAN: 18 genes
- NAN: 20 genes
- TRN: their 38-gene union
- Separate 17-gene major-cell-type `Neutrophils` signature in Table S5

TRN and the generic-neutrophil set overlap only at `CYP4F3` and `MGAM`.
TAN and NAN are disjoint.

LuCA analysis code uses `TACSTD2` with `AGR2` as an `undifferentiated`
epithelial/cancer annotation marker. That is not evidence that TACSTD2
tracks TRNs.

Census 2025-11-08, LuCA extended atlas
`1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1`, TACSTD2 only:

| cell type | n | detected fraction | mean expression |
| --- | ---: | ---: | ---: |
| malignant cell | 91528 | 0.557 | 0.000562 |
| epithelial cell of lung | 28763 | 0.551 | 0.000394 |
| neutrophil | 19368 | 0.041 | 0.000036 |
| plasma cell | 35538 | 0.030 | 0.000007 |
| CD8 T cell | 160394 | 0.012 | 0.000006 |

Malignant-cell detection was `67.1%` in LUAD annotations and `54.6%` in
LUSC annotations. Residual neutrophil TACSTD2 is low, assay-dependent, and
compatible with ambient RNA or doublets. No cell-level p-value is reported:
cells are not independent replicates.

The standardized Census object lacks LuCA's custom sample field, so tumor
versus normal-adjacent donor pseudobulk was not computed.

## TCGA bulk result

Public processed input: UCSC Xena TOIL `tcga_RSEM_gene_tpm`
(log2(TPM+0.001)), primary tumors (`-01`) labeled LUAD or LUSC in Thorsson
et al. All queried signature genes were present (`38/38` TRN, `63/63` LCAM
subtype genes).

Scoring:

- LCAM: Leader method on unlogged TPM (column fractions, `log10(1e-6+x)`,
  gene z-score, subtype mean, subtype z-score, `LCAM-hi - LCAM-lo`)
- TRN/TAN/NAN/generic neutrophils: mean gene z-score of `log1p(TPM)`
- Immune content: Thorsson leukocyte fraction
- Quartiles: computed on this `992`-sample set
- Correlations: Spearman, two-sided, from `scipy.stats.spearmanr`

| subset | n | TACSTD2 vs LCAM | TACSTD2 vs TRN | LCAM vs TRN |
| --- | ---: | --- | --- | --- |
| all primary | 992 | -0.041 (p=0.196) | 0.221 (p=2.01e-12) | -0.195 (p=5.96e-10) |
| LUAD | 508 | -0.085 (p=0.056) | 0.335 (p=8.27e-15) | -0.141 (p=0.0015) |
| LUSC | 484 | -0.137 (p=0.0026) | 0.136 (p=0.0028) | -0.249 (p=3.01e-08) |
| immune-low Q1 | 248 | 0.013 (p=0.835) | 0.192 (p=0.0024) | -0.025 (p=0.693) |
| TROP2-high Q4 | 248 | 0.167 (p=0.0084) | -0.084 (p=0.187) | -0.196 (p=0.0020) |
| immune-low and TROP2-high | 63 | 0.212 (p=0.095) | -0.155 (p=0.225) | -0.182 (p=0.154) |

TACSTD2 vs leukocyte fraction: `rho = -0.029`, `p = 0.370`, `n = 989`.
TROP2-high is not the same as immune-low.

The modest all-sample TACSTD2-TRN correlation is a bulk mixture
association. It does not make TACSTD2 a neutrophil gene, and it does not
survive as a significant association inside the immune-low TROP2-high
slice.

LCAM vs TRN `rho = -0.195` is a measured weak-to-modest anticorrelation on
TCGA. The Salcher preprint called TRN and LCAM poorly correlated; that
comparison is absent from the final Cancer Cell paper. No OAK/POPLAR
reanalysis was possible: those processed matrices are not public.

## What was not analyzed

Full LuCA H5ADs (`12.6-17.6 GB`) and Zenodo archives (`>=5.0 GB`) were not
downloaded. Protected OAK/POPLAR ICI matrices were not used. No p-value or
correlation was copied from memory or filled in when a test was not run.

A valid patient-specific TRN/LCAM call still needs that patient's tumor
expression matrix.

## Files

- `signatures.tsv`: published genes and LCAM subtype structure
- `signature_overlap.tsv`: set intersections
- `leader_immune_vs_epithelial.tsv`: extracted Leader processed rows
- `leader_lcam_de_tacstd2.tsv`: extracted Leader LCAM-hi vs lo DE row
- `tacstd2_relevant_cell_types.tsv`: LuCA Census one-gene summaries
- `tacstd2_malignant_vs_neutrophil_by_assay.tsv`: assay check
- `tacstd2_by_disease.tsv`: disease-stratified Census summaries
- `tcga_primary_nsclc_scores.tsv`: per-tumor scores
- `tcga_signature_coverage.tsv`: genes recovered from TOIL
- `tcga_spearman.tsv`: measured Spearman tests
- `tcga_quartile_strata.tsv`: leukocyte by TACSTD2 quartile means
- `tcga_stratum_counts.tsv`: sample counts and cutoffs
- `reproduce.py`: Census query
- `public_processed.py`: Leader extraction and TCGA scoring

```bash
python results/hunt_luca_sig/reproduce.py
python results/hunt_luca_sig/public_processed.py --toil /path/to/tcga_RSEM_gene_tpm.gz
```

`public_processed.py` expects local copies of the Leader `input_tables`
files listed above and the Xena TOIL matrix. Those large matrices are not
committed.

## Provenance

- Salcher et al., Cancer Cell 2022:
  [doi:10.1016/j.ccell.2022.10.008](https://doi.org/10.1016/j.ccell.2022.10.008)
- Salcher [Table S4](https://ars.els-cdn.com/content/image/1-s2.0-S1535610822004998-mmc1.pdf)
  and [Table S5](https://ars.els-cdn.com/content/image/1-s2.0-S1535610822004998-mmc3.xlsx)
- LuCA source
  [`a3c0184`](https://github.com/icbi-lab/luca/commit/a3c018451850fbb50b674bb606ea7348182420cb)
- LuCA extended atlas, CELLxGENE Census `2025-11-08`, dataset
  `1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1`
- Leader et al., Cancer Cell 2021:
  [doi:10.1016/j.ccell.2021.10.009](https://doi.org/10.1016/j.ccell.2021.10.009)
- Leader scoring source
  [`get_LCAM_scores.R`](https://github.com/effiken/Leader_et_al/blob/4a884161d50ed768963603c8a0aea38ea4c9299b/scripts/get_LCAM_scores.R)
- UCSC Xena TOIL `tcga_RSEM_gene_tpm`
- Thorsson et al. TCGA immune metrics, as bundled by Leader
