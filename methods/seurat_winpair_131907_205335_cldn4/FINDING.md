# FINDING — Seurat RPCA win-pair GSE131907 + GSE205335, CLDN4-only

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. Not a bigger merge.
Winning pair is given (PR #320: author-malignant CLDN4 %pos vs T/NK,
Spearman ρ=−0.479, N=43; Q4 vs Q1 r=−0.705, n=23). Thesis is already
correct and is **not** re-derived: CLDN4-high malignant cells sit with
lower T/NK and lower own IFN/MHC. This folder adds a **Seurat v5
`IntegrateLayers` (RPCA)** object of the same two public GEO sets.

Primary engine: **R + Seurat** (see `results/sessionInfo.txt`). Python is
used only to stream the GSE131907 UMI TSV. No Python-only primary. No
GSE148071. GSE207422 is not added.

Honest unit = **GSE131907 GEO Sample** + **GSE205335 patient**. Cells are
counts, not replicates. On the PR #320 tumor gate the 21 GSE131907 samples
have unique patient numbers, so sample n = patient n = 21 there; the unit
is still named Sample because that is the GEO label.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE131907 + GSE205335 only | Not 148071 / 207422 / 189357 / 123902 |
| Integration | Seurat v5 `IntegrateLayers` RPCA on `dataset` layers | Integration is for UMAP; scores use joined RNA log-norm |
| CLDN4 | Single gene; %pos = counts > 0 in author-malignant cells | Capped ≤120 malignant cells / unit |
| T/NK fraction | Locked author full-sample `frac_tnk` (PR #320) | Not recomputed from the 50/50-capped object |
| IFN / MHC | `AddModuleScore`; Hallmark IFNα∪IFNγ; custom MHC-I + MHC-II | CLDN4 excluded from both modules |
| Quartiles | Within-cohort Seurat CLDN4 %pos | Mid quartiles unused in Q4 vs Q1 |
| Model | Within-cohort Spearman; DerSimonian–Laird on Fisher-z; stacked MW | Two cohorts; p-values are descriptive |

## Honest n

_Filled after the Seurat run. See `results/tables/n_honest.tsv`._

## Results

_Pending `03_seurat_integrate.R`. Do not quote numbers from this stub._

## What this is not

- Not a dual-high TACSTD2×CLDN4 gate.
- Not GSE148071 or any third cohort.
- Not a cell-level Wilcoxon sold as n.
- Not evidence that CLDN4 *causes* T/NK loss or IFN/MHC drop.
- Not a Python Scanpy primary.

## Files

- `scripts/03_seurat_integrate.R` — Seurat merge, RPCA, module scores, unit tests
- `results/tables/patient_scores.tsv`
- `results/tables/patient_level_tests.tsv`
- `results/figures/`

## Reproduce

```bash
bash methods/seurat_winpair_131907_205335_cldn4/scripts/00_download.sh
Rscript methods/seurat_winpair_131907_205335_cldn4/scripts/01_keep_cells.R
python3 methods/seurat_winpair_131907_205335_cldn4/scripts/02_extract_gse131907.py
Rscript methods/seurat_winpair_131907_205335_cldn4/scripts/03_seurat_integrate.R
```
