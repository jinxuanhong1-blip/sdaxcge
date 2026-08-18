# FINDING — Seurat RPCA win-pair GSE131907 + GSE205335, CLDN4-only

ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. Not a bigger merge.
Winning pair is given (PR #320: author-malignant CLDN4 %pos vs T/NK,
Spearman ρ=−0.479, N=43; Q4 vs Q1 r=−0.705, n=23). Thesis is already
correct and is **not** re-derived: CLDN4-high malignant cells sit with
lower T/NK and lower own IFN/MHC. This folder adds a **Seurat v5
`IntegrateLayers` (RPCA)** object of the same two public GEO sets.

Primary engine: **R 4.6.1 + Seurat 5.5.1** (`results/sessionInfo.txt`).
Python only streams the GSE131907 UMI TSV. No Python-only primary. No
GSE148071. GSE207422 is not added.

Honest unit: **GSE131907 GEO Sample** (locked PR #320) and **GSE205335
patient**. Cells are counts, not replicates. The 21 gated GSE131907
samples collapse to **17** unique patient numbers (EBUS+NS share 6, 12,
13, 19). Patient-collapsed n=39 is a sensitivity, not a replacement of
the locked n=43.

## Design (locked)

| Piece | Choice | Honest limitation |
|---|---|---|
| Cohorts | GSE131907 + GSE205335 only | Not 148071 / 207422 / 189357 / 123902 |
| Integration | Seurat v5 `IntegrateLayers(RPCAIntegration)` on dataset layers | UMAP only; scores use joined RNA log-norm |
| CLDN4 | Single gene; %pos = counts > 0 in author-malignant cells | ≤120 malignant cells / unit |
| T/NK fraction | Locked author full-sample `frac_tnk` (PR #320) | Not recomputed from the capped object |
| IFN / MHC | `AddModuleScore`; Hallmark IFNα∪IFNγ (222 genes); custom MHC-I + MHC-II (29) | CLDN4 excluded from both |
| Quartiles | Within-cohort; Seurat %pos and locked %pos both reported | Seurat Q4 on GSE131907 is 4 (ties), locked Q4 is 5 |
| Model | Within-cohort Spearman; DerSimonian–Laird on Fisher-z | Two cohorts; p-values are descriptive |

## Honest n

| Item | n | Note |
|---|---:|---|
| GSE131907 gated samples (n_mal≥20) | **21** | GEO Sample = locked unit |
| GSE131907 unique patient numbers | **17** | EBUS+NS pairs for 6/12/13/19 |
| GSE205335 patients | **22** | tumor GSMs collapsed; P4001 n_mal=27 kept |
| Combined locked units | **43** | 21 + 22; cells are not n |
| Patient-collapsed units | **39** | 17 + 22 |
| Cells in integrated object | 10045 | cap ≤120 mal + ≤120 T/NK / unit |
| Author-malignant / T/NK in object | 4983 / 5062 | labels only |

Seurat CLDN4 %pos tracks locked %pos (n=43, ρ=0.974, p=4.8e-28).

## CLDN4 vs T/NK (holds)

T/NK = locked full-sample fraction. CLDN4 = Seurat malignant %pos unless
said otherwise. Combined = DerSimonian–Laird on Fisher-z of the two
cohort Spearmans.

| test | n | ρ | p | I² | stacked Q4 vs Q1 |
|---|---:|---:|---:|---:|---|
| Seurat %pos vs T/NK (DL) | 43 | **−0.452** | 0.00306 | 0% | r=−0.467, 12/10, p=0.0698 (Seurat Q) |
| Locked %pos vs T/NK (DL) | 43 | **−0.479** | 0.00152 | 0% | **r=−0.697, 12/11, p=0.00511** (locked Q) |
| Patient-collapse Seurat %pos vs T/NK | 39 | −0.528 | 0.00489 | 30% | — |
| Patient-collapse locked %pos vs T/NK | 39 | −0.613 | 0.00690 | 56% | — |

Within-cohort Seurat %pos vs T/NK: GSE131907 n=21 ρ=−0.522 p=0.0153;
GSE205335 n=22 ρ=−0.380 p=0.0813.

Locked %pos vs T/NK recovers PR #320 (ρ=−0.479, N=43). Locked Q4 vs Q1
recovers the given r=−0.705 on n=23 (here r=−0.697, n_Q1/n_Q4=12/11).
Seurat Q4 on GSE131907 is 4 because of ties in the capped %pos — do not
quote that tail as the PR #320 n=5.

## Malignant IFN / MHC (same-sign, thin on this object)

Scores = mean `AddModuleScore` in capped author-malignant cells.
Positive ρ = CLDN4-high with higher IFN/MHC. Thesis direction is
**negative**.

| test | n | ρ | p | I² |
|---|---:|---:|---:|---:|
| Seurat %pos vs IFN (DL) | 43 | −0.244 | 0.129 | 0% |
| Seurat %pos vs MHC (DL) | 43 | −0.139 | 0.534 | 46% |
| Locked %pos vs Seurat IFN (DL) | 43 | −0.309 | 0.0886 | 23% |
| Locked %pos vs Seurat MHC (DL) | 43 | −0.205 | 0.497 | 71% |

GSE205335 alone is the cleaner IFN/MHC arm (locked %pos vs IFN ρ=−0.465
p=0.0293; vs MHC ρ=−0.470 p=0.0272). GSE131907 IFN/MHC on the 120-cell
cap is null (IFN ρ=−0.102; MHC ρ=+0.088). Combined MHC I²=46–71% — do
not quote a pooled MHC win from this Seurat object.

This does **not** reopen the thesis. Genome-wide patient-pseudobulk DE
on the same pair (muscat-folder extra) already has malignant IFN/MHC
down in CLDN4-high. Here the new method is Seurat RPCA + module scores
on a memory cap. IFN/MHC is same-sign and under-powered, not a flip.

## What this is not

- Not a dual-high TACSTD2×CLDN4 gate.
- Not GSE148071 or any third cohort.
- Not a cell-level Wilcoxon sold as n.
- Not evidence that CLDN4 *causes* T/NK loss or IFN/MHC drop.
- Not a Python Scanpy primary. Seurat installed and ran.
- Not a claim that MHC is significant on this capped object.

## Files

- `scripts/03_seurat_integrate.R` — merge, RPCA, module scores, unit tests
- `results/sessionInfo.txt` — R 4.6.1, Seurat 5.5.1
- `results/tables/patient_scores.tsv` — 43 units
- `results/tables/patient_level_tests.tsv`
- `results/tables/extra_unit_tests.tsv` — patient collapse + locked Q
- `results/tables/n_honest.tsv`
- `results/figures/` — UMAP (dataset / compartment / CLDN4 / IFN), unit scatters, Q4 vs Q1 boxes

## Reproduce

```bash
bash methods/seurat_winpair_131907_205335_cldn4/scripts/00_download.sh
Rscript methods/seurat_winpair_131907_205335_cldn4/scripts/01_keep_cells.R
python3 methods/seurat_winpair_131907_205335_cldn4/scripts/02_extract_gse131907.py
Rscript methods/seurat_winpair_131907_205335_cldn4/scripts/03_seurat_integrate.R
```
