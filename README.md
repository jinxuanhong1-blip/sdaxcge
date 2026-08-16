# Hunt: CLDN / tight-junction family vs TACSTD2 (Trop-2) in public lung RNA

A small, fully reproducible co-expression analysis. It asks a focused question:

> In human lung, is the tight-junction / apical-junction "CLDN family" —
> **CLDN1, CLDN7, F11R (JAM-A), PARD3, and CLDN4** — co-expressed with
> **TACSTD2** (Trop-2, the target of the antibody-drug conjugate sacituzumab
> govitecan)?

To avoid a single-dataset fluke, it uses **three independent public per-sample
lung RNA-seq datasets** and takes the **triple intersection** of genes robustly
positively co-expressed with TACSTD2 in all three.

## Datasets (all public, re-downloadable)

| key | source | tissue | unit |
|---|---|---|---|
| `gtex_lung` | GTEx v8 (bulk-gex, tpms-by-tissue) | normal lung | raw TPM |
| `tcga_luad` | UCSC Xena / GDC hub, STAR TPM | lung adenocarcinoma | log2(TPM+1) |
| `tcga_lusc` | UCSC Xena / GDC hub, STAR TPM | lung squamous cell | log2(TPM+1) |

Exact URLs, byte sizes and SHA-256 hashes are recorded in `data/provenance.json`
after download. Raw data is **not** committed (see `.gitignore`); it is fetched
by `src/download_data.py`.

## Method (short)

Genome-wide **Spearman** rank correlation of TACSTD2 against every gene, computed
independently per dataset. Spearman is monotone-invariant, so the raw-TPM vs
log2-TPM difference does not matter and the three datasets are directly
comparable. Genes are matched across datasets by base Ensembl gene id. The
"triple-intersect" is the set of genes with BH-FDR < 0.05 **and** rho >= 0.30 in
**all three** datasets. See `results/hunt_cldn_family/REPORT.md` for full methods.

## Headline result (honest)

All five panel genes are positively and FDR-significantly co-expressed with
TACSTD2 in every dataset, but at the strict `rho >= 0.30`-in-all-three bar only
**CLDN4** clears it (mean rho ≈ 0.57; it is the #2 gene genome-wide in TCGA-LUAD).
CLDN1, CLDN7 and F11R are borderline (they clear 0.30 in two of three datasets);
PARD3 is the weakest, particularly in the tumour datasets.

Importantly, the TACSTD2 co-expression core is dominated by epithelial genes
(KRT19, CDH1, SDC1, PERP, EVPL, MAL2, NECTIN4, ...), so these positive
correlations largely reflect a shared **epithelial-content** program in bulk
tissue rather than a TACSTD2-specific mechanism. Read the caveats in the report.

## Reproduce

```bash
pip install -r requirements.txt
python src/download_data.py   # ~380 MB -> data/ (git-ignored)
python src/analyze.py         # -> results/hunt_cldn_family/
```

## Outputs (`results/hunt_cldn_family/`)

- `REPORT.md` — full write-up, generated from the numbers (nothing hand-edited)
- `panel_summary.csv` — the 5 panel genes + TACSTD2, per dataset (rho, FDR, rank)
- `triple_intersection.csv` — the robust positive core common to all 3 datasets
- `set_sizes.csv` — robust-set sizes and pairwise/triple overlaps
- `summary.json` — machine-readable headline numbers
- `per_dataset_correlations/*.csv.gz` — full genome-wide ranked tables per dataset
- `fig_panel_rho.png`, `fig_gtex_scatter.png` — figures
