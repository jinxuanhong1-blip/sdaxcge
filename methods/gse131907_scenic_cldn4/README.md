# methods/gse131907_scenic_cldn4 — malignant CLDN4-high regulons (GSE131907)

Additive public scRNA GRN. **A10 ELF3–CLDN4 bulk RNA is taken as given** and is not re-tested as a discovery.

[`methods/scrna_scenic`](../scrna_scenic/) asked whether ELF3/GRHL1/KLF4 AUCell tracks TACSTD2/CLDN4 on GSE207422. [`methods/scrna_scenic_cldn4`](../scrna_scenic_cldn4/) inferred CLDN4-high regulons on GSE207422 (marker-proxy malignant-like). This folder asks the same regulon question on **GSE131907 author-labeled malignant cells**.

**CLDN4-high regulons only.** CLDN4-low edges are used only to mark high-specific targets and are not published as a regulon table.

**No ChIP peaks are invented.** Full pySCENIC cisTarget is not run.

Finding: [`FINDING.md`](FINDING.md). Tables: [`results/gse131907_scenic_cldn4/`](../../results/gse131907_scenic_cldn4/).

## Question (pre-specified)

1. Gate **author-malignant** epithelium on public **GSE131907** UMI (`Cell_subtype` in `{Malignant cells, tS1, tS2, tS3}`).
2. Split **within-sample** CLDN4-high vs low (median `log1p(CP10k)` among malignant; sample kept if ≥20 such cells).
3. Infer a lightweight GRN **in CLDN4-high cells only** (Pearson + public priors).
4. Report the **regulon table**. Honest n = samples / cells as labeled.

## Run

```bash
python3 methods/gse131907_scenic_cldn4/scripts/download_gse131907.py
python3 methods/gse131907_scenic_cldn4/scripts/download_priors.py   # optional; snapshot is committed
python3 methods/gse131907_scenic_cldn4/scripts/analyze.py
```
