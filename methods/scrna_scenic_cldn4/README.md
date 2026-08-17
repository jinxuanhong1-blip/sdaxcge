# methods/scrna_scenic_cldn4 — CLDN4-high malignant regulons (GSE207422)

Additive public scRNA GRN. **A10 ELF3–CLDN4 bulk RNA is taken as given** and is not re-tested as a discovery. [`methods/scrna_scenic`](../scrna_scenic/) asked whether ELF3/GRHL1/KLF4 AUCell tracks TACSTD2/CLDN4. This folder asks a different question: **which TF regulons are recovered from CLDN4-high malignant-like cells?**

**CLDN4-high regulons only.** CLDN4-low edges are used only to mark high-specific targets and are not published as a regulon table.

**No ChIP peaks are invented.** Full pySCENIC cisTarget is not run.

Finding: [`FINDING.md`](FINDING.md). Tables: [`results/scrna_scenic_cldn4/`](../../results/scrna_scenic_cldn4/).

## Question (pre-specified)

1. Gate malignant-like epithelium on public **GSE207422** UMI (same marker proxy as `scrna_scenic`; CopyKAT IDs are not on GEO).
2. Split **within-sample** CLDN4-high vs low (median `log1p(CP10k)` among malignant-like; sample kept if ≥20 such cells).
3. Infer a lightweight GRN **in CLDN4-high cells only** (Pearson + public priors).
4. Report the **regulon table**. Honest n = samples / cells as labeled.

## Run

```bash
python3 methods/scrna_scenic_cldn4/scripts/download_gse207422.py
python3 methods/scrna_scenic_cldn4/scripts/download_priors.py   # optional; snapshot is committed
python3 methods/scrna_scenic_cldn4/scripts/analyze.py
```
