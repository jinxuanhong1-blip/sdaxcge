# Additive Stereo-seq / BGI spatial, CLDN4-only

Human lung cancer Stereo-seq only. No Visium fallback. No private 8-KL. No TACSTD2. No claim-failed language.

## Spatial unit

GSE328481 deposits author cell-bin AnnData (`uns/bin_type = cell_bins`, `resolution = 500` nm). Analyses are run at:

1. **cell-bin** (native)
2. **bin50** (sum counts / first `bin50_x,y` per `bin50_location_id`)

QC: cell-bin ≥50 genes and ≥20 UMI; bin50 ≥100 genes and ≥50 UMI.

Expression for correlations is log1p(CP10K).

## Epithelial / CD8

- Epithelial: author `anno` in `{cancer_cell, epi}`. Bin50 is epithelial if ≥50% of cells in that bin have those labels.
- CD8A+ unit: raw `CD8A` count > 0.
- Nearest-CD8 reference: CD8A+ and **not** epithelial (mixed bins cannot self-match).
- CLDN4-high epithelial: epithelial with CLDN4>0; if ≥40 detected, keep the top half of those log1p values.
- CLDN4-low epithelial: epithelial with CLDN4==0.

Distance is Euclidean on `(x, y) * 0.5` µm.

## Statistics (locked)

1. Spearman(CLDN4, CD8A) on all QC units, and again inside epithelial units.
2. KRT8 residual: OLS residual of log1p(CLDN4) ~ log1p(KRT8), then Spearman vs CD8A.
3. Median nearest CD8A+ distance for CLDN4-high vs CLDN4-low epithelial units; Mann–Whitney per sample.
4. Sample is the unit of inference (n = 11). Report median ρ / median ∆µm and two-sided Wilcoxon signed-rank on the 11 values.

A nonzero-union Spearman (units with CLDN4>0 or CD8A>0) is recorded as a **detection-exclusivity** check, not as an expression-level neighborhood claim.

Reproduce:

```bash
bash scripts/stereo_seq_luad/download_gse328481.sh data/GSE328481
GSE328481_DIR=data/GSE328481 python3 scripts/stereo_seq_luad/analyze_gse328481.py
python3 scripts/stereo_seq_luad/make_figures.py
```
