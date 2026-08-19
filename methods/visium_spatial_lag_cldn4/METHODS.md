# Methods — Visium spatial lag (CLDN4-only)

Public Visium LUAD/NSCLC only. No private KL. Section is the unit.

## Why lag, not same-spot

A standard Visium spot is 55 µm and routinely mixes epithelium and T cells.
Same-spot Spearman(CLDN4, CD8A) and nearest-CD8A-high Euclidean µm therefore
collapse toward 0 even when high-CLDN4 epithelium excludes CD8 from *adjacent*
spots. Those readouts are not the exclusion test and are not reported as the punch.

## Input

Space Ranger filtered feature-barcode matrix plus `tissue_positions`
(`scalefactors_json.json` when present). H&E, BAM, and FASTQ were not kept.

Sources: 10x CytAssist lung cancer demos; GEO GSE189487, GSE273378, GSE300676,
GSE307534 (invasive + precursor; normal omitted).

## QC and normalization

Keep spots with ≥200 detected genes and ≥100 UMI. Convert genes of interest
to log1p(CPM to 10,000) using that spot’s UMI total.

## Epithelial-like

Mean of whichever of `KRT8`, `KRT18`, `KRT19`, `EPCAM` are present.
Epithelial-like = score ≥ section median.

## Neighbors

Primary: Visium hex ring-1 (self excluded; ~100 µm center-to-center).
Sensitivity: Euclidean annulus 80–150 µm.

## Contrast filter (drop, do not test)

Drop a section if any of:

- CD8A+ spots < 50 or < 5% of QC spots
- CLDN4+ among epithelial-like < 10%, or CLDN4 IQR among epi = 0
- < 80 epithelial-like spots with a ring-1 neighbor
- neighbor-CD8A SD among those spots = 0

## Tests (each kept section)

1. **Lag-ρ**: Spearman(index CLDN4, mean ring-1 CD8A) on epithelial-like spots.
   95% CI by Fisher z.
2. **Q4 vs Q1 Δ**: among those spots, equal-sized bottom/top 25% of CLDN4;
   Δ = mean(neighbor CD8A | Q4) − mean(neighbor CD8A | Q1). Welch 95% CI.
   Negative Δ = CLDN4-high epithelium has less CD8A next door.
3. **Morisita–Horn** overlap on a 200 µm grid between CLDN4-Q4 epi spots and
   CD8A-high spots (CD8A ≥ q75 if q75>0, else CD8A>0). 0 = segregated, 1 = mixed.
4. **KRT8 residual on the lag**: OLS residual of index CLDN4 ~ KRT8, then
   Spearman vs neighbor CD8A (not same-spot CD8A).
5. **Optional interface**: restrict index spots to epi-like spots with ≥1
   non-epi ring-1 neighbor; repeat lag-ρ and Q4−Q1 Δ.

## Cross-section

Wilcoxon signed-rank vs 0 on the vector of section lag-ρ and on the vector
of section Δ. Report n, n_negative / n_positive.

## Distance scale

If `spot_diameter_fullres` is present: µm/pixel = 55 / diameter.
Otherwise calibrate so median nearest-neighbor pixel distance = 100 µm.
