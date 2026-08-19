# Methods

Pre-specified. CLDN4-only. Section is the unit.

## Input

Space Ranger filtered feature-barcode matrix plus `tissue_positions` (and `scalefactors_json.json` when present). Huge H&E TIFFs, BAM, FASTQ, and Loupe files were not downloaded.

## QC and normalization

Keep spots with ≥200 detected genes and ≥100 UMI. Convert the gene of interest to log1p(CPM to 10,000) using that spot’s UMI total.

## Epithelial-like spots

Mean of whichever of `KRT8`, `KRT18`, `KRT19`, `EPCAM` are present in the matrix. Epithelial-like = score ≥ section median.

## CD8A-high

Spots with CD8A ≥ the 75th percentile of all QC spots when that cutoff is >0; otherwise CD8A > 0 (sparse CD8A). Distance to the nearest CD8A-high spot excludes self.

## Tests (each section)

1. Spearman correlation of CLDN4 vs CD8A (all QC spots, and epithelial-like only).
2. Among epithelial-like spots, CLDN4 Q4 vs Q1:
   - Euclidean distance (µm) to the nearest CD8A-high spot.
   - Hex ring-1 mean neighbor CD8A.
   Mann–Whitney U, two-sided.
3. KRT8 residual: OLS residual of CLDN4 on KRT8 (same log1p-CPM scale). Spearman of residual vs CD8A; Q4 vs Q1 of the residual with the same distance and neighbor tests.

## Distance scale

If `spot_diameter_fullres` is present: µm/pixel = 55 / spot_diameter_fullres. Otherwise calibrate so the median nearest-neighbor pixel distance equals 100 µm (Visium center-to-center).

## Not done

No TACSTD2 dual-high filter. No private 8-KL. No deconvolution. No claim-failed framing.
