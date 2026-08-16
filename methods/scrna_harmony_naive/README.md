# Methods · naive LUAD/NSCLC Harmony joint object

Additive public scRNA only. Extra atlas *n* for epithelial / malignant-like
TACSTD2 and CLDN4 versus per-donor T/NK and B/TLS-like fractions. **Not an
ICI, MPR, or response test.**

## Inclusion

| Accession | Paper | What was kept | Both genes |
|---|---|---|---|
| GSE131907 | Kim *Nat Commun* 2020 | Tumor sites only: tLung, tL/B, mLN, PE, mBrain | yes (streamed UMI) |
| GSE253013 | Xiang/Sze *Cancer Res* 2024 | 9 treatment-naïve LUAD; tumor lanes | yes (GEO RDS) |
| GSE148071 | Wu *Nat Commun* 2021 | 42 advanced NSCLC diagnostic biopsies | checked on extract |
| GSE127465 | Zilionis *Immunity* 2019 | Tumor (not blood); epithelium exists | yes |
| GSE154826 | Leader *Cancer Cell* 2021 | **Omitted** | — |

Leader/GSE154826 is CD45-bead CITE-seq. Epithelium sits in the authors’
`epi_endo_fibro_doublet` gate (two digest patients without CD45 enrichment).
That is not a whole-tumor epithelial compartment, so it is not in the joint
object.

GSE148071 GEO has no prior-therapy field. The series is advanced diagnostic
biopsies (Wu 2021); it is included as an open tumor atlas with both genes,
not as an ICI cohort. Zilionis: 6/7 patients untreated (Table S1); the treated
ID is not in the GEO metadata, so all seven tumor donors are kept and the
limitation is stated.

## Joint object

1. Stream a shared lineage-plus-target gene panel from each public matrix
   (full transcriptomes are not held in RAM).
2. Restrict to tumor / tumor-site cells.
3. Assign a **joint** lineage by argmax of mean log1p marker scores
   (epithelial, T, NK, B, myeloid, fibroblast, endothelial).
4. Malignant-like = epithelial AND near-zero normal-lung markers
   (`SFTPA2`, `SFTPC`, `AGER`, `SCGB1A1`, `SCGB3A1`, `TPPP3`, `FOXJ1`).
   This is a marker proxy, not CopyKAT.
5. T/NK = joint T + NK. B/TLS-like = B/plasma, or non-epithelial cells with
   high CXCL13/CCL19/CCL21. This is a dissociated co-occurrence proxy, not
   a scored follicle.
6. Harmony (`harmonypy`) on dataset-wise z-scored log1p panel genes → PCA
   (≤20 PCs) → Harmony batch = `dataset`. Embedding uses ≤2,500 cells per
   donor. Scores use **all** tumor cells, not the subsample.

## Statistics

Unit = tumor donor (Kim tumor **site**, Xiang/Wu/Zilionis **patient**).
Eligible: ≥20 epithelial and ≥20 T/NK cells (B/TLS contrasts also require
≥10 B cells). Tests are Spearman ρ and two-sided *p* at that unit.
Dataset-stratified and dataset-residualized rank Spearmans are reported
alongside the pooled joint numbers. Primary joint contrasts are BH-adjusted
together. Cell-level correlations are not reported.

## Reproduce

```bash
python3 methods/scrna_harmony_naive/download.py
python3 methods/scrna_harmony_naive/extract.py
python3 methods/scrna_harmony_naive/extract_gse253013.py
python3 methods/scrna_harmony_naive/assemble_gse253013.py
python3 methods/scrna_harmony_naive/integrate.py
```

Outputs: `results/scrna_harmony_naive/`.

## Honest joint numbers

Filled after the public-matrix run (see `results/scrna_harmony_naive/summary.json`).
