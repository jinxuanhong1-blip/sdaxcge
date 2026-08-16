# B6 leftover — public Visium lung: CLDN4-high vs immune niches

**Output dir:** `results/w200/B6_visium_extra/`
**Script:** `scripts/b6_visium_extra_cldn4_immune.py`
**Fetch (small files only):** `scripts/fetch_b6_visium_extra.sh`

## Dataset (not E-MTAB-13530)

Public **10x Genomics Visium CytAssist** section:

- Name: Human Lung Cancer (FFPE), squamous cell carcinoma
- Page: https://www.10xgenomics.com/datasets/human-lung-cancer-ffpe-2-standard
- Space Ranger 2.0.0 outputs, license **CC BY 4.0**
- Published metrics: 3,858 spots under tissue; median 6,174 genes / 13,991 UMI per spot
- **Not** ArrayExpress E-MTAB-13530 (that accession is a different lung spatial study)

Huge raw files were **not** downloaded. See [SKIPPED_FILES.md](SKIPPED_FILES.md).

## Question

Do **CLDN4-high** epithelial/tumor spots occupy the same space as **immune niches**,
or are they spatially exclusive?

## Methods (short)

1. Load filtered H5 + spatial folder (`sc.read_visium`).
2. QC: drop spots with <500 counts or <250 genes; genes in <5 spots.
   Kept **3,816 / 3,858** spots, 18,053 genes.
3. Normalize (10k), log1p, HVG, PCA, neighbors, Leiden (res=1.0).
4. Score spots:
   - **CLDN4** (log-norm expression)
   - **Pan-immune** (`PTPRC` + T / B-plasma / myeloid / NK markers)
   - Lineage scores and a small epithelial set (`EPCAM`, `KRT5/6A/17`, `SFTPC/B`, `NAPSA`)
5. **CLDN4-high:** top quartile of CLDN4.
6. **Immune niche:** top-quartile immune score **and** Leiden cluster whose
   mean immune score is above the 75th percentile of cluster means
   (clusters 5, 7, 13, 16).
7. Tests: Spearman, Fisher exact, Mann–Whitney, squidpy neighborhood
   enrichment (6-neighbor grid), nearest-immune distance, Wilcoxon DE.

## Niche counts

| Niche | Spots |
|---|---|
| CLDN4-high | 859 |
| Immune niche | 402 |
| CLDN4-high & immune (overlap) | 95 |
| Other | 2,460 |

Overlap is small: 95 / 954 CLDN4-high spots also sit in an immune niche.

## Findings — CLDN4 and immune niches exclude each other

1. **Anti-correlation.** CLDN4 vs pan-immune score: Spearman **ρ = −0.25**,
   p = 5.6×10⁻⁵⁶.
2. **Co-membership is depleted.** Fisher exact on CLDN4-high × immune-niche:
   **OR = 0.68**, p = 1.0×10⁻³ (2,460 / 402 / 859 / 95).
3. **Immune niches are CLDN4-low.** CLDN4 in immune-niche vs other spots:
   Mann–Whitney p = 4.1×10⁻²² (immune-niche median CLDN4 = 0).
4. **Immediate neighbors avoid each other.** Neighborhood enrichment
   z(CLDN4-high, Immune niche) = **−14**. Self-enrichment is strong for
   immune niches (z = 51) and present for CLDN4-high (z = 9).
5. **Tissue-scale distance is not farther than background.** Median distance
   to the nearest immune-niche spot is 283 µm (CLDN4-high) vs 284 µm (other),
   MWU p = 0.35. Exclusion is **local** (they are not 6-neighbor adjacent),
   not a global “opposite sides of the section” effect.
6. **DE matches the labels.** CLDN4-high vs immune niche: up `CLDN4`, `CES1`,
   `FXYD3`, `CD9`, `PERP`, `KRT8`, `SOX2`, `KRT19`. Immune niche: `IGKC`,
   `IGHG1`, `IGHA1`, `JCHAIN`, `MZB1` (plasma), plus stroma (`VIM`, `DCN`, `MMP2`).
7. **Lineage scores.** Immune niches are B/plasma- and myeloid-high; CLDN4-high
   spots are epithelial-high and immune-low. The 95 mixed spots look like
   immune-niche edges that still express CLDN4.

## Figures

- `figures/01_spatial_overview.png` — CLDN4, immune score, niche map, Leiden
- `figures/02_cldn4_vs_immune_stats.png` — scatter, violins, signature heatmap
- `figures/03_spatial_relationship.png` — nhood z-scores and distance histograms

## Tables

- `tables/spot_annotations.csv`
- `tables/contingency_cldn4high_x_immuneniche.csv`
- `tables/signature_means_by_niche.csv`
- `tables/nhood_enrichment_zscores.csv`
- `tables/DE_CLDN4high_vs_immune_niche.csv`
- `stats_summary.json`

## Reproduce

```bash
python3 -m venv .venv && .venv/bin/pip install scanpy squidpy leidenalg igraph matplotlib pandas numpy scipy anndata
bash scripts/fetch_b6_visium_extra.sh
.venv/bin/python scripts/b6_visium_extra_cldn4_immune.py
```

Data stay under `data/` (gitignored). Results are committed here.
