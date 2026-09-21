# Squidpy graphs: CLDN4-high tumor vs CD8 and NK (He 2022 CosMx)

Additive layer on the CellCharter clustered CosMx NSCLC object. It does not replace the locked exclusion result.

Locked result, left as locked: around CLDN4-high tumor cells the cytotoxic neighborhood ratio is 0.36 at 50 µm and 0.52 at 100 µm, in 8/8 sections and 5/5 donors (one-sided sign P = 0.031). Effector transcripts in the cells that are nearby are not lower (GZMB / PRF1 / NKG7 / IFNG high/low 1.11–1.22, 0/8 sections decreased). The short name for that pair of facts is exclusion, not muzzling.

## Data

- Figshare [25976224.v2](https://doi.org/10.6084/m9.figshare.25976224.v2), file `cosmx_human_nsclc_clustered.h5ad` (figshare file id 46841842).
- MD5 `71a84bdcf3cc625431aded47bad9a109`.
- 765,771 cells × 960 genes. Eight sections, five donors: Lung5 (LUAD-5 R1/R2/R3), Lung6 (LUSC-6), Lung9 (LUAD-9 R1/R2), Lung12, Lung13.
- Author `cell_type` in the object. Tumor = `tumor 5/6/9/12/13` (302,313 cells). CD8 = `T CD8 memory` + `T CD8 naive` (15,955). NK = `NK` (7,814). Epithelial cells stay in `other`.
- `obsm['spatial']` is the release’s global pixel coordinate. Radii use 0.18 µm per pixel (CosMx SMI NSCLC image scale). One FOV (LUAD-5 R1, fov 6) spans 5,445 × 3,626 px, about 0.98 × 0.65 mm at that scale. The h5ad does not store the scale factor.
- CLDN4-high / low is a within-section median split of the object’s normalized `X[:, CLDN4]` among tumor cells (`>` median vs `<=` median). In LUSC-6 and LUAD-12 the tumor median is 0 (75.5% and 60.6% zeros), so “high” there means CLDN4 detected. Thresholds: `results/cosmx_squidpy_fov/tables/cldn4_thresholds.csv`.

## Methods

Squidpy 1.8.3, seed 20260921.

1. **Neighborhood enrichment.** `gr.spatial_neighbors_radius` at 50 µm and 100 µm, then `gr.nhood_enrichment` (1,000 label permutations). The z-score shuffles every label, so a large negative z for both tumor arms is tumor–immune geometry. The CLDN4 contrast is z(high, effector) − z(low, effector).
2. **Interaction matrix.** `gr.interaction_matrix` on the same radius graph. Reported number is the mean count of CD8 neighbors, NK neighbors, or CD8+NK neighbors per tumor cell (edge count / cells in the arm).
3. **Co-occurrence.** `gr.co_occurrence` with thresholds 50 µm and 100 µm. In Squidpy 1.8.3 those thresholds are cumulative (pairs at distance ≤ t). The score is the Squidpy co-occurrence ratio; 1 is the distance-matched base rate inside that FOV. A two-block toy (cross-type score 0 at short range, about 1 once the threshold covers both blocks) was used to check the installed function. Co-occurrence is all-pairs, so it is run per FOV, not on a whole section.

Section graphs use every cell in the section (FOV borders stay connected). That is the paired test below. FOV graphs (sample × fov; 232 FOVs) are the distribution. A FOV is graphed when it has at least 30 CLDN4-high tumor cells, 30 CLDN4-low tumor cells, and 10 CD8+NK cells: 202/232. CD8-specific FOV summaries also require ≥10 CD8 (193 FOVs). NK-specific summaries require ≥10 NK (150 FOVs).

Donor summaries give each section equal weight, then each donor one value, so Lung5’s three sections are not three donors. The one-sided sign test is a binomial test on the number of negative high−low differences. Section p-values also include an exact Wilcoxon signed-rank test. A separate 999-draw shuffle reassigns the CLDN4-high/low labels among tumor cells on the fixed Squidpy radius graph and tests the CD8+NK neighbor-count difference. That shuffle is not `nhood_enrichment`.

FOV-level tests are descriptive. FOVs are nested in sections.

## Interaction matrix (section graph)

Mean CD8+NK neighbors per tumor cell. Ratio = high / low. Permutation p is the within-tumor label shuffle, one-sided for high < low (999 draws; smallest p is 0.001).

| Section | 50 µm high | 50 µm low | 50 µm ratio | 50 µm perm p | 100 µm high | 100 µm low | 100 µm ratio | 100 µm perm p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| LUAD-5 R1 | 0.217 | 0.223 | 0.972 | 0.249 | 1.412 | 1.327 | 1.064 | 0.997 |
| LUAD-5 R2 | 0.221 | 0.257 | 0.859 | 0.001 | 1.453 | 1.423 | 1.021 | 0.873 |
| LUAD-5 R3 | 0.228 | 0.247 | 0.925 | 0.034 | 1.393 | 1.268 | 1.098 | 1.000 |
| LUSC-6 | 0.153 | 0.141 | 1.083 | 0.990 | 0.873 | 0.816 | 1.069 | 1.000 |
| LUAD-9 R1 | 0.430 | 0.664 | 0.648 | 0.001 | 2.651 | 3.416 | 0.776 | 0.001 |
| LUAD-9 R2 | 0.299 | 0.562 | 0.533 | 0.001 | 1.607 | 2.417 | 0.665 | 0.001 |
| LUAD-12 | 0.263 | 0.569 | 0.462 | 0.001 | 1.560 | 2.452 | 0.636 | 0.001 |
| LUAD-13 | 1.552 | 1.711 | 0.907 | 0.001 | 7.392 | 7.749 | 0.954 | 0.001 |

Unweighted across 8 sections:

| Radius | Pair | Mean count high | Mean count low | Median ratio | Sections high < low | Sign P | Wilcoxon P (less) | Wilcoxon P (two-sided) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 50 µm | CD8 | 0.322 | 0.408 | 0.801 | 7/8 | 0.0352 | 0.0078 | 0.0156 |
| 50 µm | NK | 0.098 | 0.139 | 0.876 | 5/8 | 0.363 | 0.0742 | 0.148 |
| 50 µm | CD8+NK | 0.420 | 0.547 | 0.883 | 7/8 | 0.0352 | 0.0117 | 0.0234 |
| 100 µm | CD8 | 1.711 | 1.942 | 0.948 | 7/8 | 0.0352 | 0.0273 | 0.0547 |
| 100 µm | NK | 0.582 | 0.666 | 0.996 | 4/8 | 0.637 | 0.273 | 0.547 |
| 100 µm | CD8+NK | 2.293 | 2.609 | 0.987 | 4/8 | 0.637 | 0.156 | 0.312 |

The section that is higher, rather than lower, around CLDN4-high is LUSC-6 at both radii for CD8, NK, and CD8+NK. LUSC-6 has 184 NK cells. LUAD-13 has 105 NK cells; its CD8 count (3,251) is the stable arm there.

At 100 µm the Lung5 serial sections have slightly more CD8+NK neighbors around CLDN4-high (ratios 1.06, 1.02, 1.10). The 100 µm CD8+NK deficit is carried by LUAD-9, LUAD-12, and LUAD-13.

Figure: `results/cosmx_squidpy_fov/figures/section_paired_cd8nk_neighbors.png`.

## Donor summary (equal weight per section)

CD8+NK neighbor-count ratio (mean of section means, high/low):

| Donor | 50 µm | 100 µm |
| --- | ---: | ---: |
| Lung5 | 0.916 | 1.060 |
| Lung6 | 1.083 | 1.069 |
| Lung9 | 0.595 | 0.730 |
| Lung12 | 0.462 | 0.636 |
| Lung13 | 0.907 | 0.954 |

CD8 alone is lower around CLDN4-high in 4/5 donors at both radii (Lung6 is higher). One-sided sign P = 0.1875. CD8+NK is lower in 4/5 donors at 50 µm (same sign P) and in 3/5 donors at 100 µm (Lung5 and Lung6 higher; sign P = 0.50).

## Neighborhood enrichment

Section-mean z versus CD8 at 50 µm is −65.5 (CLDN4-high) and −60.7 (CLDN4-low). At 100 µm the means are −76.8 and −70.0. Both arms are far below 0. That shared deficit is the all-label null on a radius graph, which mixes tumor–stroma separation into the z-score. It is not a CLDN4-specific result.

The paired CLDN4 contrast (z high minus z low) is negative in 7/8 sections for CD8 at 50 µm and at 100 µm. The positive section is LUSC-6 (z high −62.4 vs z low −87.6 at 50 µm), the same section whose neighbor counts are higher around CLDN4-high. NK z is more negative for CLDN4-high in 5/8 sections at 50 µm and 4/8 at 100 µm.

## Co-occurrence (FOV means)

Score < 1 means that effector label is less common within the radius than the FOV base rate. Both tumor arms sit below 1 for CD8 at ≤50 µm in every section (FOV-mean scores from 0.17 to 0.85). The CLDN4 contrast is the paired difference.

CD8 co-occurrence is lower around CLDN4-high than around CLDN4-low in 7/8 sections at ≤50 µm and at ≤100 µm (sign P = 0.0352). NK co-occurrence is lower in 7/8 sections at ≤50 µm and in 5/8 at ≤100 µm. LUSC-6 is again the section with a higher score around CLDN4-high. Donor sign for the CD8 difference is 4/5 (sign P = 0.1875) at both radii.

Among passing FOVs, the CD8+NK neighbor count is lower around CLDN4-high in 158/202 FOVs at 50 µm and 135/202 at 100 µm. Those fractions describe FOVs nested inside the 8 sections.

Figure: `results/cosmx_squidpy_fov/figures/section_cooccurrence_50um.png`.

## Effector transcripts in CD8/NK next to tumor

This block is the comparison to the muzzling half of the locked result. It is not a Squidpy graph score. For each author CD8 or NK cell, the nearest author tumor cell within 50 µm or 100 µm assigns the cell to CLDN4-high or CLDN4-low. The table is the ratio of mean CPM (10⁴ × count / `n_counts`).

| Gene | 50 µm sections with ratio > 1 | 50 µm median ratio | 100 µm sections with ratio > 1 | 100 µm median ratio (min–max) |
| --- | ---: | ---: | ---: | --- |
| GZMB | 7/8 | 1.189 | 8/8 | 1.184 (1.092–1.458) |
| PRF1 | 5/8 | 1.121 | 5/8 | 1.117 (0.741–1.381) |
| NKG7 | 5/8 | 1.053 | 4/8 | 1.010 (0.921–1.208) |
| IFNG | 4/8 | 1.017 | 3/8 | 0.929 (0.363–1.185) |

GZMB is higher next to CLDN4-high tumor in 8/8 sections at 100 µm. The other three genes are mixed. IFNG is detected in only a few percent of these neighbors (often under 10%), and the low ratios in LUSC-6 (0.48 at 50 µm) and LUAD-12 (0.38 at 50 µm) sit on that sparse count. Sample sizes and detection fractions are in `results/cosmx_squidpy_fov/tables/muzzling_ratios.csv`.

Figure: `results/cosmx_squidpy_fov/figures/muzzling_cpm_ratio_50um.png`.

## Comparison with the locked result

Same object (765,771 cells, 8 sections, 5 donors) and the same author tumor and CD8/NK labels. Different statistic.

**Exclusion.** The Squidpy radius-graph counts show fewer CD8 neighbors around CLDN4-high tumor in 7/8 sections at 50 µm and at 100 µm (median high/low 0.801 and 0.948). Co-occurrence shows the same 7/8. Combined CD8+NK counts are lower in 7/8 sections at 50 µm (median ratio 0.883) and in 4/8 sections at 100 µm (median ratio 0.987). Donors are 4/5 for CD8 at both radii (Lung6 is higher around CLDN4-high; one-sided sign P = 0.1875). The locked donor result remains 5/5 with sign P = 0.031, from its own neighbor-ratio definition. Median ratios in this file are milder than the locked 0.36 and 0.52. Those locked ratios stay as locked.

**Muzzling.** GZMB CPM next to CLDN4-high tumor is higher in 8/8 sections at 100 µm (median ratio 1.18), in the same direction as the locked high/low range 1.11–1.22. PRF1, NKG7, and IFNG are not higher in every section under this nearest-tumor CPM definition, so this file does not restate “0/8 decreased” for all four genes. No gene is lower in 8/8 sections. The graph scores themselves only measure proximity.

## What this file does not say

No private 8-KL matrix, no ICI label, no TACSTD2 gate, no CellChat probability, no claim that a Visium spot correlation is spatial exclusion. Absolute neighborhood z-scores and absolute co-occurrence scores below 1 are reported as tissue architecture shared by CLDN4-low tumor.

## Rerun

The h5ad is not in git (about 2.6 GB). Download figshare file 46841842 to `data/cosmx_human_nsclc_clustered.h5ad`, then:

```bash
python3 scripts/cosmx_squidpy_cldn4_cd8nk.py
```

Tables: `results/cosmx_squidpy_fov/tables/`. Full test dump: `summary.json`.
