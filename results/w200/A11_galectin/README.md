# A11 — galectin genes vs TACSTD2-high public lung

## Honest result (≤200 words)
Public TCGA primary tumors were analyzed separately: LUAD n=510 and LUSC n=484. TACSTD2-high means the upper within-cohort half; it is not a validated biological subtype.

- **LGALS3:** LUAD ρ=0.38, high/low log2 ratio=+0.41, q=2.8e-17; LUSC ρ=0.30, high/low log2 ratio=+0.39, q=1.1e-10
- **LGALS9C:** LUAD ρ=0.22, high/low log2 ratio=+0.60, q=2.4e-06; LUSC ρ=0.13, high/low log2 ratio=+0.30, q=0.01
- **LGALS9B:** LUAD ρ=0.19, high/low log2 ratio=+0.90, q=4.8e-05; LUSC ρ=0.14, high/low log2 ratio=+0.62, q=0.006

Discordant significant genes: LGALS1, LGALS2, LGALS9, LGALS12. These are associations, not evidence that TROP2 regulates galectins. Bulk RNA mixes malignant, stromal, and immune cells; TACSTD2 is epithelial, whereas several galectins are also abundant in non-malignant compartments. Histology-stratified agreement is therefore the strongest claim supported here. Protein abundance, spatial co-expression, treatment response, and single-cell malignant-cell effects were not tested.
One-cohort signals with matching direction: LGALS7 (LUSC), LGALS7B (LUSC).

## Methods
Values are cBioPortal PanCancer Atlas batch-normalized RNA-seq RSEM. Spearman correlation is primary; high/low median ratios and Mann–Whitney tests are descriptive. BH correction covers 26 testable gene-by-cohort comparisons; LGALS16 was unavailable.

## Sources
- [LUAD study](https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018)
- [LUSC study](https://www.cbioportal.org/study/summary?id=lusc_tcga_pan_can_atlas_2018)
- [cBioPortal API](https://www.cbioportal.org/api/swagger-ui/index.html)
