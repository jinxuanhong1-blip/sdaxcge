# A11 — galectin genes vs TACSTD2-high public lung

## Honest result (≤200 words)
TCGA primary tumors, analyzed separately: LUAD n=510, LUSC n=484. TACSTD2-high is the upper half, not a validated subtype.

- **LGALS3:** TCGA_LUAD ρ=0.38, high/low log2 ratio=+0.41, q=2.8e-17; TCGA_LUSC ρ=0.30, high/low log2 ratio=+0.39, q=1.1e-10
- **LGALS9C:** TCGA_LUAD ρ=0.22, high/low log2 ratio=+0.60, q=2.4e-06; TCGA_LUSC ρ=0.13, high/low log2 ratio=+0.30, q=0.01
- **LGALS9B:** TCGA_LUAD ρ=0.19, high/low log2 ratio=+0.90, q=4.8e-05; TCGA_LUSC ρ=0.14, high/low log2 ratio=+0.62, q=0.006

**LGALS3** is the only gene that is TCGA-concordant, replicated in independent LUAD bulk (OncoSG ρ=0.41; CPTAC LUAD ρ=0.39), retained after CPTAC LUAD ESTIMATE-purity partial Spearman (ρ=0.38), and present in tumor-cell-only CCLE 2025 (lung n=206 ρ=0.68; NSCLC n=142 ρ=0.58). It did **not** replicate in CPTAC LUSC bulk (q>0.05).

**LGALS9B/9C** stay positive in TCGA and CCLE, but independent bulk is weak: 9C none; 9B only CPTAC LUSC. Both are low-abundance in several matrices.

Discordant TCGA genes: LGALS1, LGALS2, LGALS9, LGALS12. LUSC-only: LGALS7 (TCGA_LUSC), LGALS7B (TCGA_LUSC). Unavailable: CAS_LUAD (no TACSTD2 in the public matrix). Associations only; not proof that TROP2 regulates galectins. Protein, spatial, ICI, and single-cell malignant-cell effects were not tested.

## Methods
Spearman is primary. High/low ratios and Mann–Whitney tests are descriptive. BH correction is within family (TCGA; other bulk; cell lines). Partial Spearman residualizes ranks on ESTIMATE purity. Scales are not pooled.

## Sources
- [TCGA LUAD](https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018)
- [TCGA LUSC](https://www.cbioportal.org/study/summary?id=lusc_tcga_pan_can_atlas_2018)
- [OncoSG LUAD](https://www.cbioportal.org/study/summary?id=luad_oncosg_2020)
- [CPTAC LUAD](https://www.cbioportal.org/study/summary?id=luad_cptac_2020)
- [CPTAC LUSC](https://www.cbioportal.org/study/summary?id=lusc_cptac_2021)
- [CAS LUAD](https://www.cbioportal.org/study/summary?id=luad_cas_2020)
- [CCLE Broad 2025](https://www.cbioportal.org/study/summary?id=ccle_broad_2025)
