# LUAD-only recut: CLDN4 vs CD8 / TNK / Immune

**Additive histology split of existing public sets.** Mixed LUAD+LUSC often dilutes or blends two diseases. This file re-scores **CLDN4** against CD8 / TNK / Immune on **LUAD-only** subsets where histology is a deposited label. **LUSC is kept as a separate extra row**, not dropped and not averaged into the LUAD claim.

No new private data. Patient (or primary-tumor sample) is the unit. Honest n in every cell.

## Verdict

On **TCGA-LUAD** (n=515 primaries) CLDN4 is only weakly CD8-low / CYT-low. The mixed TCGA NSCLC CD8 number is a blend: LUAD -0.118 (p=0.00721, n=515) vs LUSC extra -0.017 (p=0.699, n=502) vs mixed -0.022 (p=0.489, n=1017). After ABSOLUTE, LUAD CD8 is -0.091 (p=0.0416, n=502) and LUSC is -0.031 (p=0.488, n=493). **Do not quote the mixed TCGA rho as a LUAD finding.**

**CPTAC-LUAD** is already LUAD-only: CLDN4 RNA vs CIBERSORT CD8 is weakly negative and not a strong ImmuneScore hit. **CPTAC-LSCC protein** (LUSC extra) is the stronger anti-CD8 / anti-ImmuneScore protein row.

**GSE218989** and **GSE285029** stay mixed. Histology is not on GEO and not in the public supplements checked here. **GSE148071** official Supplementary Table 1 is aggregate only (18 ADC / 18 SQ / 6 NSCLC); GEO has age/sex only. Fig S2 has per-patient subtypes but no machine-readable table was deposited — **not split**.

scRNA that is already LUAD: **GSE131907 tLung** and **GSE253013**. Both are small-n and not a ρ ≈ −0.4 to −0.5 T/NK package.

## Split table (primary)

| cohort | histology | layer | CLDN4 vs | n | ρ | p | adjust | why this row |
|---|---|---|---|---:|---:|---:|---|---|
| TCGA | **LUAD** | bulk_RNA | CD8A | 515 | -0.118 | 0.00721 | none | Xena HiSeqV2 primary -01; ESTIMATE ImmuneScore; ABSOLUTE purity |
| TCGA | **LUAD** | bulk_RNA | CYT | 515 | -0.167 | 0.000146 | none | Xena HiSeqV2 primary -01; ESTIMATE ImmuneScore; ABSOLUTE purity |
| TCGA | **LUAD** | bulk_RNA | ImmuneScore | 515 | -0.065 | 0.143 | none | Xena HiSeqV2 primary -01; ESTIMATE ImmuneScore; ABSOLUTE purity |
| TCGA | **LUAD** | bulk_RNA | CD8A | 502 | -0.091 | 0.0416 | ABSOLUTE | Xena HiSeqV2 primary -01; ESTIMATE ImmuneScore; ABSOLUTE purity |
| TCGA | **LUSC** | bulk_RNA | CD8A | 502 | -0.017 | 0.699 | none | Xena HiSeqV2 primary -01; ESTIMATE ImmuneScore; ABSOLUTE purity |
| TCGA | **LUSC** | bulk_RNA | CYT | 502 | -0.055 | 0.22 | none | Xena HiSeqV2 primary -01; ESTIMATE ImmuneScore; ABSOLUTE purity |
| TCGA | **LUSC** | bulk_RNA | ImmuneScore | 501 | -0.004 | 0.925 | none | Xena HiSeqV2 primary -01; ESTIMATE ImmuneScore; ABSOLUTE purity |
| TCGA | **LUSC** | bulk_RNA | CD8A | 493 | -0.031 | 0.488 | ABSOLUTE | Xena HiSeqV2 primary -01; ESTIMATE ImmuneScore; ABSOLUTE purity |
| TCGA | **mixed** | bulk_RNA | CD8A | 1017 | -0.022 | 0.489 | none | Xena HiSeqV2 primary -01; ESTIMATE ImmuneScore; ABSOLUTE purity; pooled LUAD+LUSC (dilution / Simpson risk) |
| CPTAC-LUAD | **LUAD** | bulk_RNA | CIBERSORT_CD8 | 110 | -0.225 | 0.018 | none | CPTAC LUAD Gillette 2020; already LUAD-only |
| CPTAC-LUAD | **LUAD** | bulk_RNA | ImmuneScore | 110 | -0.079 | 0.413 | none | CPTAC LUAD Gillette 2020; already LUAD-only |
| CPTAC-LSCC | **LUSC** | protein | CD8A_RNA | 78 | -0.437 | 6.4e-05 | none | CPTAC LSCC Satpathy 2021; LUSC extra row |
| CPTAC-LSCC | **LUSC** | protein | ImmuneScore | 78 | -0.432 | 7.9e-05 | none | CPTAC LSCC Satpathy 2021; LUSC extra row |
| GSE218989 | **mixed_unsplit** | bulk_RNA | CD8A | 355 | -0.172 | 0.00113 | none | Kang 2024 SMC-KAIST ICI TPM; GEO+Supp8 have no histology column; cannot LUAD-split |
| GSE285029 | **mixed_unsplit** | bulk_RNA | CD8A | 234 | +0.067 | 0.305 | none | Koh 2025 JITC pre-ICI WTS n=234; GEO characteristics = tissue/cell/genotype only; no histology |
| GSE131907_tLung | **LUAD** | scRNA | TNK_frac | 11 | -0.045 | 0.894 | none | Kim 2020 LUAD atlas; primary tLung; author epithelial mean vs T/NK fraction |
| GSE131907_tLung | **LUAD** | scRNA | CD8_frac | 11 | +0.336 | 0.312 | none | Kim 2020 LUAD atlas; primary tLung; author epithelial mean vs CD8 fraction |
| GSE253013 | **LUAD** | scRNA | TNK_frac | 9 | -0.333 | 0.381 | none | Sze/Xiang 2024 treatment-naive LUAD tumors only (ANT held out); malignant-like mean vs T/NK |
| GSE148071 | **mixed_unsplit** | scRNA | T_frac_author | 42 | +0.081 | 0.611 | none | Wu 2021 advanced NSCLC; GEO age/sex only; Supp Table 1 is aggregate 18 ADC/18 SQ/6 NSCLC — no per-patient histology table |

Full numeric table: `tables/split_table.tsv`. Every test: `tables/all_tests.tsv`.

## Histology inventory (nothing invented)

| cohort | histology source | n LUAD | n LUSC | n mixed / unknown | split? |
|---|---|---:|---:|---:|---|
| TCGA-LUAD | project is LUAD | 515 | 0 | 0 | yes |
| TCGA-LUSC | project is LUSC | 0 | 502 | 0 | yes |
| CPTAC-LUAD | CPTAC LUAD proteogenomic cohort | 110 | 0 | 0 | yes |
| CPTAC-LSCC | CPTAC LSCC / LUSC proteogenomic cohort | 0 | 108 | 0 | yes |
| GSE218989 | none on GEO or Supp Data 8 | 0 | 0 | 355 | **no** |
| GSE285029 | none on GEO series matrix | 0 | 0 | 234 | **no** |
| GSE131907 | Kim 2020 LUAD atlas (all LUAD) | 11 | 0 | 0 | yes |
| GSE253013 | Sze/Xiang 2024 LUAD | 9 | 0 | 0 | yes |
| GSE148071 | Supp Table 1 aggregate only; GEO age/sex only | 0 | 0 | 42 | **no** |

## TCGA (recomputed here)

Xena `TCGA.LUAD.sampleMap/HiSeqV2` and `TCGA.LUSC.sampleMap/HiSeqV2`, log2(RSEM+1). Primary tumors only (`-01`), one sample per patient. ImmuneScore = MD Anderson ESTIMATE RNAseqV2. Purity = PanCanAtlas ABSOLUTE. CYT = mean(GZMA, PRF1). GEP18 = unweighted mean of within-cohort z-scores of the 18 Ayers genes (not NanoString TIS weights).

| endpoint | LUAD unadj | LUSC extra unadj | mixed unadj | LUAD \| ABSOLUTE | LUSC \| ABSOLUTE |
|---|---|---|---|---|---|
| CD8A | -0.118 (p=0.00721, n=515) | -0.017 (p=0.699, n=502) | -0.022 (p=0.489, n=1017) | -0.091 (p=0.0416, n=502) | -0.031 (p=0.488, n=493) |
| CYT | -0.167 (p=0.000146, n=515) | -0.055 (p=0.22, n=502) | -0.071 (p=0.0236, n=1017) | -0.152 (p=0.000634, n=502) | -0.074 (p=0.0989, n=493) |
| ImmuneScore | -0.065 (p=0.143, n=515) | -0.004 (p=0.925, n=501) | +0.109 (p=0.000486, n=1016) | -0.018 (p=0.681, n=502) | -0.019 (p=0.675, n=493) |
| GEP18 | -0.122 (p=0.00567, n=515) | -0.048 (p=0.283, n=502) | -0.082 (p=0.00856, n=1017) | -0.096 (p=0.031, n=502) | -0.070 (p=0.12, n=493) |

LUAD n=515 primaries (502 with ABSOLUTE). LUSC n=502 primaries (493 with ABSOLUTE). Mixed n=1017.

The mixed CD8 rho sits between the two histology clouds and is **closer to the LUSC null** than to the LUAD weak-negative. That is dilution. Mixed ImmuneScore even **flips sign** (+0.11) while both histologies are near zero / weakly negative — a between-histology Simpson blend. LUSC is the extra row; it is not the LUAD claim.

## CPTAC (already split by design)

CPTAC-LUAD (Gillette 2020) and CPTAC-LSCC (Satpathy 2021) are separate public proteogenomic sets. No pooling.

| cohort | histology | CLDN4 vs | n | ρ (unadj) |
|---|---|---|---:|---|
| CPTAC-LUAD | LUAD | RNA vs CIBERSORT CD8 | see table | -0.225 (p=0.018, n=110) |
| CPTAC-LUAD | LUAD | RNA vs ImmuneScore | see table | -0.079 (p=0.413, n=110) |
| CPTAC-LSCC | **LUSC extra** | protein vs CD8A RNA | see table | -0.437 (p=6.4e-05, n=78) |
| CPTAC-LSCC | **LUSC extra** | protein vs ImmuneScore | see table | -0.432 (p=7.9e-05, n=78) |

LUSC protein vs CD8/Immune is the stronger negative. LUAD RNA is weaker. Do not average them.

## GSE218989 and GSE285029 — cannot split

**GSE218989** (Kang et al., *Nat Commun* 2024; 355 ICI TPM patients): GEO series matrix has treatment / outcome / ethnicity. Published Supplementary Data 8 (`Clinical_table_v230613`) joins OS/PFS time and has **no histology column**. Mixed CLDN4–CD8A ρ = −0.172 (p=0.0011, n=355) is a mixed-NSCLC number. Not a LUAD-only recut.

**GSE285029** (Koh et al., *J Immunother Cancer* 2025; 234 pre-ICI WTS): GEO characteristics are tissue=Lung, cell type=cancer, genotype=wt. No histology, response, or purity. Mixed CLDN4–CD8A is near zero; CLDN4–GEP18 is weakly positive. Not split.

## scRNA

| cohort | histology | n | CLDN4 vs | ρ (p) |
|---|---|---:|---|---|
| GSE131907 tLung | LUAD (Kim 2020 atlas) | see table | epi mean vs T/NK frac | -0.045 (p=0.894, n=11) |
| GSE131907 tLung | LUAD | see table | epi mean vs CD8 frac | +0.336 (p=0.312, n=11) |
| GSE253013 | LUAD (Sze/Xiang 2024) | see table | malig-like mean vs T/NK | -0.333 (p=0.381, n=9) |
| GSE148071 | **mixed, unsplit** | 42 | malig mean vs author T frac | +0.081 (p=0.611, n=42) |

GSE131907 and GSE253013 do not need a histology recut — they are LUAD cohorts. GSE148071 is the mixed advanced-NSCLC atlas (Wu 2021). Official Supplementary Table 1 reports 18 adenocarcinoma / 18 squamous / 6 NSCLC **in aggregate**. GEO SOFT has age and sex only. Fig S2 colors patients by pathology vs scRNA-combined subtype, but that is a figure, not a deposited per-patient table. This recut does **not** invent P1–P42 labels.

## How to read this against mixed-NSCLC CLDN4 claims

1. **LUAD-only TCGA CLDN4–CD8 is small and real.** Mixed TCGA CD8 (ρ = −0.022) is not “the LUAD association”; it is LUAD diluted by a LUSC-null cloud.
2. **LUSC extra is not one story.** TCGA-LUSC CLDN4–CD8/CYT/Immune is near zero. CPTAC-LSCC **protein** vs CD8/Immune is strongly negative (n=78). Keep the protein row extra; do not average RNA and protein.
3. **ICI bulks that lack histology stay mixed.** GSE218989 CLDN4–CD8A ρ = −0.17 is real on 355 patients and still not a LUAD-only result.
4. **scRNA LUAD n is honest and small.** tLung n=11 and GSE253013 n=9 are not a meta-analysis substitute.

## Methods

- TCGA: stream Xena HiSeqV2; keep `-01` primaries; one row per patient. Spearman; ABSOLUTE partial = Pearson of rank residuals. Fisher-z 95% CI (`1/√(n−3)` unadj; `1/√(n−4)` partial).
- CPTAC / scRNA / GSE218989 / GSE285029: harvested public patient tables from prior additive PRs; statistics recomputed here on those tables except the two GEO mixed rows taken from the already-computed public matrices (histology absent, so no recut).
- GSE148071 T fraction: official Supplementary Data 1 T_cell / total cells (Wu 2021).
- No inferred histology. No ICI-response claim from TCGA/CPTAC.

Reproduce:

```
pip install -r methods/luad_only_cldn4/requirements.txt
python3 methods/luad_only_cldn4/analyze.py
```

TCGA raw matrices are downloaded to `/tmp/luad_only_cldn4` (gitignored). Harvested tables live in `data/`.

## Files

- `tables/split_table.tsv` — primary LUAD / LUSC-extra / mixed-unsplit rows
- `tables/all_tests.tsv` — every Spearman computed here
- `tables/histology_inventory.tsv`
- `tables/tcga_luad_samples.tsv`, `tables/tcga_lusc_samples.tsv`
- `figures/fig1_forest_split.png`
- `figures/fig2_tcga_luad_vs_lusc.png`
- `figures/fig3_cptac_luad_vs_lusc.png`
- `figures/fig4_scrna_luad_and_mixed.png`
- `figures/fig5_tcga_dilution_cd8a.png`
