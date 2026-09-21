# HTAN open lung scRNA: CLDN4 in malignant cells vs DNA-damage and STING scores

**Slice:** HTAN portal lung scRNA-seq inventory, then scores on the two combined Level-4 objects that are anonymously downloadable.  
**Skipped:** dbGaP `phs002371` raw FASTQ/BAM. Synapse per-sample matrices that return HTTP 403 without a session.  
**No fabricated statistics.** ρ / p / q below were computed by `scripts/htan_cldn4_sting/02_analyze.py`.

This does not revise the locked CosMx exclusion result or the concordant-4 patient-level immune correlation. These objects have no ICI response labels.

---

## What was tested

Within open HTAN lung single-cell objects, do **CLDN4-higher malignant cells** carry higher **DNA-repair / DSB** scores and higher **cGAS–STING** scores than CLDN4-lower cells from the same donor?

The inferential unit is the **donor**. For each donor with at least 40 cells in the stratum, a within-donor Spearman ρ was computed. The reported test is a two-sided Wilcoxon signed-rank test on those donor ρ values. BH q is computed inside the 18-test family (3 strata × 6 modules). Cell-pooled Spearman is saved only as a descriptive table: with tens of thousands of cells the p-values are not an independent-sample test.

A second, separate question is the **between-donor** Spearman of donor-mean CLDN4 vs donor-mean score. Those two questions can point different ways.

## Data used

| Object | Portal file | CELLxGENE h5ad | What was scored |
|---|---|---|---|
| Chan et al. *Cancer Cell* 2021 combined | `syn23626795` `adata.combined.mnnc.010920.h5ad` | `a9d92e38-9a6e-401b-8484-74bb15122341.h5ad` (1.40 GB, 147,137 cells) | Author labels `SCLC-A/N/P` (malignant SCLC) and `NSCLC` in LUAD histology |
| Glasner et al. *Nat Immunol* 2023 human LUAD | `syn51033592` `glasner_etal_globalAnndata_20230112.vHTA.h5ad` | `3e1c5296-5c2d-49df-8423-268bcb01f175.h5ad` (598 MB, 82,991 cells) | `cell_lineage == Epithelial` in LUAD specimens |

Expression in both objects is already on a compressed scale (99th percentile ≈ 2.3). No extra log1p was applied. Scores are the mean of per-gene z-scores **inside the stratum**, so immune cells do not set the scale.

Glasner epithelium is **not** an author malignant call and **not** a CNV call. Every specimen in that object is LUAD (primary, metastasis, or recurrence), and the epithelial compartment still mixes tumor cells with residual lung epithelium. CLDN4 median in that compartment is 0 (30.3% of cells > 0), so the high/low split is `CLDN4 > 0` vs `CLDN4 = 0` when the donor median is 0.

## Inventory (portal, 2026-09-21)

Lung `scRNA-seq` files in `htan_2026_932`: **2,390**. Phase-2 database `htan2_2026_926` returned **0** lung RNA rows.

| Atlas | Open on Synapse (Level 3–4) | Controlled (Level 1–2, dbGaP) | Scored here |
|---|---|---|---|
| HTAN MSK | L3 csv 338, mtx 110; L4 csv 191, txt 76, hdf5 13 | fastq 680, bam 118 | Chan combined + Glasner human LUAD h5ads (CELLxGENE copies of two L4 hdf5 files) |
| HTAN BU (Lung PCA, Linh/UCLA GGO–tumor–normal) | L3 csv 45 (logNorm / raw / batchAdj); L4 csv 32 (**cluster IDs only**) | fastq 504, bam 16 | No. Anonymous file API returned HTTP 403. Raw is `phs002371` |
| HTAN HTAPP lung | L3 mtx 34 + tsv 68; L4 tsv 35 (LUAD / mixed / one squamous) | fastq 78, bam 34 | No. Same HTTP 403 |
| HTAN DFCI lung | none | fastq 18 | No processed matrix |

Full file list: `notes/htan_cldn4_sting/lung_scrna_files.tsv`. Later Chan `adata.scvi.h5ad` (`syn68258699`) is on Synapse and was not on the anonymous CELLxGENE collection used here.

Mouse Glasner Visium slides are a different assay and were not scored.

## Modules

| Score | Definition | Coverage |
|---|---|---|
| `HALLMARK_DNA_REPAIR` | MSigDB Hallmark 2020, 150 genes | 148/150 Chan (missing `ZNRD1`, `MRPL40`); 149/150 Glasner (missing `ZNRD1`) |
| `DDR_DSB` | 24 HGNC DSB/checkpoint genes (`ATM`, `ATR`, `CHEK1/2`, `H2AX`, `TP53BP1`, `BRCA1/2`, `RAD51`, …) | 24/24. `H2AFX` accepted as an alias of `H2AX` |
| `PROLIF` | 10-gene cycle control (`MKI67`, `TOP2A`, `PCNA`, `MCM2/4/6`, `CDK1`, `CCNB1`, `BIRC5`, `TYMS`) | 10/10 |
| `STING_SENSOR` | `CGAS`, `STING1`, `TBK1`, `IKBKE`, `IRF3`, `IRF7` | 6/6. `TMEM173` accepted as an alias of `STING1` |
| `STING_ISG` | 14 type-I ISGs (`ISG15`, `IFIT1/2/3`, `MX1/2`, `OAS1/2/3`, `CXCL10`, `CCL5`, `IFI44`, `IFI44L`, `RSAD2`) | 14/14 |
| `REACTOME_STING` | R-HSA-1834941 (16 genes) | 15/16 (missing `TREX1`). Overlaps DSB genes (`MRE11`, `PRKDC`, `XRCC5/6`), so it is not an independent STING test |

Hallmark DNA repair is heavy on `POLR2` subunits. `DDR_DSB` is the tighter damage set. Both are reported.

## Strata

| Stratum | Cells | Donors with any cells | Donors in the within-donor test | CLDN4 mean | CLDN4 > 0 | CLDN4 median |
|---|---:|---:|---:|---:|---:|---:|
| Chan `SCLC-A/N/P` | 54,313 | 19 | 19 | 0.799 | 62.1% | 0.639 |
| Chan `NSCLC` and histo LUAD | 3,121 | 20 | 8 (≥40 cells and non-constant CLDN4) | 1.070 | 73.2% | 0.987 |
| Glasner LUAD epithelium | 21,855 | 23 | 21 | 0.571 | 30.3% | 0 |
| Chan normal-histology epithelium | 754 | 4 | not tested (n donors < 6) | 1.326 | 65.0% | 1.265 |
| Chan `NSCLC` label inside SCLC histology | 220 | 18 | not tested (almost all donors < 40 cells) | 1.660 | 50.5% | 0.592 |

## Results

### 1. Within a donor, CLDN4-higher cells score higher for DNA repair and STING

Median within-donor Spearman ρ. Fraction is donors with ρ > 0. q is BH across the 18 primary tests.

**Chan malignant SCLC (19 donors)**

| Score | Median ρ | Donors ρ>0 | Wilcoxon p | q |
|---|---:|---:|---:|---:|
| Hallmark DNA repair | 0.125 | 19/19 | 3.8×10⁻⁶ | 1.0×10⁻⁵ |
| DDR / DSB | 0.042 | 15/19 | 1.2×10⁻³ | 0.0021 |
| Proliferation | 0.002 | 10/19 | 0.57 | 0.57 |
| STING sensor | 0.138 | 19/19 | 3.8×10⁻⁶ | 1.0×10⁻⁵ |
| STING ISG | 0.122 | 18/19 | 7.6×10⁻⁶ | 1.7×10⁻⁵ |
| Reactome STING | 0.107 | 17/19 | 5.3×10⁻⁵ | 1.1×10⁻⁴ |

SCLC proliferation is flat, so the DNA-repair and STING ρ values are not a cell-cycle echo. After residualizing each module on the proliferation score inside the donor, median partial ρ stays positive: DNA repair 0.156 (19/19), DDR/DSB 0.058 (16/19), STING sensor 0.139 (19/19), STING ISG 0.126 (18/19); all Wilcoxon p ≤ 1.6×10⁻⁴.

The SCLC effect is small. Pooled across cells, DDR/DSB vs CLDN4 is ρ = −0.037. That pooled sign is a mixture of donors and is not the donor test.

**Chan malignant NSCLC in LUAD histology (8 donors with ≥40 cells)**

All six modules are positive in 8/8 donors. Wilcoxon p = 0.0078, q = 0.0083 for each (the smallest two-sided Wilcoxon p at n = 8 when every ρ is positive).

| Score | Median ρ |
|---|---:|
| Hallmark DNA repair | 0.352 |
| DDR / DSB | 0.265 |
| Proliferation | 0.204 |
| STING sensor | 0.274 |
| STING ISG | 0.324 |
| Reactome STING | 0.239 |

Partial ρ given proliferation remains positive in 8/8 donors (DNA repair 0.323, DDR/DSB 0.220, STING sensor 0.237, STING ISG 0.291). Twelve of the 20 LUAD donors with an `NSCLC` label have fewer than 40 such cells and are not in this test.

**Glasner LUAD epithelium (21 donors; not a malignant call)**

| Score | Median ρ | Donors ρ>0 | Wilcoxon p | q |
|---|---:|---:|---:|---:|
| Hallmark DNA repair | 0.277 | 20/21 | 1.9×10⁻⁶ | 9×10⁻⁶ |
| DDR / DSB | 0.102 | 20/21 | 1.9×10⁻⁶ | 9×10⁻⁶ |
| Proliferation | 0.084 | 16/21 | 0.0055 | 0.0083 |
| STING sensor | 0.194 | 21/21 | 9.5×10⁻⁷ | 9×10⁻⁶ |
| STING ISG | 0.204 | 20/21 | 2.9×10⁻⁶ | 1.0×10⁻⁵ |
| Reactome STING | 0.219 | 20/21 | 1.9×10⁻⁶ | 9×10⁻⁶ |

Partials given proliferation stay in the same direction (DNA repair 0.264, 20/21; STING sensor 0.193, 21/21). The same direction is present in the primary-specimen epithelium subset (17 donors; Hallmark DNA repair median ρ 0.280, STING sensor 0.194). Metastasis-labeled epithelium has only 4 donors and was not tested.

Single genes, median within-donor ρ: **STING1** is positive (SCLC 0.083, 17/19; NSCLC 0.202, 7/8; Glasner 0.223, 21/21). **CGAS** is weak (SCLC 0.020, 13/18; NSCLC 0.159, 8/8; Glasner 0.065, 21/21). **H2AX** is positive (SCLC 0.090, 18/19; NSCLC 0.238, 7/7; Glasner 0.205, 21/21). The sensor module is not a CGAS-only effect.

### 2. Between donors, the same means mostly do not reproduce that within-donor pattern

Donor-mean Spearman, donors with ≥40 cells. q is BH inside the same 18 between-donor tests.

| Stratum | Score | n | ρ | p | q |
|---|---|---:|---:|---:|---:|
| Chan SCLC | Hallmark DNA repair | 19 | −0.344 | 0.15 | 0.21 |
| Chan SCLC | DDR / DSB | 19 | −0.282 | 0.24 | 0.29 |
| Chan SCLC | STING sensor | 19 | +0.384 | 0.10 | 0.21 |
| Chan SCLC | STING ISG | 19 | +0.070 | 0.78 | 0.78 |
| Chan NSCLC | STING sensor | 8 | +0.762 | 0.028 | 0.17 |
| Chan NSCLC | Hallmark DNA repair | 8 | +0.571 | 0.14 | 0.21 |
| Glasner epithelium | STING sensor | 21 | **+0.675** | **7.8×10⁻⁴** | **0.014** |
| Glasner epithelium | STING ISG | 21 | +0.562 | 0.0080 | 0.072 |
| Glasner epithelium | Hallmark DNA repair | 21 | +0.414 | 0.062 | 0.19 |
| Glasner epithelium | DDR / DSB | 21 | −0.149 | 0.52 | 0.55 |

The one between-donor result that survives this FDR is Glasner epithelial **STING sensor** (and it remains after residualizing donor-mean proliferation: partial ρ = 0.729, p = 1.8×10⁻⁴, n = 21). Glasner is the unlabeled epithelial compartment, so that row is not a malignant-cell claim.

Chan SCLC between-donor DNA repair points the other way from the within-donor test and is not significant. Patients with higher average malignant CLDN4 are not a DNA-repair-high group in this atlas.

## Reading

Inside a tumor, CLDN4-higher malignant SCLC cells and CLDN4-higher NSCLC cells have higher DNA-repair and STING-module scores than CLDN4-lower cells from the same donor. Proliferation does not account for that within-donor pattern. The SCLC DNA-repair correlation is small, and the tighter DSB gene set is weaker than the Hallmark set.

Across donors, average CLDN4 does not mark DNA-repair-high SCLC. These tables do not support a claim that CLDN4-high malignant cells are STING-low.

## Reproduce

```bash
pip install -r scripts/htan_cldn4_sting/requirements.txt
mkdir -p /tmp/htan_sting/data
curl -fL -o /tmp/htan_sting/data/chan_combined.h5ad \
  https://datasets.cellxgene.cziscience.com/a9d92e38-9a6e-401b-8484-74bb15122341.h5ad
curl -fL -o /tmp/htan_sting/data/glasner_luad.h5ad \
  https://datasets.cellxgene.cziscience.com/3e1c5296-5c2d-49df-8423-268bcb01f175.h5ad
python3 scripts/htan_cldn4_sting/00_summarize_inventory.py
python3 scripts/htan_cldn4_sting/02_analyze.py
```

Tables: `results/htan_cldn4_sting/tables/`. Figures: `donor_mean_cldn4_vs_scores.png`, `median_within_donor_rho.png`.
