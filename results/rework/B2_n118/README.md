# B2 rework · find the n≈118 CCLE NSCLC **protein** table

**Honest verdict: n=118 cannot be found on the public CCLE proteomics table.**

The exact public table is Nusinow et al. 2020 / Gygi lab `protein_quant_current_normalized.csv.gz` (375 unique CCLE lines). No natural lung or NSCLC filter of that table has n≈118.

| Filter on that table | n lines | n pairwise TACSTD2+CLDN4 | Spearman ρ | vs user 0.69 |
|---|---:|---:|---:|---|
| User claim | **118** | **118** | **0.69** | — |
| Table S1 Tissue = Lung | 77 | 45 | **0.693 → 0.69** | ρ matches; **n does not** |
| Oncotree NSCLC | 63 | 35 | 0.727 → 0.73 | nearby ρ; **n does not** |
| Oncotree Lung | 76 | 44 | 0.704 → 0.70 | nearby ρ; **n does not** |

n=118 is a published CCLE **mRNA** NSCLC count (Augustyn et al. *PNAS* 2014), not a protein-table n. We did not mash tissues together to force n=118.

---

## 中文

### 一句话结论

**公开 CCLE 蛋白组表上找不到 n=118 的 NSCLC 切片。** 能对上的精确公开表是 Nusinow 2020 / Gygi 的 `protein_quant_current_normalized.csv.gz`（375 个独立细胞系）。该表的肺系是 **77** 株，Oncotree NSCLC 是 **63** 株；TACSTD2+CLDN4 成对完整病例分别是 **45** 和 **35**。用户的 **ρ=0.69** 对得上 S1 肺系成对完整病例（n=45，ρ=0.693），但对不上 n=118。

n=118 来自 CCLE **转录组**文献（Augustyn *PNAS* 2014：NSCLC mRNA n=118），不是蛋白表。

---

## Exact public table

| Item | Value |
|---|---|
| File | `protein_quant_current_normalized.csv.gz` |
| URL | https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz |
| Excel twin (manuscript Table S2) | https://gygi.hms.harvard.edu/data/ccle/Table_S2_Protein_Quant_Normalized.xlsx |
| Sample info (Table S1) | https://gygi.hms.harvard.edu/data/ccle/Table_S1_Sample_Information.xlsx |
| Gygi page | https://gygi.hms.harvard.edu/publications/ccle.html |
| Paper | Nusinow DP et al. *Cell* 2020;180:387-402.e16 |
| DOI | 10.1016/j.cell.2019.12.023 |
| SHA256 (this run) | `b72a9ff39b123c801249ea212844ac1f500c3f38fa0f399ec6ad08ee1be2a7ed` |
| Protein rows | 12,755 |
| Sample columns | 378 (3 lines run in two TMT plexes) |
| Unique CCLE codes | **375** |
| TACSTD2 | UniProt P09758, quantified in all 375 lines |
| CLDN4 | UniProt O14493, quantified in 197 / 375 lines (many NA; low peptide coverage) |

DepMap’s “Proteomics” download of `protein_quant_current_normalized.csv` is this same Gygi file (the portal link redirects to the Gygi page). cBioPortal `ccle_broad_2019` `data_protein_quantification.txt` is a gene-collapsed copy of the same matrix, not a different n=118 table.

Lineage labels used only to *name* NSCLC vs SCLC, not to invent a new matrix:

| File | URL | Role |
|---|---|---|
| `Table_S1_Sample_Information.xlsx` | Gygi, above | CCLE 2012 `Tissue of Origin` (Lung = 77 unique lines) |
| `Model.csv` | DepMap Public 24Q4 Figshare+ file 51065297 | `OncotreePrimaryDisease` / `OncotreeSubtype` |

---

## Why n=118 is not in this table

Natural filters, unique cell lines:

| Cohort | n lines in table | n with TACSTD2 | n with CLDN4 | n both |
|---|---:|---:|---:|---:|
| All unique lines | 375 | 375 | 197 | 197 |
| S1 Lung | **77** | 77 | 45 | **45** |
| S1 Lung minus NCI-H226 | 76 | 76 | 44 | 44 |
| Oncotree Lung | 76 | 76 | 44 | 44 |
| Oncotree NSCLC | **63** | 63 | 35 | **35** |
| Oncotree LUAD | 37 | 37 | 22 | 22 |
| Oncotree LUSC | 11 | 11 | 7 | 7 |
| Oncotree NET / SCLC | 10 | 10 | 7 | 7 |
| S1 solid organ | 335 | 335 | 186 | 186 |

Closest line counts to 118 are 77 (all S1 lung) and 63 (Oncotree NSCLC). Closest pairwise counts are 45 and 35. **None is ≈118** (we pre-declared |n−118| ≤ 5 as “approximately”).

NCI-H226 is `Tissue of Origin = Lung` in Table S1 / CCLE 2012 naming (`NCIH226_LUNG`) but `OncotreeLineage = Pleura` (epithelioid mesothelioma) in DepMap 24Q4. Dropping it changes 77→76, not toward 118.

Arbitrary unions of 2–3 S1 tissues can be forced to 118 lines (e.g. Lung + Breast + Kidney). Those mash-ups are **not** a published NSCLC definition and were not used for the ρ claim. See `n118_tissue_combo_hunt.csv`.

### Where n=118 actually comes from

Augustyn et al., *PNAS* 2014 (doi:10.1073/pnas.1410419111):

> genome-wide **mRNA** expression data collected from 206 human lung cell lines, including **NSCLC (n = 118)**, SCLC (n = 29), and … HBECs/HSAECs (n = 59)

The same RNA n=118 was reused later (e.g. Nilsson et al., *Cancer Cell* 2023, CD70 / EMT in NSCLC cell lines). That is Affymetrix / CCLE expression, **not** Nusinow MS protein.

Other public protein matrices checked previously and not an n=118 TACSTD2–CLDN4 table:

- ProCan–DepMapSanger DIA-MS (Gonçalves 2022, 949 lines): CLDN4 absent
- CCLE RPPA500 / TCPA: TACSTD2 and CLDN4 not on the antibody panel

---

## Recompute on the exact table (not on n=118)

Pre-specified: Spearman primary (`scipy.stats.spearmanr`), Pearson secondary, two-sided, pairwise complete cases. Unique CCLE codes; the three plex-duplicate lines (CAL120, SW948, HCT15) are averaged. 5,000-resample bootstrap 95% CI, seed 0. **No filter was tuned to hit 0.69 or 118.**

**S1 Lung, pairwise complete (figure `fig_scatter_S1_lung.png`)**

- n = **45**
- Spearman ρ = **0.6931** (p = 1.31×10⁻⁷) → **0.69** at 2 d.p.
- Pearson r = 0.687
- Bootstrap 95% CI [0.490, 0.824]; 0.69 is inside the CI
- This is the slice that reproduces the user’s **ρ**, on the wrong **n**

**Oncotree NSCLC, pairwise complete (figure `fig_scatter_Oncotree_NSCLC.png`)**

- n = **35**
- Spearman ρ = **0.7266** (p = 7.73×10⁻⁷) → **0.73**
- Pearson r = 0.744
- Nearby to 0.69 (|Δ|=0.037) but not a 2-d.p. match; n is not 118

CLDN4 missingness is the reason pairwise n is far below line n: 32/77 S1 lung lines have TACSTD2 but no CLDN4 protein value in the Gygi table.

---

## Methods (pre-specified)

1. Download the Gygi TSV (official “Protein Quantitation (TSV Format)” of Table S2) and Table S1. Do not use a derived portal extract as the primary matrix.
2. Identify sample columns by the `_TenPxNN` suffix. Collapse the three duplicated CCLE codes by mean.
3. Take the single `Gene_Symbol == TACSTD2` row and the single `Gene_Symbol == CLDN4` row.
4. Define lung / NSCLC from Table S1 and from DepMap 24Q4 `Model.csv` joined on `CCLEName`. Report both. Do not invent a third histology list to hit 118.
5. Correlate with pairwise complete cases only. Do not impute CLDN4 NAs to inflate n toward 118.

---

## Reproduce

```bash
python3 -m pip install -r scripts/rework/B2_n118/requirements.txt
python3 scripts/rework/B2_n118/download.py --outdir data/B2_n118
python3 scripts/rework/B2_n118/analyze.py --data data/B2_n118 --outdir results/rework/B2_n118
```

Machine verdict: `verdict.txt`. Full numbers: `summary.json`, `correlations.csv`.

---

## Files

| File | Role |
|---|---|
| `verdict.txt` | One-screen verdict |
| `summary.json` | Machine-readable claim, table IDs, all cohorts |
| `correlations.csv` | Every natural cohort |
| `S1_lung_protein_TACSTD2_CLDN4.csv` | 77 S1 lung lines, protein values + labels |
| `Oncotree_NSCLC_protein_TACSTD2_CLDN4.csv` | 63 Oncotree NSCLC lines |
| `all_375_protein_TACSTD2_CLDN4.csv` | Full unique-line extract |
| `n118_tissue_combo_hunt.csv` | Tissue unions with n near 118 (not used) |
| `download_manifest.json` | URLs, bytes, SHA256 |
| `fig_scatter_S1_lung.png` | Lung pairwise scatter (ρ=0.69, n=45) |
| `fig_scatter_Oncotree_NSCLC.png` | NSCLC pairwise scatter (ρ=0.73, n=35) |
| `fig_n_vs_118.png` | Line/pairwise n vs the claimed 118 |
