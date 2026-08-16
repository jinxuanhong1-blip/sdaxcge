# Methods — CPTAC LSCC protein vs xCell, purity residual

## Scope
Treatment-naive CPTAC LSCC/LUSC tumor **protein** only. Prespecified immune
endpoints are freeze xCell **CD8** (`xCell_T_cell_CD8+`) and **immune score**
(`xCell_immune_score`). This is a recompute of that PR23 slice plus a purity
residual. No ICI labels. No RNA-target tests. No LUAD.

## Data
Open CPTAC pan-cancer `data_freeze_v1.2_reorganized` on
`cptac-pancancer-data` (us-west-2). Protein tumor matrix, phenotype, and meta
only. Filenames match the LinkedOmics CPTAC-pancan-LSCC table / PR23.
Standalone `LSCC_xcell.txt` returns HTTP 403; xCell lives in `LSCC_phenotype.txt`.

Satpathy et al., *Cell* 2021, PMID 34358469: newly diagnosed resected LSCC,
no prior chemotherapy or radiotherapy. PDC proteome: PDC000234.

## Identifiers
- TACSTD2 / TROP2: `ENSG00000184292` (row `ENSG00000184292.7`)
- CLDN4: `ENSG00000189143` (row `ENSG00000189143.9`)
Match by Ensembl prefix. Pairwise-complete n is reported per test. CLDN4
protein NA=30/108 is TMT dropout.

## Purity covariates
Primary: **WES_purity** (DNA; 107/108). Sensitivity: **WGS_purity** (104/108)
and ESTIMATE TumorPurity = `cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)`
(Yoshihara et al. 2013). ESTIMATE purity is RNA-derived and circular with
immune content; it is not the honest residual.

## Statistics
- Marginal Spearman, pairwise-complete; n < 4 → NA
- Partial Spearman: rank-transform protein, xCell, and purity; OLS-residualize
  the protein and xCell ranks on the purity rank (intercept included); Pearson
  of those residuals. Two-sided t, df = n − 3. Fisher-z 95% CI uses n − 4.
- Literal purity residual: Spearman of OLS residuals of the *raw* protein and
  *raw* xCell values on raw purity (same complete cases).
- BH-FDR across the four prespecified primary partial-Spearman tests on WES
  purity. CD8 subsets and non-WES covariates are exploratory / sensitivity.

## How to rerun
```bash
pip install -r scripts/w200/CPTAC_LSCC_xcell/requirements.txt
python3 scripts/w200/CPTAC_LSCC_xcell/download.py
python3 scripts/w200/CPTAC_LSCC_xcell/analyze.py
```
