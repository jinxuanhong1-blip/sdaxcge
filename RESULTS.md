# CosMx NSCLC: is CLDN4 in a tight-junction module that anti-correlates with CD8?

Official He et al. 2022 CosMx SMI 960-plex (NanoString public flat files). 766,334 QC cells, 232 FOVs, 8 sections, 5 patients. QC and the RNA-epithelial / CD8 definitions match the contact analysis (225 FOVs have at least 40 epithelial cells and 5 CD8 T cells). This does not replace the locked 20 µm contact odds ratio or the 50/100 µm cytotoxic ratios.

## Answer

CLDN4 is not in a classical tight-junction module on this panel. It is the only claudin present. OCLN, TJP1, TJP2, TJP3, F11R, CLDN1/3/5/7/18, CGN, MARVELD2/3, CRB3, PARD3, JAM2, and JAM3 are absent. CDH1 is an adherens-junction gene. EPCAM and TACSTD2 are epithelial adhesion genes.

CLDN4’s spatial partners are an epithelial program that mixes those adhesion genes with keratins. The strongest Hotspot local-correlation partners (median Z across 232 FOVs) are KRT19 (11.14), EPCAM (10.81), TACSTD2 (10.81), KRT8 (10.16), EZR (9.79), and KRT7 (9.35). CDH1 is weaker (4.90). KRT18 is 5.93. Under Hotspot’s default program cut (≥20 genes, core genes only), CLDN4 falls in a module with at least two adhesion partners and at least two keratins in 121/232 FOVs. A keratin-free adhesion module occurs in 2/232 FOVs.

That epithelial field anti-correlates with the CD8 neighborhood, modestly for CLDN4 itself and more strongly for the module score. Among RNA-epithelial cells, Spearman correlation of CLDN4 with the count of CD8 T cells inside 50 µm has a five-patient mean of **−0.065** (median FOV −0.071; 225 FOVs). Four of five patients and seven of eight sections are negative. Lung6 is positive (+0.022). A toroidal shift of the CD8 pattern gives one-sided p = 0.005 (none of 199 shifts ≤ the observed patient mean). A patient sign test does not reject at n = 5 (4/5, p = 0.19). The Hotspot fine-module score, on the 209 FOVs where CLDN4 is assigned, has patient-mean **−0.248** (4/5 patients; Lung6 +0.029; same shift p = 0.005).

The association is the epithelial program, not a CLDN4 residual past keratins. A KRT7/8/18/19 score is negative in 5/5 patients and 8/8 sections (patient-mean −0.110). An on-panel adhesion score (CLDN4, CDH1, EPCAM, TACSTD2) is −0.100, also 5/5 and 8/8. After rank-residualizing CLDN4 on the keratin score, the patient-mean is **−0.017**. Keratin and the adhesion score carry the consistent sign; CLDN4 beyond keratin does not add a meaningful neighborhood effect.

Same-cell CLDN4 versus CD8A/CD8B transcript correlation is slightly positive (median FOV ρ +0.038; every patient mean is positive). The neighborhood result is a spatial count association, not same-cell anti-correlation. Local-correlation Z with CD8A is only −2.22. The strongest negative partners of CLDN4 are stromal genes (IGFBP7 −7.56, BGN −6.98, COL3A1 −6.72, COL1A1 −6.49), which is epithelial versus stroma geography.

## Module calls

Adhesion-like means the module contains CLDN4 and at least two of {CDH1, EPCAM, TACSTD2}. Keratin-mixed means it also contains at least two of {KRT7, KRT8, KRT18, KRT19}. Counts are FOVs (232 with ≥80 QC cells).

| Call | Hotspot program (≥20, core) | Hotspot fine (≥4) | HMRF fine (Pearson ≥ 0.40, ≥4) |
|---|---:|---:|---:|
| Adhesion and keratin | 121 | 35 | 54 |
| Adhesion-like, keratin-free | 2 | 14 | 2 |
| Keratin without adhesion | 37 | 55 | 30 |
| Other module | 34 | 111 | 23 |
| CLDN4 unassigned | 38 | 17 | 123 |

Hotspot program co-membership, among FOVs where CLDN4 was assigned and the partner passed the autocorrelation filter: TACSTD2 159/190, KRT19 154/193, EZR 152/185, KRT7 147/194, EPCAM 141/187, CDH1 106/142, KRT8 139/193. CD8A, CD8B, PTPRC, and CD3D share that module in 0 FOVs. ESAM, an endothelial junction gene, shares it in 1/130 eligible FOVs.

The fine module’s most frequent partners, out of 215 assigned FOVs, are TACSTD2 (107), EZR (97), KRT7 (86), and EPCAM (85). CDH1 joins that fine module in 42/157 eligible FOVs. The fine cut was included so a four-gene adhesion set could be called; it usually is not. The stable object is the larger epithelial program.

HMRF (2-state Potts, symmetrized 8-NN, β = 0.5) agrees when it assigns a module: 54/232 FOVs are adhesion-and-keratin, and a keratin-free adhesion module occurs twice. CLDN4 is unassigned in 123/232 FOVs, including every Lung6 FOV, so the five-patient HMRF score test is undefined. Median HMRF high-state correlation with CLDN4 is 0.39 for EPCAM, 0.35 for KRT19, 0.34 for KRT8, 0.27 for TACSTD2, and 0.20 for CDH1. A Pearson cut of 0.50 sits above the 99th percentile of these sparse fields and forms no modules; 0.40 was set from that scale check on one FOV before the cohort run.

## CD8 neighborhood field

Primary endpoint: epithelial CLDN4 (library-size log1p) versus CD8 T cell count within 50 µm (0.18 µm/pixel). CD8 T cells are CD8A or CD8B positive, CD3D/E/G positive, and not RNA-epithelial. The statistic is the unweighted mean of five patient means. The null shifts CD8 coordinates on a torus inside the FOV box (199 shifts, seed 20260921).

| Patient | FOVs | CLDN4 ρ | Hotspot fine-module ρ | Keratin ρ |
|---|---:|---:|---:|---:|
| Lung5 | 89 | −0.152 | −0.458 | −0.259 |
| Lung6 | 23 | +0.022 | +0.029 | −0.027 |
| Lung9 | 65 | −0.042 | −0.184 | −0.063 |
| Lung12 | 28 | −0.121 | −0.386 | −0.188 |
| Lung13 | 20 | −0.034 | −0.241 | −0.014 |
| Patient mean |  | −0.065 | −0.248 | −0.110 |

Module FOV counts are smaller where CLDN4 was unassigned (Lung5 80, Lung6 20, Lung9 64, Lung12 25, Lung13 20). At 100 µm, CLDN4’s patient-mean is −0.041 (3/5 patients; Lung6 and Lung13 positive; shift p = 0.005).

A 5-domain Potts HMRF on 10 expression PCs puts the highest-CLDN4 domain at median epithelial purity 0.95. Its CD8 fraction is lower than the other domains by a median of 0.030 (215/225 spatial FOVs). That contrast is the composition of an epithelial domain. Lung6 and Lung9 mean deltas are about −0.007.

## What this does not say

The continuous CLDN4 correlation is small, misses Lung6, and is weaker than the keratin score. It is a different statistic from the locked cytotoxic ratios and the contact odds ratio. Same-cell transcript ranks are not an exclusion result. Visium same-spot correlations are not part of this analysis. No private 8-KL data were used.

## Files

- Script: `scripts/cosmx_hotspot_hmrf_modules.py`
- FOV table: `results/cosmx_hotspot_hmrf/tables/fov_metrics.tsv`
- Partners: `results/cosmx_hotspot_hmrf/tables/cldn4_spatial_coupling.tsv`
- Members: `results/cosmx_hotspot_hmrf/tables/cldn4_module_members.tsv`
- Figures: `results/cosmx_hotspot_hmrf/figures/`
