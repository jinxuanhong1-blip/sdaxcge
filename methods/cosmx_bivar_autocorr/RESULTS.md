# RESULTS — CosMx bivariate Moran / Lee / Gi* of CLDN4 vs CD8A

**Dataset.** Figshare 25976224, `cosmx_human_nsclc_clustered.h5ad` (He et al. 2022 CosMx NSCLC; CellCharter object).
Cells in the object: **765771**. Sample–FOV units analyzed: **232**. Sections: **8**. Patients: **5**.
Expression: log1p(counts / n_counts × 10,000) from `layers/counts`. Coordinates: `obsm['spatial']` global pixels, 0.18 µm/px.
No private 8-KL data. CLDN4 only (no TACSTD2 gate).

## Estimand

Primary: partial bivariate Moran's I, CLDN4 with the row-standardized lag of CD8A, both genes residualized by within-FOV OLS on KRT8 and EPCAM. Graph: cells within 50 µm. Co-primary graph: kNN k=15.
Lee's L uses the same residualized fields. Gi* uses binary star weights on the same neighbor sets.
Within-tumor cross-correlation residualizes CLDN4 on KRT8 and EPCAM inside the tumor-cell set (author labels tumor 5/6/9/12/13) and correlates that residual with the lag of CD8A. The null shuffles the residual among tumor cells.
KRT8–CD8A and EPCAM–CD8A Moran use the same CD8A permutations as raw CLDN4–CD8A Moran and are the epithelial baselines.

Permutation p-values use 199 within-unit shuffles. Section Wilcoxon and the patient sign test are the confirmatory summaries. FOV inverse-variance p-values assume independent FOVs and are reported as descriptive effect sizes (RE = DerSimonian–Laird). I² in the tables is a fraction from 0 to 1.

## Headline

Primary partial Moran I (radius 50 µm, both genes residualized on KRT8+EPCAM): median FOV estimate +0.0002 (112/232 FOVs negative; FOV Wilcoxon p=0.5494). Section medians negative in 3/8 (Wilcoxon p=0.7266). Patients negative in 2/5 (one-sided sign p=0.8125).

Raw Moran I on the same graph (no residualization): median FOV estimate +0.0019 (99/232 FOVs negative; FOV Wilcoxon p=0.9995). Section medians negative in 4/8 (Wilcoxon p=0.8086). Patients negative in 2/5 (one-sided sign p=0.8125).

Within-tumor partial correlation on the same graph: median FOV estimate +0.0196 (80/217 FOVs negative; FOV Wilcoxon p=1.0000). Section medians negative in 0/8 (Wilcoxon p=1.0000). Patients negative in 0/5 (one-sided sign p=1.0000).

On the primary 50 µm graph the partial Moran I is centered at zero (median +0.0002, FOV two-sided Wilcoxon p=0.9012). Raw Moran I is a very small positive (median +0.0019, FOV two-sided Wilcoxon p=0.0011). That FOV p-value treats FOVs as independent. Section medians of raw I are split evenly, and the patient sign test is not negative. The raw shift has the same direction as the KRT8 and EPCAM baselines. Co-primary kNN k=15 partial I is also centered at zero (see the graph table). No graph gives a patient-level negative sign test.

Within tumor cells, CLDN4 and the CD8A lag are weakly positively correlated (median r +0.0196, FOV two-sided Wilcoxon p=4.22e-06; section medians are positive in every section that has enough tumor cells). The top-minus-bottom CLDN4 quartile difference in mean neighbor CD8A is a small positive, not a deficit. CD8A is detected in a similar fraction of tumor cells and other cells, so this transcript field is not a CD8 T-cell mask. The locked result counted cytotoxic cell types around CLDN4-high tumor cells. That count is a different estimand and is not recomputed here.

## Partial Moran I across neighbor graphs

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 232 | 110 | +0.0005 | 2/8 | 0.7695 | 2/5 | 0.8125 | +0.0007 | +0.117 |
| knn_k15 | 232 | 123 | -0.0005 | 6/8 | 0.1250 | 4/5 | 0.1875 | -0.0000 | +0.174 |
| knn_k30 | 232 | 121 | -0.0003 | 5/8 | 0.3711 | 2/5 | 0.8125 | +0.0000 | +0.156 |
| radius_20um | 232 | 118 | -0.0001 | 5/8 | 0.5273 | 3/5 | 0.5000 | +0.0008 | +0.175 |
| radius_50um | 232 | 112 | +0.0002 | 3/8 | 0.7266 | 2/5 | 0.8125 | +0.0001 | +0.127 |
| radius_100um | 232 | 118 | -0.0001 | 3/8 | 0.5273 | 1/5 | 0.9688 | -0.0001 | +0.254 |
| delaunay | 232 | 116 | +0.0000 | 3/8 | 0.7266 | 2/5 | 0.8125 | +0.0009 | +0.149 |

## Raw Moran I (CLDN4, lag CD8A)

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 232 | 86 | +0.0029 | 3/8 | 0.8750 | 2/5 | 0.8125 | +0.0045 | +0.539 |
| knn_k15 | 232 | 96 | +0.0017 | 3/8 | 0.8086 | 2/5 | 0.8125 | +0.0030 | +0.562 |
| knn_k30 | 232 | 96 | +0.0018 | 3/8 | 0.8438 | 1/5 | 0.9688 | +0.0022 | +0.559 |
| radius_20um | 232 | 99 | +0.0038 | 4/8 | 0.8086 | 2/5 | 0.8125 | +0.0047 | +0.533 |
| radius_50um | 232 | 99 | +0.0019 | 4/8 | 0.8086 | 2/5 | 0.8125 | +0.0023 | +0.540 |
| radius_100um | 232 | 104 | +0.0008 | 3/8 | 0.8086 | 1/5 | 0.9688 | +0.0009 | +0.573 |
| delaunay | 232 | 84 | +0.0033 | 4/8 | 0.8438 | 3/5 | 0.5000 | +0.0048 | +0.526 |

## Lee's L, partial

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 232 | 101 | +0.0011 | 3/8 | 0.6797 | 2/5 | 0.8125 | +0.0012 | +0.222 |
| knn_k15 | 232 | 107 | +0.0004 | 2/8 | 0.9023 | 1/5 | 0.9688 | +0.0005 | +0.242 |
| knn_k30 | 232 | 114 | +0.0001 | 4/8 | 0.5781 | 2/5 | 0.8125 | +0.0001 | +0.256 |
| radius_20um | 232 | 92 | +0.0019 | 1/8 | 0.9609 | 1/5 | 0.9688 | +0.0019 | +0.294 |
| radius_50um | 232 | 115 | +0.0001 | 4/8 | 0.7266 | 2/5 | 0.8125 | +0.0002 | +0.246 |
| radius_100um | 232 | 114 | +0.0001 | 3/8 | 0.5273 | 1/5 | 0.9688 | +0.0000 | +0.337 |
| delaunay | 232 | 95 | +0.0013 | 2/8 | 0.9609 | 1/5 | 0.9688 | +0.0016 | +0.234 |

## CLDN4 residual only (CD8A not residualized)

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 232 | 89 | +0.0034 | 2/8 | 0.9727 | 2/5 | 0.8125 | +0.0032 | +0.365 |
| knn_k15 | 232 | 96 | +0.0018 | 3/8 | 0.8750 | 2/5 | 0.8125 | +0.0020 | +0.440 |
| knn_k30 | 232 | 92 | +0.0013 | 3/8 | 0.8750 | 2/5 | 0.8125 | +0.0016 | +0.437 |
| radius_20um | 232 | 93 | +0.0024 | 1/8 | 0.9805 | 0/5 | 1.0000 | +0.0031 | +0.377 |
| radius_50um | 232 | 93 | +0.0015 | 2/8 | 0.9609 | 2/5 | 0.8125 | +0.0016 | +0.425 |
| radius_100um | 232 | 102 | +0.0004 | 4/8 | 0.8438 | 2/5 | 0.8125 | +0.0008 | +0.466 |
| delaunay | 232 | 91 | +0.0029 | 3/8 | 0.9609 | 2/5 | 0.8125 | +0.0032 | +0.355 |

## Epithelial baselines (raw Moran)

KRT8 with lag CD8A:

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 232 | 90 | +0.0055 | 3/8 | 0.9453 | 1/5 | 0.9688 | +0.0060 | +0.585 |
| knn_k15 | 232 | 92 | +0.0035 | 3/8 | 0.9023 | 1/5 | 0.9688 | +0.0041 | +0.601 |
| knn_k30 | 232 | 89 | +0.0035 | 3/8 | 0.9258 | 1/5 | 0.9688 | +0.0035 | +0.606 |
| radius_20um | 232 | 90 | +0.0048 | 3/8 | 0.9453 | 1/5 | 0.9688 | +0.0057 | +0.537 |
| radius_50um | 232 | 89 | +0.0032 | 3/8 | 0.8750 | 1/5 | 0.9688 | +0.0034 | +0.618 |
| radius_100um | 232 | 96 | +0.0013 | 3/8 | 0.8438 | 1/5 | 0.9688 | +0.0016 | +0.600 |
| delaunay | 232 | 87 | +0.0053 | 3/8 | 0.9453 | 1/5 | 0.9688 | +0.0065 | +0.568 |

EPCAM with lag CD8A:

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 232 | 98 | +0.0027 | 4/8 | 0.8750 | 2/5 | 0.8125 | +0.0039 | +0.586 |
| knn_k15 | 232 | 102 | +0.0018 | 4/8 | 0.8438 | 2/5 | 0.8125 | +0.0029 | +0.622 |
| knn_k30 | 232 | 103 | +0.0016 | 3/8 | 0.7695 | 1/5 | 0.9688 | +0.0023 | +0.645 |
| radius_20um | 232 | 97 | +0.0033 | 3/8 | 0.9023 | 1/5 | 0.9688 | +0.0044 | +0.559 |
| radius_50um | 232 | 101 | +0.0018 | 3/8 | 0.7695 | 1/5 | 0.9688 | +0.0023 | +0.640 |
| radius_100um | 232 | 109 | +0.0006 | 3/8 | 0.7266 | 1/5 | 0.9688 | +0.0009 | +0.622 |
| delaunay | 232 | 96 | +0.0024 | 3/8 | 0.9258 | 1/5 | 0.9688 | +0.0043 | +0.566 |

## Within-tumor CLDN4 quartile vs neighbor CD8A

Among tumor cells with at least one neighbor, mean neighbor CD8A in the top CLDN4 quartile minus the bottom quartile. Negative means CLDN4-high tumor cells have less CD8A around them. `tumor_q_diff_partial` uses CLDN4 residualized on KRT8 and EPCAM inside the tumor set.

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 216 | 94 | +0.0033 | 2/8 | 0.9023 | 1/5 | 0.9688 | +0.0045 | +0.325 |
| knn_k15 | 216 | 87 | +0.0038 | 1/8 | 0.9961 | 1/5 | 0.9688 | +0.0041 | +0.514 |
| knn_k30 | 216 | 80 | +0.0049 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0049 | +0.631 |
| radius_20um | 216 | 85 | +0.0056 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0048 | +0.324 |
| radius_50um | 216 | 77 | +0.0035 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0051 | +0.650 |
| radius_100um | 216 | 74 | +0.0020 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0031 | +0.772 |
| delaunay | 216 | 87 | +0.0075 | 3/8 | 0.8750 | 1/5 | 0.9688 | +0.0043 | +0.304 |

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 216 | 102 | +0.0016 | 3/8 | 0.7695 | 2/5 | 0.8125 | +0.0043 | +0.300 |
| knn_k15 | 216 | 99 | +0.0045 | 1/8 | 0.9727 | 1/5 | 0.9688 | +0.0036 | +0.491 |
| knn_k30 | 216 | 91 | +0.0037 | 1/8 | 0.9805 | 1/5 | 0.9688 | +0.0039 | +0.574 |
| radius_20um | 216 | 102 | +0.0038 | 3/8 | 0.8750 | 3/5 | 0.5000 | +0.0039 | +0.313 |
| radius_50um | 216 | 96 | +0.0030 | 1/8 | 0.9883 | 1/5 | 0.9688 | +0.0037 | +0.589 |
| radius_100um | 216 | 92 | +0.0029 | 2/8 | 0.9805 | 1/5 | 0.9688 | +0.0028 | +0.721 |
| delaunay | 216 | 98 | +0.0043 | 2/8 | 0.8438 | 2/5 | 0.8125 | +0.0034 | +0.209 |

CD8A detection (expression > 0), median across primary-radius FOVs: **0.082** of tumor cells and **0.070** of other cells. Median fraction of cells positive for both CLDN4 and CD8A: **0.026**.

## Within-tumor partial cross-correlation

Pearson r of tumor-cell CLDN4 with the CD8A lag. This is not a Moran I.

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 217 | 92 | +0.0071 | 2/8 | 0.9883 | 1/5 | 0.9688 | +0.0087 | +0.237 |
| knn_k15 | 217 | 85 | +0.0120 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0116 | +0.428 |
| knn_k30 | 217 | 79 | +0.0166 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0162 | +0.521 |
| radius_20um | 217 | 90 | +0.0058 | 1/8 | 0.9922 | 0/5 | 1.0000 | +0.0094 | +0.257 |
| radius_50um | 217 | 80 | +0.0196 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0179 | +0.545 |
| radius_100um | 218 | 70 | +0.0246 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0219 | +0.719 |
| delaunay | 217 | 88 | +0.0063 | 2/8 | 0.9727 | 1/5 | 0.9688 | +0.0086 | +0.201 |

Within-tumor raw CLDN4 (no residualization) is in `tables/fov_graph_stats.tsv` (`tumor_raw_r`). The within-tumor number is a Pearson correlation, not a Moran I.

## Primary graph by section and patient

Values are medians of FOV-level estimates.

| unit | partial_I | raw_I | partial_L | tumor_partial_r | krt8_I |
| --- | --- | --- | --- | --- | --- |
| LUAD-12 | +0.0012 | +0.0089 | +0.0017 | +0.0247 | +0.0105 |
| LUAD-13 | -0.0008 | -0.0007 | -0.0004 | +0.0271 | +0.0039 |
| LUAD-5 R1 | -0.0000 | -0.0009 | -0.0007 | +0.0320 | -0.0050 |
| LUAD-5 R2 | -0.0010 | -0.0050 | -0.0008 | +0.0203 | -0.0077 |
| LUAD-5 R3 | +0.0004 | -0.0019 | -0.0009 | +0.0244 | -0.0022 |
| LUAD-9 R1 | +0.0001 | +0.0074 | +0.0002 | +0.0298 | +0.0070 |
| LUAD-9 R2 | +0.0005 | +0.0076 | +0.0009 | +0.0105 | +0.0091 |
| LUSC-6 | +0.0003 | +0.0016 | +0.0015 | +0.0102 | +0.0053 |

Patient means of those section medians:

| unit | partial_I | raw_I | partial_L | tumor_partial_r | krt8_I |
| --- | --- | --- | --- | --- | --- |
| Lung12 | +0.0012 | +0.0089 | +0.0017 | +0.0247 | +0.0105 |
| Lung13 | -0.0008 | -0.0007 | -0.0004 | +0.0271 | +0.0039 |
| Lung5 | -0.0002 | -0.0026 | -0.0008 | +0.0256 | -0.0050 |
| Lung6 | +0.0003 | +0.0016 | +0.0015 | +0.0102 | +0.0053 |
| Lung9 | +0.0003 | +0.0075 | +0.0006 | +0.0202 | +0.0080 |

On the primary radius, median within-FOV R² of CLDN4 on KRT8+EPCAM is **0.126**, and of CD8A on KRT8+EPCAM is **0.001**. Median same-cell Spearman CLDN4 vs CD8A is **+0.029** (raw) and **+0.089** (residuals). Same-cell Spearman is co-expression inside a cell. Moran's I is the neighbor lag.

## Gi* overlap

**radius_50um.** Residual Gi* (CLDN4 z>1.96 and CD8A z<−1.96): median overlap fraction 0.0000 versus independence 0.0006 (median OR 0.325). FOVs with any analytic CD8 cold cells: 187/232. Fisher p<0.05 in 30/232 FOVs. Wilcoxon on (overlap − expected), alternative greater: p=0.9937. Median minimum residual CD8 Gi* z=-2.26. Raw-expression rank cold (CD8 Gi* at or below the FOV 10th percentile) overlap median 0.0132 versus independence 0.0193 (Wilcoxon greater p=1.0000).

**knn_k15.** Residual Gi* (CLDN4 z>1.96 and CD8A z<−1.96): median overlap fraction 0.0000 versus independence 0.0000 (median OR omitted because analytic cold cells are rare). FOVs with any analytic CD8 cold cells: 6/232. Fisher p<0.05 in 3/232 FOVs. Wilcoxon on (overlap − expected), alternative greater: p=0.0579. Median minimum residual CD8 Gi* z=-1.30. Raw-expression rank cold (CD8 Gi* at or below the FOV 10th percentile) overlap median 0.0323 versus independence 0.0370 (Wilcoxon greater p=1.0000).

Analytic CD8 cold spots can be scarce when CD8A is sparse, because the local mean cannot fall far below a near-zero background. The rank-based cold set is the check for that limit. Overlap is a hotspot coincidence, not a contact probability.

## Whole-section graphs

The same statistics on one graph per tissue section (global coordinates, so neighbors can cross FOV borders). Each row is a section, not a FOV.

Partial Moran I:

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 8 | 4 | +0.0001 | 4/8 | 0.8086 | 2/5 | 0.8125 | +0.0003 | +0.424 |
| knn_k15 | 8 | 5 | -0.0010 | 5/8 | 0.3203 | 4/5 | 0.1875 | -0.0004 | +0.532 |
| knn_k30 | 8 | 4 | -0.0004 | 4/8 | 0.3203 | 3/5 | 0.5000 | -0.0003 | +0.500 |
| radius_20um | 8 | 5 | -0.0006 | 5/8 | 0.6797 | 2/5 | 0.8125 | +0.0002 | +0.551 |
| radius_50um | 8 | 4 | -0.0000 | 4/8 | 0.4727 | 4/5 | 0.1875 | -0.0002 | +0.565 |
| radius_100um | 8 | 4 | -0.0001 | 4/8 | 0.3203 | 2/5 | 0.8125 | -0.0005 | +0.764 |
| delaunay | 8 | 4 | +0.0005 | 4/8 | 0.6797 | 2/5 | 0.8125 | +0.0003 | +0.539 |

Raw Moran I:

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 8 | 3 | +0.0025 | 3/8 | 0.8086 | 1/5 | 0.9688 | +0.0044 | +0.954 |
| knn_k15 | 8 | 3 | +0.0026 | 3/8 | 0.8086 | 1/5 | 0.9688 | +0.0031 | +0.963 |
| knn_k30 | 8 | 3 | +0.0031 | 3/8 | 0.8086 | 1/5 | 0.9688 | +0.0026 | +0.964 |
| radius_20um | 8 | 3 | +0.0033 | 3/8 | 0.8750 | 1/5 | 0.9688 | +0.0046 | +0.954 |
| radius_50um | 8 | 3 | +0.0033 | 3/8 | 0.8086 | 1/5 | 0.9688 | +0.0026 | +0.965 |
| radius_100um | 8 | 3 | +0.0033 | 3/8 | 0.6289 | 1/5 | 0.9688 | +0.0009 | +0.969 |
| delaunay | 8 | 3 | +0.0025 | 3/8 | 0.8438 | 1/5 | 0.9688 | +0.0048 | +0.955 |

Within-tumor partial Pearson r:

| graph | n | n_neg | median | section_neg | section_p | patient_neg | patient_sign_p | RE | I2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| knn_k6 | 8 | 0 | +0.0124 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0149 | +0.737 |
| knn_k15 | 8 | 0 | +0.0176 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0217 | +0.847 |
| knn_k30 | 8 | 0 | +0.0250 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0308 | +0.909 |
| radius_20um | 8 | 0 | +0.0141 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0140 | +0.460 |
| radius_50um | 8 | 0 | +0.0284 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0332 | +0.917 |
| radius_100um | 8 | 0 | +0.0373 | 0/8 | 1.0000 | 0/5 | 1.0000 | +0.0448 | +0.951 |
| delaunay | 8 | 1 | +0.0146 | 1/8 | 0.9961 | 0/5 | 1.0000 | +0.0144 | +0.697 |

On the whole-section 50 µm graph, partial Moran I section-as-unit: median FOV estimate -0.0000 (4/8 FOVs negative; FOV Wilcoxon p=0.4727). Section medians negative in 4/8 (Wilcoxon p=0.4727). Patients negative in 4/5 (one-sided sign p=0.1875). Because each section is a single test, the FOV count in that sentence is the section count.

## What this layer is

An expression-field autocorrelation layer on the public He et al. 2022 CosMx cohort. The locked cell-type result (CLDN4-high tumor neighborhoods contain fewer cytotoxic cells at 50/100 µm; exclusion without loss of GZMB/PRF1/NKG7/IFNG) is a different estimand and is not recomputed here.
Lung6 is the LUSC section. Lung5 and Lung9 contribute serial sections; patient summaries average those sections so technical replicates do not count as extra patients.

## Reproduce

```bash
bash methods/cosmx_bivar_autocorr/download.sh
python3 methods/cosmx_bivar_autocorr/analyze.py --h5ad data/cosmx_human_nsclc_clustered.h5ad
```

Outputs: `methods/cosmx_bivar_autocorr/tables/` and `methods/cosmx_bivar_autocorr/figures/`.
