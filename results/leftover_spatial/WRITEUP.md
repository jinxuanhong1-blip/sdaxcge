# Leftover public lung spatial series: CLDN4 / TACSTD2 vs immune neighborhood

Additive slice only. User B6 is taken as given. This folder does **not** re-run GSE221322, the 10x FFPE Visium demo (PR #126), E-MTAB-13530 (PR #134), or the PR #164 CosMx leftover.

Question: in **other** 2023–2026 public GeoMx / Visium / CosMx lung series that actually measure TACSTD2 or CLDN4, how do those genes sit relative to a T/B neighborhood?

Every ρ below is computed. n and p are reported. Series that lacked both genes after a real matrix check were skipped.

---

## English

### Hunt

GEO e-utils returned 110 human lung spatial-ish series (2023–2026). ArrayExpress leftover of interest was E-MTAB-14560/14566 (Visium + GeoMx lung/breast/DLBCL); public processed files there are images/FASTQ, not count matrices.

`tables/hunt_catalog.tsv` is the decision log. Analyzed leftovers (both genes present):

| Accession | Platform | Tissue | Spatial unit |
|---|---|---|---|
| [GSE263196](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263196) | Visium | SCLC, 5 sections | hex rings |
| [GSE273378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273378) | Visium | stage I LUAD, 16 sections | hex rings |
| [GSE265899](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE265899) | GeoMx WTA | lung cancer, paired tumor/immune AOIs | same-AOI + pair |
| [GSE289483](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289483) | GeoMx WTA | pulmonary pleomorphic carcinoma, 9 patients | same-AOI + CD45± |
| [GSE334014](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334014) | GeoMx WTA | asthma / control airway | PanCK vs stroma |
| [GSE326968](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326968) | GeoMx WTA | chronic lung allograft dysfunction | same-AOI |
| [GSE276083](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276083) | CosMx 979-plex | mycobacterial lung | FOV (no public x/y) |

Real gene-check skips: GSE261345 CTA (TACSTD2 yes, **CLDN4 no**); GSE186213 COVID CTA (TACSTD2 yes, **CLDN4 no**); GSE246011 LUAD-14C (both genes, **no coordinates**); GSE268850 IPF (both genes, almost no T-cell genes); GSE307534 / GSE288758 / GSE301973 / GSE299786 (archive-only or multi-GB).

### Methods (short)

Signatures locked in `scripts/leftover_spatial/common.py`. Visium: in-tissue, ≥200 genes, log1p(CP10K), hex rings 1/2/3, partial Spearman controlling for EPCAM/KRT/CDH1. GeoMx: depositor Q3 matrices. CosMx: signature rows only; public `spatial_coordinates` file is Width/Height, so the neighborhood is FOV-mean. Full text: `methods/leftover_spatial/playbook.md`.

### Results

**Visium GSE273378 (stage I LUAD, 16 sections).** Same-spot CLDN4 vs T+B: median ρ = **−0.095**, 13/16 negative, Wilcoxon p = **0.0076**. TACSTD2 vs T+B: median ρ = **−0.094**, 12/16 negative, p = **0.025**. Hex ring 1, CLDN4 vs neighbor T+B: median ρ = **−0.116**, 11/16 negative, p = **0.039**. After partialling the broad epithelial score, ring-1 median partial ρ = **−0.029**, p = **0.78**. The neighborhood signal is small and is not independent of epithelial content.

**Visium GSE263196 (SCLC, 5 sections).** Same-spot CLDN4 vs T+B: median ρ = **−0.059**, 4/5 negative, p = **0.625**. TACSTD2 vs T+B: median ρ = **+0.154**, 5/5 positive, p = **0.0625**. Ring 1 matches the same-spot signs. n = 5 is too small for a stable signed-rank test; the TACSTD2 direction is the opposite of T/B exclusion.

**GeoMx GSE265899 (PD-L1 lung cancer; 48 tumor + 47 immune AOIs).** Antibody-cut pairs: CLDN4 higher in tumor than immune in 46/47 pairs (median 103 vs 18, Wilcoxon p = 2.8×10⁻¹⁴); T+B higher in immune in 47/47 pairs (p = 1.4×10⁻¹⁴). All-AOI ρ(CLDN4, T+B) = **−0.821** (n = 95, p = 2.4×10⁻²⁴). **Within tumor AOIs** the same contrast is still negative: ρ = **−0.399** (n = 48, p = 0.0049); TACSTD2 vs T+B ρ = **−0.646** (n = 48, p = 7.2×10⁻⁷). Within immune AOIs CLDN4 vs T+B is not significant (ρ = −0.180, n = 47, p = 0.23).

**GeoMx GSE289483 (pleomorphic carcinoma, 114 AOIs / 9 patients).** All-AOI ρ(CLDN4, T+B) = **−0.480** (n = 114, p = 6.3×10⁻⁸); TACSTD2 vs T+B ρ = **−0.511** (p = 6.2×10⁻⁹). The adenocarcinoma component keeps the negative (CLDN4 ρ = **−0.518**, n = 39, p = 7.3×10⁻⁴). The **sarcomatoid** component does not (CLDN4 ρ = **+0.115**, n = 45, p = 0.45). CD45− AOIs remain negative (CLDN4 ρ = **−0.337**, n = 37, p = 0.041). CLDN4 level itself does not differ by CD45 segment (MW p = 0.56).

**GeoMx GSE334014 (asthma/control, 90 PanCK + 94 stroma).** All-AOI ρ(CLDN4, T+B) = **−0.683** (n = 184, p = 1.3×10⁻²⁶). Within PanCK: ρ = **−0.503** (n = 90, p = 4.5×10⁻⁷). Within stroma: ρ = **−0.071** (n = 94, p = 0.49). CLDN4 and TACSTD2 are higher in PanCK than stroma (MW p = 3.1×10⁻²³ and 2.0×10⁻²⁷).

**GeoMx GSE326968 (CLAD, 80 AOIs).** ρ(CLDN4, T+B) = **−0.419** (n = 80, p = 1.1×10⁻⁴); TACSTD2 vs T+B ρ = **−0.440** (p = 4.4×10⁻⁵). No public CD45 key.

**CosMx GSE276083 (mycobacterial lung; leftover CosMx, not PR #164).** 300,841 QC cells, CLDN4 detected in 7.0%, TACSTD2 in 4.2%. Same-cell ρ(CLDN4, T+B) = **+0.049** (n = 300,841, p = 4.3×10⁻¹⁶¹). FOV means (n = 47 FOVs, median 5,942 cells): ρ(CLDN4, T+B) = **+0.309** (p = 0.035); ρ(TACSTD2, T+B) = **+0.347** (p = 0.017); ρ(CLDN4, T) = **+0.482** (p = 6.1×10⁻⁴). Public coordinates are Width/Height only — no micron ball.

### Extra figures

- `figures/visium_samespot_rho.png` — per-section same-spot ρ
- `figures/visium_neighborhood_rings.png` — median hex-ring ρ and partial ρ
- `figures/geomx_target_vs_TB.png` — leftover GeoMx ρ
- `figures/cosmx_fov_scatter.png` — FOV-mean CLDN4/TACSTD2 vs T+B

### What the leftover numbers say

Tumor / epithelial GeoMx compartments and stage I LUAD Visium spots show a **negative** CLDN4/TACSTD2 vs T+B association. In LUAD Visium the hex-ring version shrinks to near zero after epithelial content is removed. SCLC Visium (n = 5) does not give a CLDN4 neighborhood, and TACSTD2 is weakly **positive**. CosMx FOVs in infected lung are **positive**. Pleomorphic **sarcomatoid** AOIs are null. These are extra public series, not a re-score of B6.

---

## 中文

只补 2023–2026 年尚未做过的公开肺空间系列。不重做 GSE221322、PR #126 Visium、E-MTAB-13530、PR #164 CosMx。矩阵里两个基因都没有的系列已跳过。

**Visium GSE273378（I 期肺腺癌，16 张切片）** 同点 CLDN4 对 T+B 中位 ρ = −0.095（13/16 为负，p = 0.0076）；六边形 1 环中位 ρ = −0.116（p = 0.039）。用宽上皮评分做偏相关后 1 环中位偏 ρ = −0.029（p = 0.78），邻域信号不独立于上皮含量。

**Visium GSE263196（SCLC，5 张）** CLDN4 同点中位 ρ = −0.059（p = 0.625）；TACSTD2 中位 ρ = +0.154（5/5 为正，p = 0.0625）。样本量只有 5，符号检验不稳定。

**GeoMx GSE265899** 肿瘤 AOI 的 CLDN4 高于配对免疫 AOI（46/47，p = 2.8×10⁻¹⁴）。全部 AOI ρ(CLDN4, T+B) = −0.821；**肿瘤内部**仍为负（ρ = −0.399，n = 48，p = 0.0049）。

**GeoMx GSE289483** 全部 AOI ρ = −0.480（n = 114）。腺癌成分 ρ = −0.518（n = 39，p = 7.3×10⁻⁴）；**肉瘤样成分不负**（ρ = +0.115，n = 45，p = 0.45）。

**GeoMx GSE334014** 全部 ρ = −0.683；PanCK 内 ρ = −0.503（n = 90）；间质内 ρ = −0.071（n = 94，p = 0.49）。

**CosMx GSE276083** 同细胞 ρ = +0.049；FOV 水平 ρ(CLDN4, T+B) = +0.309（n = 47，p = 0.035）。公开坐标文件没有 XY，不能做微米邻域。

目录与数字：`tables/`。复现见 `methods/leftover_spatial/playbook.md`。
