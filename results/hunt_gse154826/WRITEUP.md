# GSE154826 CITE-seq: TACSTD2 vs LCAM / CD8

**Verdict: the expected negative is not supported.**

User direction was negative (higher TACSTD2 → lower LCAM / CD8). In the authors’ own LCAM cohort this is null, and the LCAM point estimates are weakly *positive*. The only clean negative is a composition artifact: more epithelial leak in a CD45+ library → lower mean CD8 protein. That is not a Trop-2–LCAM biology result.

中文先看：**预期负相关不成立。** Leader 自己的 V2 beads 肿瘤 LCAM 队列（26 例）里，门控内 TACSTD2 对 LCAM 是 ρ=+0.25、p=0.22；对 CD8 Trm/T 是 ρ=−0.02、p=0.91。CITE-seq 上全细胞 TACSTD2 / 门控比例对 CD8 蛋白的负相关，是 CD45+ 文库里上皮漏检的成分效应，不是肿瘤 Trop-2 高则 LCAM/CD8 低。

---

## Question

Does TACSTD2 anti-correlate with the Leader LCAM score, or with CD8, in GSE154826?

## Data (what this dataset actually is)

| Item | Fact |
|---|---|
| Accession | [GSE154826](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154826) |
| Paper | Leader, Grout et al. *Cancer Cell* 2021; metadata + LCAM code from [effiken/Leader_et_al](https://github.com/effiken/Leader_et_al) |
| Cells used | 361,929 author-annotated barcodes (91 samples, 35 patients) |
| LCAM unit | patient × tissue, after the paper’s within-lineage normalisation |
| Primary cohort | V2 chemistry + CD45+ beads + Tumor, ≥80 immune cells: **n=26** (22 LUAD, 4 LUSC) |
| Digest / no CD45 | **2 patients** (695 LUAD, 706 LUSC). Correlation not computed. |
| Trop-2 / TACSTD2 ADT | **Not on the panel** |
| ADTs that exist | EPCAM, CD8, CD45, CD3, plus the rest of the hashed CITE panel |
| CD45 enrichment | 61 GEO libraries CD45+ beads, 10 CD45+ FACS. Epithelium was gated as `epi_endo_fibro_doublet`. |

This is an immune CITE-seq study. TACSTD2 RNA lives in the authors’ leftover epi/endo/fibro/doublet gate (mean 0.76 UMI/cell, 20.5% detect), not in T cells (0.009 UMI, 0.8% detect). Measuring “sample TACSTD2” on all annotated cells mostly measures how many epithelial barcodes leaked through CD45 enrichment.

LCAM is taken verbatim from `scripts/figure_5abcd_s5a.R`:

- hi = `T_activated` + `IgG` + `MoMac-II`
- lo = `B` + `AM` + `cDC2` + `AZU1_mac` + `Tcm/naive_II` + `cDC1`
- frequencies re-normalised inside T / B&plasma / MNP / lin_neg
- score = Σ log(hi + 0.01) − Σ log(lo + 0.01)

A second score uses the bulk-RNA gene list in `get_LCAM_scores.R` on immune-cell pseudobulk. Same conclusion.

TACSTD2 was scored three ways: all annotated cells; the epi/endo/fibro gate; epi-like gate cells (EPCAM/KRT UMI > endo and fibro). CD8 was scored as CD8 Trm / T, a broader CD8 set (Trm+GZMK+NKlike) / T, immune-cell CD8A CPM, and CD8 ADT where present.

## Primary result (paper LCAM set, n=26)

| x | y | ρ | p |
|---|---|---|---|
| gate TACSTD2 CPM | LCAM difference | **+0.25** | 0.22 |
| epi-like TACSTD2 CPM | LCAM difference | +0.08 | 0.68 |
| all-cell TACSTD2 CPM | LCAM difference | +0.02 | 0.94 |
| immune TACSTD2 CPM | LCAM difference | −0.08 | 0.71 |
| gate fraction | LCAM difference | −0.07 | 0.73 |
| gate TACSTD2 CPM | CD8 Trm / T | **−0.02** | 0.91 |
| all-cell TACSTD2 CPM | CD8 Trm / T | +0.02 | 0.92 |
| gate TACSTD2 CPM | broad CD8 / T | +0.15 | 0.46 |
| gate TACSTD2 CPM | immune CD8A CPM | +0.21 | 0.31 |
| gate TACSTD2 CPM | bulk-gene LCAM | +0.21 | 0.31 |

No test in the paper’s LCAM set is significant. LCAM point estimates are positive, opposite the requested sign. Residualising gate fraction does not create a negative (partial ρ for gate TACSTD2 vs LCAM = +0.25, p=0.22).

V2 beads LUAD only (n=22): gate TACSTD2 vs LCAM ρ=+0.25, p=0.27. LUSC beads n=4 — no Spearman.

## All tumor patients (n=35), including digest

Gate TACSTD2 vs LCAM: ρ=**+0.39**, p=0.019. Same sign as the primary, now nominally significant, still the **wrong sign** for the claim. BH q across the hunt table is 0.48 — do not treat this as a confirmed positive either. CD8 Trm stays null (ρ=−0.04, p=0.82).

The two digest patients (real epithelium, no CD45):

| patient | histology | gate TACSTD2 CPM | LCAM | CD8 Trm / T |
|---|---|---|---|---|
| 695 | LUAD | 1.61 | 10.67 | 0.052 |
| 706 | LUSC | 3.08 | 6.22 | 0.090 |

Two points. Not a correlation.

## CITE-seq protein (no Trop-2 antibody)

CD8 ADT is the only protein CD8 readout. EPCAM ADT exists on the full panel (mainly the two digest patients). There is no TACSTD2/Trop-2 ADT.

| cohort | x | y | n | ρ | p |
|---|---|---|---|---|---|
| patient Tumor | gate TACSTD2 | CD8 ADT (immune) | 8 | −0.14 | 0.74 |
| patient Tumor | all-cell TACSTD2 | CD8 ADT (immune) | 8 | −0.67 | 0.071 |
| patient Tumor | gate fraction | CD8 ADT (immune) | 8 | −0.69 | 0.058 |
| sample Tumor | gate TACSTD2 | CD8 ADT (immune) | 20 | −0.16 | 0.51 |
| sample Tumor | all-cell TACSTD2 | CD8 ADT (immune) | 20 | **−0.52** | 0.018 |
| sample Tumor | gate fraction | CD8 ADT (immune) | 20 | **−0.80** | 2.6×10⁻⁵ |

The negative CD8-protein signal tracks **how many non-immune barcodes are in the library**, not TACSTD2 inside those barcodes. Once TACSTD2 is measured in the gate, the CD8 correlation collapses. Cell-level TACSTD2 UMI vs CD8 ADT is ρ=−0.06 (n=72,463, p~10⁻⁶³) — different lineages, not a sample-level test.

## What this does not say

- It does not say Trop-2 protein is unrelated to LCAM. Trop-2 protein was not measured.
- It does not say a bulk-tumor TACSTD2–CD8 anti-correlation (seen in some CPTAC LUAD protein slices) is false. This experiment depleted CD45− cells.
- It does not say LCAM is poorly defined. LCAM is computed from the authors’ clusters and code.
- It does say: **in GSE154826, with the authors’ LCAM and the TACSTD2 that actually exists here, the expected negative is not observed.**

## Outputs

- `tables/key_stats.json`, `primary_focus.tsv`, `spearman.tsv`, `partial_spearman.tsv`
- `tables/patient_tissue_metrics.tsv`, `sample_metrics.tsv`, `digest_two_patients.tsv`
- `figures/fig1_primary_tacstd2_vs_lcam.png` … `fig7_cite_composition.png`
- scripts: `scripts/lib_gse154826.py`, `01_extract_batches.py`, `02_analyze.py`

---

## 中文

GSE154826 是 Leader 2021 的 NSCLC **CD45+ CITE-seq**，不是肿瘤上皮图谱。作者把上皮/内皮/成纤维/双胞打进 `epi_endo_fibro_doublet` 门。TACSTD2 几乎只在这扇门里（检测率 20.5%），T 细胞里接近零。CITE 面板没有 Trop-2 抗体，只有 EPCAM / CD8 / CD45。

LCAM 按论文脚本：T_activated + IgG + MoMac-II 对 B + AM + cDC2 + AZU1_mac + Tcm/naive_II + cDC1，先做谱系内归一。主队列是论文自己的 V2 beads 肿瘤，26 个病人。

主结果：门控 TACSTD2 对 LCAM ρ=+0.25（p=0.22），对 CD8 Trm/T ρ=−0.02（p=0.91）。不是负相关。把消化法两个病人加进去，门控 TACSTD2 对 LCAM 变成 ρ=+0.39（p=0.019），符号仍反。CITE 上全细胞 TACSTD2 / 门控比例对 CD8 蛋白的负相关（样本水平 ρ=−0.52 / −0.80）是上皮漏检的成分效应；门控内 TACSTD2 对 CD8 蛋白是 ρ=−0.16、p=0.51。

结论：在这份公开 CITE-seq 里，用作者的 LCAM 和实际测到的 TACSTD2，**用户预期的负号不成立**。不能用它支持 “TROP2 高则 LCAM/CD8 低”。
