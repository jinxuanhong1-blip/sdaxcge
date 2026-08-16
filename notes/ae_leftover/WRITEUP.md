# ArrayExpress leftover lung ICI / TROP2 / CLDN4

中文在后。

## English

This slice asked for ArrayExpress leftovers: lung ICI / TACSTD2 / CLDN4 **not already covered**, skip E-MTAB-13530 if already hunted, write to `results/w200/AE_leftover/`, and be honest.

### Search

BioStudies ArrayExpress was re-queried with lung × ICI-drug / checkpoint synonyms plus TACSTD2 / TROP2 / CLDN4 / claudin-4. Sixty-eight unique accessions. GEO mirrors (`E-GEOD-*` or any GSE in the record) were excluded. First-wave decisions from PR #3 were tagged `already_covered_*` and not re-claimed.

### Honest result

**No leftover lung ICI cohort with processed expression and ICI response exists on ArrayExpress.**

The first wave already had the only ICI-relevant processed studies (E-MTAB-13704 GEMM aPD-L1 combinations; E-MTAB-15883 anti-PD-1 model). E-MTAB-9451 is checkpoint-motivated NSCLC immune profiling without ICI-treated samples. E-MTAB-10633 and E-MTAB-8867 are ICI-adjacent but metadata-only.

Two first-wave includes (E-MTAB-13704, E-MTAB-10633) did not even re-hit this leftover synonym list because the records say `aPDL1` / TGF-β trap, not nivolumab/pembrolizumab/atezolizumab. That is a wording gap, not a new dataset.

`lung AND TROP2` = 0 hits. `lung AND atezolizumab` = 0 hits. Every TACSTD2/TROP2 ArrayExpress record is intestine or CRC. CLDN4 records are GEO mirrors.

### E-MTAB-13530

Skipped as requested. It is Visium of human NSCLC + adjacent lung, not an ICI experiment. Processed Visium files are present and under 2 GB. The companion scRNA atlas E-MTAB-13526 has annotated objects of 45–58 GB and was not downloaded.

### What we did download

E-MTAB-15784 is the only leftover with open processed matrices under 2 GB that is lung and immune. It is **autologous TIL cytotoxicity against patient-derived NSCLC organoids**. There is no checkpoint-drug arm and no responder/non-responder label. Calling this an ICI leftover would be dishonest.

14 MTX libraries / 126,201 cells / 5 patients (4 LUAD, 1 LUSC). TACSTD2 and CLDN4 are epithelial: ~95–98% and ~88–92% detection in tumor/normal organoids (mean UMI ~14–23 and ~4–8). TIL libraries are ~96–98% CD3D+ and essentially TACSTD2/CLDN4-negative (0–1%). Co-cultures are TIL-heavy; LCP94_CO is a failed/empty epithelial co-culture (TACSTD2 0.2%) and must not be averaged with the others. LCP90 tumor-organoid SDRF points at LCP89_T_ORG files — depositor copy-paste; no unique LCP90 tumor matrix exists.

This supports only a narrow claim: in these NSCLC organoids, TACSTD2 is nearly ubiquitous and higher-UMI than CLDN4, and both are epithelial rather than TIL genes. It does **not** speak to ICI response.

### Metadata-only leftovers worth knowing

- E-MTAB-12508: KP lung + FLT3L/αCD40. The unpublished title says the model is *refractory to checkpoint inhibitors*. The deposited experiment is DC therapy. No matrix.
- E-MTAB-13710: human NSCLC FRC/TLS niches. Companion to excluded vector-IO mouse study E-MTAB-13708. No ICI. No matrix.

Human clinical lung ICI RNA-seq that people actually want (POPLAR/OAK, CheckMate, KEYNOTE, IMpower) is not in open ArrayExpress. Those live in EGA, GEO, or trial-restricted portals.

## 中文

本切片要求：ArrayExpress 上尚未覆盖的肺 ICI / TACSTD2 / CLDN4；若 E-MTAB-13530 已搜过则跳过；结果写到 `results/w200/AE_leftover/`；如实写。

### 结论（先说清楚）

**ArrayExpress 上没有“漏网”的、带处理后表达矩阵和 ICI 疗效标签的肺癌队列。** 真正的 ICI 命中（E-MTAB-13704、E-MTAB-15883，以及无治疗样本的 E-MTAB-9451、仅元数据的 E-MTAB-10633 / 8867）已在第一轮（PR #3）登记。本轮同义词检索甚至重新打不到 13704 / 10633（记录写的是 aPDL1 / TGF-β，不是药名）。

`lung AND TROP2`、`lung AND atezolizumab` 均为 0 条。ArrayExpress 里的 TROP2 全是肠/结直肠；CLDN4 命中是 GEO 镜像。

E-MTAB-13530（NSCLC Visium）按要求跳过，不是 ICI 实验。配套 scRNA E-MTAB-13526 的注释 h5ad 为 45–58 GB，未下载。

唯一下载并分析的遗漏是 **E-MTAB-15784**：NSCLC 类器官 + 自体 TIL 共培养。这是 TIL 疗法平台，**没有 PD-1/PD-L1 臂，没有疗效标签**。类器官里 TACSTD2 几乎全阳性（约 95–98%，UMI 高于 CLDN4）；TIL 库里两者接近 0。不能据此谈 ICI 耐药。

仅元数据的相邻记录：E-MTAB-12508（KP 肺 + FLT3L/αCD40，标题写对 checkpoint 耐药，实验不是 ICI，无矩阵）；E-MTAB-13710（人 NSCLC FRC/TLS，无 ICI，无矩阵）。

临床肺 ICI 转录组（POPLAR/OAK、CheckMate、KEYNOTE、IMpower）不在开放 ArrayExpress，而在 EGA / GEO / 试验门户。
