# Public GeoMx only / 仅使用公开 GeoMx 数据

This playbook is written against **public** GeoMx WTA accessions. Do not add dbGaP,
restricted, or unpublished counts to this repository.

**中文:** 本手册只针对**公开**的 GeoMx WTA 登录号。不要把 dbGaP、受限访问或未发表的计数加入本仓库。

## Accessions / 登录号

| Accession | What it is / 内容 | Paper / 文献 |
|---|---|---|
| [GSE271689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271689) | Yale NSCLC GeoMx WTA DCC (supplementary `GSE271689_RAW.tar`) | Aung et al., *Nat Genet* 2025;57:2482–2493. PMID 41073787 |
| [GSE292098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292098) | Companion GeoMx DSP accession from the same study | same |

Fetch:

```bash
# from the repo root / 在仓库根目录
Rscript methods/geomx/templates/R/08_fetch_public_geo.R
# or
GEOMX_GEO=GSE292098 Rscript methods/geomx/templates/R/08_fetch_public_geo.R
```

DCCs land in `methods/geomx/data/public/<GSE>/dcc/`. An empty phenotype scaffold is written next to them.
**You** join public phenotype from the GEO SOFT / series matrix / the paper supplement.
Overall survival tables are not always in the series matrix; if they are not public, run
everything through script 04 (compartment biology) and stop. Do not import a restricted
clinical table to finish script 06.

**中文:** DCC 落在 `methods/geomx/data/public/<GSE>/dcc/`。旁边会生成一张空的表型脚手架。**由你**从 GEO SOFT /
series matrix / 论文补充材料拼接公开表型。总生存表不一定在 series matrix 里；若未公开，分析做到
脚本 04（分区生物学）即止。不要为了跑完脚本 06 而导入受限临床表。

## PKC / 探针试剂盒

The Human WTA PKC (`Hs_R_NGS_WTA_v1.0.pkc` or the current public equivalent) is distributed
by Bruker/NanoString with the GeoMx NGS toolchain. It is not redistributed here. Point
`paths.pkc_files` at your local copy.

**中文:** Human WTA PKC 由 Bruker/NanoString 随 GeoMx NGS 工具链分发，本仓库不二次分发。
把 `paths.pkc_files` 指向本地副本。

## Morphology-marker mismatch / 形态学标记不一致

The GSE271689 GEO protocol text lists **S100B** as a tumor-compartment marker (with CD45 and
CD68). The *Nature Genetics* methods for the NSCLC cohorts describe **PanCK**. Reconcile
against the published methods and the actual AOI names in the DCC annotation **before**
modelling. Do not silently recode S100B → Tumor.

**中文:** GSE271689 的 GEO 协议文本把肿瘤分区标记写作 **S100B**（与 CD45、CD68 并列），而
*Nature Genetics* 对 NSCLC 队列的方法学写的是 **PanCK**。建模前对照正式发表的方法学和 DCC
注释中的实际 AOI 名称进行核对，不要悄悄把 S100B 改名为 Tumor。

## TACSTD2 / CLDN4

Both targets are on the public Human WTA panel (`TACSTD2`, `CLDN4`). No custom PKC and no
private probe annotation are required. Lung-relevant priors (epithelium-restricted TROP2;
tight-junction CLDN4) are cited in `references.md` A8–A10.

**中文:** 两个靶点都在公开的 Human WTA panel 上，不需要定制 PKC 或非公开探针注释。与肺相关的先验
（上皮限定的 TROP2；紧密连接蛋白 CLDN4）见 `references.md` A8–A10。

## What this repo will never contain / 本仓库永远不会包含

- Patient-identifiable clinical tables beyond what GEO already released
- FASTQ / BAM (controlled-access sequencing)
- Unpublished in-house GeoMx runs
- Vendor contract data

**中文:** 不会包含 GEO 已发布范围之外的可识别临床表、受控测序 FASTQ/BAM、未发表的内部 GeoMx 实验、
或厂商合同数据。
