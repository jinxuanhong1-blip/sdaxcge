# methods/seurat_scanpy

Methods-only playbook and parameterized templates for a 2025–2026 single-cell
RNA-seq (scRNA-seq) workflow targeting **lung cancer immune-checkpoint-inhibitor
(ICI)** cohorts, spanning **Seurat v5**, **Scanpy**, and **OSCA / Bioconductor**.

仅方法学的手册与可参数化模板，面向 2025–2026 年**肺癌免疫检查点抑制剂 (ICI)**
队列的单细胞 RNA 测序流程，覆盖 **Seurat v5**、**Scanpy** 与 **OSCA / Bioconductor**。

## No data policy / 无数据原则

This directory contains **methods and code templates only**. It ships **no datasets,
no results, and no invented numbers, gene signatures, thresholds, or citations of
specific findings**. Every cohort-specific value is a `# TODO` you must supply.

本目录**仅含方法与代码模板**，**不含任何数据集、结果，也不编造任何数字、基因签名、
阈值或对具体研究结论的引用**。所有队列特定的值均为需你填写的 `# TODO`。

## Contents / 目录

| File | What it covers |
|---|---|
| [`playbook.md`](playbook.md) | Bilingual (中文 + English) decisions + mechanics: layers/assays, QC, BPCells, sketching, integration APIs, annotation (CellTypist/Azimuth/SingleR), TACSTD2/CLDN4 module scoring, reproducibility. |
| [`templates/seurat_v5_template.R`](templates/seurat_v5_template.R) | Seurat v5 end-to-end: BPCells on-disk counts, split layers, `SketchData`/`ProjectData`, `IntegrateLayers` (CCA/RPCA/Harmony/FastMNN/scVI), Azimuth `lungref`, `AddModuleScore`. |
| [`templates/osca_bioc_template.R`](templates/osca_bioc_template.R) | OSCA/Bioconductor: `SingleCellExperiment`, `emptyDrops`, MAD QC, `scDblFinder`, `scran` norm/HVG/`denoisePCA`, `batchelor::fastMNN`, `SingleR`, `AUCell`. |
| [`templates/scanpy_template.py`](templates/scanpy_template.py) | Scanpy: AnnData layers, per-sample QC + `scrublet`, batch-aware HVGs, PCA/neighbors/UMAP/Leiden (igraph flavor), Harmony, `score_genes`. |
| [`templates/celltypist_scvi_template.py`](templates/celltypist_scvi_template.py) | scVI/scANVI integration + label transfer, CellTypist annotation (`Immune_All_Low`, `Human_Lung_Atlas`). |
| [`signatures/tacstd2_cldn4_modules.md`](signatures/tacstd2_cldn4_modules.md) | How to build & validate the TACSTD2 / CLDN4 modules, with a provenance log (empty by design). |
| [`signatures/modules.template.json`](signatures/modules.template.json) | Machine-readable module template you fill with sourced genes. |

## Suggested reading order / 建议阅读顺序

1. `playbook.md` — understand the decisions and the end-to-end order of operations.
2. The template for your stack (R: Seurat v5 or OSCA; Python: Scanpy + CellTypist/scVI).
3. `signatures/` — curate the TACSTD2 / CLDN4 gene sets before scoring.

## Verified API references (2025–2026) / 已核实 API 参考

- Seurat v5 integration (`IntegrateLayers`, 5 methods, `JoinLayers`), `SketchData`/`ProjectData`, BPCells — Satija Lab Seurat v5 vignettes.
- Azimuth lung reference: `lungref` v2.0.0 (Human Lung Cell Atlas, HLCA v2) via SeuratData.
- CellTypist ≥ 1.7 (`celltypist.annotate`, `models.download_models`); models incl. `Immune_All_Low.pkl`, `Immune_All_High.pkl`, `Human_Lung_Atlas.pkl`.
- Scanpy ≥ 1.10 (Leiden `flavor="igraph"`), scvi-tools (SCVI/SCANVI), OSCA (`scran`/`scater`/`batchelor`/`SingleR`/`AUCell`).
