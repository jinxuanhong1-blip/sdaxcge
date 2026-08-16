# methods/geomx — GeoMx WTA playbook (public data only)

Bilingual (zh+en) methods playbook for **GSE271689-class** GeoMx Whole Transcriptome Atlas
studies: ROI compartments (PanCK/CK, CD45, CD68), mixed-effects models for multiple ROIs
per patient, overall-survival models, and compartment-specific questions for **TACSTD2**
and **CLDN4** in lung.

**中文:** 面向 **GSE271689 类** GeoMx 全转录组研究的中英双语方法学手册：ROI 分区（PanCK/CK、
CD45、CD68）、每位患者多个 ROI 的混合效应模型、总生存模型，以及肺癌中 **TACSTD2** / **CLDN4**
的分区特异性问题。

This directory is methods-only. It contains no counts. Templates are scaffolds, not a
validated pipeline. **Public GeoMx only** — see [`public_data.md`](public_data.md).

**中文:** 本目录仅含方法学，不含计数。模板是脚手架，不是经过验证的流水线。**仅使用公开 GeoMx
数据** — 见 [`public_data.md`](public_data.md)。

## Read first / 先读

| File | What |
|---|---|
| [`playbook.md`](playbook.md) | Full methods (zh+en): units, QC, normalization, LMM, targets, OS, power, pitfalls |
| [`references.md`](references.md) | Annotated 2024–2026 DSP stats papers + foundational methods (PMID/DOI verified) |
| [`public_data.md`](public_data.md) | GSE271689 / GSE292098 fetch rules; S100B vs PanCK note |
| [`checklists/reporting_checklist.md`](checklists/reporting_checklist.md) | Submission checklist (zh+en) |
| [`templates/methods_text/`](templates/methods_text/) | Fill-in methods paragraphs (en, zh) |
| [`templates/config/study_config.example.yml`](templates/config/study_config.example.yml) | All analysis choices in one file |
| [`templates/R/`](templates/R/) | Runnable scaffolds, `00` → `08` |

## Run order / 运行顺序

From the **repository root**, after copying the example config:

```bash
cp methods/geomx/templates/config/study_config.example.yml \
   methods/geomx/templates/config/study_config.yml
# edit paths, then:
export GEOMX_CONFIG=methods/geomx/templates/config/study_config.yml

Rscript methods/geomx/templates/R/08_fetch_public_geo.R   # public DCC only
# join public phenotype into methods/geomx/data/public/GSE271689/annotation_scaffold.csv
Rscript methods/geomx/templates/R/00_setup.R
Rscript methods/geomx/templates/R/01_qc_dcc_to_spe.R
Rscript methods/geomx/templates/R/02_normalization_batch.R
Rscript methods/geomx/templates/R/03_de_compartment.R
Rscript methods/geomx/templates/R/04_target_genes.R        # TACSTD2 / CLDN4
Rscript methods/geomx/templates/R/05_patient_aggregation.R
Rscript methods/geomx/templates/R/06_survival_os.R         # skip if OS is not public
Rscript methods/geomx/templates/R/07_power_simulation.R    # no data required
```

Scripts 00, 07 and the helpers in `utils_geomx.R` run with base R + `yaml` + `survival`.
Scripts 01–03 need Bioconductor (`GeomxTools`, `standR`, `edgeR`, `limma`, `SpatialExperiment`).

**中文:** 脚本 00、07 与 `utils_geomx.R` 只需 base R + `yaml` + `survival`。01–03 需要
Bioconductor。若公开数据没有 OS 表，做到 04 即止，不要为了 06 导入受限临床信息。
