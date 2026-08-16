# Demo data · 演示数据

Two demo datasets are referenced by the playbook. Both are real, public
accessions; **no accessions or statistics are invented anywhere in this repo.**
本手册引用两个演示数据集，均为真实公开登录号；**本仓库任何处都不编造登录号或统计量。**

Downloaded data and any derived objects stay under `data/` and are git-ignored
(see `data/.gitignore`); only scripts are committed.
下载的数据与派生对象保留在 `data/` 且被 git 忽略（见 `data/.gitignore`），仅提交脚本。

---

## 1. GeoMx DSP-WTA — `GSE271689` (primary demo · 主演示)

- NanoString GeoMx DSP **Whole Transcriptome Atlas**, NSCLC, PD-1 immunotherapy;
  has per-patient outcomes suitable for the ICI OS template.
  GeoMx DSP **全转录组**，NSCLC，PD-1 免疫治疗；含患者级结局，适配 ICI OS 模板。
- Supplementary archive `GSE271689_RAW.tar` ≈ 36 MB, 586 per-AOI `.dcc.gz`
  files — **well under the 2 GB processed-data budget**, so this is the demo we
  drive end-to-end.
  补充档约 36 MB，586 个 AOI 的 `.dcc.gz`——**远小于 2 GB 预算**，作为端到端演示。

```bash
bash methods/spatial/demo/download_geomx_gse271689.sh
bash methods/spatial/demo/run_demo.sh      # guarded runner (see below)
```

Two files are **required but NOT part of the GEO archive** — you must supply them:
两个文件**必需但不在 GEO 档案中**，需自行提供：

1. **WTA PKC** `Hs_R_NGS_WTA_v1.0.pkc` — the probe→gene map, distributed by
   NanoString. Place at `data/pkc/Hs_R_NGS_WTA_v1.0.pkc` (or set `$PKC`).
   WTA 探针→基因映射文件，由 NanoString 提供。
2. **AOI annotation sheet** `geomx_annotation.csv` — one row per AOI. Reconstruct
   from the GEO series matrix / sample characteristics. Minimum columns:
   每个 AOI 一行，由 GEO series matrix / 样本特征重建。最少列：

   | column | meaning · 含义 |
   | --- | --- |
   | `Sample_ID` | DCC id, matches the DCC filename · 与 DCC 文件名对应 |
   | `slide` | slide/scan id · 切片/扫描号 |
   | `roi` | region of interest id · ROI 号 |
   | `segment` | compartment, e.g. Tumor / Immune / Stroma · 区室 |
   | `patient` | patient/subject id (random effect) · 患者号（随机效应） |
   | `response`, `os_months`, `os_event` | optional clinical cols for template 09 · 模板09可选临床列 |

The runner checks every prerequisite and, if anything is missing, prints exactly
what to supply and exits **without producing any results**.
运行器逐项检查前置条件，若有缺失则明确提示并退出，**不产生任何结果**。

---

## 2. Visium — `E-MTAB-13530` (secondary demo · 次演示)

- 10x Genomics Visium of human NSCLC lesions and non-involved lung tissue
  (36 sections); the source study deconvolved it with cell2location.
  10x Visium，人 NSCLC 病灶与非受累肺组织（36 张切片）；原研究用 cell2location 去卷积。
- **Size caveat:** the full raw dataset (all sections + H&E images) can exceed
  the 2 GB budget. For a quick demo, fetch a **single processed section / .h5ad**.
  **体积提醒：** 完整原始数据（全部切片 + H&E 图）可能超过 2 GB。快速演示只取**单张
  处理后切片 / .h5ad**。

```bash
# lists guidance; pass one relative file path to fetch just that file
bash methods/spatial/demo/download_visium_emtab13530.sh
```

Then run the Visium templates 01→05 (see the playbook for the full command list).
随后运行 Visium 模板 01→05（完整命令见手册）。
