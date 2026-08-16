# methods/spatial · Spatial transcriptomics methods · 空间转录组方法

Playbook + runnable templates for **10x Visium** and **NanoString GeoMx DSP-WTA**
tumour-immune microenvironment analysis, centred on the niche/neighbourhood of
**`TACSTD2` (TROP2) / `CLDN4`** tumour spots vs **T / B / TLS** spots, spatially
variable genes, GeoMx compartment models, and ICI overall survival.
面向 **10x Visium** 与 **NanoString GeoMx DSP-WTA** 的肿瘤免疫微环境分析手册与可运行
模板，聚焦 **`TACSTD2` (TROP2) / `CLDN4`** 肿瘤 spot 相对 **T/B/TLS** spot 的生态位/
邻域、空间可变基因、GeoMx 区室模型与 ICI 总生存。

**Scope · 范围:** whole-transcriptome only. CosMx-1K / Xenium-IO are excluded
unless both `TACSTD2` and `CLDN4` are on the panel. **No invented accessions or
statistics.**
仅全转录组。CosMx-1K / Xenium-IO 除非面板含这两个基因否则排除。**不编造登录号或统计量。**

## Start here · 从这里开始

- **[`playbook.md`](playbook.md)** — the full bilingual playbook (zh + en).
  完整双语手册。
- **[`config/config.yaml`](config/config.yaml)** — all thresholds, gene panels,
  and paths. 所有阈值、基因面板与路径。
- **[`demo/README.md`](demo/README.md)** — datasets & how to run the demo.
  数据集与演示运行方式。

## Layout · 结构

```
methods/spatial/
├── playbook.md                        # bilingual playbook · 双语手册
├── config/config.yaml                 # central config · 中央配置
├── environment/
│   ├── requirements-visium.txt        # Python (templates 01-05)
│   └── install-geomx.R                # R/Bioconductor (templates 06-09)
├── templates/
│   ├── _utils.py                      # config + panel gate · 配置与面板门控
│   ├── 01_visium_qc.py
│   ├── 02_visium_normalization.py
│   ├── 03_visium_deconvolution_cell2location.py
│   ├── 04_visium_niche_neighborhood.py   # TACSTD2/CLDN4 vs T/B/TLS · 核心
│   ├── 05_visium_svg_spatialde2.py
│   ├── 06_geomx_qc_normalization.R
│   ├── 07_geomx_compartment_models.R     # linear mixed models · 混合模型
│   ├── 08_geomx_deconvolution_spatialdecon.R
│   └── 09_ici_os_survival.R              # runs only if labels exist · 有标签才运行
└── demo/
    ├── download_geomx_gse271689.sh       # ~36 MB, verified · 已验证
    ├── download_visium_emtab13530.sh
    ├── run_demo.sh                       # guarded runner · 守卫式运行器
    └── data/.gitignore                   # data never committed · 数据不入库
```
