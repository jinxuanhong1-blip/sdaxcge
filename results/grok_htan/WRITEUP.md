# TACSTD2 / CLDN4 vs immune features in open HTAN/HCA lung objects

**Slice:** HTAN/HCA open processed objects **&lt;2 GB** only.  
**Skipped:** dbGaP **phs002371 raw**; any single file ≥2 GB.  
**No fabricated statistics.** All ρ / p / FDR below were computed by `scripts/grok_htan/`.

---

## English

### What was tested

Whether **TACSTD2** (TROP2; ENSG00000184292) and **CLDN4** (ENSG00000189143) mark an epithelial barrier program that **spatially or compositionally excludes immune cells** in open Human Tumor Atlas Network (HTAN) and Human Cell Atlas (HCA) lung objects.

This tests, in **open processed** data only, a claim in the spirit of Bessede et al. *Clin Cancer Res* 2024 (PMID 38048058) and TROP2/claudin immune-exclusion work: high TACSTD2/CLDN4 ↔ less T-cell infiltration. These HTAN/HCA objects **do not carry ICI response labels**, so the test is immune composition / neighborhood, not ORR/PFS.

### Data used (verified CELLxGENE assets, each &lt;2 GB)

| Key | Source | Object | n | Size |
|-----|--------|--------|---|------|
| `hca_travaglini_10x` | HCA | Travaglini et al. *Nature* 2020 10x lung | 65,662 cells | 0.60 GB |
| `hca_madissoon_lung` | HCA | Madissoon et al. *Genome Biol* 2020 lung parenchyma | 57,019 cells | 0.42 GB |
| `htan_sclc_combined` | HTAN MSK | Chan et al. *Cancer Cell* 2021 combined SCLC/NSCLC | 147,137 cells | 1.46 GB |
| `luad_histology_scrna` | open CXG | LUAD histologic-subtype scRNA (10.1186/s40164-025-00740-6) | 117,266 cells | 0.67 GB |
| `htan_visium_*` (8 slides) | HTAN MSK | Glasner et al. *Nat Immunol* 2023 Visium **mouse** lung | 4,992 spots/slide | 63–78 MB |
| `hca_visium_LngSP10193347` | HCA | Madissoon et al. *Nat Genet* 2023 Visium lung parenchyma | 4,992 spots | 0.75 GB |

Downloads live under `/tmp/grok_htan_data/` (not committed). Manifest: `notes/grok_htan/download_manifest.json`. Full catalog: `notes/grok_htan/catalog.tsv`.

### Explicitly skipped

- **dbGaP phs002371 raw** (HTAN controlled Level 1/2 FASTQ/BAM). Catalogued only.
- **HLCA full** 21.86 GB and **HLCA core** 5.87 GB (Sikkema *Nat Med* 2023).
- **LuCA** core/extended 12.9 / 17.6 GB (Salcher *Cancer Cell* 2022).
- HTAN SCLC epithelial-only 8.95 GB and other Chan subsets ≥2 GB.
- Madissoon “all cells and nuclei” 2.01 GB (just over cap).
- TsankovLab HTAN_Lung Zenodo zip **18.86 GB** (10.5281/zenodo.16546233).
- CELLxGENE “HTAN VUMC” collection `a48f5033-…` is **colorectal polyps**, not lung.

Travaglini + Madissoon **are** HLCA source studies available as &lt;2 GB objects; the integrated HLCA h5ads are not.

### Methods (short)

- Genes resolved from Ensembl IDs or `feature_name` (human + mouse `Tacstd2`/`Cldn4`).
- Signature score = mean per-gene z-score of genes **present**; coverage recorded; missing genes not imputed.
- Signatures: CD8_Tcell, CYT (GZMA/PRF1), Cytotoxic_effector, IFNG_6gene, TLS_12chemokine, TGFB_exclusion, Immune_general, B_cell, Myeloid, NK.
- Associations: Spearman ρ (two-sided). Sample-level contrasts require n≥4. Multiple testing: BH FDR within the concatenated table.
- Spatial: in-tissue spots; 8-NN neighbor mean of immune scores; Stouffer combination across 8 HTAN Visium slides using sign(ρ).
- Median-split MWU on Visium was **unreliable** when TACSTD2 median = 0 (all spots called “high”); those MWU rows are NA and are **not interpreted**.

### Results

#### 1. TACSTD2 and CLDN4 are epithelial and co-expressed

Epithelial Spearman TACSTD2 vs CLDN4:

| Dataset | n epithelial | ρ | p | q |
|---------|-------------:|--:|--:|--:|
| Travaglini 10x (healthy) | 23,626 | **0.683** | ~0 | ~0 |
| Madissoon lung | 5,906 | **0.308** | 1.0×10⁻¹²⁹ | 4.6×10⁻¹²⁹ |
| HTAN SCLC combined | 64,302 | **0.202** | ~0 | ~0 |
| LUAD histology | 37,048 | **0.670** | ~0 | ~0 |

Immune-compartment co-expression is near zero (ρ 0.041–0.068) in all four objects.

Highest mean TACSTD2: basal / AT1 / goblet / club / multiciliated epithelium (Travaglini). LUAD “epithelial cell of lung”: TACSTD2 mean 8.88 (85% &gt;0), CLDN4 mean 11.89 (91% &gt;0). HTAN SCLC epithelial: TACSTD2 sparse (mean 0.14, 8.7% &gt;0) while CLDN4 is common (61% &gt;0) — expected for neuroendocrine-rich SCLC vs LUAD.

#### 2. Sample-level immune fraction: composition vs epithelial-intrinsic expression

**All-cell mean** TACSTD2/CLDN4 vs immune fraction is often **negative** (more epithelium in the library → less immune). Examples:

- LUAD `TACSTD2_all` vs `frac_immune`: ρ=−0.562, p=0.015, q=0.028, **n=18**
- SCLC `CLDN4_all` vs `frac_immune`: ρ=−0.687, p=5.2×10⁻⁷, q=1.2×10⁻⁶, **n=42**
- Madissoon `CLDN4_all` vs `frac_immune`: ρ=−0.614, p=0.005, q=0.010, **n=19**
- Travaglini `TACSTD2_all` vs `frac_cd8`: ρ=−0.821, p=0.023, q=0.042, **n=7** (small n)

This is **not** evidence of exclusion; it is a mixing fraction.

**Epithelial-restricted** mean TACSTD2 vs immune features is **not negative**. Where n is adequate:

- SCLC `TACSTD2_epi` vs `frac_immune`: ρ=**+0.563**, p=1.0×10⁻⁴, q=2.3×10⁻⁴, n=42
- SCLC `TACSTD2_epi` vs TLS_12chemokine: ρ=**+0.641**, p=4.9×10⁻⁶, q=1.1×10⁻⁵, n=42
- SCLC `TACSTD2_epi` vs IFNG_6gene: ρ=**+0.517**, p=4.5×10⁻⁴, q=9.8×10⁻⁴, n=42
- LUAD `TACSTD2_epi` vs sample CD8_Tcell score: ρ=**+0.748**, p=3.6×10⁻⁴, q=7.8×10⁻⁴, n=18
- LUAD `TACSTD2_epi` vs IFNG_6gene: ρ=**+0.829**, p=2.1×10⁻⁵, q=4.9×10⁻⁵, n=18
- Madissoon `TACSTD2_epi` vs `frac_immune`: ρ=**+0.691**, p=0.001, q=0.002, n=19

Travaglini epithelial-vs-fraction tests are non-significant (n=7). SCLC `frac_cd8` is NA because the published coarse label is “T cell”, not CD8.

**Interpretation:** in these open objects, tumors/samples whose **epithelial cells** express more TACSTD2 have **equal or higher**, not lower, immune / IFN / TLS scores. The Bessede-style “high TACSTD2 = T-cell excluded” rule is **not supported** here.

#### 3. Within-epithelial cell-level immune gene scores

These scores measure immune **transcripts inside cells labeled epithelial** (ambient RNA, doublets, or epithelial IFN programs) — not infiltration.

- Healthy Travaglini: TACSTD2 vs Immune_general ρ=−0.605 (n=23,626); vs CD8_Tcell ρ=−0.182. Barrier epithelial cells simply lack lymphocyte genes.
- LUAD / SCLC / Madissoon: weak-to-moderate **positive** correlations (e.g. LUAD TACSTD2 vs IFNG_6gene ρ=0.513; vs TGFB_exclusion ρ=0.581; SCLC TACSTD2 vs IFNG_6gene ρ=0.274). Do not treat as infiltration.

#### 4. Spatial (Visium)

**HTAN Glasner Visium is mouse** (`ENSMUSG*`, `Tacstd2`/`Cldn4`). Ctrl vs DT slides are a Treg-depletion model, not human ICI.

Stouffer meta across **8** HTAN slides (median ρ, two-sided Stouffer p):

| Test | median ρ | Stouffer p | q |
|------|---------:|-----------:|--:|
| Spot CLDN4 vs TGFB_exclusion | +0.173 | 7.5×10⁻²¹⁸ | 2.2×10⁻²¹⁶ |
| Spot TACSTD2 vs IFNG_6gene | +0.091 | 3.4×10⁻⁷³ | 2.5×10⁻⁷² |
| Spot TACSTD2 vs CLDN4 | +0.067 | 3.6×10⁻²⁴ | 9.4×10⁻²⁴ |
| Spot CLDN4 vs CD8_Tcell | −0.025 | 4.3×10⁻⁴ | 5.0×10⁻⁴ |
| Neighbor CLDN4 vs Immune_general | **−0.135** | 2.6×10⁻¹⁰⁶ | 2.5×10⁻¹⁰⁵ |
| Neighbor CLDN4 vs CD8_Tcell | **−0.090** | 3.7×10⁻³⁵ | 1.3×10⁻³⁴ |
| Neighbor TACSTD2 vs Immune_general | +0.067 | 1.1×10⁻⁵⁰ | 6.1×10⁻⁵⁰ |
| Neighbor TACSTD2 vs CD8_Tcell | +0.020 | 1.2×10⁻¹⁶ | 2.4×10⁻¹⁶ |

**Cldn4** spots have slightly fewer immune / CD8 neighbors (direction consistent with local exclusion, small effect). **Tacstd2** neighbors are weakly **enriched**, not depleted, for immune / IFN scores.

Human HCA Visium parenchyma (n=2,081 in-tissue): TACSTD2–CLDN4 ρ=0.248, p=1.5×10⁻³⁰; both genes **positively** correlate with IFNG_6gene (TACSTD2 ρ=0.244; CLDN4 ρ=0.229) and Immune_general. Compatible with airway/gland structure, not tumor exclusion.

### Caveats

1. No ICI outcome labels in these objects.
2. HLCA / LuCA / Tsankov integrated atlases were skipped as &gt;2 GB; results are study-level, not the full integrated lung cancer atlas.
3. phs002371 raw (including unpublished HTAN Level 1/2 lung) was not used.
4. SCLC TACSTD2 is sparse; LUAD is the better human tumor test.
5. Glasner Visium is **mouse**, multiplex-deconvolution scores in `obs` were not used as ground truth.
6. Visium spots are mixed; neighbor correlations are not single-cell distances.
7. Cell-level immune scores in epithelial cells are confounded by ambient RNA.
8. Sample n is modest (7–42). Travaglini n=7 is underpowered.
9. LUAD histology object is open CELLxGENE lung cancer, not an HTAN/HCA atlas; included as the only &lt;2 GB LUAD scRNA with both compartments.

### How to rerun

```bash
pip install -r scripts/grok_htan/requirements.txt
python3 scripts/grok_htan/00_catalog.py
python3 scripts/grok_htan/01_download.py          # writes /tmp/grok_htan_data
PYTHONPATH=scripts/grok_htan python3 scripts/grok_htan/02_analyze_scrna.py
PYTHONPATH=scripts/grok_htan python3 scripts/grok_htan/03_analyze_spatial.py
```

Tables: `results/grok_htan/tables/`. Figures: `results/grok_htan/figures/`.

---

## 中文

### 检测了什么

在 **HTAN / HCA 开放、单文件 &lt;2 GB** 的肺脏处理后对象中，检验 **TACSTD2（TROP2）** 与 **CLDN4** 是否构成上皮屏障程序，并在空间或样本组成上 **排斥免疫细胞**。  
不使用 dbGaP **phs002371 原始数据**。这些对象 **没有 ICI 疗效标签**，因此只检验免疫组成 / 邻域，不检验 ORR/PFS。不编造统计量。

### 使用的数据

- HCA Travaglini 10x 健康肺（65,662 细胞，0.60 GB）
- HCA Madissoon 肺实质（57,019 细胞，0.42 GB）
- HTAN MSK SCLC 合并对象 Chan 2021（147,137 细胞，1.46 GB）
- 开放 CELLxGENE LUAD 组织学亚型 scRNA（117,266 细胞，0.67 GB）
- HTAN MSK Glasner 2023 **小鼠** 肺 Visium 8 张切片（各约 70 MB）
- HCA Madissoon 2023 人肺实质 Visium 1 张（0.75 GB）

跳过：phs002371 raw；HLCA 整合 5.87/21.86 GB；LuCA 12.9/17.6 GB；Tsankov HTAN_Lung zip 18.86 GB；Chan 上皮-only 8.95 GB；CXG「HTAN VUMC」实为结直肠息肉。

### 主要结果

1. **共表达（上皮）**：Travaglini ρ=0.683（n=23,626）；Madissoon ρ=0.308（n=5,906）；SCLC ρ=0.202（n=64,302）；LUAD ρ=0.670（n=37,048）。免疫区室共表达接近 0。TACSTD2/CLDN4 最高在基底/AT1/杯状/club/纤毛上皮。SCLC 上皮 TACSTD2 稀疏（阳性率 8.7%），CLDN4 较常见（61%）；LUAD 上皮两者均高（85% / 91%）。

2. **样本水平（关键）**：全细胞均值 TACSTD2/CLDN4 与免疫比例常为负，这是 **上皮占比** 造成的混合效应（如 SCLC CLDN4_all vs frac_immune ρ=−0.687，n=42，p=5.2×10⁻⁷）。  
   **限定上皮细胞** 后，TACSTD2 与免疫/IFN/TLS **不为负**：SCLC TACSTD2_epi vs frac_immune ρ=+0.563（n=42，p=1.0×10⁻⁴）；vs TLS ρ=+0.641；LUAD TACSTD2_epi vs 样本 CD8 评分 ρ=+0.748（n=18，p=3.6×10⁻⁴）；vs IFNG ρ=+0.829。  
   **不支持**「TACSTD2 高 = T 细胞排斥」的简单规则。

3. **上皮细胞内免疫基因评分**：健康肺为负（Travaglini TACSTD2 vs Immune_general ρ=−0.605），反映细胞身份而非浸润；肿瘤对象为弱–中等正相关，受环境 RNA/双细胞干扰，不能当作浸润。

4. **空间**：Glasner Visium 为 **小鼠**。8 张切片 Stouffer：Cldn4 与 TGFB_exclusion 正（中位 ρ=+0.173）；Cldn4 的 8 邻域 Immune_general 为负（中位 ρ=−0.135），方向符合局部排斥、效应小。Tacstd2 邻域免疫/IFN 为弱正，**不是**排斥。人 HCA Visium 中两基因与 IFN/免疫均为正，符合气道/腺体结构。

### 限制

无 ICI 结局；整合图谱因 &gt;2 GB 未用；phs002371 raw 未用；SCLC 的 TACSTD2 稀疏；小鼠 Visium 不能外推人肿瘤 ICI；Visium spot 为混合细胞；上皮内免疫评分受环境 RNA 影响；样本数 7–42。中位拆分 MWU 在 TACSTD2 中位数为 0 时失效，未解释。

### 复现

见上文 English “How to rerun”。输出仅在 `notes/grok_htan/`、`scripts/grok_htan/`、`results/grok_htan/`。
