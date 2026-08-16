# B6_radius — honest rework of CLDN4 vs immune neighborhoods

**Verdict: DOES_NOT_SUPPORT_CLAIM.** Tumor-only Visium hex rings 1/2/3, with immune-rich spots removed from the CLDN4 score, still give **|median partial ρ| ≤ 0.038** (Wilcoxon across 20 tumour sections all p ≥ 0.22). This does **not** rescue the PR67 mismatch (|ρ| ≤ 0.06) and does **not** support the user claim that CLDN4-high tumour domains spatially avoid immune niches.

GeoMx GSE271689 was restricted to the **tumor / CK compartment**. That test is same-AOI co-expression, not spatial rings. Yale (n=37) CLDN4 vs CD8A ρ = −0.39 (p=0.018) did not replicate in the Greek tumor matrix (n=61, CLDN4 vs T-eff ρ = −0.05, p=0.68).

Thresholds were pre-specified. We did not tune filters to enlarge |ρ|.

---

## 中文结论

**结论：不支持用户声称。** 仅肿瘤切片、仅上皮高且非免疫富集的 spot、Visium 六边形 1/2/3 环邻域，CLDN4 与广义免疫邻域的截面中位偏 Spearman **|ρ| ≤ 0.038**，跨 20 张切片的 Wilcoxon 均不显著。PR67 的 |ρ|≤0.06 错配**没有被救回来**。

GeoMx 只分析肿瘤（CK）区室，这是**同一 AOI 内共表达**，不是空间环。Yale 发现集 CLDN4–CD8A ρ=−0.39 在 Greek 验证集未复现。

---

## 1. Why this rework exists / 为何重做

| Item | Value |
|---|---|
| User claim (B6) | CLDN4-high regions spatially avoid immune niches |
| PR67 (`fable_spatial`) | Visium CLDN4 vs immune neighborhood \|median partial ρ\| ≤ 0.06; no strong exclusion |
| This folder | `results/rework/B6_radius/` only |

Requested changes vs PR67:

1. **Tumor-spot only** — E-MTAB-13530 `P*_T*` sections; index spots = epithelial-high.
2. **Radii 1 / 2 / 3 rings** — Visium hex-grid graph distance, not “≤6 nearest pixels”.
3. **Exclude immune-rich spots from the CLDN4 score** — immune score ≥ section 75th percentile.
4. **GeoMx tumor compartment only** — GSE271689 CK / PanCK AOIs; no stromal pairing.

GSE189487 (in PR67) was **not** requested here and was not re-run.

---

## 2. Data (processed only) / 数据

| Dataset | Platform | What was used | Decision |
|---|---|---|---|
| [E-MTAB-13530](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530) | 10x Visium WTA | 20 tumour sections (`P*_T*`); `filtered_feature_bc_matrix.h5` + `spatial.tar` (~358 MB) | **Ran** |
| [GSE271689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271689) | GeoMx DSP WTA | Tumor-compartment sheets from [Nat Genet 2025](https://www.nature.com/articles/s41588-025-02351-7) MOESM10 (Fig. 6b Yale n=37; Fig. 6d Greek n=61); GEO DCC + series matrix as sensitivity | **Ran** |

FASTQ/SRA were not downloaded. SHA-256 and bytes: `download_manifest.json`. Cache: `/tmp/b6_radius_data` (not committed).

GEO series-matrix cell types: CK 212, CD68 192, CD45 175, NA 7. Tumor compartment = **CK only**.

---

## 3. Methods (pre-specified) / 预指定方法

### Visium

Per tumour section:

- In-tissue spots, UMI ≥ 250; CP10K + log1p.
- Epithelial score: EPCAM, KRT7/8/18/19/17, CDH1, ELF3, SFN (excludes TACSTD2/CLDN4/CLDN3/7).
- Broad immune and T-effector scores: same marker lists as PR67 (`analyze.py`).
- **Tumor spots:** epithelial score ≥ section median.
- **Immune-rich:** broad immune ≥ section 75th percentile. These spots are **removed from the CLDN4/TACSTD2 index**. Neighborhood immune still uses every QC spot in the ring.
- **Rings:** Visium 6-neighbor lattice (row, col±2 and row±1, col±1). Ring *k* = graph distance *k*. Self excluded. Require ≥ 3 neighbors.
- Primary stat: **partial Spearman** of index-spot CLDN4 vs ring immune, controlling for the index epithelial score. Cross-section inference: Wilcoxon signed-rank on the 20 ρ values.
- Stouffer combination of per-section p-values is reported but is **anti-conservative** (spots are spatially dependent). Wilcoxon is the primary p.
- Sensitivities (reported, not used to chase a story): keep immune-rich spots in the CLDN4 score; also control for own immune; PR67-like k=6 pixel neighbors; T-effector neighborhood; TACSTD2.

### GeoMx

- **Primary:** published tumor-compartment matrices (Yale Q3-like counts, log2(x+1); Greek batch-corrected → rank inverse-normal). Stromal sheets 6f/6h were **not** used.
- Correlate tumor-AOI CLDN4/TACSTD2 with a T-eff / broad-immune score **in the same AOI**.
- **This is not a ring neighborhood test.** GeoMx AOIs have no hex grid.
- **Yale T-eff score is thin:** only CD8A and CXCL9 are on the Fig. 6b gene list. Single-gene CD8A is the cleaner Yale number.
- DCC raw CPM (212 CK AOIs) is a sensitivity only. Positive ρ there is **not** interpreted as biology (library-size / cellularity).

`did_we_tune_thresholds`: **false**.

---

## 4. Visium results / Visium 结果

20/20 tumour sections ran (8 patients). Index spots after immune-rich exclusion: 461–1442 per section (`visium_section_qc.tsv`).

**Primary: CLDN4 vs broad-immune neighborhood, index = tumor ∩ not immune-rich, partial ρ (epithelial-controlled)**

| Ring | n sections | median ρ | n negative | Wilcoxon p | Stouffer p (anti-conservative) |
|---|---:|---:|---:|---:|---:|
| 1 | 20 | **−0.038** | 13/20 | 0.22 | 0.0038 |
| 2 | 20 | **+0.016** | 8/20 | 0.78 | 0.70 |
| 3 | 20 | **+0.013** | 8/20 | 0.60 | 0.12 |

|median ρ| at every ring is **≤ 0.038**, still inside the PR67 |ρ| ≤ 0.06 band. Ring 1 is weakly negative and not significant by Wilcoxon. Rings 2–3 flip positive.

Figure: `figures/visium_cldn4_rho_by_ring.png`.

**T-effector neighborhood (same index; not primary)** — direction is **positive**, opposite of exclusion:

| Ring | median partial ρ | Wilcoxon p |
|---|---:|---:|
| 1 | +0.025 | 0.070 |
| 2 | +0.032 | 0.036 |
| 3 | +0.058 | 0.027 |

**PR67-like k=6 pixel neighbors** (immune-rich still excluded): CLDN4 vs broad immune median partial ρ = −0.035, Wilcoxon p = 0.29. Including immune-rich spots in the CLDN4 score: −0.030, p = 0.15. Same small-effect regime as PR67.

TACSTD2 vs broad-immune rings is also small (median partial ρ −0.040 / −0.063 / −0.030; Wilcoxon p 0.23 / 0.45 / 0.52).

Full tables: `visium_per_section.tsv`, `visium_meta.tsv`.

---

## 5. GeoMx tumor compartment / GeoMx 肿瘤区室

**Published tumor matrices (primary GeoMx)** — `geomx_tumor_compartment_corr.tsv`

| Cohort | Gene | vs | n | Spearman ρ | p |
|---|---|---|---:|---:|---:|
| Yale tumor | CLDN4 | T-eff (CD8A+CXCL9 only) | 37 | −0.361 | 0.028 |
| Yale tumor | CLDN4 | CD8A | 37 | −0.386 | 0.018 |
| Yale tumor | CLDN4 | broad immune (9 genes present) | 37 | −0.312 | 0.060 |
| Yale tumor | TACSTD2 | T-eff | 37 | −0.095 | 0.58 |
| Greek tumor | CLDN4 | T-eff (8 genes) | 61 | −0.053 | 0.68 |
| Greek tumor | CLDN4 | CD8A | 61 | −0.035 | 0.79 |
| Greek tumor | CLDN4 | broad immune | 61 | +0.034 | 0.80 |
| Greek tumor | TACSTD2 | T-eff | 61 | −0.309 | 0.015 |

Yale CLDN4–CD8A is a discovery-only same-AOI anti-correlation. It is **not** a spatial-ring result and **failed** in the Greek tumor compartment. We do not treat it as confirmation of B6.

**GEO DCC CPM (sensitivity, 212 CK AOIs):** CLDN4 vs T-eff ρ = **+0.52**. Opposite sign to the published Q3/batch-corrected matrices. This is the expected artefact of raw CPM (depth / cellularity). **Do not use DCC CPM as evidence for or against exclusion.**

Figure: `figures/geomx_tumor_cldn4_vs_immune.png`.

---

## 6. Honest reading / 诚实解读

1. The user claim is **spatial exclusion of meaningful size**. After the requested rework, Visium still shows **tiny** |ρ| and no Wilcoxon-significant negative ring effect. **Claim not supported.**
2. PR67’s |ρ| ≤ 0.06 mismatch is **reproduced**, not fixed, by tumor-only + rings + immune-rich exclusion.
3. The weak positive CLDN4–T-eff ring signal (rings 2–3) is the opposite of “avoid immune niches.” Effect size remains small.
4. GeoMx tumor-compartment CLDN4–CD8A in Yale is a **different question** (same AOI, not neighborhood) and does not replicate in Greek.
5. Stouffer p-values look small for some Visium tests because each section has hundreds of spots. They are **not** the primary inference.

---

## 7. Limitations / 局限

- Visium spots are mixtures. Partial correlation and immune-rich exclusion reduce, but do not remove, composition confounding.
- Epithelial-high = top 50% of an 9-gene score, not pathologist tumour masks or CNV calls.
- Immune-rich = 75th percentile of a marker score, not deconvolution.
- No multiple-testing correction across genes × rings × immune axes; all p-values are nominal.
- Yale Fig. 6b is a gene subset (T-eff reduced to CD8A+CXCL9).
- GeoMx cannot test hex rings.
- E-MTAB-13530 has no ICI labels. GSE271689 ICI OS was **not** re-analysed (out of scope for this radius rework).

---

## 8. Reproduction / 复现

```bash
pip install -r results/rework/B6_radius/requirements.txt
python3 results/rework/B6_radius/download.py
python3 results/rework/B6_radius/analyze.py
```

Outputs stay in `results/rework/B6_radius/`. Machine summary: `summary.json` (`honest_verdict.label = DOES_NOT_SUPPORT_CLAIM`).
