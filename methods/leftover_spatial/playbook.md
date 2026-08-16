# Leftover public lung spatial series (2023–2026)

Bilingual methods for an **additive** leftover slice. User B6 (same-spot / compartment anti-colocalization of CLDN4/TACSTD2 with T/B in tumor Visium and GeoMx) is taken as given and is not re-run.

**中文：** 只补 2023–2026 年尚未做过的公开肺 Visium / GeoMx / CosMx。不重做 GSE221322、10x FFPE Visium（PR #126）、E-MTAB-13530（PR #134）、CosMx（PR #164）。系列必须在真实矩阵里测到 TACSTD2 或 CLDN4；两个基因都没有则跳过。

---

## Scope

| Include | Exclude |
|---|---|
| Human lung GeoMx / Visium / CosMx, public 2023–2026 | GSE221322, PR #126 Visium, E-MTAB-13530, PR #164 CosMx |
| Processed count / Q3 matrices with TACSTD2 **or** CLDN4 | FASTQ-only, dbGaP, files used only as 10 GB+ archives |
| Extra figures: target vs immune **neighborhood** (hex ring, FOV, or antibody compartment) | Claim-audit language; re-analysis of B6 folders |

Neighborhood is defined by the platform:

- **Visium:** honeycomb rings 1/2/3 from `(row, col)` using neighbors `(row, col±2)` and `(row±1, col±1)`.
- **GeoMx:** same-AOI immune score; antibody-cut compartments (tumor/immune, CD45±, PanCK/stroma) when the deposit labels them.
- **CosMx:** FOV means if public files lack CenterX/Y.

---

## 范围（中文）

只做公开、已处理、含 TACSTD2 或 CLDN4 的肺空间系列。Visium 用六边形环；GeoMx 用同一 AOI / 分割室；CosMx 若无细胞坐标则退回 FOV。不重做已指定系列，不写 claim-audit 措辞。

---

## Signatures (locked before looking at leftover ρ)

| Score | Genes |
|---|---|
| Target epithelial | `CLDN4`, `TACSTD2` |
| Broad epithelial (control) | `EPCAM`, `KRT8`, `KRT18`, `KRT19`, `CDH1`, `KRT7` |
| T | `CD3D`, `CD3E`, `CD3G`, `CD2`, `CD8A`, `CD8B`, `TRAC`, `CD247`, `IL7R` |
| B | `MS4A1`, `CD79A`, `CD79B`, `CD19`, `BANK1`, `CD22` |

A gene is used only if it is present in that matrix. Scores are the mean of available members.

---

## Visium

1. Keep in-tissue spots with ≥200 detected genes.
2. Log1p(CP10K).
3. **Same-spot:** Spearman of CLDN4, TACSTD2, or the two-gene mean versus T, B, and T+B.
4. **Neighborhood:** Spearman of spot CLDN4 (or TACSTD2) versus mean T+B of hex-ring *k* neighbors, *k* = 1, 2, 3; require ≥3 neighbors.
5. **Partial:** rank residual Spearman controlling for the broad epithelial score.
6. Section is the unit of inference. Report median ρ, how many sections are negative, and a two-sided Wilcoxon signed-rank *p* on the per-section ρ values.

---

## GeoMx

1. Use the depositor Q3 / log2-Q3 gene × AOI matrix (no DCC re-quantification unless that is the only public file).
2. Same signatures as Visium.
3. Spearman across AOIs, then within labeled compartments.
4. When paired tumor/immune or PanCK/stroma AOIs exist, add a paired Wilcoxon or Mann–Whitney on target levels (compartment geography, not a hex grid).

---

## CosMx

1. Stream only the signature gene rows from the public count matrix.
2. QC: require a matching annotation row and `nCount ≥ 20` when that column exists.
3. Same-cell Spearman, plus FOV-mean Spearman if x/y are not in the public coordinate file.
4. Do not invent micron radii when CenterX/Y are absent.

---

## Reporting

Every ρ is written with **n** and **p**. No pooled claim across platforms. Catalog every hunted series, including skips after a real gene check.

Reproduce:

```bash
python3 scripts/leftover_spatial/write_catalog.py
python3 scripts/leftover_spatial/analyze_visium.py
python3 scripts/leftover_spatial/analyze_geomx.py
python3 scripts/leftover_spatial/analyze_cosmx.py
python3 scripts/leftover_spatial/make_figures.py
```

Processed public files are downloaded to `data/` (gitignored). Outputs live under `results/leftover_spatial/`.
