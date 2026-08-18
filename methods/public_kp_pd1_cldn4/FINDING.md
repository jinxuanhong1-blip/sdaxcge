# Finding — public mouse KP / Kras PD-1 Cldn4 (vs KL-resistant leftovers)

Additive public **MOUSE**. **Cldn4-only.** Thesis already correct. No dual-high. **No private 8-KL.**

Question: among public mouse **lung tumors** with anti-PD-1 / ICB **and** an ICI-sensitive genotype (KP, Kras, KP GEMM, MC38-lung if any), does mouse-level **Cldn4** track T/NK, and does PD-1 vs control move T/NK and Cldn4?

This page is the sensitive-model contrast to KL-resistant leftovers (GSE182228 / GSE232730 / TISMO LLC). Series are scored **separately**. No mega-merge.

Primary table: `tables/GSE114601_pd1_primary.tsv`.

---

## Hunt verdict

Example accessions from the request were checked and are **not** mouse KP lung PD-1:

| Accession | What it actually is |
|---|---|
| GSE131928 | Human GBM |
| GSE139914 | Human memory encoding |
| GSE146256 | *C. elegans* |
| GSE161084 | Mouse psoriasis epidermis |
| GSE179994 | Human anti-PD-1 T cells |
| GSE203353 | Human placenta / GDM |
| GSE207422 | Human neoadjuvant NSCLC (no mouse) |
| GSE155972 | LLC Setdb1 ICB — leftover; Cldn4 at floor. Not re-scored as the only analysis |

**TISMO KP / MC38 / CMT-167:** TISMO has **no KP lung line**. **CMT-167** is lung but has **no ICB arm** (ICB gene export HTTP 500). **MC38** ICB Cldn4 is at the floor (TPM mean 0.021; 22/22 TPM < 1) and is **colorectal**, not MC38-lung. **KPB25L** has high Cldn4 but is mammary. LLC leftover stands.

Usable public processed matrices with a sensitive genotype **and** PD-1 / ICB:

| Series | Model | Assay | PD-1 contrast | Honest n | Cldn4 | Scored |
|---|---|---|---|---|---|---|
| **GSE114601** | KP GEMM (KrasTrp53) nodules | bulk | Vehicle vs anti-PD-1 | **2 vs 2 mice** | **not floor** (Vehicle mean 70.9) | **yes — primary** |
| **GSE157880** | HKP1 (Kras/p53) orthotopic lung | bulk | IgG vs PD-1 at 0 / 4 / 8 Gy | 3 vs 2 / 2 vs 3 / 3 vs 3 **pooled libraries** | present | yes |
| **GSE133604** | KP GEMM lung 10x ± Asf1a | scRNA | Ctrl vs Ctrl+PD-1 | **1 vs 1 library** | sparse in epi (~1.6% pos) | yes; direction only |
| **GSE169194** | KPM (Kras/p53/Msh2) total viable | bulk | IgG vs A2V+aPD-1 (**no PD-1 mono**) | **3 vs 3 mice** | present | yes; combo only |
| GSE246922 | KP CD45− parental vs ICB-relapse | bulk | not on-treatment PD-1 vs control | 3 vs 3 (1st relapse) | present | Cldn4 only; T/NK no-go |
| GSE275877 | LKR13-H Msh2-KO ± ICI | scRNA | 1 library / arm | n=1 | not opened | unfiltered 535 MB MTX / ~34M barcodes; skip |

Catalog: `catalog.tsv` (example misses + leftovers + scored rows).

---

## 1. GSE114601 — KP GEMM anti-PD-1 (primary; n=2 vs 2)

Adeegbe et al. GEO [GSE114601](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114601): lung tumor **nodules**, `gemm model: KrasTrp53`, Vehicle / JQ1 / anti-PD-1 / combo, **n=2 mice / arm**. Processed `GSE114601_counts.normalized.csv.gz`. Cldn4 + 16/16 T/NK genes present.

| arm | n mice | mean Cldn4 (norm) | mean Cldn4 log2p1 | mean T/NK log2p1 |
|---|---:|---:|---:|---:|
| Vehicle | 2 | 70.90 | 6.167 | 0.423 |
| anti-PD-1 | 2 | 31.56 | 4.845 | 0.605 |
| JQ1 | 2 | 89.45 | 6.078 | 0.409 |
| anti-PD-1+JQ1 | 2 | 163.10 | 7.299 | 0.383 |

Cldn4 is **not** at the LLC / TISMO floor.

### Cldn4 vs T/NK (mouse unit)

| subset | n | Spearman ρ | p |
|---|---:|---:|---:|
| all 8 nodules | 8 | −0.40 | 0.32 |
| Vehicle + anti-PD-1 only | 4 | −1.00 | *n=4 exact rank inversion; scipy p=0 — not a powered test* |

Direction is **higher Cldn4, lower T/NK**. n is too small to call a law.

### PD-1 vs Vehicle

| metric | n | Δ (PD-1 − Vehicle) | Welch p | MWU p |
|---|---|---:|---:|---:|
| Cldn4 (norm) | 2 / 2 | **−39.3** | 0.23 | 0.33 |
| Cldn4 log2p1 | 2 / 2 | −1.32 | 0.32 | 0.33 |
| T/NK log2p1 | 2 / 2 | **+0.18** | 0.38 | 0.33 |

Same sign as an ICI-sensitive PD-1 benefit (Cldn4 down, T/NK up). **Honest n=2 vs 2; not significant.** JQ1 alone does not move either score. Combo Cldn4 goes **up** (opposite of PD-1 mono) — do not pool with PD-1 mono.

---

## 2. GSE157880 — HKP1 lung, PD-1 ± RT (pooled libraries)

Ban et al. child of GSE157883. HKP1 = Kras/p53 orthotopic lung. Each GEO sample is **3–4 mice pooled**. Unit = library, not mouse.

Primary (no RT): **IgG 0 Gy n=3 vs PD-1 0 Gy n=2**.

| dose | n IgG / PD-1 | Cldn4 Δ | Cldn4 Welch p | T/NK Δ | T/NK Welch p |
|---|---|---:|---:|---:|---:|
| 0 Gy (primary) | 3 / 2 | −1.52 | 0.39 | −0.40 | 0.34 |
| 4 Gy | 2 / 3 | +4.76 | 0.18 | +0.54 | 0.073 |
| 8 Gy | 3 / 3 | +1.43 | 0.27 | −0.02 | 0.95 |

0 Gy is **not** a T/NK-up PD-1 effect (both Cldn4 and T/NK slightly down; ns). 4 Gy PD-1 libraries are Cldn4-higher and T/NK-higher than IgG+RT — RT-confounded. Spearman Cldn4 vs T/NK at 0 Gy: ρ=−0.70, p=0.19, n=5.

---

## 3. GSE133604 — KP GEMM 10x ± anti-PD-1 (n=1 vs 1)

Li / Wong Asf1a paper. Shared `GSE133604_genes.tsv.gz` + four MTX. Cldn4 present. Whole-tumor 10x (not CD45−).

| library | n cells (UMI≥200) | n epi | n T/NK | T/NK fraction | epi Cldn4 mean log1p | epi Cldn4 %pos |
|---|---:|---:|---:|---:|---:|---:|
| Ctrl | 6074 | 320 | 2863 | 0.471 | 0.0143 | 1.56% |
| Ctrl+PD-1 | 6683 | 341 | 3486 | 0.522 | 0.0041 | 0.59% |
| Asf1a-KO | 6362 | 398 | 1847 | 0.290 | 0.0052 | 0.50% |
| KO+PD-1 | 6094 | 351 | 2785 | 0.457 | 0.0105 | 0.85% |

WT KP: PD-1 library is Cldn4-lower and T/NK-higher than Ctrl. **Direction only (n=1).** Epithelial Cldn4 is sparse — this is not the GSE114601 nodule abundance. Asf1a-KO is a second genotype; do not average with WT.

---

## 4. GSE169194 — KPM total viable, A2V ± aPD-1 (combo)

Martinez-Usatorre / De Palma. `KrasLSL-G12D/+; p53fl/fl; Msh2-/-` — **more ICI-sensitive than KP**. Total viable cells (not TAM/CD4 sorts): IgG n=3, A2V n=3, A2V+aPD-1 n=3 (`total_A2V_aPD1_1` not deposited). **No PD-1 monotherapy.**

| contrast | n | Cldn4 Δ (log2-norm) | Cldn4 Welch p | T/NK Δ | T/NK Welch p |
|---|---|---:|---:|---:|---:|
| A2V+aPD-1 vs IgG | 3 / 3 | −0.71 | 0.44 | +0.15 | 0.84 |
| A2V vs IgG | 3 / 3 | +0.10 | 0.90 | −0.69 | 0.40 |

IgG-only Cldn4 vs T/NK: ρ=−0.50, p=0.67, n=3. All nine total-cell samples: ρ=−0.82, p=0.007 — **treatment-mixed**, not a within-control law. Do not quote p=0.007 as the PD-1 result.

---

## 5. GSE246922 — KP CD45− ICB relapse (T/NK no-go)

Parental / chronic IFNG / ICB-relapsed **CD45−** KP cells, n=3 / group. Cldn4 present (VST ~9.2). T/NK genes are leak, not infiltrate. Parental vs 1st relapse Cldn4 Δ = +0.33 VST, Welch p=0.40. Not an on-treatment PD-1 vs control table.

---

## Combinatorial (not a merge)

`tables/combinatorial_pd1_vs_control.tsv` stacks **per-series** PD-1 (or nearest ICB) contrasts. Expression matrices were **not** concatenated.

| series | clean PD-1 mono? | n | Cldn4 | T/NK |
|---|---|---|---|---|
| GSE114601 | yes | 2 vs 2 mice | down (ns) | up (ns) |
| GSE157880 0 Gy | yes (pooled) | 3 vs 2 libraries | down (ns) | down (ns) |
| GSE133604 WT | yes | 1 vs 1 | down (direction) | up (direction) |
| GSE169194 | **no** (A2V+PD-1) | 3 vs 3 | down (ns) | up (ns) |

Two of three PD-1-mono KP/Kras series move Cldn4 down and T/NK up (GSE114601, GSE133604). HKP1 0 Gy does not. **No series is powered.** This is the public sensitive-model contrast: Cldn4 is **measurable** in KP nodules (GSE114601), unlike TISMO LLC / MC38 floor. It does **not** prove a Cldn4-high / T-excluded KP law.

---

## What this is not

- Not LLC as the only analysis (GSE155972 / GSE239485 / GSE297630 / GSE297632 cataloged, not re-scored).
- Not private 8 KL.
- Not TISMO Tacstd2 49/64.
- Not a mega-merged KP+KPM+HKP1 matrix.
- Not MC38-lung (none found; TISMO MC38 is colon and floor).
- GSE275877 LKR13 ICI scRNA exists as 1-library arms but the public MTX is an unfiltered droplet dump — not scored.

---

## Methods (short)

- Inclusion: public mouse **lung tumor** + processed matrix + Cldn4 + immune or T/NK genes + anti-PD-1/ICB **and** KP / Kras / KPM / HKP1 / MC38-lung. FASTQ not used.
- Cldn4 only. T/NK = unweighted mean of present Cd3d/e/g, Cd2, Cd8a/b1, Cd4, Nkg7, Gzma/b, Prf1, Klrb1c, Ncr1, Klrd1, Klrc1, Ifng (log2p1 on count-like matrices; deposited log2 on GSE169194).
- Unit = biological mouse when deposited (GSE114601, GSE169194 totals). GSE157880 = pooled library. GSE133604 = 10x library. GSE246922 = CD45− tumor.
- scRNA (GSE133604): UMI≥200; epithelial = Epcam+ or (Cdh1+ and Krt8+); T/NK gate = Cd3d/e, Cd8a, Nkg7, Ncr1; Epcam wins over ambient T.
- Tests: Spearman; Welch t and two-sided MWU only if both arms n≥2. n=1 reports Δ only.
- Combinatorial = stacked per-series rows. No joint normalization.

Reproduce:

```bash
# matrices in /tmp/geo_dl (GEO FTP URLs in catalog / analyze.py)
python3 methods/public_kp_pd1_cldn4/analyze.py
```

---

## Files

| File | Role |
|---|---|
| `catalog.tsv` | Hunt: example misses, leftovers, scored series |
| `tables/GSE114601_pd1_primary.tsv` | Sensitive-model PD-1 table (n=2 vs 2) |
| `tables/GSE114601_mouse_level.tsv` | 8 KP nodules |
| `tables/GSE157880_sample_level.tsv` | 16 HKP1 pooled libraries |
| `tables/GSE169194_total_mouse_level.tsv` | 9 KPM total-viable mice |
| `tables/GSE133604_library_level.tsv` | 4 KP 10x libraries |
| `tables/GSE246922_KP_CD45neg.tsv` | KP CD45−; T/NK marked leak |
| `tables/spearman.tsv` | Cldn4 vs T/NK per series / subset |
| `tables/contrasts.tsv` | All arm contrasts |
| `tables/combinatorial_pd1_vs_control.tsv` | Stacked PD-1 rows, not a merge |
| `tables/summary.json` | Machine-readable n |

---

## 中文摘要

公开小鼠、只看 **Cldn4**。不用私有 8 KL。不把 LLC 当唯一分析。

点名的 GSE131928 / 139914 / 146256 / 161084 / 179994 / 203353 / 207422 **都不是**小鼠 KP 肺 PD-1。GSE155972 是 LLC Setdb1，Cldn4 贴地，不重做。TISMO 无 KP 肺系；CMT-167 无 ICB；MC38 Cldn4 贴地且是肠癌，不是 MC38-lung。

能打分的敏感模型：**GSE114601 KP GEMM 结节 anti-PD-1 vs Vehicle，n=2 vs 2**（Cldn4 不在地板；PD-1 后 Cldn4 降、T/NK 升，Welch p=0.23 / 0.38）。**GSE157880 HKP1** 0 Gy 为混合文库 n=3 vs 2，方向不稳。**GSE133604 KP 10x** n=1 vs 1，上皮 Cldn4 很稀，只报方向（Cldn4 降、T/NK 升）。**GSE169194 KPM** 只有 A2V+PD-1 联合，无单药 PD-1。各系列分开报，不合并矩阵。
