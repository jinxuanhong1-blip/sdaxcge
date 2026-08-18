# Finding — TISMO LLC Cldn4 (mouse / line), public processed only

Additive public **mouse**. **Cldn4-only.** Tacstd2 TISMO 49/64 and the human CLDN4 thesis are taken as given and are not re-scored.

Question: in TISMO **LLC** (and any other TISMO mouse lung line that has both Cldn4 and an ICI outcome), does mouse-level **Cldn4** track T/NK or IFN/MHC, and is there an ICI-response split?

**No user-private 8 KL.** Values are TISMO’s uniformly processed `log2(TPM+1)` gene-module export plus the public vivo metadata table (Zeng et al., *NAR* 2022, PMID 34534350). No FASTQ.

Primary tables: `tables/line_level.tsv`, `tables/mouse_level.tsv`.

## Line-level table

TISMO `cancerType == Lung carcinoma` lines: **CMT-167, LLC** (9 studies, 19 design groups).

Lung lines that appear in the Cldn4 ICB gene export: **LLC**.

| line | TISMO cancer | ICI outcome in TISMO | Cldn4 in ICB export | n mice (ICB export) | mean Cldn4 log2(TPM+1) | mean Cldn4 TPM |
|---|---|---|---|---:|---:|---:|
| **CMT-167** | Lung carcinoma | absent | no | 0 | NA | NA |
| **LLC** | Lung carcinoma | GSE155972 anti-PD1+anti-CTLA4; R/NR = Setdb1 genotype | yes | 33 | 0.285 | 0.401 |

**CMT-167** is a TISMO lung line (GSE100412, orthotopic, untreated, n=3) but has **no ICB arm and no response label**. The Cldn4 ICB gene-module query for CMT-167 returns HTTP 500 (empty ICI set). **KPB25L** appears in the Cldn4 ICB export (GSE124821) but TISMO labels it **Mammary cancer, NOS**, not lung — it is not added. MLE12 is a TISMO lung-adenocarcinoma *cell line* annotation only; it has no in-vivo ICI rows.

The only ICI-outcome lung design is **GSE155972 LLC** (Griffin et al., *Nature* 2021): subcutaneous flank, anti-PD1 + anti-CTLA4, WT vs Setdb1_KO. Response is **confounded with genotype**: WT ICB = Non-responders (n=6); Setdb1_KO ICB = Responders (n=7). Baseline is untreated (n=10 / genotype). Honest ICI n = **1 study, 1 line, 2 genotype arms, 33 mice**.

## Mouse-level table (GSE155972 LLC, n=33)

TISMO values are `log2(TPM+1)`. T/NK = mean of 9/9 genes (Cd8a, Cd3e, Cd3d, Cd2, Nkg7, Gzmb, Prf1, Ifng, Ncr1). IFN = mean-z of 6/6 leftover genes. MHC-I = mean-z of 6/6 leftover genes. z is across these 33 mice.

| arm | group | n mice | Cldn4 log2p1 | Cldn4 TPM | T/NK mean | IFN mean-z | MHC-I mean-z |
|---|---|---:|---:|---:|---:|---:|---:|
| Setdb1_KO | Baseline | 10 | 0.267 | 0.224 | 3.196 | 0.027 | 0.155 |
| Setdb1_KO | ICB | 7 | 0.601 | 1.299 | 4.448 | 0.330 | 0.750 |
| WT | Baseline | 10 | 0.128 | 0.097 | 2.422 | -0.651 | -0.870 |
| WT | ICB | 6 | 0.209 | 0.157 | 3.968 | 0.656 | 0.317 |

Cldn4 sits near the detection floor: mean TPM = **0.40**; **32/33** mice have TPM < 1. One mouse (**SRX8918393**, Setdb1_KO ICB) is Cldn4 = 3.23 log2p1 (TPM ≈ 8.4) and is marked in `mouse_level.tsv`.

### Cldn4 vs T/NK or IFN/MHC (Spearman, mouse unit)

| subset | n | Cldn4 vs T/NK genes ρ (p) | Cldn4 vs TISMO T/NK infil ρ (p) | Cldn4 vs IFN ρ (p) | Cldn4 vs MHC-I ρ (p) |
|---|---:|---|---|---|---|
| all GSE155972 mice | 33 | 0.122 (0.50) | 0.023 (0.90) | 0.120 (0.51) | 0.141 (0.43) |
| drop SRX8918393 | 32 | 0.063 (0.73) | -0.055 (0.76) | 0.201 (0.27) | 0.089 (0.63) |

Per-arm correlations are in `tables/spearman.tsv`. A positive Cldn4–T/NK ρ on near-floor Cldn4 is **not** a Cldn4-high / T-low exclusion pattern.

### ICI (honest)

| contrast | n | Cldn4 Δ log2p1 | Cldn4 Welch p | T/NK Δ | T/NK Welch p |
|---|---|---:|---:|---:|---:|
| WT ICB vs baseline (NR arm) | 10 / 6 | 0.081 | 0.17 | 1.546 | 0.0012 |
| Setdb1_KO ICB vs baseline (R arm) | 10 / 7 | 0.334 | 0.48 | 1.252 | 0.01 |
| treated R vs NR (genotype-confounded) | 6 / 7 | 0.392 | 0.41 | 0.480 | 0.29 |

There is **no within-genotype responder vs non-responder** contrast in TISMO LLC. Do not read the R vs NR row as an ICI-outcome test of Cldn4.

T/NK / IFN genes **do** move with ICB (same design can detect an immune shift). Cldn4 does not at a usable effect size; it stays near floor.

## Verdict

TISMO LLC Cldn4 is **present** (not absent) but **near the expression floor**. The mouse/line-level table exists.

- **Lines with Cldn4 + ICI outcome:** LLC only (GSE155972). CMT-167 has Cldn4-capable metadata as lung but **no ICI outcome**.
- **Mouse n:** 33 (10 WT baseline, 6 WT ICB/NR, 10 Setdb1_KO baseline, 7 Setdb1_KO ICB/R).
- **Cldn4 vs T/NK or IFN/MHC:** see Spearman table; Cldn4 is too low to support a Cldn4-high / immune-low or Cldn4-high / IFN-high LLC state.
- **ICI response:** labels exist but are the Setdb1 genotype. Cldn4 ICB vs baseline is not significant on either arm.

This does **not** reopen the Tacstd2 49/64 tally (2/64 of those cohorts are this LLC study). It does not use private 8 KL.

## Methods (short)

- Inclusion: TISMO `cancerType == Lung carcinoma` AND (Cldn4 in the public ICB gene export) AND (ICI treatment or response field). That intersection is LLC GSE155972.
- Expression: TISMO gene-module `downVivoExprn` type=3, `icbList` = the six ICB treatments, `tumorList` = LLC (and All for the lung-ICI census). Values = `log2(TPM+1)` as deposited.
- T/NK score = unweighted mean of present T/NK genes (already log2p1). Extra T/NK = mean-z of TISMO public infiltrate rows (CD8 T mMCPcounter, T CD8 TIMER/CIBERSORT, T mMCPcounter, T NK xCell, NK mMCPcounter). IFN / MHC-I = leftover 6-gene lists from the GSE239485 Cldn4 page, mean of gene-wise z on these mice.
- Tests: Spearman on the mouse; Welch t and two-sided MWU when both arms have n≥2. Outlier SRX8918393 is kept in the primary table and dropped only in the sensitivity rows.
- Public processed only. No GEO FASTQ. No user-private 8 KL.

Reproduce:

```bash
python3 methods/tismo_llc_cldn4/fetch.py
python3 methods/tismo_llc_cldn4/analyze.py
```

## Files

| File | Role |
|---|---|
| `tables/line_level.tsv` | TISMO lung lines: Cldn4 + ICI present/absent |
| `tables/mouse_level.tsv` | one row per GSE155972 mouse |
| `tables/spearman.tsv` | Cldn4 vs T/NK, IFN, MHC-I |
| `tables/contrasts.tsv` | ICB vs baseline and confounded R vs NR |
| `tables/lung_inventory.tsv` | all TISMO lung-carcinoma design groups |
| `tables/gene_coverage.tsv` | which genes were in the LLC export |
| `tables/summary.json` | machine-readable verdict |
| `raw/` | live TISMO exports used here |

## 中文摘要

公开 TISMO 小鼠、只看 **Cldn4**。不做 Tacstd2 49/64，不用私有 8 KL。

TISMO 肺癌系只有 **LLC** 和 **CMT-167**。有 Cldn4 **且** 有 ICI 结局的只有 **LLC GSE155972**（皮下，抗 PD-1+CTLA4；WT ICB=NR n=6，Setdb1_KO ICB=R n=7，对照各 n=10）。CMT-167 无 ICB。Cldn4 接近检测下限（多数 TPM<1）。小鼠水平 Cldn4 与 T/NK、IFN、MHC-I 的相关见 `tables/spearman.tsv`；不能支持 Cldn4 高 / 免疫低。R vs NR 与基因型完全重叠，不是独立 ICI 结局检验。
