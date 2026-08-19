# RESULTS — TISMO / syngeneic cell lines: TROP2 and the tight-junction program

Additive public-data analysis. Thesis treated as given: TROP2-high / immune-resistant tumors run a CLDN4-family tight-junction (TJ) barrier program. **LLC is labeled LLC, not KL.** No private 8-KL matrices. GSE76628 was not used.

## 1. Sources

- **TISMO** (Zeng et al., *Nucleic Acids Res* 2022, PMID 34534350). Official Data Download from https://tismo.pku-genomics.org/: in-vitro expression RDS (605 samples × 21,729 genes), in-vivo expression RDS (1,518 samples), plus cell-line / vitro / vivo annotation tables.
- **GSE274352** (Fernández-García / Barbacid, PMID 39186651). KP and KL murine LUAD **cell lines**, empty-vector vs IFN-β / STING. Primary genotype contrast = empty-vector KL vs KP.
- **GSE167381** (Pierce et al., *Nature* 2021, PMID 33981036). Immortalized LUAD lines LU1/LU2 (LKB1-unrestorable) and LR1/LR2 (LKB1-restorable); vehicle = LKB1-off, 4-OHT restores LKB1 in LR lines.
- **GSE295685** (Tango / STK11-null KL cells ± TNG260). Deposited matrix is **KL only** (no KP arm).

## 2. TISMO inventory — true KL lines are absent

- In-vitro annotation: **605** samples from **49** cell lines (TISMO's 49-line vitro collection).
- In-vivo annotation: **1518** samples from **68** models.
- Lung vitro lines: **LLC, MLE12**.
- Lung vivo models: **CMT-167, LLC**.
- Cell-genotype strings containing Lkb1/Stk11: **0** (two `Alkbh5_KO` B16 rows match a naive `lkb` substring and were excluded). Strings containing Kras+Trp53 as a KL/KP GEMM label: **0**.
- **No TISMO line is a true KL (Kras/Lkb1/Stk11) or KP (Kras/Trp53) GEMM cell line.** KPC/KPB25L are pancreatic / mammary, not lung KP. **LLC is not KL** (spontaneous C57BL/6 lung carcinoma, 1951; TISMO genotype = WT). CMT-167 is a spontaneous CMT-64 derivative. MLE12 is an FVB/N SV40-immortalized alveolar line (1992), not KL/KP.
- Therefore the TISMO tests below use: (i) all 49 vitro lines ranked by Tacstd2, (ii) TISMO ICB R vs NR labels, (iii) public GEO KL/KP cell-line RNA-seq.

## 3. TISMO vitro — TROP2-high lines carry the TJ program

Unit = **cell line**. Expression = median of TISMO-processed baseline samples (`Baseline==1`). **n = 49** lines (2 lung). Tacstd2 is exactly 0 in a large fraction of lines, so the primary split is **TROP2-positive (Tacstd2 > 0; n = 21) vs TROP2-floor (Tacstd2 = 0; n = 28)**. Tertiles and a median split among expressed lines are in the tables.

TJ-core genes (present / requested): `Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, Ocln, Tjp1, Nectin2`.

| contrast | n | statistic | value | p |
|---|---:|---|---:|---:|
| Tacstd2 vs TJ-core (49 lines) | 49 | Spearman ρ | 0.634 | 9.84e-07 |
| Tacstd2 vs Cldn4 | 49 | Spearman ρ | 0.700 | 2.14e-08 |
| TROP2-positive vs floor, TJ-core | 21 vs 28 | MWU Δ median | 0.718 | 0.0003 |
| TROP2-positive vs floor, Cldn4 | 21 vs 28 | MWU Δ median | 4.349 | 1.21e-05 |
| GSEA GO TJ on Tacstd2 ranking | 89 | NES | 1.905 | 0.0010 |
| GSEA Reactome keratinization on Tacstd2 ranking | 145 | NES | 2.345 | 0.0010 |

TROP2-positive vs TROP2-floor, TJ-core score: median 0.390 vs -0.328 (Δ = 0.718, two-sided Mann–Whitney p = 0.0003).
Cldn4: median 4.354 vs 0.005 (Δ = 4.349, p = 1.21e-05).
GO tight junction (GO:0070160) score vs Tacstd2: Spearman ρ = 0.677, p = 9.48e-08, n = 49.
Reactome keratinization score vs Tacstd2: Spearman ρ = 0.684, p = 6.15e-08, n = 49.

Lung vitro ranking (Tacstd2, high → low):

| line | class | n_baseline | Tacstd2 | Cldn4 | TJ-core |
|---|---|---:|---:|---:|---:|
| MLE12 (vitro) | lung_other | 4 | 5.917 | 0.000 | 0.373 |
| LLC (vitro) | LLC | 19 | 0.000 | 0.000 | -0.401 |
| CMT-167 (vivo naive) | CMT-167 | 3 | 0.527 | 1.863 | 0.765 |
| LLC (vivo naive) | LLC | 31 | 0.183 | 0.057 | -0.611 |

## 4. TISMO ICB-resistant vs sensitive

ICB-treated samples with TISMO `ICB_study` ∈ {R, NR}: **n_NR = 169**, **n_R = 265**.
On those samples, Tacstd2 NR vs R: median 0.417 vs 0.479 (Δ = -0.062, p = 0.0444). TJ-core: Δ = 0.136, p = 0.0048.
Model-level (cell line with ≥3 ICB-labeled samples; resistant-leaning = NR fraction ≥ 0.6; sensitive-leaning = NR fraction ≤ 0.4): **n_resistant = 3**, **n_sensitive = 9**, mixed/underpowered held out.
Baseline TJ-core, resistant- vs sensitive-leaning models: Δ = 0.866, p = 0.2091 (n = 3 vs 9).

## 5. Public KL vs KP cell lines — GSE274352

Empty-vector libraries: **KL n = 3** (KL1/KL2/KL3) vs **KP n = 3** (KP1/KP2/KP3). Values = log2(normalized count + 1). Unit = line (one empty library per line).

| gene / score | KL median | KP median | Δ (KL−KP) | MWU p |
|---|---:|---:|---:|---:|
| Tacstd2 | 5.223 | 6.868 | -1.645 | 0.1000 |
| Cldn4 | 9.709 | 9.472 | 0.237 | 1.0000 |
| tj_core_score | -0.055 | -0.589 | 0.535 | 0.7000 |
| go_tj_score | 0.062 | -0.141 | 0.203 | 0.1000 |
| reactome_keratin_score | 0.172 | -0.295 | 0.467 | 0.1000 |

GSEA on Welch *t* (KL empty − KP empty), GO tight junction: NES = 1.57, p = 0.0090, hits = 77. Reactome keratinization: NES = 2.06, p = 0.0020, hits = 57.

## 6. Additional public GEMM cell-line RNA-seq

**GSE167381** LKB1-restorable LR1/LR2 (vehicle = LKB1-off, 4-OHT = restored). n_off = 4, n_on = 4 libraries (2 lines × 2 technical replicates). Tacstd2 LKB1-off vs on: Δ = 0.122, p = 1.0000. Cldn4: Δ = -0.781, p = 0.3429. TJ-core: Δ = -0.126, p = 0.1143. LU1/LU2 are constitutively unrestorable LKB1-null lines and are tabulated separately.

**GSE295685** is KL-only (KRAS-G12D / STK11-null), n = 2 DMSO + 2 TNG260. Mean log2(TPM+1) Tacstd2 (DMSO) = 3.850; Cldn4 = 10.408. TNG260 − DMSO: Tacstd2 Δ = 1.220, p = 0.3333; Cldn4 Δ = 0.362, p = 0.3333; TJ-core Δ = 0.686, p = 0.3333. No KP arm is deposited, so this series does not support a KL vs KP test.

## 7. GSEA / ssGSEA (TISMO 49 lines)

Ranking = Spearman ρ of each gene with Tacstd2 across the 49 line medians (positive = co-expressed with TROP2).

| set | n mapped | hits | ES | NES | p |
|---|---:|---:|---:|---:|---:|
| TJ_CORE | 9 | 9 | 0.791 | 1.889 | 0.0010 |
| GO_CC_TIGHT_JUNCTION | 89 | 89 | 0.562 | 1.905 | 0.0010 |
| GO_BP_TIGHT_JUNCTION_ASSEMBLY | 50 | 50 | 0.572 | 1.859 | 0.0010 |
| REACTOME_KERATINIZATION | 145 | 145 | 0.670 | 2.345 | 0.0010 |
| GO_BP_KERATINOCYTE_DIFFERENTIATION | 37 | 37 | 0.636 | 1.994 | 0.0010 |

Leading-edge genes are in `tables/gsea_tismo_tacstd2_prerank.tsv`.

## 8. Gene lists

- **TJ-core (pre-specified):** Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, Ocln, Tjp1, Nectin2.
- **GO Cellular Component Tight Junction (GO:0070160)** and **GO BP Tight Junction Assembly (GO:0120192):** Enrichr GO 2023, human symbols mapped to mouse by case-insensitive match (`tables/gene_sets_mapped.tsv`).
- **Reactome Keratinization (R-HSA-6805567)** used because Enrichr GO BP 2021/2023/2025 do not contain a `Keratinization` term; keratinocyte-differentiation (GO:0030216) is reported alongside it.

## 9. Methods (one paragraph)

TISMO RDS matrices were read with pyreadr. Vitro tests use `Baseline==1` sample medians per cell line (n=49). Scores are the mean of per-gene z-scores across lines (TJ-core; GO:0070160; Reactome keratinization); a rank-percentile ssGSEA-like mean is reported as a companion. Primary TROP2 split is Tacstd2 > 0 versus Tacstd2 = 0 (the median is 0 because 28/49 lines sit at the floor); tertiles and a median split among expressed lines are in the tables. Group tests are two-sided Mann–Whitney; correlations are Spearman. GSEA is a weighted Kolmogorov–Smirnov enrichment on the Tacstd2-correlation ranking (or KL−KP Welch *t*), 1,000 permutations. GSE274352 / GSE295685 / GSE167381 use depositor processed matrices; counts/TPM were log2(x+1). Honest n is cell lines (TISMO, GSE274352) or libraries (GSE167381 technical replicates, GSE295685 biological replicates).

## 10. Figures

- `figures/fig1_tismo_vitro_trop2_tj.png` — 49-line Tacstd2 rank, Tacstd2 vs TJ-core, high/low boxes, TJ-core heatmap.
- `figures/fig2_lung_and_icb.png` — lung syngeneic ranking; ICB model and sample contrasts.
- `figures/fig3_gse274352_kl_vs_kp.png` — KL vs KP empty-vector cell lines.
- `figures/fig4_extra_public_kl_lines.png` — GSE167381 LKB1 restoration; GSE295685 KL ± TNG260.
- `figures/fig5_gsea_tismo.png` — prerank barcode plots on the 49-line Tacstd2 ranking.

## 11. What a paper can use

In TISMO's 49 untreated syngeneic cell lines, Tacstd2 correlates with the pre-specified TJ-core score (ρ=0.63, p=9.84e-07, n=49); TROP2-positive lines (n=21) have higher TJ-core than TROP2-floor lines (n=28; Δ=0.72, p=0.0003). GO tight-junction and Reactome keratinization GSEA on the Tacstd2 ranking are NES=1.90 and 2.35 (both permutation p=0.001). True KL/KP lines are absent from TISMO; lung coverage is LLC (not KL; Tacstd2/Cldn4 at floor), MLE12 (Tacstd2-high, Cldn4 floor), and CMT-167 (vivo naive; Tacstd2 0.53, Cldn4 1.86). ICB-treated TISMO samples: NR (n=169) have higher TJ-core than R (n=265; Δ=+0.14, p=0.0048). In GSE274352 GEMM cell lines (empty vector; n=3 KL vs 3 KP) single-gene MWU is underpowered (Tacstd2 Δ=−1.65, p=0.10; Cldn4 Δ=+0.24, p=1.0; TJ-core Δ=+0.53, p=0.70), while prerank GSEA of KL−KP *t* enriches GO tight junction (NES=1.57, p=0.009) and Reactome keratinization (NES=2.06, p=0.002).

