# Mouse lung ICI scRNA pool — Tacstd2 / Cldn4

**Additive.** User A4 TISMO (49/64 Tacstd2 up after ICB, Wilcoxon p=5.8e-5) is **taken as given** and is not re-audited. This slice is public **mouse lung ICI scRNA**, not TISMO bulk.

## Verdict

No public mouse-lung ICI scRNA series in this pool has **n≥2 biological libraries per arm** of epithelial/tumor Tacstd2 or Cldn4 after ICB vs a true no-ICB control. Several series **lack epithelial/tumor cells** (CD45+/CD3+ sorts). GSE283827 (SCLC ± aPD-1) cannot be scored: MTX has 32,589 genes and **no features.tsv**.

Paired same-library epithelial vs T/NK: Tacstd2 is higher in epithelial in **10/13** samples (Wilcoxon p=0.068 — not a firm restriction claim). Cldn4 is higher in epithelial in **11/13** (p=9.8e-4). GSE157881 deposited CD45− Tacstd2 is **lower** than CD45+ T/NK. This is not an ICB-response claim and does **not** reproduce TISMO 49/64.

## Catalog (scored vs documented)

| Series | Model | ICB in *these* libraries | Epithelial/tumor | T/NK | Tacstd2 / Cldn4 |
|---|---|---|---|---|---|
| GSE157881 | HKP1 lung, 0 vs 4 Gy | **No** (RT only) | CD45− present | CD45+ present | scored |
| GSE157882 | HKP1, club-cell DT after 4 Gy | **No** | **ABSENT** (CD45+ only) | present | scored in T/NK |
| GSE176091 | VC-LUAD, CA170 vs PBS | CA170 (VISTA/PD-L1) | **ABSENT** (CD45+ TIL) | present; n=2 vs 2 | scored in T/NK |
| GSE267557 | MC38-bearing lung, young vs aged | both aPD-1 | **ABSENT** (CD45+) | present | not downloaded (1.2 GB); no control arm |
| GSE268525 | LLC1-sgLkb1 lung met | ICI vs RT+ICI (**no untreated**) | scored if Epcam/Cdh1+Krt8 | scored | n=1 vs 1 |
| GSE283827 | Rb/p53 SCLC ± aPD-1 ± ERBB2i | yes (Ly16 aPD-1, Ly29 vehicle) | authors: epithelial | authors: immune | **not scored — no gene names** |
| GSE303943 | CMT167R subcutaneous | **No** (PKCi vs solvent) | scored | scored | not ICB; subcutaneous |
| GSE133604 leftover | KP ± aPD-1 | yes | scored | scored | n=1 vs 1; features = mm10-3.0 31053 |
| GSE129297 leftover | SCLC GEMM ± aPD-1 | yes | scored | scored | n=1 vs 1; unfiltered 10x, UMI≥500 |
| GSE232730 leftover | KPL-3M CD45+ | yes | **ABSENT** | present | n=1 vs 1 |
| GSE222158 leftover | FVB lung CD45+ | yes | **ABSENT** | present | n=1 vs 1 |
| GSE275877 leftover | LKR13 ICI R vs acquired-NR | yes | unknown | unknown | **skip**: 31,053 × 33.97M unfiltered barcodes |

## Sample-level epithelial / T/NK Tacstd2

### Epithelial / tumor

| Series | Sample | Arm | n cells | Tacstd2 mean | Tacstd2 %pos | Cldn4 mean | Cldn4 %pos |
|---|---|---|---|---|---|---|---|
| GSE157881 | 0Gy_CD45neg | control_0Gy | 5979 | 0.043 | 0.5 | 0.008 | 0.1 |
| GSE157881 | 4Gy_CD45neg | RT_4Gy | 4313 | 0.058 | 0.8 | 0.022 | 0.3 |
| GSE268525 | Flox_ICI | ICI | 345 | 0.126 | 14.2 | 0.012 | 1.4 |
| GSE268525 | Flox_RI | RT_ICI | 173 | 0.125 | 15.0 | 0.016 | 2.3 |
| GSE268525 | CKO_RI | Sfrp2cKO_RT_ICI | 238 | 0.103 | 11.8 | 0.025 | 2.9 |
| GSE303943 | CMT167R_CON | solvent | 45 | 0.000 | 0.0 | 0.000 | 0.0 |
| GSE303943 | CMT167R_PKCi | PKCi | 94 | 0.026 | 3.2 | 0.000 | 0.0 |
| GSE133604 | KP_Ctrl_IgG | control | 327 | 0.020 | 2.4 | 0.014 | 1.5 |
| GSE133604 | KP_Ctrl_aPD1 | aPD1 | 351 | 0.019 | 2.3 | 0.004 | 0.6 |
| GSE133604 | KP_Asf1aKO_IgG | Asf1aKO | 404 | 0.025 | 3.0 | 0.005 | 0.5 |
| GSE133604 | KP_Asf1aKO_aPD1 | Asf1aKO_aPD1 | 361 | 0.013 | 1.9 | 0.010 | 0.8 |
| GSE129297 | SCLC_Ctrl | control | 1063 | 0.003 | 0.2 | 0.124 | 15.1 |
| GSE129297 | SCLC_aPD1 | aPD1 | 868 | 0.004 | 0.5 | 0.106 | 12.6 |
| GSE129297 | SCLC_YKL | YKL | 1001 | 0.014 | 1.5 | 0.119 | 13.9 |
| GSE129297 | SCLC_combo | YKL_aPD1 | 918 | 0.013 | 1.6 | 0.113 | 13.6 |

### T/NK

| Series | Sample | n cells | Tacstd2 mean | Tacstd2 %pos |
|---|---|---|---|---|
| GSE157881 | 0Gy_CD45pos | 2073 | 0.077 | 1.1 |
| GSE157881 | 4Gy_CD45pos | 4020 | 0.099 | 1.4 |
| GSE157882 | PBS_CD45pos | 1998 | 0.192 | 2.7 |
| GSE157882 | DT_CD45pos | 1952 | 0.165 | 2.2 |
| GSE176091 | CD45_Control_1 | 3057 | 0.017 | 1.9 |
| GSE176091 | CD45_Control_2 | 3515 | 0.013 | 1.5 |
| GSE176091 | CD45_CA170_1 | 3682 | 0.008 | 1.0 |
| GSE176091 | CD45_CA170_2 | 1944 | 0.007 | 0.9 |
| GSE268525 | Flox_ICI | 1433 | 0.004 | 0.6 |
| GSE268525 | Flox_RI | 946 | 0.000 | 0.0 |
| GSE268525 | CKO_RI | 1335 | 0.004 | 0.6 |
| GSE303943 | CMT167R_CON | 172 | 0.012 | 1.7 |
| GSE303943 | CMT167R_PKCi | 280 | 0.005 | 0.4 |
| GSE133604 | KP_Ctrl_IgG | 2861 | 0.013 | 1.7 |
| GSE133604 | KP_Ctrl_aPD1 | 3484 | 0.015 | 1.9 |
| GSE133604 | KP_Asf1aKO_IgG | 1844 | 0.013 | 1.8 |
| GSE133604 | KP_Asf1aKO_aPD1 | 2778 | 0.017 | 2.4 |
| GSE129297 | SCLC_Ctrl | 1289 | 0.003 | 0.4 |
| GSE129297 | SCLC_aPD1 | 956 | 0.003 | 0.4 |
| GSE129297 | SCLC_YKL | 1210 | 0.021 | 2.8 |
| GSE129297 | SCLC_combo | 1652 | 0.007 | 0.9 |
| GSE232730 | KPL_CD45_Ctrl | 2254 | 0.069 | 7.9 |
| GSE232730 | KPL_CD45_aPD1 | 2551 | 0.060 | 7.3 |
| GSE222158 | CD45_Ctl | 1305 | 0.078 | 8.7 |
| GSE222158 | CD45_PD1 | 1319 | 0.064 | 7.6 |

## Contrasts (honest n / p)

| Series | Contrast | Compartment | Gene | n_a | n_b | Δ | Welch p | MWU p | Note |
|---|---|---|---|---|---|---|---|---|---|
| GSE157881 | 4Gy_vs_0Gy | epithelial | Tacstd2 | 1 | 1 | 0.015 | NA (n<2) | NA (n<2) | RT not ICB; n=1 vs 1 — p not computed |
| GSE157881 | 4Gy_vs_0Gy | epithelial | Cldn4 | 1 | 1 | 0.014 | NA (n<2) | NA (n<2) | RT not ICB; n=1 vs 1 |
| GSE157881 | 4Gy_vs_0Gy_CD45pos | tnk | Tacstd2 | 1 | 1 | 0.022 | NA (n<2) | NA (n<2) | immune sort |
| GSE157881 | 4Gy_vs_0Gy_CD45pos | tnk | Cldn4 | 1 | 1 | -0.006 | NA (n<2) | NA (n<2) | immune sort |
| GSE157881 | CD45neg_vs_CD45pos_0Gy | epithelial_minus_tnk | Tacstd2 | 1 | 1 | -0.034 | NA (n<2) | NA (n<2) | author CD45 sort; n=1 library each |
| GSE157881 | CD45neg_vs_CD45pos_0Gy | epithelial_minus_tnk | Cldn4 | 1 | 1 | 0.001 | NA (n<2) | NA (n<2) | author CD45 sort; n=1 library each |
| GSE157881 | CD45neg_vs_CD45pos_4Gy | epithelial_minus_tnk | Tacstd2 | 1 | 1 | -0.041 | NA (n<2) | NA (n<2) | author CD45 sort; n=1 library each |
| GSE157881 | CD45neg_vs_CD45pos_4Gy | epithelial_minus_tnk | Cldn4 | 1 | 1 | 0.022 | NA (n<2) | NA (n<2) | author CD45 sort; n=1 library each |
| GSE176091 | CA170_vs_control | tnk | Tacstd2 | 2 | 2 | -0.007 | 0.17 | 0.333 | CD45+ TILs; epithelial ABSENT |
| GSE176091 | CA170_vs_control | tnk | Cldn4 | 2 | 2 | -0.000 | 0.977 | 1 | CD45+ TILs; epithelial ABSENT |
| GSE268525 | RTICI_vs_ICI | epithelial | Tacstd2 | 1 | 1 | -0.001 | NA (n<2) | NA (n<2) | no untreated arm; n=1 vs 1 |
| GSE268525 | RTICI_vs_ICI | epithelial | Cldn4 | 1 | 1 | 0.004 | NA (n<2) | NA (n<2) | no untreated arm; n=1 vs 1 |
| GSE268525 | RTICI_vs_ICI | tnk | Tacstd2 | 1 | 1 | -0.004 | NA (n<2) | NA (n<2) |  |
| GSE268525 | RTICI_vs_ICI | tnk | Cldn4 | 1 | 1 | 0.000 | NA (n<2) | NA (n<2) |  |
| GSE303943 | PKCi_vs_solvent | epithelial | Tacstd2 | 1 | 1 | 0.026 | NA (n<2) | NA (n<2) | subcutaneous; not ICB vs control; n=1 vs 1 |
| GSE303943 | PKCi_vs_solvent | epithelial | Cldn4 | 1 | 1 | 0.000 | NA (n<2) | NA (n<2) | subcutaneous; n=1 vs 1 |
| GSE303943 | PKCi_vs_solvent | tnk | Tacstd2 | 1 | 1 | -0.007 | NA (n<2) | NA (n<2) |  |
| GSE303943 | PKCi_vs_solvent | tnk | Cldn4 | 1 | 1 | 0.002 | NA (n<2) | NA (n<2) |  |
| GSE133604 | KP_aPD1_vs_IgG | epithelial | Tacstd2 | 1 | 1 | -0.001 | NA (n<2) | NA (n<2) | leftover; n=1 vs 1 |
| GSE133604 | KP_aPD1_vs_IgG | epithelial | Cldn4 | 1 | 1 | -0.010 | NA (n<2) | NA (n<2) | leftover; n=1 vs 1 |
| GSE133604 | KP_aPD1_vs_IgG | tnk | Tacstd2 | 1 | 1 | 0.002 | NA (n<2) | NA (n<2) | leftover |
| GSE133604 | KP_aPD1_vs_IgG | tnk | Cldn4 | 1 | 1 | 0.000 | NA (n<2) | NA (n<2) | leftover |
| GSE129297 | SCLC_aPD1_vs_Ctrl | epithelial | Tacstd2 | 1 | 1 | 0.001 | NA (n<2) | NA (n<2) | leftover; n=1 vs 1 |
| GSE129297 | SCLC_aPD1_vs_Ctrl | epithelial | Cldn4 | 1 | 1 | -0.018 | NA (n<2) | NA (n<2) | leftover; n=1 vs 1 |
| GSE129297 | SCLC_aPD1_vs_Ctrl | tnk | Tacstd2 | 1 | 1 | 0.000 | NA (n<2) | NA (n<2) | leftover |
| GSE129297 | SCLC_aPD1_vs_Ctrl | tnk | Cldn4 | 1 | 1 | 0.000 | NA (n<2) | NA (n<2) | leftover |
| GSE232730 | aPD1_vs_control | tnk | Tacstd2 | 1 | 1 | -0.009 | NA (n<2) | NA (n<2) | CD45+ only; leftover; n=1 vs 1 |
| GSE232730 | aPD1_vs_control | tnk | Cldn4 | 1 | 1 | 0.000 | NA (n<2) | NA (n<2) | CD45+ only |
| GSE222158 | aPD1_vs_control | tnk | Tacstd2 | 1 | 1 | -0.014 | NA (n<2) | NA (n<2) | CD45+ only; leftover; n=1 vs 1 |
| GSE222158 | aPD1_vs_control | tnk | Cldn4 | 1 | 1 | 0.000 | NA (n<2) | NA (n<2) | CD45+ only |
| POOL | epithelial_vs_tnk_paired | epithelial_minus_tnk | Cldn4 | 13 | 13 | 0.042 | NA (n<2) | 0.000977 | Wilcoxon signed-rank W=1.0; samples with ≥20 cells each compartment |
| POOL | epithelial_vs_tnk_paired | epithelial_minus_tnk | Tacstd2 | 13 | 13 | 0.029 | NA (n<2) | 0.0681 | Wilcoxon signed-rank W=19.0; samples with ≥20 cells each compartment |

## Combined direction (epithelial ICB-ish + epi vs T/NK)

Only **GSE176091** has n=2 vs 2, and that is **T/NK in a CD45+ sort** (epithelial absent). All epithelial ICB-vs-control tests are **n=1 vs 1** (leftover GSE133604, GSE129297). Sign is reported; p is not.

| Series | Contrast | Compartment | Gene | Direction (treated − ref) | n | p |
|---|---|---|---|---|---|---|
| GSE133604 leftover | KP aPD-1 vs IgG | epithelial | Tacstd2 | down (tiny) | 1 vs 1 | NA |
| GSE133604 leftover | KP aPD-1 vs IgG | epithelial | Cldn4 | down | 1 vs 1 | NA |
| GSE129297 leftover | SCLC aPD-1 vs Ctrl | epithelial | Tacstd2 | up (tiny) | 1 vs 1 | NA |
| GSE129297 leftover | SCLC aPD-1 vs Ctrl | epithelial | Cldn4 | down | 1 vs 1 | NA |
| GSE176091 | CA170 vs PBS | T/NK (CD45+; epi **ABSENT**) | Tacstd2 | down | 2 vs 2 | Welch 0.17; MWU 0.33 |
| GSE176091 | CA170 vs PBS | T/NK | Cldn4 | tie | 2 vs 2 | Welch 0.98; MWU 1 |
| GSE268525 | RT+ICI vs ICI (no untreated) | epithelial | Tacstd2 | tie/down | 1 vs 1 | NA |
| GSE268525 | RT+ICI vs ICI | epithelial | Cldn4 | up | 1 vs 1 | NA |
| GSE157881 | CD45− vs CD45+ (RT, not ICB) | epi − T/NK | Tacstd2 | **down** (epi < T/NK) | 1 vs 1 | NA |
| POOL | same-library epi vs T/NK | epithelial − T/NK | Tacstd2 | 10/13 epi > T/NK | 13 paired | Wilcoxon 0.068 |
| POOL | same-library epi vs T/NK | epithelial − T/NK | Cldn4 | 11/13 epi > T/NK | 13 paired | Wilcoxon 9.8e-4 |

Paired epithelial vs T/NK Tacstd2: **10/13** samples epithelial > T/NK (samples with ≥20 cells in both); Wilcoxon signed-rank p=0.068. Cldn4: **11/13**, p=9.8e-4.

## Methods (short)

- Public GEO processed MTX/TSV only. Tacstd2 = ENSMUSG00000051397; Cldn4 = ENSMUSG00000047501.
- Unsorted 10x: epithelial/tumor = Epcam+ or (Cdh1+ and Krt8+) or Ascl1/Chga/Insm1+; T/NK = Cd3d/e, Cd8a, Nkg7, Ncr1 and not epithelial.
- GSE157881/882 deposited log-like values (used as-is). 10x MTX: log1p(UMI).
- GSE129297/GSE133604 gene names: Cell Ranger mm10-3.0 31,053-gene table from GSE275877 features (same n_genes).
- Welch / MWU only if n≥2 samples/arm. No invented R vs NR labels.

## 中文摘要

TISMO 用户 A4（49/64 Tacstd2 在 ICB 后升高，p=5.8e-5）**视为已知，不重算**。本切片是公开小鼠肺 ICI **单细胞**，不是 TISMO bulk。
候选集中：GSE157881/882 是放疗/club 细胞清除而非 ICB；GSE176091 / GSE267557 / GSE232730 / GSE222158 为 CD45+ 或 CD3+，**无上皮/肿瘤细胞**；GSE283827 有 aPD-1 但 **未提供 features.tsv**，无法对 Tacstd2/Cldn4 诚实计分；GSE268525 无未治疗对照；GSE303943 为皮下 PKCi。
有上皮且有 ICB vs 对照的只剩 leftover GSE133604（KP）与 GSE129297（SCLC），均为 **n=1 vs 1**，只报方向，不算 p。配对上皮 vs T/NK：Tacstd2 10/13 上皮更高（Wilcoxon p=0.068）；Cldn4 11/13（p=9.8e-4）。GSE157881 CD45− Tacstd2 反而低于 CD45+ T/NK。不能外推 TISMO 49/64。
