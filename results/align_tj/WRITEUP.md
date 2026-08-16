# USER-ALIGN: TROP2 / tight-junction public replication

All products live under `results/align_tj/`. Scripts in `scripts/`, numbers in `tables/`. Large GEO/Xena downloads stay in `data/` and are not committed.

Thresholds were fixed before scoring the user list: Spearman **r ≥ 0.10** (lenient) or **r ≥ 0.20** (strict), BH-FDR `< 0.05`.

## Question

Does TROP2-high NSCLC carry a keratinization / skin-barrier / tight-junction program with EMT down, and do **CLDN1, CLDN4, CLDN7, F11R, PARD3** plus **ELF3, GRHL1, KLF4, TFAP2A** co-express with **TACSTD2 / CLDN4** in independent public human and mouse data?

## Datasets

| Resource | Use | n | Note |
|---|---|---|---|
| TCGA LUAD+LUSC HiSeqV2 (UCSC Xena) | bulk DE, GSEA, Spearman | 1017 primary tumors (515 LUAD, 502 LUSC) | high tertile is LUSC-enriched (222/117) vs low (136/203) |
| GSE207422 scRNA UMI | marker-called malignant/epithelial | 10,183 EPCAM+ / PTPRC-low epithelial cells (of 92,330) | no author per-cell labels; no CNV |
| GSE131907 scRNA (Kim et al. 2020) | author-annotated tumor epithelial | 13,852 cells (tLung tS1–tS3 + tL/B malignant) | third NSCLC layer |
| GSE207422 bulk log2 TPM | same-study bulk check | 24 tumors | **not** independent |
| COXPRESdb v8 Hsa-u / Mmu-u | corpus-wide mutual rank | all tissues | **not lung-restricted** |
| Tabula Muris FACS lung (GSE109774) | mouse lung epithelial Smart-seq2 | 138 cells | Cldn4 in only 12% of cells |

## 1. TROP2-high vs low GSEA (keratinization, barrier, TJ, EMT)

Ranking: Welch *t*, high vs low TACSTD2 tertile. Libraries: GO BP/CC 2021, Hallmark 2020, KEGG 2021.

| Theme | Cohort | Lead term | NES | FDR |
|---|---|---|---:|---:|
| Skin barrier | NSCLC pooled | epidermis development | 3.35 | 0 |
| Skin barrier | **LUAD only** | skin development | 2.64 | 0 |
| Skin barrier | **LUSC only** | epidermis development | 3.17 | 0 |
| Keratinization | NSCLC pooled | keratinocyte differentiation | 2.85 | 0 |
| Keratinization | **LUAD only** | keratinocyte differentiation | 2.38 | 0 |
| Keratinization | **LUSC only** | keratinocyte differentiation | 2.92 | 0 |
| Tight junction | NSCLC pooled | cell-cell junction | 1.98 | 0.032 |
| Tight junction | **LUAD only** | Hallmark Apical Junction / KEGG Tight junction | 2.10 / 2.09 | **0.004 / 0.005** |
| Tight junction | **LUSC only** | tight junction assembly | 2.18 | **0.004** |
| EMT | NSCLC pooled | Hallmark EMT | **+1.30** | 0.50 |
| EMT | LUAD only | Hallmark EMT | **+1.83** | 0.056 |
| EMT | LUSC only | Hallmark EMT | **−1.52** | 0.24 |

Single-cell module scores (mean log CP10K+1 of the same GO/Hallmark sets) in TROP2-high vs low malignant cells:

| Theme | GSE207422 Δ (high−low) | GSE131907 Δ |
|---|---:|---:|
| keratinization | +0.14 (p ~ 10^−296) | +0.06 |
| skin barrier (cornified envelope) | +0.28 | +0.04 |
| tight junction (GO CC / KEGG) | +0.09 / +0.07 | +0.06 / +0.06 |
| Hallmark EMT | **+0.025** (up) | **+0.037** (up) |

**GSEA call.** Keratinization and epidermal-barrier programs are up in TROP2-high tumors **inside LUAD and inside LUSC**, so this is not a squamous-mix artifact. Tight-junction / apical-junction sets are FDR-significant **within each histology** (stronger than the pooled NSCLC run). EMT-down does **not** replicate as a pan-NSCLC finding: pooled and LUAD NES are positive; LUSC is directionally down but FDR 0.24. Both scRNA malignant compartments also score EMT **up**, not down, in TROP2-high cells.

## 2. Triple NSCLC intersection (TCGA ∩ GSE207422 malig ∩ GSE131907 tumor epithelial)

Pre-specified: gene is TROP2-positive in a dataset if Spearman r ≥ threshold and BH-FDR `< 0.05`.

| Call | TCGA | GSE207422 | GSE131907 | pairwise TCGA∩207422 | **triple** | ≥2 of 3 |
|---|---:|---:|---:|---:|---:|---:|
| r > 0, FDR<0.05 | 4512 | 4838 | 11569 | 1911 | 1629 | 5042 |
| **r ≥ 0.10, FDR<0.05** | 3329 | 1458 | 1158 | 754 | **228** | 1085 |
| **r ≥ 0.20, FDR<0.05** | 1212 | 362 | 60 | 190 | **20** | 211 |

The strict 20-gene triple (TACSTD2 itself plus 19 partners) is an epithelial/barrier set, not a random draw:

`KRT19, PERP, RND3, JUP, CLDN4, EMP1, SERINC2, SPINT1, ACTN4, SPINT2, ELF3, CD9, HES1, H3F3B, NET1, TM4SF1, LIPH, MAFF, LGALS3`

**User genes in that strict triple: CLDN4 and ELF3 only.**

### User-list membership (honest)

| Gene | TCGA r | GSE207422 r | GSE131907 r | triple r≥0.10 | triple r≥0.20 | COXPRESdb human / mouse (in TACSTD2 top-1000) |
|---|---:|---:|---:|---|---|---|
| **CLDN4** | 0.25 (0.46 LUAD / 0.41 LUSC) | **0.43** | **0.32** | **yes** | **yes** | rank 3 / **rank 1** |
| **ELF3** | 0.27 (0.29 / 0.47) | 0.24 | **0.32** | **yes** | **yes** | rank 63 / rank 181 |
| CLDN7 | 0.11 (LUAD ns / LUSC 0.40) | 0.25 | 0.20 | yes | no | rank 7 / rank 58 |
| PARD3 | 0.16 (LUAD ns / LUSC 0.15) | 0.16 | 0.11 | yes | no | **absent / absent** |
| GRHL1 | 0.44 | 0.15 | 0.17 | yes | no | rank 40 / **absent** |
| KLF4 | 0.32 | **0.41** | 0.12 | yes | no | rank 852 / **absent** |
| TFAP2A | 0.30 | 0.19 | 0.12 | yes | no | rank 171 / rank 85 |
| CLDN1 | 0.44 | 0.34 | **0.097** | **no** (2/3; 131907 just under 0.10) | no | rank 79 / rank 244 |
| F11R | 0.30 | **0.044** | 0.20 | **no** (2/3) | no | rank 391 / rank 911 |

**How to read this.** At the pre-specified lenient cut, seven of nine user genes land in the 228-gene triple. That cut is easy to meet in 10k-cell Spearman tests (tiny r still gives tiny p). The **strict** cut, and the COXPRESdb neighborhood, shrink the claim to **CLDN4 + ELF3**. PARD3 scrapes r≥0.10 in all three NSCLC sets but is a weak partner (r 0.11–0.16), fails r≥0.20 everywhere, and is missing from both COXPRESdb top-1000 lists — do not treat it as a TROP2 intersection gene. CLDN1 is biologically strong in TCGA and GSE207422 and in COXPRESdb; it misses the GSE131907 r≥0.10 line by 0.003. CLDN7 is LUSC-linked in bulk. F11R fails GSE207422 malignant cells (r=0.044).

## 3. ELF3 / GRHL1 / KLF4 / TFAP2A vs TACSTD2 and vs CLDN4

| TF | vs TACSTD2 (TCGA LUAD / LUSC / 207422 / 131907) | vs CLDN4 (same four) |
|---|---|---|
| **ELF3** | 0.29 / **0.47** / 0.24 / **0.32** | **0.46 / 0.48 / 0.23 / 0.39** |
| GRHL1 | 0.40 / 0.39 / 0.15 / 0.17 | 0.39 / 0.24 / 0.15 / 0.12 |
| KLF4 | 0.14 / 0.37 / **0.41** / 0.12 | −0.12 / 0.16 / 0.31 / 0.28 |
| TFAP2A | 0.24 / 0.15 / 0.19 / 0.12 | 0.15 / −0.16 / 0.01 ns / 0.08 |

**ELF3 is the only TF that is positive vs both TACSTD2 and CLDN4 in every human NSCLC layer**, including both TCGA histologies. ELF3–CLDN4 is also COXPRESdb human rank 2 (MR 84) and mouse rank 3 (MR 55).

GRHL1 tracks TACSTD2 in human (and is in the r≥0.10 triple) but is absent from the mouse Tacstd2 neighborhood. KLF4 is strong in GSE207422 malignant cells vs both TACSTD2 and CLDN4, weaker and histology-unstable in TCGA vs CLDN4, and absent in mouse COXPRESdb. TFAP2A does not track CLDN4.

GSE207422 bulk (n=24, same study) gives ELF3–TACSTD2 r=0.87 and ELF3–CLDN4 r=0.89. Those numbers are not independent replication.

## 4. Mouse public coexpression

### COXPRESdb mouse unified (all tissues)

| Query → partner | Rank | Mutual rank | logit |
|---|---:|---:|---:|
| **Tacstd2 → Cldn4** | **1** | 66 | 8.03 |
| Tacstd2 → Cldn7 | 58 | 791 | 4.38 |
| Tacstd2 → Tfap2a | 85 | 1049 | 3.95 |
| Tacstd2 → Elf3 | 181 | 1758 | 3.14 |
| Tacstd2 → Cldn1 | 244 | 2207 | 2.77 |
| **Elf3 → Cldn4** | **3** | 55 | 8.28 |
| Elf3 → Tacstd2 | 170 | 1758 | 3.14 |
| Grhl1 / Klf4 vs Tacstd2 or Cldn4 | not in top 1000 | — | — |
| Pard3 vs Tacstd2 | not in top 1000 | — | — |

### Tabula Muris lung epithelial (n=138)

Underpowered. Tacstd2 in 32% of cells; Cldn4 in 12%.

| Pair | r | p |
|---|---:|---:|
| Elf3–Tacstd2 | 0.18 | 0.035 |
| Elf3–Cldn4 | 0.03 | 0.69 |
| Cldn4–Tacstd2 | 0.15 | 0.074 |
| Cldn7–Tacstd2 | 0.27 | 0.001 |
| F11r–Tacstd2 | 0.41 | 5×10^−7 |
| Pard3 / Grhl1 vs Tacstd2 | ~0 | ns |

Lung-restricted mouse data support Tacstd2 with Cldn7/F11r and a weak Elf3–Tacstd2 link. They cannot confirm Cldn4 (too sparse). The mouse claim for the TACSTD2–CLDN4–ELF3 triangle rests on COXPRESdb, not on Tabula Muris.

## 5. What we would claim vs what we would not

**Supported.**

- TROP2-high NSCLC is a keratinization / epidermal-barrier state in LUAD, in LUSC, and in two malignant-cell scRNA cohorts.
- Tight-junction / apical-junction programs are up with TROP2; within-histology GSEA is FDR-significant.
- The **strict triple public intersection** of TROP2-correlated genes is a 20-gene epithelial/barrier set that includes **CLDN4** and **ELF3**.
- **ELF3–CLDN4–TACSTD2** is the cross-species neighborhood (TCGA both histologies, both scRNA cohorts, COXPRESdb human and mouse).

**Not supported, or only partial.**

- **EMT down** as a pan-NSCLC TROP2-high property. Only LUSC is directionally down (FDR 0.24); LUAD and both scRNA sets go the other way.
- **PARD3** as a TROP2 intersection gene. Lenient triple only; weak r; missing from COXPRESdb in both species.
- **CLDN1** as a *triple* member (misses GSE131907 r≥0.10 by 0.003). Still a strong 2/3 + COXPRESdb gene.
- **F11R** in GSE207422 malignant cells (r=0.044).
- **CLDN7** as pan-NSCLC (LUSC-linked in bulk).
- **KLF4 / TFAP2A / GRHL1 vs CLDN4** as a conserved axis. Human-only or sign-unstable; mouse COXPRESdb negative for Grhl1/Klf4.
- GSE207422 bulk r values (0.66–0.89) as independent replication.

## Methods (short)

- TCGA: primary tumors only (`-01A/01`); TACSTD2 tertiles; Welch *t* + Spearman; gseapy prerank, 1000 permutations. Histology-split GSEA uses the same libraries.
- GSE207422 scRNA: lineage = max z-scored marker set; malignant proxy = epithelial ∩ EPCAM+ ∩ PTPRC < EPCAM. Theme score = mean log(CP10K+1).
- GSE131907: `Cell_type == Epithelial cells` and `Sample_Origin ∈ {tLung, tL/B}`.
- Intersection thresholds are hard-coded in `scripts/09_intersection.py`.
- COXPRESdb API v4, `hsa-u` and `mmu-u`, top 1000 neighbors.

## Files

- `WRITEUP.md` — this document
- `tables/tcga_gsea_themes.tsv`, `tables/tcga_luad_gsea_themes.tsv`, `tables/tcga_lusc_gsea_themes.tsv`
- `tables/tcga_coexpr_targets.tsv`, `tables/tcga_coexpr_by_histology.tsv`
- `tables/gse207422_malig_*.tsv`, `tables/gse131907_malig_*.tsv`
- `tables/intersection_user_genes.tsv`, `tables/intersection_sizes.tsv`
- `tables/intersection_genes_r10.tsv` (228 genes), `tables/intersection_genes_r20.tsv` (20 genes)
- `tables/coxpresdb_target_hits.tsv`, `tables/tabula_muris_lung_coexpr.tsv`
