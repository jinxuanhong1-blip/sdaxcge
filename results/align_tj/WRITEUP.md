# USER-ALIGN: TROP2 / tight-junction public replication

All products live under `results/align_tj/`. Scripts are in `scripts/`, numbers in `tables/`. Large GEO/Xena downloads stay in `data/` and are not committed.

This write-up states what public data do and do not support. Thresholds were fixed before looking at the user gene list in the intersection step (`r ≥ 0.10` and `r ≥ 0.20`, BH-FDR `< 0.05`).

## Question

Does TROP2-high NSCLC carry a keratinization / skin-barrier / tight-junction program with EMT down, and do the user genes **CLDN1, CLDN4, CLDN7, F11R, PARD3** plus the TFs **ELF3, GRHL1, KLF4, TFAP2A** co-express with **TACSTD2 / CLDN4** in independent public human and mouse data?

## Datasets

| Resource | Use | n | Note |
|---|---|---|---|
| TCGA LUAD+LUSC HiSeqV2 (UCSC Xena, log2 norm) | bulk DE, preranked GSEA, Spearman vs TACSTD2 | 1017 primary tumors (515 LUAD, 502 LUSC) | TROP2-high tertile is LUSC-enriched (222 LUSC / 117 LUAD) vs low (136 / 203) |
| GSE207422 bulk log2 TPM | same-study bulk coexpression | 24 tumors | **not** an independent cohort |
| GSE207422 scRNA UMI | malignant/epithelial cells (marker-called; no author per-cell labels) | running | epithelial = EPCAM+ / PTPRC-low |
| GSE131907 scRNA (Kim et al. 2020) | author-annotated tLung + tL/B epithelial cells | running | third NSCLC scRNA layer |
| COXPRESdb v8 Hsa-u / Mmu-u | corpus-wide coexpression (mutual rank) | all-tissue RNA-seq+array | **not lung-restricted** |
| Tabula Muris FACS lung (GSE109774) | mouse lung epithelial Smart-seq2 | 138 epithelial cells | small n; Cldn4 sparse (12%) |

## 1. TCGA TROP2-high vs low GSEA

Ranking: Welch *t* (high tertile minus low tertile) on 1017 tumors. Libraries: GO BP/CC 2021, Hallmark 2020, KEGG 2021.

| Theme | Lead term | NES | FDR | Direction in TROP2-high |
|---|---|---|---|---|
| Skin barrier | epidermis development (GO:0008544) | 3.35 | 0 | **up** |
| Skin barrier | skin development (GO:0043588) | 3.18 | 0 | **up** |
| Skin barrier | cornified envelope (GO:0001533) | 3.18 | 0 | **up** |
| Skin barrier | establishment of skin barrier (GO:0061436) | 2.54 | 0 | **up** |
| Keratinization | keratinocyte differentiation (GO:0030216) | 2.85 | 0 | **up** |
| Tight junction | cell-cell junction (GO:0005911) | 1.98 | 0.032 | **up** |
| Tight junction | Hallmark Apical Junction | 1.83 | 0.083 | up, FDR borderline |
| Tight junction | tight junction assembly (GO:0120192) | 1.78 | 0.108 | up, not FDR-significant |
| Tight junction | KEGG Tight junction | 1.67 | 0.185 | up, not FDR-significant |
| EMT | Hallmark EMT | +1.30 | 0.50 | **not down**; weak up, FDR ns |

**Honest GSEA call.** Keratinization and skin/epidermal barrier are strongly and specifically enriched in TROP2-high tumors. Tight-junction / apical-junction terms are directionally up; only the broad cell-cell junction set is FDR `< 0.05`. EMT-down does **not** replicate: Hallmark EMT is a weak positive NES with FDR 0.50.

**Confounder.** The high tertile is squamous-enriched, and squamous tumors express keratinization programs. Within-histology coexpression (below) is the cleaner test of TACSTD2 coupling.

## 2. User TJ genes vs TACSTD2 (human)

Spearman *r* (BH-FDR on Spearman *p* for TCGA).

| Gene | TCGA NSCLC r | LUAD-only r | LUSC-only r | GSE207422 bulk r (n=24) |
|---|---:|---:|---:|---:|
| CLDN1 | 0.44 | 0.37 | 0.42 | 0.88 |
| CLDN4 | 0.25 | **0.46** | **0.41** | 0.86 |
| CLDN7 | 0.11 | 0.05 (ns) | **0.40** | 0.83 |
| F11R | 0.30 | 0.18 | 0.25 | 0.87 |
| PARD3 | 0.16 | −0.02 (ns) | 0.15 | 0.66 |

Pooled CLDN4–TACSTD2 (r=0.25) is **weaker than the within-histology values** because TACSTD2 is higher in LUSC while CLDN4 is higher in LUAD. That is a histology-mix artifact, not a weak biology.

**User-list overlap in TCGA (r ≥ 0.10, FDR `< 0.05`):** CLDN1, CLDN4, CLDN7, F11R, PARD3 all pass in pooled NSCLC. **Honest histology split:** CLDN7 and PARD3 fail in LUAD; PARD3 is only weakly positive in LUSC. The robust TJ members are **CLDN1, CLDN4, F11R**. CLDN7 is LUSC-restricted. PARD3 is the weakest and should not be claimed as a pan-NSCLC TROP2 partner.

GSE207422 bulk agrees in sign for every user gene but is the same study as the scRNA layer and n=24, so it is not an independent vote.

## 3. ELF3 / GRHL1 / KLF4 / TFAP2A vs TACSTD2 and CLDN4 (human)

| TF | vs TACSTD2 NSCLC / LUAD / LUSC | vs CLDN4 NSCLC / LUAD / LUSC |
|---|---|---|
| **ELF3** | 0.27 / 0.29 / **0.47** | **0.56 / 0.46 / 0.48** |
| GRHL1 | **0.44 / 0.40 / 0.39** | 0.02 / **0.39 / 0.24** |
| KLF4 | 0.32 / 0.14 / 0.37 | −0.21 / −0.12 / 0.16 |
| TFAP2A | 0.30 / 0.24 / 0.15 | −0.33 / 0.15 / −0.16 |

**ELF3 is the cleanest TF.** It co-expresses with both TACSTD2 and CLDN4 **inside LUAD and inside LUSC**. The ELF3–CLDN4 pair is the strongest TF–TJ correlation in TCGA (r≈0.46–0.56).

**GRHL1** tracks TACSTD2 in both histologies. The pooled GRHL1–CLDN4 r≈0 is a histology-mix artifact; within LUAD/LUSC the pair is positive.

**KLF4 and TFAP2A** track TACSTD2 but **do not** stably track CLDN4 (sign flips across histologies). Do not claim a KLF4/TFAP2A–CLDN4 axis from TCGA.

COXPRESdb human unified network (all tissues, mutual rank; top-1000 neighbors):

| Query → partner | Rank | Mutual rank | logit |
|---|---:|---:|---:|
| TACSTD2 → CLDN4 | 3 | 212 | 6.29 |
| TACSTD2 → CLDN7 | 7 | 460 | 5.15 |
| TACSTD2 → GRHL1 | 40 | 1055 | 3.90 |
| TACSTD2 → ELF3 | 63 | 1445 | 3.41 |
| TACSTD2 → CLDN1 | 79 | 1734 | 3.12 |
| ELF3 → CLDN4 | **2** | 84 | 7.64 |
| ELF3 → TACSTD2 | 75 | 1445 | 3.41 |
| CLDN4 → ELF3 | 3 | 84 | 7.64 |
| CLDN4 → TACSTD2 | 9 | 212 | 6.29 |
| TACSTD2 → PARD3 | not in top 1000 | — | — |
| TFAP2A → CLDN4 | not in top 1000 | — | — |

Corpus-wide coexpression independently ranks **CLDN4 among the nearest neighbors of TACSTD2** and **ELF3 among the nearest neighbors of CLDN4**.

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
| Cldn4 → Elf3 | 4 | 55 | 8.28 |
| Cldn4 → Tacstd2 | 6 | 66 | 8.03 |
| Grhl1 / Klf4 vs Tacstd2 or Cldn4 | not in top 1000 | — | — |
| Pard3 vs Tacstd2 | not in top 1000 | — | — |

Mouse corpus-wide data **replicate the TACSTD2–CLDN4–ELF3 triangle**. They do **not** support Grhl1 or Klf4 as close Tacstd2/Cldn4 neighbors.

### Tabula Muris lung epithelial (n=138 Smart-seq2)

Tissue-matched but underpowered. Tacstd2 in 32% of cells; Cldn4 in 12%; Grhl1/Tfap2a almost off.

| Pair | Spearman r | p |
|---|---:|---:|
| Elf3–Tacstd2 | 0.18 | 0.035 |
| Elf3–Cldn4 | 0.03 | 0.69 |
| Cldn4–Tacstd2 | 0.15 | 0.074 |
| Cldn7–Tacstd2 | 0.27 | 0.001 |
| F11r–Tacstd2 | 0.41 | 5e-7 |
| Pard3–Tacstd2 | −0.01 | 0.91 |
| Grhl1–Tacstd2 | 0.01 | 0.94 |

Lung-restricted mouse data support Tacstd2 with **Cldn7 / F11r** and a weak Elf3–Tacstd2 link. They do **not** have the power to confirm Cldn4 (too sparse) or Grhl1.

## 5. Triple public intersection (NSCLC)

Pre-specified: gene is “TROP2-positive” in a dataset if Spearman r ≥ 0.10 and BH-FDR `< 0.05`.

**Status.** TCGA genome-wide table is done (3,329 genes at r≥0.10 FDR<0.05; 1,212 at r≥0.20). GSE207422 malignant-cell and GSE131907 tumor-epithelial genome-wide Spearman tables are still running. The intersection gene list and the per-user-gene 2/3 and 3/3 calls will be written to `tables/intersection_user_genes.tsv` when both scRNA passes finish.

Until then, the **already-complete independent public layers** that can vote on the user list are:

1. TCGA NSCLC (bulk, histology-split)
2. COXPRESdb human
3. COXPRESdb mouse

**Provisional 3-layer vote (r≥0.10 FDR<0.05 in TCGA **and** partner in COXPRESdb top-1000 of TACSTD2 in human **and** of Tacstd2 in mouse):**

| Gene | TCGA | COXPRESdb human | COXPRESdb mouse | 3-layer |
|---|---|---|---|---|
| CLDN1 | yes (both histologies) | yes (rank 79) | yes (rank 244) | **yes** |
| CLDN4 | yes (both; r≈0.41–0.46) | yes (rank 3) | yes (**rank 1**) | **yes** |
| CLDN7 | LUSC only | yes (rank 7) | yes (rank 58) | partial (fails LUAD) |
| F11R | yes (weaker) | yes (rank 391) | yes (rank 911) | weak yes |
| PARD3 | LUAD no; LUSC weak | **no** | **no** | **no** |
| ELF3 | yes (both; vs CLDN4 r≈0.5) | yes (rank 63; CLDN4 rank 2) | yes (rank 181; Cldn4 rank 3) | **yes** |
| GRHL1 | yes vs TACSTD2 | yes (rank 40) | **no** | human-only |
| KLF4 | yes vs TACSTD2 | yes (rank 852) | **no** | human-only, not vs CLDN4 |
| TFAP2A | yes vs TACSTD2 | yes (rank 171) | yes (rank 85) | yes vs TACSTD2; **not** vs CLDN4 |

GSE207422 / GSE131907 malignant-cell Spearman will be the NSCLC-specific third (and fourth) vote. They can only **narrow** this list (dropout will weaken r); they cannot invent PARD3 support that COXPRESdb already lacks.

## 6. What we would claim vs what we would not

**Supported (public, multi-resource).**

- TROP2-high NSCLC is a keratinization / epidermal-barrier state.
- Tight-junction / apical-junction programs trend up with TROP2; the specific KEGG/GO “tight junction” sets are directionally consistent but only partly FDR-significant in bulk GSEA.
- **CLDN1, CLDN4, ELF3** form a cross-species TACSTD2 neighborhood (TCGA within-histology + COXPRESdb human + COXPRESdb mouse).
- **ELF3–CLDN4** is the strongest TF–TJ edge (TCGA r≈0.5 in both histologies; COXPRESdb human rank 2; mouse rank 3).

**Not supported, or only partial.**

- **EMT down in TROP2-high:** not replicated (Hallmark EMT NES +1.30, FDR 0.50).
- **PARD3** as a TROP2/TJ intersection gene: fails LUAD, absent from COXPRESdb top-1000 in both species.
- **CLDN7** is LUSC-linked, not pan-NSCLC.
- **KLF4 / TFAP2A vs CLDN4:** unstable or absent.
- **GRHL1 vs CLDN4:** real within histology in TCGA, not recovered in mouse COXPRESdb.
- GSE207422 bulk r values (0.66–0.89) should not be quoted as independent replication.

## Methods (short)

- TCGA: primary tumors only (`-01`); tertiles of TACSTD2; Welch *t* + Spearman; gseapy prerank 1000 permutations.
- GSE207422 scRNA: no author cell types; lineage = max z-scored marker set; malignant proxy = epithelial ∩ EPCAM+ ∩ PTPRC < EPCAM. Theme scores = mean log(CP10K+1) of MSigDB/GO sets used in the GSEA.
- GSE131907: author `Cell_type == Epithelial cells` and `Sample_Origin ∈ {tLung, tL/B}`.
- Intersection thresholds fixed in `scripts/09_intersection.py`.
- COXPRESdb API v4, databases `hsa-u` and `mmu-u`, top 1000 neighbors.

## Files

- `tables/tcga_gsea_themes.tsv` — theme-filtered GSEA
- `tables/tcga_coexpr_targets.tsv`, `tables/tcga_coexpr_by_histology.tsv`
- `tables/tcga_tf_cldn4_trop2_spearman_matrix.tsv`
- `tables/coxpresdb_target_hits.tsv`
- `tables/tabula_muris_lung_coexpr.tsv`
- `tables/gse207422_bulk_coexpr_targets.tsv`
- `tables/intersection_*.tsv` — filled after scRNA genome-wide passes
