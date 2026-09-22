# Master PPT table — story-fit locked results

Public PRs only. No new fits. No private 8KL matrices. Every coefficient is copied from the cited PR or from the public handoff.

**41 limbs: 29 YES, 12 PARTIAL.** Part 1 is Tacstd2 / TROP2 immune-cold plus Tacstd2 tracking tight junction. Part 2 is the CLDN4 pin. NO, NULL, and opposite results are in the skip list and are not slide support.

**Largest effect** is the largest eligible story-direction number the source PR reported for that limb. **Locked primary** is the pre-specified number when a later sweep is larger. A searched maximum is not a second confirmatory test. Searched p-values are descriptive unless that PR reported a permutation that included the search.

Paste `MASTER_STORY_FIT.pptx` if you need slides. This markdown file is the full table, including caveats that do not fit on a slide.

## Part 1 — Tacstd2 cold

| ID | Limb | Fit | Largest effect | Locked primary | n | PR |
|---|---|---|---|---|---|---|
| P1-C1 | CosMx He2022: TACSTD2-high tumor vs CD8+NK neighbors | YES | CD8+NK count hi/lo 0.455 at 10 µm and 0.670 at 20 µm; 8/8 sections and 5/5 patients; sign P=0.031. Immune-fraction ratios 0.519 / 0.643. CD8-only count ratios 0.502 / 0.726. | Same ratios. This is the short-range TACSTD2 result, not a replacement for locked CLDN4 0.36 / 0.52. | 295,877 tumor cells; 8 sections; 5 patients | #726 (table in #735, #738, #739) |
| P1-C2 | OncoSG LUAD: TACSTD2 vs T/NK-class scores after published PURITY | YES | IMSIG NK partial ρ=−0.421; IMSIG T partial ρ=−0.404 (n=169). | CD8A partial ρ=−0.309 (95% CI −0.440 to −0.165, p=4.69e-5). GEP18 partial ρ=−0.349. A1 immune partial ρ=−0.318. Unadjusted CD8A ρ=−0.380. | n=169 (portal lists 181; public z-score matrix has 169) | #139 (restated #735, #736, #739) |
| P1-C3 | TCGA LUAD+LUSC: TACSTD2 vs CD8A after ABSOLUTE purity | YES | LUSC partial ρ=−0.244 (p=4.3e-8, n=493). | Pooled partial ρ=−0.191 (p=1.31e-9, n=995). LUAD partial ρ=−0.104 (p=0.020, n=502). CYT pooled partial ρ=−0.151. | 995 tumors with ABSOLUTE purity | #107 (restated #739) |
| P1-C4 | TCGA ESTIMATE ImmuneScore after ABSOLUTE purity | PARTIAL | LUSC partial ρ=−0.131 (p=0.003, n=493). GEP18 is also LUSC-only (partial ρ=−0.226). | Pooled ImmuneScore partial ρ=−0.107 (p=7.4e-4, n=995). LUAD partial ρ=+0.072 (p=0.109). | 995 tumors | #107 |
| P1-C5 | TCGA eight primary cohorts: TACSTD2 vs CD8 after keratin | PARTIAL | All-patient KRT5/6 sweep, CD3/CD8 z-mean: pooled ρ=−0.119 (95% CI −0.176 to −0.062, p=4.98e-5, 7/8 negative). | Pre-specified KRT8+KRT18+KRT19, CD8A+CD8B mean: pooled ρ=−0.069 (95% CI −0.135 to −0.002, p=0.043, I²=76%, 6/8). | LUAD 516, LUSC 501, BRCA 1095, CESC 304, KIRC 533, STAD 412, BLCA 406, PAAD 178 | #593 primary; #705 largest all-patient sweep |
| P1-C6 | CPTAC LUAD: TROP2 protein vs xCell CD8 / immune | PARTIAL | LUAD TROP2 protein vs xCell immune ρ=−0.309 (n=110, p=0.0010); vs xCell CD8 ρ=−0.289 (p=0.0022). | After WES purity, CD8 partial ρ=−0.261 (n=108, p=0.0067); immune partial ρ=−0.272. | LUAD n=110; LSCC n=108 does not replicate | #99 (restated #735, #739) |
| P1-C7 | GSE31210 primary LUAD: TACSTD2 vs CD8A / ESTIMATE ImmuneScore | PARTIAL | Unadjusted ImmuneScore ρ=−0.306 (p=2.75e-6); unadjusted CD8A ρ=−0.289 (p=9.76e-6). | After ESTIMATE, partial ρ=−0.023 (ImmuneScore, p=0.735) and −0.107 (CD8A, p=0.108). | n=226 primary tumors | #736 |
| P1-C8 | Public mouse: Tacstd2 % epithelial vs T fraction | PARTIAL | GSE295824 Spearman ρ=−0.612 (n=16, p=0.012). Partial given Cldn4 ρ=−0.468 (p=0.082). | Five-study pool, within-study rank z: total c=−0.277, permutation p=0.137, n=34. Indirect path through Cldn4 includes 0. | 16 mice in the one negative study; pool n=34, k=5 | #724 (study call in #739) |
| P1-C9 | GSE137244 KL vs KP: IFN-compact score | PARTIAL | IFN compact (15 genes) Δ=−0.871, exact MW p=0.0317. Libraries overlap. | Same. APM MHC-I Δ=−0.500, p=0.31 (not significant). STING core is higher in KL (Δ=+0.475), not lower. | 5 KL vs 5 KP libraries | #685 |
| P1-C10 | GSE131907: TACSTD2 %pos is higher in malignant / epithelial cells than in T/NK | YES | MPE carcinoma-like TACSTD2 %pos 83.8 vs T 1.4 / NK 3.1 / B 1.3 / myeloid 5.6. | Tumor-site epithelial TACSTD2 %pos 75.0 vs T 2.8 / NK 3.2 / B 2.6 / myeloid 9.7 (#230). | MPE carcinoma-like cells n=259; tumor-site comparison is cell-compartment %pos, not a patient ρ | #745; tumor-site gap #230 |

| ID | Caveat (do not drop) |
|---|---|
| P1-C1 | At 50 µm the count ratio is 0.854 (7/8); at 100 µm 0.947 (4/8). Inside CLDN4-low the 10 µm immune-fraction ratio stays 0.507 (8/8 and 5/5); rank attenuation after CLDN4 is only 0.124. Not accounted for by the cell's own CLDN4. |
| P1-C2 | East-Asian surgical LUAD, not ICI. Do not quote IMSIG neutrophil partial ρ=−0.456 as CD8. |
| P1-C3 | LUSC is stronger than LUAD (Fisher p=0.023). Single-gene CD8A, not a deconvolution fraction. Not ICI. |
| P1-C4 | LUSC-only. Do not quote the pooled ImmuneScore as a LUAD finding. |
| P1-C5 | The sweep maximum inflates \|ρ\|. CESC partial ρ=+0.118 on the pre-specified keratin control. Bulk association is not spatial exclusion. KRT5/6-low tertile grid max ρ=−0.152 is a subset, not the all-patient row. |
| P1-C6 | LUAD-only. LSCC CD8 ρ=−0.080 (p=0.41). Treatment-naive surgery. RNA TACSTD2 vs xCell is null (\|ρ\|≤0.12). TROP2 vs CD8A protein CIs cross 0 in both histologies (#721). |
| P1-C7 | The negative sign does not survive the ESTIMATE partial. East-Asian surgical LUAD. |
| P1-C8 | Not a KL-vs-KP experiment. GSE264739 is opposite (ρ=+0.829, p=0.058). Do not merge with private 8KL. The n=4 ceiling ρ=−1.000 (exact p=0.083) is exploratory (#743), not this row. |
| P1-C9 | The Tacstd2 split is completely genotype-confounded (every KL library is above every KP library). A 192-cell sweep found 0 cells with Tacstd2 and Cldn4 up together with NHEJ down and IFN up. |
| P1-C10 | This is a dissociated detection gap, not patient-level exclusion. MPE sample Spearman is underpowered (3 PE samples with ≥20 carcinoma-like cells). Tumor-site TACSTD2 vs CD8 ρ=+0.11 (n=36, p=0.54, #230). |

## Part 1 — Tacstd2 → TJ

| ID | Limb | Fit | Largest effect | Locked primary | n | PR |
|---|---|---|---|---|---|---|
| P1-T1 | GSE137244 KL vs KP: Tacstd2, Cldn4, and TJ scores | YES | Cldn4 Δ=+5.570. Cldn4-edge (11 genes) Δ=+2.449. TJ_TISMO (7 genes) Δ=+3.269. Tacstd2 Δ=+3.238. TJ epithelial (18) Δ=+1.627. | Tacstd2 +3.238, Cldn4 +5.570, and every TJ score above except the 107-gene within-arm neighborhood, all with exact MW p=0.00794 and complete KL>KP separation. | 5 vs 5 libraries | #685 (deltas restated #746, #743, #724) |
| P1-T2 | Concordant-4 malignant TACSTD2 Q4 vs Q1: junction and keratin up | YES | ORA: Hallmark apical junction rank 1, enrichment 5.89, FDR 1.3e-5. Keratinization rank 3, enrichment 13.98, same FDR. | GSEA keratinization NES=+2.83 (rank 2 of positive NES, FDR 0.001). KEGG tight junction NES=+2.04 (FDR 0.001) but NES-rank 22. TJ family logFC=+0.223, FDR 0.0042 (206 genes, TACSTD2 held out). | DE contrast 19 vs 15; expression units 64 (not n=65) | #741 |
| P1-T3 | GSE31210: TACSTD2-high GSEA, tight junction and keratin up, EMT down | YES | KEGG tight junction NES=+2.00 (primary FDR 0.00183). Keratinization NES=+1.91. Hallmark EMT NES=−2.15. | Same three calls are the pre-specified A8 pass. Focal CLDN4 MAS5 Δ=+704, FDR=2e-9. | n=226 | #736 |
| P1-T4 | OncoSG: TACSTD2-high keratin up and EMT down; KEGG TJ not up | PARTIAL | Keratinization NES=+1.92 (FDR 0.002). Hallmark EMT NES=−2.33 (FDR 0.002). | KEGG tight junction NES=−1.00 (FDR 0.32). Hallmark apical junction NES=−1.48 (FDR 0.0045), down rather than up. Focal CLDN4 z Δ=+1.28, FDR=2e-6. | n=169 | #736 |
| P1-T5 | TCGA histology-split: TACSTD2-high tight-junction GSEA | YES | LUSC tight-junction assembly NES=2.18, FDR=0.004. | LUAD KEGG tight junction NES=2.09, FDR=0.005. Both are within-histology, not a squamous-mix artifact. | TCGA-LUAD and TCGA-LUSC, TACSTD2-high vs low | #200 |
| P1-T6 | GSE131907 tumor-site epithelium: TACSTD2 prerank tight-junction GSEA | YES | CUSTOM_TJ_CORE NES=2.21, FDR=0. | Still NES=2.17, FDR=0, after removing CLDN4 and TACSTD2. KEGG TJ NES=1.74, FDR=0.0075. | n=36 tumor-site samples | #230 |
| P1-T7 | GSE137244 and GSE165641: is TJ the top Tacstd2-high pathway? | PARTIAL | Focused epithelial TJ-core NES=2.34, rank 1 on GSE137244; GSE165641 KL GEMM NES=2.17, rank 1 (#746). | Broad 589-set universe (#744): no TJ, adhesion, or claudin-family set is in ranks 1–3. GSE137244 best TJ rank is 21 (GOBP TJ organization NES=+1.99, FDR 0.086). | GSE137244 is 5 vs 5 and genotype-confounded | #746 focused rank; #744 broad rank |
| P1-T8 | Public mouse epithelium: Tacstd2-high vs low TJ_TISMO module | YES | 7/7 mice Tacstd2-high > low; mean Δ=+0.261; binomial p=0.0078; all 7 also MW p<0.05. Rank enrichment p=3.6e-9. | Same frozen TJ_TISMO module. | 7 mice | #743 |
| P1-T9 | GSE131907 MPE: TACSTD2 Q4 vs Q1 TJ-core | YES | TJ-core Δ=+0.16, p=9.5e-8. | Holds after dropping CLDN4 (p=4.7e-8). CLDN4 co-detection OR=5.00 (88% vs 60%). Primary tS sample ρ(TACSTD2, CLDN4)=0.77 (n=10, p=0.0092). | MPE carcinoma-like cells (author Malignant label is 0 in PE) | #745 |
| P1-T10 | DepMap / Gygi CCLE lung lines: TROP2 with CLDN4 and TJ scores | YES | RNA NSCLC, CLDN4-edge (11 genes) ρ=+0.696 (n=143, p=4.92e-22). Protein grid maximum ρ=0.846 (NSCLC primary lines, partial on subtype, n=17, p=7.1e-5). | Gygi S1 lung TenPx TACSTD2–CLDN4 protein ρ=0.693 (n=45, p=1.31e-7). RNA CLDN4 ρ=+0.667; TJ epithelial ρ=+0.668; TJ TISMO ρ=+0.663 (all n=143). | RNA n=143 NSCLC lines; protein n=45 lung lines | #740 RNA and n=45 protein; #697 protein grid maximum |
| P1-T11 | TCGA: TACSTD2 tracks CLDN4 after KRT8/18/19 | YES | TACSTD2–CLDN4 keratin-partial pooled ρ=0.315 (19/21 carcinomas, p=3.82e-20). | Same. TACSTD2–CLDN7 pooled ρ=0.242 (19/21, p=4.35e-13). LUAD/LUSC plus the keratin funnel are 8/8 for both. | 21 epithelial carcinomas, primary tumors | #593 |
| P1-T12 | GSE131907: CLDN4-positive cells co-detect TACSTD2 | YES | Primary tS cells, CLDN4+ vs CLDN4−: TACSTD2 detection 92.0% vs 60.7%, OR=7.45. | Sample-mean Spearman ρ=0.770 (n=10 primary tS samples, p=0.0092). MPE carcinoma-like OR=5.00 (88.4% vs 60.5%, p=3.7e-5). | 6,352 primary tS cells; 10 samples; MPE carcinoma-like n=259 | #635 (MPE OR also in #745) |
| P1-T13 | TCGA-LUAD: structural 15-gene TJ vs CD8 / GEP | YES | TJ-15 vs CD8 Spearman ρ=−0.29 (p=1.6e-11); vs GEP ρ=−0.28 (p=9.0e-11). | ESTIMATE-purity partial ρ=−0.26 (p=1.7e-9, n=511). | n=515 tumors (partial n=511) | #90 |
| P1-T14 | CosMx: epithelial program that contains CLDN4 and TACSTD2 vs CD8 at 50 µm | YES | Hotspot-smoothed 10-gene program, outer quartiles, patient-mean Spearman ρ=−0.228; 5/5 patients and 8/8 sections; toroidal-shift p=0.005. High/low CD8-neighbor count ≈0.61×. | Keratin-free adhesion score (CLDN4, CDH1, EPCAM, TACSTD2) ρ=−0.213, still 5/5 and 8/8. Keratin-residualized program ρ=−0.081. | 225 FOVs | #643 (restated #738) |

| ID | Caveat (do not drop) |
|---|---|
| P1-T1 | Handoff TJ +3.03 is a different average. Do not relabel TJ_TISMO +3.269 as +3.03. In vivo Tacstd2, not TJ: GSE6135 Δ=+1.795 (7 vs 5, p=0.00253); GSE164758 Δ=+0.858 (9 vs 8, p=8.2e-5) (#743). |
| P1-T2 | Broad TJ vs T/NK ρ=+0.137 (p=0.60); TJ-core vs T/NK ρ=−0.191 (p=0.41). Do not write that TJ-high excludes T/NK. IFN and MHC are not down on this TACSTD2 split. Locked CLDN4 %pos vs T/NK stays ρ=−0.531. |
| P1-T3 | Hallmark apical junction NES=−1.16 (FDR 0.099), so not every junction set is up. Surgical LUAD, not ICI. |
| P1-T4 | Keratin and EMT match the slide. The KEGG tight-junction set does not. |
| P1-T5 | EMT-down does not replicate as pan-NSCLC (pooled NES +1.30, FDR 0.50). Strict triple intersection (TCGA NSCLC ∩ GSE207422 ∩ GSE131907, r≥0.20, FDR<0.05) is 20 genes and includes CLDN4 and ELF3. |
| P1-T6 | Keratinization NES=1.48 (FDR 0.031) is the weaker companion. Dissociated 10x, treatment-naive, not spatial. |
| P1-T7 | Do not write “TJ is pathway #1” from the broad ranking. Closest broad adhesion call is KEGG CAMs rank 10 in GSE164758 (NES=+2.06). |
| P1-T8 | Public integrate cohorts only. Not merged with private 8KL. |
| P1-T9 | Does not replace locked author-malignant CLDN4–T/NK results. Sample-level TACSTD2 vs T/NK is not supported (3 PE samples). |
| P1-T10 | The 0.846 interval still covers 0.69. Do not quote an RPPA n≈118 (no TROP2/CLDN4 antibodies). Hallmark IFN-γ ρ=+0.339 is positive, so lines are not immune-cold. No T cells and no ICI labels. |
| P1-T11 | COAD and READ are the two non-positive cohorts. COXPRESdb: mouse Tacstd2→Cldn4 is rank 1; human TACSTD2→CLDN4 is rank 3 (#200). |
| P1-T12 | MPE double-positive expression levels are null (ρ=−0.025). ELF3 co-detection OR=22.6 is a different gene. |
| P1-T13 | Broad KEGG tight junction vs CD8 is ρ≈+0.07 (p≈0.10) and is not this limb. This is bulk RNA, not scRNA exclusion. |
| P1-T14 | Continuous CLDN4 alone is ρ=−0.065 (4/5; Lung6 positive). The spec is the most negative 5/5 result in a 2,754-score grid. Only CLDN4 is on the 960-plex. Does not replace locked 0.36 / 0.52. |

## Part 2 — CLDN4 pin

| ID | Limb | Fit | Largest effect | Locked primary | n | PR |
|---|---|---|---|---|---|---|
| P2-1 | Concordant-4: malignant CLDN4 vs patient T/NK fraction | YES | scVI CLDN4-high fraction vs T/NK-among-malignant, within-cohort rank: DL ρ=−0.736 (I²=22%); LR χ²=36.25. Cohorts −0.764 / −0.517 / −0.802 / −0.867. I²=0 neighbor: χ²=41.03, ρ=−0.702 (#703). | %pos (UMI>0) DL ρ=−0.531 (95% CI −0.697 to −0.312, p=1.65e-5, I²=0%, N=65). Stacked Q4 vs Q1 Cliff δ=−0.724 (19/16, p=2.88e-4). Joint maximum of 6,435 count panels on both \|ρ\| and \|Cliff δ\| (#712). I²=0 partial on KRT18 is ρ=−0.533. | N=65 (13+21+22+9). Not cell counts. | #539 / #503 locked; #712 joint max; #703 largest \|ρ\| |
| P2-2 | GSE131907 author-malignant CLDN4 vs patient T/NK | YES | Mean log-normalized CLDN4 ρ=−0.532 (p=0.013, n=21). | %pos ρ=−0.478 (p=0.028). Q4 vs Q1 rank-biserial r=−0.600 (6 vs 5, p=0.12). | n=21 patients | #541 |
| P2-3 | CosMx He2022: CLDN4-high tumor has fewer immune neighbors | YES | Immune-cell fraction 0.229× at 9 µm (CLDN4 at or above the section 75th percentile of positive counts vs count 0); 8/8 sections and 5/5 patients (#696). Prior count≥1 at 10 µm is 0.377×. | Cytotoxic neighbor ratio 0.36 at 50 µm and 0.52 at 100 µm; 8/8 and 5/5; sign P=0.031. Not recomputed. | 8 sections; 5 patients; 765,771 cells in the atlas | Locked ratios in the public handoff and #698; largest 8/8 fold #696 |
| P2-4 | CosMx: exclusion, not muzzling, of nearby effectors | YES | GZMB / PRF1 / NKG7 / IFNG high/low 1.11–1.22; decreased in 0/8 sections. | Same. Squidpy nearest-tumor GZMB CPM is higher in 8/8 at 100 µm (median ratio 1.18, #629). | 8 sections; 5 patients | Public handoff; restated #629 and #738 |
| P2-5 | CosMx: CLDN4-high tumor is less likely to contact a CD8 cell | YES | FOV-stratified Mantel–Haenszel OR 0.560 (0.547–0.572). | Pooled OR 0.584 (0.567–0.602). 5,000 FOV-restricted permutations p=0.00020. OR<1 in 8/8 slides and 5/5 tissues. | Official 960-plex; overall contact rate 8.9% | #567 |
| P2-6 | CosMx: CLDN4-specific cross-type g(r) vs CD8 | YES | Mean Δg_inhomogeneous = −0.043 at r = 22.5 µm (Wilcoxon p=1.1e-12; Stouffer p=2.0e-10). | Same headline. 40/159 FOVs have label-permutation p<0.05; 17 survive BH-FDR. | 159/233 FOVs; 8 samples; 5 patients | #561 |
| P2-7 | Concordant-4 malignant CLDN4-high: IFN and MHC pathways down | YES | Hallmark IFN-γ NES=−3.87 (patient-bootstrap 95% −4.00 to −2.58). No eligible IFN test in a 310-test grid was more negative (#709). | Prior fgsea IFN-γ NES=−3.85 (bootstrap −4.01 to −2.70); IFN-α −3.42; MHC-I −2.75 (#691). Re-fit MHC-I NES=−2.71. Family OLS logFC: chemokine −0.964, IFN −0.584, MHC −0.779 (18 vs 16, #503). | Expression n=64; Q4 vs Q1 = 18 vs 16 | #709 largest NES; #691 fgsea; #503 family OLS |
| P2-8 | CPTAC LSCC: CLDN4 protein vs immune RNA and protein; structural TJ-15 is null | YES | Eligible composite HLA-A + HLA-C + HLA-F + CD8A protein ρ=−0.531 (n=78, CI −0.680 to −0.333; WES partial −0.471; search permutation p=0.001) (#708). | CLDN4 protein vs ImmuneScore ρ=−0.432 (n=78, p=7.9e-5; WES partial −0.298). GEP18 ρ=−0.461. CD8A RNA ρ=−0.437. CD8A protein ρ=−0.444. MHC-I (HLA-A/B/C) ρ=−0.410. TJ-15 protein vs ImmuneScore ρ=−0.082 (n=108, p=0.40). | CLDN4 quantified in 78/108 LSCC tumors | #708 largest composite; #289 ImmuneScore and TJ-15; #596 MHC-I and CD8A protein |
| P2-9 | CLDN4-high malignant signature vs CD8A across bulk LUAD cohorts | YES | Same winning spec, ImmuneScore meta ρ=−0.574 (I²=0%). That is also the ImmuneScore maximum. | CD8A meta ρ=−0.533 (fixed effect, I²=0%, n sum=420, k=4). ssGSEA α=0.75, first 163 genes of the locked 221-gene list, unadjusted. The 221-gene z-mean baseline is ρ=−0.409 (I²=71%). | OncoSG 169; GSE273377 103 and 60; GSE282774 58; GSE233774 30 | #693 |
| P2-10 | Concordant-4: CLDN4-high malignant cells send more barrier ligands to T/NK | YES | CellPhoneDB expression proportion, outer 10% CLDN4: +27.8 percentage points (mean of F11R, NECTIN2, CDH1, LGALS9); n=60; Wilcoxon p=3.69e-11; 4/4 cohorts positive (#713). | Q4 vs Q1 mean +25.6 percentage points (n=64, p=4.91e-12). CellPhoneDB family Δ=+0.253 (expr_prop 0.10, n=53, 52/53 patients positive, #711). LIANA consensus barrier family BH q from 1.3e-9 to 4.9e-8 (#579). | n=60 patients for the decile; n=64 for Q4 vs Q1 | #713 largest proportion; #711 score Δ; #579 consensus q |
| P2-11 | Concordant-4 Milo: T/NK neighbourhoods down where malignant CLDN4 is high | YES | Largest median \|log2FC\| among specs with ≥20 down-hits: 3.102 (24 down, 0 up; k=40, d=30, top 20% vs rest). | Primary Milo model is unchanged. Pareto count maximum: 111 T/NK neighbourhoods down and 4 up out of 708 (median \|log2FC\| 1.596; k=12, d=15, top 40%). | 65 units in every fit | #701 |
| P2-12 | OncoSG: CLDN4 vs CD8A / ImmuneScore after PURITY | YES | Unadjusted ImmuneScore ρ=−0.432 (p=4.59e-9); unadjusted CD8A ρ=−0.416 (p=1.85e-8). | PURITY partial ImmuneScore ρ=−0.308 (p=4.90e-5); partial CD8A ρ=−0.285 (p=1.81e-4). | n=169 | #736 |
| P2-13 | GSE31210: CLDN4 vs CD8A / ImmuneScore | PARTIAL | Unadjusted ImmuneScore ρ=−0.366 (p=1.50e-8); unadjusted CD8A ρ=−0.341 (p=1.45e-7). | ESTIMATE partial ImmuneScore ρ=−0.031 (p=0.639); partial CD8A ρ=−0.127 (p=0.057). | n=226 | #736 |
| P2-14 | TCGA: CLDN4 vs CD8 after keratin adjustment | PARTIAL | All-patient KRT5/6 sweep: CLDN4 pooled ρ=−0.134 (95% CI −0.194 to −0.073, p=1.68e-5, 8/8 negative) (#705). | Pre-specified KRT8/18/19 CD8 mean: CLDN4 pooled ρ=−0.081, beside TACSTD2 −0.069 and CLDN7 −0.109 (#593). Best shared continuous panel (CD3, KRT5/6, no purity): CLDN4 ρ=−0.127 (8/8). | Eight primary cohorts | #593 primary; #705 all-patient maximum |
| P2-15 | Concordant-4: CLDN4–T/NK association remains after TACSTD2 | YES | CLDN4 %pos vs T/NK, partial on TACSTD2: ρ=−0.543 (p=2.0e-5, I²=0%, negative in 4/4). | Unadjusted ρ=−0.531. Two-stage residual after TACSTD2 is r=−0.460 (p=1.2e-4, #725). SEM direct path −0.525 (95% CI −0.726 to −0.264). | N=65 | #733 and #718; residual #725 |
| P2-16 | Public scRNA: CLDN4 vs CLDN7 / EPCAM on the same T/NK endpoint | PARTIAL | After keratin and CLDN7 are both partialled, CLDN4 unique ρ=−0.465 (permutation p=0.0002). CLDN7 given CLDN4 is ρ=−0.184 (p=0.15). EPCAM's unique partial collapses (ρ=−0.013). | Keratin-partial %pos vs T/NK: CLDN4 −0.478 (p=2.7e-4, I²=0), EPCAM −0.368, CLDN7 −0.352 (p=0.10, I²=56%), CLDN1 +0.225. Head-to-head CLDN4 vs CLDN7 Δ=−0.206, permutation p=0.17; vs EPCAM p=0.18. | N=65 concordant-4 units | #644 |
| P2-17 | SEM: graph TACSTD2 → CLDN4 → T/NK on concordant-4 | PARTIAL | Concordant-4 BIC prefers TACSTD2 → CLDN4 → immune by ΔBIC 14.66 versus the reverse order and 17.95 versus independence. Indirect path −0.263 (95% CI −0.439 to −0.123). | Direct CLDN4 path −0.525 (CI −0.726 to −0.264). Total TACSTD2 effect −0.151 (CI −0.394 to +0.136) includes zero. TCGA seven-cohort BIC sum also prefers this graph (ΔBIC 11.60) with pooled indirect −0.022 and a mediation ratio 42.1% whose CI is 19.5 to 114.9. | Concordant-4 n=65; TCGA seven-cohort sum n=3,444 | #725 |

| ID | Caveat (do not drop) |
|---|---|
| P2-1 | The scVI \|ρ\| is a searched maximum (permutation p≤0.0015). Sign-only ρ=−0.552 has I²=67% and is not promoted. Largest single locked cohort is GSE123902 ρ=−0.659. Patient-bootstrap median Δ of T/NK = −0.339 (interval −0.444 to −0.093, permutation p=0.003, #714). Do not add GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526. |
| P2-2 | IFN family logFC is about 0 in this cohort alone. Four PE captures with n_malignant=0 are why the patient ρ is softer than the sample row (ρ=−0.522). |
| P2-3 | 0.229× does not replace 0.36 / 0.52 (different radius, cut, and immune definition). A raw 0.203× failed the absent-arm floor. Patient-bootstrap winner is median corrected log2 ratio −1.307 (equal-patient ratio of means 0.413) at 10 µm, interval entirely below 0 (#714). At 40 µm, median section means are 1.52 vs 2.69 (median Δ=−1.20), lower in 8/8 (#551). Squidpy median-split ratio 0.883 is 7/8 only (#629). |
| P2-4 | Write exclusion, not muzzling. No effector gene is lower in 8/8 sections. |
| P2-5 | A PanCK-high / CD45-low protein index does not show this deficit. F11R and NECTIN2 are absent from the panel. CDH1 is enriched in the contacts that remain (Δ log-norm +0.053). |
| P2-6 | Not detected in Lung5 serial sections. Homogeneous g≪1 versus complete spatial randomness is tumor–stroma geometry and is not a CLDN4-specific claim. |
| P2-7 | GSE131907 alone is flat for IFN. Dropping GSE205335, IFN-γ NES falls to −2.40. KEGG tight junction NES=+1.13 and its bootstrap crosses zero — this ranking is not “TJ pathway up.” Largest signature logFC is −0.971; largest canonical gene is CCL5 −2.451. |
| P2-8 | Treatment-naive surgery, not ICI. Filling the 30 undetected CLDN4 values weakens \|ρ\| (best filled row −0.474, n=108). The composite is a searched sum; ImmuneScore −0.432 is the pre-specified single-protein row. TJ-15 null on the same freeze is the within-TJ pin. |
| P2-9 | Purity or stromal partial meta ρ=−0.416. The single largest cell (GSE282774 ImmuneScore ρ=−0.690, size 22) is not the cross-study spec. Prefix length was the only size knob; genes were not picked on their CD8 correlation. |
| P2-10 | Do not quote CellChat probability +0.0037 (#616) as the effect size. The probability-scale searched maximum is +0.191 on n=19. CXCL9/10/11 are not barrier-sized. This is expression, not spatial exclusion. |
| P2-11 | No specification wins both the count and the \|log2FC\|. The earlier 65-down row (median \|log2FC\| 1.991) is inside the grid and off the front. |
| P2-12 | East-Asian surgical LUAD, not ICI. |
| P2-13 | The partial correlations are not significant. |
| P2-14 | Effects are small, and the sweep inflates them. The surface-gene screen does not pin CLDN4: LUAD rank 7, partial +0.122, with CLDN7 / EPCAM / MUC1 stronger (public handoff). BRCA keratin caveat on the TACSTD2 screen is −0.071. |
| P2-15 | TACSTD2 %pos vs T/NK is itself null (ρ=−0.112, p=0.46). The immune association sits on CLDN4. This is not a claim that TACSTD2 cold is mediated by CLDN4: the concordant-4 and CosMx mediation matrix is null (#718), and no grid row kept a negative direct TACSTD2 coefficient with a mediation proportion inside (0, 1) (#733). |
| P2-16 | Public data do not separate CLDN4 from CLDN7 by a difference test. Inside GSE131907 the CLDN7 keratin partial is stronger (−0.658 vs −0.514). Dropping GSE205335 reverses the meta rank. TCGA keratin CD8 ρ moves from −0.069 to −0.040 after CLDN4 (share 41.3%), but the CLDN4-minus-CLDN7 share interval contains 0, and CLDN4's further share after CLDN7 and EPCAM is 1.8% (#728). Where TACSTD2–T/NK is negative, CLDN7 attenuates it about as much as CLDN4 (#732). The stable pin is the private KD co-culture, which is not in this repository. |
| P2-17 | Do not quote the 174% path ratio as a mediation share. GSE205335 alone prefers independence. TCGA cohorts do not vote as one graph, and adding LUSC flips the sum. A lower BIC is a covariance description, not a knockdown. |

## Do not paste as support

| Item | Call | Number | Source |
|---|---|---|---|
| Concordant-4 TACSTD2 %pos vs T/NK | NULL | ρ=−0.112, p=0.46, I²=15%, N=65. Do not quote CLDN4 ρ=−0.531 as TACSTD2. | #718 / #739 |
| CosMx TACSTD2 at 50 and 100 µm | Not 8/8 | Count ratios 0.854 (7/8) and 0.947 (4/8). Short-range 0.455 / 0.670 is the YES row. | #726 |
| CPTAC TROP2 protein vs CD8A protein | NO | LUAD ρ=−0.038 and LSCC ρ=−0.103; both CIs cross 0. | #721 |
| CPTAC LSCC TROP2 protein vs xCell CD8 | NULL | ρ=−0.080, n=108, p=0.41. | #99 |
| TISMO Tacstd2 after ICB | Not fewer T/NK | 49/64 groups up, Wilcoxon p=5.8e-5. Cldn4 is 34/64, not 49/64. Baseline CD8 link is positive. No true KL/KP lung line. | #152 / #173 |
| Public GEMM Cldn4 %pos vs T/NK | Opposite | REML pooled ρ=+0.477, p=0.302, I²=72%, 34 mice. IFN pooled ρ=+0.092. | #677 |
| Mouse n=4 Cldn4 or Tacstd2 vs T fraction | Exploratory | ρ=−1.000, exact p=0.083. Not confirmatory. IFN/APM on the same grid can be ρ=+1. | #700 / #743 |
| DepMap NSCLC TROP2 vs Hallmark IFN-γ | Opposite of cold | ρ=+0.339 (n=143). TJ does not explain it away (partial +0.277). | #740 |
| Concordant-4 broad TJ score vs T/NK | NULL | ρ=+0.137, p=0.60. TJ-core ρ=−0.191, p=0.41. Immune exclusion stays on CLDN4 %pos. | #741 |
| CLDN4 Q4 vs Q1 KEGG tight junction | Not up | NES=+1.13; patient bootstrap crosses zero. | #691 |
| TACSTD2 → CLDN4 → immune mediation | NULL as a proportion | CosMx and concordant-4 primary mediation rows are null (#718). No spec kept mediation proportion inside (0, 1) with a negative direct coefficient (#733). CosMx 10–20 µm TACSTD2 exclusion is not accounted for by CLDN4 (rank attenuation 0.124 / 0.103, #726). | #718 / #733 / #726 |
| TCGA surface ranking as a CLDN4 pin | Does not pin | LUAD CLDN4 rank 7, partial +0.122. CLDN7, EPCAM, and MUC1 are stronger. Head-to-head vs CLDN7 on concordant-4 is p=0.17 (#644). | Handoff; #644; #728 |
| GSE126044 non-responders have higher TJ | Fragile | Recovered only for one 7-gene z-mean (p=0.019, n=5 vs 11). CLDN4 alone p=0.115. Response is aliased with FFPE. | #90 |
| Visium / Stereo-seq / GeoMx / official Xenium lung panel | Do not upgrade | Mixed spots, missing CLDN4, or opposite sign. Do not write a same-spot correlation as spatial exclusion. | Handoff; #544–#566 except #551 and #561 |
| Human ICI bulk CLDN4 rises after resistance | Not supported | GSE126044, GSE135222, GSE248249, GSE248378 and related sets are small, mixed, or opposite. | Handoff |
| Public mouse merged with private 8KL | Do not merge | GSE165641, GSE180963, GSE154977 stay separate from the private eight KL matrices. | Handoff |
| Discordant human scRNA added to concordant-4 | Do not add | GSE148071, GSE127465, GSE154826, GSE200563, E-MTAB-13526 dilute or flip the sign. | Handoff |
| Broad Tacstd2-high GSEA, TJ in ranks 1–3 | NO | any_cohort_highlight_in_top3 = false on a 589-set universe (#744). The focused-set rank 1 (#746) stays PARTIAL in P1-T7. | #744 |

## What this table does not do

- It does not recompute CosMx ratios 0.36 / 0.52 or concordant-4 ρ=−0.531.
- It does not treat a same-spot Visium correlation as spatial exclusion.
- It does not merge public mouse Harmony objects with the private eight KL matrices.
- It does not claim TISMO 49/64 is immune exclusion. That lock is Tacstd2 up after ICB.
- It does not claim public data uniquely rank CLDN4 over CLDN7. P2-16 is PARTIAL on purpose.
