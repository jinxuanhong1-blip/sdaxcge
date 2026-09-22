# Part 1 paste sheet — Tacstd2 cold and/or TJ

One existing PR, one number. Copied from that PR’s text. Null and opposite headlines are left out. Nothing here was refit.

**Cold** = TACSTD2 / TROP2 tracks fewer T/NK, CD8, or immune score, or a lower neighbor fraction.
**TJ** = TACSTD2 / TROP2 tracks tight junction, CLDN4, keratin, or barrier up, or CLDN4 / TJ tracks fewer immune cells.

Do not read a Visium same-spot correlation as spatial exclusion. Do not read TISMO 49/64 (Tacstd2 up after ICB) as immune-cold. Concordant-4 TACSTD2 vs T/NK (ρ=−0.112, p=0.46) is null and is not in this list. Locked concordant-4 is GSE123902 + GSE131907 + GSE205335 + GSE189357 only. A row that names GSE148071 is that PR’s own merge, not the locked n=65.

**28 cold PRs. 109 TJ PRs. 137 total.**

## Tacstd2 cold

| PR | Number | What |
|---:|---|---|
| #6 | ρ=−0.309 | CPTAC LUAD TACSTD2 protein vs xCell immune |
| #8 | ρ=−0.309 | CPTAC LUAD TACSTD2 protein vs xCell immune |
| #39 | −0.456 | OncoSG LUAD TACSTD2 vs neutrophils, purity partial |
| #48 | −0.18 | TCGA LUAD+LUSC TACSTD2 vs CD8 (19/19 immune features negative) |
| #58 | -0.665 | Post-durvalumab NSCLC TACSTD2 vs immune score |
| #73 | ρ=−0.71 | Post-durvalumab NSCLC TACSTD2 vs CYT |
| #81 | −0.16 | TCGA-LUAD TACSTD2 vs GZMB, ABSOLUTE partial; weak, not total immune |
| #82 | −0.456 | OncoSG LUAD TACSTD2 vs neutrophils, purity partial |
| #89 | −0.187 | TCGA plus OncoSG TACSTD2 vs immune effector |
| #99 | −0.309 | CPTAC LUAD TACSTD2 protein vs xCell immune |
| #107 | −0.244 | TCGA-LUSC TACSTD2 vs CD8 partial Spearman |
| #112 | −0.665 | Post-durvalumab TACSTD2 vs ESTIMATE ImmuneScore |
| #139 | −0.309 | OncoSG LUAD TACSTD2 vs CD8A |
| #143 | −0.244 | TCGA-LUSC TACSTD2 vs CD8A |
| #157 | −0.166 | TCGA-LUSC TACSTD2 vs CD8 after ESTIMATE (ImmuneScore is opposite) |
| #235 | −0.284 | CPTAC LUAD TACSTD2 vs CD8 after ESTIMATE |
| #238 | ρ=−0.781 | GSE248378 bulk TACSTD2 versus T/NK |
| #242 | r = −0.229 | TCGA NSCLC TACSTD2 versus xCell CD8 |
| #243 | ρ = −0.577 | GSE248378 TACSTD2 versus CD8A residual |
| #256 | −0.83 | GSE253013 malignant TACSTD2 %pos vs T/NK (n=9) |
| #263 | −0.33 | Naive NSCLC malignant TACSTD2 versus T/NK |
| #292 | ρ=−0.262 | GSE218989 TACSTD2 versus CD8A |
| #705 | −0.119 | TCGA TACSTD2 versus CD8 after KRT5/6 |
| #719 | −0.217 | OncoSG TACSTD2 versus CD8A given CLDN4 |
| #726 | 0.519 | CosMx TACSTD2-high immune-neighbor fraction |
| #735 | 0.455 / 0.670 | CosMx TACSTD2 CD8+NK neighbors at 10/20µm |
| #738 | 0.455 / 0.670 | CosMx TACSTD2 CD8+NK neighbors at 10/20µm |
| #739 | 0.455 / 0.670 | CosMx TACSTD2 CD8+NK neighbors at 10/20µm |

## TJ

| PR | Number | What |
|---:|---|---|
| #23 | ρ=−0.432 | CPTAC LSCC CLDN4 protein vs ImmuneScore |
| #54 | 0.693 | Gygi CCLE lung TACSTD2 vs CLDN4 protein |
| #69 | −0.291 | TCGA-LUAD tight-junction score vs CD8 |
| #71 | 0.709 | TCGA-PAAD TACSTD2 vs CLDN4 co-expression |
| #85 | ρ=0.693 | Gygi lung protein TACSTD2 vs CLDN4 |
| #90 | −0.29 | TCGA-LUAD 15-gene TJ vs CD8 |
| #101 | 0.714 | CCLE RNA TACSTD2 vs CLDN4 (rank 2/24, not #1) |
| #111 | 0.693 | Gygi lung protein TACSTD2 vs CLDN4 |
| #117 | +3.18 | TCGA-LUSC TACSTD2-high keratinization GSEA |
| #126 | ρ = −0.25 | Visium lung CLDN4 vs immune score, same-spot (not micron exclusion) |
| #200 | NES 2.18 | LUSC TROP2-high tight-junction assembly GSEA |
| #230 | NES=2.21 | GSE131907 TROP2-high epithelium TJ GSEA |
| #241 | ρ=−0.31 | TCGA-LUAD TJ 7-gene versus CD8 |
| #245 | −0.296 | CPTAC LUAD TJ-15 protein versus ImmuneScore |
| #264 | ρ=+0.833 | GSE207422 TACSTD2 versus keratin score |
| #277 | ρ=0.854 | GSE131907 LUAD TACSTD2 versus barrier/keratin |
| #289 | −0.432 | CPTAC-LSCC CLDN4 protein versus ImmuneScore |
| #296 | ρ=−0.57 | GSE283829 CLDN4 versus ImmuneScore |
| #299 | +0.471 | GSE31210 CLDN4 vs TACSTD2 (CD8 partial is NS) |
| #311 | δ_rb=−0.48 | GSE31210 CLDN4 Q4 versus Q1 CD8A |
| #313 | +0.366 | GSE4573 LUSC CLDN4 vs TACSTD2 |
| #315 | Δ=−2189 | CPTAC LSCC CLDN4 protein Q4 versus ImmuneScore |
| #317 | δ=−0.620 | CPTAC LSCC CLDN4 Q4 versus ImmuneScore |
| #319 | ρ ≈ −0.43 | CPTAC-LSCC protein CLDN4 versus CD8 |
| #329 | −0.57 | GSE182328 CLDN4 Q4 versus Q1 CD8A |
| #333 | −0.416 | OncoSG LUAD CLDN4 versus CD8A |
| #336 | −0.578 | GSE248378 CLDN4 versus T/NK score |
| #352 | −0.20 | TCGA-LUAD CLDN4 Q4 versus Q1 CD8A |
| #364 | −0.183 | GSE8894 NSCLC CLDN4 versus CD8A partial |
| #371 | −0.328 | GSE72094 CLDN4 Q4 versus Q1 CD8A |
| #372 | −0.219 | GSE42127 LUAD CLDN4 versus CD8A |
| #378 | −0.292 | GSE37745 LUSC CLDN4 versus CD8A partial |
| #382 | −0.403 | GSE11969 LUAD CLDN4 versus CD8A |
| #391 | −0.425 | GSE179351 baseline CLDN4 vs T/NK |
| #405 | −0.481 | GSE10072 LUAD CLDN4 vs CD8A |
| #423 | −0.484 | GSE253564+GSE148071 CLDN4 vs T/NK |
| #448 | ρ=−0.530 | GSE131907+GSE205335 CLDN4 vs T/NK |
| #451 | −0.420 | Three-cohort CLDN4 percent-positive vs T/NK |
| #452 | ρ=−0.483 | Triple scRNA CLDN4 vs T/NK |
| #453 | −0.412 | LUAD-only scRNA CLDN4 vs T/NK |
| #458 | −0.524 | GSE207422+GSE205335 CLDN4 vs T/NK |
| #466 | −2.09 | Winning-pair CLDN4 Q4 vs T/NK ALR |
| #467 | ρ=−0.362 | GSE131907+GSE205335 CLDN4 vs T/NK |
| #503 | −0.531 | Concordant-four CLDN4 vs T/NK |
| #506 | +0.099 | Concordant-4 CLDN4-high barrier-family ΔP |
| #525 | ρ=−0.452 | Seurat winning-pair CLDN4 vs T/NK |
| #536 | −0.435 | GSE123902+GSE205335 CLDN4 vs T/NK |
| #538 | ρ=−0.461 | Triple Seurat CLDN4 vs T/NK |
| #539 | ρ = −0.531 | Concordant-four Seurat CLDN4 vs T/NK |
| #541 | ρ=−0.532 | GSE131907 patient CLDN4 vs T/NK |
| #543 | ρ = 0.63 | TISMO vitro lines Tacstd2 vs TJ-core |
| #551 | −1.20 | CosMx NSCLC CLDN4-high vs CD8+NK counts |
| #561 | −0.043 | CosMx CLDN4–CD8 Δg at 22.5 µm |
| #567 | 0.584 | CosMx CLDN4-high CD8 contact odds ratio |
| #569 | 0.763 | CosMx CLDN4 cytotoxic-neighbor ratio in tumor core |
| #575 | 0.693 | DepMap lung protein TROP2 versus CLDN4 |
| #579 | 1.3e-9 | Concordant-4 CLDN4-high barrier LR, smallest BH q |
| #582 | −2.066 | GSE334497 Trop2 knockout keratinization down |
| #593 | 0.315 | TCGA TACSTD2 versus CLDN4 after keratin |
| #596 | ρ=−0.444 | CPTAC LSCC CLDN4 versus CD8A protein |
| #606 | NES=+2.14 | GSE137244 KL versus KP tight-junction NES |
| #609 | 0.55 | CosMx CLDN4-high versus low immune fraction |
| #612 | +3.269 | GSE137244 KL versus KP seven-gene TJ |
| #617 | +5.570 | GSE137244 KL vs KP Cldn4 |
| #621 | ρ=−0.73 | GSE325414 LUAD CLDN4 percent-positive versus T/NK |
| #629 | 0.883 | CosMx Squidpy CLDN4-high/low CD8+NK ratio at 50 µm |
| #630 | −1.92 | SKOV3 TACSTD2 shRNA CLDN4 co-drop |
| #635 | ρ=0.770 | GSE131907 primary CLDN4 versus TACSTD2 |
| #636 | +0.164 | CosMx CLDN4 ΔR² for immune neighborhood beyond keratin |
| #638 | 0.746 | Concordant-4 malignant TACSTD2 versus CLDN4 |
| #643 | −0.228 | CosMx epithelial program versus CD8 neighbors |
| #644 | −0.478 | Concordant-4 CLDN4 versus T/NK after keratin |
| #666 | −0.516 | Concordant-4 CLDN4 versus T-cell fraction |
| #668 | −0.175 | CosMx CLDN4 field versus CD8+NK density |
| #669 | −0.976 | CosMx donor-mean CLDN4 versus CD8+NK count |
| #673 | 0.377 | CosMx CLDN4-positive immune fraction at 10µm |
| #682 | −0.531 | Concordant-4 CLDN4 %pos vs T/NK |
| #685 | +3.269 | GSE137244 KL versus KP TISMO TJ score |
| #691 | −3.85 | Concordant-4 CLDN4-high Hallmark IFN-γ NES |
| #692 | +0.259 | Concordant-4 CLDN4-high core-barrier sender Δ |
| #693 | −0.533 | CLDN4-high signature versus LUAD CD8A |
| #694 | +5.570 | GSE137244 KL versus KP Cldn4 |
| #696 | 0.229 | CosMx CLDN4-high immune fraction at 9µm |
| #697 | 0.846 | Gygi primary NSCLC TROP2 versus CLDN4 protein |
| #698 | 0.573 | CosMx CLDN4 cytotoxic-neighbor ratio, 7/8 sensitivity |
| #701 | −0.531 | Concordant-4 Milo keeps CLDN4 %pos vs T/NK |
| #703 | −0.736 | Concordant-4 scVI CLDN4 versus T/NK |
| #706 | +0.245 | Concordant-4 CLDN4-high barrier LR-gated Δ |
| #708 | −0.531 | CPTAC LSCC CLDN4 versus HLA plus CD8A |
| #711 | +0.2533 | Concordant-4 CellPhoneDB barrier-family Δ |
| #712 | ρ=−0.531 | Concordant-4 CLDN4 percent-positive versus T/NK |
| #714 | −1.307 | CosMx CLDN4 immune-fraction log2 ratio |
| #715 | −0.448 | Concordant-4 four-gene module versus T/NK |
| #716 | +25.19 | Concordant-4 TACSTD2-high barrier ligand proportion |
| #724 | 0.502 | Public mouse Tacstd2 percent versus Cldn4 |
| #725 | −0.525 | Concordant-4 CLDN4 direct path to T/NK |
| #731 | −1.915 | SKOV3 TACSTD2 shRNA CLDN4 log2FC |
| #732 | −0.543 | Concordant-4 CLDN4 vs T/NK after TACSTD2 |
| #733 | ρ=−0.531 | Concordant-4 CLDN4 percent-positive versus T/NK |
| #736 | NES=+2.00 | GSE31210 TACSTD2-high KEGG tight junction |
| #737 | +3.269 | GSE137244 KL versus KP TJ_TISMO score |
| #740 | 0.693 | DepMap lung TROP2 versus CLDN4 protein |
| #741 | NES=+2.83 | Concordant-4 TACSTD2-high keratinization NES |
| #742 | ρ=−0.531 | Concordant-4 CLDN4 percent-positive versus T/NK |
| #743 | +0.261 | Public-mouse Tacstd2-high TJ module mean Δ (7/7) |
| #745 | +0.16 | GSE131907 MPE TACSTD2-high TJ-core |
| #746 | 2.34 | GSE137244 epithelial TJ-core NES |
| #747 | −0.531 | Allowed opening: concordant-4 CLDN4 %pos vs T/NK |
| #748 | −0.478 | Concordant-4 keratin-partial CLDN4 versus T/NK |

