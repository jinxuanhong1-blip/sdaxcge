# Evidence table — SKB264 opening, part 1

Locked public numbers only. Three arms: Tacstd2 KL>KP, TACSTD2 immune-negative, DEG → tight junction.

No new tests. GSE137244 is bulk cell-line RNA-seq. Claudin screening is not in this table. TISMO post-ICB Tacstd2 induction is not an immune-negative row.

Machine-readable copy: `evidence_table.tsv`.

## 1. Tacstd2, KL > KP

Scale is the deposit’s own expression scale. Delta is mean(KL) − mean(KP). Test is two-sided exact Mann–Whitney. GSE137244 is log2(FPKM+1), libraries not mice.

| Cohort | Material | n KL vs KP | Tacstd2 Δ | MW p | Role |
|---|---|---:|---:|---:|---|
| GSE137244 | bulk cell-line RNA-seq | 5 vs 5 | **+3.238** (handoff +3.24) | **0.00794** | locked reference; complete separation |
| GSE164758 | primary tumor bulk | 9 vs 8 | +0.858 | 8.227×10⁻⁵ | in vivo, separates |
| GSE6135 | mouse array, all histologies | 7 vs 5 | +1.795 | 0.002525 | separates |
| GSE6135 | adenocarcinoma only | 4 vs 5 | +0.989 | 0.01587 | sensitivity |
| GSE137396 | in vivo nodule | 5 vs 5 | +1.144 | 0.09524 | same sign, not significant |
| GSE244452 | syngeneic bulk | 3 vs 3 | +7.228 | 0.1 | n=3 floor |
| GSE274352 | cultured cells, empty vector | 3 vs 3 | −1.543 | 0.1 | n=3 floor, opposite sign |
| GSE274351 | LCM early adenoma | 5 vs 5 | −0.860 | 0.5476 | opposite sign |

Sources: PR #606, #685, #743. Handoff rounding of the GSE137244 delta is +3.24. The p-value is the 5-vs-5 floor (2/252).

## 2. Immune-negative (TACSTD2 / TROP2 only)

| Cohort | Endpoint | Number | Call |
|---|---|---|---|
| CosMx He 2022 | CD8+NK neighbors, TACSTD2 high vs low | 10 µm ratio 0.455; 20 µm 0.670; 8/8 and 5/5; sign P=0.03125 | YES, short range |
| CosMx, same | immune-neighbor fraction | 10 µm 0.519; 20 µm 0.643; 8/8 and 5/5; sign P=0.03125 | YES |
| concordant-4 | malignant TACSTD2 % detected vs T/NK | ρ=−0.112 (95% CI −0.388 to +0.183); p=0.460; I²=15%; N=65 | NULL |
| TCGA LUAD+LUSC | TACSTD2 vs CD8A ‖ ABSOLUTE | pooled ρ=−0.191 (p=1.31×10⁻⁹, n=995); LUAD −0.104 (n=502); LUSC −0.244 (n=493) | YES |
| TCGA, same | TACSTD2 vs ImmuneScore ‖ ABSOLUTE | pooled ρ=−0.107; LUAD +0.072 (p=0.109); LUSC −0.131 (p=0.003) | PARTIAL, LUSC only |
| TCGA, 8 cohorts | TACSTD2 vs CD8 ‖ KRT8/18/19 | pooled ρ=−0.069; p=0.043; I²=76%; 6/8 negative | WEAK |
| OncoSG LUAD | TACSTD2 vs CD8A ‖ purity | ρ=−0.309 (95% CI −0.440 to −0.165); p=4.69×10⁻⁵; n=169 | YES |
| OncoSG, same | immune / IMSIG T / IMSIG NK ‖ purity | −0.318 / −0.404 / −0.421 | YES |
| CPTAC LUAD protein | TROP2 vs xCell CD8 | ρ=−0.289 (n=110, p=0.0022); WES partial −0.261 (n=108) | YES |
| CPTAC LSCC protein | TROP2 vs xCell CD8 | ρ=−0.080 (n=108, p=0.41) | NULL |
| Public mouse pool | Tacstd2 % epithelium vs T fraction | c=−0.277; perm p=0.137; n=34; k=5 | WEAK / ns |
| GSE295824 alone | Tacstd2 % vs T | ρ=−0.612; n=16; p=0.012 | one study; not the pool |

Sources: PR #739, which copies #726, #718, #107, #593, #139, #99, #724. CosMx 50 µm and 100 µm TACSTD2 ratios are 0.854 (7/8) and 0.947 (4/8); the slide uses 10–20 µm.

## 3. DEG → tight junction

| Cohort | Contrast | Readout | Number | How to say it |
|---|---|---|---|---|
| concordant-4 malignant | TACSTD2 Q4 vs Q1, stacked | ORA, Hallmark apical junction | rank **1**; overlap 13; enrichment **5.89**; p=3.79×10⁻⁷; FDR **1.28×10⁻⁵** | top pathway in this ORA |
| same | same | ORA, Hallmark EMT | rank 2; enrichment 5.74; FDR 1.28×10⁻⁵ | next term |
| same | same | ORA, GO keratinization | rank 3; enrichment 13.98; FDR 1.28×10⁻⁵ | next epithelial term |
| same | same | family score, TJ 206 genes | logFC **+0.223**; FDR 0.004187 | program up |
| same | same | family score, epithelial adhesion | logFC +0.734; FDR 0.004187 | program up |
| same | same | GSEA, Hallmark apical junction | NES +2.555; FDR 0.00146; positive-NES rank **8** | enriched, not NES #1 |
| same | same | GSEA, KEGG tight junction | NES **+2.040**; FDR 0.00146; positive-NES rank **22** | enriched, not NES #1 |
| same | same | GSEA, GO keratinization | NES +2.825; FDR 0.00146; positive-NES rank 2 | near-top NES |
| same | same | GSEA, leading positive term | Hallmark TNFα/NF-κB, NES +3.200 | do not call this TJ |
| same | TJ score vs T/NK | broad TJ, 206 genes | ρ=+0.137; p=0.605; N=64 | not immune exclusion |
| Public mouse epithelium | Tacstd2 count >0 vs =0, within mouse | TJ_TISMO, 7 genes | **7/7** Δ>0; mean Δ **+0.2608**; binomial p **0.007812**; rank p **3.601×10⁻⁹** | module direction |
| GSE31210 bulk LUAD | TACSTD2 Q4 vs Q1 | KEGG tight junction | NES +2.00; FDR 0.0018 | bulk corroboration, not the malignant ORA |
| OncoSG bulk LUAD | TACSTD2 Q4 vs Q1 | KEGG tight junction | NES −1.00; FDR 0.32 | bulk, not enriched |

DE n for concordant-4 is 19 vs 15 (PR #741). Do not quote n=65 as that n. Mouse module: PR #743. Bulk LUAD GSEA: PR #736.

## Not used as a top-pathway number

| Item | Locked number | Why it stays out of the slide sentence |
|---|---|---|
| GSE137244 bulk, Hallmark apical junction on the KL vs KP rank | NES −1.15; FDR 0.31 | not higher in KL (PR #606) |
| GSE137244, epithelial TJ core NES | +2.14 inside a 5-set junction collection | not a genome-wide rank (PR #606) |
| GSE137244 / GEMM bulk, 589-set GSEA | TJ not in ranks 1–3 | PR #744 |
| Handoff “TJ +3.03” | different average from the frozen lists | PR #606 and #685 do not return +3.03; TJ7 on the same matrix is +3.269 and is not relabeled |
