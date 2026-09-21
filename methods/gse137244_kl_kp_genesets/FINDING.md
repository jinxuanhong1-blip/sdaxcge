# FINDING — GSE137244 KL vs KP gene-set scores

**ADDITIVE public cell-line RNA-seq.** Deng et al., *Nat Cancer* 2021 ([PMID 34142094](https://pubmed.ncbi.nlm.nih.gov/34142094/); [GSE137244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137244)). Five KP libraries and five KL libraries from nodule-derived lines, plus one normal lung held out of the test. No private 8-KL matrices. No TISMO LLC standing in for KL.

Score = unweighted mean of log2(FPKM+1) from `GSE137244_counts.fpkm.csv.gz`. Delta = mean(KL) − mean(KP). Positive means higher in KL. Two-sided exact Mann–Whitney treats the 10 tumor libraries as the units (252 equally likely assignments). Benjamini–Hochberg q is computed on the four pre-specified scores only: Cldn4/Tacstd2, NHEJ, STING core, IFN compact.

---

## Library table

| library | GSM | genotype | line name on the GEO title |
|---|---|---|---|
| B6AL10-1-RNA … B6AL10-5-RNA | GSM4073816–820 | KP (`KrasLSL-G12D/+; Trp53`) | all five titles use **B6AL10** |
| KL155mix-control-2-RNA | GSM4073821 | KL (`KrasLSL-G12D/+; Lkb1`) | KL155 |
| KL47-1-untreated-1-RNA | GSM4073822 | KL | KL47 |
| KLC-RNA, KLD-RNA, KLE-RNA | GSM4073823–825 | KL | KLC, KLD, KLE |
| normal-lung-RNA | GSM4073826 | normal lung | held out |

Honest n is **5 vs 5 libraries**. The KP arm is one GEO line prefix (B6AL10), not five independently named KP lines, and not five mice. The rank p below is the library test. That is the test that returns p = 2/252 = 0.00794 when every KL library sits above every KP library.

---

## Scores

Full table: `tables/set_contrasts.tsv`. Per library: `tables/sample_scores.tsv`. Per gene: `tables/gene_deltas.tsv`.

| set | genes used | Δ KL−KP | separation | exact p | q (primary 4) |
|---|---:|---:|---|---:|---:|
| **Cldn4** | 1 | **+5.570** | KL>KP | 0.00794 | locked single gene |
| **Tacstd2** | 1 | **+3.238** | KL>KP | 0.00794 | locked single gene |
| **Cldn4/Tacstd2** | 2 | **+4.404** | KL>KP | 0.00794 | 0.016 |
| TJ (Cldn3/4/6/7, Cdh1, F11r, Ocln) | 7 | **+3.269** | KL>KP | 0.00794 | continuity row |
| **NHEJ** (KEGG 2019 mouse) | 13/13 | **+0.302** | overlap (U=24) | 0.0159 | 0.021 |
| **STING core** | 5/5 | **+0.475** | KL>KP | 0.00794 | 0.016 |
| **IFN compact** | 15/15 | **−0.871** | overlap (U=2) | 0.0317 | 0.032 |
| Hallmark IFN-α | 89/97 | −0.634 | overlap (U=2) | 0.0317 | secondary |
| Hallmark IFN-γ | 185/200 | −0.600 | overlap (U=2) | 0.0317 | secondary |
| KEGG cytosolic DNA-sensing | 60/61 | −0.201 | overlap (U=2) | 0.0317 | secondary |

Cldn4 +5.57 and Tacstd2 +3.24 on this FPKM file match the locked single-gene deltas, and both are complete library separation (p = 0.00794). The two-gene mean is +4.40, also complete separation (rank-biserial +1).

The 7-gene TJ mean used for the TISMO TJ score is **+3.27** on this same scale, also complete separation. The locked TJ summary remains +3.03. This folder records the gene list that produces +3.27. Cldn3, Cldn4, Cldn6, Cldn7, and Ocln are each higher in every KL library; Cdh1 (−0.48) and F11r (−0.10) overlap.

---

## IFN is lower in the KL libraries

IFN compact Δ = **−0.87** (p = 0.032, q = 0.032, rank-biserial −0.84). Four of five KL libraries sit below every KP library. KLD is the KL library that overlaps the KP range (IFN compact 2.02; the other KL libraries are 0.94–1.21).

The same direction is in both Hallmark lists matched by symbol to this matrix (IFN-α −0.63, 89/97 symbols; IFN-γ −0.60, 185/200 symbols). Unmatched human symbols, including HLA genes and OAS1/OASL, are listed in `tables/genes_unmatched.tsv` and were left unmatched.

Members with the largest KL-low deltas: Ifit3 −3.32 (KP>KL, complete), Irf7 −2.04, Ifit1 −1.98, Isg15 −1.76, Ifit2 −1.53 (KP>KL, complete). Dropping Ifit3 leaves the set at −0.70. Ifnb1 FPKM is 0 in every library, including normal lung, so it does not move the contrast. These are tumor-cell-line transcripts.

---

## NHEJ mRNA is slightly higher in KL

KEGG 2019 Mouse non-homologous end-joining (13/13 genes) Δ = **+0.30** (U = 24, p = 0.016, q = 0.021, rank-biserial +0.92). One KP library (B6AL10-3) sits inside the KL range, so the sets overlap. Leave-one-out stays positive for every dropped gene (lowest is +0.25 after dropping Rad50, the largest single-gene delta at +0.96 and the only NHEJ gene with complete KL>KP separation). Xrcc6 is the main gene in the other direction (−0.38).

This score is mRNA abundance of the KEGG NHEJ gene list. It is the transcript pattern on these libraries: a small KL-high shift, spread across the set.

---

## STING core mean is Ikbke

STING core is Mb21d1 (cGAS; the symbol Cgas is absent), Tmem173 (STING; the symbol Sting1 is absent), Tbk1, Ikbke, and Irf3. ISGs are kept in the IFN scores so the two sets are not the same list.

The five-gene mean is **+0.48**, complete separation, q = 0.016. **Ikbke is +2.16 and is 91% of that set delta** (complete KL>KP). Dropping Ikbke leaves **+0.053**.

| STING core gene | Δ KL−KP | separation |
|---|---:|---|
| Mb21d1 (cGAS) | +0.039 | overlap |
| Tmem173 (STING) | **−0.127** | overlap |
| Tbk1 | +0.525 | KL>KP |
| Ikbke | **+2.163** | KL>KP |
| Irf3 | −0.223 | overlap |

Tmem173 itself does not separate. The KEGG cytosolic DNA-sensing pathway (60/61 symbols; CGAS unmatched) is slightly KL-low (Δ = −0.20, p = 0.032), in the same direction as the IFN scores, because that pathway list contains the ISGs and interferon genes.

---

## What to quote

- Quote **Cldn4 +5.57, Tacstd2 +3.24, Cldn4/Tacstd2 +4.40**, library complete separation, p = 0.00794. This is the KL>KP epithelial result on these lines.
- Quote **IFN compact −0.87** (and Hallmark IFN-α/−0.63, IFN-γ/−0.60) as lower interferon-stimulated transcripts in the KL libraries.
- Quote **NHEJ +0.30** as a small KL-high mRNA shift of the KEGG NHEJ list.
- Quote the STING result as **Ikbke +2.16 with Tmem173 −0.13**, and the five-gene mean only together with the leave-one-out (+0.05 without Ikbke).

Figures: `figures/scores_kl_vs_kp.png`, `figures/gene_deltas.png`.
