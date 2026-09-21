# FINDING — concordant-4 CLDN4 quartile GSEA is an observational proxy, not a knockdown

The purification sweep (T/B/NK markers removed, quintiles, continuous CLDN4, malignant-cell AUCell, NHEJ/STING) is in `FINDING_SWEEP.md`. It does not replace the numbers below.

**NOT a true KD.** No CRISPR, no siRNA, no CLDN4-loss culture. The split is the endogenous lowest versus highest malignant CLDN4 %pos quartile in the concordant-4 cohorts (GSE123902, GSE131907, GSE205335, GSE189357). Patient, donor, or sample is the unit. Positive NES means the gene set is higher in the **lowest** CLDN4 quartile.

The private comparator is a direction supplied for this wave, not a statistic recomputed here: **IFN up after CLDN4 loss**. The private knockdown matrix is not in this repository.

Pre-specified match rule, IFN-γ only: on this lowest-versus-highest rank, NES > 0 and nominal p < 0.05 is the **same direction** as that private result. NES < 0 and nominal p < 0.05 is the **opposite direction**. Otherwise the proxy does not resolve the direction. APM and tight junction are reported beside it. They were not given a private direction.

---

## Call

| readout | this proxy | private KD |
|---|---|---|
| IFN-γ | NES +3.789, nominal p 0.001, headline FDR 0.001 | IFN up after loss |
| direction | **same direction** | — |
| IFN-α | NES +3.326, p 0.001, FDR 0.001 | not supplied |
| MHC-I / APM | NES +2.596, p 0.001, FDR 0.001 | not supplied |
| KEGG TJ, CLDN4 removed | NES -1.061, p 0.336, FDR 0.336 — not resolved | not supplied |

The IFN leading edge is not a clean tumor-cell ISG list. The top of the lowest-versus-highest rank is `TRAC, SAMSN1, CD48, SRGN, CTSW, GYPC, BCL2A1, TRBC1`. IFN-γ's leading edge includes `GZMA`, `CD69`, `CD86`, and `LCP2` along with `STAT1`, `TAP1`, and `PSMB9`. That mix fits immune transcripts inside the malignant gate, in line with the locked concordant-4 result that CLDN4-high tumors carry fewer T/NK cells. A tumor-cell knockdown is a different experiment. The IFN sign matches the private KD. This proxy does not identify the mechanism.

Honest GSEA n is **18 lowest vs 16 highest** malignant pseudobulks (not the T/NK N=65). Genes ranked: **15779**. CLDN4 itself, on the Q4-minus-Q1 coefficient this rank negates, is logFC +1.678, p 0.012. That is the split-gene check, not a pathway result.

Same direction means the low-CLDN4 quartile is the IFN-high side, which is the sign the private knockdown has. It does not mean this quartile split replicates the knockdown, estimates a knockdown effect, or replaces it.

---

## Locked design

| Item | Choice |
|---|---|
| What this is | Observational proxy |
| What this is not | A CLDN4 knockdown, knockout, or siRNA |
| Cohorts | GSE123902 + GSE131907 + GSE205335 + GSE189357 |
| Not included | GSE148071, GSE127465, GSE207422, GSE154826, CD45+ or T-only sets |
| Malignant matrix | Existing UMI-sums (marker-malignant in 123902 and 189357; author-malignant in 131907 and 205335) |
| Quartile | Within-cohort malignant CLDN4 %pos, rank then qcut. Q2 and Q3 are not in the binary contrast |
| Model | log2(TMM-CPM+1), OLS, cohort covariates, coefficient CLDN4_Q4 |
| Rank | −t of that coefficient. Positive = higher in Q1 |
| Engine | Weighted KS p=1, 1000 gene-set permutations, seed 42. Same `gsea_core.py` as the GSE68465 prerank |
| Headline sets | Hallmark IFN-γ, Hallmark IFN-α, custom MHC-I/APM (21 genes), KEGG tight junction with CLDN4 removed |
| FDR | BH inside those four sets, per rank |
| GSE189357 | Q4 n=2, so no single-cohort GSEA. TD6 and TD9 stay in the stack |
| P4001 | Not in the GSE205335 malignant UMI-sum |

---

## Headline NES (lowest vs highest)

| set | ES | NES | nom p | FDR | n in rank |
|---|---:|---:|---:|---:|---:|
| IFN-γ | +0.673 | +3.789 | 0.001 | 0.001 | 198 |
| IFN-α | +0.678 | +3.326 | 0.001 | 0.001 | 95 |
| MHC-I / APM | +0.761 | +2.596 | 0.001 | 0.001 | 21 |
| KEGG TJ (CLDN4 out) | -0.254 | -1.061 | 0.336 | 0.336 | 147 |

Leading edge (first genes on the low-CLDN4 side of the rank):

- IFN-γ: `CCL5,GZMA,CSF2RB,GBP4,CD69,TAP1,IL10RA,EPSTI1,IL15RA,STAT1,SAMHD1,FGL2,CD86,JAK2,LCP2,SAMD9L,PSME2,PSMB9,PFKP,STAT4,OAS2,IFI35,LAP3,RTP4,TNFAIP3`
- MHC-I / APM: `TAP2,TAP1,PSMB9,HLA-E,HLA-B,B2M,HLA-F,HLA-C,ERAP1,PSMB10,PSMB8,IRF1,HLA-A,TAPBPL,TAPBP,NLRC5,ERAP2,HLA-G`
- KEGG TJ without CLDN4: `DLG2,NEDD4L,MAP3K1,ROCK1,PCNA,MICALL2,ROCK2,PPP2R2D,STK11,CLDN9,MARVELD3,PPP2CB,CLDN18,SCRIB,PRKAA1,CDC42,CRB3,MAPK9,IGSF5,AMOTL1,MAPK8,JUN,EPB41L4B,CACNA1D,TJP2`

---

## Per cohort

GSE189357 is not tested alone. The other three are thin; they are context for the stack, not separate claims. IFN-γ stays positive in each testable cohort. GSE205335 has the largest IFN-γ NES. Dropping it leaves IFN-γ NES positive (see the leave-one-out table). GSE205335 is also the cohort where KEGG TJ flips positive on this rank (higher TJ in the low-CLDN4 arm). The stacked TJ result does not survive that, and it is not a TJ claim.

| cohort | n Q1 / n Q4 | IFN-γ NES (p) | APM NES (p) | TJ NES (p) |
|---|---|---:|---:|---:|
| GSE123902 | 4 / 3 | +2.520 (0.001) | +2.284 (0.001) | -1.207 (0.093) |
| GSE131907 | 6 / 5 | +1.291 (0.002) | +1.104 (0.086) | -1.197 (0.136) |
| GSE205335 | 5 / 6 | +3.939 (0.001) | +2.836 (0.001) | +1.343 (0.012) |
| GSE189357 | Q4 n<3 | skipped | skipped | skipped |

Leave-one-cohort-out of the stacked rank:

| dropped | IFN-γ NES (p) | APM NES (p) | TJ NES (p) |
|---|---:|---:|---:|
| GSE123902 | +3.712 (0.001) | +2.706 (0.001) | -0.950 (0.561) |
| GSE131907 | +3.811 (0.001) | +2.807 (0.001) | +1.174 (0.040) |
| GSE205335 | +2.413 (0.001) | +1.956 (0.002) | -1.208 (0.110) |
| GSE189357 | +3.627 (0.001) | +2.723 (0.001) | -1.037 (0.321) |

---

## Sensitivities

CLDN4 removed from the ranked list. The primary TJ set already excluded CLDN4; this also stops the split gene from occupying the bottom of the rank.

| set | ES | NES | nom p | FDR | n in rank |
|---|---:|---:|---:|---:|---:|
| IFN-γ | +0.673 | +3.623 | 0.001 | 0.001 | 198 |
| IFN-α | +0.678 | +3.372 | 0.001 | 0.001 | 95 |
| MHC-I / APM | +0.761 | +2.685 | 0.001 | 0.001 | 21 |
| KEGG TJ (CLDN4 out) | -0.254 | -1.060 | 0.313 | 0.313 | 147 |
| KEGG TJ (CLDN4 in) | -0.254 | -1.063 | 0.336 | — | 147 |
| GO TJ organization (CLDN4 out) | -0.354 | -1.333 | 0.058 | — | 69 |

Same primary rank, with the circular KEGG TJ set that still contains CLDN4, and with GO tight-junction organization (CLDN4 removed). GO organization is not in the four-set FDR. On this quartile rank its nominal p is 0.054. That is not a TJ hit.

| set | ES | NES | nom p | FDR | n in rank |
|---|---:|---:|---:|---:|---:|
| IFN-γ | +0.673 | +3.789 | 0.001 | 0.001 | 198 |
| IFN-α | +0.678 | +3.326 | 0.001 | 0.001 | 95 |
| MHC-I / APM | +0.761 | +2.596 | 0.001 | 0.001 | 21 |
| KEGG TJ (CLDN4 out) | -0.254 | -1.061 | 0.336 | 0.336 | 147 |
| KEGG TJ (CLDN4 in) | -0.261 | -1.105 | 0.251 | — | 148 |
| GO TJ organization (CLDN4 out) | -0.354 | -1.340 | 0.054 | — | 69 |

Continuous malignant CLDN4 %pos (n=64 units in the count matrices, cohort covariates). Rank = −t of the %pos z-score, so positive still means higher when CLDN4 is lower. IFN-γ and APM stay positive. KEGG TJ stays unresolved. GO organization reaches nominal p < 0.05 on this continuous rank only; it was not a pre-specified headline set.

| set | ES | NES | nom p | FDR | n in rank |
|---|---:|---:|---:|---:|---:|
| IFN-γ | +0.567 | +2.985 | 0.001 | 0.001 | 198 |
| IFN-α | +0.555 | +2.602 | 0.001 | 0.001 | 95 |
| MHC-I / APM | +0.710 | +2.397 | 0.001 | 0.001 | 21 |
| KEGG TJ (CLDN4 out) | -0.259 | -1.149 | 0.147 | 0.147 | 151 |
| KEGG TJ (CLDN4 in) | -0.278 | -1.225 | 0.092 | — | 152 |
| GO TJ organization (CLDN4 out) | -0.364 | -1.423 | 0.029 | — | 71 |

---

## What this is not

- Not a true CLDN4 knockdown, and not a substitute for the private KD.
- Not the public CLDN4/TACSTD2 KD experiments (those remain a separate result: they do not cleanly open IFN/APM).
- Not a re-audit of the concordant-4 T/NK ρ = −0.53.
- Not evidence that low CLDN4 causes IFN. The quartile is endogenous. Malignant UMI-sums can still carry immune transcripts. Histology and cohort are only partly held by the cohort covariate.
- Not a claim about APM or tight junction in the private KD. Only the IFN direction was supplied.
- Not GSE148071, GSE127465, GSE207422, or GSE154826.

## Files

- `tables/gsea_headline.tsv` — primary four sets
- `tables/gsea_all.tsv` — primary, per cohort, leave-one-out, drop-CLDN4, continuous
- `tables/direction_vs_private_kd.tsv`
- `tables/rank_low_vs_high.tsv` — gene, Q4-vs-Q1 t, lowest-vs-highest stat
- `tables/quartile_membership.tsv`
- `figures/fig_headline_nes.png` — NES bar
- `figures/fig_enrichment_curves.png`
- `figures/fig_per_cohort_nes.png`

Reproduce:

```bash
python3 methods/concordant4_cldn4_kd_match_gsea/analyze.py
```
