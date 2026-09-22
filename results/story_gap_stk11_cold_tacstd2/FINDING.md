# FINDING — Disease context story gap (public only)

## Question

Can a single disease-context slide honestly say **KL/STK11 cold tumor + Tacstd2/Cldn4 high** using only public evidence, without fabricating human antigen elevation?

## Answer (one paragraph)

**Yes for the cold limb. Yes for Tacstd2/Cldn4-high only if the antigen limb is annotated as mouse KL cell-line bulk (GSE137244 lock).** Human STK11 public bulk (TCGA-LUAD, PR #49) shows the **opposite** antigen direction (lower TACSTD2/CLDN4). Public human STK11 scRNA (GSE280232) cannot power a genotype→TACSTD2 test. Mouse scRNA supports epithelial Tacstd2/Cldn4 biology with small-n caveats and must not be blurred with bulk locks.

## Locked antigen numbers (do not recompute)

| Gene | KL − KP (mean log2 FPKM+1) | Exact MW (5 vs 5) | Source |
|---|---:|---:|---|
| Tacstd2 | **+3.238** | **0.00794** | Handoff; PR #606 |
| Cldn4 | **+5.570** | **0.00794** | Handoff; PR #606 |
| TJ (handoff) | **+3.03** | 0.00794 | Handoff |

GSE137244 = Deng et al. *Nat Cancer* 2021, five Kras;Lkb1 vs five Kras;Trp53 cell-line libraries. Not a tumor immune microenvironment assay.

## Cold limb (human)

| Evidence | Number | Class |
|---|---|---|
| Skoulidis 2018 SU2C ORR | KL **7.4%** / KP **35.7%** / K-only **28.6%** | Literature |
| Skoulidis PROSPECT IHC | CD3 P=**0.0019**; CD8 P=**0.0072** | Literature |
| Public MSK KRAS-restricted DCB (PR #49) | Hellmann OR **0.24**; Rizvi2018 OR **0.38**; writeup pooled MH **0.36** | Public reanalysis |

## Human STK11 antigen — honest reverse

TCGA-LUAD (PR #49):

- STK11mut vs WT: TACSTD2 log2FC **−0.51** (FDR **1.0×10⁻⁴**); CLDN4 **−0.47** (FDR **~3×10⁻⁵**)
- Within KRAS-mut: STK11 co-mut TACSTD2 **−0.71** (p=1.0×10⁻⁴); CLDN4 **−0.68** (p=5.7×10⁻⁵)

**Do not put “human STK11 = Tacstd2/Cldn4 high” on the slide.** Put the reverse as a caveat strip.

## Bulk vs scRNA (required honesty)

1. **GSE137244** = mouse **cell-line bulk** → antigen lock.  
2. **GSE179500** = mouse **tumor bulk** → modest Tacstd2↑ when LKB1 off.  
3. **GSE179502 / GSE180963** = mouse **scRNA** → epithelial expression / restore direction; n small; GSE180963 does **not** show cold as lower immune fraction.  
4. **TCGA** = human **bulk** → STK11 cold is literature; antigen is **lower**.  
5. **GSE280232** = human **scRNA** → TACSTD2 epithelial; no powered STK11mut>wt.

## What fills the story gap

The gap was “disease context slide needs both cold and antigen-high.” Public data fill it by **splitting species/assay**:

- Cold = human STK11/KL clinical + IHC literature (+ public MSK direction).  
- Antigen-high = **mouse GSE137244 lock**, not human bulk.  
- Human bulk antigen = **honest negative** for the high claim.

## Files

- `PPT_EVIDENCE_TABLE.md` — paste into slide  
- `tables/evidence_table.tsv` — machine table  
- `provenance.json` — PR citations  
- `figures/honesty_matrix.png` — bulk/scRNA × cold/antigen grid  
- `RESULTS.md` — short PR summary  

No private 8KL matrices. No TISMO LLC-as-KL. No new GEO downloads for this PR.
