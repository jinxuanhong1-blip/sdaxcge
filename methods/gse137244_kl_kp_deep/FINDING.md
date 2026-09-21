# FINDING — GSE137244 KL vs KP, beyond mean Tacstd2 / Cldn4

Public cell-line RNA-seq only. Deng et al., Cancer Discovery 2021 ([PMID 34142094](https://pubmed.ncbi.nlm.nih.gov/34142094/)), [GSE137244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137244). No private 8-KL matrices. No TISMO LLC. Normal lung (`normal-lung-RNA`) is excluded.

This note does not replace the locked single-gene means. It adds gene-set structure, per-library scores, and the exact Mann–Whitney floor those means sit on.

## Design

Ten polyA libraries from nodule-derived lines. KP = KrasG12D; Trp53 (B6AL10-1..5). KL = KrasG12D; Lkb1 (KL155, KL47, KLC, KLD, KLE). Names that say control or untreated are the baseline libraries in this series.

Primary scale is **mean log2(FPKM+1)**, KL minus KP. That is the scale of the locked Tacstd2 and Cldn4 deltas. Inference for a 5-vs-5 score is the two-sided exact Mann–Whitney test. Preranked GSEA permutes genes, not libraries, so its p-values are not bound by that floor.

Genotype check on the same scale: Stk11 is lower in every KL library than every KP library (Δ = −2.523, exact p = 0.00794). Trp53 is lower in every KP library (Δ = +2.433, same floor). The arm labels match the genotypes.

## Locked genes reproduce

| gene | Δ log2(FPKM+1) | exact MW p | separation |
|---|---:|---:|---|
| Tacstd2 | **+3.238** | 0.00794 | all 5 KL > all 5 KP |
| Cldn4 | **+5.570** | 0.00794 | all 5 KL > all 5 KP |

Handoff rounding is Tacstd2 +3.24 and Cldn4 +5.57. Both p-values are the 5-vs-5 floor (below), not a finer tail.

## Tight junction is a claudin leading edge, not one mean

A single “TJ +3.03” is not what the pre-specified sets return on this scale.

| set | genes in matrix | mean Δ | exact p | ssGSEA direction |
|---|---:|---:|---:|---|
| Epithelial TJ core (18 structural genes) | 18 | **+1.627** | 0.00794, complete | KL higher, complete |
| GO bicellular tight junction | 76 | +0.430 | 0.00794, complete | KL not higher (p = 0.15) |
| KEGG tight junction | 167 | +0.219 | 0.00794, complete | KL not higher (p = 0.15) |
| Hallmark apical junction | 196 | −0.116 | 0.15 | KL lower, p = 0.0317 |
| GO gap junction | 32 | −0.277 | 0.15 | KL lower, p = 0.0317 |
| GO adherens junction | 150 | +0.123 | 0.00794, complete | flat (p = 0.55) |

The broad means completely separate because a few high-abundance claudins pull every KL library above every KP library. The shift itself is small (+0.22 to +0.43). Rank-based ssGSEA, which is not dominated by those few genes, does **not** place the broad GO or KEGG lists higher in KL. Hallmark apical junction is a cytoskeletal mix and is not KL-high.

Preranked GSEA (Welch t on log2(FPKM+1); 15,254 genes with FPKM ≥ 1 in at least 2 libraries; 1,000 permutations) does recover a KL-high structural program:

| set | null | NES (Welch) | FDR | NES (log2FC rank) |
|---|---|---:|---:|---:|
| Epithelial TJ core | 5 junction sets | **+2.14** | 0 | +2.21 |
| KEGG tight junction | 5 junction sets | **+1.76** | 0.011 | +1.60 |
| GO bicellular TJ | 5 junction sets | **+1.72** | 0.014 | +1.86 |
| GO adherens | 5 junction sets | +1.59 | 0.026 | +1.39 (FDR 0.060) |
| GO gap junction | 5 junction sets | +0.92 | 0.60 | −1.12 |
| Hallmark apical junction | 50 hallmarks | −1.15 | 0.31 | −1.09 |
| Hallmark EMT | 50 hallmarks | **−1.92** | 0.012 | −2.13 |
| Hallmark IFN-γ | 50 hallmarks | **−1.65** | 0.029 | −1.91 |
| Hallmark IFN-α | 50 hallmarks | −1.41 | 0.094 | −1.82 |

FDR = 0 means no null NES in the 1,000 permutations exceeded the observed one. Junction FDRs are inside a 5-set collection; hallmark FDRs are inside all 50 hallmarks. Gap junction, tested in the same 5-set collection, is null. The epithelial-core leading edge (11/18) is Cldn4, Cldn1, Cldn7, Cgnl1, Marveld2, Marveld3, Tjp1, Tjp2, Ildr1, Cldn3, Ocln.

Volcano callouts on the same log2 scale (`figures/volcano_tj_callouts.png`). Positive = higher in KL. p = 0.00794 is the ceiling, so these genes are not ranked by p.

| higher in KL | Δ | p | lower in KL | Δ | p |
|---|---:|---:|---|---:|---:|
| Cldn6 | +7.99 | 0.00794 | Cldn2 | −6.62 | 0.00794 |
| Cldn4 | +5.57 | 0.00794 | Cldn18 | −5.20 | 0.00794 |
| Cldn7 | +4.28 | 0.00794 | Cldn5 | −4.86 | 0.00794 |
| Cldn3 | +4.16 | 0.00794 | Esam | −4.21 | 0.00794 |
| Cldn1 | +3.57 | 0.00794 | | | |
| Tacstd2 | +3.24 | 0.00794 | | | |
| Marveld2 | +1.60 | 0.00794 | | | |
| Ocln | +1.47 | 0.00794 | | | |
| Tjp1 | +1.35 | 0.00794 | | | |

Cldn6 is the largest claudin shift. It is on the volcano and in the GO/KEGG leading edges. It is not in the 18-gene epithelial core (that core is the simple-epithelium set: Cldn1/3/4/7 plus occludin, MARVEL, ZO, JAM, cingulin, crumbs, angulins). F11r (JAM-A) is abundant and flat (Δ = −0.10, p = 0.69). Cldn5 and Cldn18 going down is the endothelial/alveolar claudin pattern leaving these lines, not a loss of Cldn4.

## EMT is lower in KL

Hallmark EMT ssGSEA is completely lower in the five KL libraries (exact p = 0.00794). The mean log2 score is Δ = −0.843, also complete. Preranked NES = −1.92 (FDR 0.012 within 50 hallmarks). The leading edge is matrix and mesenchymal genes (Bgn, Vcam1, Vim, collagens, Postn, Fbln, Mmp2), not a junction program. On this panel, KL lines are the more epithelial, less EMT-like arm.

## KL lines are not IFN-cold

The question is cell-intrinsic ISG. These cultures have no immune infiltrate, so a low score is not an immune-excluded tumor.

Rule used here: IFN-cold only if **both** hallmark IFN-α and IFN-γ are completely lower in every KL library than every KP library.

That rule is not met.

- Mean log2, IFN-α and IFN-γ: both Δ = −0.66, exact p = 0.0317 (U = 2, two pairwise inversions). Not the floor.
- Compact ISG (15 genes): Δ = −0.871, same p = 0.0317.
- ssGSEA IFN-γ: completely lower, p = 0.00794, but the gap is thin (KL max 15,887 vs KP min 16,030).
- ssGSEA IFN-α: p = 0.0317, not complete.

The inversions are the same two libraries. KP **B6AL10-3** is IFN-low (IFN-α mean 2.37, against 3.31–3.40 in the other four KP libraries). KL **KLD** is the least low KL line (IFN-α mean 3.11) and sits above B6AL10-3. KL47 also clears B6AL10-3 by a small margin. The other four KP libraries sit above every KL library.

So the 5-vs-5 mean leans IFN-lower, one discrete step above the floor, and it is not a cold wall. After near-replicate libraries are averaged (next section), the IFN mean contrast is not separated (cluster p = 0.53 for both hallmark IFN sets and for the compact ISG).

## Soft power: the exact MW floor

For 5 vs 5 with no ties there are C(10,5) = 252 equally likely rank assignments. The two-sided exact p-values are a ladder, not a continuous tail (`tables/mw_exact_ladder.tsv`):

| KL wins out of 25 | Cliff’s δ | exact p | α = 0.05 |
|---:|---:|---:|---|
| 25 or 0 | ±1.00 | **2/252 = 0.00794** | yes, the floor |
| 24 or 1 | ±0.92 | 4/252 = 0.0159 | yes |
| 23 or 2 | ±0.84 | 8/252 = 0.0317 | yes, last step below 0.05 |
| 22 or 3 | ±0.76 | 14/252 = 0.0556 | no |

Nothing can be reported as more significant than 0.00794. A gene with Δ = +5.57 (Cldn4) and a gene with a tiny complete separation share that p. Rejecting at 0.05 requires Cliff’s δ at least 0.84 (at most 2 of 25 pairwise orderings inverted).

Among 16,218 genes with max FPKM ≥ 1, **4,544** hit that floor (28%). Benjamini–Hochberg at FDR 0.05 needs only 2,575 genes at p = 0.00794 before the floor itself passes, so the floor genes are BH-significant and 5,778 genes pass overall. That list does not pick out Tacstd2 or Cldn4. They are tied with thousands of other complete separations. Effect size and the pre-specified GSEA are what rank them. The volcano draws the floor as a ceiling for that reason.

### The 5-vs-5 floor assumes independent libraries

Pearson r of log2(FPKM+1) on genes with mean FPKM ≥ 1 (`figures/library_correlation.png`):

- B6AL10-1 vs B6AL10-2 / -4 / -5 = 0.989 / 0.988 / 0.986
- B6AL10-3 vs B6AL10-1 = 0.822
- KLC vs KLE = 0.986
- KL155 vs KL47 = 0.940

Single-linkage at r ≥ 0.98 collapses the ten libraries to **4 KL clusters vs 2 KP clusters** (the quadruplicate B6AL10 block, B6AL10-3 alone, and KLC+KLE). The exact two-sided floor for 4 vs 2 is 2/C(6,2) = 2/15 = **0.133**. No contrast at that coarser n can reach 0.05, one-sided floor included (1/15 = 0.067).

At that coarser n, Tacstd2, Cldn4, the epithelial TJ mean, and the EMT mean are still completely separated (every KL cluster vs both KP clusters). Their exact p is 0.133, which is the new floor, not a reversal. IFN means are not separated (p = 0.53). The locked p = 0.00794 is the library-level floor. It is not an independent-line p-value.

## What this adds, and what it does not

KL libraries are higher for Tacstd2 and for the epithelial claudins Cldn1/3/4/7, lower for EMT, and only lean lower for IFN. Broad TJ averages look separated because those claudins dominate the mean; ssGSEA of the full GO/KEGG lists does not. The 5-vs-5 p-value cannot go below 0.00794, and the KP arm is mostly one correlated block.

Do not quote 0.00794 as genome-wide discovery. Do not read cell-line ISG as immune exclusion. Do not merge this series with private 8-KL or with TISMO.
