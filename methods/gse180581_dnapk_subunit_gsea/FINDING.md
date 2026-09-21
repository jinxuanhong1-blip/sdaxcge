# GSE180581: Ku70 knockdown versus parental 293T raises IFN signaling

The contrast that shows interferon and DNA-sensing genes up is the full knockdown versus parental 293T + siControl, ranked by mean log2 fold-change. DNA-PK subunit genes (XRCC5, XRCC6, PRKDC) are removed from every set. Positive NES means the set is higher in the knockdown. BH-FDR is within the 12 IFN / STING / DNA-sensing sets of that contrast.

On-target RNA still drops in the matched siRNA contrast: PRKDC log2FC −1.68, XRCC6 −2.01, XRCC5 −2.02.

## Ku70 versus parental siControl

Ku70+/- + siKu70 (n=3) versus parental 293T + siControl (n=3).

| Set | NES | nominal p | FDR | mean log2FC |
|---|---:|---:|---:|---:|
| Reactome IFN signaling | **+1.45** | 0.004 | **0.048** | +0.186 |
| Reactome IFN-α/β | +1.48 | 0.017 | 0.096 | +0.173 |
| Reactome antiviral ISGs | +1.41 | 0.024 | 0.096 | +0.180 |
| KEGG cytosolic DNA-sensing | +1.31 | 0.069 | 0.21 | +0.228 |
| Reactome STING, subunits removed | +0.67 | 0.53 | 0.53 | +0.166 |
| Hallmark IFN-α | +0.98 | 0.29 | 0.31 | +0.048 |
| Hallmark IFN-γ | +1.02 | 0.23 | 0.31 | +0.048 |

Reactome interferon signaling is the FDR < 0.05 call (145 genes in the rank). The leading edge includes BST2, ISG20, IRF5, IRF7, RIGI, IFITM1, IFI35, HLA-C, and HLA-E.

The same Ku70 samples sit above parental on the mean log2 expression of the 32 KEGG cytosolic DNA-sensing genes: knockdown scores 7.392, 7.429, 7.408 versus parental 7.230, 7.184, 7.130 (Δ = +0.228). One-sided Welch p for an increase is 0.0044 (two-sided p = 0.0089). Among the 12 Ku70 score tests, BH-FDR of that one-sided p is 0.026. The DNA-sensing leading edge includes IRF7, STING1, IRF3, MAVS, and NFKB1. The Reactome STING set itself is up (mean log2FC +0.17) and does not pass FDR (score one-sided p = 0.075).

## DNA-PKcs and Ku80, same parental contrast

DNA-PKcs knockdown versus parental: KEGG cytosolic DNA-sensing mean log2FC +0.109, GSEA NES +1.15, nominal p = 0.12, FDR 0.38. The sample-score one-sided p is 0.024 (two-sided 0.048) and does not pass BH within the DNA-PKcs score tests (FDR 0.29). Reactome IFN signaling NES is +1.09 (FDR 0.38).

Ku80 knockdown versus parental does not pass FDR on this panel. Its strongest nominal increase is GO type I IFN signaling, NES +0.96, nominal p = 0.23.

Across all three subunits and all 12 score tests, the strictest BH-FDR on the one-sided p-values is 0.079 (Ku70 DNA-sensing and Ku70 IFN signaling). The FDR < 0.05 results above are the within-Ku70 families, not a 36-test family.

## Matched siControl contrast

Inside each monoallelic line, siRNA versus siControl, ranked by Welch t, does not reproduce the parental-contrast increase. Hallmark IFN-α NES is −1.20 / −1.00 / −0.85 for DNA-PKcs / Ku70 / Ku80, none with FDR < 0.05. A published STING “down” call on that rank is the depleted subunit gene; removing XRCC5, XRCC6, and PRKDC removes it. That within-line rank is a different contrast: it subtracts the monoallelic background that already differs from parental 293T.

cGAS, IFNB1, and the type III interferons remain essentially uncounted. There is no exogenous DNA in the protocol. The upregulation above is a basal shift of the expressed IFN-signaling and cytosolic DNA-sensing genes in Ku70-depleted cells relative to parental 293T.
