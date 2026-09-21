# GSE180581: DNA-PKcs vs Ku70 vs Ku80, IFN / STING / APM / DNA-sensing

Unstimulated HEK293T. Each contrast is the monoallelic line plus the matching siRNA versus the same line plus siControl (3 vs 3, 72 h, 50 nM). No exogenous DNA is in the protocol. Positive NES means the set sits at the knockdown end of a Welch *t* rank on log2(median-of-ratios + 1).

On-target depletion is in the RNA: PRKDC log2FC **−1.68** (Welch *p* = 2.9×10⁻⁵), XRCC6 **−2.01** (*p* = 0.0011), XRCC5 **−2.02** (*p* = 0.0015). The other two subunits move by at most 0.18 log2. Replicates on genes with mean count ≥ 20 correlate at Pearson ≥ 0.964. PC1 of variable genes separates the parental line from the monoallelic backgrounds; the siRNA arm does not form its own cloud, which is why the contrasts stay inside each background.

## What this matrix can score

cGAS, IFI16, ZBP1, AIM2, TREX1, IFNB1, IFNA1, IFNG, IFNL1/2/3, and CXCL10 are at or under a few raw counts and are not in the rank (count ≥ 10 in at least 3 of 21 libraries). STING1 is low (about 5–20 counts). Hallmark IFN and MHC-I sets are carried by the genes that are actually present (STAT1/2, ISG15, IRF3/7, HLA, TAP), not by an induced interferon transcriptome.

Four small GO sets (type III IFN signaling, cellular response to type III IFN, response to type I IFN, type II IFN signaling) have 6–7 genes in the rank universe and are not scored. Twenty-six sets are scored. BH-FDR is within those 26, separately for each subunit.

## Subunit comparison

| Set | DNA-PKcs NES (FDR) | Ku70 NES (FDR) | Ku80 NES (FDR) |
|---|---:|---:|---:|
| Hallmark IFN-α | −1.20 (0.23) | −1.00 (0.37) | −0.85 (0.62) |
| Hallmark IFN-γ | −1.11 (0.25) | −1.10 (0.35) | −0.98 (0.55) |
| KEGG cytosolic DNA-sensing | +0.93 (0.37) | +0.99 (0.37) | −1.16 (0.44) |
| Reactome STING, as published | −1.56 (0.15) | −1.49 (0.35) | −1.16 (0.44) |
| Reactome STING, subunit genes removed | +1.12 (0.40) | +1.16 (0.40) | +1.71 (0.087) |
| Custom MHC-I / APM (21 genes) | +0.73 (0.51) | −1.11 (0.35) | −0.72 (0.55) |
| Reactome class I MHC processing (332 genes) | −1.19 (0.15) | −1.13 (0.35) | **−1.50 (0.026)** |

One panel call passes FDR 0.05: Ku80 Reactome class I MHC-mediated antigen processing, NES −1.50, nominal *p* = 0.001, FDR 0.026, 332 genes in the rank. The leading edge is ubiquitin ligases and constitutive proteasome subunits (UBE2, FBX, PSMB2, PSMD), not HLA or TAP. The same set is negative in DNA-PKcs (NES −1.19, nominal *p* = 0.019, FDR 0.15) and Ku70 (NES −1.13, nominal *p* = 0.091, FDR 0.35). The compact MHC-I/APM list does not pass FDR in any subunit. HLA-A/B/C log2FC stays within ±0.16. Ku80 B2M is a gene-level decrease (log2FC −0.67, Welch *p* = 0.0018).

KEGG antigen processing is nominally up in all three (NES +1.44 / +1.36 / +1.52; nominal *p* = 0.028 / 0.028 / 0.008; FDR 0.15 / 0.35 / 0.10). Its leading edge is HSPA1B, HSPA5, LGMN, and CALR, with some MHC-II genes, so this is not an MHC-I program opening.

The published Reactome STING set looks down because its leading edge is the depleted subunit: PRKDC in the DNA-PKcs knockdown, XRCC6 and XRCC5 in the Ku70 knockdown. Removing XRCC5, XRCC6, and PRKDC flips that NES positive. The Ku80 residual is nominal (*p* = 0.010, sensitivity-panel FDR 0.087, 8 genes). The custom DNA-sensor list does the same thing: primary NES −1.47 / −1.55 / −1.41, and after the three subunit genes are removed +0.92 / −0.95 / −0.79, none with FDR < 0.05. KEGG cytosolic DNA-sensing does not contain the DNA-PK genes and is not significant in any subunit.

Hallmark IFN-α and IFN-γ are negative in all three subunits. DNA-PKcs Hallmark IFN-α is NES −1.20, nominal *p* = 0.053, FDR 0.23. The strongest nominal IFN decrease is DNA-PKcs GO positive regulation of type I IFN production, NES −1.36, nominal *p* = 0.023, FDR 0.15. Ku70 is the only subunit whose median IFN-set NES is positive (+0.96 versus −0.95 and −0.72), and none of its 13 IFN sets has nominal *p* < 0.05. Ku70 does not open an interferon program on this matrix. ISG15 in the Ku70 knockdown is log2FC −0.67, Welch *p* = 0.041.

Genome-wide Welch *t* Spearman between subunits is 0.32–0.34 (n ≈ 14,550). Panel NES Spearman is 0.47–0.58 (26 sets). Genes with |log2FC| > 1 and Welch *p* < 0.05: DNA-PKcs 35, Ku70 66, Ku80 63. Anisenko et al., Biochimie 2022 (PMID 35430316) reported 7, 219, and 29 genes at a >2-fold change for DNA-PKcs, Ku70, and Ku80 with their own differential-expression pipeline. Those counts are theirs. This Welch reanalysis is a different statistic and does not reproduce 219 Ku70-only genes.

## Reading

Matched depletion of DNA-PKcs, Ku70, or Ku80 in this unstimulated HEK293T series does not open IFN, STING, MHC-I/APM, or cytosolic DNA-sensing transcription. The STING and custom DNA-sensor enrichments that look down are the knockdown of the subunit gene inside the set. The one FDR < 0.05 pathway result is a Ku80 decrease of the broad Reactome class I MHC / ubiquitin-proteasome set, shared in sign with the other two subunits and not seen in the compact MHC-I gene list. DNA-PKcs is the quietest arm by the Welch |log2FC| > 1 count. HEK293T here has essentially no cGAS or type I/III interferon transcript, and the experiment has no DNA stimulus, so this is a basal-transcriptome comparison rather than a test of DNA-triggered interferon.
