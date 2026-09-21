# Results: no CRISPRi count matrix

Villagomez et al., Scientific Reports, 10 Nov 2025. [doi:10.1038/s41598-025-23137-1](https://doi.org/10.1038/s41598-025-23137-1). PMID 41214101.

**Final.** This paper deposited no raw CRISPRi RNA-seq or proteome. GEO, SRA, BioProject, GSA, PRIDE, and OmicsDI have no matrix for the OVCAR3 / OVCA429 CLDN4 knockdown. The pSTING and ISRE results below are the authors’ prose and figure legends. No fold-change was printed, and none is inferred from the images.

Perturbation: CRISPRi with dCas9-KRAB-MeCP2 in OVCAR3 and OVCA429 (guide `GCTGGCTTGCGCATCAGGAC`). OVCAR8, which lacks endogenous CLDN4, carries the overexpression arm.

## What each assay says

| Assay | Condition | Direction in the text | Where |
|---|---|---|---|
| pSTING, confocal | OVCAR8 CLDN4 overexpression | Apparent decrease | Fig. 2a |
| pSTING, confocal | OVCAR3 and OVCA429 CLDN4 knockdown | Noticeable increase | Fig. 2b, c |
| pSTING, immunoblot | CLDN4 overexpression | Apparently elevated | Fig. 2e |
| pSTING, immunoblot | OVCAR3 and OVCA429 knockdown | Further increased, described as consistent with the knockdown confocal images | Fig. 2e |
| Total STING and pTBK1, immunoblot | Overexpression | Decreased total STING and decreased pTBK1 | Supplementary Fig. 6a |
| Total STING and pTBK1, immunoblot | OVCAR3 knockdown | Reduced STING and reduced pTBK1 | Supplementary Fig. 6b |
| STING and pTBK1, immunoblot | OVCA429 knockdown | “Distinct profile”; no direction is stated | Supplementary Fig. 6c |
| pSTING, immunoblot, CMP 48 h | OVCAR8 | Reduced | Fig. 2e |
| pSTING, immunoblot, CMP 48 h | OVCAR3 | Increased | Fig. 2e |
| ISRE luciferase, basal | OVCAR8 overexpression | Significant increase (7 experiments) | Fig. 3a, b |
| ISRE luciferase, basal | OVCAR3 knockdown | Dramatic reduction (6 experiments; unpaired t-test and Mann–Whitney) | Fig. 3c |
| ISRE luciferase, cGAMP | OVCAR8 overexpression versus wild type | No significant further increase | Fig. 3d |
| ISRE luciferase, cGAMP | OVCAR3 knockdown | Failed to reach the wild-type level (3 experiments) | Fig. 3e |

The knockdown arm is the conflict. Confocal and immunoblot both say pSTING is higher after CLDN4 loss in OVCAR3 and OVCA429. The ISRE reporter in OVCAR3 says the type I interferon response is lower, and cGAMP does not restore it to the wild-type level. The authors call claudin-4 a positive regulator of that interferon response on the basis of the reporter.

The overexpression arm is not one direction. Fig. 2a says the pSTING confocal signal falls. The Fig. 2e sentence says phosphorylated STING is apparently elevated on the immunoblot of the same overexpression. Both sentences are in the paper. The text gives no densitometry that chooses between them. Basal ISRE still rises with overexpression (Fig. 3a).

pTBK1 does not track the ISRE split. The paper reports lower pTBK1 with overexpression (Supplementary Fig. 6a), where ISRE is higher, and lower pTBK1 in the OVCAR3 knockdown (Supplementary Fig. 6b), where ISRE is lower.

## How the paper places the two readouts

The lysosome paragraph is the authors’ statement that a pSTING signal need not mean an active interferon program. After the confocal result they write that pSTING co-localizes with LAMP1 and appears enclosed in lysosomes, “suggesting that pSTING may be inactivated through lysosomal sequestration.” They then write that claudin-4 co-localizes with LAMP1 and that “the effect of claudin-4 on pSTING is related to its inactivation through lysosomes.” That sentence covers the co-localization. It does not report a count of sequestered versus signaling pSTING in the knockdown, and it does not quote an ISRE number.

The same results section keeps the interferon claim on the reporter: knockdown causes a dramatic reduction of the type I interferon response, and the authors treat claudin-4 as a positive regulator of that response.

Rab7 is their proposed partner, from the earlier OVCAR3 BioID. In the knockdown, transient Rab7 changed the pSTING image relative to wild type: in OVCAR3 CLDN4-knockdown cells, pSTING “appeared reduced compared to WT cells, opposite to the phenotype previously observed (Fig. 2b).” On the reporter, Rab7 transfection lowered ISRE when claudin-4 was present (significant in the overexpressing cells; p = 0.0538 in OVCAR3 wild type) and “no such effect was observed in OVCAR3 claudin-4 KD cells” (Fig. 3f, g).

The discussion states the branch they consider unresolved. cGAMP-activated STING can drive type I interferon through TBK1/IRF3 or autophagy through LC3 independent of TBK1, and micronuclear DNA could push one while suppressing the other. They write that claudin-4’s control of genome stability and micronuclei clearance “may represent a critical link,” and that “further investigation is needed to fully elucidate its impact on the cGAS-STING pathway.”

## What this file does not contain

No new RNA-seq, proteome, or densitometry. Open-cohort Spearman tables remain in `INVENTORY.md` and are tumor correlations, not the CRISPRi contrast.
