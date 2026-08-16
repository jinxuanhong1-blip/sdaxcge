# Public lung omics hunt: TROP2 ADCs

Search date: 2026-08-16

## Bottom line

| Agent | Direct public lung treatment omics? | Tight junction (TJ) down? | Interferon (IFN) up? |
|---|---|---|---|
| sacituzumab tirumotecan (sac-TMT; SKB264; MK-2870) | **No analyzable ADC-treated dataset found.** HRA011882 profiles osimertinib-induced persisters used to motivate sac-TMT, not cells or tumors exposed to sac-TMT. | Not testable from a direct public lung perturbation dataset. | Not testable from a direct public lung perturbation dataset. |
| sacituzumab govitecan (SG; IMMU-132; Trodelvy) | **No analyzable ADC-treated lung omics dataset found.** GSE302284 profiles the osimertinib persister state; SG was tested for tumor control, but the deposited omics samples are not an SG perturbation series. | Not testable in lung. | Not testable in lung. |
| datopotamab deruxtecan (Dato-DXd; DS-1062/a) | **Yes, but controlled access:** ICARUS-LUNG01, EGAD50000001014, paired patient NSCLC biopsies. Published aggregate results are public; sample-level sequencing requires EGA approval. | **No evidence reported.** The paper did not report a TJ gene set. Tumor-compartment EMT increased, which is not equivalent to TJ loss. | **Yes, preliminary.** In seven paired spatial-transcriptomic cases, immune-compartment Hallmark IFN-alpha response increased on treatment (NES 1.69, FDR-adjusted p=0.01). |

The honest answer is therefore: **IFN up has direct, treatment-linked support for Dato-DXd in lung; TJ down does not.** There is no public direct lung omics test of the joint “TJ down / IFN up” hypothesis for sac-TMT or SG.

## Best direct evidence: ICARUS-LUNG01

Dato-DXd was given to 100 patients with pretreated advanced NSCLC. The molecular subsets were small:

- Bulk RNA-seq: 20 paired baseline/on-treatment biopsies (3 responders, 17 non-responders).
- GeoMx spatial transcriptomics: 7 paired patients, split into tumor and immune compartments.
- The bulk comparison found treatment-associated immune pathways enriched in responders relative to non-responders (FDR-adjusted p ≤ 0.05), but only three responder RNA-seq pairs were available.
- Spatial GSEA found IFN-alpha response up in the immune compartment (NES 1.69; adjusted p=0.01).
- Spatial GSEA also found EMT up in the tumor compartment (NES 1.56; adjusted p=0.0001) and KRAS signaling up (NES 1.42; adjusted p=0.006).
- The authors explicitly call these transcriptomic findings preliminary. The study is single-arm, samples were selected by tissue adequacy, on-treatment timing mixed cycle 1 day 3 and cycle 2 day 3, and there were few responders.

The raw/de-identified omics are catalogued at EGA, not open-download data. Reanalysis requires a data-access application. The public analysis-code repository does not include the analysis-ready expression tables.

## Why the TJ claim remains unproven in lung

The closest mechanistic evidence is outside lung:

- GSE334497 (mouse TNBC) compares TROP2 wild-type with TROP2 knockout tumors. TROP2 loss reduced claudin/tight-junction programs and increased inflammatory, antitumor-immune, and T-cell cytotoxicity programs.
- The same breast study used naked hRS7 (the antibody component of SG), not ADC-treated transcriptomics, and showed reduced claudin 7 plus improved T-cell accessibility with anti-PD-1.

That makes “TROP2 blockade can weaken a claudin barrier” biologically plausible, but it cannot be relabeled as evidence that sac-TMT, SG, or Dato-DXd causes TJ down in lung. TROP2 genetic knockout, naked antibody, ADC payload effects, tissue context, and combination with checkpoint blockade are different perturbations.

## Dataset notes

See `datasets.tsv` for a machine-readable inventory and `claims.tsv` for the claim-level verdicts.

Important naming distinction: “sacituzumab” is ambiguous. Sacituzumab tirumotecan (SKB264/MK-2870; belotecan-derived payload) and sacituzumab govitecan (IMMU-132/Trodelvy; SN-38 payload) share an anti-TROP2 antibody lineage but are different ADCs.

## Sources

- ICARUS-LUNG01 article: https://doi.org/10.1016/j.ccell.2026.03.017
- ICARUS-LUNG01 EGA dataset: https://ega-archive.org/datasets/EGAD50000001014
- ICARUS-LUNG01 bulk RNA-seq study: https://ega-archive.org/studies/EGAS50000000732
- ICARUS-LUNG01 GeoMx study: https://ega-archive.org/studies/EGAS50000001679
- Public ICARUS analysis code and access instructions: https://github.com/gustaveroussy/ICARUS-01_Public
- Sac-TMT persister-state dataset: https://ngdc.cncb.ac.cn/gsa-human/browse/HRA011882
- Sac-TMT lung paper: https://pubmed.ncbi.nlm.nih.gov/42314664/
- SG lung persister-state dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302284
- SG lung paper: https://doi.org/10.1158/2159-8290.CD-24-1515
- TROP2/claudin breast mechanism dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334497
- TROP2/claudin paper: https://doi.org/10.1136/jitc-2025-012265
- Direct SG ovarian dataset (excluded by tissue): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278664
- Dato-DXd breast PDX dataset (excluded by tissue): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235812

## Search scope

Aliases were searched across GEO/NCBI Entrez, SRA/BioProject, OmicsDI, ArrayExpress/BioStudies, EGA, GSA-Human/NGDC, Europe PMC/PubMed, ClinicalTrials.gov, article data-availability statements, and general web indexing. Search terms included SKB264, SKB-264, sac-TMT, sacituzumab tirumotecan, MK-2870, sacituzumab govitecan, IMMU-132, Trodelvy, Dato-DXd, datopotamab deruxtecan, DS-1062, and DS-1062a.

This is a documented search, not proof that an unindexed sponsor dataset does not exist. “None found” means no public sample-level dataset matching both lung and direct exposure was identifiable on the search date.
