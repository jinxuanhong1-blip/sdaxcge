# Sources and provenance

Accessed 2026-08-16. “Canonical” means the authors' published software/assay,
not merely the same named gene set.

| Asset | Primary source | Implementation provenance |
|---|---|---|
| TIDE, Dysfunction, Exclusion, MDSC, CAF, TAM M2 | Jiang et al., *Nature Medicine* 2018, DOI [10.1038/s41591-018-0136-1](https://doi.org/10.1038/s41591-018-0136-1); Fu et al., *Genome Medicine* 2020, DOI [10.1186/s13073-020-0721-z](https://doi.org/10.1186/s13073-020-0721-z) | `run_tide.py` delegates to the authors' [TIDEpy v1.3](https://github.com/liulab-dfci/TIDEpy). Whole-transcriptome fitted references are not replaced by marker averages. TIDEpy documents limited NSCLC validation and warns that unreferenced `log2(RPKM+1)` is not sufficient. |
| 18-gene T-cell-inflamed GEP / TIS | Ayers et al., *JCI* 2017, DOI [10.1172/JCI91190](https://doi.org/10.1172/JCI91190); Danaher et al., *JITC* 2018, DOI [10.1186/s40425-018-0367-1](https://doi.org/10.1186/s40425-018-0367-1) | Public genes are in `signatures.tsv`. The clinical score is a housekeeping-normalized weighted sum. Published validation papers state the weights are NanoString intellectual property; therefore `TIS18_public_mean` is explicitly a coefficient-free research proxy, not TIS. |
| IPS | Charoentong et al., *Cell Reports* 2017, DOI [10.1016/j.celrep.2016.12.019](https://doi.org/10.1016/j.celrep.2016.12.019) | Port of the authors' BSD-3-Clause [`IPS.R` and `IPS_genes.txt`](https://github.com/icbi-lab/Immunophenogram), retrieved 2026-08-16. HGNC aliases `IARS→IARS1` and `SDPR→CAVIN2` were updated. Per-sample z scores, factor means, signed class means and 0–10 mapping are preserved. |
| IFNG6 | Ayers et al., *JCI* 2017, DOI [10.1172/JCI91190](https://doi.org/10.1172/JCI91190); benchmarked by Lee & Ruppin, *Nature Communications* 2020, DOI [10.1038/s41467-020-18546-x](https://doi.org/10.1038/s41467-020-18546-x) | Mean log expression of `IFNG, STAT1, IDO1, CXCL9, CXCL10, HLA-DRA`, matching the public benchmark definition. |
| Higgs IFNG4 | Higgs et al., *Clinical Cancer Research* 2018, DOI [10.1158/1078-0432.CCR-17-3451](https://doi.org/10.1158/1078-0432.CCR-17-3451) | Mean of `IFNG, CD274, LAG3, CXCL9`; upper-tertile positive in the original 97-case durvalumab NSCLC cohort. The continuous score is implemented; derive a threshold only within a declared validation design. |
| CTL5 | Jiang et al., *Nature Medicine* 2018 | Mean log expression of `CD8A, CD8B, GZMA, GZMB, PRF1`; this is the CTL level, not TIDE Exclusion. |
| CAF collagen-3 | Chakravarthy et al., *Nature Communications* 2018, DOI [10.1038/s41467-018-06654-8](https://doi.org/10.1038/s41467-018-06654-8); Tyler & Tirosh, *Scientific Reports* 2023, DOI [10.1038/s41598-023-28480-9](https://doi.org/10.1038/s41598-023-28480-9) | Mean log expression of `COL1A1, COL1A2, COL3A1`. It is a transparent bulk stromal proxy and is not the TIDE CAF genome-wide score. |
| Hallmark EMT | Liberzon et al., *Cell Systems* 2015, DOI [10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004); MSigDB set [M5930](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION) | All 200 human symbols downloaded from the official GRP endpoint on 2026-08-16. `HALLMARK_EMT_mean` is the arithmetic mean on log-scale input, not ssGSEA/GSVA. |
| Thompson EMT/inflammation | Thompson et al., *Lung Cancer* 2020, DOI [10.1016/j.lungcan.2019.10.012](https://doi.org/10.1016/j.lungcan.2019.10.012) | `score_thompson_emt.py` reproduces cohort-wise gene z-scores, signed 12-gene EMT, 27-gene inflammation, `inflammation−EMT`, and `−0.60×EMT+0.19×inflammation`. The small retrospective cohort internally fit and evaluated the weights. |
| Pan-F-TBRS | Mariathasan et al., *Nature* 2018, DOI [10.1038/nature25501](https://doi.org/10.1038/nature25501) | Public 19-gene stromal set. `score_ssgsea.R` provides ssGSEA; the local mean is explicitly labeled. Historical aliases `CTGF→CCN2` and `FAM101B→RFLNB` are updated. Originally developed in urothelial cancer, not NSCLC. |
| hMENA–TGF-β CAF9 | Mercurio et al., *JITC* 2026, DOI [10.1136/jitc-2025-013098](https://doi.org/10.1136/jitc-2025-013098) | Nine public genes, exact ssGSEA route and reported `TGFB1` > cohort 10th-percentile gate in `score_ssgsea.R`; validated in OAK and SU2C-MARK. The arithmetic mean is only a QC proxy. |
| Bindea cytotoxic cells | Bindea et al., *Immunity* 2013, DOI [10.1016/j.immuni.2013.10.003](https://doi.org/10.1016/j.immuni.2013.10.003) | Public marker mean for quick QC. Faithful Bessede-style estimates require the Bindea set through `ConsensusTME::geneSetEnrichment`, as described below. |
| TACSTD2 / TROP2 NSCLC analysis | Bessede et al., *Clinical Cancer Research* 2024, DOI [10.1158/1078-0432.CCR-23-2566](https://doi.org/10.1158/1078-0432.CCR-23-2566) | POPLAR/OAK pretreatment RNA-seq; ConsensusTME/Bindea ssGSEA; optimized TACSTD2 survival split; Cox models; randomized atezolizumab-versus-docetaxel comparison. The paper did **not** calculate TIDE Exclusion. |
| EMT in PD-L1-high NSCLC ICI | Kim et al., *British Journal of Cancer* 2024, DOI [10.1038/s41416-024-02698-4](https://doi.org/10.1038/s41416-024-02698-4) | Supports EMT as an NSCLC ICI resistance analysis; use the exact paper-specific signature for strict figure replication, or declare Hallmark EMT as a prespecified substitute. |
| CAF NSCLC ICI | Hu et al., *Journal of Translational Medicine* 2023, PMID [37010552](https://pubmed.ncbi.nlm.nih.gov/37010552/) | Paper-specific fitted CAF risk model is cohort-derived. Do not relabel the portable collagen-3 or TIDE CAF score as that risk model. |

## Representative NSCLC ICI applications, 2020–2026

This is a method catalog, not a claim that every application independently
validated treatment prediction.

| Year | Method used | NSCLC ICI paper | Evidence boundary |
|---|---|---|---|
| 2020 | Immune/IFN-related signatures | Hwang et al., *Scientific Reports*, DOI [10.1038/s41598-019-57218-9](https://doi.org/10.1038/s41598-019-57218-9) | Anti-PD-1 NSCLC cohort; paper-specific M1 and peripheral-T-cell signatures, not IFNG6/TIS. |
| 2021 | T-cell-inflamed GEP | Poma et al., *Cancers*, DOI [10.3390/cancers13153828](https://doi.org/10.3390/cancers13153828) | Pembrolizumab NSCLC; used a distinct 16-of-18 NanoString-panel average. Do not treat it as Ayers' weighted TIS. |
| 2023 | NanoString TIS and IFN Gamma GEPs | Gil et al., *Translational Oncology*, DOI [10.1016/j.tranon.2023.101818](https://doi.org/10.1016/j.tranon.2023.101818) | Real-world advanced NSCLC; 66 profiled patients. Proprietary IO360 weighted scores; favorable survival associations, limited response association. |
| 2023 | CAF risk signature | Hu et al., *Cancer Immunology, Immunotherapy*, DOI [10.1007/s00262-023-03428-0](https://doi.org/10.1007/s00262-023-03428-0) | Cohort-fitted classifier evaluated in two ICB NSCLC cohorts; not interchangeable with TIDE CAF/collagen-3. |
| 2024 | TACSTD2 and Bindea/ConsensusTME | Bessede et al., *Clinical Cancer Research*, DOI [10.1158/1078-0432.CCR-23-2566](https://doi.org/10.1158/1078-0432.CCR-23-2566) | Randomized atezolizumab versus docetaxel trials plus independent ICI cohort; no TIDE analysis. |
| 2024 | EMT | Kim et al., *British Journal of Cancer*, DOI [10.1038/s41416-024-02698-4](https://doi.org/10.1038/s41416-024-02698-4) | Two ICI NSCLC cohorts; adverse association concentrated in PD-L1-high disease. |
| 2025 | TIDE Exclusion/CAF | Li et al., *Communications Biology*, DOI [10.1038/s42003-025-09087-4](https://doi.org/10.1038/s42003-025-09087-4) | Spatial fibroblast ecotypes evaluated in ICB-treated NSCLC datasets; TIDE-based secondary characterization. |
| 2025 | TIDE and Charoentong IPS | PRAS40/AKT1S1 NSCLC study, PMCID [PMC12514079](https://pmc.ncbi.nlm.nih.gov/articles/PMC12514079/) | TIDE/IPS were in-silico response proxies; a separate plasma ICB cohort did not turn those scores into randomized treatment-predictive validation. |
| 2026 | hMENA–TGF-β CAF9 | Mercurio et al., *JITC*, DOI [10.1136/jitc-2025-013098](https://doi.org/10.1136/jitc-2025-013098) | Nine-gene GSVA/ssGSEA activity developed in OAK atezolizumab and validated in SU2C-MARK anti-PD-1; retain the reported TGFB1 gate. |

The 2025 Tempus “Immune Profile Score,” DOI
[10.1136/jitc-2025-011363](https://doi.org/10.1136/jitc-2025-011363), is a
different proprietary multiomic model (TMB plus RNA features). It is **not**
Charoentong's Immunophenoscore implemented here, despite sharing the acronym
IPS.

## Scope note for 2020–2026

Many 2020–2026 NSCLC papers *apply* TIDE, IPS, Hallmark EMT, or newly fitted CAF
models to TCGA prognosis rather than validate an ICI-predictive biomarker.
This catalog prioritizes methods with a public definition and papers with an
ICI-treated NSCLC cohort. A TCGA-only association is not evidence of treatment
prediction; treatment prediction requires an interaction or comparable treated
control arm.
