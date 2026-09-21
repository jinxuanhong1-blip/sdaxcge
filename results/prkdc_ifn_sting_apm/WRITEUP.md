# DNA-PKcs inhibitor and PRKDC knockdown RNA-seq: IFN, STING, APM

## Result

Across eight epithelial-cancer studies, loss of DNA-PKcs activity does not move IFN or the antigen-presentation machinery as a shared program. The study-level random-effects mean log2 change (positive = higher after inhibitor or knockdown) is:

| signature | studies up | mean log2 | 95% CI | P | BH q | I² |
|---|---|---|---|---|---|---|
| IFN (Hallmark α+γ, MHC/APM genes removed; 206 genes) | 4/8 | +0.010 | −0.079 to +0.099 | 0.83 | 0.83 | 0.99 |
| STING transcriptional readout (16 genes) | 5/8 | +0.131 | −0.039 to +0.300 | 0.13 | 0.39 | 0.98 |
| APM (21 genes: HLA, TAP, immunoproteasome, NLRC5) | 3/8 | +0.008 | −0.065 to +0.081 | 0.83 | 0.83 | 0.96 |

I² near 1 means these averages are not a common effect. Individual studies go both ways. IFN and APM stay compatible with zero if each cell line is entered separately (13 lines), if gene log2 changes are winsorized at ±2, and if the meta is split into selective inhibitors versus PRKDC knockdown or knockout. STING is the exception at cell-line level: +0.15 (0.01 to 0.30), P=0.035, I²=0.97, 10/13 lines up. That estimate counts the six SCLC lines as six votes. Pooling them first, which is the study-level contrast, puts the STING interval back across zero.

What does move, in some models only:

- **HCT116 DNA-PK knockout** (GSE285698, normoxia) raises APM (+0.33) and the classical MHC-I quartet HLA-A/B/C/B2M (+0.58). IFN is +0.18. The STING score is +0.48 but rests on 4 detected genes (IRF7 +1.34, CCL5 +0.40). IFNB1 and CXCL10 are undetected.
- **8505C thyroid, KU-57788, 24 h** (GSE319513) raises HLA-A +0.68, HLA-B +0.66, B2M +1.05 (quartet +0.75) while the full APM mean is −0.04 and IFN is −0.11. STING is +0.22. PRKDC mRNA itself falls −0.80, more than a pure catalytic block usually does.
- **C4-2, NU7441** (GSE116765) is flat for IFN (+0.01) and STING (0.00). APM is +0.14 and the MHC-I quartet is +0.25. The TORKi specificity control CC-223 moves the quartet by +0.29, so this bump is not DNA-PK-specific.
- **CWR22Rv1-AR-EK siDNA-PKcs** (GSE242255; PRKDC −3.32) gives a small IFN rise (+0.11) and a clear MHC-I drop (quartet −0.78, HLA-A −1.58). AZD7648 on the same cells is near flat (IFN +0.05, STING +0.04, quartet +0.21).
- **SUM159 shDNA-PKcs, normoxia** (GSE167956; PRKDC −2.26) is slightly down (IFN −0.07, STING −0.03, APM −0.10). Hypoxia is more negative.
- **Six SCLC lines, NU7441, 5 days** (GSE273409), pooled as one study: IFN +0.01, STING +0.13, APM −0.09. H196 is the STING outlier (CXCL10 +1.01, CCL5 +1.29) without an APM increase. H146 is the largest IFN change (+0.23). The six-line APM mean does not rise.

Selective-inhibitor studies only (5 studies): IFN −0.028 (P=0.51), STING +0.103 (P=0.19), APM +0.011 (P=0.78). The MHC-I quartet is the one consistently positive slice in that drug set: 5/5 studies have a positive point estimate, random-effects mean +0.32 (0.03 to 0.62), P=0.032, I²=0.97. Genetic knockdown does not reproduce it (mean +0.01, P=0.93): the prostate siRNA and the SUM159 hairpin lower HLA, and the HCT116 knockout raises it.

IFNB1 is essentially absent from these epithelial lines, so a STING score here is not an interferon-beta induction. It is the mean of whichever of CCL5, CXCL10, IRF7, NFKBIA and related genes are detected (often 4–11 of 16).

## Non-epithelial cancer, and HEK293T

Non-epithelial cancer (liposarcoma peposertib, HOS osteosarcoma siPRKDC, SK-N-MC NU7441, U937 NU7441; TC32 siPRKDC excluded because PRKDC only fell −0.23) does not show induction. The selective-inhibitor subset is negative for IFN (−0.28, 0/4 up, P=0.008). That signal is the 40 µM NU7441 Ewing arrays (GSE85202) plus U937, not peposertib. TC32’s genome-wide median log2 change is −0.24, and MX1 falls by several log2 units, so the IFN mean is partly a high-dose downward shift plus outlier probes. Winsorizing gene log2 changes at ±2 still leaves the drug-subset IFN mean at −0.22. LPS853 and T778 peposertib (0.5 µM, no doxorubicin) are within ±0.05 for IFN and APM.

HEK293T DNA-PKcs knockdown (GSE180581, matched heterozygous background; PRKDC −1.72) does not induce the programs (IFN −0.09, STING −0.05, APM −0.13). It is not in the cancer meta.

U937 plus transfected DNA, with or without NU7441 at 8 h, changes IFN by −0.05. That is a different question (does the inhibitor block a DNA stimulus) and is not in the meta.

## How studies were chosen

GEO was searched on 2026-09-21 for DNA-PKcs, PRKDC, NU7441, NU7026, AZD7648, M3814/peposertib, KU-57788 and knockdown/knockout language, restricted to expression profiling by high-throughput sequencing, plus the matching array series. Inclusion required a DNA-PKcs catalytic inhibitor (NU7441/KU-57788, NU5455, AZD7648, peposertib) or a PRKDC/DNA-PKcs knockdown, knockout, or shRNA, in an epithelial or cancer cell model, with at least two samples per arm and a genome-wide matrix.

The epithelial-cancer meta uses one contrast per study. Where a study has both a genetic perturbation and a drug, the genetic arm is the meta contrast and the drug is reported beside it. Six SCLC lines are pooled to one study estimate before the across-study meta, so one deposit does not cast six votes. CC-115 is dual DNA-PK/mTOR and stays out of the selective-inhibitor meta. CC-223 is a TORKi specificity control.

Not scored, with the reason in `tables/screen.tsv`: osteoblasts, primary B cells, Ku-only depletion, PTK7 or CHK2 or PRMT1 perturbations, CRISPR screens, CUT&Tag, 3-hydroxyflavone in AML, and C4-2 AZD7648 SLAM-seq (GSE287819; RAW only, n=2, no gene-count matrix).

## Score

For each contrast, a gene in the signature is kept if it is detected in at least half the samples of one arm. The sample score is the mean log2 expression of those genes. The contrast effect is the difference of arm means. Raw counts use log2(CPM+1). Normalized counts, salmon values, FPKM and array intensities use log2(x+1). GSE129436 is already a signed log matrix and is not logged again. Genetic arms must lower PRKDC by at least 0.3 log2 or they are dropped (TC32 siPRKDC). Catalytic inhibitors are not required to lower PRKDC mRNA. The genome-wide median log2 change is near zero for the RNA-seq contrasts that enter the epithelial meta, so those signature shifts are not a global expression offset.

Across-study pooling is DerSimonian-Laird random effects on the log2 changes. BH q is across the three pre-specified signatures (IFN, STING, APM) in the epithelial study-level meta.

## Files

- `tables/contrasts.tsv` — every scored arm
- `tables/meta.tsv` — random-effects summaries
- `tables/epithelial_study_estimates.tsv` — one row per study
- `tables/leave_one_study_out.tsv`
- `tables/gene_lfc.tsv` — PRKDC, ISGs, HLA, STING machinery
- `tables/screen.tsv` — include / exclude log
- `figures/01_forest_epithelial_study.png`
- `figures/02_forest_epithelial_cell_lines.png`
- `figures/03_signature_heatmap.png`
- `figures/04_prkdc_qc.png`

Reproduce: `bash results/prkdc_ifn_sting_apm/scripts/download.sh && python3 results/prkdc_ifn_sting_apm/scripts/score_meta.py`
