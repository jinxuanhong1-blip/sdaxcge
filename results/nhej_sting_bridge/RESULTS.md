# NHEJ-component loss does not raise IFN / STING / APM transcripts

## Question

The middle step of CLDN4 → NHEJ → IFN is: loss of an NHEJ component raises
interferon, STING-axis, and MHC-I antigen-presentation (APM) transcripts.
These two public RNA-seq series perturb NHEJ components and do not perturb
CLDN4. They test only that middle step.

## What the counts show

They do not show that upregulation.

In HEK293T, Ku70, Ku80, and DNA-PKcs depletion are on target, and the IFN
program is not there to raise. cGAS, IFNβ, CXCL10, and the classic ISGs sit
at the count floor. Detected APM genes do not go up. The only STING-axis gene
with padj < 0.05 is MAVS after Ku70 depletion (log2FC +0.28).

In HeLa, the samples deposited as XRCC4(−/−) have a **lower** IFN and APM
program than scrambled-gRNA cells, in both experiments. Hallmark interferon-α
and interferon-γ are enriched among genes that go down (NES −2.86 and −2.85).
That is the opposite direction from an upward middle step. The XRCC4 mRNA
itself drops in only one of the two experiments.

## Design

**GSE180581** (Anisenko et al., Data in Brief 2021, PMID 34849385; Biochimie
2022, PMID 35430316). HEK293T with a monoallelic knockout of Ku70, Ku80, or
DNA-PKcs, then 50 nM siRNA against the remaining allele for 72 h, versus
siControl in the **same** heterozygous line. n = 3. Parental 293T + siControl
is not the control. STAR counts on hg38, tested here with PyDESeq2 (Wald,
BH). Genes need ≥10 counts in ≥3 of the 6 samples in that contrast.

**GSE135274** (Benjamin & Schiller, Int. J. Mol. Sci. 2022, PMID 35054780).
HeLa XRCC4(−/−) (guide g2G3) versus scrambled gRNA, without mirin. Two
experiments. CLC read-1 and read-2 columns are technical and were summed
before testing. Every culture, including the controls, was transfected with a
TALEN NHEJ reporter, so this is not an unchallenged basal knockout. n = 2;
PyDESeq2 warns that the dispersion prior is poorly estimated. Direction was
therefore also checked as a per-experiment library-size log2FC, which does
not use that fit.

Mirin arms (alternative-NHEJ block, and the double block the paper emphasizes)
are secondary. They are in the tables. The XRCC4 question uses the no-mirin
contrast.

## On-target check

| Contrast | Target | log2FC | padj | Normalized mean, depleted vs control | Genes at padj < 0.05 |
|---|---|---:|---:|---|---:|
| Ku70 siRNA | XRCC6 | −2.01 | 3.1×10⁻⁵² | 135 vs 539 | 785 |
| Ku80 siRNA | XRCC5 | −1.98 | 8.4×10⁻⁵⁰ | 134 vs 541 | 574 |
| DNA-PKcs siRNA | PRKDC | −1.70 | 3.1×10⁻⁶¹ | 1,126 vs 3,628 | 46 |
| XRCC4 KO, no mirin | XRCC4 | −0.54 | 0.55 | 46 vs 66 | 1,005 |

Replicate Pearson correlations on log2(normalized counts) are 0.97–0.99 in
every arm. The HEK rank order (Ku70 ≫ Ku80 ≫ DNA-PKcs by number of FDR hits)
matches the published description of these lines. The exact “>2-fold” gene
counts in that paper are not re-derived here.

XRCC4 raw counts are 59 vs 59 in experiment 1 and 33 vs 71 in experiment 2
(size-factor log2FC −0.10 and −1.11). The series metadata labels g2G3 as
XRCC4(−/−). These counts do not confirm that label in experiment 1.

## HEK293T: the IFN program is at the floor

Count-filter medians (positive would be up in the depleted cells):

| Panel | Ku70 | Ku80 | DNA-PKcs |
|---|---|---|---|
| STING axis (detected / 13) | +0.10 (6/13, 6 up) | +0.14 (6/13, 6 up) | +0.07 (6/13, 5 up) |
| IFN core (29) | −0.01 (13 detected) | +0.14 (13) | −0.12 (13) |
| APM (16) | −0.18 (14; down p = 0.029) | +0.10 (14) | −0.09 (14) |
| Hallmark IFNα median | +0.01 | +0.05 | −0.04 |
| Hallmark IFNγ median | −0.03 | −0.02 | −0.05 |

Pre-ranked GSEA on the Wald statistic (2,000 permutations; BH across the 12
Hallmark tests in this run):

- Ku70 and Ku80 Hallmark IFNα/γ: NES between −1.03 and +0.62, BH q ≥ 0.27.
- DNA-PKcs Hallmark IFNγ NES −1.48 (nominal p = 0.0055, BH q = 0.0094) and
  IFNα NES −1.44 (nominal p = 0.015, BH q = 0.022). The per-gene median is
  only about −0.04 to −0.05. That is a downward rank shift, not an induced
  program, and it is the wrong direction for the middle step.

Why the STING-axis sign test is not an induced program: the six genes that
pass the filter in the Ku depletions are STING1, TBK1, IKBKE, IRF3, IRF7, and
MAVS, each with |log2FC| mostly under 0.2. None of STING1, TBK1, IRF3, or
IRF7 has padj < 0.05. The one gene that does is **MAVS** after Ku70 depletion
(log2FC +0.28, padj = 0.0034). cGAS, IFI16, ZBP1, DDX58, IFIH1, IFNB1, and
CXCL10 do not pass the count filter. IFNB1, CXCL10, OAS2, and RSAD2 are at
zero counts in these libraries. cGAS normalized means are below 1.

APM does not go up. The clearest single APM gene is B2M after Ku80 depletion,
and it goes **down** (log2FC −0.64, padj = 0.0015). HLA-A stays near 0
(Ku70 +0.03, Ku80 −0.06, DNA-PKcs −0.08).

## HeLa XRCC4-labeled cells: IFN and APM go down

No-mirin contrast, XRCC4(−/−) versus scramble. IFN-core median log2FC
**−1.18** (1 up / 26 down; Wilcoxon p for a downward shift = 1.9×10⁻⁷). APM
median **−0.42** (4 up / 11 down; downward p = 0.0062). STING-axis median
**−0.26** (3 up / 7 down; downward p = 0.014).

Hallmark interferon-α NES **−2.86**, interferon-γ NES **−2.85**. Nominal
p = 0.0005, which is the resolution floor of 2,000 permutations (1/2001).
BH q = 0.001 across the 12 Hallmark tests.

Genes that move, and that move in the same direction in both experiments:

| Gene | DESeq2 log2FC | padj | Experiment 1 log2FC | Experiment 2 log2FC |
|---|---:|---:|---:|---:|
| ISG15 | −1.90 | 8.3×10⁻¹⁰ | −1.52 | −2.22 |
| IFI27 | −1.70 | 5.6×10⁻⁷ | −1.25 | −2.12 |
| OAS2 | −1.71 | 0.010 | −0.62 | −2.58 |
| OAS1 | −1.27 | 2.9×10⁻⁴ | −0.66 | −1.80 |
| RSAD2 | −1.13 | 0.035 | −0.39 | −1.77 |
| IFIT1 | −1.18 | 0.065 | −0.58 | −1.78 |
| STING1 | −0.69 | 0.012 | −0.78 | −0.59 |
| HLA-A | −0.65 | 1.6×10⁻⁸ | −0.76 | −0.52 |
| HLA-B | −0.87 | 5.1×10⁻¹⁵ | −0.95 | −0.77 |
| B2M | −0.45 | 8.1×10⁻⁵ | −0.45 | −0.44 |
| TAP1 | −0.42 | 0.013 | −0.22 | −0.59 |

MX1 is the exception among the core ISGs (log2FC +0.71, baseMean 16, padj not
significant; raw counts 23 vs 12 and 18 vs 12). STAT1 is mixed by experiment
(+0.19 / −0.42) and not significant in the joint fit (log2FC −0.14, padj
0.75). IFNB1 is absent (0, 0, 1, 0). cGAS is expressed (baseMean 168) and does
not move (log2FC −0.26, padj 0.62).

The mirin arms do not flip the sign. XRCC4 versus scramble inside mirin:
IFN-core median −1.41, Hallmark NES −3.12 / −3.13. XRCC4+mirin versus
untreated scramble: IFN-core median −0.36, Hallmark NES −2.08 / −2.11.

## Where this leaves CLDN4 → NHEJ → IFN

The middle arrow, as an upregulation of IFN / STING / APM transcripts after
NHEJ-component loss, is not what these two transcriptomes show.

- HEK293T can test the knockdown (the targets fall about 4-fold) and cannot
  test an IFN transcriptional response. The sensors and the ISG output are
  at the floor, and the APM genes that are detected stay flat or move down.
- HeLa expresses the program, and the XRCC4-labeled samples have less of it.
  Both experiments agree on ISG15, OAS1/2, IFI27, HLA-A/B, B2M, and STING1.
  The XRCC4 transcript drop is limited to experiment 2, and every sample
  carries a TALEN break, so this is a lower IFN/APM state in the deposited
  knockout arm, not a clean basal demonstration that removing XRCC4 turns
  the program on.

These series are not lung, and they are not a CLDN4 perturbation. They do not
speak to the CLDN4 → NHEJ step.

## Figures and tables

- `fig_nhej_ifn_panels.png` — panel medians and the on-target log2FC.
- `fig_focus_log2fc.png` — gene-level DESeq2 log2FC. A dot means the gene
  failed the count filter.
- `panel_summary.tsv`, `gsea_prerank.tsv`, `on_target.tsv`, `focus_genes.tsv`,
  `nhej_genes.tsv`, `gse135274_experiment_concordance.tsv`, `sample_sheet.tsv`.
- Full Wald tables: `*_de.tsv.gz`.

Reproduce: `python3 methods/nhej_sting_bridge/analyze.py`.
