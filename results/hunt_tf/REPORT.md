# Hunt: ELF3 / GRHL1 / KLF4 / TFAP2A with TACSTD2 and CLDN4, NKX2-1 opposite

Honest public-data search. Co-expression is not evidence of direct transcriptional regulation.

## Claim being tested

In public **human and mouse lung RNA**, the user TF set **ELF3, GRHL1, KLF4, TFAP2A** co-correlates with **TACSTD2** and **CLDN4**, and **NKX2-1** anti-correlates.

## What would count as support (fixed before opening expression matrices)

Pair (raw Spearman on log2(CPM+1), n ≥ 12 after filters):

- positive pair: ρ ≥ 0.30 and two-sided p < 0.05
- NKX2-1 pair: ρ ≤ −0.20 and two-sided p < 0.05

Cohort call:

- **SUPPORTED** — ≥ 6 of 8 positive pairs **and** both NKX2-1 pairs anti
- **PARTIAL** — ≥ 4 of 8 positive pairs (NKX2-1 not required)
- **NOT_SUPPORTED** — otherwise
- **UNINFORMATIVE** — targets/TFs mostly undetected, or n < 12

The same rules are applied after rank-partialling **EPCAM** (epithelial content) and, separately, **SFTPC** (AT2 / alveolar content). A raw SUPPORTED call that becomes NOT_SUPPORTED after EPCAM residualization is **composition-confounded**, not a validated network.

Score (continuous, used for ranking only):
`mean ρ(4 TFs × 2 targets) − mean ρ(NKX2-1 × 2 targets)`.

A 300-shuffle permutation of TF sample labels gives an empirical p for that score. It is a within-cohort sanity check, not a study-wide FDR.

## Search space

Source: [recount3](https://rna.recount.bio/) open data (human Gencode v26 / G026, mouse Gencode vM23 / M023). Same pipeline for both species.

1. All SRA studies with ≥ 6 runs in the recount3 project index (6,831 human, 8,195 mouse).
2. Per-study SRA metadata text is searched for lung / airway / pulmonary disease terms.
3. Runs are labelled tissue / culture-organoid / cell line / sorted-or-fluid. Single-cell studies are flagged and not downloaded.
4. Download list (pre-specified, not chosen by looking at correlations):
   - **Anchors (always):** GTEx lung, TCGA LUAD, TCGA LUSC
   - **SRA tissue:** ≥ 12 tissue lung RNA-seq runs, single-cell fraction < 0.25, median spots ≥ 5×10^5, and either ≥ 40% of runs are lung-tagged or the lung tag is run-level
   - **SRA culture / cell line:** culture/organoid with the same n/depth filters; cell-line panels (n ≥ 24 or a multi-line title). Single-line nanoparticle time courses are mostly excluded.
   - Sorted cells, blood, BAL, platelets, and isolated immune studies are excluded

This is a screen of public lung RNA that recount3 already processed. It is not every GEO series, and it is not single-cell ATAC or spatial data.

## Why bulk lung can fake this network

TACSTD2 and CLDN4 are high in airway / regenerating epithelium. NKX2-1 and SFTPC are high in AT2 / alveolar lineage. ELF3, GRHL1, KLF4 and TFAP2A are epithelial transcription factors with airway-leaning published roles. In mixed lung tissue the same sample-to-sample swing in **airway vs alveolar fraction** will raise ELF3/GRHL1/KLF4/TFAP2A/TACSTD2/CLDN4 together and lower NKX2-1, without any of those TFs writing TACSTD2 or CLDN4.

That is why EPCAM- and SFTPC-partial correlations, and the culture/cell-line slice, are required before calling the network supported.


## Results

Scored **248** cohort slices (135 human, 113 mouse). **240** were informative under the pre-specified filters.

- Raw calls: NOT_SUPPORTED 154, PARTIAL 83, SUPPORTED 3
- EPCAM-residual calls: NOT_SUPPORTED 190, PARTIAL 47, SUPPORTED 3
- SFTPC-residual calls: NOT_SUPPORTED 171, PARTIAL 68, SUPPORTED 1

### Verdict

**The pre-specified 4-TF + NKX2-1-anti package is not supported** in public human or mouse lung RNA.

3 of 240 informative slices meet the raw SUPPORTED rule; 1 still meets it after EPCAM residualization. That is not a validation. The two largest clean human matrices fail the full rule, mouse has zero SUPPORTED calls, and the few raw hits are small epithelial-cell / EMT / single-lab tumor series, not independent bulk parenchyma.

What actually holds, and what does not:

- **ELF3** is the most consistent positive correlate of TACSTD2 and especially CLDN4 (GTEx, TCGA LUAD/LUSC, CCLE, LUAD line catalogue).
- **GRHL1** tracks both targets in TCGA tumors and cell-line panels, but is weakly *negative* in GTEx lung. It is not a universal lung co-correlate.
- **KLF4** is mixed: present in GTEx and some line panels, absent or reversed in TCGA LUAD vs CLDN4.
- **TFAP2A** is not part of this module. It is near-zero in GTEx and often fails in tumors and CCLE.
- **NKX2-1 anti-correlation is rejected.** In GTEx lung NKX2-1 is *positively* correlated with TACSTD2 (ρ=0.56) and CLDN4 (ρ=0.71). In a 160-line lung cancer panel it tracks the same genes (ρ≈0.90). TCGA tumors are near zero / slightly positive. Apparent anti-correlation shows up in infection, EMT, and injury time courses where AT2 / NKX2-1 programs drop — that is biology of damage, not a TF network.
- A smaller claim — ELF3 ± GRHL1 ± KLF4 co-vary with TACSTD2/CLDN4 as an epithelial/barrier state — is visible in tumors and cell-line panels. That is still co-expression, not regulation.

### Anchor cohorts (always downloaded)

| cohort | n | raw call | EPCAM call | mean ρ TFs | mean ρ NKX2-1 | score | p_perm |
|---|---:|---|---|---:|---:|---:|---:|
| GTEx_LUNG | 655 | PARTIAL | NOT_SUPPORTED | 0.20 | 0.64 | -0.43 | 1.00 |
| TCGA_LUAD_tumor | 542 | PARTIAL | PARTIAL | 0.27 | 0.10 | 0.17 | 0.00 |
| TCGA_LUAD_normal | 59 | NOT_SUPPORTED | NOT_SUPPORTED | 0.27 | 0.12 | 0.16 | 0.13 |
| TCGA_LUSC_tumor | 504 | PARTIAL | PARTIAL | 0.27 | 0.03 | 0.24 | 0.00 |
| TCGA_LUSC_normal | 51 | NOT_SUPPORTED | PARTIAL | 0.25 | 0.31 | -0.06 | 0.66 |

#### GTEx_LUNG pairwise raw Spearman

| TF | TACSTD2 ρ (p) | CLDN4 ρ (p) |
|---|---|---|
| ELF3 | 0.41 (0.0e+00) | 0.65 (0.0e+00) |
| GRHL1 | -0.11 (4.6e-03) | -0.10 (9.9e-03) |
| KLF4 | 0.32 (0.0e+00) | 0.38 (0.0e+00) |
| TFAP2A | -0.02 (6.6e-01) | 0.08 (5.3e-02) |
| NKX2-1 | 0.56 (0.0e+00) | 0.71 (0.0e+00) |

#### TCGA_LUAD_tumor pairwise raw Spearman

| TF | TACSTD2 ρ (p) | CLDN4 ρ (p) |
|---|---|---|
| ELF3 | 0.36 (0.0e+00) | 0.52 (0.0e+00) |
| GRHL1 | 0.39 (0.0e+00) | 0.37 (0.0e+00) |
| KLF4 | 0.17 (1.0e-04) | -0.07 (9.0e-02) |
| TFAP2A | 0.25 (0.0e+00) | 0.18 (0.0e+00) |
| NKX2-1 | -0.05 (2.9e-01) | 0.25 (0.0e+00) |

#### TCGA_LUSC_tumor pairwise raw Spearman

| TF | TACSTD2 ρ (p) | CLDN4 ρ (p) |
|---|---|---|
| ELF3 | 0.48 (0.0e+00) | 0.49 (0.0e+00) |
| GRHL1 | 0.39 (0.0e+00) | 0.27 (0.0e+00) |
| KLF4 | 0.36 (0.0e+00) | 0.17 (1.0e-04) |
| TFAP2A | 0.13 (2.9e-03) | -0.13 (3.1e-03) |
| NKX2-1 | -0.04 (3.4e-01) | 0.11 (1.4e-02) |

#### SRP186687 pairwise raw Spearman

| TF | TACSTD2 ρ (p) | CLDN4 ρ (p) |
|---|---|---|
| ELF3 | 0.60 (0.0e+00) | 0.75 (0.0e+00) |
| GRHL1 | 0.48 (0.0e+00) | 0.65 (0.0e+00) |
| KLF4 | 0.46 (0.0e+00) | 0.42 (0.0e+00) |
| TFAP2A | 0.25 (4.0e-04) | 0.05 (4.7e-01) |
| NKX2-1 | 0.05 (4.5e-01) | 0.29 (0.0e+00) |

#### DRP001919 pairwise raw Spearman

| TF | TACSTD2 ρ (p) | CLDN4 ρ (p) |
|---|---|---|
| ELF3 | 0.73 (0.0e+00) | 0.68 (2.0e-04) |
| GRHL1 | 0.69 (1.0e-04) | 0.69 (1.0e-04) |
| KLF4 | 0.37 (6.5e-02) | 0.42 (3.2e-02) |
| TFAP2A | 0.48 (1.3e-02) | 0.23 (2.5e-01) |
| NKX2-1 | -0.06 (7.7e-01) | 0.26 (2.0e-01) |

### Per-TF pair pass rates (informative tissue, tumor, cell-line)

| TF | target | n pairs | median ρ | fraction ρ≥0.30, p<0.05 | fraction ρ≤−0.20, p<0.05 |
|---|---|---:|---:|---:|---:|
| ELF3 | TACSTD2 | 194 | 0.39 | 0.49 | 0.08 |
| ELF3 | CLDN4 | 209 | 0.50 | 0.56 | 0.03 |
| GRHL1 | TACSTD2 | 196 | 0.26 | 0.38 | 0.09 |
| GRHL1 | CLDN4 | 211 | 0.19 | 0.35 | 0.11 |
| KLF4 | TACSTD2 | 197 | 0.08 | 0.28 | 0.13 |
| KLF4 | CLDN4 | 212 | -0.00 | 0.24 | 0.20 |
| TFAP2A | TACSTD2 | 194 | 0.14 | 0.24 | 0.06 |
| TFAP2A | CLDN4 | 209 | 0.18 | 0.25 | 0.09 |
| NKX2-1 | TACSTD2 | 175 | 0.16 | 0.32 | 0.17 |
| NKX2-1 | CLDN4 | 183 | 0.19 | 0.33 | 0.18 |

### Human vs mouse tissue

- **human tissue/tumor:** n=70, raw NOT_SUPPORTED 38, PARTIAL 29, SUPPORTED 3; EPCAM NOT_SUPPORTED 50, PARTIAL 18, SUPPORTED 2; median score 0.04
- **mouse tissue/tumor:** n=105, raw NOT_SUPPORTED 78, PARTIAL 27; EPCAM NOT_SUPPORTED 91, PARTIAL 14; median score 0.01

### Culture / cell-line slice (less composition)

65 informative culture/cell-line slices. Raw: NOT_SUPPORTED 38, PARTIAL 27. EPCAM: NOT_SUPPORTED 49, PARTIAL 15, SUPPORTED 1.

| cohort | material | n | raw | EPCAM | score | title |
|---|---|---:|---|---|---:|---|
| SRP101626 | culture_or_organoid | 16 | PARTIAL | PARTIAL | 1.05 | Expression data from fresh human embryonic lung epithelial tip and stalk cells a |
| ERP009282 | culture_or_organoid | 18 | PARTIAL | SUPPORTED | 1.04 | Systems biology approaches reveal low-dose effects of nanoparticles |
| SRP046226 | culture_or_organoid | 24 | PARTIAL | PARTIAL | 0.92 | Phenotypic responses of differentiated asthmatic human airway epithelial culture |
| SRP048565 | culture_or_organoid | 18 | NOT_SUPPORTED | NOT_SUPPORTED | 0.63 | Host transcriptome analysis of Aspergillus fumigatus infection in Airway Epithel |
| SRP106050 | culture_or_organoid | 48 | PARTIAL | PARTIAL | 0.50 | RNA-seq of TRRAP knock-down HBEC ALI cultures |
| SRP096589 | culture_or_organoid | 44 | PARTIAL | PARTIAL | 0.47 | Gene signature profiles in respiratory epithelium infected with nontuberculous m |
| SRP091771 | cell_line | 15 | NOT_SUPPORTED | NOT_SUPPORTED | 0.46 | Gene expression profiling associated with knockdown of RNF20 in human normal and |
| DRP001919 | cell_line | 26 | PARTIAL | NOT_SUPPORTED | 0.44 | Omics catalogue of lung adenocarcinoma cell lines |

### Composition diagnostic

Mean raw Spearman across informative cohorts:

| pair | mean ρ |
|---|---:|
| TACSTD2 vs SCGB1A1 (club) | 0.15 |
| TACSTD2 vs SFTPC (AT2) | 0.11 |
| CLDN4 vs SCGB1A1 | 0.09 |
| CLDN4 vs SFTPC | 0.05 |
| NKX2-1 vs SFTPC | 0.42 |
| NKX2-1 vs SCGB1A1 | 0.21 |

If TACSTD2/CLDN4 track SCGB1A1 and anti-track SFTPC, while NKX2-1 tracks SFTPC, the hunt recovered the airway-vs-alveolar axis, not a new TF circuit.

### The three raw SUPPORTED hits (do not over-read them)

| cohort | n | EPCAM | note |
|---|---:|---|---|
| SRP066794 | 24 | PARTIAL | EMT time course in culture (mislabelled tissue). TFAP2A goes the *wrong* way (ρ≈−0.8). ELF3/GRHL1/KLF4 collapse together during EMT. |
| SRP076732 | 18 | SUPPORTED | Title is human lung epithelial cells, n=18. ELF3 vs TACSTD2 is 0.01. Only double-SUPPORTED hit; too small and not parenchyma. |
| SRP223534 | 29 | PARTIAL | 29 lung tumor RNAs from one lab, multiple aliquots per patient. Batch/patient structure can inflate ρ. |

### EPCAM-residual SUPPORTED (includes cohorts that were only PARTIAL on raw)

| cohort | raw | n | score | title |
|---|---|---:|---:|---|
| SRP076732 | SUPPORTED | 18 | 1.28 | Transcriptome analysis of human lung epithelial cells |
| ERP009282 | PARTIAL | 18 | 1.04 | Systems biology approaches reveal low-dose effects of nanoparticles |
| SRP157975 | PARTIAL | 48 | 0.91 | Ex Vivo Lung Perfusion as a Human Platform for Preclinical Small Molecule Testin |

None of these is a large independent bulk-lung parenchyma series. They do not rescue the 4-TF + NKX2-1 package.

## What this is not

- Not ChIP, CUT&RUN, motif, or reporter evidence.
- Not single-cell (those studies were excluded; recount3 gene sums of scRNA-seq are the wrong object).
- Not a claim that GRHL1 or ELF3 is *the* TACSTD2/CLDN4 factor. Sister PR A10 tests GRHL1 in TCGA with a different matrix.
- recount3 ends in 2020-era SRA. Later GEO lung series are not in this hunt.


## How to rerun

```bash
# metadata cache (once): scripts/00_fetch_metadata.sh
python3 scripts/01_select_lung_studies.py
python3 scripts/02_choose_cohorts.py
python3 scripts/03_fetch_matrices.py
python3 scripts/04_analyze.py
python3 scripts/05_figures.py
python3 scripts/06_write_report.py
```

Gene IDs are in `scripts/genes.json`.
