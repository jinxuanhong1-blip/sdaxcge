# 2024–2026 NSCLC scRNA hunt (CLDN4 and T cells, clinical or ICI)

Search date 2026-09-21. GEO DataSets query:

`(NSCLC OR LUAD OR "non-small cell lung" OR "lung adenocarcinoma") AND ("single cell" OR scRNA OR scRNA-seq) AND gse[Entry Type] AND 2024:2026[PDAT]`

224 series. Human series were read for design, sample characteristics, and supplementary-file size. A series is usable for malignant CLDN4 vs T/NK only if the same samples contain tumor epithelium and T/NK cells on a whole-transcriptome assay, with clinical or ICI labels, and an open count matrix.

Locked cohorts were not reopened: GSE123902, GSE131907, GSE205335, GSE189357, GSE207422. Already scored in this repo and not rerun: GSE291670, GSE241934, GSE253013.

## Ranked open candidates

| Rank | Accession | What it is | n | CLDN4 + T/NK | Clinical / ICI | Access | Decision |
|---|---|---|---|---|---|---|---|
| 1 | [GSE325414](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325414) | Jimenez et al. PEF treat-and-resect, BD Rhapsody WTA + author labels | 25 donors, 156,467 cells | Yes. Tumor epithelium, T, CD8, NK | Clinical stage IA2–IB, histology. Not ICI. All deposited libraries are the PEF arm | Open matrix 774 MB + cell metadata | **Run.** Primary = RESRU untreated tumor |
| 2 | [GSE233203](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233203) | Stage IV LUAD pleural fluid, ABCP (atezolizumab + bevacizumab + chemo) | 7 (3 response, 4 non-response) | Likely unsorted effusion (tumor + immune). Not solid tumor | ICI response | Open 10x MTX, 519 MB | Best unused **ICI-response** series. Not run: effusion, n=7 |
| 3 | [GSE253718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253718) | EGFR-mutant LUAD tumor, all cells | 6 (3 naive, 3 TKI-resistant) | Yes, unsorted tumor | Treatment group, not ICI | Open 10x MTX, 760 MB | Smaller n than GSE325414 |
| 4 | [GSE333596](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333596) | Stage I vs III LUAD tumor + matched normal | 6 tumors | Fresh tumor 10x, expected mixed | Stage | Open RAW.tar, 3.6 GB | Clinical stage, small n |
| 5 | [GSE293360](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293360) | T1N0 vs T1N2 NSCLC plus N2 nodes | 7 primary tumors | Tumor scRNA | N stage | Open RAW.tar, 947 MB | Clinical N, modest n |
| 6 | [GSE164789](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164789) | Precursor LUAD atlas (AIS/MIA/IAC era), tumor + adjacent | ~30 tumor libraries | Mixed tumor scRNA | Lesion stage, not ICI | Open RAW.tar, 2.0 GB | Large, but precursor and no deposited patient labels |
| 7 | [GSE194070](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194070) | Untreated NSCLC ± COPD, tumor + adjacent | 8 tumors | Count matrix of tumor scRNA | COPD, untreated | Open counts, 69 MB | Clinical comorbidity, not ICI |
| 8 | [GSE267108](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267108) | LUAD grouped by PD-L1 IHC | 8 tumors | Processed single-cell object | PD-L1, not response | Open processed tar, 249 MB | ICI-adjacent biomarker only |
| 9 | [GSE274934](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274934) | Total cells, WT vs ALK+ LUAD | 9 scRNA tumors | Unsorted lung | Title is immunotherapy; GEO labels are genotype, not response | Open RAW.tar, 14 GB | Too large for the gain in n; response not on the samples |
| 10 | [GSE337519](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337519) | Paired pre/post neoadjuvant chemo-immunotherapy | 1 patient | 10x tumor | ICI, paired | Open RAW.tar, 33 MB | n=1 cannot support a Spearman |

## Open but cannot score malignant CLDN4

| Accession | Why it stops | n | Annotation |
|---|---|---|---|
| [GSE243013](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013) | CD45+ immune atlas. No malignant epithelium | 234 post-neoadjuvant | ICI / pathologic response |
| [GSE280232](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280232) | Sorted T cells from tumor and adjacent lung | 13 patients (KRAS/STK11) | Neoadjuvant nivolumab + ipilimumab |

GSE280232 is the largest unused neoadjuvant ICI T-cell series in this window. It has no tumor transcriptome, so CLDN4 in malignant cells is not measurable.

## Not used as the run

Spatial, cell-line, blood-only, mouse, and bulk series in the same GEO hit list were skipped. GSE311609 remains a Xenium panel without CLDN4 (already noted in the handoff).
