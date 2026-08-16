# Hunt log — public paired pre/post ICI lung RNA (TACSTD2 / CLDN4)

Goal: find **same-patient** pre- vs post-ICI **lung tumor** RNA (or IHC) that can test
whether TACSTD2/TROP2 **rises after IO**, as claimed from TISMO and the Zhejiang IHC
H-score 94 → 121 (claim **A7**).

Search date: 2026-08-16. Sources: GEO FTP, OmicsDI, PubMed/PMC, TISMO portal,
sibling PR #102 (`cursor/claim-a6-a11-a7-8834`).

## Verdict (A7 analog)

| Claim | Public analog | Direction match? | Stats |
|---|---|---|---|
| A7: TROP2 **rises** after ICI (Zhejiang IHC 94→121; TISMO) | **GSE248249** — 13 same-patient NSCLC pairs, PD-(L)1, pre vs **acquired resistance** | **No.** TACSTD2 is flat. | Wilcoxon p=0.95; median Δlog2=−0.16; 6↑ / 7↓ |
| A7 (CLDN4, same pairs) | GSE248249 | **Falls** (not hypothesized) | Wilcoxon p=0.017; median Δlog2=−0.89; 2↑ / 11↓ |
| A7 mouse TISMO Tacstd2 rise | TISMO portal has no programmatic dump; closest public mouse **lung** ICB RNA = **GSE246922** (KP, LLC1) | **No.** Tacstd2 slightly down, n.s. | KP ICB-res vs parental log2FC=−0.20, p=0.21 |
| A7 public paired lung TROP2 **IHC** | None with reusable numeric pairs | — | Inoue 2025 mixed treatments (not ICI-only); JCO 8591 abstract not deposited |

**A7 analog claim: NOT SUPPORTED** in public paired ICI lung RNA. The hypothesized rise
is not reproduced. GSE248249 TACSTD2 is flat across 13 pairs (p=0.95) and across the
4 same-site pairs (p=1.00). **Only 1 pair is lung-to-lung** (Patient 01, Δlog2 +0.60,
descriptive). CLDN4 falls in 11/13 pairs (p=0.017).

## Datasets inspected

| Accession | Tissue / assay | Timepoints | Same-patient pairs? | TACSTD2 measurable? | Used? | Why / why not |
|---|---|---|---|---|---|---|
| **GSE248249** | NSCLC tumor, Clariom D microarray | 13 pre + 29 post (acquired resistance) | **Yes, 13 pairs** | Yes (TC0100014340.hg.1) | **Yes — primary A7 analog** | Only public same-patient pre/post ICI *lung tumor* transcriptome |
| GSE246922 | Mouse KP / LLC1 lung, RNA-seq VST | parental / IFNγ / ICB-resistant | cell-line replicates, not patients | Yes (ENSMUSG00000051397) | Yes — mouse analog | Same paper; TISMO dump unavailable |
| GSE207422 scRNA | NSCLC tumor, BD Rhapsody | 3 pre + 12 post neoadjuvant PD-1+chemo | **0 pairs** (different patients) | Yes | Yes — POST axis | Cross-patient; TACSTD2 slightly down (p=0.29) |
| GSE207422 bulk | NSCLC, log2 TPM | 24 **pre only** | No | Yes | Yes — PRE vs response | No post RNA deposited |
| GSE126044 | NSCLC, counts | 16 **pre only** | No | Yes | Yes — PRE | No on/post |
| GSE135222 | NSCLC, TPM | 27 **pre only** | No | Yes (Ensembl) | Yes — PRE | No on/post |
| GSE91061 | **Melanoma**, FPKM | Pre + On nivolumab | Yes, 43 pairs | Yes (Entrez 4070) | Orthogonal only | Not lung. Responders ↑ TACSTD2 on-Rx (p=0.019) |
| GSE179994 | NSCLC, scRNA+TCR | pre / post PD-1±chemo | some site-matched | **No** (T cells only) | No | No epithelial transcriptome |
| GSE176021 | NSCLC neoadjuvant nivo, scRNA | TIL / blood | T-cell focused | No (T cells) | No | Cannot score TACSTD2 |
| GSE111414 | NSCLC **PBMC** CD8 RNA | pre + week 4 nivo | Yes (blood) | Not tumor TACSTD2 | No | Wrong compartment |
| GSE154286 | NSCLC, targeted 201-gene | pre / post **chemo** | 29 pairs | Panel may lack TACSTD2; not ICI | No | Chemotherapy, not ICI |
| GSE136961 | NSCLC, Oncomine 395-gene | **pre only** | No | Immune panel | No | Baseline only |
| GSE249000 | Mouse CT26 (colon) RNA | IFNγ / ICB-resistant | cell line | Yes (mouse) | No | Not lung |
| GSE145281 | ICI blood/tumor scRNA | mixed | not lung-paired bulk | T/myeloid focused | No | Not paired lung bulk TACSTD2 |
| POPLAR/OAK (EGAS00001005013) | NSCLC atezolizumab | **pre only** | No | Yes (controlled EGA) | No | DAC/controlled; pretreatment |
| Tempus CRC-24-0605 | NSCLC RNA | naive / post-ICI / 56 pairs | Yes, but **licensed** | Yes (TROP2 “maintained”) | No | Not public |
| TISMO (tismo.pku-genomics.org) | mouse syngeneic ICB | pre / ICB | many models | Yes (UI) | **Not downloadable** | Shiny UI only; GitHub `zexian/TISMO_data` is scripts, no matrices |
| Inoue 2025 (PMC11801068) | NSCLC TROP2 IHC | paired pre/post **mixed** Rx | some | Protein | No | 5 ICI patients: **no TROP2 change**; not ICI-only series |
| Zhejiang IHC (in-house) | NSCLC TROP2 IHC | paired ICI | yes (94→121) | Protein | **Out of scope** | Not public; this slice does not use it |

## What would have counted as a direction match
Same-patient lung tumor, ICI between biopsies, TACSTD2/TROP2 **up** with a
pre-specified paired test (Wilcoxon) at nominal p<0.05. None of the public
resources met that. GSE248249 is the correct design and the direction is null/down.

## Files produced from this hunt
- `results/fable_paired/tables/gse248249_paired_stats.csv`
- `results/fable_paired/tables/gse248249_patient_pairs.csv`
- `results/fable_paired/figures/fig5_gse248249_paired_pre_post.png`
- `results/fable_paired/tables/gse246922_mouse_stats.csv`
- `results/fable_paired/tables/a7_analog_verdict.json`
