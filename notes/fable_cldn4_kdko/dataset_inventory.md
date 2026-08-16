# CLDN4 knockdown/knockout dataset inventory (exhaustive)

Search executed 2026-08-16 via NCBI E-utilities (`esearch`/`esummary` on the
`gds` and `sra` databases) plus per-series verification through the GEO
`acc.cgi` text interface. Scripts: `scripts/fable_cldn4_kdko/`.

Search terms combined the gene (CLDN4 / claudin-4 / "claudin 4") with the
perturbation vocabulary (knockdown, knock-down, knockout, knock-out, shRNA,
siRNA, CRISPR, sgRNA, silencing, silenced, depletion, depleted, deficient,
loss of function, gene editing, ablation) restricted to GSE series, plus
looser title/description and lung/pulmonary nets, and an SRA sample-level
sweep. Raw hit list: `results/fable_cldn4_kdko/geo_candidates_raw.json`,
`search_broaden.json`; per-sample verification: `geo_verified.json`.

## INCLUDED — true CLDN4 loss-of-function, processed data available

| Accession | Organism | Tissue / model | Perturbation | Platform | Processed file (md5) | Notes |
|---|---|---|---|---|---|---|
| **GSE207704** | Human | Breast cancer lines T47D & MCF7 | CRISPR **CLDN4-/-** vs WT (2 reps each) | RNA-seq (Cufflinks FPKM) | `GSE207704_CLDN4_RNAseq.txt.gz` (d3097646…) | Primary. Per-group FPKM (no per-replicate values). |
| **GSE50927** | Mouse | Whole **lung** | genetic **Cldn4 KO** vs WT | RNA-seq (author edgeR table) | `GSE50927_Cldn4lungWTvsKOgenes.csv.gz` (28b682e5…) | Lung-first. Baseline WT-vs-KO contrast; low biological replication. |
| **GSE22493** | Human | Ovarian cancer SKOV-3-IP-Luc | **CLDN4 siRNA knockdown** vs CLDN4-overexpression control | 2-colour Operon oligo array (GPL10555), 3 arrays | `GSE22493_RAW.tar` (9bae4611…) | Low confidence: control channel is CLDN4-overexpressing (not scramble); CLDN4 probe direction inconsistent across arrays (see `GSE22493_CLDN4_kd_diagnostic.csv`); TACSTD2 not on array. |

CLDN4 orientation controls (sanity checks that the perturbation is real):
- GSE207704: CLDN4 FPKM WT→KO drops MCF7 86.6→52.3, T47D 43.4→20.6.
- GSE50927: `Cldn4` logFC = -6.06, FDR = 4.1e-26 (defines logFC as KO-vs-WT).
- GSE22493: CLDN4 log2(KD/ctrl) = -1.64 / +1.22 / NA across the 3 arrays (inconsistent).

## EXCLUDED — CLDN4 mentioned but NOT a CLDN4 perturbation

| Accession | Why excluded |
|---|---|
| GSE99415 / GSE99416 / GSE99417 | "CLDN4 ceRNA network" study, but samples are gastric tumour vs adjacent normal **tissue** — no CLDN4 KD/KO. |
| GSE60885 | Placenta DNA-methylation; CLDN4 only correlative. |
| GSE65107 | Lymphocytic-colitis biopsies; CLDN4/5/8 observational, no perturbation. |
| GSE26055 | **CLDN7** siRNA (wrong gene). |
| GSE48443 / GDS4961 | **CLDN18** KO lung (wrong gene). |
| GSE84742 | Colonic crypt claudin differentiation; WT surface vs crypt cells, no perturbation. |
| GSE11759 | HNF4α mutant colon; CLDN4 downstream only. |
| GSE108417, GSE244820, GSE247130, GSE247271, GSE264098, GSE268909, GSE309751, GSE309894, GSE310539 | Lung studies where CLDN4 appears only as a marker/AT2 gene; GSE268909 is a CFTR/CF model. No CLDN4 perturbation arm. |
| SRP263109 / PRJNA634686 | 4C-seq of chromatin conformation at the CLDN4 locus in an MSH2-KO line; not a CLDN4 perturbation, and >processed-expression scope. |

All included processed files are < 14 MB; nothing exceeded the 2 GB skip
threshold, so no dataset was dropped for size.
