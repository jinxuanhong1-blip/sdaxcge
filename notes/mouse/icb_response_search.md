# Search for public mouse *lung* ICI RNA with ICB response labels

Question asked: does any public mouse lung ICI RNA matrix have **per-mouse ICB responder vs non-responder** labels for Tacstd2 / Cldn4?

## Explicit exclusion: GSE76628 is not lung ICI

Verified from GEO SOFT (`notes/mouse/raw_meta/GSE76628.soft.txt`):

- Title: *Stromal-Based Signatures for the Classification of Gastric Cancer [part II]*
- Organism: Mus musculus
- Design: Ad-VEGF-A164 angiogenesis model; **angiogenic sites in flanks**; anti-VEGFR antibodies (DC101 / G6)
- Assay: microarray (GPL1261), 78 samples, CEL only
- Not lung tumor, not anti-PD-1/PD-L1 ICB, not ICB response

**Do not use GSE76628 for mouse lung ICI Tacstd2/Cldn4.**

## Original 11 accessions (no per-mouse ICB R vs NR)

| Accession | What it actually labels | Per-mouse ICB R/NR? |
|---|---|---|
| E-MTAB-13704 | treatment arm (vehicle / aPD-L1 / combos) | No |
| GSE239485 | treatment arm (vehicle / polyI:C+aPD-1 / +anti-C5aR1) | No |
| GSE297630 / GSE297632 | control vs anti-PD-1 *tolerant surviving cells* | No (not R vs NR) |
| GSE197260 | treatment/time (n=1/arm) | No |
| GSE330658 | PTX / anti-VEGF (anti-PD-L1 arm not deposited) | No |
| GSE129297 / GSE133604 / GSE222158 | 1 scRNA sample per treatment | No |
| GSE241978 | AhR-KO DE table | No |
| PXD059688 | raw MS >2GB; no usable protein table | No |

Treatment ≠ ICB response. "Tolerant surviving LLC cells" ≠ RECIST-style responder.

## GEO search for mouse ICB responder / non-responder RNA

Query (2026-08-16): `GSE[Entry Type] AND Mus musculus AND (responder OR non-responder OR nonresponder) AND (PD-1 OR PD-L1 OR checkpoint)`.

Public **R vs NR RNA exists**, but the lung-relevant hits are **not lung tumors**:

| Accession | Model / tissue | Why not used |
|---|---|---|
| GSE139475 | CT26 colon, IgG / aPD-1 responder / non-responder | colon, not lung |
| GSE214348 / GSE227404 | SCC, aPD-L1 R vs NR | skin SCC, not lung |
| GSE274421 | MC38 colorectal in CC mice, R vs NR lines | colon, not lung |
| GSE166309 | ICB R vs NR scRNA (MDSC/PIM1 paper) | not lung tumor RNA |
| GSE293819 | breast ICB systemic immune | breast |
| GSE299686 | HNSCC liquid biopsy | not lung tumor |
| GSE255484 | melanoma aCTLA-4 R vs NR | melanoma |

No GEO series in this search is **mouse lung tumor RNA with per-mouse ICB R vs NR**.

## Closest public *lung* ICI RNA added this round

| Accession | What it is | What it is not | Processed file |
|---|---|---|---|
| GSE309199 | RPM SCLC *lung* tumors, Ctrl / aPD-1 / entinostat / combo, n=3 | not per-mouse R vs NR | `GSE309199_Mouse_Azam.TPMcalculator.raw_counts.tsv.gz` (0.5 MB) |
| GSE330941 | LLC WT (ICI-refractory model) vs Ago2KO (ICI-sensitized), n=4, **day 12, no ICI on these RNA samples** | not on-treatment R vs NR | `GSE330941_filtered_tablecounts_tpm.csv.gz` (0.7 MB) |
| GSE261890 | spatial RNA, immunotherapy-sensitive vs resistant NSCLC | spatial only; RAW.tar 1.8 GB; rds 730 MB; n=2 slides | cataloged, not analyzed (spatial, not a bulk response matrix) |

TISMO (http://tismo.cistrome.org → https://tismo.pku-genomics.org) annotates ICB response for syngeneic models, but the public download SPA did not expose a usable processed matrix in this run (`/datadownload` 404; GitHub `zexian/TISMO_data` is processing scripts only). TISMO was **not** used as a fabricated source.

## Honest answer

**There is no public mouse lung ICI RNA matrix with per-mouse ICB responder vs non-responder labels that we could download as a processed file <2 GB.**

Tacstd2 / Cldn4 can be tested against:

1. **ICI treatment vs control** (replicated: E-MTAB-13704, GSE239485, GSE297630, plus GSE309199) — this is *not* ICB response.
2. **Model-level ICI sensitivity** (GSE330941 Ago2KO-sensitized vs WT-refractory LLC; RNA taken without ICI) — this is a proxy, not R vs NR.
3. **Immune-score correlation** as a transcriptomic stand-in for a "hot" TME — not clinical response.
