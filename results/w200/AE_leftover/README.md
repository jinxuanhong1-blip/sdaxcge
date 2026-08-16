# ArrayExpress leftover: lung ICI / TACSTD2 / CLDN4

Honest bottom line: **there is no leftover ArrayExpress lung ICI cohort with processed expression and ICI response labels.** The first-wave hunt (PR #3) already took the real ICI hits. What remains is one lung + immune study that is **not ICI** (TIL–organoid co-culture), two metadata-only immuno-adjacent studies, a too-large NSCLC atlas, and the requested skip of E-MTAB-13530.

## What was already covered (do not re-claim)

| Accession | First-wave call | Why it is not leftover |
|---|---|---|
| E-MTAB-13704 | include | GEMM in-situ lung aPD-L1 combinations; processed counts |
| E-MTAB-15883 | include | lung-cancer-model anti-PD-1; processed data |
| E-MTAB-9451 | include_context | NSCLC ICOS/Treg profiling; no ICI-treated samples |
| E-MTAB-10633 | metadata_only | aPD-L1 / TGF-β trap; no processed files |
| E-MTAB-8867 | metadata_only | ICI-myocarditis, 3 NSCLC patients; no processed files |
| E-MTAB-13708 | exclude | coronavirus-vector immunotherapy, not ICI |
| E-MTAB-10027, 13770, 3218, 3732, 16855 | exclude | not lung ICI |

E-MTAB-13704 and E-MTAB-10633 did **not** re-hit this leftover synonym search (`aPDL1` / TGF-β wording). They are still first-wave covered, not new leftovers.

## Requested skip

**E-MTAB-13530** (10x Visium of human NSCLC lesions and non-involved lung) was skipped as requested. It is not ICI-labeled. Processed Visium files exist and are <2 GB each; a parallel Visium leftover slice owns that atlas.

Companion **E-MTAB-13526** (scRNA atlas, ~900k cells) is leftover only in the atlas sense. Annotated h5ads are 45–58 GB (over the <2 GB rule). Not ICI-treated. Not downloaded.

## Leftover that was actually analyzed

**E-MTAB-15784** — NSCLC patient-derived organoids + autologous TILs + co-culture (Chen / Peking University Shenzhen Hospital). TIL therapy platform, **no PD-1/PD-L1/CTLA-4 arm, no ICI response labels**.

14 processed MTX libraries, 126,201 cells, all files <2 GB. TACSTD2 and CLDN4 are in the features tables (except TACSTD2 missing from LCP90_TIL).

Mean detection rate (raw UMI, no cell-type annotation):

| sample type | n libraries | TACSTD2 | CLDN4 | EPCAM | CD3D |
|---|---:|---:|---:|---:|---:|
| Tumor organoid | 2 | 95.8% | 88.9% | 64.0% | 0.8% |
| Normal organoid | 2 | 98.3% | 88.3% | 63.0% | 1.2% |
| TILs | 5 | 0.2% | 0.4% | 0.3% | 96.7% |
| Co-culture | 5 | 40.5% | 43.5% | 4.4% | 96.1% |

Per-library numbers and caveats are in `E-MTAB-15784_summary.json` and `E-MTAB-15784_key_genes.tsv`.

Do not over-read the co-culture means: LCP94_CO is TIL-dominated (TACSTD2 0.2%). LCP90 tumor-organoid SDRF reuses LCP89_T_ORG files (depositor error; no unique LCP90 tumor MTX).

## Leftover metadata only (no matrix)

- **E-MTAB-12508** — KP lung + FLT3L/αCD40. Paper title: *Type 1 DCs promote immunity to neoantigens in mutated lung cancer refractory to checkpoint inhibitors*. Experiment is DC therapy, not ICI. No processed files.
- **E-MTAB-13710** — human NSCLC FRC / TLS scRNA companion to excluded E-MTAB-13708. No ICI treatment. No processed files.

## TROP2 / CLDN4 leftovers that are not lung ICI

ArrayExpress TACSTD2/TROP2 hits are intestine/CRC (E-MTAB-11377, 11382, 11466, 16433, 16843, 16849). No lung TROP2 study.

CLDN4 hits are GEO mirrors (E-GEOD-22493 ovarian; E-GEOD-50927 ventilator lung injury). Excluded as GEO mirrors.

`lung AND TROP2` and `lung AND atezolizumab` returned **0** ArrayExpress hits.

## Reproduce

```bash
python3 scripts/ae_leftover/search_leftover.py
python3 scripts/ae_leftover/download_leftover.py
python3 scripts/ae_leftover/analyze_15784.py
python3 scripts/ae_leftover/build_catalog.py
```

Processed MTX files are not in git (`downloads/` is gitignored). SHA256s are in `download_manifest.json`.
