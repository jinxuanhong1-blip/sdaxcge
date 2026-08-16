# Extra public mouse-lung ICI RNA beyond TISMO

Additive slice for claim A4. **TISMO is taken as given** (49/64 models Tacstd2 up after ICB, p=5.8e-5) and is **not** re-audited. GSE155972 (TISMO LLC ICB) is excluded.

## Question

Do **public** mouse lung / NSCLC RNA series **outside TISMO** (LLC, KP, KL, orthotopic NSCLC; PD-1 / PD-L1 / CTLA4) show Tacstd2 or Cldn4 changing **pre vs post ICB**, or vs **response** if labeled?

GSE179500 (LKB1/KL; Tacstd2 higher when LKB1 is off) is cited as a known genotype series and is not treated as ICB.

## Inclusion

- Organism: *Mus musculus*
- Tissue / model: LLC / LL/2, KP (`Kras`+`Trp53`), KL (`Kras`+`Lkb1`/`Stk11`), 344SQ, CMT-167, HKP1, or orthotopic / GEMM NSCLC
- Treatment text: PD-1, PD-L1, CTLA4, or ICB (combo arms kept but flagged)
- Open processed gene table on GEO, compact enough to download (< ~20 MB)

## Exclusion

- TISMO source GSE155972
- SCLC GEMM ICB (scored only as out-of-scope leftovers)
- Human xenografts (e.g. NCI-H460)
- Immune-only sorts without a tumor-epithelial matrix (CD8/Treg/neutrophil/BMDM-only)
- CRISPR sgRNA counts, ChIP-seq, WGBS, TCR-seq
- scRNA 10x matrices / RDS (no compact gene × sample table)

## Search (2026-08-16)

NCBI E-utilities `db=gds`, six queries (LLC, KP, KL, orthotopic, plus a broad mouse-lung PD-1/PD-L1/CTLA4 expression search). Union = **283** GSE summaries (`notes/a4_extra_mouse_lung_ici/raw/geo_search.json`). Sixty-four candidates were page-fetched; compact matrices were downloaded for scoring (`02_download_matrices.py`).

ArrayExpress leftover noted: **E-MTAB-15883** is PD-1 fate-mapped *immune* cells from subcutaneous LLC + CTX ± aPD-1 — not tumor Tacstd2/Cldn4.

## Gene IDs

| Gene | Symbol | Ensembl (mouse) |
|---|---|---|
| Tacstd2 | Tacstd2 / Trop2 | **ENSMUSG00000051397** |
| Cldn4 | Cldn4 | **ENSMUSG00000047501** |

`ENSMUSG00000030798` is **Cd37** and was not used.

## Statistics

- Two-sided Welch *t* and Mann–Whitney U when both arms have n≥2
- log2FC = mean difference if the deposited table is already log-like (VST, voom, Clariom log2, median-centered log); otherwise log2((mean_b+1)/(mean_a+1))
- n=1 arms: values reported, *p* not computed
- No R vs NR meta-analysis: that label is absent

## Reproduce

```bash
python3 methods/a4_extra_mouse_lung_ici/01_fetch_geo_suppl.py
python3 methods/a4_extra_mouse_lung_ici/02_download_matrices.py
python3 methods/a4_extra_mouse_lung_ici/03_score_tacstd2_cldn4.py
python3 methods/a4_extra_mouse_lung_ici/04_summarize.py
```

Requires `pandas` is not required; `numpy`, `scipy`, `openpyxl`, `matplotlib`.
