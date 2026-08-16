# Mouse lung ICI scRNA pool (additive; TISMO taken as given)

Public **mouse** lung / lung-cancer **single-cell** RNA after PD-1 / PD-L1 / VISTA / ICB.
**Not TISMO bulk.** User A4 TISMO (49/64 Tacstd2 up after ICB, Wilcoxon p=5.8e-5) is **taken as given** and is not re-audited.

## Question

In public mouse-lung ICI scRNA, is **Tacstd2 / Cldn4** higher in epithelial/tumor cells after ICB vs control, and is it higher in epithelial/tumor than in T/NK when both compartments exist?

## Inclusion

- *Mus musculus*, lung or lung-cancer model (LLC, KP, KL, HKP1, SCLC GEMM, carcinogen LUAD, CMT167)
- ICB / PD-1 / PD-L1 / VISTA (CA170) in the study or in the libraries
- Open processed gene × cell table on GEO (txt or 10x MTX)

## Honest skips / flags

| Series | Why |
|---|---|
| TISMO / GSE155972 | Bulk; taken as given |
| GSE267557 | CD45+ only (young vs aged, **both aPD-1**); epithelial **ABSENT**; RAW 1.2 GB |
| GSE283827 | Autochthonous SCLC ± aPD-1 ± ERBB2i; MTX deposited **without features.tsv** (32,589 genes). Cannot map Tacstd2/Cldn4 honestly |
| GSE275877 | Leftover LKR13 ICI R vs acquired-NR; unfiltered MTX is 31,053 × **33.97M** barcodes — not scored |
| GSE303943 | CMT167R **subcutaneous**; PKCi vs solvent, **not** ICB vs control |
| GSE157881/882 | Study is RT ± ICB; **these libraries are RT / club-cell depletion, not ICB** |
| GSE176091 | CD3+/CD45+ TIL sorts — epithelial **ABSENT** |

## Reproduce

```bash
pip install -r methods/mouse_scrna_ici_pool/requirements.txt
python3 methods/mouse_scrna_ici_pool/02_download.py
python3 methods/mouse_scrna_ici_pool/03_score.py
python3 methods/mouse_scrna_ici_pool/04_summarize.py
```

Large 10x files stay under `notes/mouse_scrna_ici_pool/raw/data/` (gitignored).

## Statistics

- Per sample: mean log1p (or deposited log-like) Tacstd2/Cldn4 in epithelial/tumor vs T/NK
- Epithelial/tumor = author CD45− sort, or Epcam+ or (Cdh1+ and Krt8+) or Ascl1/Chga/Insm1+. Do **not** use Sftpc/Scgb1a1 (ambient / normal AT2-club). Do not require Ptprc==0 (ambient CD45).
- T/NK = CD3 sort, or Cd3d/e, Cd8a, Nkg7, Ncr1; Epcam/tumor markers win over ambient T UMIs
- ICB vs control: Welch + MWU only if **n≥2 samples per arm**; n=1 vs 1 reports delta only
- Combined direction table: sign of epithelial Tacstd2/Cldn4 on ICB-vs-control contrasts that exist
