# Extra accession hunt (2026-08-16)

Queried NCBI GEO for mouse lung + anti-PD-1/PD-L1 RNA series. Processed files <2 GB only.

## Included in this extension
- GSE114601 GEMM NSCLC Vehicle vs anti-PD1 (n=2)
- GSE157880 KP lung IgG vs PD-1 ± RT (0 Gy pair is clean ICB)
- GSE309199 RPM SCLC Ctrl vs aPD1 (n=3)
- GSE169196 KPM total-tumour IgG vs A2V+aPD1 (no aPD1-mono arm)

## Explicitly skipped (honest)
- GSE190264: in-vitro LLC1 chemo/MEK, no ICB antibody
- GSE114300: sorted CD4 T cells, not tumour epithelium
- GSE277610: sorted CD8 T cells
- GSE309192: control vs tumour tissue, no ICB
- GSE184000: whole-lung irAE, not tumour
- GSE315010: RAS(ON) inhibitors, no ICB antibody arm in the mouse matrix
- GSE303940: CMT167R PKC inhibitor, n=1, no ICB
- GSE194166: scRNA/TCR, 238 MB, n=1/arm
- GSE129298: scRNA n=1/arm (companion of already-included GSE129297)
- GSE157883: parent SuperSeries of GSE157880 (used the processed Bulk048)

## TISMO re-try (2026-08-16, second pass)
Probed `https://tismo.pku-genomics.org/` (SPA), frontend routes `/datadownload` `/gene` `/internaldata`, and backends `/tismo` `/rtismo` — all API paths returned 404. GitHub `zexian/TISMO_data` is processing scripts only. **49/64 Tacstd2-up after ICB was not independently recomputed.**

## Radiation confounder (keep 0 Gy as primary)
GSE157880 Tacstd2 PD-1 4Gy vs IgG 0Gy = +1.42 (p=0.003) but IgG 4Gy vs IgG 0Gy = +1.78 (p=0.013). Rise is radiation, not ICB.
