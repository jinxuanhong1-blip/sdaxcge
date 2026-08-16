# Hunt: public neoadjuvant lung IO scRNA/bulk for TACSTD2

Question: which **public** lung neoadjuvant IO datasets (scRNA or bulk) can recover

1. **Finding A** — malignant `TACSTD2` higher in NMPR than MPR
2. **Finding B** — per-patient `TACSTD2` vs T/NK Spearman **ρ < 0**

Seeds specified: `GSE207422`, `GSE243013`, `GSE146100`, `GSE229353`. Additional series from GEO E-utilities, 2023–2026. **No invented accessions.** Public processed matrices only (no FASTQ).

## Honest answer

**Single-cell malignant TACSTD2 does not recover either claim at p < 0.05.**

| Dataset | Assay | Finding A (NMPR > MPR) | Finding B (ρ vs T/NK) |
|---|---|---|---|
| GSE207422 | bulk log2TPM, pre-tx, n=24 | direction only (median 6.22 vs 5.00; MW p=0.34; 15 vs 9) | **ρ=−0.506, p=0.012** |
| GSE207422 | scRNA, marker Epi, n=14/15 | direction only (1.35 vs 1.21; p=0.64; 10 vs 4) | ρ=−0.264, p=0.34 |
| GSE241934 | scRNA, author Epi, n=35 | **no** (1.50 vs 1.67; p=0.96; 25 vs 10) | **null** (ρ=−0.014, p=0.94) |
| GSE291670 | scRNA, marker Epi, n=6 | **no** (opposite; 3 vs 3; p=0.40) | ρ=−1.0 on n=4 after filter; not interpretable |
| GSE248378 | bulk FPKM, post-tx, n=29 | **not testable** (no MPR on GEO) | **ρ=−0.781, p=5.8×10⁻⁷** |
| GSE329813 | GeoMx tumor-bed, n=22 | **yes, p=0.002** (4.80 vs 3.94; 11 vs 11) | **ρ=−0.512, p=0.015** |

GSE329813 is **not scRNA**. It is GeoMx DSP of primary-tumor-bed ROIs (mixture). It is the only public matrix here that recovers **both** directions at p<0.05, and it does **not** isolate malignant cells.

The best-powered scRNA test with author epithelial labels and public MPR is **GSE241934** (NEOTIDE + real-world, 45 patients). It does not recover A or B.

## Seeds

**GSE207422** (Hu et al., PMID 36869384; public 2023). Only seed that can test both findings on tumor cells.

- Bulk: 24 pre-treatment biopsies, 1/patient, NMPR 15 / MPR+pCR 9. TACSTD2 present.
- scRNA: 92,330 cells, 15 patients (12 post-tx surgery, 3 pre-tx). Marker rule: epithelial = (EPCAM/KRT19/KRT18/KRT8 > 0) AND PTPRC==0 (12,018 cells); T/NK = (CD3D/CD3E/NKG7/GNLY > 0) AND PTPRC>0 (41,841 cells). Overlap 0.

**GSE243013** (Liu et al., PMID 40147443; public 2025). 243 patients, public MPR (non-MPR 112, MPR 45, pCR 85). Processed immune MTX is 7.1 GB. Author `major_cell_type` is only T/NK, B, Myeloid (1,254,749 cells). **No epithelial/malignant compartment.** Cannot test malignant TACSTD2. Matrix not downloaded.

**GSE146100** (PMID 33820821; public 2021). One multiple-primary LUAD patient, 3 nodules, pembrolizumab. Per-nodule response, not MPR/NMPR. n=1 patient. Not analyzed.

**GSE229353** (PMID 37231145; public 2023). 7 CD45+ 10x MTX in `GSE229353_RAW.tar`. GEO labels are Chemo vs anti-PD1+Chemo, not MPR. Cannot test A. Immune-sorted, cannot test malignant TACSTD2.

## 2023–2026 GEO search

Queries and raw hits:

- `metadata/geo_search_gse.tsv` — title-restricted (16 GSE)
- `metadata/geo_search_broader.tsv` — broader GDS (61 GSE)

Additional series that were opened on GEO (SOFT + supplementary listing), not invented:

| Accession | Why opened | Why it can/cannot test the claims |
|---|---|---|
| GSE241934 | neoadjuvant sintilimab+chemo scRNA, processed MTX | **Yes, tested.** Author Epi + MPR. Does not recover A or B. |
| GSE248378 | neoadjuvant durvalumab bulk FPKM | Finding B only. No public MPR. |
| GSE291670 | 6-patient camrelizumab+anlotinib scRNA MTX | Tested. n=3 vs 3. A opposite. B n=4 unusable. |
| GSE329813 | GeoMx neoadjuvant pembro+chemo, MPR in titles | Tested as spatial-bulk. Recovers A and B. Not scRNA. |
| GSE225620 | neoadjuvant PD-1 blood RNA | Blood. No malignant TACSTD2. |
| GSE260770 | sintilimab GGO, exosomal FPKM | Blood exosome. |
| GSE280232 | KRAS neoadjuvant nivo+ipi 10x MTX | **Sorted T cells.** No malignant cells. No MPR on GEO. |
| GSE337519 | “sc landscape of neoadjuvant chemo-IO” | n=1 patient. No group stats. |
| GSE274934 | ALK+ IO sc multiome | 15 GB RAW.tar; not standard neoadjuvant chemo-IO. Not pulled. |
| GSE299111 | lung Tregs after chemo | Chemo only, sorted immune. |
| GSE300685 | TIL-NK in resectable smokers | NK-focused; not used. |
| GSE249568 / GSE250509 | neoadjuvant tepotinib GeoMx | MET TKI, not IO. n=1. |
| GSE240033 | canine osteosarcoma | Not human lung (false positive). |

## Methods (so the numbers can be re-run)

Scripts in `scripts/`. Python 3.12, packages in repo `requirements.txt`.

- Finding A: Mann–Whitney U on per-sample (or per-patient) malignant/epithelial mean `log1p(TACSTD2)` (scRNA) or bulk/GeoMx TACSTD2. pCR counted as MPR. Two-sided p and one-sided NMPR>MPR both stored in the JSON files.
- Finding B: Spearman (and Pearson) of that TACSTD2 value vs T/NK fraction (scRNA) or mean T/NK marker expression (bulk/GeoMx). One row per patient when multiple samples exist (post-tx preferred in GSE207422).
- scRNA marker rule is written in `analyze_scrna_GSE207422.py` and reused in `analyze_scrna_GSE291670.py`. GSE241934 uses the authors’ `major.cell.type == Epi` / `T` / `NK`.
- Minimum 20 epithelial cells (and 20 T/NK cells for B) unless a sensitivity table says otherwise.

GEO download URLs are the NCBI FTP supplementary paths recorded in each series SOFT file. Large matrices were kept local under `data/` and are **not** in git.

## What this hunt does *not* claim

- It does not claim that TACSTD2 is unrelated to neoadjuvant IO. Bulk and GeoMx mixtures show a negative TACSTD2–T/NK correlation and, in GeoMx tumor-bed, higher TACSTD2 in NMPR.
- It does not claim a confirmed **malignant-cell** NMPR>MPR effect. The two scRNA datasets that can isolate epithelium either lack power (GSE207422) or go the other way / null (GSE241934, GSE291670).
- GSE243013 cannot rescue Finding A no matter how large it is: there are no malignant cells in the public matrix.

## Outputs

- `tables/inventory.tsv` — every accession considered
- `tables/summary.json` — numeric digest
- `tables/GSE*_results.json` and `*_per_sample.csv` — per-dataset stats
- `metadata/` — SOFT, GSM titles, E-utilities search dumps
