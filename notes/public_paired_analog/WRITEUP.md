# Extra public analog: paired pre/post ICI lung RNA for TACSTD2 and CLDN4

Additive methods and results only. The private Zhejiang 25-pair IHC series (User A7: TROP2 H-score 94→121; NMPR 147 vs MPR 88) is taken as given and is not re-scored here. No public paired pre/post ICI TROP2 IHC matrix was found; the tables below are RNA (or array) analogs.

## Methods

**Search.** GEO (NCBI E-utilities, 2023–2026 publication dates) and ArrayExpress/BioStudies were queried for human lung/NSCLC series with immunotherapy or PD-1/PD-L1 terms plus paired, pre/post, or neoadjuvant language (105 GEO hits reviewed). Named leftovers were opened directly: GSE253564, GSE248378, GSE166449, GSE207422, GSE248249, GSE291670, GSE179994, GSE202417, GSE241934, GSE243013. Controlled-access EGA/GSA records were inventoried but not downloaded.

**What counted as a pair.** Same numeric patient ID, one pretreatment tumor and one on/post-ICI (or chemo-IO) tumor, with a downloadable expression matrix containing `TACSTD2` and `CLDN4`. Blood, T-cell-only, or DAC-locked series were excluded from scoring.

**GSE253564 leftover timepoints.** GSE253564 is the Altorki NCT02904954 pretreatment FPKM matrix (neoadjuvant durvalumab ± SBRT; *Nat Commun* 2023 / *Cancer Res Commun* 2024), not the Hu *Genome Med* 2023 leftover. The matching post-resection FPKM matrix is GSE248378 from the same trial. Patient numbers (`durva001` ↔ `pod01` / `36-M-PO`) join 20 patients. Major pathologic response (MPR) and arm were taken from the paper source data (Figure 2e “Major”; Figure 1b ARM). Tests: two-sided Wilcoxon signed-rank on post−pre FPKM; two-sided Mann–Whitney U for unpaired MPR vs NMPR. No log transform for the paired FPKM test (primary); log2((x+0.1)/(y+0.1)) is in the pair table as a sensitivity column.

**Other scored public sets.** GSE248249 (Memon *Cancer Cell* 2024): Affymetrix Clariom D RMA; probes `TC0100014340.hg.1` (`TACSTD2`) and `TC0700007993.hg.1` (`CLDN4`); 13 same-patient pre vs acquired-resistance pairs. GSE207422 bulk: deposited log2TPM, pretreatment only, MPR includes pCR. GSE207422 scRNA: sample-level mean log1p(UMI) in `EPCAM>0` cells; 3 pretreatment vs 12 post-treatment patients (0 pairs). GSE166449: deposited TPM as log2(TPM+1); 7 RECIST responders vs 15 non-responders; pretreatment only.

**Not done.** Re-alignment of SRA, purity adjustment, cutpoint search, or pooling across platforms. Zhejiang IHC was not used.

## Results

### Public paired lung ICI RNA exists; it is not TROP2 IHC

Two open same-patient lung-tumor ICI expression sets could be scored. Neither is TROP2 IHC. GSE253564+GSE248378 is neoadjuvant PD-L1 ± SBRT with MPR labels. GSE248249 is metastatic PD-(L)1 with acquired-resistance posts and no MPR.

**Table 1. Paired pre vs post `TACSTD2` / `CLDN4` (public analog).**

| Cohort | Gene | n pairs | Up / down | Median pre → post | Median Δ | p |
|---|---|---:|---:|---|---:|---:|
| Altorki joinable (all NMPR) | TACSTD2 | 20 | 14 / 6 | 47.8 → 45.4 FPKM | +16.2 | 0.11 |
| Altorki Arm1 durva alone | TACSTD2 | 14 | 8 / 6 | 63.0 → 44.3 | +2.0 | 0.71 |
| Altorki Arm2 durva+SBRT | TACSTD2 | 6 | 6 / 0 | 27.6 → 101.3 | +40.3 | 0.031 |
| Altorki joinable (all NMPR) | CLDN4 | 20 | 14 / 6 | 50.7 → 63.7 | +7.5 | 0.12 |
| Altorki Arm1 | CLDN4 | 14 | 8 / 6 | 50.7 → 48.5 | +2.5 | 0.81 |
| Altorki Arm2 | CLDN4 | 6 | 6 / 0 | 42.4 → 121.9 | +37.3 | 0.031 |
| GSE248249 all pairs | TACSTD2 | 13 | 6 / 7 | 4.84 → 4.75 RMA | −0.16 | 0.95 |
| GSE248249 same site | TACSTD2 | 4 | 2 / 2 | 4.43 → 4.08 | −0.02 | 1.00 |
| GSE248249 lung→lung | TACSTD2 | 1 | 1 / 0 | 5.81 → 6.41 | +0.60 | — |
| GSE248249 all pairs | CLDN4 | 13 | 2 / 11 | 5.66 → 4.83 | −0.89 | 0.017 |
| GSE248249 same site | CLDN4 | 4 | 0 / 4 | 5.20 → 4.74 | −0.74 | 0.13 |

GSE253564 leftover posts join only NMPR patients (0/20 MPR). The 11 pretreatment MPR biopsies in GSE253564 have no matching post column in GSE248378. The Arm2 n=6 rise (Wilcoxon p=0.031) is the smallest two-sided p possible when all six deltas are positive and is leveraged by one pair (durva007: TACSTD2 1.2→123 FPKM). Arm1 (durvalumab without SBRT) is flat. GSE248249 `TACSTD2` is flat; `CLDN4` falls in 11/13 pairs (p=0.017), but 12/13 pairs change anatomic site.

**Table 2. MPR / response analogs (unpaired or pretreatment).**

| Cohort | Contrast | Gene | n | Median hi vs lo | p |
|---|---|---|---:|---|---:|
| GSE253564 pre | MPR vs NMPR | TACSTD2 | 11 vs 21 | 157 vs 58 FPKM | 0.016 |
| GSE253564 pre | MPR vs NMPR | CLDN4 | 11 vs 21 | 62 vs 50 FPKM | 0.34 |
| GSE207422 bulk pre | MPR vs NMPR | TACSTD2 | 9 vs 15 | 5.00 vs 6.22 log2TPM | 0.34 |
| GSE207422 bulk pre | MPR vs NMPR | CLDN4 | 9 vs 15 | 5.85 vs 6.27 log2TPM | 0.26 |
| GSE207422 scRNA post | MPR vs NMPR (EPCAM+) | TACSTD2 | 4 vs 8 | 1.27 vs 1.41 | 0.57 |
| GSE207422 scRNA post | MPR vs NMPR (EPCAM+) | CLDN4 | 4 vs 8 | 1.44 vs 1.15 | 0.46 |
| GSE207422 scRNA | post vs pre (EPCAM+) | TACSTD2 | 12 vs 3 | 1.33 vs 1.87 | 0.36 |
| GSE207422 scRNA | post vs pre (EPCAM+) | CLDN4 | 12 vs 3 | 1.27 vs 1.25 | 0.73 |
| GSE166449 pre | RECIST R vs NR | TACSTD2 | 7 vs 15 | 2.12 vs 1.78 log2(TPM+1) | 0.41 |
| GSE166449 pre | RECIST R vs NR | CLDN4 | 7 vs 15 | 1.21 vs 1.33 | 0.95 |

GSE166449 is not paired (Lee *Cell* 2021 Samsung pembrolizumab; pretreatment titles only). Hu *Genome Med* leftover on-treatment / post-treatment bulk tumor RNA was not deposited: GSE207422 bulk is pretreatment-only (n=24); scRNA is 3 pre + 12 post from 15 different patients.

### Closest public sets that are not scorable pairs

ArrayExpress had no open 2023–2026 E-MTAB paired lung ICI RNA matrix. Closest leftovers: GSE91061 melanoma pre/on nivolumab (wrong tissue); GSE202417 paired PBMC CD8 (wrong compartment); GSE179994 T-cell/TCR (wrong compartment); GSE291670 post-only scRNA MPR 3 vs NMPR 3; GSE241934/HRA007419 and HRA006493 paired bulk/scRNA under GSA DAC; EGAD00001011302 (24 pre / 12 on nivo RNA) and NEOPREDICT-Lung EGAS00001007753 (17 NanoString pairs) under EGA DAC; a 2025 *Cellular Oncology* 13-pair neoadjuvant PD-1+chemo bulk series with no GEO accession.

## Reading

These are extra public RNA analogs, not a protein H-score replication set. The only open neoadjuvant lung ICI pairs with MPR labels (Altorki) are NMPR-only after the leftover join, and the only nominal pre→post rise is in the SBRT-containing arm (n=6). Pretreatment `TACSTD2` is higher in MPR than NMPR in GSE253564 (p=0.016) and not different in GSE207422 bulk (p=0.34). GSE248249 `TACSTD2` does not change at acquired resistance (p=0.95).

## Files

- `results/public_paired_analog/tables/paper_table_paired_pre_post.csv`
- `results/public_paired_analog/tables/paper_table_mpr_or_response.csv`
- `results/public_paired_analog/tables/altorki_paired_pre_post.csv`
- `results/public_paired_analog/tables/gse248249_pairs.csv`
- `results/public_paired_analog/tables/public_inventory_2023_2026.csv`
- `scripts/public_paired_analog/analyze.py`
