# Rework: real public KRAS/LKB1 (STK11) lung datasets with Tacstd2

**This folder is self-contained.** It replaces any analysis that treated **GSE76628 as KL lung**. That accession is nude-mouse flank skin. The real open KL / KRAS+STK11 lung series that actually contain `Tacstd2`/`TACSTD2` are listed below and were re-downloaded from GEO for this run.

Script: `analysis/kl_real/run_kl_real.py`.  
Machine-readable summary: `metrics.json`. Checksums: `provenance.json`.  
Catalog (including series not re-run here): `real_KL_lung_catalog.tsv`.

---

## Verdict in one page

| Accession | Open? | Is it KL lung? | Tacstd2 / TACSTD2 | Honest n | What you may claim |
|---|---|---|---|---|---|
| **GSE76628** | yes | **NO** — nude mouse **skin** (Ad-VEGF stroma) | high, because skin epithelium expresses Trop2 | 78 arrays | nothing about lung or KL |
| **GSE180963** | yes | **YES** — K vs KL GEMM nodules | epithelial-restricted (KL epi 65%+ vs immune 4%+) | **n = 1 mouse / genotype** | tumor-vs-immune selectivity. **Not** a KL>K claim |
| **GSE179502** | yes | **YES** — KT;Lkb1XTR sorted tumor cells | 54.7% of cells Tacstd2+; ~14% high subset | **n = 3 vs 3 mice** | a TROP2-high KL tumor subset exists; LKB1-off means look higher (mouse MW *p* = 0.10) |
| **GSE280232** | yes | **YES** — human KRAS±STK11 after neoadjuvant ICB | ~0% in real sorted T-cell libs; 72–85% in EPCAM+ mixed-tumor cells | mixed-tumor epithelium usable in **2 STK11mut + ~2 STK11wt** patients | TACSTD2 is an epithelial gene here, not a TIL gene; no powered KL-vs-K human tumor-cell test |
| **GSE179500** (companion bulk, same XTR paper) | yes | **YES** | author DE Restored vs Non-Restored **log2FC = −0.46, padj = 0.032** | p53-WT: 10 Non-Restored, 7 Restored, 4 KT | **best-powered** Tacstd2 KL contrast in this rework: higher when LKB1 is off |

**Bottom line.** A TROP2-high state exists inside real KL lung tumor epithelium (GSE179502). Bulk XTR data with biological replicates say Tacstd2 is modestly higher when LKB1 is off (GSE179500). GSE180963 can only support “Trop2 is epithelial, not immune,” and even that sits on **n = 1**. GSE280232 is the right *human* KL series but is mostly sorted T cells, so it cannot test tumor-cell TACSTD2 except in a handful of mixed libraries.

No private mouse cohort was used.

---

## GSE76628 is not KL lung

GEO title: *Stromal-Based Signatures for the Classification of Gastric Cancer [part II]* (Uhlik et al., *Cancer Res* 2016, PMID 27197264; superseries GSE76630).

| Field | Value |
|---|---|
| Tissue | Athymic **nude (Foxn1nu)** female **flank skin** |
| Perturbation | Ad-VEGF-A164 tumor-*surrogate* angiogenesis / wound-stroma model |
| Assay | Affymetrix Mouse Genome 430 2.0 **bulk microarray**, 78 samples |
| Kras / Lkb1 / lung / tumor cells / T cells | **none** (nude mice have no mature T cells) |

Tacstd2 is high on that array because Trop2 is a normal skin-epithelial gene, not because anyone sequenced a KL lung tumor. Using GSE76628 as “public KL mice” is a mismatch.

---

## What was analyzed (all three requested series are open)

Processed GEO files were downloaded 2026-08-16 from NCBI FTP (see `provenance.json`). Normalization for scRNA: `log1p(CP10k)`; a cell is called positive if that value is > 0 (at least one UMI). **The inferential unit is the mouse or the patient.** Cell-level *p*-values are reported only as descriptive and labeled pseudoreplication.

### 1. GSE180963 — true KL GEMM, **n = 1**

*Single cell RNA sequencing of tumor sections from GEMM harboring KrasG12D/+ or KrasG12D/+Lkb1fl/fl (KL) mutation.* Public 2022. Two author-filtered 10x matrices: K = 6,696 cells, KL = 7,564 cells.

Epithelial call (clustering-free): `Epcam+` AND (`Krt8`/`Krt18`/`Krt19`/`Cldn18`)+ AND `Ptprc−`.

| Genotype | Mice | Epithelial cells | Tacstd2+ in epi | Tacstd2+ in immune | Immune fraction |
|---|---|---|---|---|---|
| K (KrasG12D) | **1** | 28 | 17.9% (5/28) | 2.3% | 82.8% |
| KL (KrasG12D;Lkb1fl/fl) | **1** | 170 | 64.7% (110/170) | 4.3% | 82.3% |

**What holds:** inside the KL sample, Tacstd2 is ~epithelial and near-absent from immune cells (mean log-norm 0.61 vs 0.04). Same direction in K, on 28 epithelial cells.

**What does not hold:**

1. **n = 1 mouse per genotype.** Genotype is fully confounded with the single 10x library. The cell-level KL-vs-K *p* = 3×10−5 is pseudoreplication. It is not evidence that LKB1 loss raises Trop2.
2. This dataset does **not** reproduce an immune-desert as a lower CD45 fraction (both samples ~82% `Ptprc+`). “Cold” here would have to be T-cell state / exclusion, which this pass did not test.
3. Tumor cells are rare after dissociation. K epithelium is 28 cells. Ambient `Sftpc` is a known problem in this series; we did not use surfactant genes to call epithelium.

Figures: `figures/GSE180963_Tacstd2_compartment_genotype.png`.  
Tables: `tables/GSE180963_*.csv`.

### 2. GSE179502 — true KL tumor cells, n = 3 vs 3

Murray et al., *Cancer Discov* 2022 (PMID 35228570). Sorted Lin− tdTomato+ neoplastic cells from KT;Lkb1XTR lung tumors. Six mice, barcodes prefixed with mouse ID (so sample mapping is not guessed):

| Mouse | Cohort | Treatment | n cells | % Epcam+ | % Ptprc+ | mean Stk11 | mean Tacstd2 | % Tacstd2+ |
|---|---|---|---|---|---|---|---|---|
| CM0875 | NonRestored | Tamoxifen (no FLPo) | 3011 | 73.9 | 0.10 | 0.005 | 0.647 | 57.1 |
| ZR1932 | NonRestored | Vehicle | 1103 | 82.3 | 0.18 | 0.013 | 0.867 | 71.4 |
| ZR1966 | NonRestored | Vehicle (FLPo, no tam) | 3652 | 75.4 | 1.40 | 0.003 | 0.678 | 59.9 |
| CM0879 | Restored | Tamoxifen + FLPo | 4663 | 63.4 | 0.30 | 0.254 | 0.432 | 40.7 |
| CM0884 | Restored | Tamoxifen + FLPo | 952 | 86.8 | 0.21 | 0.294 | 0.614 | 66.9 |
| ZR1969 | Restored | Tamoxifen + FLPo | 2636 | 77.4 | 1.37 | 0.274 | 0.547 | 58.1 |

Sanity check: `Stk11` is ~50× higher in Restored mice (0.25–0.29 vs 0.003–0.013). These really are LKB1-off vs LKB1-on tumor cells.

**TROP2-high subset (all 16,017 cells):** 54.7% Tacstd2+; cells at or above the 75th percentile of positive cells = **13.7%** of the dataset. That is a real high subset inside KL / Lkb1-mutant tumor epithelium.

**LKB1-off vs restored, mouse means:** NonRestored 0.731 vs Restored 0.531 (log2FC +0.46). Mann–Whitney on **n = 3 vs 3 mice: U = 9, p = 0.10**. That is the smallest two-sided *p* you can get with 3 vs 3 and no ties — every NonRestored mouse mean is ≥ every Restored mean except the CM0884 overlap with CM0875. **Underpowered. Directionally consistent with the bulk companion. Not a definitive LKB1→Trop2 causal claim from scRNA alone.**

Immune-low cannot be tested: immune cells were FACS-depleted (0.7% `Ptprc+`).

Figures: `figures/GSE179502_Tacstd2_by_cohort.png`, `GSE179502_Tacstd2_per_mouse.png`.

### 3. GSE280232 — true human KRAS/STK11, mostly the wrong compartment for TACSTD2

Rosner et al., *Clin Cancer Res* 2024 (PMID 39545922). Open. 5′ GEX + TCR from 6 KRAS<sup>mut</sup>/STK11<sup>mut</sup> and 7 KRAS<sup>mut</sup>/STK11<sup>wt</sup> patients after neoadjuvant ICB.

GEO characteristics say most libraries are **sorted CD3+ T cells**. TACSTD2 is an epithelial gene. In those T-cell libraries TACSTD2 positivity is **0.01–0.24%** (noise / ambient). That is a negative control, not a KL tumor-cell result.

**Exception / QC flag:** patient **01-130** is labeled “Sorted T cells” but is 18–31% `EPCAM+` and 24–35% `TACSTD2+`, with only ~20–32% `CD3E+`. Treat 01-130 as a contaminated / mis-sorted library, not as T cells.

**Mixed CD45+ and CD45− tumor libraries** (the only place tumor-cell TACSTD2 is testable):

| Patient | Genotype | EPCAM+ PTPRC− cells | TACSTD2+ in that epi gate | TACSTD2+ in PTPRC+ |
|---|---|---|---|---|
| 11318-141 | STK11mut | 2928 | **85.3%** | 12.6% |
| 11318-139 | STK11mut | 1546 | **85.1%** | 6.4% |
| 11318-145 | STK11wt | 1392 | **72.1%** | 4.9% |
| 11318-143 | STK11wt | **18** | 72.2% (n=18) | 1.5% |
| MD043-132 (2 libs) | STK11wt | 69 + 45 | 87% / 82% | <0.5% |

TACSTD2 is epithelial in human KRAS-mutant tumors too. STK11mut vs wt in that epi gate is 85% vs 72% on **two vs one** epithelium-rich tumors (11318-143 and MD043-132 recovered almost no epithelium). **Do not call that a KL effect.**

Figures: `figures/GSE280232_TACSTD2_pct_by_sample.png`, `GSE280232_mixed_TACSTD2_compartment.png`.

### 4. GSE179500 — companion bulk (same paper as GSE179502)

Not in the “analyze these three” list, but it is the only open KL lung series here with a **real biological-replicate DE call** for Tacstd2. Author file `RestoredvNon-Restored`:

- Tacstd2 log2FC = **−0.463** (higher in Non-Restored / LKB1-off)
- *p* = 0.0016, **padj = 0.032**

Restricted to p53-WT samples (this reanalysis):

| Cohort | n samples | mean DESeq2-norm Tacstd2 |
|---|---|---|
| Non-Restored (LKB1-off) | 10 | 319 |
| Restored (LKB1-on) | 7 | 220 |
| KT wild-type (Lkb1 intact) | 4 | 258 |

Non-Restored vs Restored, sample-level Mann–Whitney **p = 0.003**, log2FC +0.54. Non-Restored vs KT wild-type is **not** significant (*p* = 0.45, n=4 KT). So the cleanest supported statement is: **in this XTR system, restoring LKB1 lowers Tacstd2; a KL-vs-K (never-lost) difference is not shown.**

Figure: `figures/GSE179500_Tacstd2_bulk_p53WT.png`.

---

## n = 1 caveats (read before quoting any number)

- **GSE180963 is n = 1 vs 1.** Every KL-vs-K statistic in that series is one mouse against one mouse. Cell-level tests inflate *n* to thousands and are invalid for genotype inference.
- **GSE179502 is n = 3 vs 3.** That is better, still small. Mouse-level *p* = 0.10 is the floor of a 3-vs-3 two-sided rank test. The cell-level *p* ~ 10−84 in `GSE179502_Tacstd2_contrasts.csv` is **pseudoreplication** and must not be quoted as evidence.
- **GSE280232 mixed-tumor epithelium is n = 2 STK11mut patients** with enough EPCAM+ cells. The other “STK11wt mixed” libraries are 18–69 epithelial cells.
- **GSE179500** is the only contrast in this folder whose *p*-value is a sample-level test with n ≥ 7 per arm. Even there, KT wild-type is n = 4 and does not differ from Non-Restored.
- No clustering, doublet removal, or ambient RNA correction (SoupX/CellBender) was run. Positivity is UMI ≥ 1. Ambient RNA can create low-level Tacstd2 in immune cells (GSE180963 ~4%, GSE280232 mixed immune 0.3–13%).
- GSE180963 and GSE280232 mixed libraries are immune-dominated after dissociation; epithelial percentages are not tumor purity.

---

## Other real public KL lung sets (not re-run here)

See `real_KL_lung_catalog.tsv`. Worth knowing:

- **GSE179501** — same XTR mice as GSE179502 but total-viable (tumor + immune). Right design for “TROP2-high / immune-low,” but prior crude checks showed heavy ambient/doublet mixing (`Epcam+` often `Ptprc+`). Needs a proper QC pass.
- **GSE194166** — KP vs KPL **immune-sorted** scRNA. Immune-cold KL TME, no tumor epithelium, so no Tacstd2 tumor subset.
- **GSE274351, GSE175479, GSE277929** — KL bulk / LCM context, not single-cell.

Human bulk LUAD with KRAS/STK11 calls (TCGA, etc.) can test TACSTD2 vs genotype at n ≫ 10; that is outside this GEO-scRNA rework.

---

## Reproduce

```bash
pip install -r analysis/kl_real/requirements.txt
# downloads live in data/kl_real/ (gitignored). Re-run the script; it will
# reuse files already present and re-extract GSE280232 GEX matrices.
python3 analysis/kl_real/run_kl_real.py
```

GEO URLs (processed data only; no SRA):

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179502
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280232
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179500
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76628  (mismatch record)
