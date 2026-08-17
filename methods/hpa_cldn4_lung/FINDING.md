# HPA public CLDN4 (ENSG00000189143) — lung-cancer IHC + tissue nTPM

Official Human Protein Atlas **v25.1** only ([gene](https://www.proteinatlas.org/ENSG00000189143-CLDN4), [lung cancer](https://www.proteinatlas.org/ENSG00000189143-CLDN4/pathology/lung+cancer), [tissue RNA](https://www.proteinatlas.org/ENSG00000189143-CLDN4/tissue)). Retrieved 2026-08-17.

**Semi-quantitative ordinal IHC only.** HPA publishes High / Medium / Low / Not detected. Intensity and quantity are also ordinal. **No H-score is published and none was calculated here.**

Antibody for the lung-cancer TMA: **CAB002610** (Reliability IH = Enhanced). HPA does not print the strings LUAD / LUSC; those buckets are SNOMED histology on the same patient records (`Adenocarcinoma, NOS` → LUAD; `Squamous cell carcinoma, NOS` → LUSC).

Re-run: `python3 methods/hpa_cldn4_lung/download_and_tabulate.py`

## Count table — lung-cancer TMA IHC by LUAD vs LUSC

Source: official gene XML patient staining, recounted in `tables/04_ihc_counts_luad_lusc.tsv`. Matches the knowledge table total (`cancer_data.tsv`: High 0 / Medium 9 / Low 2 / Not detected 0, n=11).

| Histology (HPA SNOMED → bucket) | High | Medium | Low | Not detected | n | % Medium | % Low |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Adenocarcinoma, NOS → **LUAD** | 0 | 4 | 2 | 0 | 6 | 66.7 | 33.3 |
| Squamous cell carcinoma, NOS → **LUSC** | 0 | 3 | 0 | 0 | 3 | 100.0 | 0.0 |
| Other labeled (malignant NOS + carcinoid) | 0 | 2 | 0 | 0 | 2 | 100.0 | 0.0 |
| Unspecified | 0 | 0 | 0 | 0 | 0 | — | — |
| **All lung-cancer TMA** | **0** | **9** | **2** | **0** | **11** | **81.8** | **18.2** |

Knowledge table (`tables/02_knowledge_lung_cancer_counts.tsv`) is **not** histology-split. The LUAD / LUSC split exists only because every XML patient has a SNOMED histology label.

**Read as a small TMA, not a histology test.** n=6 vs 3. Both LUAD Low cases are adenocarcinoma. All three LUSC cores are Medium. That is a count, not a claim that LUSC is higher. No Fisher / chi-square is reported.

Other labeled (kept out of LUAD/LUSC):

- patient 496: Neoplasm, malignant, NOS — Medium
- patient 763: Carcinoid, malignant, NOS — Medium

The TMA is therefore **not NSCLC-only**.

## Patient-level ordinal scores (audit)

`tables/03_ihc_patients.tsv`. Staining is the HPA category used in the count table. Intensity and quantity are shown as published and were **not** multiplied into an H-score.

| Patient | Sex | Age | SNOMED histology | Bucket | Staining | Intensity | Quantity | Location |
| ---: | --- | ---: | --- | --- | --- | --- | --- | --- |
| 426 | M | 48 | Adenocarcinoma, NOS | LUAD | Medium | Moderate | >75% | cytoplasmic/membranous |
| 461 | M | 60 | Adenocarcinoma, NOS | LUAD | Low | Weak | 75%–25% | cytoplasmic/membranous |
| 537 | F | 50 | Adenocarcinoma, NOS | LUAD | Medium | Moderate | 75%–25% | cytoplasmic/membranous |
| 1249 | F | 44 | Adenocarcinoma, NOS | LUAD | Medium | Moderate | 75%–25% | cytoplasmic/membranous |
| 1303 | M | 68 | Adenocarcinoma, NOS | LUAD | Medium | Moderate | 75%–25% | cytoplasmic/membranous |
| 1327 | M | 64 | Adenocarcinoma, NOS | LUAD | Low | Moderate | <25% | cytoplasmic/membranous |
| 376 | M | 64 | Squamous cell carcinoma, NOS | LUSC | Medium | Moderate | >75% | cytoplasmic/membranous |
| 385 | M | 75 | Squamous cell carcinoma, NOS | LUSC | Medium | Moderate | 75%–25% | cytoplasmic/membranous |
| 1199 | M | 71 | Squamous cell carcinoma, NOS | LUSC | Medium | Moderate | 75%–25% | cytoplasmic/membranous |
| 496 | F | 75 | Neoplasm, malignant, NOS | other | Medium | Moderate | 75%–25% | cytoplasmic/membranous |
| 763 | M | 60 | Carcinoid, malignant, NOS | other | Medium | Moderate | >75% | cytoplasmic/membranous |

Public image URLs are in the patient TSV. No ICI, stage, or survival fields exist on these records.

## RNA nTPM — lung vs other tissues

Companion **RNA**, not IHC. Official `rna_tissue_consensus.tsv` (HPA consensus = max of HPA and GTEx nTPM per tissue). Full ranks: `tables/05_ntpm_tissue_consensus.tsv`. Summary: `tables/06_ntpm_lung_vs_other.tsv`.

HPA labels CLDN4 RNA as **tissue enhanced**, **detected in many**. RNA tissue-specific nTPM callouts on the gene TSV are esophagus 317.9, urinary bladder 282.7, intestine 265.3.

| Comparison | Value |
| --- | ---: |
| Lung consensus nTPM | **66.7** |
| Rank among 51 consensus tissues (high → low) | **16 / 51** |
| Median of the other 50 tissues | 14.0 |
| Mean of the other 50 tissues | 55.6 |
| Lung / other-tissue median | **4.76×** |
| Highest other tissue | esophagus **317.9** (4.8× lung) |

| Tissue | nTPM | Rank | vs lung |
| --- | ---: | ---: | --- |
| esophagus | 317.9 | 1 | 4.8× lung |
| urinary bladder | 282.7 | 2 | 4.2× |
| colon | 265.3 | 3 | 4.0× |
| salivary gland | 215.1 | 4 | 3.2× |
| skin | 187.8 | 5 | 2.8× |
| **lung** | **66.7** | **16** | **1.0×** |
| cervix | 66.5 | 17 | ≈ lung |
| stomach | 14.4 | 26 | 0.22× |
| liver | 4.5 | 31 | 0.067× |
| spleen | 0.7 | 39 | 0.010× |
| skeletal muscle | 0.5 | 42 | 0.007× |
| heart muscle | 0.4 | 46 | 0.006× |
| bone marrow | 0.2 | 49 | 0.003× |

Lung RNA is mid-high among barrier / secretory epithelia and far above muscle, heart, marrow, and most CNS tissues. It is **not** the HPA RNA peak (esophagus / bladder / colon). Do not treat nTPM as alveolar protein IHC: those are different assays.

## Not in HPA / not done

- No H-score, IRS, or percent-positive continuous score.
- Knowledge `cancer_data.tsv` has no LUAD vs LUSC columns.
- No paired TACSTD2+CLDN4 cores, no ICI labels, no survival on the TMA.
- n=11 total; histology split is descriptive only.

## Files

| File | Content |
| --- | --- |
| `tables/04_ihc_counts_luad_lusc.tsv` | Count table (this finding) |
| `tables/03_ihc_patients.tsv` | Patient ordinal IHC + SNOMED + image URLs |
| `tables/02_knowledge_lung_cancer_counts.tsv` | Official unsplit knowledge counts |
| `tables/05_ntpm_tissue_consensus.tsv` | All 51 consensus tissues |
| `tables/06_ntpm_lung_vs_other.tsv` | Lung vs other-tissue summary |
| `tables/01_download_manifest.tsv` | URLs, bytes, SHA-256 |
| `tables/summary.json` | Machine-readable copy of the counts |

## Citation

Human Protein Atlas version 25.1, https://www.proteinatlas.org (CC BY 4.0). Uhlén et al. Tissue-based map of the human proteome. *Science* (2015).
