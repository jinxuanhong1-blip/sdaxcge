# GSE117570 — malignant CLDN4 vs same-patient T/NK (CLDN4-only)

**Additive only.** TISCH NSCLC catalog already listed this object in `methods/tisch_nsclc_pool/` (sample-level, TACSTD2+CLDN4, Spearman skipped at n=2). This folder does **not** rewrite that pool. It re-scores **CLDN4 only**, with **patient as the unit**, on the public TISCH extract.

Song et al., *Cancer Medicine* 2019 ([PMID 31033233](https://pubmed.ncbi.nlm.nih.gov/31033233/); GEO [GSE117570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117570); TISCH `NSCLC_GSE117570`). Four treatment-naïve early NSCLC cases, tumor + adjacent normal, 10x 3′. No ICI / RECIST / MPR column on GEO or in TISCH CellMetainfo. **No dual-high.** TACSTD2 is an audit gene only (P1 library is alive).

Matrix: TISCH2 `expression.h5` (150.4 MB) + `CellMetainfo_table.tsv` (1.20 MB). GEO series matrix is metadata (0 expression rows). GEO `GSE117570_RAW.tar` UMI tables exist (~19 MB) and were **not** required; TISCH already carries lineage + patient.

## Verdict

CLDN4 is epithelial/malignant-restricted in the three patients where the gene is present. Same-patient T/NK cells are near floor. That is a **compartment** contrast, not an ICI or dual-high claim.

A patient-level Spearman of malignant CLDN4 vs T/NK **fraction** is **empty**. Honest n does not reach 5.

- Catalog patients: **4** (P1–P4).
- Patients with any CLDN4 > 0 in the TISCH h5: **3** (P2, P3, P4).
- **P1 is gene-empty for CLDN4** (0/4314 cells, tumor+NAT). EPCAM is also 0/4314. KRT19 / KRT7 / NAPSA / TACSTD2 are present, so this is not a dead library. P1 cannot be used as a CLDN4-low biological point.
- Tumor malignant ≥20: all 4. Tumor T/NK ≥20 **and** CLDN4 detectable: **1** (P2 only).
- Spearman (need n≥5): **not computed**.
- Paired malignant > T/NK CLDN4, detectable patients with ≥1 T/NK: **3/3** (P2, P3, P4). P3 has 3 T/NK cells; P4 has 8.

## Honest n

| item | n | source |
|---|---:|---|
| GEO / TISCH patients | **4** | P1–P4; 8 samples (tumor+NAT) |
| TISCH cells aligned | **11453** | CellMetainfo ∩ h5 barcodes |
| Tumor cells | 4872 | TISCH `Source=Tumor` |
| Tumor malignant | 1505 | major-lineage `Malignant` |
| Tumor T/NK | 640 | CD4Tconv + CD8T + NK (only T/NK labels present) |
| Patients with CLDN4 > 0 anywhere | **3** | P2 / P3 / P4 |
| P1 CLDN4-positive cells | **0** | empty row in this extract |
| Eligible for Spearman (mal≥20, T/NK≥20, CLDN4 detectable) | **1** | P2 only |
| Eligible paired T/NK≥5 + CLDN4 detectable | **2** | P2 (52 T/NK) and P4 (8 T/NK) |
| Spearman computed | **0** | n&lt;5 |
| ICI / RECIST / MPR labels | **0** | treatment-naïve; none on GEO or TISCH |

Do not write this as n=4 CLDN4-vs-infiltration. Do not recycle P1 as CLDN4-low / T/NK-high: P1 T/NK fraction is high (0.316) but CLDN4 was not measured.

## One-row table

| dataset | unit | n catalog | n CLDN4 detectable | n Spearman-eligible | CLDN4 vs frac T/NK | paired mal&gt;T/NK | dual-high | ICI |
|---|---|---:|---:|---:|---|---|---|---|
| GSE117570 TISCH NSCLC | patient | 4 | **3** | **1** | empty (n&lt;5) | 3/3 (P3 T/NK=3) | no | no |

## Per-patient tumor (unit = patient)

TISCH values are MAESTRO `log2(TPM/10+1)`. T/NK fraction = n_TNK / n_tumor_cells.

| Patient | Histology | TNM | Sex | n tumor | n mal | CLDN4 mal mean | CLDN4 mal %pos | n T/NK | frac T/NK | CD8/CD4/NK | CLDN4 T/NK mean | CLDN4 T/NK %pos | note |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---|
| P1 | LUAD | pT1b pN1 | Male | 1824 | 128 | 0.000 | 0.0 | 577 | 0.316 | 19/545/13 | 0.000 | 0.0 | CLDN4 empty in all patient cells (TISCH h5) |
| P2 | LUSC | pT2a pN0 | Female | 1297 | 199 | 1.440 | 87.4 | 52 | 0.040 | 13/38/1 | 0.041 | 3.8 | eligible for Spearman (n still <5 globally) |
| P3 | LUAD | pT2a pNx | Female | 328 | 208 | 2.226 | 98.1 | 3 | 0.009 | 0/3/0 | 0.000 | 0.0 | T/NK n<20 in tumor; fraction underpowered |
| P4 | LUAD | pT1a pN0 | Female | 1423 | 970 | 1.989 | 97.5 | 8 | 0.006 | 1/7/0 | 0.144 | 12.5 | T/NK n<20 in tumor; fraction underpowered |

P2 is the only LUSC. P3 and P4 tumors are T-poor after dissociation (3 and 8 T/NK cells). Paper text already notes CD8/NK are low in these tumors versus NAT.

## Paired compartment (malignant CLDN4 vs same-patient T/NK CLDN4)

Primary question as asked: malignant CLDN4 versus T/NK **in the same patient**, not a dual-high gate.

| floor | n | # mal &gt; T/NK | median mal | median T/NK | Wilcoxon p (greater) |
|---|---:|---:|---:|---:|---|
| CLDN4 detectable, T/NK ≥1 | 3 | 3 | 1.989 | 0.041 | 0.125 (n=3; descriptive) |
| CLDN4 detectable, T/NK ≥5 | 2 | 2 | 1.715 | 0.092 | n=2, not tested |

P2: 1.440 vs 0.041 (52 T/NK). P3: 2.226 vs 0 (3 T/NK). P4: 1.989 vs 0.144 (8 T/NK). P1 both 0 because the gene is missing.

Tumor-lineage restriction (all four patients pooled, so P1 zeros dilute malignant %pos): malignant CLDN4 mean 1.780 (88.0% pos, n=1505); T/NK lineages sit at the floor (CD4Tconv / CD8T / NK each &lt;1% pos in the full object). Extra figure `fig3` drops P1 so the restriction is visible.

## vs T/NK fraction (secondary; empty)

| contrast | n | ρ | p | note |
|---|---:|---|---|---|
| Malignant CLDN4 mean vs frac T/NK, Spearman-eligible | 1 | — | — | n=1 (P2) |
| Same, all CLDN4-detectable tumors | 3 | — | — | n=3&lt;5; P3/P4 T/NK-thin |

Descriptive ranks among P2/P3/P4: P2 has the lowest malignant CLDN4 (1.44) and the highest T/NK fraction (0.040); P3/P4 are CLDN4-high and almost T-empty. That is three points. It is not a Spearman.

## Extra figures

| file | what |
|---|---|
| `figures/fig1_paired_malig_vs_tnk.png` | same-patient malignant vs T/NK CLDN4 means |
| `figures/fig2_malig_cldn4_vs_frac_tnk.png` | patient scatter; P1 marked gene-empty |
| `figures/fig3_cldn4_by_lineage.png` | tumor violin, P2–P4 |
| `figures/fig4_tumor_composition.png` | TISCH lineage stack |
| `figures/fig5_pctpos_malig_vs_tnk.png` | % positive, both compartments |
| `figures/fig6_nat_vs_tumor_tnk.png` | paired NAT vs tumor T/NK fraction |
| `figures/fig7_cldn4_detectability.png` | cells with CLDN4&gt;0 (P1 = 0) |
| `figures/fig8_cldn4_mean_and_tnk_n.png` | malignant mean and raw T/NK counts |

## Methods (short)

1. Inputs: TISCH2 `NSCLC_GSE117570_expression.h5` and `CellMetainfo_table.tsv` only. GEO series matrix for diagnosis / TNM / sex / smoking. No FASTQ, no dual-high, no TACSTD2 test.
2. Values are TISCH MAESTRO `log2(TPM/10+1)`. Positive = value &gt; 0.
3. Malignant = TISCH major-lineage `Malignant`. T/NK = `CD4Tconv`, `CD8T`, `NK` (the only T/NK labels in this object).
4. Unit = **patient**. Primary tissue = `Source=Tumor`. NAT is extra (composition / T/NK drop), not a CLDN4 score.
5. Detectable = at least one cell in that patient with CLDN4 &gt; 0. P1 fails.
6. Spearman only if ≥5 patients with ≥20 malignant and ≥20 T/NK and CLDN4 detectable. That filter leaves **P2**.
7. Paired compartment = patient mean CLDN4 in malignant vs T/NK. Wilcoxon is reported only as a label; n=3 is not a powered test.
8. TISCH called 1,066 “Malignant” cells in P2 NAT. Those cells are not used for the tumor score. They are a TISCH over-call risk, not a second tumor.

## How to read this

- **Public processed object exists** (TISCH h5 &lt; 2 GB). This is not an empty accession.
- **Honest n for CLDN4 vs T/NK fraction is 1 eligible / 3 detectable / 4 catalog.** Spearman is empty.
- **Honest paired compartment is 3/3** among patients where CLDN4 exists, with P3/P4 T/NK counts of 3 and 8.
- **P1 cannot support a CLDN4-low story.** The gene is absent from the extract (and so is EPCAM).
- **No dual-high.** No ICI endpoint. No claim that CLDN4 tracks infiltration in this series.

This slice adds a patient-level CLDN4-only table and extra figures. It does not add a testable correlation.

## Files

- `analyze.py` — download TISCH + GEO matrix, score, write this page
- `tables/per_patient_tumor.tsv`, `per_sample.tsv`, `detectability.tsv`, `cldn4_by_lineage.tsv`, `stats.tsv`, `one_row.tsv`, `geo_sample_traits.tsv`, `summary.json`
- `figures/fig1_*.png` … `fig8_*.png` (and PDF)

```bash
python3 -m pip install -r methods/gse117570_cldn4/requirements.txt
python3 methods/gse117570_cldn4/analyze.py
```
