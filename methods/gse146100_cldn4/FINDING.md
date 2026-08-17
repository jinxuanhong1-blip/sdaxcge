# GSE146100 — CLDN4-only: epithelial (malignant proxy) vs T/NK

**Additive only.** Prior A3 TACSTD2 slice (PR 180) is not re-run and is not
re-argued. This folder is **CLDN4 only**. TACSTD2 is never a gate. Dual-high
is not run.

Zhang et al., *J Immunother Cancer* 2021 ([PMID 33820821](https://pubmed.ncbi.nlm.nih.gov/33820821/);
GEO [GSE146100](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146100);
TISCH2 `NSCLC_GSE146100`). One 72-year-old multi-primary LUAD after
pembrolizumab. Three surgically resected nodules: W1 (NR, EGFR L858R),
W2 (R, KRAS G12C), W3 (NR, EGFR L858R/R77H).

## Honest n

The task framed this as **n=3 patients**. That would overclaim.

| item | n | source |
|---|---:|---|
| Patients | **1** | GEO / TISCH `Patient` = Patient 1 / paper |
| Nodules (samples) | **3** | W1, W2, W3 |
| TISCH cells | 10,996 | h5 + CellMetainfo |
| GEO NormData cells | 11,612 | `GSE146100_NormData.txt.gz` (229 MB) |
| TISCH **Malignant** cells | **0** | Celltype (malignancy) is Immune / Stromal only |
| TISCH Epithelial (malignant proxy) | 975 | 112 / 211 / 652 in W1 / W2 / W3 |
| TISCH T/NK | 6,854 | CD4Tconv + CD8T + NK + Treg |
| Pretreatment samples | **0** | all three nodules are post-pembrolizumab |
| Independent patients for Spearman | **1** | not computed |
| Nodules for Spearman | **3** | n&lt;5; not computed |

n=3 nodules from 1 patient is **descriptive only**. No p-hacking: Spearman is
not computed, cell-level MWU is not used as evidence, and no cutoff / dual-high
/ R-vs-NR test is searched to manufacture a p-value. Response, driver genotype,
and clone identity are the same contrast (R = KRAS W2; both NR = EGFR).

## Verdict

- **Not a 3-patient cohort.** Public record is 1 patient / 3 nodules. Do not
  write this as n=3 patients.
- **No TISCH Malignant call.** “Malignant CLDN4” here is TISCH Epithelial
  (EPCAM/KRT-high lineage). Marker-rule sensitivity is reported separately.
- **CLDN4 is epithelial-restricted** (sanity, not a claim about T/NK):
  epithelial mean 1.63 / 86.7%+ vs T/NK mean 0.03 / 1.7%+.
- **Three descriptive points, not a correlation.** W3 has the highest
  epithelial CLDN4 and the lowest T/NK fraction. W1 vs W2 is not the same
  direction (W2 is slightly lower CLDN4 and also lower T/NK than W1). Do not
  upgrade W3 into a CLDN4–T/NK law.
- **No response claim.** n=1 patient; R vs NR is genotype-confounded; all
  samples are on-treatment.

## Primary table — TISCH epithelial CLDN4 vs T/NK (n=3 nodules)

TISCH2 `expression.h5` (191 MB, &lt;2 GB) + `CellMetainfo_table.tsv`. Values are
MAESTRO `log2(TPM/10+1)`. T/NK = CD4Tconv + CD8T + NK + Treg. Dual-high was
not run.

| Nodule | Response | Genotype | n epi | CLDN4 mean | CLDN4 %pos | n T/NK | T/NK % |
|---|---|---|---:|---:|---:|---:|---:|
| W1 | NR | EGFR L858R | 112 | 1.54 | 85.7 | 2890 | 75.2 |
| W2 | R | KRAS G12C | 211 | 1.44 | 89.1 | 2314 | 65.7 |
| W3 | NR | EGFR L858R/R77H | 652 | 1.71 | 86.0 | 1650 | 45.4 |

Pooled epithelial n=975. W3 holds most epithelial cells (652/975). %CLDN4+ is
similar across nodules (86–89%); the mean shift is W3-high, W2-low.

Spearman of epithelial CLDN4 vs T/NK fraction: **not computed** (n=3 nodules,
n=1 patient). The three points are in `tables/per_nodule_tisch.tsv` and
`figures/fig_scatter_cldn4_vs_tnk.png`.

T/NK mix (not tested):

| Nodule | CD4Tconv | CD8T | NK | Treg |
|---|---:|---:|---:|---:|
| W1 | 1550 | 1054 | 192 | 94 |
| W2 | 952 | 1222 | 107 | 33 |
| W3 | 1037 | 394 | 213 | 6 |

## Sensitivity — marker-rule epithelial (not a second claim)

Any of EPCAM/KRT8/KRT18/KRT19 &gt; 0 and PTPRC = 0. Same direction as the
TISCH epithelial table: W3 highest CLDN4 / lowest T/NK. Different
normalization, so means are not pooled with TISCH.

| Matrix | Nodule | n marker-epi | CLDN4 mean | CLDN4 %pos | T/NK fraction |
|---|---|---:|---:|---:|---:|
| TISCH h5 | W1 | 152 | 0.90 | 50.7 | 75.2 (TISCH labels) |
| TISCH h5 | W2 | 220 | 1.11 | 69.1 | 65.7 (TISCH labels) |
| TISCH h5 | W3 | 572 | 1.61 | 80.9 | 45.4 (TISCH labels) |
| GEO NormData | W1 | 159 | 0.80 | 49.1 | 69.9 (CD3D/CD8A/NKG7+, EPCAM=0) |
| GEO NormData | W2 | 240 | 0.92 | 63.7 | 66.4 (CD3D/CD8A/NKG7+, EPCAM=0) |
| GEO NormData | W3 | 616 | 1.45 | 76.3 | 46.6 (CD3D/CD8A/NKG7+, EPCAM=0) |

GEO has 11,612 cells (TISCH kept 10,996). Marker-rule T/NK on GEO is a
marker proxy, not the primary TISCH lineage fraction.

## Lineage sanity (CLDN4, TISCH, all nodules)

| Lineage | n | CLDN4 mean | CLDN4 %pos |
|---|---:|---:|---:|
| Epithelial | 975 | 1.63 | 86.7 |
| Treg | 133 | 0.06 | 3.0 |
| Mono/Macro | 1588 | 0.03 | 7.0 |
| Fibroblasts | 132 | 0.04 | 5.3 |
| DC | 594 | 0.05 | 4.5 |
| NK | 512 | 0.04 | 2.5 |
| Endothelial | 134 | 0.03 | 3.0 |
| Mast | 227 | 0.04 | 2.6 |
| B | 492 | 0.03 | 2.4 |
| CD4Tconv | 3539 | 0.03 | 1.8 |
| CD8T | 2670 | 0.03 | 1.6 |

CLDN4 is not a T/NK gene in this matrix. That is why the contrast is
epithelial CLDN4 vs T/NK **fraction**, not CLDN4 inside T/NK cells.

## What was not done

- No dual-high TACSTD2+CLDN4 gate.
- No TACSTD2 high/low split.
- No Spearman, Pearson, MWU, or Q4-vs-Q1 on n=3.
- No R-vs-NR test (1 vs 2 nodules, 1 patient, genotype-confounded).
- No SRA/FASTQ and no inferCNV (no matched normal; TISCH already has no
  Malignant call).
- Prior TACSTD2 A3 writeup is left as-is.

## How to read this

This is a public 1-patient, 3-nodule LUAD scRNA object. It can show that
CLDN4 is epithelial-restricted and that the three nodules differ in
epithelial CLDN4 mean and T/NK fraction. It cannot support a patient-level
CLDN4–T/NK association, an ICI-response claim, or a dual-high story.

## Files

| File | Role |
|---|---|
| `tables/one_row.tsv` | Honest-n one-liner |
| `tables/per_nodule_tisch.tsv` | Primary 3-nodule table |
| `tables/per_nodule_geo_marker.tsv` | GEO marker-rule sensitivity |
| `tables/cldn4_by_lineage.tsv` | Compartment sanity |
| `tables/lineage_counts.tsv` | TISCH mix per nodule |
| `tables/vs_tnk.tsv` | Spearman explicitly not computed |
| `tables/sample_metadata.tsv` | Nodule / GSM / genotype |
| `summary.json` / `sanity.json` / `audit.json` | Verdict |
| `figures/fig_bars_cldn4_tnk.png` | CLDN4 and T/NK bars |
| `figures/fig_scatter_cldn4_vs_tnk.png` | 3-point scatter (no fit) |

```bash
python3 methods/gse146100_cldn4/scripts/download.py
python3 methods/gse146100_cldn4/scripts/analyze.py
```
