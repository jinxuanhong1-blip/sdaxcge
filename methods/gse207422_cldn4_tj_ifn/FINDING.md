# GSE207422 — malignant CLDN4-high vs low: TJ/keratin and IFN in the same cells

**Additive slice.** This is not dual-high. CLDN4 is the splitter. TACSTD2 is a companion gene and is **never a gate**. CLDN4 is **excluded** from the TJ module (circular otherwise). The given A3 TACSTD2 / dual-high CLDN4-vs-immune write-ups are not re-argued.

**Verdict (honest n=6 paired patients, not 12):** in A3-malignant cells, CLDN4-high (within-tumor Q4) is higher than CLDN4-low (Q1) for **both** the TJ/keratin module and the IFN ISG module **in the same cells**. All 6 eligible post patients go the same way on both modules (exact Wilcoxon p=0.031, the floor at n=6). The two deltas move **together** (ρ=+0.83, p=0.042), not as a TJ-up / IFN-down tradeoff. Keratin alone is weaker (4/6, p=0.16). An OXPHOS control does not move (4/6, p=0.22).

Do **not** cite n=12. Six of twelve post patients have <20 A3-malignant cells (including 3 of 4 MPR). Cell-level p-values are exploratory (pseudoreplication).

## Data and n

- Public GEO UMI only: **92,330** cells × **24,292** genes. Author CopyKAT / epithelium RDS barcodes are not on GEO.
- **Unit of every primary test is the patient.** Within each patient, A3-malignant cells are split on CLDN4 `log1p(CP10k)` Q4 vs Q1. Module score = mean `log1p(CP10k)` of the gene set. Paired Wilcoxon is two-sided exact.
- Marker lineages (same argmax rule as the dual-high CLDN4 slice): epithelial **11,019**; A3-malignant-like **6,627**.
- A3-malignant-like = epithelial AND zero UMI for SFTPA2/AGER/SCGB1A1/SCGB3A1/TPPP3 (**not** CopyKAT).
- Primary floor: ≥20 A3-malignant cells and ≥8 cells in each CLDN4 tail. **Eligible post n=6:** P03 (MPR), P04, P07, P09, P10, P12 (NMPR).
- Dropped post (A3-malignant n): P02 (10), P06 (1), P11 (0), P13 (3), P14 (0), P15 (4). Three of four MPR are empty or one cell.
- Missing from the public matrix: KRT14, KRT6A, KRT6B, SFTPC (same SFTPC hole as the CytoTRACE-like slice). Keratin uses KRT7/8/18/19/5/17.

## Definitions

| Item | Rule |
|---|---|
| CLDN4 split | within-patient Q4 vs Q1 of malignant `log1p(CP10k)` |
| TJ/keratin | CLDN1/3/7, OCLN, TJP1/2, F11R, PARD3, MARVELD2, CGN, CRB3, JAM3 + KRT7/8/18/19/5/17. **No CLDN4.** |
| TJ only | the 12 TJ genes above (no CLDN4, no keratins) |
| IFN ISG | 40-gene type-I ISG core (ISG15, MX1, OAS*, IFIT*, STAT1/2, IRF7/9, …) |
| IFN core-6 | IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A |
| MHC-I APM | HLA-A/B/C/E/F, B2M, TAP1/2, TAPBP, PSMB8/9/10, … |
| OXPHOS control | NDUFA1, COX5A, ATP5F1A/B, UQCRC1, SDHA, … |
| Primary n | post patients with ≥20 A3-malignant and ≥8 cells per tail |

## Primary table — same cells, paired patients

A3-malignant, 12 post patients attempted, **6 eligible**. High = CLDN4 Q4; low = CLDN4 Q1. Mean is the mean of patient-level module means.

| Test | n | high vs low | Δ median | high>low | W | p |
|---|---:|---|---:|---:|---:|---:|
| TJ/keratin (CLDN4 excluded) | **6** | 0.793 vs 0.597 | +0.133 | 6/6 | 0 | **0.031** |
| TJ only (CLDN4 excluded) | 6 | 0.495 vs 0.280 | +0.179 | 6/6 | 0 | **0.031** |
| keratin (simple+basal) | 6 | 1.388 vs 1.232 | +0.132 | 4/6 | 3 | 0.16 |
| IFN ISG core | **6** | 0.442 vs 0.368 | +0.054 | 6/6 | 0 | **0.031** |
| IFN user core-6 | 6 | 1.267 vs 1.141 | +0.133 | 6/6 | 0 | **0.031** |
| MHC-I APM | 6 | 1.413 vs 1.251 | +0.211 | 5/6 | 1 | 0.063 |
| OXPHOS control | 6 | 0.947 vs 0.869 | +0.081 | 4/6 | 4 | 0.22 |
| companion TACSTD2 | 6 | 2.074 vs 0.999 | +1.044 | 6/6 | 0 | 0.031 |
| library size (total UMI) | 6 | 14582 vs 6335 | +8144 | 5/6 | 1 | 0.063 |

P07 IFN ISG delta is +0.0008 (visually flat). The other five IFN deltas are +0.048 to +0.179. Exact two-sided p cannot go below 0.031 at n=6 when all signs agree.

### Same-cell deltas and sample-level Spearman (still n=6)

| Test | n | stat | p |
|---|---:|---|---:|
| paired-delta Spearman, TJ/keratin vs IFN ISG | 6 | ρ=**+0.829** | **0.042** |
| patient-mean CLDN4 vs TJ/keratin | 6 | ρ=+0.600 | 0.21 |
| patient-mean CLDN4 vs TJ only | 6 | ρ=+0.943 | 0.0048 |
| patient-mean CLDN4 vs IFN ISG | 6 | ρ=+0.943 | 0.0048 |
| patient-mean CLDN4 vs keratin | 6 | ρ=−0.086 | 0.87 |
| patient-mean CLDN4 vs TACSTD2 | 6 | ρ=+0.771 | 0.072 |

The TJ/keratin combined module is carried by the TJ genes, not by keratin. Sample-level CLDN4 vs IFN is the same ρ as CLDN4 vs TJ-only. Deltas of the two modules are positively correlated: patients with a larger TJ lift also have a larger IFN lift.

## Epithelial complete-case (sensitivity; not the primary)

All 12 post patients have ≥20 epithelial cells. Paired tails require ≥8 cells/side: **n=11** (P13 n=28 drops, tail<8). This includes residual / normal-lung epithelium that A3-malignant removes.

| Test | n | high vs low | high>low | p |
|---|---:|---|---:|---:|
| TJ/keratin | 11 | 0.701 vs 0.545 | 11/11 | 0.00098 |
| IFN ISG | 11 | 0.405 vs 0.345 | 10/11 | 0.014 |
| OXPHOS control | 11 | 0.807 vs 0.765 | 7/11 | 0.28 |
| total UMI | 11 | 14137 vs 9267 | 9/11 | 0.0068 |

Sample-level epithelial means, n=12 complete: CLDN4 vs TJ/keratin ρ=−0.16 p=0.62; CLDN4 vs IFN ISG ρ=+0.54 p=0.071. NMPR vs MPR epithelial CLDN4 is null (Δ=+0.005, p=0.93, n=8 vs 4), matching the dual-high CLDN4 slice. A3-malignant NMPR vs MPR is **5 vs 1** after the occupancy floor and is not tested.

## Cell-level (exploratory; do not cite as n)

Post A3-malignant cells n=4,683 (pseudoreplication). CLDN4 vs TJ/keratin ρ=+0.334; vs TJ-only ρ=+0.452; vs IFN ISG ρ=+0.231; vs keratin ρ=+0.014 (NS); vs OXPHOS ρ=−0.200. Sign matches the paired patient tests. These p-values treat cells as independent and are not a claim.

## Honest limits

1. **n=6, not 12.** Half the post cohort has almost no A3-malignant cells. MPR residual tumor is the main hole (P06=1, P11=0, P14=0). Header n stays 12 attempted / 6 tested.
2. Exact Wilcoxon p=0.031 is the smallest two-sided p at n=6. One IFN pair (P07) is a near-tie.
3. CLDN4-high cells have more UMI (p=0.063). Scores are CP10k-normalized. OXPHOS does not significantly rise, but depth is not fully ruled out.
4. This is not Hu et al. CopyKAT. A3-malignant-like can leak unmarked epithelium.
5. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 tracks the CLDN4 split (6/6) and is reported only as a companion.
6. KRT6A/B/14 and SFTPC are absent from the public UMI.

## Extra figures

- `figures/fig_extra_paired_tj_ifn.png` — same-cell paired TJ/keratin and IFN (n=6)
- `figures/fig_extra_delta_scatter.png` — ΔTJ/keratin vs ΔIFN in those patients
- `figures/fig_extra_honest_n.png` — A3-malignant occupancy; floor n=20 visible
- `figures/fig_sample_cldn4_vs_modules.png` — patient-mean CLDN4 vs each module

## Files

- `primary_table.tsv` — compact primary table
- `per_patient.tsv` — 15 samples; tests use post rows
- `paired_high_low.tsv` / `eligibility.tsv` / `tests.tsv`
- `malignant_epithelial_scores.tsv.gz` — epithelial + A3-malignant cell scores
- `summary.json` / `sample_metadata.tsv` / `lineage_counts.tsv`
- Scripts: `scripts/download.py`, `scripts/extract.py`, `scripts/analyze.py`, `scripts/gene_sets.py`

## Reproduce

```bash
python3 methods/gse207422_cldn4_tj_ifn/scripts/download.py
python3 methods/gse207422_cldn4_tj_ifn/scripts/extract.py
python3 methods/gse207422_cldn4_tj_ifn/scripts/analyze.py
```
