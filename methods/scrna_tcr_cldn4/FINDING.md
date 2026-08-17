# FINDING — TCR/BCR clone expansion vs malignant CLDN4

**GSE207422 has no public TCR or BCR.** The additive leftover that does is **GSE241934** (neoadjuvant IO+chemo; Hu/Zhang NEOTIDE + real-world). Same patients have Cell Ranger TCR clonotypes and author-epithelial CLDN4. On that join, clone expansion does **not** track malignant/epithelial CLDN4 (n=33, Spearman **ρ=−0.063, p=0.73**).

No public BCR was found (0 IGH/IGK/IGL contigs in the GSE241934 TCR tables).

---

## Hunt (same-patient TCR/BCR **and** CLDN4)

A series is usable only if a **public** VDJ/clonotype table and a **malignant or epithelial CLDN4** score exist for the **same patients**. Patient is the unit of n.

| Accession | Public TCR/BCR | Malignant/epi CLDN4 | Same-patient join? | Why |
|---|---|---|---|---|
| **GSE207422** | **No** | Yes (UMI, 15 pts) | **No** | BD Rhapsody WTA + bulk RNA only. GEO suppl = UMI matrix, sample xlsx, bulk TPM. No contig / clonotype / VDJ. Paper methods: Rhapsody WTA, no TCR/BCR library. Raw not deposited. |
| **GSE241934** | **Yes TCR** (43 pts; TRA/TRB). **No BCR** | Yes (author Epi on IIT+Real MTX) | **Yes** | Leftover neoadjuvant IO+chemo. `GSE241934_RAW.tar` = Cell Ranger `clonotypes.csv` + contig CSVs. Same `sampleID` as the author-Epi CLDN4 table. |
| GSE243013 | Yes TCR | No (CD45+ only) | No | Residual TACSTD2/CLDN4 in the immune MTX is not malignant epithelium. |
| GSE179994 | Yes TCR | No (T cells only) | No | Cannot score malignant CLDN4. |
| GSE176021 / GSE176022 | Yes (VDJ tars / bulk culture TCR) | No | No | Lymphocytes / MANAFEST bulk TCR. |
| GSE185204 / GSE186446 | Yes TCR (n=3) | No (sorted T) | No | CD3+ TIL. n=3. |
| GSE280232 | Yes (paired GEX+TCR GSMs) | No | No | GEO `cell type: Sorted T cells`. |
| GSE229353 | No public VDJ | No (CD45+) | No | Post-neoadjuvant immune GEX. |
| GSE291670 | No | Yes (snRNA MTX) | No | RAW.tar is MTX only. |
| GSE205335 | No | Yes (UMI RDS) | No | Identity + UMI only. |
| GSE267108 / GSE274595 / GSE337519 / E-MTAB-13526 | No | Epithelium on some | No | 2024–26 leftover matrices without public VDJ. |
| GSE308745 | Yes TCR+ADT | No (PBMC) | No | Blood only. |

Full hunt table: `results/hunt_tcr_cldn4.tsv`.

---

## GSE241934 — what was joined

TCR from GEO `GSE241934_RAW.tar` (43 patients). Clonotype = Cell Ranger `clonotype_id`. Expansion = clone size **≥2** (sensitivity ≥3). Metric = fraction of TCR cells in expanded clones. No BCR: contig chains are TRA 144,540 + TRB 164,732 + IGH/IGK/IGL **0**.

CLDN4 = mean log1p in **author epithelial** cells on the public IIT + Real matrices (same numbers as `methods/scrna_meta_mpr/results/per_patient_GSE241934.tsv`). This is **not** Hu GSE207422 CopyKAT. GSE241934 has no public CopyKAT IDs; author `major.cell.type` / Epi is the public malignant-ish label.

| Item | n |
|---|---|
| Patients with public TCR | **43** (all ≥315 TCR cells; median 3,981) |
| TCR patients with a CLDN4 value | **42** (P481 has 0 epithelial cells) |
| Primary floor: TCR ≥50 **and** ≥20 author-Epi cells | **33** (IIT 11 + Real 22; MPR 10 + NMPR 23) |
| Dropped for epi &lt;20 (TCR present) | 9 (P122, P33, P345, P348, P394, P483, P498, P579, P605) |
| CLDN4 table only, no TCR in RAW.tar | P223, P52 |

pCR counted as MPR. Cell-level p-values are not reported.

---

## Primary (pre-specified)

Spearman, two-sided, patient unit, n=33.

| Test | n | ρ | p |
|---|---|---|---|
| Author-Epi CLDN4 vs expanded-cell fraction (clone ≥2) | **33** | **−0.063** | **0.73** |
| CLDN4 vs expanded clones / 1k TCR cells | 33 | +0.19 | 0.29 |
| CLDN4 vs Shannon clonality | 33 | −0.15 | 0.39 |
| CLDN4 vs top-clone fraction | 33 | −0.018 | 0.92 |
| CLDN4 vs expanded-cell fraction (clone ≥3) | 33 | −0.10 | 0.58 |
| TACSTD2 vs expanded-cell fraction (clone ≥2) | 33 | −0.33 | 0.062 |

No CLDN4–expansion test is p&lt;0.05. Direction is near zero. TACSTD2 is a secondary gene and is still NS.

### Floors and cohorts (same CLDN4 vs clone≥2 expansion)

| Slice | n | ρ | p |
|---|---|---|---|
| epi ≥10 | 38 | +0.080 | 0.63 |
| epi ≥5 | 39 | +0.038 | 0.82 |
| IIT EGFR-mut only | 11 | +0.17 | 0.61 |
| Real-world EGFR-WT only | 22 | −0.20 | 0.36 |

Loosening the epithelial floor does not create a signal. The two cohorts point in opposite directions and both are NS.

### Secondary (not the CLDN4 test)

| Test | n | Result | p |
|---|---|---|---|
| Expanded-cell fraction, NMPR vs MPR | 23 vs 10 | median 0.499 vs 0.394; δ=+0.15 | 0.52 |
| Author-Epi CLDN4, NMPR vs MPR | 23 vs 10 | median 1.58 vs 1.51; δ=+0.19 | 0.40 |

Overall TCR expansion does not differ by MPR on this public join.

---

## What this is not

- Not GSE207422. That series has malignant CLDN4 and **no** public TCR/BCR.
- Not a CD8 Tex-CXCL13 or antigen-specific clone analysis. Public files are clonotype size only.
- Not BCR / TLS B-cell expansion. No BCR tables.
- Not Hu CopyKAT malignant IDs (those are still not on GEO for GSE207422, and GSE241934 uses author Epi).
- Not GSE243013 residual immune-compartment CLDN4. That is not malignant epithelium and is not used here.

---

## Reproduce

```bash
# GEO TCR archive (~115 MB); extract clonotypes + filtered contig CSVs
curl -L -o /tmp/GSE241934_RAW.tar \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE241nnn/GSE241934/suppl/GSE241934_RAW.tar
mkdir -p /tmp/gse241934/tcr
tar -xf /tmp/GSE241934_RAW.tar -C /tmp/gse241934/tcr \
  --wildcards '*_clonotypes.csv.gz' '*_filtered_contig_annotations.csv.gz'
TCR_DIR=/tmp/gse241934/tcr python3 methods/scrna_tcr_cldn4/scripts/analyze.py
```

## Files

- `results/patient_tcr_cldn4.tsv` — 45 patients (43 TCR + 2 CLDN4-only)
- `results/stats.tsv` — every test
- `results/hunt_tcr_cldn4.tsv` — series hunt
- `results/summary.json`
- `results/fig_cldn4_vs_tcr_expansion.png`
- `inputs/GSE241934_author_epi_cldn4.tsv` — author-Epi CLDN4 / TACSTD2
- `inputs/GSE241934_geo_filelist.txt` — GEO RAW.tar listing (clonotype CSVs)
- `scripts/analyze.py`
