# FINDING — GSE76628 mouse Cldn4 vs T/NK

**ADDITIVE. Cldn4-only. No Tacstd2∩Cldn4 dual-high.** Inferential unit = mouse.
Public processed matrices only. This page does **not** audit the existing
thesis (Cldn4-high malignant/epithelial cells have lower IFN/MHC and tumors
have fewer T/NK; Cldn4-low/KD opens IFN). The series cannot test that claim.

**Verdict: NO-GO. Stop. Mouse-level table is empty (n_mice scored = 0).**

GSE76628 is a **bulk Affymetrix** series of **Ad-VEGF-A164 flank-skin**
angiogenesis lesions in **athymic nude mice**, deposited for gastric-cancer
stromal signatures (Uhlik et al., *Cancer Res* 2016, PMID
[27197264](https://pubmed.ncbi.nlm.nih.gov/27197264/)). It is not mouse lung
Kras/Lkb1. Sister deposits from the same study are the same model. Cldn4 %pos
in epithelial/malignant cells and T (`Cd3d`/`Cd3e`) or NK (`Nkg7`/`Klrd1`)
fractions were **not** computed.

---

## Decision

| Question | Answer |
|---|---|
| Public processed matrix? | **Yes** — GEO series matrix, 45,101 probes × **78** samples |
| Cldn4 on the platform? | **Yes** — `1418283_at` (not the failure) |
| T/NK genes on the platform? | **Yes** — `Cd3d`, `Cd3e`, `Nkg7`, `Klrd1` (not the failure) |
| Lung Kras/Lkb1 GEMM? | **No** — Ad-VEGF-A164 intradermal flank/ear |
| Epithelial / malignant lung compartment? | **No** — whole-flank bulk; no tumors |
| Cell-level Cldn4 %pos or T/NK fraction? | **No** — MAS5 bulk microarray |
| T cells in the animals? | **No** — Crl:NU(NCr)-Foxn1nu athymic nude |
| Same-study KP-KL processed matrix? | **No** — GSE76588 / SuperSeries GSE76630 are the same Ad-VEGF model |
| Mouse-level Cldn4 vs T/NK table? | **Empty** — `tables/mouse_units.tsv` has 0 rows |
| Epithelial Q4 vs Q1 IFN/MHC/TJ? | **Not run** (n_mice usable = 0, not ≥6) |
| Dual-high Tacstd2 ∩ Cldn4? | **Not defined** |
| Thesis treated as claim-failed? | **No** |

Stop rule from the assignment: *If the series is not usable (no Cldn4, T-only,
no matrix), say so and stop.* Here the matrix and Cldn4 exist, but the series
is still unusable for the locked primary (mouse-level epithelial/malignant
Cldn4 vs T/NK fraction in lung Kras/Lkb1). That is this page.

---

## What GSE76628 actually is

GEO [GSE76628](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76628)
(“Stromal-Based Signatures for the Classification of Gastric Cancer [part II]”).
Subseries of [GSE76630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76630).
Platform GPL1261, Affymetrix Mouse Genome 430 2.0. Processed data are MAS5
intensities (scaling factor 1500) in the series matrix.

| Field | Value |
|---|---|
| Organism | *Mus musculus* |
| n arrays (honest) | **78** |
| Tissue | Flank skin |
| Model | 5×10^8 PFU Ad-VEGF-A164 intradermal flank (and ear) |
| Strain | Female athymic nude, 4–6 weeks (Charles River Crl:NU(NCr)-Foxn1nu) |
| Design | Day-0 control + Ad-VEGF day 5 / 20 / 60 × NT / DC101 (anti-VEGFR2) / G6 (anti-VEGF) |
| n per group | 8, except day-60 DC101 **n=6** |
| Inferential unit as deposited | One bulk array per mouse |
| Cell types | None resolved |

Sample-level inventory: `tables/sample_inventory.tsv`. Group n:
`tables/sample_groups.tsv`. Figure: `figures/fig_sample_n_by_group.png`.

---

## Same-study deposits (not a hidden KP-KL matrix)

The assignment allowed any processed public matrix **on this series** or a
**related KP-KL GEMM if the same study deposits more**. PMID 27197264 deposits
only the Ad-VEGF mouse arrays:

| Accession | Relation | n | KP-KL lung? |
|---|---|---:|---|
| GSE76628 | this series (part II) | 78 | no |
| GSE76588 | sister (part I) | 20 | no — same flank Ad-VEGF model |
| GSE76630 | SuperSeries | 98 | no |

Human gastric TCGA/ACRG analyses in the paper are not a GEO processed matrix
on this series and are not mouse lung Kras/Lkb1. Public KP-KL scRNA series
from **other** studies (for example GSE165641, GSE180963) were **not** opened
and were **not** substituted.

---

## Why Cldn4 vs T/NK was not scored

1. **Wrong organ and genotype.** Locked context is mouse lung Kras/Lkb1.
   These arrays are VEGF-driven vascular/stromal lesions in normal skin.
2. **No epithelial/malignant scoring unit.** Cldn4 cannot be restricted to
   epithelial or malignant cells. Bulk skin Cldn4 is a tissue-fraction mix.
3. **No T/NK fraction.** `%pos` and T/NK fractions need cell-level labels.
   Bulk `Cd3d`/`Nkg7` intensity is not a fraction.
4. **No T cells.** Nude mice lack mature T cells. A mouse-level T fraction
   would be biology the animals do not have.
5. **Do not invent a claim test.** Running bulk Cldn4–Cd3d correlations on
   this wound-healing model would not be additive to the lung thesis and
   would look like a claim audit. That was not done.

Cldn4 and the requested T/NK genes **are** on GPL1261
(`tables/probe_presence.tsv`). Probe presence does not create a usable
compartment.

---

## Public files used

| File | Role | Used |
|---|---|---|
| `GSE76628_series_matrix.txt.gz` (10.8 MB) | processed MAS5 + sample metadata | **yes** — inventory only |
| `GSE76588_series_matrix.txt.gz` (3.5 MB) | sister metadata | **yes** — confirm same model |
| `GPL1261.annot.gz` (8.4 MB) | probe symbols | **yes** — Cldn4 / T/NK presence |
| `GSE76628_RAW.tar` (292 MB CEL) | raw | **no** |

No CEL reprocessing. No scRNA matrix exists on the series. Machine table:
`tables/public_files.tsv`. Gate: `tables/nogo_gate.tsv` and
`figures/fig_nogo_gate.png`.

---

## Mouse-level table (done condition)

`tables/mouse_units.tsv` exists and has **0 data rows**.

| Item | n |
|---|---:|
| Bulk arrays in GSE76628 | 78 |
| Mice with epithelial/malignant Cldn4 %pos | **0** |
| Mice with T or NK fraction | **0** |
| Mice entering Q4 vs Q1 IFN/MHC/TJ | **0** |

Primary contrast not run. Empty table is the honest n.

---

## Reproduce

```bash
python3 methods/gse76628_mouse_cldn4_tnk/scripts/inventory.py
```

Downloads the three public processed/annotation files into
`methods/gse76628_mouse_cldn4_tnk/data/` if absent and rewrites tables and
figures. Does not score expression.
