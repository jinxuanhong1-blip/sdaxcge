# Public lung neoadjuvant/ICI scRNA meta (MPR / RECIST)

Additive methods package. Patient is the unit. Public processed matrices only.

## Question

In every public lung neoadjuvant or ICI scRNA series that has **MPR or RECIST** labels **and** malignant/epithelial cells, is malignant/epithelial **TACSTD2** or **CLDN4** higher in NMPR than MPR (or NR than R)?

Pool with a simple standardized-mean-difference meta (Hedges *g*) and an inverse-variance forest. Total **n_patients** is the number that matters.

## Inclusion

Entered in the malignant/epithelial MPR table only if all of:

1. Public processed counts (GEO supplementary MTX / UMI / RDS). No FASTQ, no DAC.
2. Human lung / NSCLC.
3. Neoadjuvant or ICI.
4. Public MPR (or pCR) labels, or public RECIST (held in a separate table).
5. Malignant or epithelial cells in the public object (author label preferred; marker epithelium if that is all GEO gives).

CD45-only series stay in a footnote.

## Must-try series

| Series | Labels | Compartment used | Role |
|---|---|---|---|
| GSE207422 (A3) | MPR 4 vs NMPR 8 post-tx; pCR with MPR | All epithelial (primary). Malignant-like = sensitivity. | A3 MPR contrast taken as given |
| GSE241934 IIT | Author Pathological Response | Author `major.cell.type == Epi` | Kept separate |
| GSE241934 Real | Author Pathological Response | Author Epi | Kept separate, then also pooled with IIT |
| GSE291670 | GEO titles MPR-1/2/3 vs Non-MPR-1/2/3 | Marker malignant | 3 vs 3 |
| GSE205335 | RECIST PR vs SD/PD | Author `lineage.sub == Malignant cells` | Separate RECIST table |

GSE207422 two MPR samples have **zero** malignant-like cells (normal-lung epithelium only). The A3-given primary row therefore uses all epithelial cells so that all four MPR patients remain. The malignant-like 8 vs 2 row is a sensitivity, not the pooled primary.

## Metric

Per patient, mean log1p expression in malignant/epithelial cells.

- GSE207422: mean log1p(UMI) in epithelial cells; ≥20 epithelial cells.
- GSE241934: mean log1p(UMI) in author Epi; ≥20 Epi cells.
- GSE291670: mean log1p(CP10k) in marker-malignant cells (all 6 patients).
- GSE205335: mean log1p(CP10k) in author malignant; ≥20 malignant cells.

Within-cohort standardization (Hedges *g*) is what is pooled. Raw means are not mixed across normalizations.

Contrast = NMPR − MPR (or NR − R). Positive *g* = higher in non-responders.

## Tests (per cohort)

- Two-sided Mann–Whitney U on patient means. Exact at 3 vs 3 (minimum two-sided *p* = 0.10).
- Hedges *g* with small-sample correction *J* = 1 − 3/(4*df* − 1); variance *(n1+n2)/(n1 n2) + g² / (2(n1+n2))*.
- Glass rank-biserial *r* = 2*U*/(*n1 n2*) − 1, with null SE from the MWU variance.

Cells are never the replicate.

## Meta

Primary: four independent MPR rows — GSE207422, GSE241934 IIT, GSE241934 Real, GSE291670.

Inverse-variance fixed effect and DerSimonian–Laird random effect. Report *Q*, *I²*, *τ²*, 95% CI, two-sided normal *p*.

Sensitivity: collapse IIT+Real into one patient-pooled row (same 53 patients, *k* = 3).

GSE205335 is **not** in the MPR pool.

## Leftover 2023–2026

GEO E-utilities (title/summary MPR or neoadjuvant scRNA, 2023–2026). Hits that fail the inclusion rule are listed in `results/leftover_inventory.tsv`. Notable exclusions:

- **GSE243013** — public MPR, *n* = 243, CD45-only. Footnote.
- **PRJNA1068179** (Molecular Cancer 2025) — paper has 3 NMPR / 2 MPR / 1 pCR; raw SRA only.
- GSE229353, GSE280232 — immune-sorted.
- GSE337519 — *n* = 1.
- GSE329813 — GeoMx, not scRNA.

## How to rerun

```bash
# leftover GEO titles (network)
python3 methods/scrna_meta_mpr/scripts/geo_leftover_search.py

# GSE241934 CLDN4 from public MTX (optional; table already stored)
bash methods/scrna_meta_mpr/scripts/download_gse241934.sh
python3 methods/scrna_meta_mpr/scripts/extract_gse241934_cldn4.py

# GSE205335 CLDN4 from public RDS (optional; table already stored)
# gunzip the GEO file twice, then:
Rscript methods/scrna_meta_mpr/scripts/extract_gse205335_cldn4.R

python3 methods/scrna_meta_mpr/scripts/run_meta.py
```

Large GEO objects stay under `/tmp` and are not committed.
