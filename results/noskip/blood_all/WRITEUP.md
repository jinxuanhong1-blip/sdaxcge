# No-skip re-collection: TACSTD2 / CLDN4 in lung ICI blood, PBMC, plasma, CTC

Scope: every **usable public** human lung ICI series in blood / PBMC / plasma / CTC **where TACSTD2 or CLDN4 is actually in the public matrix**. Blood was **not** used as a drop reason. Absences are catalogued in `ABSENCES.tsv`, not silently omitted.

Search: NCBI GEO GDS eutils, 151 unique GSE (queries in `geo_search_candidates.tsv`). Probe table: `probe_catalog.tsv`.

## What was computed vs response

| Series | Compartment | n (groups) | TACSTD2 | CLDN4 | Honest read |
|--------|-------------|------------|---------|-------|-------------|
| **GSE249262** | CTC (stage III CRT → durvalumab) | 15 stable / 9 progression, baseline | p = 0.065, med 3.49 vs 3.91 | p = 0.15, med 5.12 vs 5.51 | Only biologically plausible blood compartment. Direction: slightly higher in progression. **Not significant.** Small n. |
| **GSE285888** | PBMC scRNA, 33 pts | 7 CR / 13 PD (also 8 DR) | 0.081% cells; p = 0.96 | 0.18% cells; p = 0.54 | Ambient / EPCAM-level. **No response association.** |
| **GSE111414** | sorted CD8 PBMC | 5 R / 5 NR, baseline | median log2CPM = 0 / 0; p = 0.80 | all zero | Genes present, **not expressed.** |
| **GSE202417** | sorted CD8 PBMC array | 6 R / 8 NR, pre-tx | p = 0.49 | p = 0.18 | Probes exist. Intensities similar to EPCAM on the same array — likely background in CD8. **No association.** |

No series gives a significant TACSTD2 or CLDN4 vs ICI response result. The closest is GSE249262 TACSTD2 (p = 0.065), which is a CTC array after CRT+durvalumab, not PBMC, and is underpowered.

## Genes measured but no public response label

These were **kept and quantified**, not dropped for being blood:

| Series | Compartment | Genes | What we could test |
|--------|-------------|-------|--------------------|
| GSE225620 | whole blood RNA-seq | TACSTD2 yes; CLDN4 ~0 | pre vs post neoadjuvant PD-1. TACSTD2 p = 0.72. Paper mentions responders; **GEO has no responder field.** |
| GSE235048 | PBMC RNA-seq TPM | TACSTD2 yes (TPM ~0.04–0.13); **CLDN4 absent from matrix** | pre vs post ICI. No RECIST. TACSTD2 p = 0.59. |
| GSE305086 | whole-blood GPL570 | both probes present | Detectability only. TACSTD2 mean-intensity percentile ~42 (near background). CLDN4 canonical probe `201428_at` ~26th percentile; `1569421_at` higher — **probe-discordant.** No public response/PFS. |

## Important absences (not dropped as blood — catalogued)

Full table: `ABSENCES.tsv`.

- **GSE216297** (n = 286 platelet RNA, nivolumab, Responder/nonResponder **is** in GEO): the public matrix is a **3805-gene low-abundance-filtered** count table. TACSTD2 and CLDN4 Ensembl IDs are **not in it**. This is the largest labeled blood RNA series, and it cannot be used for these genes from the public processed file.
- **GSE266219** (PBMC/CD8 scRNA, NR / aPD1.R / JAKi.R in sample titles): public files are Seurat RDS. Not extracted in this environment (no Seurat). **Not skipped for being blood.**
- **GSE152590**: CD8 TPM matrix lacks TACSTD2/CLDN4/EPCAM. No response labels.
- Plasma **miRNA** (GSE310370, GSE207715), **cfDNA 5hmC** (GSE237087), **CyTOF immunomap** (GSE295969): TACSTD2/CLDN4 mRNA not measured.
- Tumor series that matched the search (GSE136961, GSE93157, GSE161537, GSE311200): **not blood**.

## Method

- Patient (or baseline sample) is the statistical unit. Two-sided Mann–Whitney U.
- RNA-seq: log2 CPM or TPM as deposited. Arrays: deposited intensity; Clariom probes `TC0100014340.hg.1` (TACSTD2), `TC0700007993.hg.1` (CLDN4).
- scRNA (GSE285888): stream the 487 MB gene×cell UMI matrix; patient pseudobulk CP10K and detection rate.
- Positive-control genes (PTPRC, CD3D, CD8A, GZMB) are high in blood/PBMC as expected. Epithelial controls (EPCAM, KRT8/18) track TACSTD2/CLDN4 near zero except in CTC.

## Honest conclusion

Public lung ICI **blood/PBMC/plasma** data do **not** support TACSTD2 or CLDN4 as circulating response biomarkers.

1. In PBMC / CD8 / whole blood, both genes are at or below epithelial-ambient / array-background. Tests vs response are null (GSE285888, GSE111414, GSE202417).
2. The one compartment where the genes *should* be real — **CTC RNA (GSE249262)** — shows a non-significant trend toward higher baseline TACSTD2 in patients who later progressed (p = 0.065, n = 24). That is not a claim.
3. The largest labeled blood RNA cohort (GSE216297, n = 286) **does not publish these genes** after abundance filtering.
4. Several response-labeled blood series (GSE266219, GSE306542, GSE295969) are not usable as a TACSTD2/CLDN4 matrix from public files (Seurat-only, CyTOF, or incomplete GEO deposit).

Tumor tissue remains the appropriate public compartment for these two epithelial ADC targets.

## Reproduce

```bash
python3 scripts/noskip/blood_all/01_search_geo.py
python3 scripts/noskip/blood_all/02_probe_series.py
# downloads live under /tmp/geo_blood (not committed)
python3 scripts/noskip/blood_all/03_analyze.py
python3 scripts/noskip/blood_all/04_analyze_gse285888.py
```

Outputs in `results/noskip/blood_all/`: `TACSTD2_CLDN4_vs_response.tsv`, `response_stats.tsv`, `per_sample_expression.tsv`, `ABSENCES.tsv`, figures, GEO search/probe catalogs.
