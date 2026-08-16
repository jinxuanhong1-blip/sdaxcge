# TACSTD2 / TROP2 loss-of-function transcriptomic datasets — verified inventory

Search date: 2026-08-16. Databases: NCBI GEO DataSets (`db=gds`, E-utilities) and
EBI BioStudies/ArrayExpress. Search scripts: `scripts/fable_tacstd2_kdko/01_search_geo.py`
plus manual E-utilities probes (see `results/fable_tacstd2_kdko/geo_search_*`).

## Search strategy (exhaustive)

1. Field-restricted queries on `[Title]`/`[Description]` combining `TACSTD2`/`TROP2`
   with perturbation terms (knockdown, knockout, shRNA, siRNA, sgRNA, CRISPR,
   silencing, depletion, deletion, loss, KO, KD, deficient, null, mutant).
2. Sample-tag probes: `shTACSTD2`, `siTACSTD2`, `sgTACSTD2`, `shTROP2`, `siTROP2`,
   `sgTROP2`, `TACSTD2 KO`, `TROP2 KO`, `Trop2 knockout`, `TACSTD2 shRNA`.
3. Cross-check against EBI ArrayExpress/BioStudies.

Bare `TACSTD2`/`TROP2` return ~33k / ~2k GEO records because they match any series
that merely *measures* the gene; those are not perturbation experiments and were
excluded. Broad `TACSTD2 shRNA` (129) / `TACSTD2 KO` (533) hits are CRISPR/RNAi
*library screens* or gene-list matches — none of them name TACSTD2/TROP2 as the
perturbed gene in title or summary (verified by scanning esummary docs).

ArrayExpress-only hits (E-MTAB-11466/11382/11377) are TROP2+ **cell-sorting**
(marker) studies, not perturbation, and E-GEOD-42112 == GSE42112 (marker). No
non-GEO perturbation datasets found.

## True TROP2/TACSTD2 loss-of-function datasets (all verified, processed data, <2GB)

| Accession | Organism | System | Perturbation | Design (n) | Processed file | Size | Axes |
|-----------|----------|--------|--------------|-----------|----------------|------|------|
| **GSE334497** | Mouse | 4T1 TNBC tumors | Trop2 **KO** (CRISPR) vs WT | 5 KO / 5 WT | `GSE334497_normalized_counts.csv.gz` | 1.2 MB | CLDN4, tight junction, IFN/immune |
| **GSE289287** | Human | T-47D breast, xenografts | Trop-2 **KO** vs WT (also DSG2 KO) | 4 KO / 3 WT (Trop2, xenograft) | author DESeq2 tables `GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz` (+DSG2) | ~2 MB each | desmosome/junction |
| **GSE245459** | Human | Ovarian cancer cells | sh**TACSTD2** vs shNC (± cisplatin) | 3 sh / 3 shNC (no drug) | `GSE245459_fpkm.anno.txt.gz` | 14 MB | KD context |
| **GSE15212** | Human | Colorectal (RKO) cells | **TACSTD2** siRNA vs neg-ctrl siRNA | 6 siTACSTD2 / 12 neg-ctrl (72h) | per-sample Agilent tables (GPL4133) | small | KD, independent |

### Excluded (marker / not perturbation)
- GSE42112 / GSM10326xx — TROP2 used to *sort* trophoblast cells (marker).
- GSE75748, GSE114326, GSE227698, GSE280834/835/836/876, GSE33348 — TACSTD2/TROP2
  measured or mentioned, but the perturbed gene is something else (FOXA1, CXCR2,
  G0S2, Vav2/3, etc.).

## Notes on statistics
- GSE289287 ships **author DESeq2 results** (log2FC, p, padj) → used directly for the
  TROP2-loss contrast; the DSG2 KO tables are used as a desmosome-partner comparison.
- GSE334497 provides normalized counts only (no raw) → Welch t-test + Benjamini-Hochberg
  FDR on log2 normalized counts, with Mann-Whitney U as a rank-based robustness check.
- GSE245459 provides FPKM → t-test + BH FDR on log2(FPKM+1) (n=3/3, treated as
  exploratory; caveated).
- GSE15212 is a two-colour Agilent array → per-sample processed VALUE assembled into a
  matrix; Welch t-test + BH FDR on the TACSTD2-siRNA vs neg-ctrl (72h) contrast.
