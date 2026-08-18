# GSE50927 mouse whole lung — Cldn4-high vs low IFN / MHC / TJ

**Additive public mouse. Cldn4-only. No dual-high.** Kage / Borok *Am J Physiol Lung Cell Mol Physiol* 2014, PMID [25106430](https://pubmed.ncbi.nlm.nih.gov/25106430/); GEO [GSE50927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50927). Platform is Illumina HiSeq 2000 (GPL13112), mm9 EdgeR. Unit is the **deposited condition**. This is **ventilator-induced lung injury**, not a tumour series. Thesis is taken as given: **Cldn4-high epithelium is barrier / IFN-low**. This folder does **not** audit that as failed.

Cldn4-high = WT (Cldn4-intact). Cldn4-low = Cldn4 KO (Cldn4-null). **KOlow / KOhigh are BAL-protein injury strata, not Cldn4-expression strata.** No slide was re-scored. No ICI arm. Public processed EdgeR tables only.

## Honest n

| Item | Public? | n | Note |
|---|---|---:|---|
| GEO series linked from PMID 25106430 | yes | **1** | GSE50927 only |
| Related processed mouse series | **no** | **0** | E-GEOD-50927 is the same five GSM; Wray 2009 (PMID 19376879) has no GEO; hyperoxia arm is qRT-PCR only |
| GSM / BioSample / SRX | yes | **5** | one record per condition. Do **not** use 5 as a pairwise n |
| SRA runs | yes | **10** | 2 FASTQ runs per GSM. Not a processed count matrix |
| Series-matrix expression rows | **no** | **0** | `Sample_data_row_count=0` |
| Author EdgeR tables | yes | **4** | naive KO vs WT; WT VILI vs naive WT; KO VILIhigh vs WT VILI; KO VILIlow vs WT VILI |
| Cldn4 in those tables | yes | **4/4** | Entrez 12740 / `Cldn4` |
| Design-text “all in duplicates” | text only | claimed 2 | not recoverable as two processed columns |
| **Primary genotype n (naive KO vs WT)** | yes | **1 vs 1** | this is the n used below |
| Sample-level Cldn4 Q4 vs Q1 / Spearman | **no** | **0** | no per-sample matrix |
| Tumour / LUAD / ICI | no | 0 | mixed-background whole lung |
| Stage / ESTIMATE / pathologist TJ IHC | no | 0 | injury physiology series |

Do not write n=5, n=10, or n=2. The computable public genotype n is **1 vs 1 GSM** per contrast. Author EdgeR *P* / FDR assume a dispersion on an unreplicated design; they confirm the Cldn4 KO and the paper’s VILI genes. They are **not** a sample-level test. Prerank GSEA uses gene-set permutation on the deposited logFC rank (same engine as `methods/cldn4_ko_gsea`; 1000 perms, seed 42).

## One-row table

| dataset | contrast | n high vs low | Cldn4 logFC (deposited) | IFN-γ NES (FDR) | IFN-α NES (FDR) | MHC-I NES (FDR) | TJ NES (FDR) | IFN compact mean (high−low) | verdict | ICI |
|---|---|---|---:|---|---|---|---|---:|---|---|
| GSE50927 naive lung | KO vs WT (**primary**) | **1 vs 1** | **−6.061** | **+1.508 (0.005)** | **+1.587 (0.010)** | **+1.423 (0.050)** | −0.963 (0.28) | **−0.797** | **SUPPORTS** (not a fail) | no |
| GSE50927 WT VILI | VILI vs naive WT | 1 vs 1 | **+3.959** | −0.898 (0.28) | **−1.568 (0.010)** | **−1.580 (0.023)** | **+1.247 (0.050)** | — | **COMPANION** (induction, not genotype) | no |
| GSE50927 VILI KOhigh | KO VILIhigh vs WT VILI | **1 vs 1** | **−10.16** | **+1.679 (0.005)** | +1.219 (0.096) | +1.284 (0.096) | −0.772 (0.40) | **−0.322** | **SUPPORTS** | no |
| GSE50927 VILI KOlow | KO VILIlow vs WT VILI | **1 vs 1** | **−12.69** | +1.104 (0.25) | +0.612 (0.48) | +0.664 (0.48) | −0.922 (0.48) | +0.046 | **DIRECTION SUPPORTS** (FDR≥0.05) | no |

Full numbers: `tables/one_row.tsv`, `tables/gsea_headline.tsv`, `tables/compact_panel_scores.tsv`, `tables/label_inventory.tsv`.

Positive NES on a genotype contrast = the set is enriched at the **Cldn4-loss** end (KO minus WT). That is the IFN-low-in-Cldn4-high direction.

## Primary: naive Cldn4-high (WT) vs Cldn4-low (KO)

Author table `GSE50927_Cldn4lungWTvsKOgenes.csv.gz`. Mixed-cell whole lung, no ventilator. Cldn4 logFC = **−6.061** (logCPM 2.25; submitter FDR 4.1×10⁻²⁶) confirms the KO and the sign (KO minus WT). Tacstd2 is unchanged (logFC −0.145, FDR 1) and is **not** used.

| set | n in rank | NES | FDR | compact mean logFC (high−low) | compact Wilcoxon p | call |
|---|---:|---:|---:|---:|---:|---|
| Hallmark IFN-γ | 181/188 | **+1.508** | **0.005** | — | — | UP after Cldn4 loss |
| Hallmark IFN-α | 89/94 | **+1.587** | **0.010** | — | — | UP after Cldn4 loss |
| MHC-I / APM (21 genes; all present) | 21/21 | **+1.423** | **0.050** | — | — | UP after Cldn4 loss |
| KEGG tight junction | 154/169 | −0.963 | 0.28 | — | — | NS (Cldn4 is the hit) |
| IFN compact (24/26 used; Ifng / Cxcl11 low CPM out) | 24 | — | — | **−0.797** | 5.1×10⁻⁶ | IFN lower in WT |
| MHC compact (17/17) | 17 | — | — | **−0.432** | 0.0021 | MHC lower in WT |
| TJ compact **excluding Cldn4** (10/10) | 10 | — | — | +0.049 | 0.43 | other TJ unchanged |
| Epithelial compact (8/8) | 8 | — | — | +0.179 | 0.031 | mild WT-higher keratin / AT1 |

Lead IFN genes on the naive rank include Cxcl11, Cxcl9, Ccl5, Cxcl10, Isg15, Rsad2, Oasl1, Ifi44, Gbp3, B2m. Individual deposited logFC (KO minus WT): Cxcl9 +2.88, Cxcl10 +2.54, Isg15 +1.08 (FDR 0.0031), B2m +0.66 (FDR 0.034), Stat1 +0.39. Cldn3 +0.09 and Cldn18 −0.45 are not FDR-significant, matching the paper’s AT2 qPCR (other claudins unchanged).

This is the clean Cldn4-high vs Cldn4-null split. It **supports** the thesis. It cannot HOLDS: n=1 vs 1.

## Companion: WT VILI induces Cldn4 (not a high-vs-low cut)

`GSE50927_VILIwtGenes.csv.gz`. Cldn4 logFC = **+3.959** (16-fold; FDR 7.2×10⁻⁹¹). KEGG tight junction NES **+1.247** (FDR 0.050); GO keratinization NES +1.642 (FDR 0.010). Hallmark IFN-α NES **−1.568** (FDR 0.010); MHC-I NES **−1.580** (FDR 0.023). Compact MHC mean logFC (VILI minus naive) = −0.56 (Wilcoxon p=7.6×10⁻⁵).

When the same WT lung induces Cldn4, IFN-α / MHC go **down** and the TJ rank goes **up**. That is compatible with Cldn4-high = barrier / IFN-low. It is **not** entered as the primary genotype test. Tnf / Il1b / Il6 rise with VILI itself (injury program); they are not substituted for the IFN compact.

## Companion: VILI Cldn4-null vs Cldn4-intact

KOhigh / KOlow split mice by BAL protein, not by Cldn4 RNA. Both KO arms are Cldn4-null (Cldn4 logFC −10.2 and −12.7).

| contrast | IFN-γ NES (FDR) | IFN compact high−low | MHC compact high−low | paper injury genes (KO−WT) |
|---|---|---:|---:|---|
| KO VILIhigh vs WT VILI | **+1.679 (0.005)** | −0.322 (p=0.018) | −0.335 (p=0.00038) | Tnf +2.01 FDR 3.8×10⁻⁸; Il1b +1.73 FDR 2.6×10⁻⁷; Egr1 +2.13 FDR 3.3×10⁻¹⁷ |
| KO VILIlow vs WT VILI | +1.104 (0.25) | +0.046 (p=0.43) | +0.022 (p=0.71) | Tnf / Il1b / Egr1 NS vs WT |

The high-injury Cldn4-null lung is IFN-γ–high versus Cldn4-intact VILI. The WT-like-injury Cldn4-null lung is not. That is the paper’s own subtraction and is why KOlow is **not** a Cldn4-high vs low expression class.

## Related processed series

| accession | relation | usable? |
|---|---|---|
| GSE50927 | this series (4 EdgeR tables) | **yes** |
| E-GEOD-50927 | ArrayExpress mirror | duplicate |
| Wray 2009 PMID 19376879 | Cldn4 VILI physiology (CPE peptide / siRNA) | **no matrix** |
| Hyperoxia arm, same paper | qRT-PCR only | **no matrix** |

No second public processed mouse lung Cldn4 series was added. Human in-vitro CLDN4-loss ranks (GSE207704, GSE22493) live in `methods/cldn4_ko_gsea` and are not this mouse injury unit.

## What this does not test

- A tumour, LUAD, or ICI law.
- Sample-level Spearman or Cldn4 Q4 vs Q1 (no matrix).
- n=2 biological replicates as a processed column (design text only).
- KOlow / KOhigh as Cldn4-expression quartiles.
- A broad TJ program beyond Cldn4 itself at baseline (other claudins are unchanged).
- Dual-high Cldn4+Tacstd2.
- A failed audit of the thesis. Naive IFN / MHC are up after Cldn4 loss; WT VILI induces Cldn4 and drops IFN-α / MHC. Both are **additive** and **compatible**.

## Reproduce

```bash
python3 -m pip install -r methods/gse50927_mouse_cldn4/requirements.txt
python3 methods/gse50927_mouse_cldn4/analyze.py
```

Downloads (not committed) go to `$GSE50927_MOUSE_CLDN4_DATA` (default `/tmp/gse50927_mouse_cldn4`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/matrix/GSE50927_series_matrix.txt.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtGenes.csv.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkohiGenes.csv.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkoloGenes.csv.gz`

## Files

- `analyze.py` — download, rank, prerank GSEA, compact IFN/MHC/TJ scores, figures
- `gsea_core.py` — same weighted-KS prerank engine as `methods/cldn4_ko_gsea`
- `genesets/mh.all.v2023.2.Mm.symbols.gmt` — MSigDB mouse Hallmark (IFN-α / IFN-γ)
- `genesets/custom_human_sets.json` — MHC-I / KEGG TJ / keratin lists (mouse-mapped at run time)
- `tables/label_inventory.tsv` — honest n
- `tables/related_series.tsv`
- `tables/one_row.tsv`
- `tables/gsea_headline.tsv`
- `tables/compact_panel_scores.tsv`
- `tables/gene_level.tsv`
- `tables/gene_coverage.tsv`
- `tables/cldn4_logfc.tsv`
- `tables/sample_annotation.tsv`
- `tables/summary.json`
- `figures/fig1_gsea_nes.png`
- `figures/fig2_compact_high_minus_low.png`
