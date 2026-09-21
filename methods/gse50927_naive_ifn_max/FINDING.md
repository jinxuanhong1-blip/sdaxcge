# GSE50927 naive lung — max-effect IFN/ISG panel

**Cldn4 KO minus WT, no ventilator.** Kage / Borok *Am J Physiol Lung Cell Mol Physiol* 2014, PMID [25106430](https://pubmed.ncbi.nlm.nih.gov/25106430/); GEO [GSE50927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50927). The only processed genotype column for this contrast is the author edgeR table `GSE50927_Cldn4lungWTvsKOgenes.csv.gz`.

| Item | Value |
|---|---|
| KO | GSM1232581, Cldn4 KO, no VILI, whole lung |
| WT | GSM1232580, WT, no VILI, whole lung |
| **n** | **1 vs 1** |
| logFC sign | author edgeR, **KO minus WT** |
| Cldn4 logFC | **−6.061** (logCPM 2.25; deposited FDR 4.07×10⁻²⁶) |

The series page says the samples were prepared in duplicate. Those duplicates are not two expression columns. This analysis does not call the n 2, 5, or 10.

Positive logFC means the gene is higher after Cldn4 loss. That is the IFN-up direction.

## How the maximum was chosen

Eleven locked panels were scored. Membership comes from lists already used in this repository (IFN compact 26, the fable ISG / extended / core-8 lists, the kdko type-I and type-II lists, the Ayers 2017 IFN-γ 6-gene signature) and from MSigDB mouse Hallmark IFN-α, IFN-γ, their intersection, and their union (`mh.all.v2023.2.Mm.symbols.gmt`).

No gene was added or removed because its logFC was large or small. Discordant genes stay in the panel that lists them.

Two numbers are maximized together:

- **mean logFC** of the detected genes (KO minus WT)
- **gene concordance** = fraction of those genes with logFC > 0

A panel is on the Pareto front when no other locked panel is at least as high on both and higher on one. The joint scalar is mean logFC × fraction up. A separate call keeps the largest mean among panels whose gene-level sign test (fraction up > 1/2) has p < 0.05. That gate is a sign test, not a fraction cut placed between 5/6 and 23/26.

All of those p-values are **tests across genes on one unreplicated column**. They are not mouse-level p-values.

## Pareto front

| Panel | Genes up | Mean logFC | Joint | Gene sign p | What it wins |
|---|---:|---:|---:|---:|---|
| **Ayers IFN-γ 6** | **5/6** | **+1.565** | **1.304** | 0.11 | largest mean, largest joint |
| **IFN compact 26** | **23/26** | **+1.024** | **0.906** | **4.4×10⁻⁵** | largest mean among sign-concordant panels; largest joint among panels with ≥20 genes |
| **Fable ISG** | **23/25** | **+0.542** | 0.499 | **9.7×10⁻⁶** | highest fraction up |

Full grid: `tables/panel_scores.tsv`. Hallmark IFN-α / IFN-γ means are +0.281 / +0.260 (70/89 and 131/181). They are concordant only because the sets are large, and the means are diluted by genes that are not the compact ISG module.

### Ayers IFN-γ 6 — largest mean

Mouse symbols by the existing kdko ortholog rule. H2-Ea is absent, so HLA-DRA maps to H2-Aa.

| Gene | logFC | logCPM | Deposited FDR |
|---|---:|---:|---:|
| Cxcl9 | +2.883 | 1.34 | 0.062 |
| Cxcl10 | +2.543 | 0.76 | 0.094 |
| Ifng | +2.445 | **−0.61** | 0.59 |
| Ido1 | +1.328 | 0.48 | 0.46 |
| Stat1 | +0.391 | 6.95 | 0.29 |
| H2-Aa | **−0.198** | 8.48 | 0.89 |

Mean **+1.565**. Median +1.887. 5/6 up. Gene-level Wilcoxon p = 0.031. Sign test p = 0.11. The mean is the maximum. The concordance is not: one MHC-II gene is down, and 5/6 does not reject a fair coin.

### IFN compact 26 — the concordant maximum

Mean **+1.024**. Median +0.605. **23/26 up**. Pairwise same-sign fraction 0.788. Gene-level Wilcoxon p = 6.4×10⁻⁷. Mann-Whitney versus the other 21,562 genes p = 1.4×10⁻⁹. Same unit: genes, not mice.

The three genes that stay down, and were not dropped: **Ifit2 −0.369**, **Ifit3 −0.207**, **Eif2ak2 −0.140**. All three are detected (logCPM 5.38, 6.49, 4.98).

Two genes sit at logCPM ≤ 0 and pull the mean up: **Cxcl11 +5.044** (logCPM −1.07) and **Ifng +2.445** (logCPM −0.61).

### Fable ISG — highest fraction up

25 detected mouse symbols by the existing fable first-ortholog rule (OAS1 → Oas1a, not the higher-logFC paralog Oas1g). Mean **+0.542**. **23/25 up**. The only down genes are again Ifit2 and Ifit3. Every gene has logCPM > 1, so the abundance floor does not change this row.

## Abundance floor (same lists)

logCPM > 0 removes low-count genes. It does not build a new panel.

| Panel | After the floor | Mean logFC | Genes up |
|---|---|---:|---:|
| Ayers IFN-γ 6 | drops Ifng | **+1.389** | 4/5 (sign p = 0.19) |
| IFN compact 26 | drops Cxcl11 and Ifng | **+0.797** | **21/24** (sign p = 1.4×10⁻⁴) |
| Fable ISG | nothing dropped | +0.542 | 23/25 |

The compact mean **+0.797** is the same magnitude as the earlier high-minus-low compact score of −0.797, with the sign written as KO minus WT. The unfiltered compact mean +1.024 is larger because it keeps Cxcl11 and Ifng. Both numbers are in `tables/summary.json`.

## What this does not say

- It does not say n = 2. The duplicate libraries are not in the processed table.
- It does not say the gene-level sign test is a mouse-level p-value.
- It does not drop Ifit2, Ifit3, Eif2ak2, or H2-Aa to raise concordance.
- It does not use Oas1g (+1.192) in place of the locked Oas1a ortholog (+0.762).
- It is not a tumour, LUAD, or ICI contrast. VILI tables are a different experiment and are not in this panel.
- Deposited edgeR FDR values assume a dispersion on an unreplicated design. They are copied for Cldn4 and for the per-gene table. They were not used to choose the panel.

## Reproduce

```bash
python3 -m pip install -r methods/gse50927_naive_ifn_max/requirements.txt
python3 methods/gse50927_naive_ifn_max/analyze.py
```

The download goes to `$GSE50927_NAIVE_IFN_DATA` (default `/tmp/gse50927_naive_ifn`):

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz`

## Files

- `analyze.py` — score locked panels, write tables and figures
- `genesets/locked_panels.json` — compact, fable, kdko, and Ayers lists, plus the ortholog rules
- `genesets/mh.all.v2023.2.Mm.symbols.gmt` — MSigDB mouse Hallmark, IFN-α and IFN-γ
- `tables/panel_scores.tsv` — every panel × abundance filter
- `tables/panel_membership.tsv` — gene, logFC, logCPM, deposited FDR
- `tables/headline_genes.tsv` — Ayers, compact 26, fable ISG
- `tables/one_row.tsv` — winners and the Pareto front
- `tables/summary.json`
- `figures/fig1_max_effect_lollipop.png` — Ayers 6 and compact 26
- `figures/fig2_panel_pareto.png` — mean logFC for all 11 panels; red marks the Pareto front
