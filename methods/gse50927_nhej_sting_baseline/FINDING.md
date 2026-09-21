# GSE50927 naive lung — NHEJ, Sting1/Cgas, IFN/chemokine

**Additive public mouse. Non-cancer. Baseline only.** Kage / Borok *Am J Physiol Lung Cell Mol Physiol* 2014, PMID [25106430](https://pubmed.ncbi.nlm.nih.gov/25106430/); GEO [GSE50927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50927). Whole lung, mixed 129S6 / C57BL/6 / BALB/c, Illumina HiSeq 2000, mm9 EdgeR. This wave scores the **uninjured** genotype contrast only: Cldn4 KO with no ventilator versus WT with no ventilator. VILI tables are not re-scored.

The logic under test is the tumour chain **lower NHEJ activity → cGAS/STING → IFN and CCL5/CXCL9/CXCL10**. Here that chain is asked in a **non-cancer** Cldn4-null lung, as transcripts. Cldn4-high remains the WT lung. Cldn4-low remains the KO. The locked reading from the same contrast (IFN lower when Cldn4 is intact) is taken as given and is not re-audited.

The first pass below is the c-NHEJ core, the Sting1/Cgas pair, and the IFN/chemokine output. A second pass (`sweep.py`) adds other NHEJ panels, DNA-repair and p53 sets, a mitotic-error panel used as a micronucleus proxy, cytosolic-DNA sensors, IFN-gene removal, and leave-one-out. GSVA and ssGSEA are not run: the series matrix has no expression rows.

## Honest n

| Item | n | Used |
|---|---:|---|
| Baseline GSM, WT no VILI (GSM1232580) | 1 | yes |
| Baseline GSM, Cldn4 KO no VILI (GSM1232581) | 1 | yes |
| **Deposited pairwise n** | **1 vs 1** | yes |
| Author file | 1 | `GSE50927_Cldn4lungWTvsKOgenes.csv.gz` |
| SRA runs on these two GSM | 4 | no — not a count matrix |
| VILI GSM (WT, KOlow, KOhigh) | 3 | no — out of this wave |
| Series-matrix expression rows | 0 | no per-sample matrix |
| Tumour / LUAD / ICI | 0 | whole lung |

Design text says duplicates. The deposited EdgeR table is one row per gene, so the computable n stays **1 vs 1**. Author *P* / FDR assume a dispersion on that unreplicated design. They confirm Cldn4 loss and rank genes. They are not a mouse-level test. Gene-label permutation (10,000 draws, seed 42) asks whether a module mean sits away from a random gene set of the same size. It is also not a mouse-level *P*.

logFC sign is **KO minus WT**. The anchor is Cldn4 logFC **−6.061** (FDR 4.07×10⁻²⁶). A WT-minus-KO contrast would have made Cldn4 positive.

Module means use genes with logCPM ≥ 0. The zero-count floor in this table is logCPM = −2.071. Genes below logCPM 0 are listed and plotted, and they stay out of the mean.

## One-row table

| dataset | contrast | n | Cldn4 logFC | NHEJ core mean | NHEJ FDR hits | Cgas logFC (FDR) | Sting1 logFC (FDR) | IFN/chemokine mean | IFN/chemokine FDR hits | call |
|---|---|---|---:|---:|---:|---|---|---:|---:|---|
| GSE50927 naive lung | KO vs WT, **no VILI** | **1 vs 1** | **−6.061** | **−0.029** | **0/10** | **+0.498 (0.67)** | **+0.129 (1)** | **+1.122** | **5/13** | **IFN/chemokine up; NHEJ and Sting1/Cgas mRNA flat** |

Full numbers: `tables/one_row.tsv`, `tables/module_summary.tsv`, `tables/gene_level.tsv`.

## What each limb does

| limb | genes in mean | mean logFC (KO−WT) | median | up / down | FDR < 0.05 | gene-label perm *P* |
|---|---:|---:|---:|---|---:|---:|
| c-NHEJ core | 10/10 | **−0.029** | −0.019 | 4 / 6 | **0** | 0.81 |
| NHEJ extended (53BP1 / MRN / ligation accessories) | 14/14 | +0.013 | −0.001 | 7 / 7 | **0** | 0.91 |
| Sting1 + Cgas | 2/2 | +0.313 | +0.313 | 2 / 0 | **0** | 0.21 |
| STING neighbors (expressed) | 4/6 | +0.118 | +0.081 | 3 / 1 | 1 (Zbp1) | 0.50 |
| IFN / chemokine | 13/16 | **+1.122** | **+0.627** | **11 / 2** | **5** | **< 1×10⁻⁴** |
| TJ control, Cldn4 excluded | 6/6 | −0.032 | +0.023 | 3 / 3 | 0 | 0.83 |

Mann–Whitney of the IFN/chemokine logFCs against the other expressed genes is 2.85×10⁻⁵. Wilcoxon signed-rank of those 13 genes against 0 is 0.0017. The same tests on c-NHEJ are 0.75 and 0.77. The TJ control is the check that the KO column is not shifted end to end: other claudins, Tjp1, and Ocln stay near zero while Cldn4 itself is −6.06.

### c-NHEJ

All ten core genes are expressed. None have author FDR < 0.05.

| gene | alias | logFC | logCPM | FDR |
|---|---|---:|---:|---:|
| Xrcc6 | Ku70 | +0.311 | 4.05 | 0.69 |
| Xrcc5 | Ku80 | −0.011 | 4.02 | 1 |
| Prkdc | DNA-PKcs | **+0.049** | 4.73 | 1 |
| Dclre1c | Artemis | −0.356 | 3.60 | 0.76 |
| Nhej1 | XLF | +0.208 | 1.20 | 1 |
| Xrcc4 | XRCC4 | −0.173 | 3.27 | 1 |
| Lig4 | LIG4 | +0.266 | 3.75 | 0.82 |
| Poll | Pol λ | −0.092 | 3.38 | 1 |
| Polm | Pol μ | −0.467 | 4.10 | 0.21 |
| 1110057K04Rik | Paxx | −0.028 | 5.44 | 1 |

Prkdc, the gene the cancer experiments inhibit, is +0.05. The largest core move is Polm at −0.47 (FDR 0.21). The extended set is also balanced. The most negative extended genes, Rif1 (−0.450, FDR 0.12) and Rad50 (−0.365, FDR 0.20), stay above 0.05, and H2afx goes the other way (+0.501, FDR 0.22).

### Sting1 / Cgas

mm9 names in this 2014 table: **Cgas = Mb21d1** (Entrez 214763), **Sting1 = Tmem173** (Entrez 72512).

| gene | logFC | logCPM | *P* | FDR |
|---|---:|---:|---:|---:|
| Cgas (Mb21d1) | +0.498 | 2.14 | 0.130 | 0.67 |
| Sting1 (Tmem173) | +0.129 | 4.33 | 0.479 | 1 |
| Tbk1 | +0.117 | 5.59 | 0.465 | 1 |
| Irf3 | +0.045 | 5.51 | 0.763 | 1 |
| Zbp1 | **+0.997** | 2.57 | 0.0015 | **0.046** |
| Aim2 | −0.686 | 2.54 | 0.022 | 0.26 |
| Ifnb1 | 0 | −2.07 | 1 | not detected |
| Trex1 | 0 | −2.07 | 1 | not detected |

Both asked sensors point up. Neither clears FDR. Tbk1 and Irf3, the kinase and the transcription factor immediately downstream, are flat. Zbp1 is a separate DNA sensor and is the neighbor that does clear FDR. It is reported so it is visible. It is not counted as Sting1 or Cgas.

### IFN / chemokine

Pre-specified panel: CCL5, CXCL9, CXCL10, CXCL11, Ifnb1, Ifng, and a short ISG list that keeps both directions (Ifit2 and Ifit3 are down in this table). Mean uses the 13 genes with logCPM ≥ 0.

| gene | logFC | logCPM | FDR | in mean |
|---|---:|---:|---:|---|
| Ccl5 | **+2.710** | 4.96 | **0.0011** | yes |
| Cxcl9 | +2.883 | 1.34 | 0.062 | yes |
| Cxcl10 | +2.543 | 0.76 | 0.094 | yes |
| Isg15 | **+1.075** | 2.79 | **0.0031** | yes |
| Gbp4 | **+2.385** | 6.97 | **4.8×10⁻⁷** | yes |
| Rsad2 | +1.198 | 2.62 | 0.22 | yes |
| Mx1 | +0.627 | 1.37 | 0.65 | yes |
| Ifit1 | **+0.576** | 4.62 | **0.028** | yes |
| Oasl2 | **+0.542** | 5.24 | **0.018** | yes |
| Stat1 | +0.391 | 6.95 | 0.29 | yes |
| Irf7 | +0.227 | 3.98 | 0.98 | yes |
| Ifit2 | −0.369 | 5.38 | 0.19 | yes |
| Ifit3 | −0.207 | 6.49 | 0.97 | yes |
| Cxcl11 | +5.044 | −1.07 | 0.17 | no, low count |
| Ifng | +2.445 | −0.61 | 0.59 | no, low count |
| Ifnb1 | 0 | −2.07 | 1 | no, not detected |

Ccl5 is an FDR hit. Cxcl9 and Cxcl10 are large and sit just above FDR 0.05. The median of the 13 expressed genes is +0.63, so the mean is not one chemokine. Ifit2 and Ifit3 are the two expressed genes that go down.

## Logic reading

Predicted transcript pattern if Cldn4 loss copied enzymatic NHEJ inhibition: NHEJ genes down, Sting1/Cgas up, IFN and CCL5/CXCL9/CXCL10 up.

| limb | prediction | observed at baseline | reading |
|---|---|---|---|
| c-NHEJ, including Prkdc | down | mean −0.029, 0/10 FDR < 0.05 | flat |
| Cgas mRNA | up | +0.50, FDR 0.67 | direction only |
| Sting1 mRNA | up | +0.13, FDR 1 | flat |
| IFN / chemokine | up | mean +1.12, 11/13 up, perm *P* < 1×10⁻⁴ | **up** |

The output limb is present in this uninjured Cldn4-null lung, in the same direction as the locked IFN result on this contrast (IFN higher after Cldn4 loss; IFN lower in the Cldn4-intact lung). CCL5, CXCL9, and CXCL10 are the chemokines named on the tumour side of the same logic, and they are the high end of this rank.

The upstream transcript limbs do not move. Cancer NHEJ experiments inhibit DNA-PKcs activity. This table measures Prkdc mRNA, which is unchanged. cGAS–STING activation is a phosphorylation cascade. This table measures Sting1 and Cgas mRNA. Those two facts are why a flat NHEJ / Sting1 / Cgas mRNA result can sit next to a real IFN/chemokine rise. The rise is real in the deposited rank. The table does not identify NHEJ or Sting1/Cgas mRNA as its cause.

Zbp1 (FDR 0.046) is the one DNA-sensing transcript in the neighbor list that clears the author cutoff. That is a separate gene from Cgas and Sting1.

## Sweep

Flat NHEJ / Sting1 / Cgas mRNA in the first pass is not the last test. The same deposited rank was scored again with public sets and with explicit alternate panels. Enrichment is preranked GSEA (Subramanian 2005, weight 1, 1000 gene-set permutations, seed 42, a fresh stream for each set). Positive NES means the set sits at the KO-up end. IFN-γ NES on this engine is **+1.508**, the same number as the earlier GSE50927 pass. Means, a 5,000-draw gene-label permutation, and leave-one-out use genes with logCPM ≥ 0. None of these are mouse-level tests. n stays **1 vs 1**.

GSVA and ssGSEA were not given a number. Both need a genes-by-sample matrix. `GSE50927_series_matrix.txt.gz` has `Sample_data_row_count = 0` for all five GSM, and the only data row is the `ID_REF` header. Building either score from the logFC rank would just be preranked GSEA under another name.

Decontamination removes genes that sit in Hallmark IFN-α, Hallmark IFN-γ, or the phospho-proxy ISG list (205 symbols present in this table). A set that only looks high because it contains Cxcl10 or Isg15 is not an upstream bridge. A few 2023 MSigDB symbols are absent from this 2014 mm9 table (for example `Sgo1`, `Knl1`, `Ifi27`). They are listed in `tables/sweep_set_members.tsv` and were not imputed.

### What survives

| set | class | in mean | mean logFC | NES (nominal *P*, sweep FDR) | after IFN genes removed | leave-one-out |
|---|---|---:|---:|---|---|---|
| Chromosome missegregation (explicit mitotic-error panel) | proliferation | 18/18 up | **+0.728** | **+1.73 (0.003, 0.015)** | unchanged (0 IFN genes) | min mean **+0.668** without Aurkb |
| Hallmark E2F targets | proliferation | 188 | +0.231 | **+1.67 (0.001, 0.006)** | NES +1.65 | holds |
| Hallmark G2M checkpoint | proliferation | 180 | +0.203 | **+1.53 (0.001, 0.006)** | NES +1.54 | holds |
| Proliferation markers (Mki67, Top2a, Cdk1, Ccnb1, Mcms) | proliferation | 13 up / 1 down | +0.590 | **+1.56 (0.007, 0.031)** | unchanged (0 IFN genes) | holds |
| Hallmark mitotic spindle | proliferation | 191 | +0.125 | **+1.30 (0.015, 0.046)** | NES +1.31 | holds |
| Phospho-proxy ISG (same 13 genes as the first-pass IFN/chemokine mean) | output | 11 up / 2 down | **+1.122** | **+1.58 (0.010, 0.038)** | this is the output | holds |
| Hallmark IFN-γ | output | 175 | +0.218 | **+1.508 (0.001, 0.006)** | output | holds |
| Hallmark IFN-α | output | 88 | +0.227 | **+1.60 (0.001, 0.006)** | output | holds |

Author FDR hits inside the mitotic-error panel: Aurkb **+1.753** (FDR 0.0019), Top2a **+0.938** (FDR 7.8×10⁻⁴), Aurka **+1.368** (FDR 0.040), Bub1b **+1.022** (FDR 0.031). Mki67 is **+0.757** (FDR 0.027). Every gene that is present in the 18-gene mitotic-error panel has logFC > 0. `Sgo1` and `Knl1` were on the panel and are not in the table.

This is the strongest KO-up program that is not the IFN list. It is the transcriptional pattern expected if chromosomes are mis-segregating, which is the usual source of micronuclei. It is also the pattern expected if an IFN-high whole lung simply contains more cycling cells. Bulk lung at n = 1 vs 1 cannot separate those two readings. No micronucleus was counted.

### What does not bridge

Repair, p53 / UV / arrest, cytosolic-DNA sensors, and the nuclear envelope were the pre-specified bridge family. A call required a positive decontaminated NES, a leave-one-out mean that stays positive, and BH FDR < 0.05 inside that family. **None cleared.** The top of that family is KEGG p53 signaling: decontaminated NES **+1.33**, nominal *P* 0.032, family FDR **0.52**, mean logFC **+0.17**. The gene whose removal moves that mean most is Ccnb1, a cyclin. Single genes inside it are real (Pmaip1 FDR 5.7×10⁻⁷, Sfn FDR 0.009, Bbc3 +0.781 FDR 0.024) and they do not pull the set past the family cutoff.

| set | mean logFC | NES (nominal *P*) | after IFN genes removed |
|---|---:|---|---|
| c-NHEJ core (first pass) | **−0.029** | −0.50 (0.52) | unchanged, 0/10 FDR < 0.05 |
| KEGG non-homologous end-joining | **−0.021** | −0.74 (0.42) | 6 up / 6 down |
| Ku70 / Ku80 / Prkdc only | +0.116 | n = 3, no GSEA | Prkdc +0.049, FDR 1 |
| alt-NHEJ / MMEJ | +0.109 | +0.75 (0.38) | flat |
| Hallmark DNA repair | +0.056 | +0.87 (0.39) | mean +0.050 |
| KEGG homologous recombination | +0.051 | +0.85 (0.38) | flat |
| Nuclear envelope (Lmna, Lmnb1/2, Banf1, nucleoporins, …) | +0.059 | +0.73 (0.44) | 0 FDR < 0.05. Lmnb1 is **+0.37**, not down |
| Cytosolic DNA sensors (Cgas, Sting1, Tbk1, Irf3, Zbp1, Aim2, …) | +0.100 | +0.78 (0.38) | mean **+0.029**, NES +0.49. Leave-one-out **flips** when Mb21d1 is dropped (min mean −0.014) |
| KEGG cytosolic DNA-sensing pathway | +0.231 | +1.35 (0.027, sweep FDR 0.064) | mean **+0.095**, NES +1.11, *P* 0.15, family FDR 0.52 |
| Damage-arrest transcripts | −0.022 | +0.62 (0.47) | mean **−0.102** once Cdkn1a (an ISG, +0.921, FDR 2.0×10⁻⁶) is removed |
| Hallmark p53 | +0.018 | +0.73 (0.49) | flat |
| Hallmark UV response up | +0.002 | +0.72 (0.49) | mean −0.004 |

Zbp1 (**+0.997**, FDR 0.046) is the one DNA-sensor transcript with an author FDR below 0.05. It is also in the Hallmark IFN list, so the decontamination step does not let it stand in for an independent cGAS/STING bridge. The KEGG cytosolic-DNA leading edge before that step is Ccl5, Cxcl10, Ccl4, Zbp1. After the IFN genes are removed, the pathway is no longer an enrichment. RNASEH2 transcripts are mildly up and Trex1 is undetected, which is the wrong direction for a failure to clear cytosolic DNA.

The phospho-proxy available in RNA-seq is the ISG program itself (p-STING and p-IRF3 are not in this table). That proxy is up, with mean logFC +1.12. It is the output limb, not a new upstream link.

## What this does not test

- A tumour, LUAD, or ICI law.
- VILI, KOlow, or KOhigh. Those contrasts stay in the earlier GSE50927 IFN/MHC/TJ pass.
- Mouse-level replication. n = 1 vs 1 GSM.
- Phospho-STING, micronucleus counts, or DNA-PKcs kinase activity. The sweep scores transcripts only. GSVA and ssGSEA had no sample matrix to run on.
- A sample-level Spearman. There is no per-sample matrix.
- Dual-high Cldn4 + Tacstd2. Tacstd2 is not part of this wave.
- A revision of the locked naive-lung IFN result. The IFN/chemokine numbers here are the same direction as that result.

## Reproduce

```bash
python3 -m pip install -r methods/gse50927_nhej_sting_baseline/requirements.txt
python3 methods/gse50927_nhej_sting_baseline/analyze.py
python3 methods/gse50927_nhej_sting_baseline/sweep.py
```

The download (not committed) goes to `$GSE50927_NHEJ_STING_DATA` (default `/tmp/gse50927_nhej_sting`):

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz`

## Files

- `analyze.py` — download, panel means, gene-label permutation, figure
- `sweep.py` — alternate NHEJ panels, DNA-damage / mitotic-error / cytosolic-DNA sets, prerank GSEA, IFN-gene removal, leave-one-out
- `gsea_prerank.py` — weighted prerank engine, one permutation stream per set
- `genesets/mh.all.v2023.2.Mm.symbols.gmt` — MSigDB mouse Hallmark
- `genesets/kegg_mouse_selected.gmt` — Enrichr KEGG 2019 Mouse, the ten pathways used here
- `tables/sweep_summary.tsv`
- `tables/sweep_set_members.tsv`
- `tables/sweep_leave_one_out.tsv`
- `tables/sweep_bridge_rank.tsv`
- `tables/sweep_methods.tsv`
- `tables/sweep_call.json`
- `figures/fig2_sweep_nes.png`
- `tables/one_row.tsv`
- `tables/module_summary.tsv`
- `tables/gene_level.tsv`
- `tables/label_inventory.tsv`
- `tables/sample_annotation.tsv`
- `tables/summary.json`
- `figures/fig1_baseline_nhej_sting_ifn.png`
