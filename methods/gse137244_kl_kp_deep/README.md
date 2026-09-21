# GSE137244 KL vs KP — GSEA, ssGSEA, volcano, MW floor

Public mouse cell-line RNA-seq (Deng et al., Cancer Discovery 2021, GSE137244). Five KrasG12D;Lkb1 libraries versus five KrasG12D;Trp53 libraries. Additive to the locked Tacstd2 / Cldn4 means. Numbers and the IFN-cold call are in `FINDING.md`.

## Rerun

```bash
pip install -r methods/gse137244_kl_kp_deep/requirements.txt
python3 methods/gse137244_kl_kp_deep/analyze.py
```

The script downloads `GSE137244_counts.fpkm.csv.gz` into `cache/` (gitignored) and rewrites `tables/` and `figures/`.

## Outputs

| file | content |
|---|---|
| `tables/locked_genes.tsv` | Tacstd2, Cldn4, Stk11, Trp53 |
| `tables/tj_gene_callouts.tsv` | claudins and structural TJ genes |
| `tables/gsea_prerank.tsv` | hallmark (50) and junction (5) prerank |
| `tables/gsea_leading_edge.tsv` | leading-edge genes, Welch rank |
| `tables/ssgsea_contrasts.tsv` | per-library ssGSEA and exact MW |
| `tables/mean_signature_scores.tsv` | mean log2 scores, library and cluster |
| `tables/ifn_cold.tsv` | IFN-α / IFN-γ cold rule |
| `tables/mw_exact_ladder.tsv` | exact 5-vs-5 p for every U |
| `tables/gene_stats_expressed.tsv.gz` | gene-level Δ, Welch t, exact MW, BH |
| `figures/volcano_tj_callouts.png` | volcano with the MW ceiling |
| `figures/ssgsea_per_line.png` | ssGSEA by library |
| `figures/gsea_nes.png` | NES under two rankings |
| `figures/library_correlation.png` | why the KP arm is one block |
| `tables/ifn_hunt.tsv` | Hallmark versions, Reactome ISG, leave-one and KP-block flags |
| `tables/ifn_final_call.json` | soft IFN call |
| `tables/ifn_nhej_sting_primary.tsv` | IFN, NHEJ, and STING scores |
| `figures/ifn_leaveone.png` | Δ after dropping each library |
| `figures/nhej_sting_per_line.png` | NHEJ and STING by library |

`ifn_sensitivity.py` is the IFN / NHEJ / STING pass. It does not rewrite the Tacstd2, Cldn4, or TJ tables. Hallmark GMTs and the Reactome GMT are downloaded into `cache/msig/` on first run.

## Gene sets

- `gene_sets/mh.all.v2024.1.Mm.symbols.gmt` — MSigDB mouse hallmark 2024.1 (Broad).
- `gene_sets/junction_controls.gmt` — GO:0005923 bicellular tight junction (Enrichr GO Cellular Component 2025, case-matched to this matrix; five symbols absent, listed in `unmapped_go_bicellular.txt`), KEGG 2019 mouse tight junction (167/167 matched), and MSigDB mouse GO CC gap and adherens junctions. The 2024.1 mouse GO CC release has no bicellular-tight-junction term, which is why GO:0005923 is taken from Enrichr.
- `gene_sets/epithelial_and_isg.gmt` — 18-gene epithelial TJ core (Cldn1/3/4/7, occludin, MARVEL, ZO, JAM, cingulin, crumbs, angulins) and the 15-gene compact ISG list. Cldn6 is called out on the volcano and in GO/KEGG, and is not inside the core.

GSEA ranking is the Welch t of log2(FPKM+1), KL minus KP, with a second rank by the mean log2 difference. NES > 0 means enriched toward KL. ssGSEA uses rank normalization and weight 0.25, with no sample permutation; the 5-vs-5 test on those scores is exact Mann–Whitney.
