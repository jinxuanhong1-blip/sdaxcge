# GSE180963: Tacstd2 / Cldn4 in KL vs K lung GEMM scRNA (w200)

**Question.** In the public KrasG12D/+ (K) vs KrasG12D/+;Lkb1fl/fl (KL) lung
GEMM scRNA series `GSE180963`, what is the expression of `Tacstd2` (TROP2) and
`Cldn4` in tumor epithelium vs immune cells?

**Bottom line (honest).** Both genes are **epithelial-restricted** in this
series. In the KL sample, `Tacstd2` is on in **63%** of cluster-called tumor
epithelial cells vs **3.6%** of immune cells. `Cldn4` is also epithelial-restricted
but **sparse**: **10%** of KL epithelium vs **0.03%** of immune cells (2 cells).
The two genes barely co-occur (9% of KL epithelium double-positive; Spearman
ρ = 0.16). **n = 2 (1 mouse per genotype).** We do **not** claim a KL-vs-K
difference, and we do **not** independently confirm that this KL sample is
immune-cold: immune fraction is ~82% in both samples.

---

## Dataset (public, processed)

GSE180963 — *Single cell RNA sequencing of tumor sections from GEMM harboring
KrasG12D/+ or KrasG12D/+Lkb1fl/fl (KL) mutation* (Bai et al., Southern Medical
University; public 2022-07-16). Mouse FVB, 10x Chromium 3′ v3.1, mm10. The
series summary explicitly frames LKB1 (`Stk11`) loss as producing an
immune-desert / “cold” TME.

| GSM | label | genotype | cells (author-filtered) |
|---|---|---|---|
| GSM5481386 | **K** | KrasG12D/+ | 6,696 |
| GSM5481387 | **KL** | KrasG12D/+ ; Lkb1fl/fl | 7,564 |

GEO supplementary `GSE180963_RAW.tar` (88 MB) holds author-filtered count
matrices (`K/` and `KL/` as `matrix.mtx` + `genes.tsv` + `barcodes.tsv`; 20,304
genes, no duplicate symbols). Matrices already QC’d by the authors (Seurat
3.1.5: 500–6000 features, mito < 20%, genes in ≥3 cells). The two samples were
**pooled into one 10x library** and demultiplexed by label. No unfiltered
matrices; SoupX/CellBender is not possible from what GEO hosts. No SRA FASTQ
reprocessing. No private cohort.

`Tacstd2` and `Cldn4` are both present in the gene table.

## What was done

`scripts/w200_gse180963/02_analyze.py`: load both MTX → re-apply author QC
(0 cells dropped) → CP10k + log1p → HVG(2000) / PCA(50) / neighbors / UMAP /
Leiden(res=1.0) → per-cell lineage scores → cluster majority-vote annotation →
score `Tacstd2` and `Cldn4` by cell type / compartment / genotype →
clustering-free epithelial rule as a robustness check → Tacstd2–Cldn4
co-expression in epithelium. Seed 0.

**No batch integration.** Genotype is the sample. “Correcting” batch would
erase the only contrast in the series.

`Cldn4` is **not** used to call epithelium (avoids circularity). Epithelial
markers: `Epcam`, `Krt8/18/19`, `Cldn18`, `Cdh1`, `Nkx2-1`. Marker-rule
epithelial = Epcam+ AND (keratin or Cldn18)+ AND Ptprc−. `Sftpc` is treated as
an **ambient-RNA flag**, not a caller.

## Composition (the “cold TME” is not visible as fewer immune cells)

| genotype | epi (cluster) | immune (cluster) | stromal | % Ptprc+ | % Sftpc+ (ambient) |
|---|---|---|---|---|---|
| K | 20 (0.3%) | 5,477 (81.8%) | 1,199 (17.9%) | 82.8 | 76.5 |
| KL | 231 (3.1%) | 6,180 (81.7%) | 1,153 (15.2%) | 82.3 | 77.7 |

Both samples are immune-dominated. Immune fraction is **identical**. The
LKB1-loss “cold” phenotype, if present here, is not a drop in CD45+ fraction
and was not tested as T-cell exclusion or myeloid skew. K epithelium is
**tiny** (20 cluster / 28 marker-rule cells). Only one Leiden cluster
(cluster 13) scores as epithelium.

## Tacstd2 (TROP2)

| definition | set | n | % Tacstd2+ | mean log1p CP10k |
|---|---|---|---|---|
| cluster | KL epithelium | 231 | **62.8** | 0.577 |
| cluster | KL immune | 6,180 | 3.6 | 0.039 |
| marker rule | KL epithelium | 170 | **64.7** | 0.614 |
| marker rule | KL Ptprc+ | 6,223 | 4.3 | 0.044 |
| cluster | K epithelium | 20 | 40.0 | 0.269 |
| cluster | K immune | 5,477 | 2.2 | 0.027 |

In KL, Tacstd2 is ~15× higher (mean) in epithelium than immune (descriptive
log2FC 3.87). Signal on the UMAP sits in the epithelial cluster; the ~4%
immune “positives” are low-level and consistent with ambient RNA / doublets
(`Sftpc` is “on” in ~77% of all cells, including immune). Same direction in K,
with only 20 epithelial cells to work with.

By lineage in KL: epithelium 62.8%, myeloid 6.1%, T/NK 2.6%, B 2.3%,
endothelium 2.5%, fibroblast 1.4%.

## Cldn4

| definition | set | n | % Cldn4+ | mean log1p CP10k |
|---|---|---|---|---|
| cluster | KL epithelium | 231 | **10.0** (23 cells) | 0.093 |
| cluster | KL immune | 6,180 | **0.03** (2 cells) | 0.0004 |
| marker rule | KL epithelium | 170 | 10.0 | 0.099 |
| cluster | K epithelium | 20 | 5.0 (1 cell) | 0.038 |
| cluster | K immune | 5,477 | **0** | 0 |

Cldn4 is tighter to epithelium than Tacstd2 (almost no immune counts) but
**much rarer**. This is not a Cldn4-high epithelium. Companion tight-junction
genes in the same KL epithelial cluster are common: `Cldn3` 93.5%, `Cldn7`
92.6%, `Tjp1` 65.4%, `Epcam` 87.4%. So the junction program is present;
`Cldn4` itself is the sparse member.

The two KL immune Cldn4+ cells are both in the B compartment (2 / 1,406 =
0.14%). Treat as noise.

## Tacstd2–Cldn4 co-expression in epithelium

| genotype | definition | n epi | both+ | Tacstd2 only | Cldn4 only | % both | % Cldn4 among Tacstd2+ | Spearman ρ |
|---|---|---|---|---|---|---|---|---|
| KL | cluster | 231 | 20 | 125 | 3 | 8.7 | 13.8 | 0.16 (p=0.012) |
| KL | marker | 170 | 14 | 96 | 3 | 8.2 | 12.7 | 0.09 (p=0.26) |
| K | cluster | 20 | 1 | 7 | 0 | 5.0 | 12.5 | 0.43 (p=0.06) |

Most Tacstd2+ tumor cells are **Cldn4-negative**. There is no strong
Tacstd2–Cldn4 module in this series.

## KL vs K (do not use)

Nominally, KL epithelium has more Tacstd2+ cells than K (62.8% vs 40%;
descriptive log2FC 1.10; cell-level MW p=0.010) and a non-significant Cldn4
difference (10% vs 5%; p=0.47). **This is not a genotype result.** n = 1 mouse
each, genotype = sample = library, and K contributes 20 epithelial cells.
p-values treat cells as replicates (pseudoreplication).

## Honesty caveats

1. **n = 2.** One mouse per genotype. Cross-genotype numbers are confounded
   with mouse and the single pooled library.
2. **“Cold” is not reproduced as immune-fraction loss.** Both samples ~82%
   immune. Functional coldness (T-cell exclusion, suppressive myeloid) was not
   tested.
3. **Tumor cells are rare.** Dissociation is immune-dominated. K epithelium is
   20–28 cells. KL is 170–231. Epithelial statistics, especially K and Cldn4,
   are thin.
4. **Ambient surfactant.** `Sftpc`+ in ~77% of all cells. Ambient RNA can
   explain residual Tacstd2 in immune/stroma. We did not use surfactant genes
   to call epithelium.
5. **Automated marker annotation**, not expert-reviewed. Audit:
   `tables/cluster_lineage_scores.csv`. Only cluster 13 is epithelium.
6. **Author-prefiltered input.** Cannot redo cell calling or ambient
   correction.
7. **Cell-level p-values are descriptive.** They are in the tables so the
   extract is complete; they are not inferential.
8. A sibling hunt (`results/hunt_gse180963/`, Tacstd2 only) used the same
   public matrices. Our Tacstd2 epithelial vs immune numbers independently
   match that extract. **Cldn4 is the new measurement here.**

## What is not claimed

- Not claimed: KL tumors are Tacstd2-higher than K tumors.
- Not claimed: this KL sample is immune-cold by composition.
- Not claimed: a Tacstd2-high / Cldn4-high co-positive tumor state (it is
  uncommon here).
- Not claimed: protein, ADC response, or human translation.

The supported statement is: **in this public KL lung GEMM scRNA sample,
Tacstd2 is a common epithelial transcript and Cldn4 is a rare
epithelial-restricted transcript; both are near-absent from immune cells.**

## Files

Tables (`tables/`):

- `target_contrasts.csv` — Tacstd2 and Cldn4, epithelium vs immune, plus the
  confounded KL-vs-K row.
- `targets_by_compartment_genotype.csv`, `targets_by_celltype_genotype.csv`
- `tacstd2_cldn4_coexpression_epithelial.csv`
- `tj_companion_genes_epithelial.csv` — Cldn3/Cldn7/Ocln/Tjp1/Epcam in
  cluster epithelium.
- `compartment_*`, `celltype_*`, `immune_dominance_ambient_qc.csv`,
  `cluster_lineage_scores.csv`, `qc_summary.csv`

Figures (`figures/`): `umap_genotype.png`, `umap_compartment.png`,
`umap_Tacstd2.png`, `umap_Cldn4.png`,
`violin_targets_compartment_genotype.png`,
`scatter_tacstd2_vs_cldn4_epithelial.png`, `bar_compartment_fraction.png`,
`key_numbers.png`.

Machine-readable summary: `key_stats.json`. Matrices stay out of git
(`data/`); re-download with `01_download.py`.

## Reproduce

```bash
python3 scripts/w200_gse180963/01_download.py
python3 scripts/w200_gse180963/02_analyze.py
```

scanpy, anndata, leidenalg, python-igraph, pandas, scipy, numpy, matplotlib.
