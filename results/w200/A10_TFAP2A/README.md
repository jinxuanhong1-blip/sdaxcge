# A10 — TFAP2A versus TACSTD2/CLDN4 in lung

## Honest verdict

**TFAP2A is an interesting lung-cancer biology/prognosis candidate, but public
evidence does not support treating it as an interchangeable epithelial marker
for TACSTD2 (TROP2) or CLDN4.**

- For identifying epithelial lung-tumor cells, **TACSTD2/CLDN4 remain the
  stronger pair**: they are abundant membrane/junction genes and co-vary with
  each other substantially more than either co-varies with TFAP2A.
- TFAP2A's attractive feature is different: normal bulk lung RNA is very low,
  while multiple public LUAD datasets and tissue studies report tumor
  up-regulation and adverse outcome. It is a **nuclear transcription factor**,
  not a surface marker.
- The evidence found does **not** establish TFAP2A as a direct regulator of
  TACSTD2 or CLDN4 in lung. Published TFAP2A lung mechanisms instead implicate
  targets such as PSG9, ITGB4 and ESR2.
- Histology matters. In TCGA, median TFAP2A expression is about 4.7-fold higher
  in LUSC than LUAD. Calling it a uniform “lung” marker would hide that
  difference.

### Decision

| Intended use | TFAP2A assessment |
|---|---|
| General epithelial/tumor-cell annotation | **No; weaker than TACSTD2/CLDN4** |
| Cell-surface isolation, imaging or ADC target | **No; TFAP2A is nuclear** |
| LUAD tumor-vs-normal candidate | **Promising, but validate independently** |
| LUAD progression/prognosis hypothesis | **Supported, not yet clinical-grade** |
| Upstream explanation for TACSTD2/CLDN4 | **Not demonstrated** |

## Public TCGA check

The included script queried cBioPortal's public TCGA PanCancer Atlas RSEM
profiles on 2026-08-16. It retained samples with all three genes: 510 LUAD and
484 LUSC tumors.

### Abundance

| Cohort | Gene | Median RSEM | IQR |
|---|---:|---:|---:|
| LUAD | TFAP2A | 187.1 | 83.8–416.5 |
| LUAD | TACSTD2 | 6,336.0 | 4,045.8–9,543.0 |
| LUAD | CLDN4 | 9,133.6 | 5,898.9–12,945.0 |
| LUSC | TFAP2A | 874.9 | 582.8–1,321.9 |
| LUSC | TACSTD2 | 8,830.4 | 4,998.2–15,563.4 |
| LUSC | CLDN4 | 4,298.3 | 2,567.9–6,870.5 |

Cross-gene RSEM magnitudes should not be read as protein abundance or assay
sensitivity. They do, however, show that TFAP2A RNA is much less abundant in
these bulk tumors.

### Within-cohort co-expression

| Cohort | Pair | Spearman ρ | Top-quartile intersection |
|---|---|---:|---:|
| LUAD | TFAP2A–TACSTD2 | 0.241 | 49/128 |
| LUAD | TFAP2A–CLDN4 | 0.147 | 45/128 |
| LUAD | TACSTD2–CLDN4 | **0.459** | **64/128** |
| LUSC | TFAP2A–TACSTD2 | 0.141 | 39/121 |
| LUSC | TFAP2A–CLDN4 | **−0.170** | 23/121 |
| LUSC | TACSTD2–CLDN4 | **0.392** | **56/121** |

The TACSTD2–CLDN4 relationship is the positive-control benchmark. TFAP2A has
only weak positive relationships in LUAD and does not reproduce the pair in
LUSC. Bulk correlation cannot distinguish tumor-cell regulation from tumor
purity or cell-composition effects, so even the weak positive values are not
mechanistic evidence.

## Normal lung context

Current Human Protein Atlas tissue pages report consensus normal-lung RNA of:

| Gene | Consensus nTPM | Important context |
|---|---:|---|
| TFAP2A | **0.4** | Bulk lung low; HPA nevertheless identifies respiratory ciliated cells as a lung cell type with enriched expression |
| TACSTD2 | 74.5 | Expressed by lung epithelial/progenitor compartments and inducible after lung infection |
| CLDN4 | 66.7 | Epithelial tight-junction gene; HPA assigns lung cell-type enrichment to alveolar type 2 cells |

This is the best argument for TFAP2A as a tumor/normal discriminator, but it is
not proof of tumor specificity. TFAP2A is also expressed in normal squamous
epithelia and selected normal cell types. Conversely, high normal-lung bulk RNA
makes TACSTD2/CLDN4 weaker tumor-vs-normal discriminators while still making
them useful epithelial markers.

## Literature triangulation

1. A 2021 LUAD study reported elevated TFAP2A RNA in TCGA-LUAD and three GEO
   cohorts, elevated nuclear protein in a 131-pair tissue microarray, worse
   prognosis, and experimental promotion of metastasis through
   TFAP2A→PSG9→TGF-β/EMT.
2. A later NSCLC study reported higher TFAP2A in 102 paired tumors and linked
   high expression to adverse clinicopathologic features and outcome, with
   ESR2/MAPK proposed as a mechanism. This supports biological relevance, but
   it does not validate an epithelial-marker role.
3. In malignant pleural effusion scRNA-seq, CLDN4 marked an EpCAM-positive
   metastatic cluster and correlated positively with TACSTD2. This is direct
   lung single-cell support for the benchmark pair.
4. Public and experimental infection datasets show TACSTD2 induction in lung
   epithelial cells. Thus TACSTD2 positivity is not cancer-specific.
5. HPA currently lists TFAP2A, but not TACSTD2 or CLDN4, in its lung-cancer
   prognostic summary. HPA survival analyses are retrospective, cutpoint-based
   and largely TCGA-derived; they are hypothesis-generating rather than a
   clinical validation.

## Limits

- This folder contains **no private A10 experiment or w200 matrix**; none was
  present in the repository. Results here are public-data-only.
- TCGA values are bulk primary-tumor RNA. They do not establish co-expression
  in the same cell, protein localization, diagnostic sensitivity/specificity,
  causality or performance in metastases.
- LUAD and LUSC were analyzed separately; pooling them would inflate apparent
  relationships through histology differences.
- The HPA normal-lung values and TCGA tumor RSEM values come from different
  pipelines and must not be used to calculate a fold change against each
  other.
- The literature is favorable toward TFAP2A and includes reuse of public
  cohorts. Independent, blinded protein-level validation is still needed.

## Files

- `analyze_tcga.py` — standard-library reproduction script.
- `tcga_sample_expression.tsv` — downloaded sample-level public values.
- `tcga_gene_summary.tsv` — medians and quartiles.
- `tcga_correlations.tsv` — Spearman and log-RSEM Pearson correlations.
- `tcga_top_quartile_overlap.tsv` — high-expression overlap.
- `sources.tsv` — source URLs and what each supports.

Run with:

```bash
python3 results/w200/A10_TFAP2A/analyze_tcga.py
```
