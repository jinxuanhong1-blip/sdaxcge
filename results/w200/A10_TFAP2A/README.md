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
- Public ChIP does **not** establish TFAP2A as a regulator of TACSTD2 or CLDN4
  in lung. There is **no usable lung TFAP2A ChIP-seq**. In non-lung public
  peak sets, TACSTD2 is essentially unbound and CLDN4 has no promoter-proximal
  (±1 kb) peak. Published lung TFAP2A mechanisms instead implicate PSG9,
  ITGB4 and ESR2.
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
| Upstream explanation for TACSTD2/CLDN4 | **Not supported by public ChIP** |

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

## Public ChIP check

Searched on 2026-08-16: ENCODE, ReMap 2022, and ChIP-Atlas. Overlaps used
Ensembl GRCh38 coordinates and the canonical TSS of each gene.

### Inventory

| Source | TFAP2A ChIP available | Lung? | Used? |
|---|---|---|---|
| ENCODE | ENCSR000EVP, HeLa-S3 | No | **No; revoked** |
| ReMap 2022 | MCF-7, MCF-7+E2, WA09 | No | Yes, peak overlap |
| ChIP-Atlas | 45 experiments; 0 in the Lung class | No | Yes, assembled peaks and target-gene scores |

ChIP-Atlas class counts for TFAP2A: Pluripotent stem cell 25, Digestive tract
6, Muscle 4, Epidermis 3, Others 2, Breast 1, Neural 1, Uterus 1, **Lung 0**.

There is therefore **no public lung TFAP2A ChIP-seq** to test the TACSTD2/CLDN4
pair in the relevant tissue.

### ChIP-Atlas target-gene scores (46 columns, no lung)

| Window | Gene | Non-zero experiments | Average score | Interpretation |
|---|---|---:|---:|---|
| ±1 kb | TACSTD2 | 1/46 | 2.2 | One GP5d (colon) score of 99; otherwise 0 |
| ±5 kb | TACSTD2 | 1/46 | 2.2 | Same single experiment |
| ±10 kb | TACSTD2 | 1/46 | 2.2 | Same single experiment |
| ±1 kb | CLDN4 | **0/46** | 0 | Not listed as a ±1 kb target |
| ±5 kb | CLDN4 | 7/46 | 39.0 | Distal, non-lung (colon, trophoblast, H9, 92-1) |
| ±10 kb | CLDN4 | 12/46 | 82.6 | Still no lung; HeLa/GP5d/MCF-7/HEPM/hESC |

### Peak overlap, canonical TSS

| Dataset | TACSTD2 ±1 kb | CLDN4 ±1 kb | Note |
|---|---|---|---|
| ReMap MCF-7 | No | No | CLDN4 hits start at ±5–10 kb / gene body |
| ReMap MCF-7+E2 | 1 peak | No | Isolated breast-cancer peak; not recurrent |
| ReMap WA09 | No | No | CLDN4 ±5 kb only |
| ChIP-Atlas q05 | 1 peak (GP5d) | No | Same colon experiment as the target-gene table |
| ChIP-Atlas q10 | No | No | Stricter threshold removes the TACSTD2 hit |

CLDN4's Ensembl gene span is much larger than the canonical transcript because
of alternative TSSs. Some “gene-body” hits sit several kilobases upstream of
the canonical TSS and should not be read as promoter occupancy.

### ChIP conclusion

Public ChIP does not support a TFAP2A→TACSTD2 or TFAP2A→CLDN4 model in lung.
The missing lung dataset is the decisive gap. The non-lung data that do exist
are enough to say TACSTD2 is not a recurrent TFAP2A-bound gene, and CLDN4 is
not promoter-proximal.

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
- Public ChIP is non-lung. Peak presence in MCF-7, HeLa, colon or hESC is not
  evidence of binding in lung epithelium or lung tumors.
- ChIP-Atlas assembled BEDs concatenate many experiments; q05 is more
  sensitive and noisier than q10. Target-gene scores are MACS2-derived and
  are not ChIP-qPCR validation.
- Large peak files are downloaded to `/tmp/a10_chip` at runtime and are not
  stored in the repository.

## Files

- `analyze_tcga.py` — TCGA LUAD/LUSC expression comparison.
- `analyze_chip.py` — ENCODE/ReMap/ChIP-Atlas inventory and locus overlap.
- `tcga_sample_expression.tsv` — downloaded sample-level public values.
- `tcga_gene_summary.tsv` — medians and quartiles.
- `tcga_correlations.tsv` — Spearman and log-RSEM Pearson correlations.
- `tcga_top_quartile_overlap.tsv` — high-expression overlap.
- `chip_inventory.tsv` — public TFAP2A ChIP datasets and lung status.
- `chip_locus_windows.tsv` — GRCh38 windows used for overlap.
- `chip_peak_hits.tsv` — individual overlapping peaks.
- `chip_hit_summary.tsv` — per-dataset overlap counts.
- `chip_target_gene_scores.tsv` — ChIP-Atlas target-gene scores.
- `sources.tsv` — source URLs and what each supports.

Run with:

```bash
python3 results/w200/A10_TFAP2A/analyze_tcga.py
python3 results/w200/A10_TFAP2A/analyze_chip.py
```
