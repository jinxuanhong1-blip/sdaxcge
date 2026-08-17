# CLDN4-only combo: GSE253564 + GSE148071

**This is not CellChat.** GSE253564 has no public processed scRNA, GeoMx, or bulk-with-fractions matrix. The only processed file is `GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz` (Altorki / Elemento, PMID 38401548): bulk RNA-seq FPKM, 20,187 genes × **32** pre-treatment tumors (durva001–durva061). SRA raw exists; no cell-type table. GSE148071 is scRNA (Wu et al., PMID 33953163; 42 biopsies, per-sample TXT in `GSE148071_RAW.tar`) and cannot be merged with bulk FPKM for CellChat / LIANA.

CLDN4-only. No dual-high. The two cohort rhos are **taken as given** and are not re-audited.

## Combo table — CLDN4 vs T/NK

Patient is the unit. Implied two-sided p / Fisher-z CI use the given ρ and n (Spearman *t* / *z*); they do not re-estimate ρ from either matrix.

| cohort | assay | n | CLDN4 vs T/NK ρ [95% CI] | p | I² | note |
|---|---|---:|---|---:|---:|---|
| GSE253564 | bulk FPKM, pre-treatment | 32 | −0.48 [−0.71, −0.16] | 0.00543 | — | given; not re-audited |
| GSE148071 | scRNA, partial of 42 | 25 | −0.49 [−0.74, −0.12] | 0.0129 | — | given; not re-audited |
| **RE combo** | Fisher-z DerSimonian–Laird | **57** | **−0.484 [−0.666, −0.249]** | **0.000160** | **0%** | not CellChat |

Stouffer *z* = −3.73 (p = 0.000192); Fisher combined p = 0.000741. Fixed-effect ρ = RE ρ (τ² = 0, Q = 0.002). Forest: `results/forest_CLDN4_tnk.png`. Machine table: `results/combo_cldn4_tnk.tsv`.

## Bulk chemokine / LR-proxy (GSE253564 only)

log2(FPKM+1), n=32. Pair score = mean of the two subunit logs (geometric-mean proxy of FPKM+1). BH *q* within the tested gene panel or pair list. **Bulk product proxy, not CellChat.**

### Signature scores vs CLDN4

| score | genes used | ρ vs CLDN4 | p | q |
|---|---:|---:|---:|---:|
| T-recruit receptor | 8 (CXCR3, CCR5, CCR1, CXCR6, CCR7, CXCR5, CX3CR1, CXCR4) | −0.537 | 0.00154 | 0.00617 |
| T-recruit ligand | 11 (CXCL12 absent) | −0.229 | 0.207 | 0.287 |
| IFN | 5 | −0.225 | 0.216 | 0.287 |
| checkpoint | 13 (NECTIN2→PVRL2) | −0.092 | 0.615 | 0.615 |

T-recruit receptors track the given T/NK-low direction. Ligand and IFN scores do not.

### Pair proxies with q < 0.05

| pair | pathway | ρ vs CLDN4 | p | q |
|---|---|---:|---:|---:|
| CCL21–CCR7 | T_recruit | −0.550 | 0.00111 | 0.0267 |
| CCL19–CCR7 | T_recruit | −0.506 | 0.00311 | 0.0374 |
| CXCL11–CXCR3 | T_recruit | −0.445 | 0.0108 | 0.0443 |
| CXCL13–CXCR5 | T_recruit | −0.443 | 0.0111 | 0.0443 |
| CD274–PDCD1 | checkpoint | −0.460 | 0.00801 | 0.0443 |
| PVR–CD226 | checkpoint | −0.457 | 0.00853 | 0.0443 |

CXCL12–CXCR4 dropped (CXCL12 absent). Full pair list: `results/gse253564_cldn4_lr_proxy.tsv`.

### Single genes (selected)

Negative vs CLDN4 at q < 0.05: CXCR5 −0.611, CCR7 −0.595, TIGIT −0.586, CXCR3 −0.582, PDCD1 −0.559, CXCR4 −0.535.

Opposite (higher with CLDN4, q < 0.05): CD276 +0.519, PVRL2 (NECTIN2) +0.492, CX3CL1 +0.463, PVR +0.449.

CXCL16 ligand is +0.372 (q = 0.096); CXCL9/10/11 ligands are weak once the receptor is not in the score.

T/NK marker genes are listed in `results/gse253564_cldn4_vs_genes.tsv` as documentation only and were **not** used to replace the given ρ = −0.48.

## What is not claimed

- No CellChat / LIANA / NicheNet run. One arm is bulk FPKM.
- No dual-high (TACSTD2×CLDN4) split.
- No re-audit of the two given CLDN4 vs T/NK rhos.
- Bulk pair scores mix every cell type in the biopsy; they are not sender–receiver probabilities.
