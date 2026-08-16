# Methods — CPTAC LUAD TACSTD2 protein vs MCP-counter / ESTIMATE / CIBERSORT

Self-contained rework. Outputs only under `results/rework/CPTAC_LUAD_defs/`.

## Question
PR6 found LUAD **TACSTD2 protein** vs freeze **xCell immune score** Spearman
ρ = −0.309 (p = 0.001, primary-list q = 0.024). LSCC TACSTD2 protein vs
xCell was **null** (PR23 / PR84). This folder asks whether the LUAD
association is definition-specific: does it hold for **public ESTIMATE**,
**public CIBERSORT**, and **MCP-counter** after a **WES purity residual**?

This is a treatment-naive surgical cohort (Gillette et al. *Cell* 2020,
PMID 32649874). **No ICI labels.** OS/PFS are not retested here.

## Data (open S3 freeze v1.2)
Prefix: `https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/`

| File | Role | HEAD |
|---|---|---|
| `LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt` | TMT protein | 200 |
| `LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt` | RNA for MCP-counter | 200 |
| `LUAD_phenotype.txt` | ESTIMATE, CIBERSORT, xCell, WES/WGS purity | 200 |
| `LUAD_meta.txt` | Age / Sex / Stage | 200 |
| `LUAD_{xcell,cibersort,mcpcounter,estimate}.txt` | standalone | **403** |

Targets: TACSTD2 `ENSG00000184292`, CLDN4 `ENSG00000189143` (version-stripped match).
n = 110 tumors with protein + RNA + phenotype IDs aligned.

## Immune definitions
**Primary (rework FDR family, TACSTD2 protein only):**

1. `ESTIMATE_ImmuneScore` — freeze column (Yoshihara 2013).
2. `MCPcounter_T_cells` — computed here (Becht et al. *Genome Biol* 2016).
3. `MCPcounter_Cytotoxic_lymphocytes` — computed here.
4. `CIBERSORT_T_cell_CD8+` — freeze column (immunedeconv CIBERSORT).

**Reproduction (not in the rework FDR family):** `xCell_immune_score` and
`xCell_T_cell_CD8+` must recover PR6 (ρ ≈ −0.309 / −0.289).

**Exploratory:** other MCP-counter populations, CIBERSORT Treg / M2 /
neutrophil, CIBERSORT leukocyte sum, ESTIMATE StromalScore, CLDN4 protein.

### Why MCP-counter is computed
MCP-counter is **not** a column in `LUAD_phenotype.txt`. Standalone
`LUAD_mcpcounter.txt` is HTTP 403. Scores are computed from the **public**
tumor RNA matrix using the official Becht 2016 HUGO/Ensembl marker table
(`scripts/rework/CPTAC_LUAD_defs/mcp_counter_genes.tsv`, from
https://github.com/ebecht/MCPcounter/Signatures/genes.txt).

Freeze RNA is `log2(RSEM coding UQ 1500)`. Official MCP-counter is the
geometric mean of linear expression. The arithmetic mean of log2 values
**is** that transform. Gene coverage is written to
`tables/mcp_counter_coverage.tsv`. Official “CD8 T cells” is a **single
gene** (`CD8B`) and is exploratory, not primary.

### CIBERSORT caveat
The 22 freeze CIBERSORT columns do **not** sum to 1 (median sum ≈ 1.8).
They are not relative-mode fractions. `CIBERSORT_T_cell_CD8+` is used as
published. Leukocyte sum is exploratory only.

## Statistics
- Marginal Spearman (pairwise complete).
- Primary residual: **WES_purity** (DNA; 2/110 missing).
- Partial Spearman = Pearson of rank-residuals after OLS on ranked purity
  (df = n−3). Fisher-z 95% CI uses n−4.
- Sensitivity: Spearman of OLS residuals of the raw values; WGS_purity;
  ESTIMATE cosine tumor purity (Yoshihara 2013). ESTIMATE purity is
  **circular** with ImmuneScore and is not an honest residual for that score.
- BH-FDR across the **four** TACSTD2 × primary-score WES partial tests only.
  CLDN4 gets its own 4-test FDR (comparator; TMT NA = 31/110).
- No ICI endpoint. No fabricated accessions.

## How to rerun
```bash
pip install -r scripts/rework/CPTAC_LUAD_defs/requirements.txt
python3 scripts/rework/CPTAC_LUAD_defs/download.py
python3 scripts/rework/CPTAC_LUAD_defs/analyze.py
```
