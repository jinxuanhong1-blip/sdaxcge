# scTCR / scBCR playbook — Caushi-class & GSE243013-class

> **METHODS ONLY.** Reproducible protocol. No results, no invented statistics.
> Bilingual: **English** first, **中文** second (mirrored). Chinese starts at
> [中文版](#中文版-sctcr--scbcr-方法手册).
>
> Tools: **scirpy** (Python; verified `0.25.x` API) and **immunarch** (R; `0.9`/`0.10`
> `rep*` API plus experimental `airr_*`).
>
> Scope: clone expansion versus **MPR**, and whether **TACSTD2-high** (and
> **CLDN4-high**) tumors have fewer **expanded CXCL13+ T cells**, on **public**
> processed files. Controlled raw FASTQ is **listed only**.

---

## 0. Dataset classes and the two estimands

### 0.1 Caushi-class

Caushi et al., *Nature* 2021 (PMID 34290408; NCT02259621, neoadjuvant nivolumab).
A Caushi-class object is **CD3-sorted** paired **scRNA-seq + scTCR-seq** from
tumor / adjacent lung / TDLN (plus optional blood), with:

- Pathologic response: **MPR** = residual viable tumor ≤ 10% at resection;
  **non-MPR** otherwise. Do not silently swap in RECIST.
- Optional **MANAFEST / ViraFEST** clone labels (MANA vs influenza vs EBV).
- Cell Ranger GEX + VDJ (or equivalent) per capture.

This class measures **T-cell** programs. It does **not** contain a malignant
epithelial compartment, so tumor `TACSTD2` / `CLDN4` cannot be scored from the
T-cell GEX matrix.

### 0.2 GSE243013-class

Liu et al., *Cell* 2025 (PMID 40147443). A GSE243013-class object is a
**CD45+ / immune-only** post-neoadjuvant chemo-immunotherapy atlas with:

- Paired 10x **5′ GEX + TCR** (not BCR) per surgical sample.
- Author cell types including `CD8T_Tex_CXCL13`, `CD4T_Tfh_CXCL13`,
  `CD4T_Th1-like_CXCL13`.
- Patient labels: `pathological_response` (`pCR` / `MPR` / `non-MPR`),
  `pathological_response_rate`, `radiological_response`, `cancer_type`
  (LUAD / LUSC), treatment fields.
- Author TCR table with `clonotype`, `expansion`, `TRA_cdr3`, `TRB_cdr3`.

Zero epithelial / malignant cells. Tumor `TACSTD2` / `CLDN4` is **not
estimable** from this release.

### 0.3 The two estimands (write these down before touching files)

| ID | Estimand | Unit | Contrast |
| --- | --- | --- | --- |
| **E1** | Abundance of **expanded CXCL13+ T cells** (CD8 Tex and/or CD4 Tfh/Th1-like) | **patient / sample** | MPR (pCR+MPR or MPR-only; pre-register) vs non-MPR |
| **E2** | Same endpoint | **patient / sample** | tumor **TACSTD2-high vs low** (optional parallel: **CLDN4-high vs low**) |

E2 is **not identified** on either class unless an **independent tumor-epithelial
TACSTD2 / CLDN4 score** exists for the same `patient_id` (paired bulk RNA, paired
epithelial scRNA, or IHC). Immune-cell or ambient `TACSTD2` is **not** a tumor
score. If that join key is missing, stop E2 and write “not estimable”.

### 0.4 What this playbook will not do

- Download or re-align **controlled FASTQ** (EGA / GSA-Human). List accessions only.
- Treat **cells as replicates** of a patient-level label (MPR, TACSTD2-high).
- Use GSE243013’s 6.6 GB immune count matrix as a source of **malignant** TACSTD2.
- Invent a TACSTD2 split when no tumor score is available.
- Equate RECIST with MPR, or pCR with MPR, without a pre-registered rule.

---

## 1. Public files vs controlled raw (list raw only)

Do **not** fetch the controlled rows. Public processed files are the analysis
inputs.

### 1.1 Caushi-class

| Role | Accession / DOI | What it actually is | Use |
| --- | --- | --- | --- |
| SuperSeries | [GSE173351](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173351) | Parent of the two GEO subseries | navigation |
| **Public sc GEX + VDJ** | [GSE176021](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176021) | Per-sample Cell Ranger GEX `.tar.gz` + `.vdj.tar.gz`; `GSE176021_CD3_annotations.rds.gz`, `GSE176021_CD8_annotations.rds.gz` | **primary public scTCR+GEX** |
| Public bulk TCR (MANAFEST cultures) | [GSE176022](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176022) | MiXCR clonotype tables (the Nature text called this “processed sc”; GEO title is bulk TCR-seq) | antigen-specific clone IDs |
| Public Adaptive TCRVβ | ImmuneACCESS [10.21417/JC2021N](https://doi.org/10.21417/JC2021N) | Bulk TCRVβ | blood / tissue repertoire |
| Scripts | [BKI-immuno/neoantigen-specific-T-cells-NSCLC](https://github.com/BKI-immuno/neoantigen-specific-T-cells-NSCLC) | figure / clustering code | reproduction, not a data host |
| **Controlled raw (LIST ONLY)** | EGA **EGAS00001005343** / dataset **EGAD00001007728** | 10x 5′ GEX+VDJ FASTQ (NovaSeq) | do not download here |

GEO sample metadata already carries `response status: MPR` / `Non-MPR` (e.g.
`GSM5352886` = `MD01-024_tumor_1`, Non-MPR).

### 1.2 GSE243013-class

| Role | File / accession | Size class | Use |
| --- | --- | --- | --- |
| **Public TCR table** | `GSE243013_T_with_TCR_annotation.csv.gz` | ~15 MB | E1 without the count matrix |
| **Public cell + clinical metadata** | `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz` | ~39 MB | MPR, histology, treatment |
| Public 10x TCR per sample | `GSE243013_RAW.tar` (GEO name is misleading) | ~0.5 GB | `filtered_contig_annotations.csv` for scirpy / immunarch |
| Public TIME subtype | `GSE243013_NMF_all_group_5.csv.gz` | 1 KB | NMF groups 1–5; **not** TACSTD2 |
| Public gene / barcode index | `GSE243013_genes.csv.gz`, `GSE243013_barcodes.csv.gz` | small | TACSTD2 is **in the gene index**; that does not make E2 estimable |
| Public immune counts | `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz` | **~6.6 GB gzipped** | optional CXCL13 RNA if author clusters are not enough; **immune-only** |
| **Controlled raw (LIST ONLY)** | NGDC GSA-Human / OMIX (GEO: “RAW data not provided”; raw FASTQ deposited at https://ngdc.cncb.ac.cn/ under Chinese Genetic Information Data Management). **No HRA/OMIX ID is printed on the GEO record** — do not invent one. | FASTQ | do not download here |

`GSE243013_RAW.tar` is **processed 10x TCR tarballs**, not FASTQ.

Verified metadata columns (do not assume others):
`sampleID, cellID, n_genes, n_genes_by_counts, total_counts, total_counts_mt,
pct_counts_mt, major_cell_type, sub_cell_type, total_counts_rb, pct_counts_rb,
gender, age, smoking_history, cancer_type, pre_treatment_staging,
anti-PD1_therapy, chemotherapy, targeted_therapy, cycles,
pathological_response, pathological_response_rate, radiological_response`.

Verified TCR-table columns:
`sampleID, cellID, sub_cell_type, TRA_v_gene, TRA_j_gene, TRA_c_gene, TRA_cdr3,
TRB_v_gene, TRB_j_gene, TRB_c_gene, TRB_cdr3, clonotype, expansion,
clonotype_number, T_new_name`.

### 1.3 Optional public companions (same methods, different n)

| Accession | Why it is listed | Caveat |
| --- | --- | --- |
| [GSE179994](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179994) | Public `GSE179994_all.scTCR.tsv.gz` + T-cell metadata; pre/post PD-1+chemo | RECIST-like response, not MPR; still sample-level |
| [GSE185204](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185204) | Paired scRNA/TCR after ICB, tumor / nLung / LN | n = 3 patients; descriptive only |
| [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) | Public scRNA + **bulk** RNA with MPR; **no TCR**. Controlled raw **HRA001033** (list only) | Use **only** if you need a tumor-epithelial TACSTD2 score for a *different* cohort; do not paste those scores onto GSE243013 patients |

### 1.4 BCR availability (honest)

| Class | scBCR? |
| --- | --- |
| Caushi-class | No. CD3 sort + TCR library. |
| GSE243013-class | No. Protocol is “5′ and **human TCR** library”. B cells are present in GEX; there is no VDJ-B table. |
| Generic 5′ + `vdj_b` | Yes. Same scirpy / immunarch path with `receptor_type == "BCR"`. |

If a public lung-ICI 5′+BCR release appears, reuse §4–§7 with IG loci
(`IGH`/`IGK`/`IGL`), clone definition on heavy+light, and SHM via
`ir.tl.mutational_load`. Until then, BCR endpoints are **not estimable** on
these two classes.

---

## 2. Analysis specification (freeze before EDA)

Write a version-controlled spec (`templates/config.template.yml`) containing:

1. **Population:** NSCLC, neoadjuvant PD-1 ± chemo, resection specimen (default).
2. **Primary endpoint (E1):** patient-level
   `n_expanded_CXCL13_CD8_clones / n_TCR_CD8_cells * 1000`
   (expanded CXCL13+ CD8 clones per 1,000 TCR-qualified CD8 cells).
3. **Primary contrast:** `path_response_bin` = (`pCR` or `MPR`) vs `non-MPR`.
   Sensitivity: `pCR` vs `non-pCR`; MPR-only vs non-MPR (drop pCR).
4. **Secondary endpoints:** expanded CXCL13+ fraction among CXCL13+ cells;
   expanded-cell fraction among all TCR+ T cells; Shannon / inverse-Simpson
   diversity; number of hyperexpanded clones (freq ≥ 1%).
5. **E2:** same primary endpoint vs tumor TACSTD2 tertiles (or median),
   **only if** a tumor score table keyed by `patient_id` is on disk.
   Parallel CLDN4 contrast uses the same rule.
6. **Stratification / covariates:** `cancer_type` (LUAD vs LUSC) is required
   on GSE243013-class. Also record `n_TCR_cells` (depth), treatment
   (chemo vs chemo-IO), and `radiological_response` as a separate label.
7. **Minimum depth:** drop samples with `< 50` TCR-qualified T cells (or
   `< 20` CD8) from E1. Report the exclusion count.
8. **Multiple testing:** one primary E1 test; secondaries descriptive or
   FDR-controlled as a block.

---

## 3. Ingestion

### 3.1 Directory layout

```text
data/public/
  caushi/GSE176021/<GSM>/{filtered_feature_bc_matrix/, vdj/filtered_contig_annotations.csv}
  gse243013/T_with_TCR_annotation.csv.gz
  gse243013/NSCLC_immune_scRNA_metadata.csv.gz
  gse243013/tcr/<sampleID>/filtered_contig_annotations.csv
  optional/tumor_tacstd2_cldn4.tsv    # patient_id, tacstd2, cldn4, source
data/controlled/   # empty placeholder; never commit FASTQ
```

### 3.2 Sample metadata (one row per sample)

Minimum columns:

`sample_id, patient_id, dataset, tissue, timepoint, cancer_type, ici_agent,
chemo, path_response, path_response_bin, residual_tumor_pct, recist,
n_gex_cells, n_tcr_cells, tacstd2_tumor, cldn4_tumor, tacstd2_source`.

`tacstd2_source` must be one of
`paired_bulk_rna | paired_epithelial_scrna | ihc | none`.
If `none`, leave scores NA and skip E2.

### 3.3 scirpy: 10x VDJ → MuData

Verified against scirpy `0.25`:

```python
import anndata as ad
import mudata as md
import scanpy as sc
import scirpy as ir

airr = ir.io.read_10x_vdj("filtered_contig_annotations.csv")  # AnnData, AIRR in obsm
airr.obs["sample_id"] = sample_id
# GEX on the same barcodes (optional for clone-only E1)
gex = sc.read_10x_mtx("filtered_feature_bc_matrix")
gex.obs["sample_id"] = sample_id
gex.layers["counts"] = gex.X.copy()
mdata = md.MuData({"gex": gex, "airr": airr})
```

Concatenate samples with `ad.concat` **before** `MuData` if you want one AIRR
modality. After `mdata.obs` is assembled, **do not strip** the
`gex:` / `airr:` prefixes — they collide (`gex:sample_id` and `airr:sample_id`).
Always select `airr:clone_id`, `airr:sample_id`, etc.

10x AIRR fields present after `read_10x_vdj`: `junction`, `junction_aa`,
`v_call`, `j_call`, `umi_count`, … There is **no** `duplicate_count` on 10x
contigs. Use `umi_count`.

### 3.4 immunarch: folder of 10x contig tables

```r
# metadata.txt in the same folder: Sample<TAB>path_response_bin<TAB>...
imm <- immunarch::repLoad("data/public/gse243013/tcr/", .mode = "paired")
```

`repLoad` auto-detects `10x (filt.contigs)` (`cdr3_nt`, `umis`, `barcode`,
`raw_clonotype_id`). Paired mode collapses TRA+TRB per barcode.

### 3.5 Author TCR table shortcut (GSE243013 E1)

You can compute E1 **without** Cell Ranger contig files:

```python
tcr = pd.read_csv("GSE243013_T_with_TCR_annotation.csv.gz")
meta = pd.read_csv("GSE243013_NSCLC_immune_scRNA_metadata.csv.gz")
# join clinical on sampleID (one row per sample from meta.drop_duplicates)
```

Use author `clonotype` as `clone_id` **or** rebuild clones from
`TRA_cdr3`+`TRB_cdr3` (preferred for a methods-controlled definition).
Author `expansion` is a useful QC check, not the primary definition.

---

## 4. Receptor QC and clone definition (scirpy)

```python
ir.pp.index_chains(mdata)          # default: productive + require_junction_aa
ir.tl.chain_qc(mdata)
# keep TCR (or BCR) with a usable pairing
keep = mdata.obs["airr:receptor_type"].isin(["TCR"]) & ~mdata.obs[
    "airr:chain_pairing"
].isin(["multichain", "orphan VDJ", "orphan VJ", "no IR"])
mdata = mdata[keep].copy()

ir.pp.ir_dist(mdata, metric="identity", sequence="nt")
ir.tl.define_clonotypes(
    mdata,
    receptor_arms="all",
    dual_ir="primary_only",
    key_added="clone_id",
)
# within-sample expansion labels: <= 1 / <= 2 / > 2
ir.tl.clonal_expansion(
    mdata,
    target_col="clone_id",
    expanded_in="sample_id",
    breakpoints=(1, 2),
    key_added="clonal_expansion",
)
```

**Default clone key (TCR):** nucleotide CDR3 identity on **both** primary
VJ and VDJ arms (`receptor_arms="all"`, `sequence="nt"`, `metric="identity"`).
This matches 10x `raw_clonotype_id` more closely than amino-acid identity.

**Sensitivity clone key:** amino-acid CDR3 + `normalized_hamming` (cutoff 15)
via `define_clonotype_clusters` → `cc_nhamming15` (convergent clusters).
Do not mix keys inside one contrast.

**BCR (when a `vdj_b` library exists):** same calls; filter
`receptor_type == "BCR"`. Prefer heavy+light (`receptor_arms="all"`).
Add `ir.tl.mutational_load` for SHM; do not use SHM as a clone key.

**Chain-pairing rules**

| Pairing | Default |
| --- | --- |
| `single pair` | keep |
| `extra VJ` / `extra VDJ` | keep primary only (`dual_ir="primary_only"`) |
| `multichain` | drop (droplet multiplet) |
| `orphan VJ` / `orphan VDJ` | drop from paired-clone E1; optional orphan-VDJ sensitivity |
| `no IR` | drop |

---

## 5. CXCL13+ and “expanded” definitions

### 5.1 CXCL13+ T cells

Use **one** primary rule, lock it in the spec:

| Rule | When to use |
| --- | --- |
| **A. Author cluster** (GSE243013-class default) | `sub_cell_type` ∈ {`CD8T_Tex_CXCL13`, `CD8T_terminal_Tex_LAYN`} for CD8-Tex; add `CD4T_Tfh_CXCL13`, `CD4T_Th1-like_CXCL13` only in a pre-registered CD4 secondary |
| **B. RNA gate** | log-normalized `CXCL13 > 0` (or UCell CXCL13-Tex signature) on TCR-qualified T cells |
| **C. Intersection** | A ∩ B, as a sensitivity |

Caushi-class has no author CXCL13 cluster name in the GEO RDS titles; use
rule B (and/or MANA-specific clones as a **separate** endpoint — MANA ≠ CXCL13).

Do **not** call a cell CXCL13+ from ambient contamination. If you use rule B
on a whole-tissue (not CD3-sorted) capture, run ambient correction on GEX
**before** the gate (CellBender / SoupX; see `methods/scrna/`).

### 5.2 Expanded clone

Primary: **clone size ≥ 2 cells in that `sample_id`** (scirpy
`clonal_expansion != "<= 1"`).

Sensitivities (report all, do not switch after seeing p-values):

- size ≥ 3;
- frequency ≥ 0.5% of TCR-qualified cells in the sample;
- immunarch `homeo` bin `Large`+`Hyperexpanded` (0.01–1).

Clone size is **within sample**, never across patients (public CDR3s are not
the same clone).

### 5.3 Patient-level endpoints

For each `sample_id` compute the table in `templates/03_patient_endpoints.py`:

| Column | Definition |
| --- | --- |
| `n_tcr` | TCR-qualified T cells |
| `n_cd8_tcr` | subset CD8 (author label or `CD8A`/`CD8B` RNA) |
| `n_clones` | unique `clone_id` |
| `n_exp_clones` | clones with size ≥ 2 |
| `frac_exp_cells` | cells in expanded clones / `n_tcr` |
| `n_cxcl13` | CXCL13+ T cells (rule A or B) |
| `n_exp_cxcl13_cells` | CXCL13+ **and** expanded |
| `n_exp_cxcl13_clones` | unique clones that are expanded and have ≥ 1 CXCL13+ cell |
| `exp_cxcl13_per_k` | `n_exp_cxcl13_clones / n_cd8_tcr * 1000` (**E1 primary**) |
| `frac_cxcl13_that_are_exp` | `n_exp_cxcl13_cells / n_cxcl13` |

Join `path_response_bin`, `cancer_type`, and (if present) `tacstd2_tumor`.

---

## 6. immunarch clonality / diversity (same samples)

```r
imm <- immunarch::repLoad("data/public/.../tcr/")
# volume / count
immunarch::repExplore(imm$data, .method = "volume")
# occupied repertoire space
immunarch::repClonality(imm$data, .method = "homeo")
immunarch::repClonality(imm$data, .method = "clonal.prop", .perc = 10)
# diversity (never interpret without depth)
immunarch::repDiversity(imm$data, .method = "chao1")
immunarch::repDiversity(imm$data, .method = "inv.simp")
immunarch::repDiversity(imm$data, .method = "d50")
# overlap (same patient, tumor vs nLung) — Caushi-class
immunarch::repOverlap(imm$data, .method = "jaccard", .col = "aa")
```

Subset to CXCL13+ barcodes with `select_barcodes()` (needs a `Barcode` column
from 10x load), then recompute `repClonality` / `repDiversity` on that subset.
Do **not** compare raw Chao1 across samples with 10× different `n_tcr`;
downsample (`repSample`) to the minimum passing depth as a sensitivity.

Experimental immunarch 1.0 aliases (`airr_clonality_prop`,
`airr_diversity_shannon`) may be used **in addition to**, not instead of, the
`rep*` calls above until you pin a single major version.

---

## 7. Statistics for E1 and E2

**Unit = patient.** If two captures exist per patient, collapse by summing
cells / clones or by keeping the tumor capture only (pre-register).

**E1 test (default):** two-sided Wilcoxon rank-sum of `exp_cxcl13_per_k`
between `path_response_bin` groups, **stratified by `cancer_type`**
(LUAD and LUSC separately, plus a van Elteren / blocked rank test).
Linear model sensitivity:

```text
exp_cxcl13_per_k ~ path_response_bin + cancer_type + log1p(n_cd8_tcr)
```

**E2 test (only if tumor scores exist):** same model with
`tacstd2_bin` (high = top tertile vs bottom tertile; drop middle) and a
Spearman correlation of `tacstd2_tumor` vs `exp_cxcl13_per_k` within
histology. Repeat for `cldn4_tumor`. Do **not** residualize TACSTD2 on
purity here unless a purity column is in the **same** tumor-score table
(see `methods/scrna/` and the ICI-purity playbook).

**Forbidden:** Wilcoxon on cells labeled MPR vs non-MPR; pseudoreplication
by splitting one clone across tests; using NMF TIME group as a TACSTD2
surrogate.

**Missing E2:** write a one-line skip:
“E2 not estimable: no paired tumor-epithelial TACSTD2/CLDN4 table for
these `patient_id`s.”

---

## 8. TACSTD2 / CLDN4 — how to (and how not to) get a tumor score

| Source | Allowed for E2? |
| --- | --- |
| Paired **bulk** tumor RNA (`TACSTD2`, `CLDN4` log-TPM) on the same patient | Yes |
| Paired **epithelial / malignant** scRNA mean | Yes |
| Protein / IHC H-score on the same resection | Yes |
| Immune-cell `TACSTD2` from GSE243013 counts | **No** (wrong compartment; ambient risk) |
| Caushi T-cell GEX `TACSTD2` | **No** |
| TCGA-LUAD TACSTD2 tertiles pasted onto unrelated ICI patients | **No** |
| GSE207422 malignant TACSTD2 (different trial) | **No** as a join; yes as a **separate** A3-style analysis |
| `GSE243013_NMF_all_group_5.csv` | TIME subtype, **not** TACSTD2 |

Ambient `TACSTD2`/`CLDN4` from lysed epithelium is a known soup gene. Any
RNA gate on non-epithelial barcodes must come after ambient removal and
still must not be used as a tumor score.

---

## 9. Recommended run order

1. Freeze `templates/config.template.yml`.
2. Download **only** public files in §1 (TCR table + metadata first for
   GSE243013; GSE176021 GEX+VDJ tarballs for Caushi).
3. Build `sample_metadata.tsv`. Set `tacstd2_source=none` unless a real
   tumor table is present.
4. scirpy: index → chain QC → `ir_dist` → `define_clonotypes` →
   `clonal_expansion` (`templates/01_scirpy_clonotypes.py`).
5. immunarch: `repLoad` → clonality / diversity
   (`templates/02_immunarch_clonality.R`).
6. Patient endpoints + E1 test (`templates/03_patient_endpoints.py`).
7. E2 only after a tumor-score join. Otherwise skip.
8. Sensitivities: clone key (nt vs aa-cluster), expansion threshold,
   CXCL13 rule A vs B, pCR vs MPR grouping, LUAD-only / LUSC-only.

---

## 10. Outputs (methods artifacts, not findings)

Write under a local work directory (do not commit patient-level clinical
tables if they are identifiable):

- `qc_chain_pairing.tsv` — counts of pairing categories per sample
- `clone_table.tsv` — `sample_id, clone_id, clone_size, n_cxcl13, v_call, j_call, junction_aa`
- `patient_endpoints.tsv` — §5.3 columns + labels
- `e1_spec.json` — frozen contrast and exclusion n
- `e2_status.txt` — `estimable` + source, or `not estimable`

---

## 11. Environment

See `env/environment.yml`. Pins used to validate the templates:

- Python 3.11+: `scirpy>=0.25,<0.26`, `scanpy>=1.10`, `anndata>=0.10`,
  `mudata>=0.3`, `pandas`, `scipy`
- R 4.2+: `immunarch>=0.10` (CRAN; pulls `immundata`)

---

# 中文版 · scTCR / scBCR 方法手册

> **仅方法。** 无结果、不编造统计量。受控原始 FASTQ **只列登记号、不下载**。
> 工具：scirpy、immunarch。科学问题：克隆扩增 vs **MPR**；
> **TACSTD2 高**（及 CLDN4 高）肿瘤是否更少 **扩增的 CXCL13+ T 细胞**。

---

## 0. 数据类别与两个估计量

### 0.1 Caushi 类

Caushi 等，*Nature* 2021（PMID 34290408；NCT02259621，新辅助纳武利尤单抗）。
**CD3 分选**的配对 scRNA + scTCR（肿瘤 / 癌旁 / 引流淋巴结，可含外周血）：

- 病理缓解：**MPR** = 手术时残存活肿瘤 ≤ 10%；否则 **non-MPR**。不可与 RECIST 混用。
- 可选 MANAFEST / ViraFEST 抗原标签（MANA / 流感 / EBV）。
- 每捕获一份 Cell Ranger GEX + VDJ。

本类只测 **T 细胞**，**没有**恶性上皮区室，不能从 T 细胞表达矩阵给肿瘤打
`TACSTD2` / `CLDN4` 分。

### 0.2 GSE243013 类

Liu 等，*Cell* 2025（PMID 40147443）。新辅助化疗免疫后的 **CD45+ / 仅免疫** 图谱：

- 每例手术样本 10x **5′ GEX + TCR**（无 BCR）。
- 作者分群含 `CD8T_Tex_CXCL13`、`CD4T_Tfh_CXCL13`、`CD4T_Th1-like_CXCL13`。
- 患者标签：`pathological_response`（pCR / MPR / non-MPR）、病理缓解率、影像学缓解、
  `cancer_type`（LUAD / LUSC）、治疗字段。
- 作者 TCR 表含 `clonotype`、`expansion`、TRA/TRB CDR3。

无上皮 / 恶性细胞。本释放 **不能**估计肿瘤 TACSTD2 / CLDN4。

### 0.3 两个估计量（先写后做）

| 编号 | 估计量 | 分析单位 | 对比 |
| --- | --- | --- | --- |
| **E1** | **扩增的 CXCL13+ T 细胞**丰度（CD8 Tex 和/或 CD4 Tfh） | **患者 / 样本** | MPR（预注册 pCR+MPR 或仅 MPR）vs non-MPR |
| **E2** | 同一终点 | **患者 / 样本** | 肿瘤 **TACSTD2 高 vs 低**（平行：CLDN4） |

E2 仅当同一 `patient_id` 存在**独立的肿瘤上皮** TACSTD2 / CLDN4 评分
（配对 bulk RNA、配对上皮 scRNA 或 IHC）时才可识别。免疫细胞或环境游离
`TACSTD2` **不是**肿瘤评分。缺 join 键则停止 E2，写明“不可估计”。

### 0.4 明确不做

- 下载或重比对受控 FASTQ（EGA / GSA-Human）。
- 把细胞当作患者级标签（MPR、TACSTD2-high）的重复。
- 用 GSE243013 的 6.6 GB 免疫矩阵当作**恶性** TACSTD2。
- 没有肿瘤评分时编造 TACSTD2 分组。
- 未预注册就把 RECIST 与 MPR、pCR 与 MPR 等同。

---

## 1. 公共文件 vs 受控原始（原始只列）

### 1.1 Caushi 类

| 角色 | 登记号 | 实际内容 | 用法 |
| --- | --- | --- | --- |
| SuperSeries | GSE173351 | 两个子系列的父项 | 导航 |
| **公共 sc GEX+VDJ** | **GSE176021** | 每样本 Cell Ranger GEX + `.vdj.tar.gz`；CD3/CD8 注释 RDS | **主分析输入** |
| 公共 bulk TCR | GSE176022 | MiXCR 克隆型表（Nature 正文曾写“processed sc”；GEO 标题为 bulk TCR） | 抗原特异性克隆 ID |
| Adaptive TCRVβ | ImmuneACCESS 10.21417/JC2021N | 批量 TCRVβ | 血 / 组织库 |
| **受控原始（只列）** | EGA **EGAS00001005343** / **EGAD00001007728** | 10x 5′ GEX+VDJ FASTQ | 此处不下载 |

### 1.2 GSE243013 类

| 角色 | 文件 | 用法 |
| --- | --- | --- |
| **公共 TCR 表** | `GSE243013_T_with_TCR_annotation.csv.gz`（~15 MB） | 不做表达矩阵即可做 E1 |
| **公共临床+细胞元数据** | `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz`（~39 MB） | MPR、组织学 |
| 公共 10x TCR | `GSE243013_RAW.tar`（~0.5 GB；**不是 FASTQ**） | 每样本 `filtered_contig_annotations.csv` |
| TIME 亚型 | `GSE243013_NMF_all_group_5.csv.gz` | NMF 1–5；**不是** TACSTD2 |
| 免疫 counts | `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz`（**约 6.6 GB**） | 可选 CXCL13 RNA；仍是免疫 |
| **受控原始（只列）** | NGDC GSA-Human / OMIX（GEO 写明 raw 不在 GEO，按中国人类遗传资源管理上传 ngdc.cncb.ac.cn）。**GEO 未印刷 HRA/OMIX 号，不得编造。** | 不下载 |

已核对元数据列与 TCR 表列见英文 §1.2。

### 1.3 可选公共伴随队列

GSE179994（公共 scTCR，偏 RECIST，非 MPR）；GSE185204（n=3，仅描述）；
GSE207422（有 MPR 与恶性 TACSTD2，**无 TCR**；受控原始 **HRA001033** 只列）。
GSE207422 的 TACSTD2 **不能**贴到 GSE243013 患者上。

### 1.4 BCR（如实）

Caushi 类与 GSE243013 类均**无 scBCR**。仅当存在 5′ + `vdj_b` 库时，用同一
scirpy/immunarch 路径（`receptor_type=="BCR"`，重+轻链，`mutational_load` 测 SHM）。
在此之前 BCR 终点**不可估计**。

---

## 2. 分析方案（探索前冻结）

见 `templates/config.template.yml`。

1. 人群：NSCLC，新辅助 PD-1 ± 化疗，手术标本。
2. **E1 主终点：** 每 1,000 个 TCR 合格 CD8 中的扩增 CXCL13+ CD8 克隆数
   `n_exp_cxcl13_clones / n_cd8_tcr * 1000`。
3. **主对比：** `pCR` 或 `MPR` vs `non-MPR`（敏感性：仅 pCR；去掉 pCR 的 MPR vs non-MPR）。
4. 次要：CXCL13+ 中扩增细胞比例、全体 T 中扩增细胞比例、Shannon / 逆 Simpson、
   高频克隆（≥1%）数。
5. **E2：** 仅当磁盘上有按 `patient_id` 索引的肿瘤 TACSTD2 表时，用三分位
   （或中位数）做同样终点；CLDN4 平行、同规则。
6. **分层：** GSE243013 类必须按 LUAD / LUSC；记录 TCR 深度、化疗 vs 化疗免疫、
   影像学缓解（单独标签）。
7. **最小深度：** TCR 合格 T `< 50`（或 CD8 `< 20`）的样本从 E1 排除并报告例数。
8. 多重检验：只设一个 E1 主检验；次要描述或按块 FDR。

---

## 3. 读入

目录见英文 §3.1。样本表最少字段见 §3.2。`tacstd2_source` 只能是
`paired_bulk_rna | paired_epithelial_scrna | ihc | none`。

**scirpy：** `ir.io.read_10x_vdj` → 可选 GEX → `MuData({"gex","airr"})`。
`mdata.obs` 带 `gex:` / `airr:` 前缀，**不要去掉前缀**（会重名）。
10x 字段是 `umi_count`，没有 `duplicate_count`。

**immunarch：** `repLoad(..., .mode="paired")` 自动识别
`10x (filt.contigs)`。

**GSE243013 捷径：** 只读 TCR 注释表 + 元数据即可做 E1。作者 `clonotype` 作 QC；
方法学主定义用 `TRA_cdr3`+`TRB_cdr3` 重算。

---

## 4. 受体质控与克隆定义（scirpy）

`index_chains` → `chain_qc` → 去掉 `multichain` / orphan / `no IR` →
`ir_dist(identity, nt)` → `define_clonotypes(receptor_arms="all",
dual_ir="primary_only")` → `clonal_expansion(expanded_in="sample_id")`。

**默认克隆键：** 两条主链核苷酸 CDR3 完全一致。
**敏感性：** 氨基酸 + `normalized_hamming`（cutoff 15）的
`define_clonotype_clusters`。一次对比只用一种键。

BCR（若有 `vdj_b`）：滤 `receptor_type=="BCR"`，重+轻链；SHM 用
`mutational_load`，不用 SHM 当克隆键。

配对规则：`single pair` 保留；额外 VJ/VDJ 只留 primary；`multichain` 与
orphan、`no IR` 默认剔除。

---

## 5. CXCL13+ 与“扩增”

**CXCL13+ 只选一条主规则并写入方案：**

- **A.** 作者分群（GSE243013 默认）：CD8 用 `CD8T_Tex_CXCL13`（可加
  `CD8T_terminal_Tex_LAYN`）；CD4 Tfh/Th1-like 仅作预注册次要。
- **B.** RNA 门控：标准化后 `CXCL13 > 0` 或 UCell 签名。
- **C.** A∩B 作敏感性。

Caushi 类无作者 CXCL13 分群名：用规则 B。MANA 特异性克隆是**另一终点**，≠ CXCL13。
全组织捕获做规则 B 前须先做环境 RNA 校正。

**扩增（主）：** 该 `sample_id` 内克隆大小 ≥ 2。
敏感性：≥ 3；频率 ≥ 0.5%；immunarch `Large`+`Hyperexpanded`。
克隆大小**不跨患者**。

患者级列见英文 §5.3；E1 主列是 `exp_cxcl13_per_k`。

---

## 6. immunarch 克隆性 / 多样性

`repExplore` / `repClonality(homeo, clonal.prop)` /
`repDiversity(chao1, inv.simp, d50)` / `repOverlap(jaccard)`。
CXCL13+ 子集用 `select_barcodes`。深度差一个数量级时不可直接比 Chao1，
用 `repSample` 下采样作敏感性。immunarch 1.0 的 `airr_*` 只能作附加，
不能在未钉版本时代替 `rep*`。

---

## 7. E1 / E2 统计

**单位 = 患者。** 同一患者两次捕获：预注册为只留肿瘤捕获或合并细胞/克隆。

**E1 默认：** `exp_cxcl13_per_k` 在 `path_response_bin` 间 Wilcoxon，
**按组织学分层**（LUAD、LUSC 分开 + 分层秩和）。
线性敏感性：`~ path_response_bin + cancer_type + log1p(n_cd8_tcr)`。

**E2（仅当有肿瘤评分）：** 上三分位 vs 下三分位；组内 Spearman。
CLDN4 同规则。没有纯度列不要在这里对 TACSTD2 做纯度残差。

**禁止：** 以细胞为重复的 MPR Wilcoxon；把 NMF TIME 组当作 TACSTD2。

**E2 缺失：** 写一行
“E2 不可估计：这些 `patient_id` 没有配对的肿瘤上皮 TACSTD2/CLDN4 表。”

---

## 8. 肿瘤 TACSTD2 / CLDN4 评分：能用与不能用

| 来源 | 可否用于 E2 |
| --- | --- |
| 同一患者配对肿瘤 bulk RNA | 可 |
| 配对上皮 / 恶性 scRNA 均值 | 可 |
| 同一切除标本 IHC | 可 |
| GSE243013 免疫细胞 TACSTD2 | **否** |
| Caushi T 细胞 GEX TACSTD2 | **否** |
| 把 TCGA 或 GSE207422 的 TACSTD2 贴到另一队列患者 | **否** |
| NMF 五分组 | TIME 亚型，**不是** TACSTD2 |

---

## 9. 推荐顺序

冻结配置 → 只下 §1 公共文件 → 建样本表（无肿瘤表则 `tacstd2_source=none`）→
scirpy 克隆 → immunarch 克隆性 → 患者终点与 E1 → 有肿瘤评分才做 E2 →
敏感性（克隆键、扩增阈值、CXCL13 规则、pCR/MPR、组织学）。

---

## 10–11. 产物与环境

产物：`qc_chain_pairing.tsv`、`clone_table.tsv`、`patient_endpoints.tsv`、
`e1_spec.json`、`e2_status.txt`。环境见 `env/environment.yml`。
