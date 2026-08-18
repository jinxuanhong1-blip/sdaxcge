# FINDING — GSE154989 KP GEMM Smart-seq2: Cldn4-only

**ADDITIVE public MOUSE. Cldn4-only. Thesis already correct. No dual-high.**

Public processed matrix **exists** on GEO. Epithelial / malignant Cldn4 and IFN / MHC / TJ can be scored at mouse level. **T/NK fraction cannot.** Immune cells were FACS-excluded (`tdTomato+ / CD45− / CD11b− / TER119− / CD31−`). Residual `Cd3d` / `Nkg7` / `Ptprc` is leak, not a compartment. T/NK Spearman is **empty**. This accession does **not** join the human concordant-4 pool.

---

## Decision

| Question | Answer |
|---|---|
| Processed matrix on GEO? | **Yes** — `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5` (191 MB) + rawCount H5 + sample / gene tables |
| Usable for epithelial Cldn4? | **Yes** — `Cldn4` is in the gene table; cells are FACS epithelium / tumor |
| Usable for T/NK fraction? | **No** — CD45− FACS. Design no-go, not a missing-file no-go |
| High vs low IFN / MHC / TJ? | **Yes, mouse-level** — primary = KP animals with ≥20 cells |
| Honest primary n | KP animals ≥20 cells (not 3,891 cells; not 39 deposited `mouseID` rows) |
| Join human concordant-4? | **No** — T/NK inverse not testable; mouse additive only |
| Dual-high TACSTD2 × CLDN4 | **not defined** |

T/NK stop rule: *score epithelial/malignant Cldn4 vs T/NK fraction at mouse level*. There is no T/NK fraction. Do not invent one from leftover `Cd3d`.

---

## What the public objects actually are

Marjanovic et al., *Cancer Cell* 2020 ([PMID 32707077](https://pubmed.ncbi.nlm.nih.gov/32707077/); SuperSeries [GSE152607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152607)). Series [GSE154989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154989): “Emergence of a high-plasticity cell state during lung cancer evolution, single-cell RNAseq timecourse.”

STAR Methods: live `tdTomato+ / CD45− / CD11b− / TER119− / CD31−` cells, modified SMART-Seq2. 3,891 QC full-length transcriptomes from K / KP / T lungs at defined time points. Paper text says “39 mice”; GEO `mouseID` has **39** values because **KP 30w tumors are split** (`m1_T1`…`m3_T9`). Collapsing `_T#` gives **30 biological animals**.

| Public object | Size | Used here |
|---|---|---|
| `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5` | 190.9 Mb | **yes** — sparse COO log-like TPM (`i`,`j`,`v`) |
| `GSE154989_mmLungPlate_fQC_dSp_rawCount.h5` | 75.6 Mb | no (TPM used) |
| `GSE154989_mmLungPlate_fQC_smpTable.csv.gz` | 26.1 Kb | **yes** — `mouseID`, stage |
| `GSE154989_mmLungPlate_fQC_dZ_annot_smpTable.csv.gz` | 153.9 Kb | **yes** — `clusterK12` audit |
| `GSE154989_mmLungPlate_fQC_geneTable.csv.gz` | 762.1 Kb | **yes** — mouse symbols |
| SRA FASTQ | — | **no** — public processed only |

Sister 10x series GSE154977 (3 cisplatin libraries) was **not** opened.

---

## T/NK is a design no-go

| marker | cells with v>0 / 3891 | note |
|---|---:|---|
| Epcam | 3569 | epithelium present |
| Cldn4 | 1960 | epithelial Cldn4 present |
| Ptprc (CD45) | 122 | ~3% leak |
| Cd3d | 31 | not a T compartment |
| Nkg7 | 36 | not an NK compartment |

Do not write n = 3,891 as a T/NK n. Do not write a T/NK Spearman.

---

## Honest n

| item | n | note |
|---|---:|---|
| GEO QC cells | 3891 | do not quote as analysis n |
| deposited `mouseID` | 39 | paper’s “39 mice”; KP 30w = multiple tumors |
| biological animals (strip `_T#`) | 30 | true mouse unit |
| T / normal AT2 | 5 | not KP |
| K-only | 9 | Kras; Trp53 WT |
| KP animals | 16 | named GEMM |
| KP animals with ≥20 cells (**PRIMARY**) | *filled by* `analyze.py` | Cldn4 vs IFN/MHC/TJ |
| T/NK-fraction mice | **0** | FACS CD45− |

Primary unit = biological KP mouse with ≥20 epithelial cells. Sensitivity = K+KP ≥20 cells.

---

## Mouse-level Cldn4 vs IFN / MHC / TJ

Scores are mean depositor normTPM across A8 families mapped to mouse (MHC-I/APM = curated mouse set; TJ holds `Cldn4` out). Q4 vs Q1 and Spearman require the same locks as the human additive folders (Spearman n≥4; Q4 n_units≥8).

**Numbers are written by `analyze.py` into `tables/` and copied here after the run.**

| family | cohort | n | Spearman ρ (p) | Q4−Q1 Δ median | p |
|---|---|---:|---|---:|---:|
| IFN | KP ≥20 | — | — | — | — |
| MHC-I/APM | KP ≥20 | — | — | — | — |
| TJ (Cldn4 held out) | KP ≥20 | — | — | — | — |

IFN/MHC **down** in Cldn4-high means Q4−Q1 delta < 0 or continuous ρ < 0. That arm can still be read. It does not create a T/NK fraction.

---

## Methods (locked)

- Cldn4 only. Tacstd2 is not a gate.
- Matrix: GEO `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5` only (1-based COO `i`=gene, `j`=cell, `v`=depositor normTPM). Undetected = 0.
- Mouse unit = `mouseID` with trailing `_T#` stripped. KP 30w tumors from the same `m#` are one mouse.
- Primary cohort: genotype `KP`, n_cells ≥ 20.
- Family scores: mean of mapped mouse genes (Hallmark IFNα∪IFNγ; curated mouse MHC-I/APM; KEGG∪GOBP TJ minus Cldn4).
- T/NK fraction is **not** computed.
- No SRA. No GSE154977. No human concordant-4 merge.
- Reproduce: `python3 methods/gse154989_kp_cldn4/analyze.py`

---

## How to read this

- This is **epithelial-only Smart-seq2**, not a whole-tumor 10x atlas.
- Do not quote n=3891 cells or n=39 `mouseID` as the mouse n.
- Do not treat leftover `Cd3d` as T/NK fraction.
- Thesis is unchanged. Additive mouse IFN/MHC/TJ only.
