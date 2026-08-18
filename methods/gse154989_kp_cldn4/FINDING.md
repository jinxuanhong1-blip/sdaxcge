# FINDING — GSE154989 KP GEMM Smart-seq2: Cldn4-only

**ADDITIVE public MOUSE. Cldn4-only. Thesis already correct. No dual-high.**

Public processed matrix **exists** on GEO. Epithelial Cldn4 and IFN / MHC / TJ were scored at **biological mouse** level. **T/NK fraction cannot be scored.** Cells were FACS `tdTomato+ / CD45− / CD11b− / TER119− / CD31−`. Residual `Cd3d` / `Nkg7` / `Ptprc` is leak, not a compartment. T/NK Spearman is **empty**. This accession does **not** join the human concordant-4 pool.

Primary KP (n=15 mice, ≥20 cells): IFN and MHC-I/APM are **not** down in Cldn4-high (both ρ>0, both Q4−Q1 Δ>0, both p>0.2). TJ (Cldn4 held out) tracks Cldn4 (ρ=0.757, p=0.001) — sanity that Cldn4 is a real tight-junction gene here, not a join key.

---

## Decision

| Question | Answer |
|---|---|
| Processed matrix on GEO? | **Yes** — `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5` |
| Usable for epithelial Cldn4? | **Yes** — 1,960 / 3,891 cells Cldn4+; 3,569 Epcam+ |
| Usable for T/NK fraction? | **No** — CD45− FACS. Design no-go, not a missing-file no-go |
| High vs low IFN / MHC / TJ? | **Yes, mouse-level** — primary KP n=**15** |
| IFN/MHC down in Cldn4-high (KP)? | **No** — IFN ρ=+0.339 p=0.216; MHC ρ=+0.139 p=0.621 |
| Honest primary n | **15 KP mice** (not 3,891 cells; not 39 deposited `mouseID`) |
| Join human concordant-4? | **No** — T/NK inverse not testable; KP IFN/MHC not down |
| Dual-high TACSTD2 × CLDN4 | **not defined** |

T/NK stop rule: *score epithelial/malignant Cldn4 vs T/NK fraction at mouse level*. There is no T/NK fraction. Do not invent one from leftover `Cd3d`.

---

## What the public objects actually are

Marjanovic et al., *Cancer Cell* 2020 ([PMID 32707077](https://pubmed.ncbi.nlm.nih.gov/32707077/); SuperSeries [GSE152607](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152607)). Series [GSE154989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154989): “Emergence of a high-plasticity cell state during lung cancer evolution, single-cell RNAseq timecourse.”

STAR Methods: live `tdTomato+ / CD45− / CD11b− / TER119− / CD31−` cells, modified SMART-Seq2. 3,891 QC full-length transcriptomes from K / KP / T lungs at defined time points. Paper text says “39 mice”; GEO `mouseID` has **39** values because **KP 30w tumors are split** (`m1_T1`…`m3_T9`). Collapsing `_T#` gives **30 biological animals**.

| Public object | Size | Used here |
|---|---|---|
| `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5` | 190.9 Mb | **yes** — sparse COO (`i`=gene 1…52638, `j`=cell 1…3891, `v`=depositor normTPM) |
| `GSE154989_mmLungPlate_fQC_dSp_rawCount.h5` | 75.6 Mb | no (TPM used) |
| `GSE154989_mmLungPlate_fQC_smpTable.csv.gz` | 26.1 Kb | **yes** — `mouseID`, stage |
| `GSE154989_mmLungPlate_fQC_dZ_annot_smpTable.csv.gz` | 153.9 Kb | **yes** — `clusterK12` audit |
| `GSE154989_mmLungPlate_fQC_geneTable.csv.gz` | 762.1 Kb | **yes** — mouse symbols (`Cldn4` present) |
| SRA FASTQ | — | **no** — public processed only |

Sister 10x series GSE154977 (3 cisplatin libraries) was **not** opened.

H5 is not 10x. It is a custom COO dump. Undetected genes are absent from the triplet and are scored as 0.

---

## T/NK is a design no-go

| marker | cells with v>0 / 3891 | note |
|---|---:|---|
| Epcam | 3569 | epithelium present |
| Cldn4 | 1960 | epithelial Cldn4 present |
| Ptprc (CD45) | 122 | ~3% leak |
| Cd3d | 31 | not a T compartment |
| Nkg7 | 36 | not an NK compartment |
| any T/NK audit gene | 202 | still leak; two KP 12w plates are dirtier |

`clusterK12` is 12 epithelial / tumor states. Ptprc+ rate is ≤8% in every cluster. Cluster 11 (89 cells) is Epcam-low; it is not a T/NK cluster and was not used to mint a fraction.

Do not write n=3,891 as a T/NK n. Do not write a T/NK Spearman.

---

## Honest n

| item | n | note |
|---|---:|---|
| GEO QC cells | 3891 | do not quote as analysis n |
| deposited `mouseID` | 39 | paper’s “39 mice”; KP 30w = 12 tumors from 3 mice |
| biological animals (strip `_T#`) | 30 | true mouse unit |
| T / normal AT2 | 5 | not KP; Cldn4-low, MHC-high |
| K-only | 9 | Kras; Trp53 WT |
| KP animals | 16 | named GEMM |
| KP animals with ≥20 cells (**PRIMARY**) | **15** | dropped `KP_2w_ND_m1` (4 cells) |
| K+KP animals with ≥20 cells | 24 | sensitivity only |
| T/NK-fraction mice | **0** | FACS CD45− |

Primary unit = biological KP mouse with ≥20 epithelial cells. Q4 vs Q1 tails are **4 vs 4**. That is honest and thin.

---

## Mouse-level Cldn4 vs IFN / MHC / TJ

Family score = mean depositor normTPM of mapped mouse genes. IFN = Hallmark IFNα∪IFNγ (217 / ~250 human genes mapped and present). MHC-I/APM = curated mouse set (24 genes: `H2-K1`, `H2-D1`, `H2-Q*`, `H2-T23`, `B2m`, TAP/PSMB/ERAP machinery). TJ = KEGG∪GOBP minus `Cldn4` (200 genes).

### Primary: KP, n=15

| family | n | n_Q1 / n_Q4 | Spearman ρ (p) | Q4−Q1 Δ median | rank-biserial r | MW p |
|---|---:|---|---|---:|---:|---:|
| IFN | 15 | 4 / 4 | +0.339 (0.216) | +0.228 | +0.625 | 0.200 |
| MHC-I/APM | 15 | 4 / 4 | +0.139 (0.621) | +0.193 | +0.250 | 0.686 |
| TJ (Cldn4 held out) | 15 | 4 / 4 | **+0.757 (0.001)** | +0.550 | +1.000 | 0.029 |

IFN/MHC **down** would be Δ<0 or ρ<0. Primary KP is the opposite sign, not significant. TJ-up is the expected Cldn4–junction correlation, not immune exclusion.

### Per-mouse primary units

| mouse | week | n_cells | n_tumors | Cldn4 mean | Cldn4 %pos | IFN | MHC-I/APM | TJ |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| KP_2w_ND_m2 | 2w | 27 | 1 | 1.622 | 0.296 | 1.491 | 3.340 | 1.961 |
| KP_12w_ND_m1 | 12w | 131 | 1 | 4.012 | 0.687 | 1.546 | 2.514 | 2.423 |
| KP_12w_ND_m2 | 12w | 124 | 1 | 4.942 | 0.758 | 1.808 | 3.073 | 2.600 |
| KP_12w_ND_m3 | 12w | 74 | 1 | 4.575 | 0.770 | 1.591 | 2.555 | 2.476 |
| KP_12w_ND_m4 | 12w | 62 | 1 | 3.713 | 0.435 | 1.618 | 3.091 | 2.218 |
| KP_12w_ND_m5 | 12w | 48 | 1 | 4.506 | 0.667 | 1.322 | 2.282 | 2.175 |
| KP_12w_ND_m6 | 12w | 52 | 1 | 3.452 | 0.519 | 1.535 | 2.924 | 2.154 |
| KP_20w_ND_m2 | 20w | 43 | 1 | 4.556 | 0.744 | 1.625 | 3.000 | 2.275 |
| KP_20w_ND_m3 | 20w | 60 | 1 | 5.958 | 0.817 | 1.715 | 3.146 | 2.397 |
| KP_20w_ND_m4 | 20w | 173 | 1 | 3.014 | 0.561 | 1.614 | 2.602 | 2.326 |
| KP_20w_ND_m5 | 20w | 126 | 1 | 1.495 | 0.317 | 1.359 | 2.183 | 1.750 |
| KP_20w_ND_m6 | 20w | 136 | 1 | 1.598 | 0.331 | 1.768 | 3.059 | 2.205 |
| KP_30w_ND_m1 | 30w | 520 | 4 | 2.945 | 0.448 | 1.176 | 1.910 | 1.811 |
| KP_30w_ND_m2 | 30w | 619 | 5 | 5.226 | 0.775 | 1.500 | 2.369 | 2.381 |
| KP_30w_ND_m3 | 30w | 415 | 3 | 4.392 | 0.667 | 1.609 | 2.673 | 2.386 |

Honest Spearman n = **15** (not 2,610 KP cells, not 12 KP-30w tumors, not 39 `mouseID`).

### Sensitivity (not the claim)

| contrast | family | n | ρ (p) | Q4−Q1 Δ | note |
|---|---|---:|---|---:|---|
| K+KP ≥20 | IFN | 24 | +0.270 (0.203) | +0.069 | not down |
| K+KP ≥20 | MHC-I/APM | 24 | −0.288 (0.173) | −0.171 | ns; mixes K |
| K+KP ≥20 | TJ | 24 | +0.737 (4×10⁻⁵) | +0.382 | same TJ sanity |
| all animals ≥20 (incl. T AT2) | MHC-I/APM | 28 | −0.500 (0.007) | −0.665 | **confounded** — normal AT2 are Cldn4-low / MHC-high |
| KP including 4-cell mouse | IFN | 16 | +0.103 (0.704) | +0.024 | still not down |
| KP 30w tumors as units | IFN | 12 | +0.566 (0.055) | +0.419 | not mice; same-mouse tumors |

Do **not** quote the all-animal MHC-down as support. That sign is T AT2 vs tumor, not Cldn4 vs IFN inside KP.

---

## Methods (locked)

- Cldn4 only. Tacstd2 is not a gate.
- Matrix: GEO `GSE154989_mmLungPlate_fQC_dSp_normTPM.h5` only. Undetected = 0.
- Mouse unit = `mouseID` with trailing `_T#` stripped. KP 30w tumors from the same `m#` are one mouse.
- Primary cohort: genotype `KP`, n_cells ≥ 20. Dropped `KP_2w_ND_m1` (4 cells).
- Family scores: mean of mapped mouse genes (Hallmark IFNα∪IFNγ; curated mouse MHC-I/APM; KEGG∪GOBP TJ minus Cldn4).
- Spearman requires n≥4 finite pairs. Q4 vs Q1 requires n_units≥8 and two tails.
- T/NK fraction is **not** computed.
- No SRA. No GSE154977. No human concordant-4 merge.
- Reproduce: `python3 methods/gse154989_kp_cldn4/analyze.py`

---

## How to read this

- This is **epithelial-only Smart-seq2**, not a whole-tumor 10x atlas.
- Do not quote n=3,891 cells or n=39 `mouseID` as the mouse n.
- Do not treat leftover `Cd3d` as T/NK fraction.
- Do not treat T-AT2–driven MHC-down as a KP finding.
- Thesis is unchanged. Additive mouse result: Cldn4 is present and tracks TJ; T/NK is unscorable; KP IFN/MHC is not down.
