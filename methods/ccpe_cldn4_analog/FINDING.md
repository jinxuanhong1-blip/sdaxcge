# FINDING — C-CPE pharmacologic CLDN4-axis analog (not SKB264, not lung)

**Verdict: C-CPE and CLDN4 siRNA do not move IFN / MHC-I / APM up in a defensible way.**

Public **GSE22421** (SKOV-3, C-CPE 5 µg/ml 72 h vs untreated; two-color array; **n=3**) is the pharmacologic analog. **GSE22493** (SKOV-3-IP-Luc, CLDN4 lentiviral siRNA vs CLDN4 overexpression; same platform; **n=3**) was paired because gene symbols map. **TACSTD2 is not on GPL10555** — no number is invented.

| Question | GSE22421 C-CPE / untreated (n=3) | GSE22493 CLDN4 siRNA / OE (n=3) |
|---|---|---|
| CLDN4 mRNA | deposited median **0.00** (0.11 / −0.42 / 0.00); ScanArray median **−0.56** | deposited median **−1.23** (n=2; GSM558700 missing); knockdown **not clean** |
| TACSTD2 | **absent** on GPL10555 | **absent** on GPL10555 |
| Priority 6 (IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A) | 4 up / 2 down vs 0; **none** p < 0.05 | 3 up / 3 down; **none** p < 0.05 |
| MHC-I / APM (12/16 measured) | 9 up / 3 down vs 0, but set median **+0.19 vs background +0.34**, MW p=**0.22** | 3 up / 9 down; set median **−0.45 vs bg −0.03**, MW p=**0.16** |
| IFN_IMMUNE (60/73) | set median **+0.39 vs bg +0.34**, MW p=**0.73** | set median **−0.23 vs bg −0.03**, MW p=**0.30** |
| Hallmark IFN-α / IFN-γ GSEA | NES **0.77** p=**0.86** / NES **1.01** p=**0.47** | NES **−0.68** p=**0.59** / NES **−1.02** p=**0.26** |
| Classical MHC-I | HLA-A +0.03; HLA-C **−0.40**; B2M **−0.92** | HLA-A **−1.60**; HLA-B **−1.56**; HLA-C **−1.25**; B2M **−2.56** |

GSE22421’s whole-array background median is **+0.35**. Counting genes “up vs 0” overstates IFN/TJ. Versus background, APM and the priority six sit **below** the array, and IFN_IMMUNE is at background. Hallmark IFN sets are null. The strongest Hallmark hit is **MYC_TARGETS_V1** (NES 1.48, perm p=0.000) — consistent with the paper’s metabolism/ubiquitin story, not antigen presentation.

This is **ovarian SKOV-3**, not lung, and **C-CPE**, not SKB264.

---

## What was asked

Additive public analog of a CLDN4-axis pharmacologic perturbation. Report CLDN4, TACSTD2, IFN / MHC-I / APM, and tight-junction log2FC and/or GSEA. Honest n. Do not invent numbers.

## Data (honest n)

| Item | GSE22421 | GSE22493 |
|---|---|---|
| GEO | [GSE22421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22421) | [GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493) |
| Paper | Gao et al. *Clin Cancer Res* 2011;17:1065–74. PMID [21123456](https://pubmed.ncbi.nlm.nih.gov/21123456/) | none deposited |
| Platform | [GPL10555](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10555) BWH Human Release 3.0 (Operon 60-mer, ~36k probes) | same |
| Model | SKOV-3 ovarian line | SKOV-3-IP-Luc ovarian line |
| Arrays | GSM557400 / 401 / 402 | GSM558700 / 701 / 702 |
| n | **3** two-color slides | **3** two-color slides |
| Cy3 / Cy5 | untreated / C-CPE 5 µg/ml 72 h | CLDN4 overexpression / CLDN4 lentiviral siRNA |
| Deposited VALUE | log2(treated/control) | log2(knockdown/control) |
| Control | untreated | **overexpression, not scramble/WT** |

Both records mention dye-swap. All six deposited GSMs have the same dye assignment (Cy3 = control, Cy5 = treated). Submitter is the same lab (Zhijian Gao, BWH).

**GSE22493 was paired, not skipped.** The two series share GPL10555. Symbols map by the same ORF / DESCRIPTION-prefix rule (25,102 probes; 16,834 symbols). ISG15 is present only as historical **G1P2**. TACSTD2 / TROP2 / GA733-1 / M1S1 are not annotated; only TACSTD1 is on the platform.

## Methods

1. Primary metric: GEO series-matrix deposited `VALUE`. Per gene, median across probes on each array, then median/mean across arrays.
2. One-sample t vs 0 when n≥2 (**descriptive**; n=3 has almost no power). BH-FDR within each pre-specified panel.
3. Gene-set location: two-sided Mann–Whitney of gene-mean log2 vs all genes with ≥2 arrays (background n=16,832 on GSE22421; 15,853 on GSE22493).
4. Optional GSEA: preranked weighted KS on gene-mean log2, Hallmark 2023.2 GMT, 1,000 gene-set permutations. Descriptive.
5. Sensitivity: ScanArray `log2((Ch2 Median−B)/(Ch1 Median−B))` when both channels are >0. Not the author ratio.
6. Lists frozen in `scripts/gene_sets.py` before looking at fold-changes. Missing genes stay NA.

No limma, no post-hoc gene picking, no p-values invented for absent genes.

## CLDN4 / TACSTD2 / CLDN3

`tables/axis_genes.tsv`, `tables/cldn4_diagnostic.tsv`

| Gene | GSE22421 deposited (3 arrays) | median | ScanArray median | GSE22493 deposited | median |
|---|---|---:|---:|---|---:|
| CLDN4 | +0.11 / −0.42 / 0.00 | **0.00** | **−0.56** | missing / −1.74 / −0.71 | **−1.23** (n=2) |
| CLDN3 | −0.25 / −0.67 / −2.64 | **−0.67** | **−1.81** (n=2) | −0.79 / −0.86 / −1.03 | **−0.86** (p=0.0064) |
| TACSTD2 | not on GPL10555 | — | — | not on GPL10555 | — |

C-CPE is a CLDN3/4 binder. The paper reports **protein** loss and TJ relocalization, not a transcriptional CLDN4 knockdown. On this array, CLDN4 mRNA is unchanged (deposited median 0). GSM557400 ScanArray CLDN4 is low-intensity (Cy3 Median−B = 44, Cy5 = 9) and should not be over-read. GSE22493 CLDN4 is negative on the two available deposited arrays, but GSM558701 raw Cy5/Cy3 is **+1.22** (sign flip) and the control is overexpression. **Neither series confirms a clean CLDN4 mRNA loss as the sole contrast.**

## Priority IFN / MHC-I six

`tables/priority_genes.tsv`

| Gene | GSE22421 median | direction | p (t vs 0) | GSE22493 median | direction | p |
|---|---:|---|---:|---:|---|---:|
| IFI27 | +0.54 | UP | 0.20 | −2.32 | DOWN | 0.73 |
| OAS2 | +0.08 | UP | 0.45 | +0.11 | UP | 0.46 |
| IFIT1 | −0.20 | DOWN | 0.78 | −2.12 | DOWN | 0.055 |
| MX1 | +0.18 | UP | 0.17 | +0.36 (n=2) | UP | 0.56 |
| ISG15 (G1P2) | −0.40 | DOWN | 0.17 | +0.71 | UP | 0.62 |
| HLA-A (11 probes) | +0.03 | UP | 0.59 | −1.60 | DOWN | 0.41 |

No priority gene reaches nominal p < 0.05 on C-CPE. Magnitudes are small except IFI27, and IFI27 ScanArray flips to down. HLA-A’s 11 probes do not agree (`tables/probe_level.tsv`). This is not an MHC-I-up call.

## MHC-I / APM

`tables/apm_genes.tsv`

GSE22421, 12/16 measured (NLRC5, ERAP1, ERAP2, PDIA3 absent):

| Gene | median log2 | p | direction |
|---|---:|---:|---|
| HLA-A | +0.03 | 0.59 | UP |
| HLA-B | +0.54 | 0.50 | UP |
| HLA-C | −0.40 | 0.073 | DOWN |
| B2M | −0.92 | 0.094 | DOWN |
| TAP1 | +0.57 | 0.35 | UP |
| TAP2 | −0.47 | 0.44 | DOWN |
| TAPBP | +0.26 | 0.83 | UP |
| PSMB8 | +0.23 | 0.20 | UP |
| PSMB9 | +0.28 | 0.76 | UP |
| PSMB10 | +0.93 | 0.28 | UP |
| CALR | +0.24 | 0.56 | UP |
| CANX | +0.29 | 0.13 | UP |

Nine genes are positive vs 0. That count is **not** an APM-up signature: the APM set median of gene-means is **+0.19 against a +0.34 background** (MW p=0.22, shift DOWN vs background). Classical MHC-I is split (HLA-C and B2M down). No APM gene has within-panel q < 0.05.

GSE22493 APM is worse: 9/12 down, including HLA-A/B/C, B2M, TAP1, TAP2.

## IFN set and Hallmark GSEA

`tables/geneset_stats.tsv`, `tables/hallmark_gsea_preranked.tsv`, `tables/ifn_genes.tsv`

| Set | GSE22421 measured | set median | bg median | MW p | GSEA NES (perm p) |
|---|---:|---:|---:|---:|---|
| PRIORITY | 6/6 | +0.15 | +0.34 | 0.35 | — |
| APM | 12/16 | +0.19 | +0.34 | 0.22 | — |
| IFN_IMMUNE | 60/73 | +0.39 | +0.34 | 0.73 | — |
| TJ | 25/27 | +0.19 | +0.34 | 0.61 | — |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | 70 | — | — | — | 0.77 (0.86) |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 154 | — | — | — | 1.01 (0.47) |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | 164 | — | — | — | 1.25 (0.044) |
| HALLMARK_MYC_TARGETS_V1 | 154 | — | — | — | 1.48 (0.000) |

A few inflammatory transcripts are high on C-CPE (CCL5 median +1.72, p=0.0087; IL6 +1.44, p=0.020; IL15 +0.82, p=0.0092; OAS3 +0.42, p=0.0035). Panel FDR does not pass (lowest IFN q=0.18). These are not classical MHC-I / APM opening, and they do not lift the IFN Hallmark walks.

GSE22493 IFN_IMMUNE is slightly **below** background (median −0.23 vs −0.03, p=0.30). Hallmark IFN-α/γ NES are negative.

## Tight junction

`tables/tj_genes.tsv`

C-CPE does not produce a TJ-down mRNA signature. Versus 0, 20/25 measured TJ genes are up; versus background the set is still down (median +0.19 vs +0.34, p=0.61). CLDN1 is the most consistent up gene (median +1.83, p=0.0033, q=0.083). CLDN3 is down. CLDN4 is zero. CDH1 median −0.44. GSE22493 TJ is mixed-to-down (7 up / 16 down among genes with ≥2 arrays); CLDN3 is the only TJ gene with a small p (0.0064, q=0.15).

The paper’s C-CPE effect on CLDN4 is **protein localization / barrier**, not a transcriptional TJ collapse on this chip.

## Pairing GSE22421 vs GSE22493

`tables/paired_gse22421_vs_gse22493.tsv`

Same platform, same mapping, same lists. The contrasts are **not** the same biology (ligand vs overexpression-vs-siRNA). Concordant deposited medians among priority+APM genes that are measured on both: OAS2 up/up (both ~0.1), IFIT1 down/down, MX1 up/up (tiny), PSMB8/PSMB9/CALR up/up, HLA-C/B2M/TAP2 down/down. Discordant: IFI27, ISG15, HLA-A, HLA-B, TAP1. There is **no shared IFN/MHC-I-up signature**.

## Honest ceiling

1. **C-CPE (GSE22421, n=3) does not open IFN / MHC-I / APM** on deposited log2 ratios or Hallmark GSEA. Background-adjusted APM is down; IFN is null; B2M and HLA-C trend down.
2. **CLDN4 siRNA (GSE22493, n=3) does not open IFN / MHC-I / APM.** Control is overexpression. CLDN4 knockdown is not clean on the array.
3. **TACSTD2 cannot be reported** — not annotated on GPL10555.
4. Do not export this to lung, ICI, or SKB264. C-CPE is a non-toxic CLDN3/4-binding fragment used here as a chemosensitizer analog, not a TROP2 ADC.

## What this does not claim

- It does not claim C-CPE failed to downregulate CLDN4 **protein** (the paper’s Western / IF). That assay is not in the GEO files.
- It does not claim a clean negative on every ISG. n=3 two-color Operon arrays are noisy; GSM-level sign flips are common.
- It does not treat dye-swap text as extra samples.
- It does not use SKB264, sac-TMT, or any lung ICI cohort.

## Outputs

| File | Contents |
|---|---|
| `tables/key_genes.tsv` | CLDN4, CLDN3, TACSTD2, priority 6, B2M, TAP1/2, HLA-B/C |
| `tables/axis_genes.tsv` | CLDN4 / CLDN3 / TACSTD2 |
| `tables/priority_genes.tsv` | IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A |
| `tables/apm_genes.tsv` | 16-gene MHC-I / APM |
| `tables/ifn_genes.tsv` | 73-gene IFN_IMMUNE |
| `tables/tj_genes.tsv` | 27-gene TJ panel |
| `tables/geneset_stats.tsv` | Mann–Whitney vs background |
| `tables/hallmark_gsea_preranked.tsv` | Hallmark NES / perm p |
| `tables/paired_gse22421_vs_gse22493.tsv` | same-gene deposited medians |
| `tables/cldn4_diagnostic.tsv` | deposited vs ScanArray CLDN4 |
| `tables/mapping_audit.tsv` | why 22493 was paired; why TACSTD2 is NA |
| `tables/sample_table.tsv` | channel assignment |
| `figures/fig1_gse22421_axis_priority.png` | CLDN4/CLDN3 + priority 6 |
| `figures/fig2_gse22421_apm_tj.png` | APM and TJ medians |
| `figures/fig3_paired_gse22421_gse22493.png` | paired scatter |
| `figures/fig4_geneset_medians.png` | set location |
| `figures/fig5_gse22421_ifn_panel.png` | IFN_IMMUNE per gene |
| `key_stats.json` | machine-readable summary |

## Reproduce

```bash
pip install -r methods/ccpe_cldn4_analog/scripts/requirements.txt
python3 methods/ccpe_cldn4_analog/scripts/download_data.py   # → $CCPE_CLDN4_DATA (default /tmp/ccpe_cldn4_analog)
python3 methods/ccpe_cldn4_analog/scripts/analyze.py
```

Public file md5 (this run):

| File | md5 |
|---|---|
| GSE22421_series_matrix.txt.gz | `71d4c56fd22aae0ee1fcf671aa2b9d4b` |
| GSE22421_family.soft.gz | `8678232b060a3fd4d13704f27626d409` |
| GSE22421_RAW.tar | `50fe9b2670b3d1714fa7c7c975ba6a7b` |
| GSE22493_series_matrix.txt.gz | `8f8a07cb8f3f3a90396c3feb1eb3ea0f` |
| GSE22493_family.soft.gz | `5d88cef4c873370970e6ba6be82d7b48` |
| GSE22493_RAW.tar | `9bae4611ed5e75f28b387b32f8f793db` |
