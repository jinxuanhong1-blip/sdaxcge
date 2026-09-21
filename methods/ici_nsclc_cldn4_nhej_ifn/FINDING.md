# Open ICI NSCLC bulk: CLDN4 versus response, with NHEJ and IFN scores

Patient is the unit. Cohorts are not pooled. Positive Cliff's δ means the score is higher in the benefit group. NHEJ is the mean within-cohort z of KEGG hsa03450 (13 genes) and is reported only when Ku70/80, DNA-PKcs, LIG4, XRCC4, and XLF are all present and at least 11/13 genes are present. IFN Teff is the POPLAR 8-gene set (CD8A, GZMA, GZMB, IFNG, EOMES, CXCL9, CXCL10, TBX21), reported at ≥7/8 genes. IFN Hallmark is MSigDB Hallmark interferon-gamma response (200 genes), reported at ≥180 genes.

These are pretreatment or unspecified-time associations. They are not a test of CLDN4 rising after resistance.

## What is not computable

| Source | Status |
|---|---|
| POPLAR / OAK RNA-seq (EGAS00001005013) | Closed. log2(TPM+1), counts, CPM, and clinical tables are Genentech DAC EGAC00001002120. No open per-sample CLDN4, NHEJ, or IFN score joined to atezolizumab response. |
| POPLAR published IFN summary | Fehrenbacher et al., Lancet 2016: in the Teff/IFNγ-high subgroup, atezolizumab versus docetaxel OS HR 0.43 (95% CI 0.24–0.77). That is a treatment-effect subgroup, not a within-arm CLDN4 result, and it is not an NHEJ result. |
| Bessede 2024 supplements | Open figshare 10.1158/1078-0432.c.7077754 is TACSTD2, not CLDN4 or NHEJ, and has no per-sample expression. |
| TCCIA | Zenodo 10.5281/zenodo.7969298 is circRNA (`TCCIA_Ensemble_circRNAs.tsv.gz`, `TCCIA_4_methods_circRNAs.tsv.gz`). No linear CLDN4 and no NHEJ/IFN score. The OAK/POPLAR linear matrices inside TCCIA are the same controlled EGA files. |
| GSE136961 | Oncomine Immune Response panel. CLDN4 absent. NHEJ genes absent. Not scored. |
| GSE285029 | 234 ICI lung libraries, no response/PFS/DCB field in GEO, and the deposited matrix contains negative CLDN4 values. Not scored. |
| GSE182328 | GEO label is Akkermansia, not ORR/DCB/PFS/MPR. Not used as response. |
| GSE248249 | Supplementary file is RAW.tar. Not an open bulk matrix. |

## Response contrasts

δ is Cliff's δ (benefit − no benefit). 95% intervals are percentile bootstraps (2,000 resamples, seed 1). P values are two-sided Mann–Whitney.

| Cohort | n (benefit / other) | Endpoint | CLDN4 δ (p) | NHEJ δ (p) | IFN Teff δ (p) | IFN Hallmark δ (p) |
|---|---:|---|---:|---:|---:|---:|
| GSE218989 | 355 (168 / 187) | GEO responder vs non-responder | −0.051 (0.40) | −0.008 (0.90) | **+0.194 (0.0016)** | **+0.134 (0.029)** |
| GSE190266 | 69 (17 / 52) | no PFS event by 6 months vs event | +0.317 (0.052) | not scored | +0.222 (0.17) | not scored |
| GSE135222 | 27 (7 / 20) | DCB, PFS ≥ 180 days | +0.114 (0.69) | −0.186 (0.50) | +0.457 (0.081) | +0.071 (0.81) |
| GSE207422 | 24 (9 / 15) | pathologic MPR vs NMPR | −0.289 (0.26) | **−0.556 (0.027)** | **+0.615 (0.014)** | +0.244 (0.34) |
| GSE166449 | 22 (7 / 15) | GEO title responder vs non-responder | +0.029 (0.95) | −0.181 (0.53) | +0.467 (0.091) | +0.124 (0.68) |
| GSE126044 | 16 (5 / 11) | GEO ORR responder vs non-responder | −0.527 (0.11) | −0.455 (0.18) | **+1.000 (4.6×10⁻⁴)** | **+0.673 (0.038)** |

CLDN4 does not have one direction. The largest open cohort, GSE218989, is null. GSE126044 (n=16) leans higher in non-responders and matches the earlier δ=−0.53, p=0.115. GSE190266 leans the other way (higher in 6-month benefit, p=0.052, interval crosses 0). This does not support a statement that CLDN4 is higher after ICI resistance, and it does not support a consistent pretreatment CLDN4–nonresponse association.

NHEJ versus benefit is null in GSE218989 (δ=−0.008, p=0.90). The nominal NHEJ result is GSE207422 only: higher NHEJ in NMPR (δ=−0.556, p=0.027; core c-NHEJ δ=−0.496, p=0.049). That cohort is neoadjuvant anti-PD-1 plus chemotherapy, n=9 MPR / 15 NMPR, and the same direction is not present in the large PD-1/PD-L1 cohort.

IFN Teff is higher in the benefit group in all six response cohorts. It is significant in GSE218989, GSE207422, and GSE126044 (complete separation, 5/5 responders above 11/11 non-responders). Hallmark IFN-γ agrees in sign where it is scored, and is significant in GSE218989 and GSE126044. That is the open-cohort counterpart of the POPLAR Teff/IFNγ summary, computed inside the treated arm rather than as an atezolizumab-versus-docetaxel hazard ratio.

## CLDN4 versus the scores

Spearman ρ, same samples as the row above.

| Cohort | n | CLDN4–NHEJ ρ (p) | CLDN4–NHEJ given Teff | CLDN4–NHEJ given MKI67 | CLDN4–Teff ρ (p) |
|---|---:|---:|---:|---:|---:|
| GSE218989 | 355 | **+0.404 (2.1×10⁻¹⁵)** | +0.410 | +0.466 | +0.001 (0.99) |
| GSE135222 | 27 | **+0.440 (0.022)** | +0.522 | +0.512 | +0.191 (0.34) |
| GSE126044 | 16 | +0.232 (0.39) | +0.103 | +0.281 | **−0.603 (0.013)** |
| GSE166449 | 22 | +0.081 (0.72) | +0.119 | +0.119 | −0.193 (0.39) |
| GSE207422 | 24 | +0.190 (0.38) | +0.121 | −0.127 | −0.330 (0.11) |
| GSE190266 | 69 | not scored |  |  | **+0.402 (6.3×10⁻⁴)** |
| GSE253564 pre, no response label | 32 | **−0.389 (0.028)** | −0.396 | −0.376 | **−0.425 (0.015)** |
| GSE248378 post, no response label | 29 | −0.335 (0.075) | −0.241 | +0.118 | **−0.522 (0.0037)** |

In GSE218989, CLDN4 tracks NHEJ and does not track Teff. Residualizing on Teff or on MKI67 does not remove the NHEJ correlation (CLDN4–MKI67 ρ=+0.005). The same positive CLDN4–NHEJ ρ is present in GSE135222. It is not present in the small cohorts, and the sign is negative in the early-stage durvalumab matrices (GSE253564 pretreatment; GSE248378 resection). Those two have no MPR or recurrence column in GEO, so they are correlations only, not response tests, and they are not a pre/post pair.

## Cohort notes

- GSE218989: GEO overall design is bulk RNA-seq of immunotherapy-treated lung cancer; each sample is labeled PD-1/PD-L1 inhibitor and responder or non-responder (168/187). Pre versus post is not in the sample fields. Series title belongs to the parent pan-cancer paper.
- GSE190266: Dijon anti-PD-1 monotherapy, tumor at diagnosis. One early censor excluded (PFS time 0.03, event 0). The deposited TPM stops at gene MTMR14 (16,384 columns), so Ku/DNA-PKcs/XRCC4/XLF and Hallmark genes after that point are absent. NHEJ and Hallmark are not scored. Teff is 7/8 (TBX21 absent).
- GSE135222: GEO `pfs=1` is the progression event. All 27 have PFS time evaluable at 180 days (7 DCB / 20 not). Expression is Ensembl; PRKDC is ENSG00000253729.
- GSE207422: pretreatment bulk biopsies only. BD Rhapsody libraries are excluded. Benefit is pathologic MPR, including pCR.
- GSE126044: all pretreatment. 5/16 are FFPE and were kept. NHEJ is 12/13 (MRE11 absent from the count file). Hallmark is 199/200 (`MARCHF1` absent; `WARS` counted as `WARS1`).
- GSE248378: Teff is 7/8 because IFNG is absent from the FPKM file. NHEJ is 12/13 (DNTT absent).
- GSE253564 and GSE248378 are the same durvalumab ± SBRT trial. Arm is not response.

Figure: `figures/fig_forest_cliffs.png`. GSE218989 scatter: `figures/fig_gse218989_cldn4.png`.

## Reproduce

```bash
python3 methods/ici_nsclc_cldn4_nhej_ifn/analyze.py
```

The script downloads GEO supplements into `data/cache/` (gitignored) if they are missing. Numeric rows: `tables/one_row.tsv`. Per-sample scores: `tables/per_sample.tsv`.
