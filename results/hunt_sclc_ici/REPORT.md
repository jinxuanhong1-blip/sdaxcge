# Hunt: TACSTD2 (TROP2) / CLDN4 vs the SCLC immune (SCLC-I) axis

**Verdict, one sentence:** IMpower133 cannot answer this from public data; in the two public leftovers that can (George bulk, Chan atlas), **neither ADC target is a marker of the ICI-favorable inflamed phenotype**, and TROP2 is barely expressed on SCLC tumor cells at all.

No EGA. No fabricated IMpower133 expression. No simulated patients.

---

## 1. What was asked, and what is actually public

Question: do **TACSTD2 (TROP2)** and **CLDN4** track **SCLC-I** — the inflamed subtype that Gay et al. 2021 reported as deriving the most benefit from adding atezolizumab to chemo in IMpower133?

| Source | Access | Can compute TACSTD2 × SCLC-I? | What we did |
|---|---|---|---|
| IMpower133 RNA-seq + subtype calls | **EGA-controlled** (`EGAS00001004888`) | **No** | Transcribed published cohort-level leftovers only. Template refuses to run without real files. |
| George et al. 2015 bulk SCLC (n=81) | Public (cBioPortal `sclc_ucologne_2015`) | Yes — this is Gay's discovery cohort | Full analysis |
| Chan et al. 2021 SCLC atlas | Public (cellxgene) | Partial — no SCLC-I tumor-cell label; immune = TME fraction | Full analysis |

IMpower133 leftovers that *are* public are qualitative / cohort-level (see `tables/impower133_public_leftovers.tsv`): subtype mix A 51 / N 23 / I 18 / P 7%; Ayers 18-gene GEP high in SCLC-I; SCLC-I median OS >18 mo on EP+atezo vs ~10 mo on placebo; SCLC-I vs other OS HR 0.566 (0.321–0.998) in the atezo arm only. **None of those leftovers contain TACSTD2 or CLDN4.** That absence is itself a finding about the public record.

---

## 2. Methods (transparent, not a re-implementation of Gay NMF)

- **George subtypes:** argmax over unit-variance {ASCL1, NEUROD1, POU2F3, Ayers 18-gene inflamed GEP}. This is a documented approximation, **not** the original NMF labels (those are unpublished per-sample). Sanity check against Gay's published George mix: A 36/36, N 28/31, P 16/16, I 20/17 (% this analysis / % Gay). Close enough to use; not identical.
- **Signatures:** unweighted mean z-score of published gene lists (Ayers 2017 18-gene GEP, HLA/antigen presentation, checkpoints, STING chemokines, EMT). Rooney CYT = mean log2 of GZMA+PRF1. Where the original score is weighted/proprietary, we say so and do not claim identity.
- **NE score:** mean-z(NE markers) − mean-z(non-NE markers). TACSTD2 and CLDN4 are **excluded** from the non-NE list to avoid circularity.
- **Statistics:** Spearman + BH-FDR; Kruskal–Wallis; Mann–Whitney SCLC-I vs rest. Cells are not treated as independent in Chan — primary inference is **donor-level** (n=20 SCLC donors).
- **Chan labels:** tumor cells are SCLC-A / SCLC-N / SCLC-P. There is **no SCLC-I tumor-cell class**. SCLC-I is a bulk/TME concept.

Reproduce: `bash results/hunt_sclc_ici/data/download_data.sh` then `python3 -m src.hunt_sclc_ici.run_hunt`.

---

## 3. George 2015 bulk (real, n=81)

### TACSTD2 (TROP2)

| Axis | Spearman ρ | BH q |
|---|---:|---:|
| NE score | **−0.49** | **4×10⁻⁵** |
| Epithelial score | **+0.35** | **0.008** |
| Ayers T-cell-inflamed GEP | +0.24 | 0.08 (NS) |
| HLA / antigen presentation | +0.24 | 0.08 (NS) |
| CD8A | +0.12 | 0.39 |
| SCLC-I vs rest (Mann–Whitney) | median 2.91 vs 2.72 log2(FPKM+1) | p=0.40 |

TROP2 tracks **non-neuroendocrine / epithelial** biology. Its link to the inflamed GEP is weak and does **not** survive FDR. It is **not** higher in our SCLC-I calls than in the rest.

### CLDN4

| Axis | Spearman ρ | BH q |
|---|---:|---:|
| Epithelial score | **+0.54** | **2×10⁻⁶** |
| SCLC-A axis (ASCL1) | **+0.44** | **3×10⁻⁴** |
| CD8A | **−0.29** | **0.041** |
| Ayers T-cell-inflamed GEP | −0.12 | 0.39 |
| SCLC-I vs rest | medians 6.11 vs 6.13 | p=0.97 |

CLDN4 is an **epithelial / SCLC-A** gene. The only immune association that survives FDR is **anti-CD8** — the immune-cold direction. It is indistinguishable between SCLC-I and the rest.

### Sanity that the subtype approximation is not nonsense

Independent of TACSTD2/CLDN4, known Gay subtype features recapitulate: SSTR2 highest in SCLC-N (Kruskal p=4×10⁻⁸); DLL3 highest in SCLC-A / lowest in SCLC-P (p=8×10⁻⁶); MICA highest in SCLC-P (p=0.002); EPCAM lowest in SCLC-I (q=0.039). So the I-vs-rest null for the ADC targets is not because subtypes are random.

Figures: `figures/george_corr_heatmap.png`, `george_scatter_ADC_vs_immune.png`, `george_box_TACSTD2_CLDN4_by_subtype.png`.

---

## 4. Chan 2021 atlas (real, 77,143 SCLC cells / 20 donors)

### Compartment: who actually expresses these genes?

In SCLC samples only:

| Compartment | n cells | TACSTD2 detect | CLDN4 detect |
|---|---:|---:|---:|
| Epithelial | 55,876 | **2.2%** | **61%** |
| Lymphoid | 16,461 | 0.13% | 2.2% |
| Myeloid | 3,412 | 1.4% | 3.2% |
| Mesenchymal | 1,394 | 3.6% | 3.7% |

TROP2 is **rare on SCLC tumor cells** and is **not an immune-cell transcript**. Bulk TACSTD2 is therefore not an infiltrate artifact — there is just very little of it. CLDN4 is a bona fide epithelial gene (matches the bulk epithelial correlation).

### Tumor-cell labels (A/N/P — no I)

| Chan label | n cells | TACSTD2 detect | CLDN4 detect |
|---|---:|---:|---:|
| SCLC-A | 31,865 | 0.9% | **75%** |
| SCLC-N | 19,279 | 2.9% | 50% |
| SCLC-P | 3,169 | 3.2% | **5.7%** |

CLDN4 is an SCLC-A (and to a lesser extent SCLC-N) surface gene and is essentially off in SCLC-P. TROP2 stays rare in every tumor-cell class. ASCL1/NEUROD1/POU2F3 detection rates match the labels (sanity: POU2F3 detect 86% in SCLC-P vs <0.5% in A/N).

### Patient-level: tumor ADC-target vs immune fraction (the honest test)

n=20 SCLC donors with ≥30 epithelial cells. Spearman, tumor-cell mean vs fraction immune:

| Gene | ρ vs % immune | p | BH q (all genes tested) |
|---|---:|---:|---:|
| TACSTD2 | +0.52 | 0.020 | **0.17 (NS)** |
| CLDN4 | −0.03 | 0.90 | 0.90 |

The nominal TACSTD2–immune association **does not survive FDR**, and it is **not robust**:

- 2/20 donors have dominant epithelial label **Epithelial Stroma** (not SCLC-A/N/P). They are immune-rich and TROP2-higher. Drop them (n=18 A/N/P-dominant): TACSTD2 vs % immune ρ=0.40, **p=0.10**.
- Treatment-naive only (n=7): ρ=0.57, p=0.18 — underpowered.
- CLDN4 stays null in every subset.

So Chan does **not** support “TROP2 marks the inflamed SCLC-I-like TME” once you stop treating a stroma-driven nominal p as a result.

Figures: `figures/chan_detect_by_compartment.png`, `chan_detect_by_tumor_label.png`, `chan_patient_ADC_vs_immune.png`.

---

## 5. What this does *not* say

- It does **not** say TROP2 ADCs cannot work in SCLC. It says TROP2 RNA is low and not aligned with the ICI-hot subtype. Protein, heterogeneity, and trial activity are out of scope.
- It does **not** say SCLC-I is wrong. Gay's IMpower133 OS trend is a published leftover we cannot re-test without EGA.
- It does **not** assign official Gay NMF labels. Proportions match; per-sample identity may differ.
- Chan n=20 donors is small. A null at this n is evidence against a *large* patient-level association, not proof of zero association.
- We did not digitize Gay Fig 3G heatmap pixels and pretend that is IMpower133 TACSTD2 data.

---

## 6. Bottom line

**IMpower133 public leftovers contain no TACSTD2/CLDN4 measurement.** The two public datasets that can be analyzed agree:

1. **CLDN4** = epithelial / SCLC-A, immune-cold or immune-null. Not SCLC-I.
2. **TACSTD2** = non-NE / epithelial in bulk, **almost off** in SCLC tumor cells in the atlas, and **not** a robust inflamed-TME marker.

If the therapeutic question is “should a TROP2 ADC be aimed at the ICI-hot SCLC-I subset?”, the public data say **no reason to think so**. If the question is “does IMpower133 itself show this?”, the honest answer is **the data are not public**.
