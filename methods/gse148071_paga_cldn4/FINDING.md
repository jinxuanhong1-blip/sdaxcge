# Finding — GSE148071 PAGA / diffusion scored by CLDN4

Additive public slice. **Not a TACSTD2 redo.** Primary readout is **CLDN4** on TISCH-labeled epithelium from GSE148071 (Wu et al., *Nat Commun* 2021, PMID 33953163): 42 advanced NSCLC diagnostic biopsies (Singleron GEXSCOPE). Counts are the public GEO raw UMI matrices. Labels are TISCH2 major-lineage (Malignant + leftover Alveolar / Basal / Epithelial). GEO does not ship author cell types or histology on the series matrix.

PAGA + Alveolar/AT2-rooted diffusion pseudotime. Inferential unit = **patient** (one sample each). Cell-level ρ is descriptive. Barrier/keratin score **excludes CLDN4** (no circularity). Root is TISCH Alveolar, never CLDN4-high. Eligible Spearman n happened to equal the GEO catalog (**n=42**) because every biopsy retained ≥10 analysis cells after the 500/patient cap. That is not 42 malignant tumors: leftover-epithelium–dominant **P5 / P35 / P39** stay in the mixed n=42 and are dropped from the malignant-only **n=39**.

## Verdict

Sample-level CLDN4 vs Alveolar-rooted DPT: n=42, ρ=-0.335, p=0.03. CLDN4 vs AT2 score: n=42, ρ=0.103, p=0.518. CLDN4 vs barrier/keratin (CLDN4 excluded from the score): n=42, ρ=0.590, p=3.87e-05. CLDN4 vs malignant-like: n=42, ρ=0.550, p=0.000161. CLDN4 vs TACSTD2 (comparator only): n=42, ρ=0.492, p=0.000933. Malignant-only CLDN4 vs DPT: n=39, ρ=-0.344, p=0.032 (p<0.05). Paired CLDN4-high vs low barrier/keratin: n=38, W=2.0, Δmed=0.272, p=2.18e-11. Paired CLDN4-high vs low AT2: n=38, W=332.0, Δmed=0.003, p=0.769. PAGA has 1 component(s) at connectivity>0 among 30 Leiden vertices. Leftover-epithelium–dominant patients (TISCH Malignant <10 on the object): ['P35', 'P39', 'P5']. GEO n=42 is the catalog, not the Spearman n. Mixed advanced NSCLC; no ICI labels. Not a TACSTD2 redo. No both-high gate.

## Honest n

- GEO patients / samples: **42** (Wu 2021). This is the catalog n, not the test n.
- TISCH cells: **82267**. TISCH epithelial (Malignant+Alveolar+Basal+Epithelial): **56265**.
- GEO barcodes matched to TISCH epithelium: **56265**.
- After QC (min 200 genes, min 500 UMI): **n_cells_qc = 56265** in **42** patients.
- Graph cap: max 500/patient when n>20000 (seed 0). Analysis object: **n_cells = 15474** in **n_patients = 42**.
- Lineage on the analysis object: {'Malignant': 11366, 'Alveolar': 2929, 'Basal': 679, 'Epithelial': 500}.
- Patients with ≥10 analysis cells (Spearman): **n = 42** (catalog n equals eligible n here; still not 42 malignant tumors).
- Patients with ≥10 TISCH Malignant cells on the object: **n = 39**. Leftover-epithelium–dominant (malignant <10): **P5, P35, P39**.
- Patients with ≥8 cells in both CLDN4-high and CLDN4-low arms (paired extra): **n = 38**.
- CLDN4 tertile cells: low 5302, mid 5014, high 5158.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': []}.
- No ICI / RECIST / MPR labels. Histology is not on the GEO series matrix — no LUAD/LUSC split.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- DPT root: TISCH Alveolar (median AT2 score) (root cell index 9758, patient P33).
- PAGA components at connectivity>0: **1** among 30 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).

## Primary (sample-level Spearman, BH inside this list)

| Contrast | n_samples | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs DPT | 42 | -0.335 | 0.03 | 0.0675 |
| CLDN4 vs AT2 score | 42 | 0.103 | 0.518 | 0.582 |
| CLDN4 vs club score | 42 | 0.273 | 0.0806 | 0.104 |
| CLDN4 vs basal score | 42 | 0.051 | 0.75 | 0.75 |
| CLDN4 vs barrier/keratin score | 42 | 0.590 | 3.87e-05 | 0.000349 |
| CLDN4 vs malignant-like score | 42 | 0.550 | 0.000161 | 0.000727 |
| CLDN4 vs TACSTD2 (comparator) | 42 | 0.492 | 0.000933 | 0.0028 |
| TACSTD2 vs DPT (comparator) | 42 | -0.321 | 0.0381 | 0.0686 |
| SFTPC vs DPT (control) | 42 | -0.279 | 0.0738 | 0.104 |

## Sensitivity (not in the BH family)

| Contrast | n_samples | ρ | p |
| --- | ---: | ---: | ---: |
| malignant-only CLDN4 vs DPT | 39 | -0.344 | 0.032 |
| malignant-only CLDN4 vs AT2 score | 39 | 0.170 | 0.302 |
| malignant-only CLDN4 vs barrier/keratin | 39 | 0.577 | 0.000121 |
| drop leftover-dominant CLDN4 vs DPT | 39 | -0.344 | 0.032 |

## Extra figure — CLDN4-high vs CLDN4-low (sample-paired)

Emitted: **True**. Rule: always emit (requested extra figure). Sample-level Spearman(CLDN4, barrier_keratin_no_CLDN4) n=42, ρ=0.59, p=3.87e-05. Paired tertile n=38.

| Contrast | n | W | Δmed (high−low) | p |
| --- | ---: | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 38 | 2.0 | 0.272 | 2.18e-11 |
| AT2 score high vs low | 38 | 332.0 | 0.003 | 0.769 |
| DPT high vs low | 38 | 359.0 | 0.000 | 0.875 |

## Caveats

- Cell-level p-values are not the claim. n_cells is large by construction.
- TISCH Malignant is a processed label, **not CNV re-called here**.
- Several patients are leftover-epithelium–dominant (almost no TISCH Malignant). They stay in the mixed graph and are dropped from the malignant-only sensitivity.
- Patient batch is strong (42 tumors). Per-patient cap reduces one-sample domination; it is not Harmony.
- Mixed advanced NSCLC. Do not write LUAD-only. Do not write ICI language.
- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.
- Do not write “AT2 differentiates into NSCLC because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”

## Reproduce

```bash
pip install -r methods/gse148071_paga_cldn4/requirements.txt
bash methods/gse148071_paga_cldn4/scripts/download.sh /tmp/gse148071_paga_data
python3 methods/gse148071_paga_cldn4/scripts/extract_epithelium.py \
  --data /tmp/gse148071_paga_data \
  --out /tmp/gse148071_paga_data/epithelium.h5ad
python3 methods/gse148071_paga_cldn4/scripts/analyze_paga_cldn4.py \
  --input /tmp/gse148071_paga_data/epithelium.h5ad \
  --outdir methods/gse148071_paga_cldn4 \
  --finding methods/gse148071_paga_cldn4/FINDING.md
```

Trajectory figure: `figures/fig_trajectory_cldn4.png`. Extra figure: `figures/fig_extra_cldn4_tertile.png`.
