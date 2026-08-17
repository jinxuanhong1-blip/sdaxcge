# Finding — GSE148071 REAL Slingshot/PAGA scored by CLDN4

ADDITIVE. **CLDN4 only.** Public malignant-like epithelium from GSE148071 (Wu et al., *Nat Commun* 2021, PMID 33953163): 42 advanced NSCLC diagnostic biopsies (Singleron GEXSCOPE). Counts = GEO raw UMI. Labels = TISCH2 major-lineage (Malignant + leftover Alveolar / Basal / Epithelial). **Not a TACSTD2 redo.** No dual-high gate.

Primary clock = **real Bioconductor slingshot** (Street et al. 2018; R package loaded). PAGA is geometry. DPT is a comparator only. Root = TISCH Alveolar / AT2, **never CLDN4-high**. Inferential unit = **patient**. Cell-level ρ is descriptive. Barrier/keratin **excludes CLDN4**. IFN = Hallmark IFNα ∩ IFNγ intersection (CLDN4 not in the set).

PR #394 barrier/keratin (CLDN4 excluded) patient Spearman is **given** (n=42, ρ=+0.590, p=3.87e-05) and is **not re-audited** here. This folder asks a different question: CLDN4 + barrier + IFN **along Slingshot/PAGA pseudotime**.

## Verdict

REAL Slingshot ran (17 lineages; principal=Lineage3). Patient-level CLDN4 vs Alveolar-rooted Slingshot PT: n=39, ρ=-0.441, p=0.00499. Barrier/keratin (CLDN4 excluded) vs Slingshot PT: n=39, ρ=-0.380, p=0.0171. IFN vs Slingshot PT: n=39, ρ=-0.405, p=0.0104. CLDN4 vs IFN: n=42, ρ=-0.065, p=0.682. AT2 vs Slingshot PT (root control): n=39, ρ=-0.577, p=0.00012. Paired CLDN4-high vs low barrier/keratin: n=38, W=2.0, Δmed=0.272, p=2.18e-11. Paired CLDN4-high vs low IFN: n=38, W=103.0, Δmed=0.045, p=3.95e-05. PAGA has 1 component(s) at connectivity>0 among 30 Leiden vertices. Root is TISCH Alveolar (median AT2; CLDN4-high excluded); root CLDN4 tertile=mid (not high). Leftover-epithelium–dominant patients (TISCH Malignant <10 on the object): ['P35', 'P39', 'P5']. PR #394 barrier ρ=+0.590 is given and was not re-audited. GEO n=42 is the catalog, not the Spearman n. Mixed advanced NSCLC; no ICI labels. Not a TACSTD2 redo. No both-high gate.

## Honest n

- GEO patients / samples: **42** (Wu 2021). This is the catalog n, not the test n.
- TISCH cells: **82267**. TISCH epithelial (Malignant+Alveolar+Basal+Epithelial): **56265**.
- GEO barcodes matched to TISCH epithelium: **56265**.
- After QC (min 200 genes, min 500 UMI): **n_cells_qc = 56265** in **42** patients.
- Graph cap: max 500/patient when n>20000 (seed 0). Analysis object: **n_cells = 15474** in **n_patients = 42**.
- TISCH lineage on the analysis object: {'Malignant': 11366, 'Alveolar': 2929, 'Basal': 679, 'Epithelial': 500}.
- Slingshot lineages: **17**. Principal = Lineage3 (cells on lineage: {'Lineage1': 6439, 'Lineage2': 6933, 'Lineage3': 7018, 'Lineage4': 5881, 'Lineage5': 6562, 'Lineage6': 4199, 'Lineage7': 3986, 'Lineage8': 4349, 'Lineage9': 4097, 'Lineage10': 4084, 'Lineage11': 3464, 'Lineage12': 4304, 'Lineage13': 2975, 'Lineage14': 3083, 'Lineage15': 2971, 'Lineage16': 2061, 'Lineage17': 2182}).
- Patients with ≥10 analysis cells: **n = 42** (catalog n equals this count; still not 42 malignant tumors).
- Patients with finite principal Slingshot PT (Lineage3) used for along-PT Spearman: **n = 39**. Dropped for missing principal-lineage membership: **P11, P29, P42**.
- Patients with ≥10 TISCH Malignant cells on the object: **n = 39**. Leftover-epithelium–dominant (malignant <10): **P5, P35, P39**.
- Patients with ≥8 cells in both CLDN4-high and CLDN4-low arms (paired extra): **n = 38**.
- CLDN4 tertile cells: low 5302, mid 5014, high 5158.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': [], 'IFN': []}.
- Slingshot: available=True; version=2.10.0.
- No ICI / RECIST / MPR labels. Histology is not on the GEO series matrix — no LUAD/LUSC split.
- Do not write “n=42 malignant tumors.”

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30; Slingshot PCs 10.
- Slingshot/DPT root: TISCH Alveolar (median AT2; CLDN4-high excluded) (root cell index 10534, patient P35, lineage Alveolar, Leiden 27, CLDN4 tertile mid).
- PAGA components at connectivity>0: **1** among 30 Leiden vertices.
- Barrier/keratin genes: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- IFN genes: Hallmark IFNα ∩ IFNγ (73; WARS/WARS1 alias accepted).
- PR #394 barrier ρ=+0.590 cited, not re-fit.

## Primary (patient-level Spearman along pseudotime, BH inside this list)

| Contrast | n_patients | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Slingshot PT | 39 | -0.441 | 0.00499 | 0.00997 |
| barrier/keratin (no CLDN4) vs Slingshot PT | 39 | -0.380 | 0.0171 | 0.0228 |
| IFN vs Slingshot PT | 39 | -0.405 | 0.0104 | 0.0167 |
| CLDN4 vs IFN | 42 | -0.065 | 0.682 | 0.682 |
| AT2 vs Slingshot PT (control) | 39 | -0.577 | 0.00012 | 0.00096 |
| SFTPC vs Slingshot PT (control) | 39 | -0.451 | 0.00399 | 0.00997 |
| CLDN4 vs DPT (comparator clock) | 42 | -0.256 | 0.102 | 0.117 |
| Slingshot PT vs DPT (concordance) | 39 | 0.548 | 0.000301 | 0.0012 |

## Sensitivity (not in the BH family)

| Contrast | n_patients | ρ | p |
| --- | ---: | ---: | ---: |
| malignant-only CLDN4 vs Slingshot PT | 36 | -0.608 | 8.46e-05 |
| malignant-only barrier vs Slingshot PT | 36 | -0.528 | 0.000942 |
| malignant-only IFN vs Slingshot PT | 36 | -0.429 | 0.00897 |
| malignant-only CLDN4 vs IFN | 39 | -0.064 | 0.697 |
| drop leftover-dominant CLDN4 vs Slingshot PT | 36 | -0.608 | 8.46e-05 |
| CLDN4 vs Slingshot PT (mean across lineages) | 42 | -0.460 | 0.00215 |
| IFN vs Slingshot PT (mean across lineages) | 42 | -0.146 | 0.356 |
| Lineage1 patient CLDN4 vs PT | 31 | -0.307 | 0.0927 |
| Lineage1 patient IFN vs PT | 31 | -0.328 | 0.0718 |
| Lineage2 patient CLDN4 vs PT | 30 | -0.553 | 0.00153 |
| Lineage2 patient IFN vs PT | 30 | -0.387 | 0.0347 |
| Lineage3 patient CLDN4 vs PT | 30 | -0.570 | 0.001 |
| Lineage3 patient IFN vs PT | 30 | -0.358 | 0.0518 |
| Lineage4 patient CLDN4 vs PT | 28 | -0.593 | 0.000887 |
| Lineage4 patient IFN vs PT | 28 | -0.196 | 0.316 |
| Lineage5 patient CLDN4 vs PT | 30 | -0.535 | 0.00229 |
| Lineage5 patient IFN vs PT | 30 | -0.369 | 0.045 |
| Lineage6 patient CLDN4 vs PT | 24 | -0.010 | 0.965 |
| Lineage6 patient IFN vs PT | 24 | 0.323 | 0.124 |
| Lineage7 patient CLDN4 vs PT | 24 | -0.286 | 0.175 |
| Lineage7 patient IFN vs PT | 24 | 0.268 | 0.206 |
| Lineage8 patient CLDN4 vs PT | 24 | -0.274 | 0.195 |
| Lineage8 patient IFN vs PT | 24 | 0.164 | 0.443 |
| Lineage9 patient CLDN4 vs PT | 21 | -0.513 | 0.0174 |
| Lineage9 patient IFN vs PT | 21 | 0.184 | 0.424 |
| Lineage10 patient CLDN4 vs PT | 21 | -0.481 | 0.0275 |
| Lineage10 patient IFN vs PT | 21 | 0.119 | 0.606 |
| Lineage11 patient CLDN4 vs PT | 19 | -0.414 | 0.078 |
| Lineage11 patient IFN vs PT | 19 | 0.265 | 0.273 |
| Lineage12 patient CLDN4 vs PT | 19 | -0.335 | 0.161 |
| Lineage12 patient IFN vs PT | 19 | 0.179 | 0.464 |
| Lineage13 patient CLDN4 vs PT | 14 | 0.077 | 0.794 |
| Lineage13 patient IFN vs PT | 14 | -0.363 | 0.203 |
| Lineage14 patient CLDN4 vs PT | 15 | 0.107 | 0.704 |
| Lineage14 patient IFN vs PT | 15 | -0.225 | 0.42 |
| Lineage15 patient CLDN4 vs PT | 16 | -0.003 | 0.991 |
| Lineage15 patient IFN vs PT | 16 | -0.188 | 0.485 |
| Lineage16 patient CLDN4 vs PT | 14 | -0.213 | 0.464 |
| Lineage16 patient IFN vs PT | 14 | -0.011 | 0.97 |
| Lineage17 patient CLDN4 vs PT | 13 | 0.462 | 0.112 |
| Lineage17 patient IFN vs PT | 13 | -0.379 | 0.201 |

## Slingshot lineage table (done criterion)

| Lineage | n_cells | n_patients_elig | end | CLDN4~PT ρ | barrier~PT ρ | IFN~PT ρ |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| Lineage1 | 6439 | 31 | 26 | -0.307 | -0.154 | -0.328 |
| Lineage2 | 6933 | 30 | 0 | -0.553 | -0.406 | -0.387 |
| Lineage3 | 7018 | 30 | 3 | -0.570 | -0.421 | -0.358 |
| Lineage4 | 5881 | 28 | 14 | -0.593 | -0.434 | -0.196 |
| Lineage5 | 6562 | 30 | 21 | -0.535 | -0.402 | -0.369 |
| Lineage6 | 4199 | 24 | 6 | -0.010 | -0.158 | 0.323 |
| Lineage7 | 3986 | 24 | 20 | -0.286 | 0.052 | 0.268 |
| Lineage8 | 4349 | 24 | 22 | -0.274 | -0.057 | 0.164 |
| Lineage9 | 4097 | 21 | 24 | -0.513 | -0.319 | 0.184 |
| Lineage10 | 4084 | 21 | 25 | -0.481 | -0.258 | 0.119 |
| Lineage11 | 3464 | 19 | 8 | -0.414 | -0.223 | 0.265 |
| Lineage12 | 4304 | 19 | 9 | -0.335 | -0.419 | 0.179 |
| Lineage13 | 2975 | 14 | 16 | 0.077 | 0.371 | -0.363 |
| Lineage14 | 3083 | 15 | 19 | 0.107 | 0.518 | -0.225 |
| Lineage15 | 2971 | 16 | 29 | -0.003 | 0.141 | -0.188 |
| Lineage16 | 2061 | 14 | 12 | -0.213 | -0.332 | -0.011 |
| Lineage17 | 2182 | 13 | 18 | 0.462 | 0.473 | -0.379 |

## Extra figure — CLDN4-high vs CLDN4-low (patient-paired) and along-PT

Emitted: **True**. Rule: always emit (requested extra figures). Paired tertile n=38. PR #394 barrier patient ρ=+0.590 is cited, not re-audited.

| Contrast | n | W | Δmed (high−low) | p |
| --- | ---: | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 38 | 2.0 | 0.272 | 2.18e-11 |
| IFN high vs low | 38 | 103.0 | 0.045 | 3.95e-05 |
| Slingshot PT high vs low | 32 | 144.0 | -0.007 | 0.024 |
| AT2 score high vs low | 38 | 332.0 | 0.003 | 0.769 |

## Caveats

- Cell-level p-values are not the claim. n_cells is large by construction.
- TISCH Malignant is a processed label, **not CNV re-called here**.
- Several patients are leftover-epithelium–dominant (almost no TISCH Malignant). They stay in the mixed graph and are dropped from the malignant-only sensitivity.
- Patient batch is strong (42 tumors). Per-patient cap reduces one-sample domination; it is not Harmony.
- Mixed advanced NSCLC. Do not write LUAD-only. Do not write ICI language.
- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.
- Do not write “AT2 differentiates into NSCLC because PAGA is connected.”
- Do not write “CLDN4 marks the malignant terminal.”
- Slingshot is an ordering, not a clock. Direction is the external Alveolar/AT2 root.
- PR #394 barrier ρ is given. This analysis does not re-fit that contrast as a discovery.

## Outputs (done when these exist)

- `tables/patient_level.tsv`
- `tables/lineage_level.tsv`
- `tables/tisch_lineage_level.tsv`
- `tables/patient_level_spearman.tsv`
- `figures/fig_trajectory_cldn4.png`
- `figures/fig_extra_along_pt.png`
- `figures/fig_extra_cldn4_tertile.png`
- `figures/fig_honest_n.png`

## Reproduce

```bash
bash methods/gse148071_slingshot_real_cldn4/scripts/install_tools.sh
bash methods/gse148071_slingshot_real_cldn4/scripts/download.sh /tmp/gse148071_sling_data
python3 methods/gse148071_slingshot_real_cldn4/scripts/extract_epithelium.py \
  --data /tmp/gse148071_sling_data \
  --out /tmp/gse148071_sling_data/epithelium.h5ad
python3 methods/gse148071_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse148071_sling_data/epithelium.h5ad \
  --outdir methods/gse148071_slingshot_real_cldn4 \
  --finding methods/gse148071_slingshot_real_cldn4/FINDING.md
```

Trajectory: `figures/fig_trajectory_cldn4.png`. Extra along-PT: `figures/fig_extra_along_pt.png`.
