# Finding — GSE127465 REAL Slingshot/PAGA, CLDN4 only

ADDITIVE. **CLDN4 only.** Zilionis et al., *Immunity* 2019, PMID 30979687; human NSCLC **inDrops** (GSE127465). Tumor epithelium = author Type I / Type II / club / ciliated + `PatientN-specific` malignant. Blood dropped. **No TACSTD2∩CLDN4 dual-high gate.** This folder does not merge GSE148071.

Primary clock: **REAL Slingshot** (Street et al. 2018; R `slingshot` 2.10.0) on Harmony-free PCA + Leiden, start cluster **8** (Leiden 8 max Type II fraction / AT2 (CLDN4-high cluster 7 excluded)). PAGA is geometry only. Inferential unit = **patient**. **n=7 is thin** — Spearman on 7 patients is a sign check, not a precise effect. Barrier/keratin **excludes CLDN4**. Root / start is not CLDN4-high (CLDN4-high Leiden = 7).

**What can be said (n=7 thin).** Malignant+CLDN4 exists (n_malignant=3694, n_CLDN4+=1529). Slingshot start is Leiden 8, not the CLDN4-high cluster 7. CLDN4 vs barrier/keratin (no CLDN4): n=7, ρ=0.643, p=0.119. **What cannot be said.** n=7 patient Spearman is not a precise effect. Do not read Slingshot paths as Type II → each patient's tumor.

## Verdict

REAL Slingshot produced 7 lineage(s) from start Leiden 8 (not CLDN4-high 7). Patient-level CLDN4 vs Slingshot: n=7, ρ=-0.357, p=0.432 — n=7 is thin. CLDN4 vs AT2 score: n=7, ρ=0.500, p=0.253. CLDN4 vs barrier/keratin (CLDN4 excluded): n=7, ρ=0.643, p=0.119. CLDN4 vs malignant-like: n=7, ρ=0.607, p=0.148. PAGA has 1 component(s) at connectivity>0 among 11 Leiden vertices. Author malignant is PatientN-specific; Slingshot MST may stitch discrete tumors. Not a TACSTD2 redo. No both-high gate.

## Honest n

- Human patients deposited: **7** (p1–p7). This catalog n is also the test n. **Thin.**
- Tumor cells (author metadata): **40362**. Blood cells excluded: **14411**.
- Tumor epithelium on the object: **n_cells = 4281**.
- Author malignant (`Patient*`-specific): **3694**. CLDN4>0 among them: **1529** ({'p1': 290, 'p2': 25, 'p3': 309, 'p4': 130, 'p5': 422, 'p6': 97, 'p7': 256}).
- Type II / Type I / club / ciliated: 403 / 46 / 68 / 70.
- Patients with ≥10 epithelial cells: **n = 7**.
- Patients with ≥8 cells in both CLDN4 tertile arms: **n = 7**.
- CLDN4 tertile cells: {'low': 2502, 'high': 1427, 'mid': 352}.
- Author lineage counts: {'Malignant': 3694, 'Type II cells': 403, 'Ciliated cells': 70, 'Club cells': 68, 'Type I cells': 46}.
- Genes absent from locked sets: {'AT2': [], 'AT1': [], 'club': [], 'basal': [], 'ciliated': [], 'barrier_keratin': [], 'malignant_like': []}.
- Dual-high TACSTD2∩CLDN4 gate: **not defined / not used**.
- Slingshot MST will connect patient-specific malignant clusters even when they are transcriptionally discrete. PAGA connectivity is the honesty check.

## Locked choices

- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- Batch: none (no Harmony). Author malignant is PatientN-specific by design.
- Slingshot start: Leiden 8. CLDN4-high cluster 7 was ineligible as root.
- PAGA components at connectivity>0: **1** among 11 Leiden vertices.
- Barrier/keratin: KRT8, KRT18, KRT19, KRT7, CDKN1A, PLAUR (**CLDN4 out**).
- Slingshot version: 2.10.0.

## Lineage table (done criterion)

| lineage | start | end | n_clusters | n_cells | mean_CLDN4 | mean_AT2 | frac_malignant | cluster_path |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| Lineage1 | 8 | 3 | 5 | 2343 | 0.217 | 0.502 | 0.813 | 8->5->4->0->3 |
| Lineage2 | 8 | 9 | 4 | 1109 | 0.348 | 1.095 | 0.527 | 8->5->4->9 |
| Lineage3 | 8 | 2 | 4 | 1484 | 0.345 | 1.009 | 0.706 | 8->5->4->2 |
| Lineage4 | 8 | 6 | 4 | 1318 | 0.365 | 1.222 | 0.675 | 8->5->4->6 |
| Lineage5 | 8 | 1 | 4 | 1491 | 0.483 | 0.795 | 0.713 | 8->5->4->1 |
| Lineage6 | 8 | 7 | 4 | 1292 | 0.451 | 1.050 | 0.669 | 8->5->4->7 |
| Lineage7 | 8 | 10 | 3 | 687 | 0.286 | 1.581 | 0.357 | 8->5->10 |

Machine table: `results/tables/lineages.tsv`.

## Primary (patient-level Spearman, BH inside this list)

**n=7 is thin.** p-values are descriptive.

| Contrast | n_patients | ρ | p | q |
| --- | ---: | ---: | ---: | ---: |
| CLDN4 vs Slingshot | 7 | -0.357 | 0.432 | 0.48 |
| CLDN4 vs DPT (companion) | 7 | -0.679 | 0.0938 | 0.247 |
| CLDN4 vs AT2 score | 7 | 0.500 | 0.253 | 0.316 |
| CLDN4 vs club score | 7 | 0.607 | 0.148 | 0.247 |
| CLDN4 vs basal score | 7 | -0.857 | 0.0137 | 0.0685 |
| CLDN4 vs barrier/keratin (no CLDN4) | 7 | 0.643 | 0.119 | 0.247 |
| CLDN4 vs malignant-like score | 7 | 0.607 | 0.148 | 0.247 |
| CLDN4 vs TACSTD2 (comparator) | 7 | 0.929 | 0.00252 | 0.0252 |
| SFTPC vs Slingshot (control) | 7 | -0.179 | 0.702 | 0.702 |
| AT2 score vs Slingshot (control) | 7 | -0.500 | 0.253 | 0.316 |

## Extra figure — CLDN4-high vs CLDN4-low (patient-paired)

Emitted: **True**. Paired tertile n=7.

| Paired contrast (high − low) | n | Δ median | p |
| --- | ---: | ---: | ---: |
| barrier/keratin (no CLDN4) high vs low | 7 | 0.120 | 0.0156 |
| AT2 high vs low | 7 | 0.020 | 1 |
| malignant-like high vs low | 7 | 0.095 | 0.0156 |
| Slingshot high vs low | 7 | 0.428 | 0.297 |
| DPT high vs low | 7 | 0.021 | 0.688 |

## What this does not claim

- Cell-level p-values are not the claim. n_cells is large by construction.
- **n=7 is thin.** Do not write a precise effect size from the patient Spearman.
- Author `PatientN-specific` is the paper's malignant label, **not CNV re-called here**.
- Slingshot lineages are MST paths on Leiden centers. They are not proof that Type II differentiates into each patient's tumor.
- Do not write “CLDN4 marks the malignant terminal.”
- No TACSTD2∩CLDN4 both-high gate. Not a TACSTD2 redo.
- Not ICI / MPR / RECIST. Not a merge with GSE148071.
- RNA velocity was not run (no spliced/unspliced).

## Outputs

- `results/tables/lineages.tsv` — **done criterion**
- `results/tables/sample_level_spearman.tsv`
- `results/tables/sample_means.tsv`
- `results/figures/fig_trajectory_cldn4.png`
- `results/figures/fig_extra_cldn4_tertile.png`
- `results/figures/fig_extra_lineage_cldn4.png`
- `results/figures/fig_extra_paga.png`
- `results/figures/fig_extra_patient_cldn4.png`
- `results/figures/fig_honest_n.png`
- `results/summary.json`

## Reproduce

```bash
bash methods/gse127465_slingshot_real_cldn4/scripts/install_tools.sh
python3 methods/gse127465_slingshot_real_cldn4/scripts/download.py \
  --out /tmp/gse127465_slingshot
python3 methods/gse127465_slingshot_real_cldn4/scripts/extract_epithelium.py \
  --data /tmp/gse127465_slingshot \
  --out /tmp/gse127465_slingshot/epithelium.h5ad
python3 methods/gse127465_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse127465_slingshot/epithelium.h5ad \
  --outdir methods/gse127465_slingshot_real_cldn4/results \
  --finding methods/gse127465_slingshot_real_cldn4/FINDING.md
```
