# Methods — GSE207422 Palantir destinies, CLDN4-only, MPR labels

Additive public slice. **Not dual-high.** CLDN4 is the readout. TACSTD2 is a companion gene and is never a gate. Barrier / TJ scores **exclude CLDN4**.

This is **real Palantir** (Setty et al., 2019): diffusion maps, multiscale space, Palantir pseudotime, **branch / fate probabilities (destinies)**, and entropy. PAGA is the cluster graph. Diffusion pseudotime (DPT) is not the primary fate model. A parallel MPR-traj / DPT-only agent does not replace this run.

## Object

- Accession: GSE207422 (Hu et al., *Genome Medicine* 2023; neoadjuvant anti–PD-1 NSCLC).
- Public GEO UMI matrix only. Author CopyKAT / epithelium RDS barcodes are **not** on GEO.
- Lineage: argmax of compact marker scores (epithelial vs T/NK/B/myeloid/…).
- **A3-malignant-like** = epithelial AND zero UMI for `SFTPA2`, `AGER`, `SCGB1A1`, `SCGB3A1`, `TPPP3`.
- This is not Hu et al. CopyKAT. Residual unmarked epithelium can leak.
- MPR labels: paper groups locked from the GEO sample sheet / Hu et al. Table S1 (`pCR` collapsed to MPR). Three treatment-naive (TN) samples have no MPR label.

## Palantir / PAGA

- Recompute HVG / PCA / neighbors **on the malignant-like subset**.
- Leiden resolution 0.6; HVG 3000; neighbors 30; PCs 30.
- PAGA on Leiden (display + connectivity table).
- Palantir `run_diffusion_maps` → `determine_multiscale_space` → `run_palantir`.
- **Root is not CLDN4-high.** Early cell is chosen among the CLDN4-low tertile as the cell farthest from the CLDN4-high centroid in PCA. `use_early_cell_as_start=True`.
- **Terminals are not defined as CLDN4-high or TACSTD2-high.** Palantir auto-detects terminals. Destinies are named after the fact from barrier / IFN profiles.
- MAGIC imputation is not used for tests.

## Statistics

- Inferential unit for primary tests = **patient** (post-treatment, occupancy floor ≥20 A3-malignant cells).
- Cell-level Spearman is exploratory (pseudoreplication).
- Destiny vs MPR: patient-mean fate probability, NMPR vs MPR. If MPR n&lt;2 after the floor, report as descriptive and do not claim a test.
- Programs along destinies: sample-level Spearman of CLDN4, barrier (no CLDN4), and IFN ISG vs each fate probability and vs Palantir pseudotime.
- Small n is stated in the header. Do not cite n=12 as the tested n.

## Non-claims

- Not a TACSTD2∩CLDN4 dual-high analysis.
- Not CopyKAT-malignant.
- Not “AT2 becomes tumor” (AT2 markers were zeroed by the A3 rule; there is no nLung AT2 root).
- Not RNA velocity.
- Palantir destinies on a multi-patient tumor subset can recover **patient identity**. That is reported, not hidden.
