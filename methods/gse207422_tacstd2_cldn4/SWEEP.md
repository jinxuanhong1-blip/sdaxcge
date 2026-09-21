# GSE207422 immune-endpoint method sweep

One additional sweep for TACSTD2 and CLDN4 versus CD8 and versus MPR/NMPR.
Rules were fixed in `scripts/sweep_immune.py` before these p-values were computed.
A positive finding requires the epithelial-mean test and either percent-positive or the epithelial-fraction partial correlation to agree at P<0.05.
Quartile, histology-only, and leave-one-out cuts do not reopen a closed endpoint.

Closure call: **FINAL** for GSE207422 immune endpoints (CD8 and MPR).

## Primary rows (n=12 post-treatment)

- TACSTD2 | epithelial_mean | spearman | CD8 fraction: ρ=-0.23 P=0.471 n=12
- TACSTD2 | epithelial_mean | partial_spearman_adj_epithelial_fraction | CD8 fraction: ρ=-0.27 P=0.399 n=12
- TACSTD2 | epithelial_mean | spearman | CD8+NK fraction: ρ=-0.18 P=0.572 n=12
- TACSTD2 | epithelial_mean | partial_spearman_adj_epithelial_fraction | CD8+NK fraction: ρ=-0.21 P=0.504 n=12
- TACSTD2 | epithelial_mean | spearman | NK fraction: ρ=0.29 P=0.354 n=12
- TACSTD2 | epithelial_mean | partial_spearman_adj_epithelial_fraction | NK fraction: ρ=0.29 P=0.360 n=12
- TACSTD2 | epithelial_mean | spearman | CD8 within T/NK: ρ=0.01 P=0.983 n=12
- TACSTD2 | epithelial_mean | partial_spearman_adj_epithelial_fraction | CD8 within T/NK: ρ=0.01 P=0.985 n=12
- TACSTD2 | epithelial_mean | spearman | CD8 cytotoxicity: ρ=0.22 P=0.484 n=12
- TACSTD2 | epithelial_mean | partial_spearman_adj_epithelial_fraction | CD8 cytotoxicity: ρ=0.22 P=0.491 n=12
- TACSTD2 | epithelial_mean | spearman | T/NK cytotoxicity: ρ=0.10 P=0.762 n=12
- TACSTD2 | epithelial_mean | partial_spearman_adj_epithelial_fraction | T/NK cytotoxicity: ρ=0.10 P=0.749 n=12
- TACSTD2 | epithelial_mean | NMPR vs MPR: Δ=0.366 (n=8 vs 4), exact P=0.109
- TACSTD2 | epithelial_mean | CD8 fraction top3 vs bottom3: Δ median=-0.097, P=0.700 (n=3 vs 3; minimum two-sided exact P is 0.10)
- TACSTD2 | epithelial_mean | CD8 cytotoxicity top3 vs bottom3: Δ median=0.195, P=0.700 (n=3 vs 3; minimum two-sided exact P is 0.10)
- TACSTD2 | epithelial_mean | MPR top3 vs bottom3: Δ median=NA, P=1.000 (top MPR 0/3 vs bottom 1/3; OR=0; minimum informative P is large at n=3)
- CLDN4 | epithelial_mean | spearman | CD8 fraction: ρ=-0.05 P=0.880 n=12
- CLDN4 | epithelial_mean | partial_spearman_adj_epithelial_fraction | CD8 fraction: ρ=-0.07 P=0.821 n=12
- CLDN4 | epithelial_mean | spearman | CD8+NK fraction: ρ=-0.10 P=0.746 n=12
- CLDN4 | epithelial_mean | partial_spearman_adj_epithelial_fraction | CD8+NK fraction: ρ=-0.13 P=0.685 n=12
- CLDN4 | epithelial_mean | spearman | NK fraction: ρ=-0.06 P=0.863 n=12
- CLDN4 | epithelial_mean | partial_spearman_adj_epithelial_fraction | NK fraction: ρ=-0.06 P=0.852 n=12
- CLDN4 | epithelial_mean | spearman | CD8 within T/NK: ρ=-0.36 P=0.245 n=12
- CLDN4 | epithelial_mean | partial_spearman_adj_epithelial_fraction | CD8 within T/NK: ρ=-0.37 P=0.243 n=12
- CLDN4 | epithelial_mean | spearman | CD8 cytotoxicity: ρ=-0.26 P=0.417 n=12
- CLDN4 | epithelial_mean | partial_spearman_adj_epithelial_fraction | CD8 cytotoxicity: ρ=-0.26 P=0.409 n=12
- CLDN4 | epithelial_mean | spearman | T/NK cytotoxicity: ρ=-0.29 P=0.366 n=12
- CLDN4 | epithelial_mean | partial_spearman_adj_epithelial_fraction | T/NK cytotoxicity: ρ=-0.28 P=0.372 n=12
- CLDN4 | epithelial_mean | NMPR vs MPR: Δ=0.062 (n=8 vs 4), exact P=1.000
- CLDN4 | epithelial_mean | CD8 fraction top3 vs bottom3: Δ median=-0.134, P=1.000 (n=3 vs 3; minimum two-sided exact P is 0.10)
- CLDN4 | epithelial_mean | CD8 cytotoxicity top3 vs bottom3: Δ median=-0.827, P=0.700 (n=3 vs 3; minimum two-sided exact P is 0.10)
- CLDN4 | epithelial_mean | MPR top3 vs bottom3: Δ median=NA, P=1.000 (top MPR 1/3 vs bottom 1/3; OR=1; minimum informative P is large at n=3)
- TACSTD2 | epithelial_pct | spearman | CD8 fraction: ρ=-0.27 P=0.404 n=12
- TACSTD2 | epithelial_pct | partial_spearman_adj_epithelial_fraction | CD8 fraction: ρ=-0.23 P=0.478 n=12
- TACSTD2 | epithelial_pct | spearman | CD8+NK fraction: ρ=-0.24 P=0.457 n=12
- TACSTD2 | epithelial_pct | partial_spearman_adj_epithelial_fraction | CD8+NK fraction: ρ=-0.20 P=0.536 n=12
- TACSTD2 | epithelial_pct | spearman | NK fraction: ρ=0.39 P=0.208 n=12
- TACSTD2 | epithelial_pct | partial_spearman_adj_epithelial_fraction | NK fraction: ρ=0.41 P=0.185 n=12
- TACSTD2 | epithelial_pct | spearman | CD8 within T/NK: ρ=-0.09 P=0.779 n=12
- TACSTD2 | epithelial_pct | partial_spearman_adj_epithelial_fraction | CD8 within T/NK: ρ=-0.09 P=0.781 n=12
- TACSTD2 | epithelial_pct | spearman | CD8 cytotoxicity: ρ=0.03 P=0.914 n=12
- TACSTD2 | epithelial_pct | partial_spearman_adj_epithelial_fraction | CD8 cytotoxicity: ρ=0.05 P=0.889 n=12
- TACSTD2 | epithelial_pct | spearman | T/NK cytotoxicity: ρ=0.06 P=0.863 n=12
- TACSTD2 | epithelial_pct | partial_spearman_adj_epithelial_fraction | T/NK cytotoxicity: ρ=0.04 P=0.890 n=12
- TACSTD2 | epithelial_pct | NMPR vs MPR: Δ=1.947 (n=8 vs 4), exact P=1.000
- TACSTD2 | epithelial_pct | CD8 fraction top3 vs bottom3: Δ median=-0.097, P=0.400 (n=3 vs 3; minimum two-sided exact P is 0.10)
- TACSTD2 | epithelial_pct | CD8 cytotoxicity top3 vs bottom3: Δ median=-0.077, P=1.000 (n=3 vs 3; minimum two-sided exact P is 0.10)
- TACSTD2 | epithelial_pct | MPR top3 vs bottom3: Δ median=NA, P=1.000 (top MPR 1/3 vs bottom 1/3; OR=1; minimum informative P is large at n=3)
- CLDN4 | epithelial_pct | spearman | CD8 fraction: ρ=-0.30 P=0.342 n=12
- CLDN4 | epithelial_pct | partial_spearman_adj_epithelial_fraction | CD8 fraction: ρ=-0.28 P=0.372 n=12
- CLDN4 | epithelial_pct | spearman | CD8+NK fraction: ρ=-0.30 P=0.342 n=12
- CLDN4 | epithelial_pct | partial_spearman_adj_epithelial_fraction | CD8+NK fraction: ρ=-0.28 P=0.372 n=12
- CLDN4 | epithelial_pct | spearman | NK fraction: ρ=0.32 P=0.308 n=12
- CLDN4 | epithelial_pct | partial_spearman_adj_epithelial_fraction | NK fraction: ρ=0.33 P=0.292 n=12
- CLDN4 | epithelial_pct | spearman | CD8 within T/NK: ρ=-0.30 P=0.342 n=12
- CLDN4 | epithelial_pct | partial_spearman_adj_epithelial_fraction | CD8 within T/NK: ρ=-0.30 P=0.342 n=12
- CLDN4 | epithelial_pct | spearman | CD8 cytotoxicity: ρ=-0.28 P=0.379 n=12
- CLDN4 | epithelial_pct | partial_spearman_adj_epithelial_fraction | CD8 cytotoxicity: ρ=-0.28 P=0.386 n=12
- CLDN4 | epithelial_pct | spearman | T/NK cytotoxicity: ρ=-0.16 P=0.618 n=12
- CLDN4 | epithelial_pct | partial_spearman_adj_epithelial_fraction | T/NK cytotoxicity: ρ=-0.17 P=0.598 n=12
- CLDN4 | epithelial_pct | NMPR vs MPR: Δ=-3.906 (n=8 vs 4), exact P=0.570
- CLDN4 | epithelial_pct | CD8 fraction top3 vs bottom3: Δ median=-0.247, P=0.200 (n=3 vs 3; minimum two-sided exact P is 0.10)
- CLDN4 | epithelial_pct | CD8 cytotoxicity top3 vs bottom3: Δ median=-0.827, P=0.200 (n=3 vs 3; minimum two-sided exact P is 0.10)
- CLDN4 | epithelial_pct | MPR top3 vs bottom3: Δ median=NA, P=1.000 (top MPR 1/3 vs bottom 1/3; OR=1; minimum informative P is large at n=3)

## Histology strata (n=6 each; descriptive)

- Adeno | TACSTD2 | epithelial_mean | CD8 fraction: ρ=0.20 P=0.704 n=6
- Adeno | TACSTD2 | epithelial_mean | NMPR vs MPR: Δ=-0.193 (n=4 vs 2), P=1.000
- Squamous | TACSTD2 | epithelial_mean | CD8 fraction: ρ=-0.37 P=0.468 n=6
- Squamous | TACSTD2 | epithelial_mean | NMPR vs MPR: Δ=0.823 (n=4 vs 2), P=0.133
- Adeno | CLDN4 | epithelial_mean | CD8 fraction: ρ=0.43 P=0.397 n=6
- Adeno | CLDN4 | epithelial_mean | NMPR vs MPR: Δ=-0.256 (n=4 vs 2), P=0.533
- Squamous | CLDN4 | epithelial_mean | CD8 fraction: ρ=-0.43 P=0.397 n=6
- Squamous | CLDN4 | epithelial_mean | NMPR vs MPR: Δ=0.396 (n=4 vs 2), P=0.533
- Adeno | TACSTD2 | epithelial_pct | CD8 fraction: ρ=0.14 P=0.787 n=6
- Adeno | TACSTD2 | epithelial_pct | NMPR vs MPR: Δ=-16.383 (n=4 vs 2), P=0.267
- Squamous | TACSTD2 | epithelial_pct | CD8 fraction: ρ=-0.71 P=0.111 n=6
- Squamous | TACSTD2 | epithelial_pct | NMPR vs MPR: Δ=14.914 (n=4 vs 2), P=0.267
- Adeno | CLDN4 | epithelial_pct | CD8 fraction: ρ=0.31 P=0.544 n=6
- Adeno | CLDN4 | epithelial_pct | NMPR vs MPR: Δ=-13.740 (n=4 vs 2), P=0.133
- Squamous | CLDN4 | epithelial_pct | CD8 fraction: ρ=-0.83 P=0.042 n=6
- Squamous | CLDN4 | epithelial_pct | NMPR vs MPR: Δ=8.350 (n=4 vs 2), P=0.533

## Leave-one-out, epithelial mean

- TACSTD2 | CD8 fraction Spearman ρ range -0.41 to 0.00 (P 0.212 to 1.000).
- TACSTD2 | CD8+NK fraction Spearman ρ range -0.35 to 0.06 (P 0.298 to 0.979).
- TACSTD2 | MPR Mann-Whitney P range 0.042 to 0.279 when one patient is dropped (full-sample P is the primary result).
- CLDN4 | CD8 fraction Spearman ρ range -0.20 to 0.12 (P 0.555 to 0.958).
- CLDN4 | CD8+NK fraction Spearman ρ range -0.27 to 0.05 (P 0.417 to 0.979).
- CLDN4 | MPR Mann-Whitney P range 0.630 to 1.000 when one patient is dropped (full-sample P is the primary result).

## Call

- TACSTD2: CD8 closed=True (largest |ρ| among the four primary CD8 rows 0.27); MPR closed=True (mean Δ=0.366, P=0.109; percent-positive Δ=1.947, P=1.000).
- CLDN4: CD8 closed=True (largest |ρ| among the four primary CD8 rows 0.26); MPR closed=True (mean Δ=0.062, P=1.000; percent-positive Δ=-3.906, P=0.570).

## Rows with P<0.05 that do not reopen the endpoint

- Squamous | CLDN4 | epithelial_pct | CD8 fraction: ρ=-0.83 P=0.042 n=6. Histology strata were pre-specified as descriptive.
- Same stratum, epithelial mean: ρ=-0.43 P=0.397 n=6. The percent-positive stratum does not agree with the mean at P<0.05.
- Leave-one-out TACSTD2 epithelial-mean MPR: dropping P02, P12 moves P to 0.042. The full-sample P stays the primary result.
These rows are reported because they are the only places a P falls under 0.05. They are n=6 histology or single-patient deletions. They do not meet the concordance rule.

Immune endpoints on this matrix are **FINAL**. CD8 abundance, CD8 cytotoxicity, and MPR do not meet the pre-specified concordance rule for TACSTD2 or for CLDN4. Further cuts of GSE207422 will not be treated as new evidence for these endpoints.

