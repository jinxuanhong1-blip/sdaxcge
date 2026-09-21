# Concordant-4: expression trajectory toward CLDN4-high

CLDN4-only. Locked units: GSE123902 (13) + GSE131907 (21) + GSE205335 (22) + GSE189357 (9).

The public matrices do not contain spliced and unspliced counts, so scVelo is not fit here. The analysis that did run is an expression-graph trajectory (PAGA, DPT, an MST in the Slingshot style, a CytoTRACE-like potency, cluster absorption) plus AUCell, UCell, z-mean, and AddModuleScore. It is not splicing velocity.

```bash
python3 methods/concordant4_scrna_velocity/inventory.py
python3 methods/concordant4_scrna_velocity/prepare_expression.py   # needs the GEO files under /tmp/c4data
python3 methods/concordant4_scrna_velocity/sweep_trajectory.py
python3 methods/concordant4_scrna_velocity/plot_winner.py
```

Headline numbers are in `FINDING.md`. Every grid row is in `results/tables/trajectory_grid.tsv`.
