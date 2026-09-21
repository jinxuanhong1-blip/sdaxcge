# Concordant-4 malignant trajectories (Palantir, PAGA, Slingshot)

Malignant cells from the locked four cohorts:

- GSE123902 (donor)
- GSE131907 (sample; author `Malignant cells` only)
- GSE205335 (patient; author `Malignant cells`, non-normal tissue)
- GSE189357 (patient)

The path is **AT2-like → barrier-like** inside those malignant cells. **CLDN4 is tracked along the path.** It is not used to choose the root or the barrier terminal.

**This is not an ICI clock.** Pseudotime is cell state, not time on immunotherapy. RECIST and MPR are not tested.

The locked CLDN4 %pos vs T/NK result is not recomputed. This folder tests whether **terminal barrier-fate probability** tracks the locked patient T/NK fraction.

```bash
pip install -r methods/concordant4_malignant_traj_cldn4/requirements.txt
python3 methods/concordant4_malignant_traj_cldn4/scripts/download.py --out /tmp/c4raw
python3 methods/concordant4_malignant_traj_cldn4/scripts/extract_malignant.py \
  --data /tmp/c4raw --out /tmp/c4_malignant.h5ad
python3 methods/concordant4_malignant_traj_cldn4/scripts/analyze.py \
  --input /tmp/c4_malignant.h5ad
```

Write-up: `FINDING.md`. Methods: `METHODS.md`.
