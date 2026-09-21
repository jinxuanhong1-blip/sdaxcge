NicheNet ligand-activity delta for CLDN4-high versus CLDN4-low malignant senders, concordant-4 only.

From the repository root:

```bash
bash methods/nichenet_activity_delta_max/scripts/download_prior.sh /tmp/nichenet_prior
Rscript methods/nichenet_activity_delta_max/scripts/export_prior_scores.R \
  --prior=/tmp/nichenet_prior/ligand_target_matrix_nsga2r_final.rds \
  --in=methods/nichenet_activity_delta_max/inputs \
  --out=/tmp/nichenet_work/prior_scores
PRIOR_SCORES=/tmp/nichenet_work/prior_scores \
  python3 methods/nichenet_activity_delta_max/scripts/activity_delta.py
```

Writeup: `RESULTS.md`.
