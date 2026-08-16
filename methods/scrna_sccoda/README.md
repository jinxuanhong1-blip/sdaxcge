# methods/scrna_sccoda

Additive public-only **Dirichlet-multinomial / scCODA** analysis of lung-tumor
scRNA cell-type composition versus **patient malignant TACSTD2** (and,
separately, **MPR**).

**Combinatorial question:** which **cohort × annotation** recovers **T/NK or
TLS down** in TACSTD2-high, with honest FDR and honest n (patients).

| | |
|---|---|
| Playbook | [`playbook.md`](playbook.md) |
| Cohorts | GSE207422, GSE241934 (IIT+RWC), GSE253013 |
| Engine | frequentist DM-GLM (default); scCODA HMC if importable |
| Unit | patient / sample — never cell |
| Access | GEO processed files only; 9.3 GB GSE253013 RDS catalogued, not downloaded |

## Reproduce

```bash
pip install -r requirements.txt
python3 scripts/selftest_dm.py
python3 scripts/run_all.py
```

## Results (after a run)

| File | Role |
|---|---|
| `results/sample_covariates.tsv` | per-patient TACSTD2, MPR, eligibility, n |
| `results/composition_long.tsv` | cohort × annotation × cell-type counts |
| `results/grid_effects.tsv` | DM-GLM + ALR + naive fraction tests |
| `results/recovery_tnk_tls.tsv` | primary recovery table |
| `results/summary.json` | n recovered, scCODA status, FDR contract |
| `results/fig1_recovery_heatmap.png` | ALR slopes across the grid |

See `data/public_derived/PROVENANCE.md` for the two tables that are not
re-extracted in this folder (GSE207422 DRMref; GSE253013 patient metrics).
