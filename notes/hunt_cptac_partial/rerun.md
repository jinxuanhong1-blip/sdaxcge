# hunt_cptac_partial — rerun

```bash
python3 scripts/hunt_cptac_partial/00_download.py --outdir data/hunt_cptac_partial
python3 scripts/hunt_cptac_partial/01_analyze.py --data data/hunt_cptac_partial --outdir results/hunt_cptac_partial
```

Raw freeze files stay in `data/hunt_cptac_partial/` and are not committed.
Outputs go only to `results/hunt_cptac_partial/`.
