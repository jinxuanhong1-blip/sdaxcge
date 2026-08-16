# bulk_immune

Runnable 2024–2026 playbook for bulk-RNA immune deconvolution and exclusion
scores, and their correlation with **TACSTD2** / **CLDN4**, in NSCLC ICI
cohorts plus TCGA LUAD/LUSC.

Start here: **[playbook.md](playbook.md)** (English + 中文).

```bash
cd methods/bulk_immune
pip install -r requirements.txt
bash scripts/05_run_demo.sh
```

Outputs land in `results/<cohort>/`. The numbers quoted in the playbook are
copied to `results/demo/` so they remain auditable without re-downloading GEO
or Xena.
