# GSE205335 trajectory — CLDN4 and RECIST

ADDITIVE CLDN4-only PAGA + diffusion pseudotime on public GSE205335
epithelium (advanced NSCLC ICI). Patient is the unit. RECIST contrast is
PR vs PD/SD. No GSE148071. No dual-high. Does not re-audit PR #320 T/NK.

Writeup: [`FINDING.md`](FINDING.md).

```bash
pip install -r methods/gse205335_traj_cldn4_recist/requirements.txt
python3 methods/gse205335_traj_cldn4_recist/scripts/download.py --out /tmp/gse205335_traj
python3 methods/gse205335_traj_cldn4_recist/scripts/extract_epithelium.py \
  --data /tmp/gse205335_traj --out /tmp/gse205335_traj/epithelium.h5ad
python3 methods/gse205335_traj_cldn4_recist/scripts/analyze.py \
  --input /tmp/gse205335_traj/epithelium.h5ad \
  --outdir methods/gse205335_traj_cldn4_recist/results \
  --finding methods/gse205335_traj_cldn4_recist/FINDING.md
```
