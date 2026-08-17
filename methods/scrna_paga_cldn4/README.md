# Additive PAGA / diffusion scored by CLDN4

Public GSE131907 LUAD epithelium (Kim et al. 2020). **Not a TACSTD2 redo.**

- Playbook: this folder
- Finding: [`FINDING.md`](FINDING.md)
- Trajectory figure: [`results/figures/fig_trajectory_cldn4.png`](results/figures/fig_trajectory_cldn4.png)
- Extra figure: [`results/figures/fig_extra_cldn4_tertile.png`](results/figures/fig_extra_cldn4_tertile.png)

```bash
pip install -r methods/scrna_paga_cldn4/requirements.txt
bash methods/scrna_paga_cldn4/scripts/download.sh /tmp/scrna_paga_cldn4_data
python3 methods/scrna_paga_cldn4/scripts/extract_epithelium.py \
  --data /tmp/scrna_paga_cldn4_data \
  --out /tmp/scrna_paga_cldn4_data/epithelium.h5ad
python3 methods/scrna_paga_cldn4/scripts/analyze_paga_cldn4.py \
  --input /tmp/scrna_paga_cldn4_data/epithelium.h5ad \
  --outdir methods/scrna_paga_cldn4/results \
  --finding methods/scrna_paga_cldn4/FINDING.md
```

GSE207422 is referenced in-repo but not pooled (NSCLC mixed histology; no author epithelial labels on GEO).
