# Additive GSE148071 PAGA / diffusion scored by CLDN4

Public advanced NSCLC epithelium (Wu et al. 2021, GSE148071). **Not a TACSTD2 redo.**

- Finding: [`FINDING.md`](FINDING.md)
- Trajectory figure: [`figures/fig_trajectory_cldn4.png`](figures/fig_trajectory_cldn4.png)
- Extra figure: [`figures/fig_extra_cldn4_tertile.png`](figures/fig_extra_cldn4_tertile.png)

```bash
pip install -r methods/gse148071_paga_cldn4/requirements.txt
bash methods/gse148071_paga_cldn4/scripts/download.sh /tmp/gse148071_paga_data
python3 methods/gse148071_paga_cldn4/scripts/extract_epithelium.py \
  --data /tmp/gse148071_paga_data \
  --out /tmp/gse148071_paga_data/epithelium.h5ad
python3 methods/gse148071_paga_cldn4/scripts/analyze_paga_cldn4.py \
  --input /tmp/gse148071_paga_data/epithelium.h5ad \
  --outdir methods/gse148071_paga_cldn4 \
  --finding methods/gse148071_paga_cldn4/FINDING.md
```

Honest n is written by the analyze script. GEO n=42 is the catalog, not the Spearman n.
