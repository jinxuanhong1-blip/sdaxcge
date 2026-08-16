# scripts/scrna_paga

Additive public PAGA / diffusion slice. Methods: `methods/scrna_paga/playbook.md`.

```bash
pip install -r requirements.txt
bash scripts/scrna_paga/download.sh /tmp/scrna_paga_data
python3 scripts/scrna_paga/extract_epithelium.py \
  --data /tmp/scrna_paga_data \
  --out /tmp/scrna_paga_data/epithelium.h5ad
python3 scripts/scrna_paga/analyze_paga.py \
  --input /tmp/scrna_paga_data/epithelium.h5ad \
  --outdir results/scrna_paga
```

GSE253013 (9.3 GB RDS) is over the public-file cap and is not downloaded.
GSE207422 is NSCLC (LUAD+LUSC) and is not pooled into this LUAD-only graph.
