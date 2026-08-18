# seurat_integrate_154977_180963_cldn4

ADDITIVE pair: public mouse [GSE154977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154977) (KP 30w 10x) + [GSE180963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963) (K vs KL). **R + Seurat + Harmony**. Cldn4-only. Honest n = **mouse**. Pair, not the triple. See `FINDING.md`.

```bash
bash methods/seurat_integrate_154977_180963_cldn4/scripts/install_r.sh
bash methods/seurat_integrate_154977_180963_cldn4/scripts/download.sh /tmp/pair_154977_180963
Rscript methods/seurat_integrate_154977_180963_cldn4/scripts/analyze.R \
  --data /tmp/pair_154977_180963 \
  --out methods/seurat_integrate_154977_180963_cldn4
```
