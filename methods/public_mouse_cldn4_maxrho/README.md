# Mouse-level Cldn4% on the public integrate cohorts

Maximum |Spearman ρ| of epithelial Cldn4 percent versus T fraction and versus epithelial IFN/APM. Cohorts are the three series already integrated in the public Seurat/Harmony run: GSE154977, GSE180963, GSE165641. One row is one mouse. Private 8 KL matrices are not used.

```bash
bash methods/public_mouse_cldn4_maxrho/scripts/download.sh /tmp/kpkl_10x
python3 methods/public_mouse_cldn4_maxrho/scripts/analyze.py \
  --data /tmp/kpkl_10x \
  --out methods/public_mouse_cldn4_maxrho
```

Numbers are written by the script to `FINDING.md` and `tables/`.
