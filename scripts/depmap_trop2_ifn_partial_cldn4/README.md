# TROP2 vs IFN / MHC / immune scores, partialled on CLDN4

Public DepMap 24Q4 RNA, Gygi/Nusinow CCLE TMT protein, and DepMap 24Q4 CRISPR Chronos.

The question is whether CLDN4 attenuates TROP2's association with pre-specified IFN, MHC-I, and immune-ligand scores. The TROP2–CLDN4 protein Spearman near 0.69 is recomputed only as a matrix check.

Results are written by `analyze.py` into `results/depmap_trop2_ifn_partial_cldn4/`. Read `FINDING.md` there. Do not retype the coefficients.

```bash
python3 -m pip install -r scripts/depmap_trop2_ifn_partial_cldn4/requirements.txt
python3 scripts/depmap_trop2_ifn_partial_cldn4/download.py --outdir data/depmap_trop2_ifn_partial_cldn4
python3 scripts/depmap_trop2_ifn_partial_cldn4/analyze.py --data data/depmap_trop2_ifn_partial_cldn4 --outdir results/depmap_trop2_ifn_partial_cldn4
```

`data/` is gitignored. The script drops the full matrices after extracting the genes used here.
