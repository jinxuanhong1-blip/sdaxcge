# Part1 bridge max ρ: DepMap/CCLE TROP2–CLDN4 protein + TROP2–TJ gene

Maximize Spearman ρ for (1) Gygi TACSTD2–CLDN4 protein and (2) DepMap 24Q4
TACSTD2 versus locked TJ gene scores, with pre-set histology / partial grids and
honest n. No imputation. No cell-line dropping to raise ρ.

```bash
python3 -m pip install -r scripts/depmap_trop2_tj_part1_maxrho/requirements.txt
python3 scripts/depmap_trop2_tj_part1_maxrho/download.py
python3 scripts/depmap_trop2_tj_part1_maxrho/analyze.py
```

Outputs: `results/depmap_trop2_tj_part1_maxrho/`.
