# MAX EFFECT Part1 — public mouse + GSE137244 Tacstd2

Maximize (1) bulk Tacstd2 KL>KP Δ and Welch |t|, and (2) scRNA Tacstd2-high → TJ module concordance (mice up/down). Every row is labeled bulk or scRNA. Public GEO only; no private 8KL.

```bash
# bulk (downloads into methods/.../cache/)
python3 methods/max_effect_part1_public_mouse_tacstd2/scripts/analyze_bulk_tacstd2_max.py

# scRNA
bash methods/max_effect_part1_public_mouse_tacstd2/scripts/download_scrna.sh /tmp/kpkl_10x
python3 methods/max_effect_part1_public_mouse_tacstd2/scripts/analyze_scrna_tj_concordance.py \
  --data /tmp/kpkl_10x \
  --out methods/max_effect_part1_public_mouse_tacstd2
```

Write-up: `FINDING.md`.
