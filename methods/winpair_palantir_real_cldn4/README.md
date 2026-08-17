# Winning-pair REAL Palantir (GSE131907 + GSE205335), CLDN4 only

ADDITIVE CLDN4-only Palantir on the winning pair epithelium.
Does not redo GSE131907-only PAGA (PR #325) or the Slingshot/DPT fallback (PR #449).
Does not add GSE148071 or GSE207422. No TACSTD2∩CLDN4 dual-high gate.

**Install palantir. Do not stop at empty DPT.**
Done when `results/tables/palantir_destinies.tsv` and `results/tables/palantir_pseudotime.tsv` exist.

Patient is the inferential unit. Early cell is never CLDN4-high.
Scores along destinies: CLDN4, barrier (no CLDN4), IFN (no CLDN4).

```bash
pip install -r methods/winpair_palantir_real_cldn4/requirements.txt
python3 methods/winpair_palantir_real_cldn4/scripts/download.py --out /tmp/winpair_131907_205335
python3 methods/winpair_palantir_real_cldn4/scripts/extract_epithelium.py \
  --data /tmp/winpair_131907_205335 \
  --out /tmp/winpair_131907_205335/epithelium.h5ad
python3 methods/winpair_palantir_real_cldn4/scripts/analyze.py \
  --input /tmp/winpair_131907_205335/epithelium.h5ad \
  --outdir methods/winpair_palantir_real_cldn4/results \
  --finding methods/winpair_palantir_real_cldn4/FINDING.md
```
