# GSE131907 max effect: TJ Δ and TACSTD2 vs T

Patient-level sweep on the public Kim 2020 UMI matrix. The inferential n is patients. Cell-level tests are not used.

The script first checks that author-malignant CLDN4 % positive versus T/NK matches the locked GSE131907 row (sample ρ ≈ −0.522, patient ρ ≈ −0.478, n = 21). If that check fails, it stops.

Winning rules are fixed in `scripts/analyze_maxeffect.py` and are not a confirmatory p-value. See `FINDING.md`.

```bash
bash methods/gse131907_maxeffect_tj_t/scripts/download.sh /tmp/gse131907
python3 methods/gse131907_maxeffect_tj_t/scripts/extract_genes.py
python3 methods/gse131907_maxeffect_tj_t/scripts/analyze_maxeffect.py
```

Raw matrices stay in `/tmp`. They are not committed.
