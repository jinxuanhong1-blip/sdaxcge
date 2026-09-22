# Extended Data — GSE131907 tight junction versus T cells

Patient-level figure for the GSE131907 maximum-effect sweep (PR #773). Claudin scores rise from TACSTD2 Q1 to Q4. TACSTD2 percent does not track T-cell fraction.

Regenerate the PDF, SVG and PNG, plus `CAPTION.md`, with:

```bash
python3 nature_figures/ed_gse131907_tj/plot_ed_gse131907_tj.py
```

Requires matplotlib, numpy, pandas and scipy. The typeface is Arimo, the metric-compatible Arial face available in this environment.

Source tables are copied from PR #773. See `source/PROVENANCE.md`. The script recomputes the patient-level means and correlations and stops if they disagree with the published values. Inferential n is patients.
