# methods/scrna_sccoda

Additive public scRNA **compositional** analysis for malignant TACSTD2 vs T/NK
and NMPR vs MPR epithelial / T/NK fractions.

Series (each alone, then combinations): **GSE207422**, **GSE241934**, **GSE291670**.

scCODA was not installable here (`sccoda` → rpy2 → system R). The fallback is
**CLR + ILR** on 6-part lineage counts and a **Dirichlet-multinomial** two-group
LRT with permutation p. Patient is the unit. Honest n/p only.

## Run

```bash
pip install -r methods/scrna_sccoda/requirements.txt
python3 methods/scrna_sccoda/scripts/download.py
python3 methods/scrna_sccoda/scripts/extract_gse207422_markers.py
python3 methods/scrna_sccoda/scripts/build_compositions.py
python3 methods/scrna_sccoda/scripts/run_analysis.py
```

Large GEO matrices stay under `data/` and are gitignored.

## Headline (computed, not claimed)

The only slice whose **directions** match the user trend (T/NK down when
malignant TACSTD2 is high, **and** NMPR higher epithelial / lower T/NK) is
**GSE241934 REAL** (EGFR-WT real-world neoadjuvant IO+chemo).

That H1 correlation is **ρ = −0.018, p = 0.93, n = 24** — a null. H2
epithelial one-sided p = 0.10 (21 vs 13) and is partly the residual-tumor
definition of MPR. No multi-series combination keeps both directions.

See `results/FINDING.md`.
