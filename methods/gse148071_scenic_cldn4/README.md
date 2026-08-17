# methods/gse148071_scenic_cldn4

Additive **CLDN4-only** SCENIC/GRN **proxy** on public **GSE148071** (Wu et al., *Nat Commun* 2021, PMID 33953163).

This folder does **not** use GSE207422 or GSE131907. It does **not** define dual-high (TACSTD2 is never a gate).

**Method actually run:** Aibar-style **AUCell** on public TF–target priors (TRRUST v2 + DoRothEA + CollecTRI). **Not** full pySCENIC (`grn → ctx → aucell`); cisTarget motif rankings were not downloaded.

Question (patient-level): do CLDN4-high putative-malignant cells show different **IFN / MHC / TJ** regulon activity vs CLDN4-low?

See `FINDING.md` and `results/primary_table.tsv`.

## Reproduce

```bash
python3 methods/gse148071_scenic_cldn4/scripts/download.py
python3 methods/gse148071_scenic_cldn4/scripts/download_priors.py
python3 methods/gse148071_scenic_cldn4/scripts/analyze.py
```
