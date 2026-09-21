# LIANA+ / CellPhoneDB / Connectome — concordant-4 CLDN4

Outgoing ligand–receptor scores from malignant CLDN4 Q4 vs Q1 to T/NK and myeloid.
Concordant four only (GSE123902, GSE131907, GSE205335, GSE189357). CLDN4 only. No dual-high.

```bash
bash methods/liana_consensus_concordant4_cldn4/scripts/download.sh /tmp/concordant4_raw
python3 methods/liana_consensus_concordant4_cldn4/scripts/run_consensus.py
```

Requires liana 1.10, anndata, scipy, pandas, statsmodels, matplotlib, and R with Matrix
(GSE205335 is a double-gzipped dgCMatrix). See METHODS.md and FINDING.md.
