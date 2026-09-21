# Milo/scVI: TACSTD2 T/NK neighbourhoods conditional on CLDN4

Concordant-4 only (GSE123902, GSE131907, GSE205335, GSE189357; n = 65).
The graph is the patient-batch scVI latent. The question is whether
TACSTD2-associated T/NK-down neighbourhoods shrink once CLDN4 is in the model.

Raw matrices are not committed. Download them with the concordant-4 prepare script.

```bash
python3 methods/scanvi_concordant4_cldn4/download.py --out /tmp/geo_c4
python3 methods/scanvi_concordant4_cldn4/prepare.py --raw /tmp/geo_c4
python3 methods/milo_scvi_tacstd2_cldn4/train_scvi.py
python3 methods/milo_scvi_tacstd2_cldn4/run_milo.py
```

edgeR (the Milo quasi-likelihood test) is an R dependency. Numbers are written by
`run_milo.py` into `FINDING.md` from the fitted tables.
