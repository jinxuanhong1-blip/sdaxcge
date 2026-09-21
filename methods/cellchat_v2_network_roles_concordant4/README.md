# CellChat v2 network and sender roles (concordant-4)

Additive to the locked 14-pair CellChat table. CLDN4-high vs CLDN4-low malignant cells are the senders. T/NK is the receiver context. Honest unit is the patient / locked sample.

Concordant four only: GSE123902, GSE131907, GSE205335, GSE189357. No GSE148071. No dual-high. TACSTD2 is not a gate.

```bash
Rscript methods/cellchat_v2_network_roles_concordant4/scripts/install_packages.R
bash methods/cellchat_v2_network_roles_concordant4/scripts/download.sh /tmp/concordant4_raw
Rscript methods/cellchat_v2_network_roles_concordant4/scripts/run_cellchat_v2.R --raw=/tmp/concordant4_raw
```

Primary numbers are written by the R script into `FINDING.md` and `results/tables/`.
