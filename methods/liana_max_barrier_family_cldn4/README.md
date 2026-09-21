# Max-effect barrier ligands (concordant-4)

CellPhoneDB `lr_means` and Connectome `expr_prod` for F11R, NECTIN2, CDH1, and LGALS9
from CLDN4-high versus CLDN4-low malignant cells to T/NK.
The grid and the selection rule are in `METHODS.md`. The numbers are in `FINDING.md`.

```bash
bash scripts/download.sh /tmp/concordant4_raw
python3 scripts/run_max_effect.py
```
