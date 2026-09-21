# GSE131907 pySCENIC, CLDN4-high vs low

Real pySCENIC (GRNBoost2 + cisTarget + AUCell) on author-malignant cells. Prior folders in this project that say “SCENIC proxy” did not run cisTarget.

Results: [`results/gse131907_pyscenic_cldn4/RESULTS.md`](../../results/gse131907_pyscenic_cldn4/RESULTS.md) (written by `scripts/03_contrast.py`).

Protocol: [`METHODS.md`](METHODS.md).

```bash
python3 methods/gse131907_pyscenic_cldn4/scripts/01_prepare.py
python3 methods/gse131907_pyscenic_cldn4/scripts/02_pyscenic.py
python3 methods/gse131907_pyscenic_cldn4/scripts/03_contrast.py
```

SCENIC+ is not applicable (no scATAC). Concordant-4 is not merged (15 GiB RAM).
