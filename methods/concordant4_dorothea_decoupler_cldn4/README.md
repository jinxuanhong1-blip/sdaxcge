# Concordant-4 DoRothEA / decoupler — malignant CLDN4-high vs low

Patient-level TF activity on the locked four-cohort malignant pseudobulks (PR #503 labels). Primary method is decoupler ULM on DoRothEA A+B+C.

STAT1 / STAT2 / IRF1 activity is lower in CLDN4-high (same sign as IFN-up after CLDN4 loss, and the same sign as the IFN gene score). The cell-level DoRothEA wmean in PR #454 is the opposite sign and is a different contrast. NLRC5 is not in DoRothEA. Junction TFs are not different at this n.

Writeup: `FINDING.md`. Methods: `METHODS.md`.

A specification sweep (`scripts/sweep.py`) is reported in full. The smallest joint p-value is CollecTRI GRHL2 (10 targets) up in CLDN4-high together with STAT1 down. TFAP2A in that same network is significantly the other way, and the pre-specified DoRothEA GRHL2 test stays null. Those p-values are from a search, not a second primary test.

```bash
python3 methods/concordant4_dorothea_decoupler_cldn4/scripts/analyze.py
python3 methods/concordant4_dorothea_decoupler_cldn4/scripts/sweep.py
```
