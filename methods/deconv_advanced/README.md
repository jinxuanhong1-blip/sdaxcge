# methods/deconv_advanced

METHODS ONLY. Numeric output belongs in `results/`.

仅方法学。数值产出写入 `results/`。

| File | Role |
|---|---|
| `playbook.md` | Advanced deconvolution landscape (BayesPrism, CIBERSORTx, MuSiC, EcoTyper, InstantDL scope correction), TACSTD2 × CD8/TLS protocol, comparison protocol. zh+en. |
| `bayesprism_ecotyper_vs_tacstd2_cldn4.md` | Head-to-head: BayesPrism $Z$ vs EcoTyper CS/CE against **both** `TACSTD2` and `CLDN4`. zh+en. |
| `config/barrier_genes.yaml` | Locked `BARRIER5` gene contract. |
| `run/scale_for_ecotyper.R` | Auditable log2 + per-dataset unit-variance scaling before EcoTyper recovery. |
| `analysis/e4_e5_seeds.R` | E4 concordance and E5 residuals (no immune association). |
| `analysis/bp_vs_ecotyper_pattern.R` | Pre-declared P-agree / P-gene-split / P-method-split / P-null classifier. |
