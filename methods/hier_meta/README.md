# methods/hier_meta

Bayesian / hierarchical pooling template for **TACSTD2** and **CLDN4** across
the small **open** ICI GEO cohorts already analysed in this repository.

**Methods only.** No findings live here. The script reads `results/*.csv` if
present; otherwise it prints the input contract and exits.

**Open cohorts only.** EGA / dbGaP / DAC-controlled trial transcriptomes are
never pooled. They are listed as missing and used only for MNAR sensitivity.

| Path | Role |
| --- | --- |
| [playbook.md](playbook.md) | Bilingual (EN + 中文) methods playbook |
| [templates/hier_meta.py](templates/hier_meta.py) | Primary estimator (tau quadrature, shrinkage, bias, missing EGA) |
| [templates/hier_meta.R](templates/hier_meta.R) | metafor companion |
| [templates/selftest.py](templates/selftest.py) | Unit + demo smoke tests |
| [config/config.json](config/config.json) | Default run config |
| [example/](example/) | Synthetic inputs + real-accession catalog (no invented real effects) |

```bash
python3 methods/hier_meta/templates/hier_meta.py --results-dir results
python3 methods/hier_meta/templates/hier_meta.py --demo --genes TACSTD2,CLDN4
python3 methods/hier_meta/templates/selftest.py
```
