# CPTAC LUAD TACSTD2 protein vs MCP-counter / ESTIMATE / CIBERSORT

Self-contained rework of the PR6 LUAD protein–immune slice.

PR6 reported TACSTD2 protein vs freeze **xCell immune score** ρ = −0.309
(q = 0.024). LSCC TACSTD2 protein vs xCell was null (PR23 / PR84).
This folder asks whether the LUAD association holds when immune is defined
by **public ESTIMATE**, **public CIBERSORT**, and **MCP-counter computed
from the public RNA matrix** (Becht 2016), after a **WES purity residual**.

Treatment-naive surgical LUAD (Gillette et al. *Cell* 2020). **No ICI labels.**

```bash
pip install -r scripts/rework/CPTAC_LUAD_defs/requirements.txt
python3 scripts/rework/CPTAC_LUAD_defs/download.py
python3 scripts/rework/CPTAC_LUAD_defs/analyze.py
```

Outputs: `results/rework/CPTAC_LUAD_defs/`.
Matrices (not committed): `data/rework/CPTAC_LUAD_defs/`.
