#!/usr/bin/env python3
"""Cell-level secondary view of the SG vs IgG1-SN-38 contrast in E-MTAB-16849.

Pseudobulk gives only n=9 paired units, so this reports per-timepoint effect
sizes at single-cell resolution. Cells within a tumour are not independent, so
the Mann-Whitney p-values here are descriptive, not inferential; the pseudobulk
paired tests in analyze.py remain the primary inference.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, RES = os.path.join(HERE, "data"), os.path.join(HERE, "results")

PANELS = {
    "core_APM": ["B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
                 "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5", "TAPBP"],
    "core_IFN": ["STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "GBP2",
                 "GBP4", "IRF7", "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IDO1"],
    "STING_axis": ["TMEM173", "CGAS", "IRF3"],
    "sanity": ["CLDN4", "TACSTD2"],
}
GENES = [g for p in PANELS.values() for g in p]
TP_ORDER = ["0h", "6h", "24h", "48h", "72h", "120h"]
# Replicates with an SG and a control arm at the same timepoint (matches analyze.py).
PAIRED = {"0h": ["repl1"], "6h": ["repl2"], "24h": ["repl1"], "48h": ["repl1", "repl2"],
          "72h": ["repl1", "repl2"], "120h": ["repl1", "repl2"]}

c = pd.read_csv(f"{DATA}/emtab16849_panel_per_cell.csv.gz", index_col="barcode")
c = c[[t in PAIRED and r in PAIRED[t] for t, r in zip(c["timepoint"], c["replicate"])]]
lg = np.log2(c[GENES].div(c["total_counts"], axis=0) * 1e6 + 1.0)
lg[["treatment", "timepoint"]] = c[["treatment", "timepoint"]]

rows = []
for tp in TP_ORDER:
    sub = lg[lg["timepoint"] == tp]
    a, b = sub[sub["treatment"] == "SG"], sub[sub["treatment"] == "ADC"]
    for panel, genes in PANELS.items():
        sa, sb = a[genes].mean(axis=1), b[genes].mean(axis=1)
        mw = stats.mannwhitneyu(sa, sb, alternative="two-sided")
        rows.append({
            "timepoint": tp, "panel": panel, "n_SG": len(sa), "n_control": len(sb),
            "mean_SG": sa.mean(), "mean_control": sb.mean(),
            "delta": sa.mean() - sb.mean(),
            "cohens_d": (sa.mean() - sb.mean()) /
                        np.sqrt((sa.var(ddof=1) + sb.var(ddof=1)) / 2),
            "mannwhitney_p_descriptive": mw.pvalue,
        })
t = pd.DataFrame(rows)
t.to_csv(f"{RES}/emtab16849_percell_by_timepoint.csv", index=False)

gl = []
for tp in TP_ORDER:
    sub = lg[lg["timepoint"] == tp]
    a, b = sub[sub["treatment"] == "SG"], sub[sub["treatment"] == "ADC"]
    for g in GENES:
        gl.append({"timepoint": tp, "gene": g,
                   "delta_log2": a[g].mean() - b[g].mean(),
                   "pct_detected_SG": float((a[g] > 0).mean()),
                   "pct_detected_control": float((b[g] > 0).mean())})
gt = pd.DataFrame(gl)
gt.to_csv(f"{RES}/emtab16849_percell_gene_by_timepoint.csv", index=False)

print(t.pivot(index="timepoint", columns="panel", values="delta")
       .reindex(TP_ORDER).to_string(float_format=lambda v: f"{v: .4f}"))
print("\nfull cell-level table:")
print(t.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
print("\ndetection rate (all paired cells):")
print(gt.groupby("gene")[["pct_detected_SG", "pct_detected_control"]].mean()
        .loc[GENES].to_string(float_format=lambda v: f"{v: .3f}"))
