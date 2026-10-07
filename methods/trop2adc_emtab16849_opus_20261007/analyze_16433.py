#!/usr/bin/env python3
"""Score APM / IFN / STING panels in E-MTAB-16433 (SG vs vehicle, CRC PDOX).

Design is 4 vs 4 independent mice, so the primary gene-level test is the
Wilcoxon rank-sum test and the primary score test is Welch's t-test. Batches
alternate vehicle/SG (0..7), so a batch-paired signed-rank analysis is reported
as a sensitivity check. CXCL9 is absent from this study's feature list.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, RES = os.path.join(HERE, "data"), os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

pb = pd.read_csv(f"{DATA}/emtab16433_pseudobulk_cpm.csv", index_col="mouse_ID")
present = [c[:-5] for c in pb.columns if c.endswith("__cpm")]

PANELS = {
    "core_APM": ["B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
                 "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5", "TAPBP"],
    "core_IFN": ["STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "GBP2",
                 "GBP4", "IRF7", "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IDO1"],
    "STING_axis": ["TMEM173", "CGAS", "IRF3"],
    "sanity": ["CLDN4", "TACSTD2"],
}
PANELS = {k: [g for g in v if g in present] for k, v in PANELS.items()}
GENES = [g for v in PANELS.values() for g in v]

cpm = pb[[f"{g}__cpm" for g in GENES]]
cpm.columns = GENES
lg = np.log2(cpm + 1.0)
arm = pb["treatment"]
sg_i, ve_i = arm[arm == "trodelvy"].index, arm[arm == "vehicle"].index

recs = []
for panel, genes in PANELS.items():
    for g in genes:
        a, b = lg.loc[sg_i, g].to_numpy(float), lg.loc[ve_i, g].to_numpy(float)
        mw = stats.mannwhitneyu(a, b, alternative="two-sided")
        tt = stats.ttest_ind(a, b, equal_var=False)
        recs.append({
            "panel": panel, "gene": g, "n_SG": len(a), "n_vehicle": len(b),
            "mean_log2_SG": a.mean(), "mean_log2_vehicle": b.mean(),
            "log2FC": a.mean() - b.mean(),
            "ranksum_p": mw.pvalue, "welch_t_p": tt.pvalue,
            "mean_cpm_all": float(cpm[g].mean()),
        })
gt = pd.DataFrame(recs)
gt.to_csv(f"{RES}/emtab16433_gene_level.csv", index=False)

srecs = []
for panel, genes in PANELS.items():
    a, b = lg.loc[sg_i, genes].mean(axis=1), lg.loc[ve_i, genes].mean(axis=1)
    tt = stats.ttest_ind(a, b, equal_var=False)
    mw = stats.mannwhitneyu(a, b, alternative="two-sided")
    srecs.append({
        "panel": panel, "n_genes": len(genes),
        "mean_score_SG": a.mean(), "mean_score_vehicle": b.mean(),
        "delta_score": a.mean() - b.mean(),
        "sd_SG": a.std(ddof=1), "sd_vehicle": b.std(ddof=1),
        "welch_t": tt.statistic, "welch_t_p": tt.pvalue, "ranksum_p": mw.pvalue,
    })
st = pd.DataFrame(srecs)
st.to_csv(f"{RES}/emtab16433_panel_scores.csv", index=False)

# Sensitivity: batches alternate vehicle(even)/SG(odd); pair 0-1, 2-3, 4-5, 6-7.
bp = pb[["treatment", "batch"]].copy()
bp["blk"] = bp["batch"] // 2
prec = []
for panel, genes in PANELS.items():
    sc = lg[genes].mean(axis=1)
    d = []
    for blk, sub in bp.groupby("blk"):
        s = sub[sub["treatment"] == "trodelvy"].index
        v = sub[sub["treatment"] == "vehicle"].index
        if len(s) == 1 and len(v) == 1:
            d.append(float(sc.loc[s[0]] - sc.loc[v[0]]))
    d = np.asarray(d)
    prec.append({
        "panel": panel, "n_blocks": len(d), "mean_delta": d.mean(),
        "n_pos": int((d > 0).sum()),
        "signed_rank_p": stats.wilcoxon(d, zero_method="wilcox").pvalue,
        "paired_t_p": stats.ttest_1samp(d, 0.0).pvalue,
        "deltas": np.round(d, 4).tolist(),
    })
pt = pd.DataFrame(prec)
pt.to_csv(f"{RES}/emtab16433_batchpaired_sensitivity.csv", index=False)

print("per-mouse panel scores, log2(CPM+1):")
show = pd.DataFrame({p: lg[g].mean(axis=1) for p, g in PANELS.items()})
show["treatment"] = arm
print(show.sort_values("treatment").to_string(float_format=lambda v: f"{v: .4f}"))
print("\ngene level (SG vs vehicle):")
print(gt.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
print("\npanel scores (primary, unpaired 4v4):")
print(st.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
print("\nbatch-paired sensitivity:")
print(pt.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
