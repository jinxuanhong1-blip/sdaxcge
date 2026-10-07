#!/usr/bin/env python3
"""Score APM / type-I-II IFN / STING panels in E-MTAB-16849.

Contrast: SG (sacituzumab govitecan, TROP2-ADC) vs ADC (IgG1-SN-38 non-targeting
control carrying the same payload), paired within timepoint x replicate.

log2(x+1) on CPM pseudobulk; log2FC = treated - control; gene-level Wilcoxon
signed-rank across pairs; per-sample panel score compared with a paired t-test.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

PANELS = {
    "core_APM": ["B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
                 "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5", "TAPBP"],
    "core_IFN": ["STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "GBP2",
                 "GBP4", "IRF7", "ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "IDO1"],
    "STING_axis": ["TMEM173", "CGAS", "IRF3"],
}
SANITY = ["CLDN4", "TACSTD2"]
ALL_GENES = [g for p in PANELS.values() for g in p] + SANITY
TP_ORDER = ["0h", "6h", "24h", "48h", "72h", "120h"]

pb = pd.read_csv(f"{DATA}/emtab16849_pseudobulk_cpm.csv", index_col="group")
cpm = pb[[f"{g}__cpm" for g in ALL_GENES]]
cpm.columns = ALL_GENES
lg = np.log2(cpm + 1.0)
lg = pd.concat([pb[["treatment", "timepoint", "replicate", "n_cells"]], lg], axis=1)

# Pair SG with the non-targeting control sharing the same timepoint and replicate.
key = lg["timepoint"].astype(str) + "_" + lg["replicate"].astype(str)
sg = lg[lg["treatment"] == "SG"].set_index(key[lg["treatment"] == "SG"])
ct = lg[lg["treatment"] == "ADC"].set_index(key[lg["treatment"] == "ADC"])
pairs = [k for k in sg.index if k in ct.index]
pairs.sort(key=lambda k: (TP_ORDER.index(k.split("_")[0]), k))
unpaired = sorted(set(sg.index) ^ set(ct.index))

fc = pd.DataFrame({g: sg.loc[pairs, g].values - ct.loc[pairs, g].values
                   for g in ALL_GENES}, index=pairs)
fc.index.name = "pair_timepoint_replicate"
fc.insert(0, "timepoint", [p.split("_")[0] for p in pairs])
fc.to_csv(f"{RES}/emtab16849_log2fc_per_pair.csv")


def gene_table(sub: pd.DataFrame, tag: str) -> pd.DataFrame:
    recs = []
    for panel, genes in list(PANELS.items()) + [("sanity", SANITY)]:
        for g in genes:
            d = sub[g].to_numpy(float)
            try:
                w = stats.wilcoxon(d, alternative="two-sided",
                                   zero_method="wilcox").pvalue
            except ValueError:
                w = np.nan
            mean_cpm = float(cpm[g].mean())
            recs.append({
                "panel": panel, "gene": g, "n_pairs": len(d),
                "mean_log2FC": d.mean(), "median_log2FC": float(np.median(d)),
                "n_pos": int((d > 0).sum()), "n_neg": int((d < 0).sum()),
                "wilcoxon_p": w, "mean_cpm_all_samples": mean_cpm,
            })
    t = pd.DataFrame(recs)
    t.to_csv(f"{RES}/emtab16849_gene_level_{tag}.csv", index=False)
    return t


def score_table(sub: pd.DataFrame, tag: str) -> pd.DataFrame:
    recs = []
    for panel, genes in list(PANELS.items()) + [("sanity", SANITY)]:
        sg_s = sg.loc[sub.index, genes].mean(axis=1)
        ct_s = ct.loc[sub.index, genes].mean(axis=1)
        d = (sg_s - ct_s).to_numpy(float)
        tt = stats.ttest_rel(sg_s, ct_s)
        ws = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox")
        recs.append({
            "panel": panel, "n_pairs": len(d),
            "mean_score_SG": sg_s.mean(), "mean_score_control": ct_s.mean(),
            "mean_delta_score": d.mean(), "sd_delta": d.std(ddof=1),
            "n_pos": int((d > 0).sum()),
            "t_stat": tt.statistic, "t_p": tt.pvalue,
            "wilcoxon_p": ws.pvalue,
        })
    t = pd.DataFrame(recs)
    t.to_csv(f"{RES}/emtab16849_panel_scores_{tag}.csv", index=False)
    return t


# Per-timepoint panel score deltas, for the time course.
rows = []
for p in pairs:
    r = {"pair": p, "timepoint": p.split("_")[0]}
    for panel, genes in list(PANELS.items()) + [("sanity", SANITY)]:
        r[f"delta_{panel}"] = float(sg.loc[p, genes].mean() - ct.loc[p, genes].mean())
    for g in SANITY:
        r[f"log2FC_{g}"] = float(sg.loc[p, g] - ct.loc[p, g])
    rows.append(r)
per_tp = pd.DataFrame(rows)
per_tp.to_csv(f"{RES}/emtab16849_score_delta_per_pair.csv", index=False)

print(f"paired SG-vs-control contrasts (n={len(pairs)}): {pairs}")
print(f"unpaired groups dropped: {unpaired}\n")

for tag, sub in [("all_pairs", fc), ("post_dose", fc[fc["timepoint"] != "0h"])]:
    print("=" * 78)
    print(f"[{tag}] n_pairs={len(sub)}")
    gt = gene_table(sub, tag)
    st = score_table(sub, tag)
    print(gt.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
    print()
    print(st.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
    print()

print("=" * 78)
print("per-pair score deltas")
print(per_tp.to_string(index=False, float_format=lambda v: f"{v: .4f}"))
