#!/usr/bin/env python3
"""GSE205335: TACSTD2 (TROP2) in malignant cells, responders vs non-responders,
and malignant vs T/NK cells, aggregated per patient.

Data: author-provided UMI matrix + author-provided cell annotations
(lineage.sub == "Malignant cells" is the authors' own malignant call).
Response: RECIST from GEO sample characteristics. R = PR, NR = SD/PD, NE excluded.
"""

import gzip
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/workspace/results/hunt_gse205335"
import os

os.makedirs(OUT, exist_ok=True)

MIN_CELLS = 20  # min cells per patient per compartment to report a value

# ---------- load ----------
percell = pd.read_csv("/workspace/data/percell_tacstd2.csv.gz")
ident = pd.read_csv("/workspace/data/GSE205335_Lung_IO_CellIdentity.txt", sep="\t")
gsm = pd.read_csv("/tmp/gsm_meta.csv")

df = ident.merge(percell, on="barcode", how="inner", validate="1:1")
assert len(df) == len(ident) == len(percell), (len(df), len(ident), len(percell))

# map orig.ident -> GSM metadata via title token + platform suffix
gsm["token"] = gsm["title"].str.split().str[1].str.replace("_", "-")
gsm["suffix"] = np.where(gsm["platform"].str.contains("3'"), "3P", "5P")
gsm["orig.ident"] = gsm["token"] + "-" + gsm["suffix"]

unmapped = set(df["orig.ident"]) - set(gsm["orig.ident"])
assert not unmapped, f"unmapped samples: {unmapped}"

df = df.merge(
    gsm[["orig.ident", "gsm", "patient", "tissue", "recist", "platform", "cancer subtype"]],
    on="orig.ident",
    how="left",
)

resp_map = {"PR": "R", "CR": "R", "SD": "NR", "PD": "NR", "NE": "NE"}
df["response"] = df["recist"].map(resp_map)

# normalized expression: log1p(CP10K)
df["tacstd2_cp10k"] = df["tacstd2_umi"] / df["total_umi"] * 1e4
df["tacstd2_log1p_cp10k"] = np.log1p(df["tacstd2_cp10k"])
df["tacstd2_pos"] = df["tacstd2_umi"] > 0

# compartments
df["compartment"] = np.where(
    df["lineage.sub"] == "Malignant cells",
    "Malignant",
    np.where(df["lineage.total"] == "T/NK cells", "T/NK", "Other"),
)

# ---------- per-sample and per-patient summaries ----------
def summarize(g):
    return pd.Series(
        {
            "n_cells": len(g),
            "pct_pos": 100.0 * g["tacstd2_pos"].mean(),
            "mean_log1p_cp10k": g["tacstd2_log1p_cp10k"].mean(),
            "pseudobulk_cpm": g["tacstd2_umi"].sum() / g["total_umi"].sum() * 1e6,
        }
    )


sub = df[df["compartment"].isin(["Malignant", "T/NK"])]

per_sample = (
    sub.groupby(["patient", "orig.ident", "gsm", "tissue", "recist", "response", "platform", "compartment"])
    .apply(summarize, include_groups=False)
    .reset_index()
)
per_sample.to_csv(f"{OUT}/per_sample_tacstd2.csv", index=False)

per_patient = (
    sub.groupby(["patient", "recist", "response", "compartment"])
    .apply(summarize, include_groups=False)
    .reset_index()
)
per_patient.to_csv(f"{OUT}/per_patient_tacstd2.csv", index=False)

mal = per_patient[(per_patient["compartment"] == "Malignant") & (per_patient["n_cells"] >= MIN_CELLS)].copy()
tnk = per_patient[(per_patient["compartment"] == "T/NK") & (per_patient["n_cells"] >= MIN_CELLS)].copy()

lines = []
lines.append("GSE205335 TACSTD2 (TROP2) analysis — statistical results")
lines.append("=" * 60)
lines.append(f"Cells total: {len(df)}; malignant: {(df['compartment']=='Malignant').sum()}; "
             f"T/NK: {(df['compartment']=='T/NK').sum()}")
lines.append(f"Patients with >= {MIN_CELLS} malignant cells: {len(mal)} "
             f"(R={ (mal['response']=='R').sum() }, NR={ (mal['response']=='NR').sum() }, NE={ (mal['response']=='NE').sum() })")
all_pat = gsm[["patient", "recist"]].drop_duplicates()
mal_counts = per_patient[per_patient["compartment"] == "Malignant"][["patient", "n_cells"]]
pat_mal = all_pat.merge(mal_counts, on="patient", how="left").fillna({"n_cells": 0})
excluded = pat_mal[pat_mal["n_cells"] < MIN_CELLS]
if len(excluded):
    lines.append("Patients excluded (< min malignant cells captured): "
                 + ", ".join(f"{r.patient}({int(r.n_cells)} cells, RECIST {r.recist})" for r in excluded.itertuples()))
    n_r_drop = (excluded["recist"].isin(["PR", "CR"])).sum()
    n_nr_drop = (excluded["recist"].isin(["SD", "PD"])).sum()
    lines.append(f"    NOTE: dropout is asymmetric — {n_r_drop} responder vs {n_nr_drop} non-responder "
                 "patients had (near-)zero malignant cells captured; absence of detectable tumor cells "
                 "in responders may itself reflect treatment effect or sampling, and it limits the R-arm sample size.")

# ---------- test 1: malignant TACSTD2, R vs NR (per patient) ----------
lines.append("")
lines.append("[1] Malignant-cell TACSTD2: Responders (PR) vs Non-responders (SD/PD)")
lines.append("    Unit = patient (samples pooled per patient). NE patients excluded.")
mal_rn = mal[mal["response"].isin(["R", "NR"])]
for metric in ["mean_log1p_cp10k", "pct_pos", "pseudobulk_cpm"]:
    r_vals = mal_rn.loc[mal_rn["response"] == "R", metric].values
    n_vals = mal_rn.loc[mal_rn["response"] == "NR", metric].values
    u, p = stats.mannwhitneyu(r_vals, n_vals, alternative="two-sided")
    lines.append(
        f"    {metric}: R n={len(r_vals)} median={np.median(r_vals):.3f} | "
        f"NR n={len(n_vals)} median={np.median(n_vals):.3f} | Mann-Whitney U={u:.1f}, p={p:.4f}"
    )

# PD-only sensitivity
mal_pd = mal[mal["recist"].isin(["PR", "PD"])]
r_vals = mal_pd.loc[mal_pd["recist"] == "PR", "mean_log1p_cp10k"].values
p_vals = mal_pd.loc[mal_pd["recist"] == "PD", "mean_log1p_cp10k"].values
if len(r_vals) >= 2 and len(p_vals) >= 2:
    u, p = stats.mannwhitneyu(r_vals, p_vals, alternative="two-sided")
    lines.append(
        f"    sensitivity PR vs PD only (mean_log1p_cp10k): PR n={len(r_vals)} "
        f"median={np.median(r_vals):.3f} | PD n={len(p_vals)} median={np.median(p_vals):.3f} | p={p:.4f}"
    )

# ---------- test 2: malignant vs T/NK within patient (paired) ----------
lines.append("")
lines.append("[2] TACSTD2: Malignant vs T/NK cells, paired within patient (all patients incl. NE)")
paired = mal.merge(tnk, on="patient", suffixes=("_mal", "_tnk"))
for metric in ["mean_log1p_cp10k", "pct_pos"]:
    a = paired[f"{metric}_mal"].values
    b = paired[f"{metric}_tnk"].values
    try:
        w, p = stats.wilcoxon(a, b, alternative="two-sided")
        txt = f"Wilcoxon signed-rank W={w:.1f}, p={p:.2e}"
    except ValueError as e:
        txt = f"Wilcoxon not computable ({e})"
    lines.append(
        f"    {metric}: n={len(paired)} pairs | malignant median={np.median(a):.3f} | "
        f"T/NK median={np.median(b):.3f} | {txt}"
    )

# sanity check markers
mal_cells = df[df["compartment"] == "Malignant"]
tnk_cells = df[df["compartment"] == "T/NK"]
lines.append("")
lines.append("[sanity] EPCAM+ fraction: malignant "
             f"{100*(mal_cells['epcam_umi']>0).mean():.1f}% vs T/NK {100*(tnk_cells['epcam_umi']>0).mean():.1f}%; "
             "PTPRC(CD45)+ fraction: malignant "
             f"{100*(mal_cells['ptprc_umi']>0).mean():.1f}% vs T/NK {100*(tnk_cells['ptprc_umi']>0).mean():.1f}%")

with open(f"{OUT}/stats_results.txt", "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))

# ---------- plots ----------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# panel A: R vs NR malignant per patient
ax = axes[0]
order = ["R", "NR", "NE"]
colors = {"R": "#2166ac", "NR": "#b2182b", "NE": "#999999"}
for i, grp in enumerate(order):
    v = mal.loc[mal["response"] == grp, "mean_log1p_cp10k"].values
    if len(v) == 0:
        continue
    x = np.random.default_rng(0).normal(i, 0.06, len(v))
    ax.scatter(x, v, color=colors[grp], s=60, zorder=3, alpha=0.85)
    ax.hlines(np.median(v), i - 0.25, i + 0.25, color="black", lw=2, zorder=4)
ax.set_xticks(range(len(order)))
ax.set_xticklabels([f"{g}\n(n={ (mal['response']==g).sum() })" for g in order])
ax.set_ylabel("TACSTD2 mean log1p(CP10K), malignant cells")
ax.set_title("Malignant TACSTD2 per patient\nby ICI response (RECIST)")

# panel B: %positive R vs NR
ax = axes[1]
for i, grp in enumerate(order):
    v = mal.loc[mal["response"] == grp, "pct_pos"].values
    if len(v) == 0:
        continue
    x = np.random.default_rng(1).normal(i, 0.06, len(v))
    ax.scatter(x, v, color=colors[grp], s=60, zorder=3, alpha=0.85)
    ax.hlines(np.median(v), i - 0.25, i + 0.25, color="black", lw=2, zorder=4)
ax.set_xticks(range(len(order)))
ax.set_xticklabels([f"{g}\n(n={ (mal['response']==g).sum() })" for g in order])
ax.set_ylabel("% TACSTD2+ malignant cells")
ax.set_title("Fraction TACSTD2+ malignant cells\nper patient")

# panel C: paired malignant vs T/NK
ax = axes[2]
for r in paired.itertuples():
    ax.plot([0, 1], [r.mean_log1p_cp10k_mal, r.mean_log1p_cp10k_tnk],
            color="grey", alpha=0.5, lw=1, zorder=2)
ax.scatter([0] * len(paired), paired["mean_log1p_cp10k_mal"], color="#762a83", s=50, zorder=3, label="Malignant")
ax.scatter([1] * len(paired), paired["mean_log1p_cp10k_tnk"], color="#1b7837", s=50, zorder=3, label="T/NK")
ax.set_xticks([0, 1])
ax.set_xticklabels(["Malignant", "T/NK"])
ax.set_xlim(-0.5, 1.5)
ax.set_ylabel("TACSTD2 mean log1p(CP10K)")
ax.set_title(f"TACSTD2 malignant vs T/NK\npaired per patient (n={len(paired)})")
ax.legend(frameon=False)

plt.tight_layout()
plt.savefig(f"{OUT}/tacstd2_summary.png", dpi=200)
plt.close()

# per-patient bar detail
fig, ax = plt.subplots(figsize=(12, 5))
mo = mal.sort_values(["response", "mean_log1p_cp10k"], ascending=[True, False])
bar_colors = [colors[r] for r in mo["response"]]
ax.bar(range(len(mo)), mo["mean_log1p_cp10k"], color=bar_colors)
ax.set_xticks(range(len(mo)))
ax.set_xticklabels([f"{p}\n{r}" for p, r in zip(mo["patient"], mo["recist"])], fontsize=8)
ax.set_ylabel("TACSTD2 mean log1p(CP10K), malignant cells")
ax.set_title("Malignant TACSTD2 per patient (blue=R, red=NR, grey=NE)")
plt.tight_layout()
plt.savefig(f"{OUT}/tacstd2_per_patient_bars.png", dpi=200)
plt.close()

print(f"\nOutputs written to {OUT}")
