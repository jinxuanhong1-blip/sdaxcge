#!/usr/bin/env python3
"""
Leftover-series TACSTD2/CLDN4 vs ICI-outcome tests.

Only series that have (a) an open processed matrix containing the gene(s)
and (b) a per-sample ICI outcome that can be joined without inventing labels:

  GSE221733  NSCLC GeoMx DSP CTA, immunotherapy-treated, Responder /
             Non-responder + follow-up/vital status. TACSTD2 is on the
             panel; CLDN4 is not. Patient-level mean of PanCK+ AOIs.

  GSE248378  NSCLC neoadjuvant durvalumab +/- SBRT, post-resection FPKM.
             Both genes present. GEO has only Arm1/Arm2; published Nature
             Communications source data (41467_2023_44195 MOESM6) maps
             DurvaNNN -> numeric sample IDs. One GEO sample (45-M-PO) has
             an Arm2 vs paper-monotherapy conflict and is excluded from
             the primary test.

  GSE193049  BALF (not tumour) RNA from PD-1-blockade responders vs
             non-responders, n=7. Both genes present.

No fabricated statistics: every p-value below is computed here from the
joined tables written to results/fable_geo_2022_2023/analysis/.
"""
import gzip
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "fable_geo_2022_2023"
DATA = OUT / "data"
ANA = OUT / "analysis"
ANA.mkdir(parents=True, exist_ok=True)

RESULTS = []


def cliffs_delta(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    gt = sum(x > y for x in a for y in b)
    lt = sum(x < y for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def mw_row(cohort, gene, outcome, g1, a, g2, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    row = {
        "cohort": cohort, "gene": gene, "outcome": outcome,
        "group1": g1, "n1": int(len(a)), "median1": round(float(np.median(a)), 3),
        "group2": g2, "n2": int(len(b)), "median2": round(float(np.median(b)), 3),
        "stat": float(u), "p_value": float(p),
        "effect_cliffs_delta": round(float(cliffs_delta(a, b)), 3),
        "note": "",
    }
    RESULTS.append(row)
    return row


def boxplot(groups, data, title, fname, ylabel):
    fig, ax = plt.subplots(figsize=(4.2, 4))
    vals = [data[g] for g in groups]
    bp = ax.boxplot(vals, tick_labels=groups, patch_artist=True, widths=0.55)
    colors = ["#4C72B0", "#C44E52"]
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, g in enumerate(groups, 1):
        y = np.asarray(data[g], float)
        x = rng.normal(i, 0.05, len(y))
        ax.scatter(x, y, color="black", s=16, zorder=3, alpha=0.7)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(ANA / fname, dpi=140)
    plt.close(fig)


def km_plot(times, events, groups, title, fname):
    fig, ax = plt.subplots(figsize=(5.2, 4))
    kmf = KaplanMeierFitter()
    labels = sorted(set(groups))
    for lab, color in zip(labels, ["#4C72B0", "#C44E52"]):
        mask = groups == lab
        kmf.fit(times[mask], events[mask], label=f"{lab} (n={int(mask.sum())})")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=color)
    a, b = labels[0], labels[1]
    res = logrank_test(times[groups == a], times[groups == b],
                       events[groups == a], events[groups == b])
    ax.set_title(f"{title}\nlog-rank p={res.p_value:.3f}", fontsize=9)
    ax.set_xlabel("time")
    ax.set_ylabel("survival")
    fig.tight_layout()
    fig.savefig(ANA / fname, dpi=140)
    plt.close(fig)
    RESULTS.append({
        "cohort": title.split()[0], "gene": "TACSTD2", "outcome": "survival_logrank",
        "group1": a, "n1": int((groups == a).sum()), "median1": "",
        "group2": b, "n2": int((groups == b).sum()), "median2": "",
        "stat": float(res.test_statistic), "p_value": float(res.p_value),
        "effect_cliffs_delta": "",
        "note": "median-split TACSTD2; log-rank",
    })


# ---------------------------------------------------------------- GSE221733
def analyze_gse221733():
    expr = pd.read_excel(DATA / "GSE221733__GSE221733_4301_CTA_norm.xlsx")
    expr = expr.rename(columns={expr.columns[0]: "Sample_title"})
    if "CLDN4" not in expr.columns:
        print("GSE221733: CLDN4 not on CTA panel (expected)", file=sys.stderr)
    meta = pd.read_csv(DATA / "GSE221733_sample_meta.csv")
    j = meta.merge(expr[["Sample_title", "TACSTD2"]], on="Sample_title", how="inner")
    j = j[j["response"].isin(["Responder", "Non-responder"])].copy()
    # patient-level mean within PanCK+ (tumour epithelium)
    pos = j[j["segment"] == "PanCK pos"]
    pat = (pos.groupby(["patient id", "response", "status"], as_index=False)
             .agg(TACSTD2=("TACSTD2", "mean"),
                  followup=("followup", "first"),
                  n_aoi=("TACSTD2", "size")))
    pat.to_csv(ANA / "GSE221733_patient_PanCKpos.csv", index=False)
    by = {g: pat.loc[pat.response == g, "TACSTD2"].tolist()
          for g in ["Responder", "Non-responder"]}
    mw_row("GSE221733", "TACSTD2", "ICI_response_PanCKpos_patient",
           "Responder", by["Responder"], "Non-responder", by["Non-responder"])
    boxplot(["Responder", "Non-responder"], by,
            f"GSE221733 NSCLC DSP PanCK+ (n={len(pat)} patients)\nTACSTD2 vs ICI response",
            "GSE221733_TACSTD2_response.png", "TACSTD2 (DSP CTA, normalized)")

    # PanCK- stroma/immune as a negative-control compartment
    neg = j[j["segment"] == "PanCK neg"]
    patn = (neg.groupby(["patient id", "response"], as_index=False)
              .agg(TACSTD2=("TACSTD2", "mean")))
    byn = {g: patn.loc[patn.response == g, "TACSTD2"].tolist()
           for g in ["Responder", "Non-responder"]}
    if all(len(v) >= 2 for v in byn.values()):
        mw_row("GSE221733", "TACSTD2", "ICI_response_PanCKneg_patient",
               "Responder", byn["Responder"], "Non-responder", byn["Non-responder"])

    # OS: follow-up (days) + vital status. Median-split TACSTD2 on PanCK+.
    surv = pat.dropna(subset=["followup", "status", "TACSTD2"]).copy()
    surv["event"] = (surv["status"].astype(str).str.lower() == "deceased").astype(int)
    med = float(surv["TACSTD2"].median())
    surv["tac_group"] = np.where(surv["TACSTD2"] >= med, "TACSTD2-high", "TACSTD2-low")
    km_plot(surv["followup"].astype(float).values, surv["event"].values,
            surv["tac_group"].values,
            "GSE221733 TACSTD2 median-split OS (PanCK+)",
            "GSE221733_TACSTD2_OS.png")
    print(f"GSE221733 patients PanCK+ {len(pat)} R/NR "
          f"{pat.response.value_counts().to_dict()}", file=sys.stderr)
    return pat


# ---------------------------------------------------------------- GSE248378
def analyze_gse248378():
    expr = pd.read_csv(DATA / "GSE248378__GSE248378_Durva_Post_FPKMs.txt.gz",
                       sep="\t", index_col=0)
    meta = pd.read_csv(DATA / "GSE248378_sample_meta.csv")
    src = pd.read_excel(DATA / "41467_2023_44195_MOESM6_ESM.xlsx",
                        sheet_name="Figure 1b")
    src = src.iloc[:, :5].copy()
    src.columns = ["ref", "paper_arm", "resection", "status", "dfs_months"]
    src = src.dropna(subset=["ref"])
    src["num"] = src["ref"].map(lambda t: int(re.search(r"(\d+)", str(t)).group(1)))
    # pandas default NA list includes the string "None", which is the
    # paper's non-MPR label — keep it as a real category.
    mpr = pd.read_excel(DATA / "41467_2023_44195_MOESM6_ESM.xlsx",
                        sheet_name="Figure 2e", keep_default_na=False)
    mpr = mpr.iloc[:, :3].copy()
    mpr.columns = ["ref", "resection2", "mpr"]
    mpr = mpr[mpr["ref"].astype(str).str.startswith("Durva")]
    src = src.merge(mpr[["ref", "mpr"]], on="ref", how="left")

    meta["num"] = meta["Sample_title"].map(
        lambda t: int(re.search(r"(\d+)", str(t)).group(1)))
    j = meta.merge(src, on="num", how="left")
    j["arm_ok"] = (
        ((j["treatment"] == "Arm1") & (j["paper_arm"] == "Durvalumab")) |
        ((j["treatment"] == "Arm2") & (j["paper_arm"] == "SBRT + durva"))
    )
    for g in ["TACSTD2", "CLDN4"]:
        j[g] = j["Sample_title"].map(expr.loc[g])
    j["recurrence"] = np.where(
        j["status"].astype(str).str.contains("Recurrence", case=False, na=False),
        "recurrence", "no_recurrence")
    j["dfs_event"] = j["status"].astype(str).str.contains(
        "Recurrence|Dead", case=False, na=False).astype(int)
    j.to_csv(ANA / "GSE248378_joined_outcomes.csv", index=False)

    # Paper text: post-treatment RNAseq is 29 of 31 non-MPR tumours
    # (no MPR contrast exists in this matrix). Recurrence 9 vs 20.
    print(f"GSE248378 mapped {len(j)}; arm-concordant {int(j.arm_ok.sum())}; "
          f"discordant {j.loc[~j.arm_ok, 'Sample_title'].tolist()}; "
          f"mpr={j.mpr.value_counts(dropna=False).to_dict()}",
          file=sys.stderr)

    for g in ["TACSTD2", "CLDN4"]:
        sub2 = j.dropna(subset=[g, "recurrence"])
        by2 = {lab: sub2.loc[sub2.recurrence == lab, g].tolist()
               for lab in ["no_recurrence", "recurrence"]}
        if all(len(v) >= 2 for v in by2.values()):
            row = mw_row("GSE248378", g, "recurrence_postRx_FPKM_n29",
                         "no_recurrence", by2["no_recurrence"],
                         "recurrence", by2["recurrence"])
            row["note"] = "paper non-MPR post-Rx set; 45-M-PO kept (recurrence label OK)"
            boxplot(["no_recurrence", "recurrence"], by2,
                    f"GSE248378 NSCLC post-ICI non-MPR (n={len(sub2)})\n{g} vs recurrence",
                    f"GSE248378_{g}_recurrence.png", f"{g} (FPKM)")
        # sensitivity: drop the one arm-discordant sample
        sub3 = j[j["arm_ok"]].dropna(subset=[g, "recurrence"])
        by3 = {lab: sub3.loc[sub3.recurrence == lab, g].tolist()
               for lab in ["no_recurrence", "recurrence"]}
        if all(len(v) >= 2 for v in by3.values()):
            row = mw_row("GSE248378", g, "recurrence_postRx_FPKM_armConcordant",
                         "no_recurrence", by3["no_recurrence"],
                         "recurrence", by3["recurrence"])
            row["note"] = "excludes 45-M-PO (GEO Arm2 vs paper monotherapy)"
    surv = j.dropna(subset=["TACSTD2", "dfs_months"]).copy()
    med = float(surv["TACSTD2"].median())
    surv["tac_group"] = np.where(surv["TACSTD2"] >= med, "TACSTD2-high", "TACSTD2-low")
    km_plot(surv["dfs_months"].astype(float).values, surv["dfs_event"].values,
            surv["tac_group"].values,
            "GSE248378 TACSTD2 median-split DFS",
            "GSE248378_TACSTD2_DFS.png")
    return j


# ---------------------------------------------------------------- GSE193049
def analyze_gse193049():
    tar_path = DATA / "GSE193049_RAW.tar"
    meta = pd.read_csv(DATA / "GSE193049_sample_meta.csv")
    frames = []
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.name.endswith(".csv.gz"):
                continue
            gsm = m.name.split("_")[0]
            raw = tf.extractfile(m).read()
            df = pd.read_csv(gzip.GzipFile(fileobj=__import__("io").BytesIO(raw)))
            df = df.iloc[:, :2]
            df.columns = ["gene", gsm]
            frames.append(df.set_index("gene")[gsm])
    expr = pd.concat(frames, axis=1)
    expr.to_csv(ANA / "GSE193049_BALF_expr.csv")
    per = meta.copy()
    per["response"] = np.where(
        per["Sample_title"].str.contains("non_responder", case=False),
        "non_responder", "responder")
    for g in ["TACSTD2", "CLDN4"]:
        per[g] = per["gsm"].map(expr.loc[g])
    per.to_csv(ANA / "GSE193049_per_sample.csv", index=False)
    for g in ["TACSTD2", "CLDN4"]:
        by = {lab: per.loc[per.response == lab, g].tolist()
              for lab in ["responder", "non_responder"]}
        mw_row("GSE193049", g, "PD1_BALF_response",
               "responder", by["responder"],
               "non_responder", by["non_responder"])
        boxplot(["responder", "non_responder"], by,
                f"GSE193049 BALF (n={len(per)})\n{g} vs PD-1 response",
                f"GSE193049_{g}_response.png", f"{g} (author units)")
    print(f"GSE193049 n={len(per)} {per.response.value_counts().to_dict()}",
          file=sys.stderr)
    return per


def main():
    analyze_gse221733()
    analyze_gse248378()
    analyze_gse193049()
    res = pd.DataFrame(RESULTS)
    res["p_value"] = res["p_value"].astype(float).round(4)
    res.to_csv(ANA / "leftover_marker_outcome_tests.csv", index=False)
    print("\n=== leftover tests ===", file=sys.stderr)
    print(res.to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
