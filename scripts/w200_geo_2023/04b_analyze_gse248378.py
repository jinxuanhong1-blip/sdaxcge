#!/usr/bin/env python3
"""
GSE248378 leftover analysis.

GEO deposits post-durvalumab resected-tumour FPKMs (n=29) with treatment arm
and histology, but *no* recurrence / MPR label. The depositing paper
(Altorki et al., Nat Commun 2023, 10.1038/s41467-023-44195-x) states this
matrix is the 29 non-MPR tumours used in Fig. 5d / Fig. 6 (9 recurred, 20
did not). Figure 5d source data lists the 29 ITGAE (CD103) FPKM values by
recurrence group. Those 29 values match the GEO ITGAE column one-to-one at
2 decimal places, so recurrence can be recovered from *public* source data
without guessing sample IDs.

This is still a leftover-limited test:
  * POST-treatment, non-MPR only (MPR tumours were excluded from the matrix)
  * Twist-exome capture RNA-seq (~15k genes), not a full transcriptome
  * Recurrence, not radiographic ICI response
  * n=9 vs 20

Outputs under results/w200/GEO_2023/analysis/
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from geo_common import OUT, CACHE, fetch_series_matrix

ANA = OUT / "analysis"
ANA.mkdir(parents=True, exist_ok=True)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False

SRC = Path("/tmp/w200_cache/altorki/MOESM6.xlsx")
EXPR = CACHE / "GSE248378__GSE248378_Durva_Post_FPKMs.txt.gz"


def cliffs_delta(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    gt = sum(float(x > y) for x in a for y in b)
    lt = sum(float(x < y) for x in a for y in b)
    return (gt - lt) / (len(a) * len(b))


def parse_meta(text):
    gsm = None
    titles = None
    chars = []
    for line in text.splitlines():
        if line.startswith("!series_matrix_table_begin"):
            break
        if not line.startswith("!Sample_"):
            continue
        parts = line.rstrip("\n").split("\t")
        key = parts[0].lstrip("!")
        vals = [p.strip().strip('"') for p in parts[1:]]
        if key == "Sample_geo_accession":
            gsm = vals
        elif key == "Sample_title":
            titles = vals
        elif key.startswith("Sample_characteristics_ch"):
            chars.append(vals)
    df = pd.DataFrame({"gsm": gsm, "title": titles})
    parsed = {}
    for row in chars:
        if len(row) != len(gsm):
            continue
        for i, cell in enumerate(row):
            if ":" in cell:
                k, v = cell.split(":", 1)
                parsed.setdefault(k.strip().lower(), [None] * len(gsm))
                parsed[k.strip().lower()][i] = v.strip()
    for k, col in parsed.items():
        df[k] = col
    return df


def load_fig5d_recurrence():
    raw = pd.read_excel(SRC, sheet_name="Figure 5d", header=None)
    # row 2 is the group header; data starts at row 3
    no_rec = pd.to_numeric(raw.iloc[3:, 0], errors="coerce").dropna().astype(float)
    rec = pd.to_numeric(raw.iloc[3:, 1], errors="coerce").dropna().astype(float)
    return {round(float(v), 2): "no_recurrence" for v in no_rec} | {
        round(float(v), 2): "recurrence" for v in rec
    }


def boxplot(groups, data, title, fname, ylabel):
    if not HAVE_MPL:
        return
    fig, ax = plt.subplots(figsize=(4.4, 4))
    vals = [np.asarray(data[g], float) for g in groups]
    bp = ax.boxplot(vals, tick_labels=groups, patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], ["#4C72B0", "#C44E52"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, g in enumerate(groups, 1):
        y = np.asarray(data[g], float)
        ax.scatter(rng.normal(i, 0.05, len(y)), y, color="black", s=16, zorder=3, alpha=0.75)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(ANA / fname, dpi=140)
    plt.close(fig)


def test_row(cohort, gene, outcome, g1, a, g2, b, note=""):
    a, b = np.asarray(a, float), np.asarray(b, float)
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "cohort": cohort, "gene": gene, "outcome": outcome,
        "group1": g1, "n1": int(len(a)), "median1": round(float(np.median(a)), 3),
        "group2": g2, "n2": int(len(b)), "median2": round(float(np.median(b)), 3),
        "mannwhitney_U": float(u), "p_value": round(float(p), 4),
        "cliffs_delta": round(float(cliffs_delta(a, b)), 3),
        "note": note,
    }


def main():
    expr = pd.read_csv(EXPR, sep="\t", index_col=0)
    meta = parse_meta(fetch_series_matrix("GSE248378"))
    rec_map = load_fig5d_recurrence()

    itgae = expr.loc["ITGAE"].astype(float)
    keys = [round(float(v), 2) for v in itgae]
    if len(set(keys)) != len(keys):
        raise SystemExit("ITGAE FPKM not unique; cannot join Fig.5d")
    missing = [k for k in keys if k not in rec_map]
    if missing:
        raise SystemExit(f"ITGAE values not in Fig.5d: {missing}")

    per = meta.copy()
    per = per.set_index("title")
    per["ITGAE_FPKM"] = itgae.reindex(per.index)
    per["recurrence"] = [rec_map[round(float(v), 2)] for v in per["ITGAE_FPKM"]]
    for gene in ("TACSTD2", "CLDN4"):
        per[gene] = expr.loc[gene].reindex(per.index).astype(float)
    per["arm"] = per["treatment"].map({
        "Arm1": per["treatment"],  # keep raw; paper: Arm labels not defined in GEO
    })
    # Infer arm from title when possible (M = monotherapy, RT = dual) as a
    # *secondary* label only; GEO treatment Arm1/Arm2 is the deposited field.
    def title_arm(t):
        t = str(t)
        if "-M-" in t:
            return "monotherapy_from_title"
        if "-RT-" in t:
            return "dual_from_title"
        return ""
    per["title_arm_guess"] = [title_arm(t) for t in per.index]
    per.to_csv(ANA / "GSE248378_per_sample.csv")

    tests = []
    note = ("post-treatment non-MPR tumours only; recurrence recovered from "
            "Nat Commun source data Fig.5d via unique ITGAE FPKM")
    for gene in ("TACSTD2", "CLDN4"):
        a = per.loc[per["recurrence"] == "no_recurrence", gene]
        b = per.loc[per["recurrence"] == "recurrence", gene]
        tests.append(test_row(
            "GSE248378", gene, "recurrence_in_nonMPR",
            "no_recurrence", a, "recurrence", b, note))
        boxplot(["no_recurrence", "recurrence"],
                {"no_recurrence": a.values, "recurrence": b.values},
                f"GSE248378 {gene} vs recurrence (non-MPR, post-ICI)",
                f"GSE248378_{gene}_recurrence.png",
                f"{gene} FPKM")

        a2 = per.loc[per["treatment"] == "Arm1", gene]
        b2 = per.loc[per["treatment"] == "Arm2", gene]
        tests.append(test_row(
            "GSE248378", gene, "GEO_treatment_arm",
            "Arm1", a2, "Arm2", b2,
            "treatment assignment, NOT ICI response; Arm1/Arm2 as deposited"))

    # Spearman TACSTD2 vs CLDN4
    rho, p = stats.spearmanr(per["TACSTD2"], per["CLDN4"])
    tests.append({
        "cohort": "GSE248378", "gene": "TACSTD2_vs_CLDN4",
        "outcome": "spearman_within_cohort",
        "group1": "TACSTD2", "n1": len(per), "median1": round(float(per["TACSTD2"].median()), 3),
        "group2": "CLDN4", "n2": len(per), "median2": round(float(per["CLDN4"].median()), 3),
        "mannwhitney_U": "", "p_value": round(float(p), 4),
        "cliffs_delta": "",
        "note": f"Spearman rho={rho:.3f}",
    })

    new = pd.DataFrame(tests)
    dest = ANA / "marker_outcome_tests.csv"
    if dest.exists() and dest.stat().st_size > 10:
        try:
            old = pd.read_csv(dest)
        except pd.errors.EmptyDataError:
            old = pd.DataFrame()
        if not old.empty and "cohort" in old.columns:
            old = old[old["cohort"] != "GSE248378"]
            new = pd.concat([old, new], ignore_index=True)
    new.to_csv(dest, index=False)
    (ANA / "GSE248378_join_audit.json").write_text(json.dumps({
        "n_geo": int(len(per)),
        "n_fig5d": len(rec_map),
        "n_recurrence": int((per["recurrence"] == "recurrence").sum()),
        "n_no_recurrence": int((per["recurrence"] == "no_recurrence").sum()),
        "itgae_unique": True,
        "join": "round(ITGAE_FPKM, 2) == Fig.5d source-data value",
        "source": "41467_2023_44195_MOESM6_ESM.xlsx sheet Figure 5d",
        "doi": "10.1038/s41467-023-44195-x",
        "selection": "29/31 non-MPR resected tumours with post-treatment RNA-seq",
    }, indent=2))
    print(pd.DataFrame(tests).to_string(index=False))
    print(per["recurrence"].value_counts().to_string())


if __name__ == "__main__":
    main()
