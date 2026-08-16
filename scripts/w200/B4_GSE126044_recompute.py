#!/usr/bin/env python3
"""B4 / w200 exact recompute: GSE126044 CLDN4 and tight-junction (TJ) score, responder vs non-responder.

Context
-------
A prior batch summary claimed for GSE126044: "NR higher TJ, p=0.019".
This script honestly recomputes, from the raw count matrix, the responder (R) vs
non-responder (NR) contrast for:

  * CLDN4 alone (log2 CPM)
  * A tight-junction (TJ) signature score, computed several defensible ways so the
    reader can see how sensitive the p-value is to the gene set and scoring method.

Dataset: GSE126044 (Cho et al.), pre-treatment NSCLC tumor biopsies, anti-PD-1.
Raw counts = integer gene x sample matrix, gene symbols as row index.
n = 16 patients (5 responders, 11 non-responders).

Normalisation: log2 CPM = log2(1e6 * count / library_size + 1)  (same as the
existing fable_geo pipeline in this repo, so numbers are comparable).

Signature score (default): per-gene z-score across the 16 samples (log2 CPM),
then the unweighted mean of z-scores over genes that are present and have non-zero
variance. This is the standard "mean-z" signature score. A "mean log2 CPM"
variant (no z-scoring) is reported alongside as a sensitivity check.

Statistics: two-sided Mann-Whitney U (Wilcoxon rank-sum) is the primary test
(small n, non-parametric); Welch's t-test is reported as a secondary check.
Effect sizes: rank-biserial correlation / AUC (from U) and Cohen's d.

Nothing here is tuned to hit any particular p-value. All gene sets, all scoring
methods, and both tests are reported regardless of outcome.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
OUT = ROOT / "results" / "w200" / "B4_GSE126044"
DATA = OUT / "data"
TAB = OUT / "tables"
FIG = OUT / "figures"
for d in (TAB, FIG):
    d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Gene sets
# --------------------------------------------------------------------------
# Curated canonical / structural tight-junction genes (claudins, occludin,
# tricellulin, ZO proteins, JAMs, cingulin, and core apical polarity scaffolds
# that build the junction). This is the "core TJ" set used as the primary TJ score.
CORE_TJ = [
    # claudins
    "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "CLDN6", "CLDN7", "CLDN8",
    "CLDN9", "CLDN10", "CLDN11", "CLDN12", "CLDN14", "CLDN15", "CLDN16",
    "CLDN17", "CLDN18", "CLDN19", "CLDN20", "CLDN23",
    # occludin / MARVEL
    "OCLN", "MARVELD2", "MARVELD3",
    # ZO / scaffolds
    "TJP1", "TJP2", "TJP3", "SYMPK", "CGN", "CGNL1",
    # junctional adhesion molecules
    "F11R", "JAM2", "JAM3",
    # apical polarity complexes that anchor the TJ
    "PARD3", "PARD6A", "PARD6B", "PRKCZ", "PRKCI",
    "INADL", "PATJ", "MPP5", "CRB3", "AMOT",
]

# Repo-native TJ modules from hunt-tj-gsea (MSigDB + custom literature).
# Included as a priori sensitivity sets, not as post-hoc p-hunting.
CUSTOM_TJ_CORE = [
    "AMOT", "AMOTL1", "CGN", "CGNL1", "CLDN1", "CLDN3", "CLDN4", "CLDN7",
    "CLDN8", "CRB3", "F11R", "INADL", "JAM2", "JAM3", "MAGI1", "MAGI3",
    "MARVELD2", "MPDZ", "OCLN", "PARD3", "PARD6A", "PARD6B", "PATJ",
    "TJP1", "TJP2", "TJP3",
]
CUSTOM_CLAUDIN_PAR_FOCAL = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]


def load_kegg_tj(path: Path) -> list[str]:
    """Parse the GENE section of a KEGG `get/hsa04530` flat file into symbols."""
    if not path.exists():
        return []
    symbols: list[str] = []
    in_gene = False
    for line in path.read_text().splitlines():
        if line.startswith("GENE"):
            in_gene = True
            line = line[len("GENE"):]
        elif line and not line.startswith(" "):
            # a new top-level section began
            in_gene = False
        if not in_gene:
            continue
        m = re.search(r"\d+\s+([A-Za-z0-9\-]+);", line)
        if m:
            symbols.append(m.group(1))
    # de-duplicate, preserve order
    seen: set[str] = set()
    out = []
    for s in symbols:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------
def log2cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0)
    cpm = counts.divide(lib, axis=1) * 1e6
    return np.log2(cpm + 1)


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = pd.read_csv(DATA / "GSE126044_counts.txt.gz", sep="\t", index_col=0)
    clin = pd.read_csv(DATA / "GSE126044_clinical.csv")
    rows = []
    for _, r in clin.iterrows():
        key = r["title"].replace("RNA-seq_", "")
        rows.append({
            "sample": key,
            "gsm": r["gsm"],
            "response": r["patient_response"],
            "sample_type": r["sample"],
        })
    meta = pd.DataFrame(rows).set_index("sample")
    meta = meta.loc[counts.columns]
    return counts, meta


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------
def signature_score(logcpm: pd.DataFrame, genes: list[str], method: str):
    """Return (per-sample score Series, list of genes actually used).

    method="meanz"    -> mean of per-gene z-scores (across samples)
    method="meanexpr" -> mean log2 CPM (no z-scoring)
    """
    present = [g for g in genes if g in logcpm.index]
    sub = logcpm.loc[present]
    # drop zero-variance genes (uninformative; z-score undefined)
    var = sub.var(axis=1)
    used = var[var > 0].index.tolist()
    sub = sub.loc[used]
    if method == "meanz":
        z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1), axis=0)
        score = z.mean(axis=0)
    elif method == "meanexpr":
        score = sub.mean(axis=0)
    else:
        raise ValueError(method)
    return score, used


def compare(score: pd.Series, response: pd.Series) -> dict:
    idx_r = response[response == "responder"].index
    idx_nr = response[response == "non-responder"].index
    r = score.loc[idx_r].astype(float)
    nr = score.loc[idx_nr].astype(float)
    u, p_mw = stats.mannwhitneyu(r, nr, alternative="two-sided")
    n1, n2 = len(r), len(nr)
    auc_r = u / (n1 * n2)  # P(responder > non-responder)
    rank_biserial = 2 * auc_r - 1
    t, p_t = stats.ttest_ind(r, nr, equal_var=False)
    # Cohen's d (pooled)
    sp = np.sqrt(((n1 - 1) * r.var(ddof=1) + (n2 - 1) * nr.var(ddof=1)) / (n1 + n2 - 2))
    d = (r.mean() - nr.mean()) / sp if sp > 0 else np.nan
    higher = "responders" if r.median() > nr.median() else "non-responders"
    return {
        "n_responder": n1,
        "n_nonresponder": n2,
        "median_R": round(float(r.median()), 4),
        "median_NR": round(float(nr.median()), 4),
        "mean_R": round(float(r.mean()), 4),
        "mean_NR": round(float(nr.mean()), 4),
        "direction": f"higher in {higher}",
        "AUC_R_gt_NR": round(float(auc_r), 4),
        "rank_biserial": round(float(rank_biserial), 4),
        "cohens_d_R_minus_NR": round(float(d), 4),
        "mannwhitney_U": float(u),
        "p_mannwhitney": round(float(p_mw), 4),
        "welch_t": round(float(t), 4),
        "p_welch_t": round(float(p_t), 4),
    }


def boxplot(score: pd.Series, response: pd.Series, title: str, ylab: str, fname: Path):
    order = ["responder", "non-responder"]
    data = [score[response[response == lvl].index].values for lvl in order]
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    ax.boxplot(data, tick_labels=["R", "NR"], showfliers=False, widths=0.6)
    colors = ["#2c7fb8", "#d95f0e"]
    for i, (vals, c) in enumerate(zip(data, colors), start=1):
        x = np.random.normal(i, 0.06, size=len(vals))
        ax.scatter(x, vals, alpha=0.8, s=26, color=c, edgecolor="k", linewidth=0.3)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylab)
    fig.tight_layout()
    fig.savefig(fname, dpi=140)
    plt.close(fig)


# --------------------------------------------------------------------------
def main() -> None:
    np.random.seed(0)
    counts, meta = load()
    response = meta["response"]
    logcpm = log2cpm(counts)

    kegg_tj = load_kegg_tj(DATA / "kegg_hsa04530.txt")

    # persist the resolved gene sets for reproducibility
    (DATA / "geneset_core_tj.txt").write_text("\n".join(CORE_TJ) + "\n")
    (DATA / "geneset_custom_tj_core.txt").write_text("\n".join(CUSTOM_TJ_CORE) + "\n")
    (DATA / "geneset_custom_claudin_par_focal.txt").write_text(
        "\n".join(CUSTOM_CLAUDIN_PAR_FOCAL) + "\n"
    )
    if kegg_tj:
        (DATA / "geneset_kegg_hsa04530.txt").write_text("\n".join(kegg_tj) + "\n")

    # ---- per-sample table ----
    persample = meta.reset_index()
    if "sample" not in persample.columns:
        persample = persample.rename(columns={persample.columns[0]: "sample"})
    persample["CLDN4_log2cpm"] = logcpm.loc["CLDN4", persample["sample"]].values
    persample["TACSTD2_log2cpm"] = (
        logcpm.loc["TACSTD2", persample["sample"]].values
        if "TACSTD2" in logcpm.index else np.nan
    )

    score_defs = [
        ("core_TJ", CORE_TJ, "meanz"),
        ("core_TJ_no_CLDN4", [g for g in CORE_TJ if g != "CLDN4"], "meanz"),
        ("core_TJ_meanexpr", CORE_TJ, "meanexpr"),
        ("CUSTOM_TJ_CORE", CUSTOM_TJ_CORE, "meanz"),
        ("CUSTOM_CLAUDIN_PAR_FOCAL", CUSTOM_CLAUDIN_PAR_FOCAL, "meanz"),
    ]
    if kegg_tj:
        score_defs += [
            ("KEGG_hsa04530", kegg_tj, "meanz"),
            ("KEGG_hsa04530_meanexpr", kegg_tj, "meanexpr"),
        ]

    used_map = {}
    for name, genes, method in score_defs:
        s, used = signature_score(logcpm, genes, method)
        persample[name] = s.loc[persample["sample"]].values
        used_map[name] = used

    persample.to_csv(TAB / "per_sample_scores.csv", index=False)

    # ---- statistics (all samples; primary) ----
    rows = []

    def add_row(feature, typ, n_used, score, resp, subset):
        row = {
            "subset": subset,
            "feature": feature,
            "type": typ,
            "n_genes_used": n_used,
        }
        row.update(compare(score, resp))
        rows.append(row)

    cldn4 = logcpm.loc["CLDN4"]
    add_row("CLDN4 (log2 CPM)", "single_gene", 1, cldn4, response, "all")
    if "TACSTD2" in logcpm.index:
        add_row("TACSTD2 (log2 CPM)", "single_gene", 1,
                logcpm.loc["TACSTD2"], response, "all")
    for name, genes, method in score_defs:
        s, used = signature_score(logcpm, genes, method)
        add_row(name, f"signature_{method}", len(used), s, response, "all")

    # Fresh-only sensitivity: all 5 FFPE samples are NR, so sample type is
    # fully confounded with response. Restricting to fresh biopsies is a
    # legitimate robustness check, not a p-hunt.
    fresh_idx = meta.index[meta["sample_type"] == "fresh"]
    resp_fresh = response.loc[fresh_idx]
    add_row("CLDN4 (log2 CPM)", "single_gene", 1,
            cldn4.loc[fresh_idx], resp_fresh, "fresh_only")
    core_z, used_core = signature_score(logcpm, CORE_TJ, "meanz")
    add_row("core_TJ", "signature_meanz", len(used_core),
            core_z.loc[fresh_idx], resp_fresh, "fresh_only")
    custom_z, used_custom = signature_score(logcpm, CUSTOM_TJ_CORE, "meanz")
    add_row("CUSTOM_TJ_CORE", "signature_meanz", len(used_custom),
            custom_z.loc[fresh_idx], resp_fresh, "fresh_only")

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(TAB / "R_vs_NR_stats.csv", index=False)

    # gene-set membership / detection table
    memb_rows = []
    for name, genes, method in score_defs:
        for g in genes:
            memb_rows.append({
                "geneset": name,
                "gene": g,
                "in_matrix": g in logcpm.index,
                "used_in_score": g in used_map[name],
            })
    pd.DataFrame(memb_rows).to_csv(TAB / "geneset_membership.csv", index=False)

    crosstab = pd.crosstab(meta["response"], meta["sample_type"], margins=True)
    crosstab.to_csv(TAB / "response_by_sample_type.csv")

    # ---- figures ----
    boxplot(cldn4, response, "GSE126044 CLDN4 (R vs NR)", "log2 CPM",
            FIG / "CLDN4_R_vs_NR.png")
    boxplot(core_z, response, "GSE126044 core TJ score (mean-z, R vs NR)",
            "TJ score (mean z)", FIG / "coreTJ_R_vs_NR.png")
    boxplot(custom_z, response, "GSE126044 CUSTOM_TJ_CORE (mean-z, R vs NR)",
            "CUSTOM_TJ_CORE (mean z)", FIG / "customTJ_R_vs_NR.png")
    if kegg_tj:
        kegg_z, _ = signature_score(logcpm, kegg_tj, "meanz")
        boxplot(kegg_z, response, "GSE126044 KEGG TJ score (mean-z, R vs NR)",
                "KEGG TJ score (mean z)", FIG / "keggTJ_R_vs_NR.png")

    # ---- console report ----
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 40)
    print("GSE126044 recompute (n=%d: %d R, %d NR)" % (
        len(response), (response == "responder").sum(),
        (response == "non-responder").sum()))
    print("\nR vs NR statistics:")
    print(stats_df.to_string(index=False))

    summary = {
        "dataset": "GSE126044",
        "n_total": int(len(response)),
        "n_responder": int((response == "responder").sum()),
        "n_nonresponder": int((response == "non-responder").sum()),
        "normalisation": "log2 CPM = log2(1e6*count/lib + 1)",
        "primary_test": "two-sided Mann-Whitney U (Wilcoxon rank-sum)",
        "results": rows,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print("\nWrote:")
    for p in [TAB / "R_vs_NR_stats.csv", TAB / "per_sample_scores.csv",
              TAB / "geneset_membership.csv", OUT / "summary.json"]:
        print("  ", p.relative_to(ROOT))


if __name__ == "__main__":
    main()
