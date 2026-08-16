"""GSE248249 — the only public same-patient pre/post ICI *lung tumor* RNA we found.

Memon et al., Cancer Cell 2024. NSCLC treated with PD-(L)1 blockade.
  * 13 pre-treatment + 29 acquired-resistance (post) FFPE tumors (Clariom D).
  * 13 same-patient pairs (every pre has a matched post).
  * Post timepoint is **acquired resistance**, not an unselected on-treatment biopsy.

Probe mapping (clariomdhumantranscriptcluster.db / Entrez):
  TACSTD2 4070 -> TC0100014340.hg.1
  CLDN4   1364 -> TC0700007993.hg.1
  EPCAM   4072 -> TC0200007506.hg.1   (epithelial positive control)
  CD274  29126 -> TC0900006559.hg.1   (PD-L1, IO-axis control)
"""
from __future__ import annotations

import gzip
import re

import numpy as np
import pandas as pd

import config as C
from common import compare_groups, paired_test, savefig, PALETTE
import matplotlib.pyplot as plt

PROBES = {
    "TACSTD2": "TC0100014340.hg.1",
    "CLDN4": "TC0700007993.hg.1",
    "EPCAM": "TC0200007506.hg.1",
    "CD274": "TC0900006559.hg.1",
}


def _parse_meta(path) -> pd.DataFrame:
    titles, accs, chars = [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
    df = pd.DataFrame({"title": titles}, index=accs)
    for row in chars:
        key = None
        for cell in row:
            if ":" in cell:
                key = cell.split(":", 1)[0].strip()
                break
        if key is None:
            continue
        vals = [c.split(":", 1)[1].strip() if ":" in c else np.nan for c in row]
        col = key
        i = 1
        while col in df.columns:
            col = f"{key}_{i}"
            i += 1
        df[col] = vals
    df["patient"] = df["title"].str.extract(r"Patient (\d+)", expand=False)
    tp = df.get("timepoint", pd.Series(index=df.index, dtype=object)).astype(str).str.lower()
    df["timepoint"] = np.where(tp.str.contains("pre"), "pre",
                       np.where(tp.str.contains("post"), "post", "unknown"))
    return df


def _extract_probes(path, probes: dict[str, str]) -> pd.DataFrame:
    wanted = {v: k for k, v in probes.items()}
    rows = {}
    with gzip.open(path, "rt") as fh:
        header = None
        in_table = False
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if not in_table:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            if parts[0] == "ID_REF":
                header = parts[1:]
                continue
            if parts[0] in wanted:
                rows[wanted[parts[0]]] = [float(x) if x not in ("", "NA", "null") else np.nan
                                          for x in parts[1:]]
    if header is None:
        raise RuntimeError("no ID_REF header in series matrix")
    missing = [g for g in probes if g not in rows]
    if missing:
        raise RuntimeError(f"probes not found in matrix: {missing}")
    return pd.DataFrame(rows, index=header)


def main() -> None:
    matrix = C.raw_path("GSE248249_matrix")
    meta = _parse_meta(matrix)
    expr = _extract_probes(matrix, PROBES)
    long = expr.join(meta[["title", "patient", "timepoint", "tumor site", "gender"]],
                     how="inner")
    long.to_csv(C.TABLES_DIR / "gse248249_paired_long.csv")
    print(long.groupby("timepoint").size())
    print("patients with both timepoints:",
          sorted(set(long.loc[long.timepoint == "pre", "patient"])
                 & set(long.loc[long.timepoint == "post", "patient"])))

    rows = []
    # unpaired pre vs post (all samples)
    for gene in PROBES:
        pre = long.loc[long.timepoint == "pre", gene].values
        post = long.loc[long.timepoint == "post", gene].values
        rows.append(compare_groups(
            "GSE248249", gene, "unpaired post_vs_pre (acquired resistance)",
            post, pre, "post(AR)", "pre").as_row())

    # paired Wilcoxon
    pre_df = long[long.timepoint == "pre"].set_index("patient")
    post_df = long[long.timepoint == "post"].set_index("patient")
    common = pre_df.index.intersection(post_df.index)
    # if a patient has >1 post, take the first (titles are unique; keep first)
    post_df = post_df[~post_df.index.duplicated(keep="first")]
    pre_df = pre_df[~pre_df.index.duplicated(keep="first")]
    common = pre_df.index.intersection(post_df.index)

    pair_rows = []
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, gene in zip(axes, C.TARGET_GENES):
        pre_v = pre_df.loc[common, gene].astype(float)
        post_v = post_df.loc[common, gene].astype(float)
        stat, p, med_delta = paired_test(pre_v.values, post_v.values)
        n_up = int((post_v.values - pre_v.values > 0).sum())
        n_dn = int((post_v.values - pre_v.values < 0).sum())
        rows.append({
            "dataset": "GSE248249", "gene": gene,
            "comparison": "paired post_vs_pre (acquired resistance)",
            "n_hi": int(len(common)), "n_lo": int(len(common)),
            "log2fc_hi_vs_lo": float(np.mean(post_v - pre_v)),
            "auc": np.nan, "test": "Wilcoxon signed-rank",
            "statistic": stat, "p_value": p,
            "n_pairs": int(len(common)),
            "n_pairs_up": n_up, "n_pairs_down": n_dn,
            "median_delta_post_minus_pre": med_delta,
        })
        for pt in common:
            ax.plot([0, 1], [pre_v.loc[pt], post_v.loc[pt]], "-o",
                    color="#3182bd", alpha=0.75, markersize=5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Pre", "Post (acquired\nresistance)"])
        ax.set_ylabel(f"{gene} (SST-RMA, log2)")
        ax.set_title(f"GSE248249 {gene}: paired pre→AR\n"
                     f"median Δ={med_delta:+.3f}, {n_up}↑/{n_dn}↓, p={p:.3f}, n={len(common)}",
                     fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
        pair_rows.append(pd.DataFrame({
            "patient": list(common), "gene": gene,
            "pre": pre_v.loc[common].values, "post": post_v.loc[common].values,
            "delta": (post_v - pre_v).loc[common].values,
        }))

    fig.suptitle("GSE248249 NSCLC PD-(L)1 blockade — same-patient pre vs acquired-resistance RNA",
                 fontsize=11)
    savefig(fig, C.FIGURES_DIR / "fig5_gse248249_paired_pre_post.png")

    pd.concat(pair_rows, ignore_index=True).to_csv(
        C.TABLES_DIR / "gse248249_patient_pairs.csv", index=False)

    # also report EPCAM / CD274 paired as controls
    for gene in ["EPCAM", "CD274"]:
        pre_v = pre_df.loc[common, gene].astype(float)
        post_v = post_df.loc[common, gene].astype(float)
        stat, p, med_delta = paired_test(pre_v.values, post_v.values)
        n_up = int((post_v.values - pre_v.values > 0).sum())
        n_dn = int((post_v.values - pre_v.values < 0).sum())
        rows.append({
            "dataset": "GSE248249", "gene": gene,
            "comparison": "paired post_vs_pre (acquired resistance) [control]",
            "n_hi": int(len(common)), "n_lo": int(len(common)),
            "log2fc_hi_vs_lo": float(np.mean(post_v - pre_v)),
            "auc": np.nan, "test": "Wilcoxon signed-rank",
            "statistic": stat, "p_value": p,
            "n_pairs": int(len(common)),
            "n_pairs_up": n_up, "n_pairs_down": n_dn,
            "median_delta_post_minus_pre": med_delta,
        })

    out = pd.DataFrame(rows)
    out.to_csv(C.TABLES_DIR / "gse248249_paired_stats.csv", index=False)
    print(out.to_string())


if __name__ == "__main__":
    main()
