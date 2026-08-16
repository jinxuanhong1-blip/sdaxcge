#!/usr/bin/env python3
"""Extract TACSTD2/CLDN4 from leftover matrices that actually measure them.

Honest scope
------------
GSE141479 (Clariom D, CD8+ PBMC, nivolumab pre/post) is the only 2020 leftover
with both genes on the platform AND an ICI treatment context. GEO does **not**
deposit a per-sample response / PFS / OS label, so this is an epithelial-null
compartment check (same role as GSE111414 in the 2019–2021 wave), not an
outcome test.

GSE154286 is a 201-gene targeted panel without TACSTD2/CLDN4 (documented).
GSE99995 is IFN-γ / PD-L1 stratified LUAD, not an ICI-treated leftover.
"""
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results" / "w200" / "GEO_2020"
DL = RES / "downloads"
CLIN = RES / "clinical"
TAB = RES / "tables"
FIG = RES / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

PROBES = {
    "TACSTD2": "TC0100014340.hg.1",
    "CLDN4": "TC0700007993.hg.1",
}


def read_matrix_rows(path, wanted):
    wanted = set(wanted)
    found = {}
    with gzip.open(path, "rt") as f:
        in_table = False
        header = None
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            cols = [c.strip().strip('"') for c in line.rstrip("\n").split("\t")]
            if header is None:
                header = cols
                continue
            if cols[0] in wanted:
                found[cols[0]] = [float(x) if x not in ("", "null", "NA") else np.nan
                                  for x in cols[1:]]
    return header[1:], found


def analyze_gse141479():
    acc = "GSE141479"
    samples, expr = read_matrix_rows(
        DL / acc / f"{acc}_series_matrix.txt.gz", PROBES.values()
    )
    clin = pd.read_csv(CLIN / f"{acc}_geo_characteristics.csv")
    clin = clin.set_index("geo_accession").loc[samples]
    df = pd.DataFrame({
        "geo_accession": samples,
        "individual": clin["individual"].values,
        "treatment": clin["treatment"].values,
        "cell_type": clin["cell_type"].values,
        "tissue": clin["tissue"].values,
        "diagnosis": clin["diagnosis"].values,
        "source": clin["source_name_ch1"].values,
    })
    for gene, probe in PROBES.items():
        df[gene] = expr[probe]
        df[f"{gene}_probe"] = probe
    df.to_csv(TAB / f"{acc}_expr.csv", index=False)

    rows = []
    for gene in PROBES:
        pre = df.loc[df.treatment == "Pre-treatment", gene]
        post = df.loc[df.treatment == "Post-treatment (nivolumab)", gene]
        u, p = stats.mannwhitneyu(pre, post, alternative="two-sided")
        auc = u / (len(pre) * len(post))
        rows.append({
            "dataset": acc,
            "compartment": "PBMC CD8+",
            "gene": gene,
            "probe": PROBES[gene],
            "outcome": "pre vs post nivolumab (NOT clinical response — GEO has no R/NR/PFS/OS)",
            "n": len(df),
            "n_group1": int(len(pre)),
            "n_group2": int(len(post)),
            "median_pre": round(float(pre.median()), 3),
            "median_post": round(float(post.median()), 3),
            "mean": round(float(df[gene].mean()), 3),
            "sd": round(float(df[gene].std()), 3),
            "auc_pre_gt_post": round(float(auc), 3),
            "p_value": round(float(p), 4),
            "direction": "higher pre" if pre.median() > post.median() else "higher post",
            "interpretable_as_ICI_outcome": "NO",
            "note": "epithelial gene in sorted CD8+ T cells; no per-sample ICI outcome in GEO",
        })
        # paired patients with both timepoints
        both = df.groupby("individual").filter(lambda x: set(x.treatment) >=
                                               {"Pre-treatment", "Post-treatment (nivolumab)"})
        pairs = []
        for ind, g in both.groupby("individual"):
            a = g.loc[g.treatment == "Pre-treatment", gene]
            b = g.loc[g.treatment == "Post-treatment (nivolumab)", gene]
            if len(a) and len(b):
                pairs.append((float(a.iloc[0]), float(b.iloc[0])))
        if pairs:
            a, b = zip(*pairs)
            w, pw = stats.wilcoxon(a, b)
            rows[-1]["n_paired"] = len(pairs)
            rows[-1]["paired_wilcoxon_p"] = round(float(pw), 4)

    # boxplot
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    order = ["Pre-treatment", "Post-treatment (nivolumab)"]
    for ax, gene in zip(axes, PROBES):
        data = [df.loc[df.treatment == lvl, gene].values for lvl in order]
        ax.boxplot(data, tick_labels=["pre", "post nivo"], showfliers=False)
        for i, lvl in enumerate(order, start=1):
            y = df.loc[df.treatment == lvl, gene].values
            x = np.random.normal(i, 0.06, size=len(y))
            ax.scatter(x, y, s=16, alpha=0.7, color="#377eb8")
        ax.set_title(f"{acc}: {gene} (CD8+ PBMC)")
        ax.set_ylabel("Clariom D normalized intensity")
    fig.tight_layout()
    fig.savefig(FIG / f"{acc}_prepost_boxplot.png", dpi=130)
    plt.close(fig)
    return df, rows


def document_gse154286_panel():
    """Write the 201 panel IDs so absence of TACSTD2/CLDN4 is auditable."""
    path = DL / "GSE154286" / "GSE154286_series_matrix.txt.gz"
    ids = []
    with gzip.open(path, "rt") as f:
        in_table = False
        header = False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            if not header:
                header = True
                continue
            ids.append(line.split("\t", 1)[0].strip().strip('"'))
    out = TAB / "GSE154286_panel_ids.txt"
    out.write_text("\n".join(ids) + "\n")
    blob = "\n".join(ids).upper()
    rec = {
        "accession": "GSE154286",
        "n_panel_ids": len(ids),
        "TACSTD2_present": "TACSTD2" in blob or "TROP2" in blob,
        "CLDN4_present": "CLDN4" in blob,
        "note": "SAKK19/09 NSCLC chemo (not ICI) targeted panel; genes absent",
    }
    (TAB / "GSE154286_panel_gene_check.json").write_text(json.dumps(rec, indent=2))
    return rec


def main():
    df, rows = analyze_gse141479()
    panel = document_gse154286_panel()
    res = pd.DataFrame(rows)
    res.to_csv(TAB / "combined_TACSTD2_CLDN4_leftover_2020.csv", index=False)
    print(res.to_string(index=False))
    print("\nGSE154286 panel:", panel)
    print("Wrote tables to", TAB)


if __name__ == "__main__":
    main()
