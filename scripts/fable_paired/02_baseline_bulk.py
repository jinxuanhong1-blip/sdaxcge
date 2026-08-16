"""Baseline (pre-treatment) bulk RNA-seq: TACSTD2/CLDN4 vs ICI response.

Three independent NSCLC cohorts profiled before anti-PD-(L)1 therapy:
  * GSE207422 bulk  - 24 pre-treatment biopsies, MPR vs NMPR pathologic response.
  * GSE126044       - 16 pre-treatment tumors, responder vs non-responder.
  * GSE135222       - 27 pre-treatment tumors, durable benefit proxied by PFS>=6 mo.

Each cohort is analysed independently; consistency across them is the verification.
"""
from __future__ import annotations

import gzip
import re

import numpy as np
import pandas as pd

import config as C
from common import compare_groups, strip_box, savefig, PALETTE
import matplotlib.pyplot as plt


def parse_series_matrix(path) -> pd.DataFrame:
    """Parse a GEO series_matrix into a DataFrame indexed by GSM accession."""
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
        # each characteristic row is "key: value"; use key from first non-empty cell
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
            col = f"{key}_{i}"; i += 1
        df[col] = vals
    return df


def _target_rows_symbol(path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df.loc[[g for g in C.TARGET_GENES if g in df.index]]


# ---------------------------------------------------------------------------
def cohort_gse207422_bulk() -> list[dict]:
    expr = _target_rows_symbol(C.raw_path("GSE207422_bulk_expr"))  # log2TPM
    meta = pd.read_excel(C.raw_path("GSE207422_bulk_meta")).dropna(subset=["Sample"])
    pr = meta["Pathologic Response"].astype(str).str.upper()
    meta["response"] = np.where(pr.str.contains("NMPR"), "non_responder",
                        np.where(pr.str.contains("MPR") | pr.str.contains("PCR"),
                                 "responder", "unknown"))
    meta = meta.set_index("Sample")
    rows, plot = [], {}
    for gene in C.TARGET_GENES:
        vals = expr.loc[gene]
        r = [vals[s] for s in meta.index if s in vals and meta.loc[s, "response"] == "responder"]
        nr = [vals[s] for s in meta.index if s in vals and meta.loc[s, "response"] == "non_responder"]
        rows.append(compare_groups("GSE207422_bulk", gene,
                                    "baseline: responder_vs_nonresponder",
                                    r, nr, "responder(MPR)", "non_responder(NMPR)").as_row())
        plot[gene] = ({"responder": np.array(r), "non_responder": np.array(nr)},
                      f"{gene} (log2 TPM)")
    return rows, plot, "GSE207422 bulk (baseline)"


def cohort_gse126044() -> list[dict]:
    counts = pd.read_csv(C.raw_path("GSE126044_counts"), sep="\t", index_col=0)
    # counts -> log2 CPM
    cpm = counts / counts.sum(axis=0) * 1e6
    logcpm = np.log2(cpm + 1)
    sm = parse_series_matrix(C.raw_path("GSE126044_matrix"))
    resp_col = [c for c in sm.columns if "response" in c][0]
    sm["dis"] = sm["title"].str.replace("RNA-seq_", "", regex=False)
    sm["response"] = np.where(sm[resp_col].str.contains("non", case=False), "non_responder",
                              "responder")
    dis2resp = dict(zip(sm["dis"], sm["response"]))
    rows, plot = [], {}
    for gene in C.TARGET_GENES:
        if gene not in logcpm.index:
            continue
        vals = logcpm.loc[gene]
        r = [vals[s] for s in vals.index if dis2resp.get(s) == "responder"]
        nr = [vals[s] for s in vals.index if dis2resp.get(s) == "non_responder"]
        rows.append(compare_groups("GSE126044", gene,
                                    "baseline: responder_vs_nonresponder",
                                    r, nr, "responder", "non_responder").as_row())
        plot[gene] = ({"responder": np.array(r), "non_responder": np.array(nr)},
                      f"{gene} (log2 CPM)")
    return rows, plot, "GSE126044 (baseline)"


def cohort_gse135222(pfs_months_cut: float = 6.0) -> list[dict]:
    tpm = pd.read_csv(C.raw_path("GSE135222_tpm"), sep="\t", index_col=0)
    tpm.index = tpm.index.str.split(".").str[0]  # strip ENSG version
    logtpm = np.log2(tpm + 1)
    sm = parse_series_matrix(C.raw_path("GSE135222_matrix"))
    tcol = [c for c in sm.columns if "pfs.time" in c.lower() or "pfs time" in c.lower()]
    tcol = tcol[0] if tcol else [c for c in sm.columns if "time" in c.lower()][0]
    sm["pfs_time"] = pd.to_numeric(sm[tcol], errors="coerce")
    sm["nsclc"] = sm["title"].str.replace(" ", "", regex=False)  # "NSCLC 990"->"NSCLC990"
    sm["response"] = np.where(sm["pfs_time"] >= pfs_months_cut * 30.0,
                              "responder", "non_responder")
    name2resp = dict(zip(sm["nsclc"], sm["response"]))
    rows, plot = [], {}
    for gene in C.TARGET_GENES:
        ens = C.ENSEMBL_IDS[gene]
        if ens not in logtpm.index:
            continue
        vals = logtpm.loc[ens]
        r = [vals[s] for s in vals.index if name2resp.get(s) == "responder"]
        nr = [vals[s] for s in vals.index if name2resp.get(s) == "non_responder"]
        rows.append(compare_groups("GSE135222", gene,
                                    f"baseline: DCB(PFS>={pfs_months_cut}mo)_vs_NDB",
                                    r, nr, "responder(DCB)", "non_responder(NDB)").as_row())
        plot[gene] = ({"responder": np.array(r), "non_responder": np.array(nr)},
                      f"{gene} (log2 TPM)")
    return rows, plot, "GSE135222 (baseline)"


def main() -> None:
    all_rows = []
    cohorts = [cohort_gse207422_bulk(), cohort_gse126044(), cohort_gse135222()]
    for rows, _, _ in cohorts:
        all_rows.extend(rows)
    stats_df = pd.DataFrame(all_rows)
    stats_df.to_csv(C.TABLES_DIR / "baseline_bulk_stats.csv", index=False)
    print(stats_df.to_string())

    # figure: 3 cohorts x 2 genes
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for j, (_, plot, title) in enumerate(cohorts):
        for i, gene in enumerate(C.TARGET_GENES):
            groups, ylabel = plot[gene]
            named = {"responder": groups["responder"], "non-responder": groups["non_responder"]}
            strip_box(axes[i, j], named,
                      {"responder": PALETTE["responder"], "non-responder": PALETTE["non_responder"]},
                      ylabel, f"{title}\n{gene}")
    fig.suptitle("Baseline TACSTD2 / CLDN4 by ICI response (three NSCLC cohorts)", fontsize=12)
    savefig(fig, C.FIGURES_DIR / "fig2_baseline_bulk_cohorts.png")


if __name__ == "__main__":
    main()
