#!/usr/bin/env python3
"""
Step 6: extract TACSTD2 / CLDN4 (TROP2 = TACSTD2 protein) from the series that
carry an embedded expression table, map probes -> gene symbols via the platform
annotation, group samples heuristically, and compute descriptive + comparison
statistics.

Writes:
  results/geo_sweep/marker_values_embedded.tsv   (per gene, per sample)
  results/geo_sweep/marker_stats_embedded.tsv    (per gene, per group)
  notes/geo_sweep/analysis_embedded.json         (per-series detail + notes)
"""
import json
import re
import sys
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_annot import build_symbol_map  # noqa: E402

NOTES = Path("notes/geo_sweep")
RESULTS = Path("results/geo_sweep")
MATRIX_DIR = RESULTS / "matrices"

TACSTD2_SET = {"TACSTD2", "TROP2", "TROP-2", "GA733-1", "M1S1", "EGP-1", "EGP1"}
CLDN4_SET = {"CLDN4", "CLAUDIN-4", "CLAUDIN4", "CPETR1", "CPE-R", "CPER"}
GENE_SETS = {"TACSTD2": TACSTD2_SET, "CLDN4": CLDN4_SET}


def load_expr(acc):
    """Return (expr_df indexed by probe id, sample_titles list)."""
    path = MATRIX_DIR / f"{acc}_series_matrix.txt.gz"
    sample_titles, sample_geo = [], []
    with gzip.open(path, "rt", errors="replace") as fh:
        lines = fh.readlines()
    for ln in lines:
        if ln.startswith("!Sample_title"):
            sample_titles = [c.strip().strip('"') for c in ln.split("\t")[1:]]
        elif ln.startswith("!Sample_geo_accession"):
            sample_geo = [c.strip().strip('"') for c in ln.split("\t")[1:]]
    # locate table
    start = end = None
    for i, ln in enumerate(lines):
        if ln.startswith("!series_matrix_table_begin"):
            start = i + 1
        elif ln.startswith("!series_matrix_table_end"):
            end = i
            break
    if start is None:
        return None, sample_titles, sample_geo
    from io import StringIO
    tbl = "".join(lines[start:end])
    df = pd.read_csv(StringIO(tbl), sep="\t", index_col=0)
    df.index = df.index.astype(str).str.strip('"')
    df.columns = [c.strip().strip('"') for c in df.columns]
    df = df.apply(pd.to_numeric, errors="coerce")
    return df, sample_titles, sample_geo


def group_label(title):
    t = title.strip()
    # cut verbose ", biological rep1, Exp..." style replicate suffixes
    t = re.sub(r"[,;]?\s*(biological\s+)?rep(licate)?\.?\s*\d+.*$", "", t,
               flags=re.IGNORECASE)
    # cut trailing "_1", " 1", "-1" replicate indices
    t = re.sub(r"[ _\-.]*(rep(licate)?)?[ _\-.]*\d+$", "", t, flags=re.IGNORECASE)
    t = t.strip(" ,;_-")
    return t if t else title


def guess_scale(df):
    mx = np.nanmax(df.values)
    return "log2-like" if mx < 30 else "linear"


def analyze(acc, gpl):
    df, titles, geos = load_expr(acc)
    detail = {"accession": acc, "platform": gpl, "genes": {}}
    if df is None or df.empty:
        detail["status"] = "no_expression_table"
        return detail, [], []
    smap = build_symbol_map(gpl) if gpl else {}
    # If row IDs already look like symbols, use identity mapping too.
    ids_upper = {i: i.upper() for i in df.index}
    scale = guess_scale(df)
    detail["scale_guess"] = scale
    detail["n_probes_mapped"] = len(smap)

    if not titles or len(titles) != df.shape[1]:
        titles = list(df.columns)
    groups = [group_label(t) for t in titles]
    detail["sample_titles"] = titles
    detail["group_labels"] = groups

    values_rows, stats_rows = [], []
    for gene, aliases in GENE_SETS.items():
        # probes whose mapped symbol (annotation OR identity) is in aliases
        probes = [p for p in df.index
                  if smap.get(p, "").upper() in aliases
                  or ids_upper[p] in aliases]
        gd = {"n_probes": len(probes), "probes": probes[:20]}
        if not probes:
            gd["status"] = "gene_not_on_platform"
            detail["genes"][gene] = gd
            continue
        sub = df.loc[probes]
        collapsed = sub.mean(axis=0)  # mean across probes per sample
        gd["status"] = "ok"
        gd["per_sample"] = {}
        for col, val in collapsed.items():
            gd["per_sample"][col] = None if pd.isna(val) else float(val)
            values_rows.append({
                "accession": acc, "gene": gene, "sample": col,
                "group": groups[list(df.columns).index(col)] if col in list(df.columns) else "",
                "value": None if pd.isna(val) else float(val),
                "scale": scale,
            })
        # per-group descriptive stats
        gdf = pd.DataFrame({"value": collapsed.values, "group": groups})
        grp_stats = {}
        for gname, sub2 in gdf.groupby("group"):
            vals = sub2["value"].dropna().values
            grp_stats[gname] = {
                "n": int(len(vals)),
                "mean": float(np.mean(vals)) if len(vals) else None,
                "median": float(np.median(vals)) if len(vals) else None,
                "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else None,
            }
            stats_rows.append({
                "accession": acc, "gene": gene, "group": gname,
                "n": len(vals),
                "mean": round(float(np.mean(vals)), 4) if len(vals) else "",
                "median": round(float(np.median(vals)), 4) if len(vals) else "",
                "sd": round(float(np.std(vals, ddof=1)), 4) if len(vals) > 1 else "",
                "scale": scale,
            })
        gd["group_stats"] = grp_stats
        # pairwise comparison if exactly 2 groups
        gvals = [s["value"].dropna().values for _, s in gdf.groupby("group")]
        gnames = [n for n, _ in gdf.groupby("group")]
        if len(gvals) == 2 and all(len(v) >= 2 for v in gvals):
            try:
                u, p = stats.mannwhitneyu(gvals[0], gvals[1], alternative="two-sided")
                m0, m1 = np.mean(gvals[0]), np.mean(gvals[1])
                diff = m1 - m0
                log2fc = diff if scale == "log2-like" else (
                    float(np.log2(m1 / m0)) if m0 > 0 and m1 > 0 else None)
                gd["comparison"] = {
                    "group_a": gnames[0], "group_b": gnames[1],
                    "mean_a": float(m0), "mean_b": float(m1),
                    "log2FC_b_vs_a": None if log2fc is None else round(float(log2fc), 4),
                    "mannwhitney_p": round(float(p), 5),
                }
            except Exception as e:  # noqa: BLE001
                gd["comparison"] = {"error": str(e)}
        detail["genes"][gene] = gd
    detail["status"] = "analyzed"
    return detail, values_rows, stats_rows


# (accession, GPL) for the embedded-expression series (miRNA/panel-only noted)
TARGETS = [
    ("GSE108417", "GPL6246"),
    ("GSE33348", "GPL6246"),
    ("GSE48443", "GPL6885"),
    ("GSE206613", "GPL32069"),
    ("GSE210340", "GPL20775"),
    ("GSE266384", "GPL11202"),
    ("GSE310370", "GPL21572"),   # miRNA array -> no mRNA markers expected
    ("GSE317352", "GPL32057"),
]


def main():
    all_details, all_values, all_stats = [], [], []
    for acc, gpl in TARGETS:
        print(f"Analyzing {acc} ({gpl}) ...", flush=True)
        try:
            detail, vrows, srows = analyze(acc, gpl)
        except Exception as e:  # noqa: BLE001
            detail, vrows, srows = {"accession": acc, "platform": gpl,
                                    "status": f"error: {e}"}, [], []
        all_details.append(detail)
        all_values.extend(vrows)
        all_stats.extend(srows)
        for gene in ("TACSTD2", "CLDN4"):
            g = detail.get("genes", {}).get(gene, {})
            print(f"    {gene}: {g.get('status')} probes={g.get('n_probes')}")

    (NOTES / "analysis_embedded.json").write_text(json.dumps(all_details, indent=2))
    if all_values:
        pd.DataFrame(all_values).to_csv(
            RESULTS / "marker_values_embedded.tsv", sep="\t", index=False)
    if all_stats:
        pd.DataFrame(all_stats).to_csv(
            RESULTS / "marker_stats_embedded.tsv", sep="\t", index=False)
    print("\nWrote analysis_embedded.json, marker_values_embedded.tsv, "
          "marker_stats_embedded.tsv")


if __name__ == "__main__":
    main()
