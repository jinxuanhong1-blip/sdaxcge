#!/usr/bin/env python3
"""Analyse TACSTD2 / CLDN4 detectability in GSE305086 whole-blood bulk arrays
(Affymetrix HG-U133 Plus 2.0, GPL570) from aNSCLC patients under PD-(L)1 IO.

Samples (173): 76 baseline + 76 follow-up (after 4 cycles) + 10 non-smoker and
11 smoker age-matched controls. Values are log2 RMA-normalised intensities.

There are no public per-patient response labels in the GEO metadata (only the
treatment regimen: 1LIO / 2LIO / CHTIO), so this dataset is used for:
  1. Detectability of TACSTD2 / CLDN4 in whole blood, expressed as the
     percentile rank of each probe's mean intensity among all 26,453 probes
     (a probe near the bottom is at array background = effectively absent).
  2. Disease effect: patient baseline vs age-matched controls (Mann-Whitney U).
  3. Treatment effect: paired baseline vs follow-up (Wilcoxon signed-rank).
T-cell genes FOXP3 / GATA3 / PDCD1 serve as positive controls: the source study
reports them deregulated in patient blood, so a working pipeline should move.

Outputs (results/fable_blood_ici/):
    gse305086_sample_annotation.csv
    gse305086_target_detectability.csv
    gse305086_contrasts.csv
    gse305086_detectability.png
    gse305086_summary.json
"""
import os
import re
import gzip
import json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GEO_DIR = os.environ.get("GEO_DIR", "/tmp/geo")
RESULTS = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..",
                                       "results", "fable_blood_ici"))
EXPR = os.path.join(GEO_DIR, "GSE305086_Expression_matrix_final.csv.gz")
SERIES = os.path.join(GEO_DIR, "GSE305086_series_matrix.txt.gz")

# Probe -> gene map (from GPL570 annotation; probes confirmed present in matrix).
PROBE_GENE = {
    "202285_s_at": "TACSTD2", "202286_s_at": "TACSTD2", "202287_s_at": "TACSTD2",
    "201428_at": "CLDN4", "1569421_at": "CLDN4",
    "201839_s_at": "EPCAM",
    "207238_s_at": "PTPRC", "212587_s_at": "PTPRC",
    "210164_at": "GZMB", "214617_at": "PRF1",
    "213539_at": "CD3D", "213915_at": "NKG7", "205758_at": "CD8A",
    "209602_s_at": "GATA3", "209603_at": "GATA3", "209604_s_at": "GATA3",
    "207634_at": "PDCD1",
}
# FOXP3 probes (221333_at/221334_s_at/224211_at) are filtered out of this
# processed matrix, so GATA3 and PDCD1 are used as the T-cell positive controls.
TARGET_PROBES = ["202285_s_at", "202286_s_at", "202287_s_at",
                 "201428_at", "1569421_at"]


def parse_series(path):
    """Return DataFrame indexed by short sample name (matrix column) with
    geo_accession, sample_type, group, and control smoking status."""
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()
    def row(tag):
        for ln in lines:
            if ln.startswith(tag):
                return [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        return None
    titles = row("!Sample_title")
    gsm = row("!Sample_geo_accession")
    stype, group = None, None
    for ln in lines:
        if ln.startswith("!Sample_characteristics_ch1"):
            vals = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
            if vals and vals[0].startswith("sample type:"):
                stype = [v.split(":", 1)[1].strip() for v in vals]
            elif vals and vals[0].startswith("group:"):
                group = [v.split(":", 1)[1].strip() for v in vals]
    recs = []
    for t, g, st, gr in zip(titles, gsm, stype, group):
        short = t.split(" ", 1)[0]           # "C1 baseline" -> "C1"
        ctrl = ("non-smoker" if "non-smoker" in t else
                "smoker" if "smoker" in t else "")
        recs.append({"sample": short, "geo_accession": g,
                     "sample_type": st, "group": gr, "control_status": ctrl})
    return pd.DataFrame(recs).set_index("sample")


def main():
    os.makedirs(RESULTS, exist_ok=True)
    meta = parse_series(SERIES)
    meta.to_csv(os.path.join(RESULTS, "gse305086_sample_annotation.csv"))

    expr = pd.read_csv(EXPR, sep=";", index_col=0)
    expr.index.name = "probe"
    print(f"expression matrix: {expr.shape[0]} probes x {expr.shape[1]} samples")

    # keep only panel probes actually present in this (filtered) matrix
    global PROBE_GENE
    PROBE_GENE = {p: g for p, g in PROBE_GENE.items() if p in expr.index}
    print(f"panel probes present: {len(PROBE_GENE)}")

    # align sample metadata to matrix columns
    meta = meta.reindex(expr.columns)
    assert meta["sample_type"].notna().all(), "unmatched matrix columns"

    baseline = meta.index[meta.sample_type == "baseline"]
    followup = meta.index[meta.sample_type == "follow-up"]
    control = meta.index[meta.sample_type == "control group"]
    print(f"baseline={len(baseline)} follow-up={len(followup)} control={len(control)}")

    # ---- 1. Detectability: percentile of probe mean among all probes -------
    probe_mean = expr.mean(axis=1)
    pct_rank = probe_mean.rank(pct=True) * 100
    det_rows = []
    for probe, gene in PROBE_GENE.items():
        if probe not in expr.index:
            continue
        det_rows.append({
            "probe": probe, "gene": gene,
            "is_target": gene in ("TACSTD2", "CLDN4"),
            "mean_log2": round(float(probe_mean[probe]), 4),
            "median_log2": round(float(expr.loc[probe].median()), 4),
            "percentile_rank": round(float(pct_rank[probe]), 2),
        })
    det = pd.DataFrame(det_rows).sort_values("percentile_rank", ascending=False)
    det.to_csv(os.path.join(RESULTS, "gse305086_target_detectability.csv"),
               index=False)

    # ---- 2 & 3. Contrasts on gene-level values (mean of probes per gene) ----
    gene_expr = expr.loc[list(PROBE_GENE)].groupby(
        pd.Series(PROBE_GENE)).mean()   # gene x sample

    con_rows = []
    for gene in gene_expr.index:
        b = gene_expr.loc[gene, baseline].to_numpy(float)
        c = gene_expr.loc[gene, control].to_numpy(float)
        u, p = mannwhitneyu(b, c, alternative="two-sided")
        con_rows.append({
            "gene": gene, "contrast": "baseline_vs_control", "test": "MannWhitneyU",
            "baseline_mean": round(b.mean(), 4), "other_mean": round(c.mean(), 4),
            "n_baseline": len(b), "n_other": len(c), "stat": float(u),
            "p_value": float(p),
        })
        # paired baseline vs follow-up (match patient id: Cx vs Cx+4)
        pairs = [(s, s + "+4") for s in baseline
                 if (s + "+4") in gene_expr.columns]
        b_p = gene_expr.loc[gene, [p0 for p0, _ in pairs]].to_numpy(float)
        f_p = gene_expr.loc[gene, [p1 for _, p1 in pairs]].to_numpy(float)
        w, pw = wilcoxon(b_p, f_p)
        con_rows.append({
            "gene": gene, "contrast": "baseline_vs_followup_paired",
            "test": "WilcoxonSignedRank",
            "baseline_mean": round(b_p.mean(), 4), "other_mean": round(f_p.mean(), 4),
            "n_baseline": len(b_p), "n_other": len(f_p), "stat": float(w),
            "p_value": float(pw),
        })
    con = pd.DataFrame(con_rows)
    con.to_csv(os.path.join(RESULTS, "gse305086_contrasts.csv"), index=False)

    # ---- 4. Figure: where do targets sit in the whole array distribution ----
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(probe_mean, bins=80, color="#cccccc", edgecolor="none")
    ymax = ax.get_ylim()[1]
    colors = {"TACSTD2": "#d62728", "CLDN4": "#9467bd", "EPCAM": "#ff7f0e",
              "PTPRC": "#1f77b4", "FOXP3": "#2ca02c"}
    seen = set()
    for probe, gene in PROBE_GENE.items():
        if gene not in colors or probe not in probe_mean:
            continue
        x = probe_mean[probe]
        ax.axvline(x, color=colors[gene], lw=1.5,
                   label=gene if gene not in seen else None)
        seen.add(gene)
    ax.set_xlabel("probe mean log2 intensity (all 173 whole-blood samples)")
    ax.set_ylabel("number of probes")
    ax.set_title("GSE305086 whole-blood array: TACSTD2 / CLDN4 vs array background\n"
                 "(left tail = background/absent, right = expressed)")
    ax.legend(title="probe (gene)")
    fig.tight_layout()
    fig.savefig(os.path.join(RESULTS, "gse305086_detectability.png"), dpi=140)

    # ---- 5. Summary --------------------------------------------------------
    summary = {
        "dataset": "GSE305086", "platform": "GPL570 (HG-U133 Plus 2.0)",
        "n_samples": int(expr.shape[1]), "n_probes": int(expr.shape[0]),
        "n_baseline": len(baseline), "n_followup": len(followup),
        "n_control": len(control),
        "targets": {},
    }
    for gene in ("TACSTD2", "CLDN4"):
        sub = det[det.gene == gene]
        summary["targets"][gene] = {
            "probes": sub["probe"].tolist(),
            "max_percentile_rank": float(sub["percentile_rank"].max()),
            "mean_log2_range": [float(sub["mean_log2"].min()),
                                float(sub["mean_log2"].max())],
        }
    with open(os.path.join(RESULTS, "gse305086_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print("\nDetectability (percentile rank among all probes):")
    print(det.to_string(index=False))
    print("\nContrasts (targets + controls):")
    print(con[con.gene.isin(["TACSTD2", "CLDN4", "FOXP3", "GATA3", "PDCD1"])]
          .to_string(index=False))


if __name__ == "__main__":
    main()
