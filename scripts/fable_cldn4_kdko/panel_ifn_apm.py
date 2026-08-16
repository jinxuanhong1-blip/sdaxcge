"""Focused test of ONE directional hypothesis:

    CLDN4 knockdown / knockout OPENS (up-regulates) the interferon /
    MHC-class-I / antigen-presentation-machinery (APM) program.

The user's private RNA-seq shows this UP direction for
IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A, TAP1, TAP2. Here we test, honestly
and one-sidedly, which public CLDN4 loss-of-function datasets reproduce the
SAME (upward) direction.

Panels:
  CORE      : the 8 genes the user named.
  EXTENDED  : full type-I ISG + MHC-I + APM program.

Datasets: GSE50927 (mouse lung Cldn4 KO, edgeR table with p/FDR),
GSE207704 (human breast CRISPR CLDN4-/-, per-group FPKM),
GSE22493 (human ovarian CLDN4 siRNA, 2-colour array; low confidence).
"""
import gzip
import glob
import json
import math
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

RAW = "../../results/fable_cldn4_kdko/raw"
OUT = "../../results/fable_cldn4_kdko"

CORE = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A", "TAP1", "TAP2"]

ISG = ["ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2",
       "IFIT1", "IFIT2", "IFIT3", "IFI27", "IFI44", "IFI44L", "IFI6",
       "USP18", "HERC5", "DDX58", "IFIH1", "XAF1", "BST2", "ISG20",
       "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "GBP1", "GBP2"]
MHC1_APM = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1",
            "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "NLRC5", "ERAP1"]
EXTENDED = ISG + MHC1_APM

# human -> mouse ortholog candidates for panel genes (Title-case default).
MOUSE = {
    "HLA-A": ["H2-K1", "H2-D1", "H2-Q7"], "HLA-B": ["H2-K1", "H2-D1"],
    "HLA-C": ["H2-D1", "H2-K1"], "HLA-E": ["H2-T23"], "HLA-F": [],
    "B2M": ["B2m"], "IFI27": ["Ifi27l2a", "Ifi27l2b", "Ifi27"],
    "IFIT1": ["Ifit1", "Ifit1bl1"], "OAS1": ["Oas1a", "Oas1g", "Oas1b"],
    "OASL": ["Oasl1", "Oasl2"], "MX1": ["Mx1"], "MX2": ["Mx2"],
    "GBP1": ["Gbp2", "Gbp3"], "IFI44L": [], "ISG20": ["Isg20"],
    "PSMB10": ["Psmb10"], "ERAP1": ["Erap1"],
}


def to_mouse(sym):
    if sym in MOUSE:
        return MOUSE[sym]
    return [sym[0].upper() + sym[1:].lower()] if "-" not in sym else [sym]


# ---------------- per-dataset fold-change dictionaries ----------------
def fc_gse50927():
    df = pd.read_csv(f"{RAW}/GSE50927_Cldn4lungWTvsKOgenes.csv.gz")
    df = df.rename(columns={"Marker.Symbol": "symbol"}).dropna(subset=["symbol"])
    df["absfc"] = df["logFC"].abs()
    df = df.sort_values("absfc", ascending=False).drop_duplicates("symbol")
    d = {r.symbol: (r.logFC, r.PValue, r.FDR) for r in df.itertuples()}
    return d, dict(zip(df.symbol, df.logFC))  # (per-gene), (fc-only for bg)


def fc_gse207704():
    df = pd.read_csv(f"{RAW}/GSE207704_CLDN4_RNAseq.txt.gz", sep="\t").rename(columns={
        "MCF7_CLDN4KO_FPKM (fpkm)": "MCF7_KO", "MCF7_WT_FPKM (fpkm)": "MCF7_WT",
        "T47D_CLDN4KO_FPKM (fpkm)": "T47D_KO", "T47D_WT_FPKM (fpkm)": "T47D_WT",
        "gene_short_name": "symbol"}).dropna(subset=["symbol"])
    for c in ["MCF7_KO", "MCF7_WT", "T47D_KO", "T47D_WT"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    g = df.groupby("symbol")[["MCF7_KO", "MCF7_WT", "T47D_KO", "T47D_WT"]].sum()
    ps = 1.0
    g["MCF7"] = np.log2((g.MCF7_KO + ps) / (g.MCF7_WT + ps))
    g["T47D"] = np.log2((g.T47D_KO + ps) / (g.T47D_WT + ps))
    g["mean"] = g[["MCF7", "T47D"]].mean(axis=1)
    g["expr"] = g[["MCF7_KO", "MCF7_WT", "T47D_KO", "T47D_WT"]].max(axis=1) >= 1.0
    per = {i: (g.at[i, "mean"], g.at[i, "MCF7"], g.at[i, "T47D"],
               bool(g.at[i, "expr"]),
               (g.at[i, "MCF7_WT"], g.at[i, "MCF7_KO"], g.at[i, "T47D_WT"], g.at[i, "T47D_KO"]))
           for i in g.index}
    bg = {i: g.at[i, "mean"] for i in g.index if g.at[i, "expr"]}
    return per, bg


def _parse_scanarray(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        lines = fh.read().splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("BEGIN DATA"):
            header = lines[i + 1].split("\t"); start = i + 2; break
    c = {x.strip(): j for j, x in enumerate(header)}
    per = {}
    for ln in lines[start:]:
        f = ln.split("\t")
        if len(f) <= c["Ch2 Median - B"] or not re.match(r"^\d+$", f[0].strip()):
            continue
        sym = f[c["Name"]].strip().strip('"').split("--")[0].strip()
        try:
            a = float(f[c["Ch1 Median - B"]]); b = float(f[c["Ch2 Median - B"]])
        except ValueError:
            continue
        if a <= 0 or b <= 0:
            continue
        per.setdefault(sym, []).append(math.log2(b / a))
    gl = {g: float(np.mean(v)) for g, v in per.items()}
    med = np.median(list(gl.values()))
    return {g: v - med for g, v in gl.items()}


def fc_gse22493():
    arrays = [_parse_scanarray(p) for p in sorted(glob.glob(f"{RAW}/gse22493/GSM*.txt.gz"))]
    genes = set().union(*[set(a) for a in arrays])
    per, bg = {}, {}
    for g in genes:
        vals = [a[g] for a in arrays if g in a]
        if len(vals) < 2:
            continue
        m = float(np.mean(vals))
        p = float(stats.ttest_1samp(vals, 0.0).pvalue) if len(vals) == 3 else None
        per[g] = (m, p, len(vals))
        bg[g] = m
    return per, bg


# ---------------- honest directional statistics ----------------
def directional_stats(fcs, bg_fcs):
    """fcs: list of panel log-fold-changes (finite). Returns up-direction tests."""
    fcs = [x for x in fcs if np.isfinite(x)]
    n = len(fcs)
    n_up = sum(1 for x in fcs if x > 0)
    out = {"n_panel_detected": n, "n_up": n_up, "n_down": n - n_up,
           "median_logFC": (round(float(np.median(fcs)), 4) if n else None),
           "mean_logFC": (round(float(np.mean(fcs)), 4) if n else None)}
    # sign test: are more than half up?
    if n >= 3:
        out["sign_test_p_up"] = float(stats.binomtest(n_up, n, 0.5, alternative="greater").pvalue)
        # one-sided Wilcoxon signed-rank vs 0 (panel shifted positive)
        try:
            out["wilcoxon_p_up"] = float(stats.wilcoxon(fcs, alternative="greater").pvalue)
        except ValueError:
            out["wilcoxon_p_up"] = None
        # one-sided Mann-Whitney panel vs background (panel > background)
        bgv = [v for v in bg_fcs.values() if np.isfinite(v)]
        out["mwu_p_panel_gt_bg"] = float(
            stats.mannwhitneyu(fcs, bgv, alternative="greater").pvalue)
    return out


def evaluate(name, per, bg, panel_syms, mouse=False):
    rows = []
    fcs = []
    for h in panel_syms:
        syms = to_mouse(h) if mouse else [h]
        found = None
        for s in syms:
            if s in per:
                found = (s, per[s]); break
        if found is None:
            rows.append({"human": h, "matched": None, "logFC": None,
                         "direction": "not detected"})
            continue
        s, val = found
        fc = val[0]
        fcs.append(fc)
        row = {"human": h, "matched": s, "logFC": round(float(fc), 4),
               "direction": "UP" if fc > 0 else "down"}
        # attach p/FDR where available
        if name == "GSE50927":
            row["PValue"], row["FDR"] = float(val[1]), float(val[2])
        elif name == "GSE207704":
            row["log2FC_MCF7"], row["log2FC_T47D"] = round(val[1], 3), round(val[2], 3)
        elif name == "GSE22493":
            row["p_1samp"], row["n_arrays"] = val[1], val[2]
        rows.append(row)
    return rows, directional_stats(fcs, bg)


def main():
    d50_per, d50_bg = fc_gse50927()
    d20_per, d20_bg = fc_gse207704()
    d22_per, d22_bg = fc_gse22493()

    report = {}
    for name, per, bg, mouse in [
        ("GSE50927", d50_per, d50_bg, True),
        ("GSE207704", d20_per, d20_bg, False),
        ("GSE22493", d22_per, d22_bg, False)]:
        report[name] = {}
        for panelname, panel in [("CORE", CORE), ("EXTENDED", EXTENDED)]:
            rows, ds = evaluate(name, per, bg, panel, mouse=mouse)
            report[name][panelname] = {"per_gene": rows, "directional": ds}

    with open(f"{OUT}/panel_ifn_apm_results.json", "w") as f:
        json.dump(report, f, indent=1)

    # tidy CSV of CORE panel across datasets
    csv_rows = []
    for name in ["GSE50927", "GSE207704", "GSE22493"]:
        for r in report[name]["CORE"]["per_gene"]:
            csv_rows.append({"dataset": name, **r})
    pd.DataFrame(csv_rows).to_csv(f"{OUT}/panel_core_by_dataset.csv", index=False)

    # console summary
    for name in ["GSE50927", "GSE207704", "GSE22493"]:
        print("\n" + "=" * 68)
        print(name)
        for panelname in ["CORE", "EXTENDED"]:
            ds = report[name][panelname]["directional"]
            print(f"  [{panelname}] up {ds['n_up']}/{ds['n_panel_detected']} "
                  f"median logFC={ds['median_logFC']} "
                  f"| sign p(up)={ds.get('sign_test_p_up')} "
                  f"wilcoxon p(up)={ds.get('wilcoxon_p_up')} "
                  f"MWU p(panel>bg)={ds.get('mwu_p_panel_gt_bg')}")
        print("  CORE per-gene:")
        for r in report[name]["CORE"]["per_gene"]:
            print("     ", r)
    print("\nSaved panel_ifn_apm_results.json + panel_core_by_dataset.csv")


if __name__ == "__main__":
    main()
