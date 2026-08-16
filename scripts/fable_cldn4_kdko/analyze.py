"""CLDN4-loss differential-expression analysis across 3 verified datasets.

Datasets (all true CLDN4 perturbations, processed data only):
  GSE207704  human breast cancer, CRISPR CLDN4-/- vs WT (T47D, MCF7), RNA-seq FPKM
  GSE50927   mouse lung, Cldn4 KO vs WT, RNA-seq (author edgeR DE table)
  GSE22493   human ovarian SKOV3, CLDN4 siRNA-KD vs CLDN4-high control, 2-colour array

Questions: does CLDN4 loss change (a) TACSTD2/Trop-2, (b) junction genes,
(c) IFN/immune genes?

Outputs tidy per-gene tables + gene-set-level statistics to results/.
"""
import gzip
import glob
import io
import json
import math
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

import gene_sets as GS

RAW = "../../results/fable_cldn4_kdko/raw"
OUT = "../../results/fable_cldn4_kdko"


def mannwhitney_set_vs_bg(values_by_gene, set_symbols):
    """Two-sided Mann-Whitney U comparing set logFC vs background logFC."""
    set_syms = set(set_symbols)
    in_set = {g: v for g, v in values_by_gene.items() if g in set_syms and np.isfinite(v)}
    bg = {g: v for g, v in values_by_gene.items() if g not in set_syms and np.isfinite(v)}
    if len(in_set) < 3:
        return {"n_set": len(in_set), "n_bg": len(bg), "U": None, "p": None,
                "median_set": (np.median(list(in_set.values())) if in_set else None),
                "median_bg": (np.median(list(bg.values())) if bg else None)}
    U, p = stats.mannwhitneyu(list(in_set.values()), list(bg.values()),
                              alternative="two-sided")
    return {"n_set": len(in_set), "n_bg": len(bg), "U": float(U), "p": float(p),
            "median_set": float(np.median(list(in_set.values()))),
            "median_bg": float(np.median(list(bg.values())))}


# ---------------------------------------------------------------------------
# GSE50927 : mouse lung Cldn4 KO vs WT (author edgeR table)
# ---------------------------------------------------------------------------
def analyze_gse50927():
    df = pd.read_csv(f"{RAW}/GSE50927_Cldn4lungWTvsKOgenes.csv.gz")
    df = df.rename(columns={"Marker.Symbol": "symbol"})
    df = df.dropna(subset=["symbol"])
    # collapse duplicate symbols by strongest |logFC|
    df["absfc"] = df["logFC"].abs()
    df = df.sort_values("absfc", ascending=False).drop_duplicates("symbol")
    fc_by_gene = dict(zip(df["symbol"], df["logFC"]))

    ms = GS.mouse_symbols()
    # per-gene table for markers
    marker_rows = []
    for setname, syms in [("TACSTD2", ms["TACSTD2"]), ("JUNCTION", ms["JUNCTION"]),
                          ("IFN_IMMUNE", ms["IFN_IMMUNE"])]:
        sub = df[df["symbol"].isin(syms)]
        for _, r in sub.iterrows():
            marker_rows.append({"set": setname, "symbol": r["symbol"],
                                "logFC_KOvsWT": round(float(r["logFC"]), 4),
                                "PValue": float(r["PValue"]), "FDR": float(r["FDR"])})
    # gene-set shift stats (logFC of set vs background)
    setstats = {}
    for setname in ["JUNCTION", "IFN_IMMUNE"]:
        setstats[setname] = mannwhitney_set_vs_bg(fc_by_gene, ms[setname])
    return {"dataset": "GSE50927", "n_genes": len(df),
            "markers": marker_rows, "setstats": setstats,
            "cldn4_logFC": fc_by_gene.get("Cldn4")}


# ---------------------------------------------------------------------------
# GSE207704 : human breast CRISPR CLDN4-/- vs WT (FPKM, per group)
# ---------------------------------------------------------------------------
def analyze_gse207704():
    df = pd.read_csv(f"{RAW}/GSE207704_CLDN4_RNAseq.txt.gz", sep="\t")
    df = df.rename(columns={
        "MCF7_CLDN4KO_FPKM (fpkm)": "MCF7_KO", "MCF7_WT_FPKM (fpkm)": "MCF7_WT",
        "T47D_CLDN4KO_FPKM (fpkm)": "T47D_KO", "T47D_WT_FPKM (fpkm)": "T47D_WT",
        "gene_short_name": "symbol"})
    df = df.dropna(subset=["symbol"])
    for c in ["MCF7_KO", "MCF7_WT", "T47D_KO", "T47D_WT"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    # collapse duplicate symbols by summing FPKM (isoform loci -> gene)
    g = df.groupby("symbol")[["MCF7_KO", "MCF7_WT", "T47D_KO", "T47D_WT"]].sum()
    ps = 1.0  # FPKM pseudocount
    g["log2FC_MCF7"] = np.log2((g["MCF7_KO"] + ps) / (g["MCF7_WT"] + ps))
    g["log2FC_T47D"] = np.log2((g["T47D_KO"] + ps) / (g["T47D_WT"] + ps))
    g["log2FC_mean"] = g[["log2FC_MCF7", "log2FC_T47D"]].mean(axis=1)
    # require detectable expression in at least one condition to reduce noise
    expressed = (g[["MCF7_KO", "MCF7_WT", "T47D_KO", "T47D_WT"]].max(axis=1) >= 1.0)
    ge = g[expressed]
    fc_by_gene = dict(zip(ge.index, ge["log2FC_mean"]))

    hs = GS.human_symbols()
    marker_rows = []
    for setname, syms in [("TACSTD2", hs["TACSTD2"]), ("JUNCTION", hs["JUNCTION"]),
                          ("IFN_IMMUNE", hs["IFN_IMMUNE"])]:
        sub = g[g.index.isin(syms)]
        for sym, r in sub.iterrows():
            marker_rows.append({"set": setname, "symbol": sym,
                                "MCF7_WT": round(float(r["MCF7_WT"]), 3),
                                "MCF7_KO": round(float(r["MCF7_KO"]), 3),
                                "T47D_WT": round(float(r["T47D_WT"]), 3),
                                "T47D_KO": round(float(r["T47D_KO"]), 3),
                                "log2FC_MCF7": round(float(r["log2FC_MCF7"]), 4),
                                "log2FC_T47D": round(float(r["log2FC_T47D"]), 4),
                                "log2FC_mean": round(float(r["log2FC_mean"]), 4)})
    setstats = {}
    for setname in ["JUNCTION", "IFN_IMMUNE"]:
        setstats[setname] = mannwhitney_set_vs_bg(fc_by_gene, hs[setname])
    return {"dataset": "GSE207704", "n_genes_expressed": int(expressed.sum()),
            "markers": marker_rows, "setstats": setstats,
            "cldn4": {"MCF7_WT": float(g.loc["CLDN4", "MCF7_WT"]) if "CLDN4" in g.index else None,
                      "MCF7_KO": float(g.loc["CLDN4", "MCF7_KO"]) if "CLDN4" in g.index else None,
                      "T47D_WT": float(g.loc["CLDN4", "T47D_WT"]) if "CLDN4" in g.index else None,
                      "T47D_KO": float(g.loc["CLDN4", "T47D_KO"]) if "CLDN4" in g.index else None}}


# ---------------------------------------------------------------------------
# GSE22493 : ovarian SKOV3, two-colour array, Cy5=CLDN4 KD vs Cy3=CLDN4-high control
# ---------------------------------------------------------------------------
def parse_scanarray(path):
    """Return dict symbol -> mean log2(Ch2/Ch1) for one array (median-centered)."""
    with gzip.open(path, "rt", errors="replace") as fh:
        lines = fh.read().splitlines()
    # find data block
    start = None
    header = None
    for i, ln in enumerate(lines):
        if ln.startswith("BEGIN DATA"):
            header = lines[i + 1].split("\t")
            start = i + 2
            break
    cols = {c.strip(): j for j, c in enumerate(header)}
    ci_name = cols["Name"]
    ci_ch1 = cols["Ch1 Median - B"]
    ci_ch2 = cols["Ch2 Median - B"]
    per_gene = {}
    for ln in lines[start:]:
        if ln.startswith("END DATA") or not ln.strip():
            continue
        f = ln.split("\t")
        if len(f) <= ci_ch2:
            continue
        if not re.match(r"^\d+$", f[0].strip()):  # data rows start with integer Index
            continue
        name = f[ci_name].strip().strip('"')
        sym = name.split("--")[0].strip()
        if not sym or sym.upper().startswith("EMPTY") or sym.upper() == "BLANK":
            continue
        try:
            ch1 = float(f[ci_ch1]); ch2 = float(f[ci_ch2])
        except ValueError:
            continue
        if ch1 <= 0 or ch2 <= 0:
            continue
        lr = math.log2(ch2 / ch1)  # Cy5(KD) / Cy3(control)
        per_gene.setdefault(sym, []).append(lr)
    # collapse replicate spots within array
    gene_lr = {g: float(np.mean(v)) for g, v in per_gene.items()}
    # global median-center (dye/loading normalization)
    med = np.median(list(gene_lr.values()))
    return {g: v - med for g, v in gene_lr.items()}


def analyze_gse22493():
    arrays = [parse_scanarray(p) for p in sorted(glob.glob(f"{RAW}/gse22493/GSM*.txt.gz"))]
    all_genes = set().union(*[set(a) for a in arrays])
    rows = {}
    for g in all_genes:
        vals = [a[g] for a in arrays if g in a]
        if len(vals) < 2:
            continue
        mean = float(np.mean(vals))
        if len(vals) == 3:
            t, p = stats.ttest_1samp(vals, 0.0)
        else:
            t, p = (np.nan, np.nan)
        rows[g] = {"log2FC_KDvsCtrl": mean, "n_arrays": len(vals),
                   "t": (float(t) if np.isfinite(t) else None),
                   "p": (float(p) if p == p else None)}
    fc_by_gene = {g: r["log2FC_KDvsCtrl"] for g, r in rows.items()}

    hs = GS.human_symbols()
    marker_rows = []
    for setname, syms in [("TACSTD2", hs["TACSTD2"]), ("JUNCTION", hs["JUNCTION"]),
                          ("IFN_IMMUNE", hs["IFN_IMMUNE"])]:
        for sym in syms:
            if sym in rows:
                r = rows[sym]
                marker_rows.append({"set": setname, "symbol": sym,
                                    "log2FC_KDvsCtrl": round(r["log2FC_KDvsCtrl"], 4),
                                    "n_arrays": r["n_arrays"],
                                    "p_1samp_t": r["p"]})
    setstats = {}
    for setname in ["JUNCTION", "IFN_IMMUNE"]:
        setstats[setname] = mannwhitney_set_vs_bg(fc_by_gene, hs[setname])
    return {"dataset": "GSE22493", "n_genes": len(rows),
            "markers": marker_rows, "setstats": setstats,
            "cldn4": rows.get("CLDN4")}


def gse22493_kd_diagnostic():
    """Per-array log2(KD/control) for CLDN4 to document knockdown efficiency."""
    rows = []
    for p in sorted(glob.glob(f"{RAW}/gse22493/GSM*.txt.gz")):
        gsm = os.path.basename(p).split(".")[0]
        with gzip.open(p, "rt", errors="replace") as fh:
            lines = fh.read().splitlines()
        for i, ln in enumerate(lines):
            if ln.startswith("BEGIN DATA"):
                header = lines[i + 1].split("\t"); start = i + 2; break
        cols = {c.strip(): j for j, c in enumerate(header)}
        for ln in lines[start:]:
            f = ln.split("\t")
            if len(f) <= cols["Ch2 Median - B"] or not re.match(r"^\d+$", f[0].strip()):
                continue
            sym = f[cols["Name"]].strip().strip('"').split("--")[0].strip()
            if sym == "CLDN4":
                ch1 = f[cols["Ch1 Median - B"]]; ch2 = f[cols["Ch2 Median - B"]]
                try:
                    c1 = float(ch1); c2 = float(ch2)
                except ValueError:
                    continue
                lr = (math.log2(c2 / c1) if c1 > 0 and c2 > 0 else None)
                rows.append({"gsm": gsm, "CLDN4_ctrl_Cy3": c1, "CLDN4_KD_Cy5": c2,
                             "log2_KD_over_ctrl": (round(lr, 3) if lr is not None else "NA(nonpos)")})
    return rows


def _write_csv(path, rows, cols):
    pd.DataFrame(rows)[cols].to_csv(path, index=False)


def main():
    results = {}
    results["GSE50927"] = analyze_gse50927()
    results["GSE207704"] = analyze_gse207704()
    results["GSE22493"] = analyze_gse22493()
    results["GSE22493_kd_diagnostic"] = gse22493_kd_diagnostic()
    with open(f"{OUT}/analysis_results.json", "w") as f:
        json.dump(results, f, indent=1)

    # ---- tidy CSV exports ----
    # per-dataset marker tables
    for acc in ["GSE207704", "GSE50927", "GSE22493"]:
        md = pd.DataFrame(results[acc]["markers"])
        md.insert(0, "dataset", acc)
        md.to_csv(f"{OUT}/markers_{acc}.csv", index=False)
    # combined gene-set statistics
    ss_rows = []
    for acc in ["GSE207704", "GSE50927", "GSE22493"]:
        for setname, st in results[acc]["setstats"].items():
            ss_rows.append({"dataset": acc, "gene_set": setname, **st})
    pd.DataFrame(ss_rows).to_csv(f"{OUT}/geneset_stats.csv", index=False)
    # kd diagnostic
    pd.DataFrame(results["GSE22493_kd_diagnostic"]).to_csv(
        f"{OUT}/GSE22493_CLDN4_kd_diagnostic.csv", index=False)

    # Human-readable dump
    for acc in ["GSE207704", "GSE50927", "GSE22493"]:
        r = results[acc]
        print("\n" + "=" * 70)
        print(acc, "| CLDN4 control:", r.get("cldn4", r.get("cldn4_logFC")))
        print("-- gene-set shift (Mann-Whitney set vs background) --")
        for s, st in r["setstats"].items():
            print(f"   {s:11} n_set={st['n_set']:3} median_set={st['median_set']} "
                  f"median_bg={st['median_bg']} p={st['p']}")
        print("-- markers --")
        for m in r["markers"]:
            print("   ", {k: v for k, v in m.items()})
    print("\nSaved analysis_results.json")


if __name__ == "__main__":
    main()
