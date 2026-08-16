"""TACSTD2 (Trop-2) / CLDN4 (Claudin-4) vs response / outcome analysis.

Reality of the 2015-2018 human-lung ICI landscape (see clinical_label_index.tsv):
  * The in-window human LUNG series that carry immune-checkpoint TREATMENT-RESPONSE
    labels use TARGETED immune panels (GSE93157 nCounter 775-gene; GSE110390
    21-gene IFN-gamma) that DO NOT measure TACSTD2 or CLDN4.
  * The in-window human LUNG series that DO measure TACSTD2/CLDN4 genome-wide
    (GSE72094 LUAD array; GSE81089 NSCLC RNA-seq) are surgical cohorts with
    OVERALL-SURVIVAL labels, not ICI response.
  * GSE91061 (Riaz 2017) has genome-wide expression AND RECIST ICI response, but
    is MELANOMA (nivolumab), not lung.

So this script does exactly what the data allow, and labels each analysis by what
it actually is (ICI response vs overall survival; lung vs melanoma). Nothing is
imputed or invented.

Analyses:
  A. Gene-availability audit of TACSTD2/CLDN4/CD274 across all processed matrices.
  B. GSE91061 (MELANOMA, nivolumab): TACSTD2/CLDN4/CD274 vs RECIST response
     (pre-treatment), Mann-Whitney U + boxplots.  [ICI response benchmark]
  C. GSE72094 (LUAD array): TACSTD2/CLDN4/CD274 vs overall survival
     (median-split log-rank) + Spearman vs CD274.  [lung, survival]
  D. GSE81089 (NSCLC RNA-seq): same as C.                    [lung, survival]

Outputs: results/fable_geo_2015_2018/analysis/*
"""
import glob
import gzip
import json
import math
import os
from datetime import datetime

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "fable_geo_2015_2018")
DATA = os.path.join(OUT, "data")
CLIN = os.path.join(OUT, "clinical")
ANA = os.path.join(OUT, "analysis")
os.makedirs(ANA, exist_ok=True)

GENES = {"TACSTD2": {"symbol": "TACSTD2", "entrez": "4070", "ensembl": "ENSG00000184292"},
         "CLDN4": {"symbol": "CLDN4", "entrez": "1364", "ensembl": "ENSG00000189143"},
         "CD274": {"symbol": "CD274", "entrez": "29126", "ensembl": "ENSG00000120217"}}

results = {}


def read_clinical(acc):
    path = os.path.join(CLIN, f"{acc}_clinical.tsv")
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in f]
    return rows


def logrank(times, events, groups):
    """Two-group log-rank test. groups in {0,1}. Returns (chi2, p, expected/obs)."""
    times = np.asarray(times, float)
    events = np.asarray(events, int)
    groups = np.asarray(groups, int)
    event_times = np.unique(times[events == 1])
    O1 = E1 = V = 0.0
    for t in event_times:
        at_risk = times >= t
        n = at_risk.sum()
        n1 = (at_risk & (groups == 1)).sum()
        d = ((times == t) & (events == 1)).sum()
        d1 = ((times == t) & (events == 1) & (groups == 1)).sum()
        if n <= 1:
            continue
        E1 += d * n1 / n
        O1 += d1
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    if V == 0:
        return 0.0, 1.0, (O1, E1)
    chi2 = (O1 - E1) ** 2 / V
    p = stats.chi2.sf(chi2, 1)
    return chi2, p, (O1, E1)


# ----------------------------------------------------------------------------
# A. Gene availability audit across all processed expression matrices
# ----------------------------------------------------------------------------
def audit_gene_availability():
    audit = []
    patt = os.path.join(DATA, "*", "*")
    for path in sorted(glob.glob(patt)):
        name = os.path.basename(path)
        acc = os.path.basename(os.path.dirname(path))
        low = name.lower()
        if not (low.endswith((".txt.gz", ".tsv.gz", ".csv.gz")) and
                any(k in low for k in ("fpkm", "count", "tpm", "expr", "raw_data",
                                       "norm", "gene", "matrix", "processed",
                                       "featurecounts", "cufflinks", "values"))):
            continue
        if "series_matrix" in low:
            continue
        # scan the FULL file; a gene is "present" if its symbol / Entrez / Ensembl
        # id appears as a row identifier (first field) or as a whole token.
        hit = {g: False for g in GENES}
        toks = {g: {info["symbol"].upper(), info["entrez"], info["ensembl"]}
                for g, info in GENES.items()}
        sep = "," if low.endswith(".csv.gz") else "\t"
        try:
            with gzip.open(path, "rt", errors="replace") as f:
                for ln in f:
                    first = ln.split(sep, 1)[0].strip().strip('"').upper()
                    for g in GENES:
                        if not hit[g] and first in toks[g]:
                            hit[g] = True
                    if all(hit.values()):
                        break
        except Exception:
            continue
        audit.append({"accession": acc, "file": name,
                      "TACSTD2": hit["TACSTD2"], "CLDN4": hit["CLDN4"],
                      "CD274": hit["CD274"]})
    with open(os.path.join(ANA, "gene_availability_audit.tsv"), "w") as f:
        f.write("accession\tfile\tTACSTD2\tCLDN4\tCD274\n")
        for r in audit:
            f.write(f"{r['accession']}\t{r['file']}\t{r['TACSTD2']}\t{r['CLDN4']}\t{r['CD274']}\n")
    results["gene_availability"] = audit
    print(f"[A] gene-availability audit: {len(audit)} processed files scanned")


# ----------------------------------------------------------------------------
# B. GSE91061 melanoma nivolumab: response benchmark
# ----------------------------------------------------------------------------
def analyze_gse91061():
    fpkm = os.path.join(DATA, "GSE91061",
                        "GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz")
    with gzip.open(fpkm, "rt") as f:
        header = f.readline().rstrip("\n").split(",")
        cols = [c.strip().strip('"') for c in header][1:]
        want = {info["entrez"]: g for g, info in GENES.items()}
        expr = {}
        for ln in f:
            parts = ln.rstrip("\n").split(",")
            gid = parts[0].strip().strip('"')
            if gid in want:
                expr[want[gid]] = np.array([float(x) for x in parts[1:]], float)
    # response per sample title
    clin = read_clinical("GSE91061")
    title2resp = {r["title"]: r.get("response", "") for r in clin}
    title2visit = {r["title"]: r.get("visit_(pre_or_on_treatment)", "") for r in clin}

    def group_for(title):
        resp = title2resp.get(title, "")
        visit = title2visit.get(title, "")
        if visit.lower() != "pre":
            return None
        if resp == "PRCR":
            return "responder"
        if resp == "PD":
            return "nonresponder"
        if resp == "SD":
            return "SD"
        return None

    groups = [group_for(t) for t in cols]
    out = {"dataset": "GSE91061", "type": "MELANOMA nivolumab (RECIST response)",
           "n_pre": sum(1 for g in groups if g), "genes": {}}
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, gene in zip(axes, ["TACSTD2", "CLDN4", "CD274"]):
        vals = expr[gene]
        r = np.array([np.log2(vals[i] + 1) for i in range(len(cols))
                      if groups[i] == "responder"])
        nr = np.array([np.log2(vals[i] + 1) for i in range(len(cols))
                       if groups[i] == "nonresponder"])
        sd = np.array([np.log2(vals[i] + 1) for i in range(len(cols))
                       if groups[i] == "SD"])
        u, p = stats.mannwhitneyu(r, nr, alternative="two-sided") if len(r) and len(nr) else (float("nan"), float("nan"))
        out["genes"][gene] = {
            "n_responder": int(len(r)), "n_nonresponder": int(len(nr)), "n_SD": int(len(sd)),
            "median_responder_log2fpkm": float(np.median(r)) if len(r) else None,
            "median_nonresponder_log2fpkm": float(np.median(nr)) if len(nr) else None,
            "mannwhitney_U": float(u), "p_value": float(p),
        }
        ax.boxplot([r, sd, nr], tick_labels=["PR/CR", "SD", "PD"], showfliers=False)
        for j, arr in enumerate([r, sd, nr], start=1):
            ax.scatter(np.random.normal(j, 0.05, len(arr)), arr, s=12, alpha=0.6)
        ax.set_title(f"{gene}  (PR/CR vs PD p={p:.3g})")
        ax.set_ylabel("log2(FPKM+1)")
    fig.suptitle("GSE91061 (MELANOMA, nivolumab) pre-treatment — ICI response benchmark")
    fig.tight_layout()
    fig.savefig(os.path.join(ANA, "GSE91061_response_boxplots.png"), dpi=130)
    plt.close(fig)
    results["GSE91061"] = out
    print(f"[B] GSE91061 melanoma: {out['n_pre']} pre-tx samples grouped; "
          f"TACSTD2 p={out['genes']['TACSTD2']['p_value']:.3g}, "
          f"CLDN4 p={out['genes']['CLDN4']['p_value']:.3g}")


# ----------------------------------------------------------------------------
# Survival helper for lung datasets
# ----------------------------------------------------------------------------
def survival_analysis(acc, expr_by_gene, sample_order, times, events, kind):
    out = {"dataset": acc, "type": kind, "n": int(len(sample_order)),
           "n_events": int(np.sum(events)), "genes": {}}
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, gene in zip(axes, ["TACSTD2", "CLDN4", "CD274"]):
        vals = expr_by_gene[gene]
        med = np.median(vals)
        grp = (vals > med).astype(int)
        chi2, p, (o1, e1) = logrank(times, events, grp)
        # Spearman vs CD274
        rho, prho = stats.spearmanr(vals, expr_by_gene["CD274"])
        out["genes"][gene] = {
            "median_expr": float(med),
            "logrank_chi2": float(chi2), "logrank_p": float(p),
            "n_high": int(grp.sum()), "n_low": int((1 - grp).sum()),
            "spearman_vs_CD274_rho": float(rho), "spearman_vs_CD274_p": float(prho),
        }
        for gv, lab in [(1, "high"), (0, "low")]:
            mask = grp == gv
            t = np.sort(times[mask])
            e = events[mask][np.argsort(times[mask])]
            n = mask.sum()
            surv = [1.0]
            tt = [0.0]
            at_risk = n
            s = 1.0
            for i in range(len(t)):
                if e[i] == 1:
                    s *= (1 - 1.0 / at_risk)
                at_risk -= 1
                tt.append(t[i]); surv.append(s)
            ax.step(tt, surv, where="post", label=f"{gene} {lab} (n={n})")
        ax.set_title(f"{gene} median-split (log-rank p={p:.3g})")
        ax.set_xlabel("time"); ax.set_ylabel("survival"); ax.set_ylim(0, 1.02)
        ax.legend(fontsize=7)
    fig.suptitle(f"{acc} — {kind}")
    fig.tight_layout()
    fig.savefig(os.path.join(ANA, f"{acc}_survival_km.png"), dpi=130)
    plt.close(fig)
    results[acc] = out
    print(f"[surv] {acc}: n={out['n']} events={out['n_events']} | "
          + ", ".join(f"{g} logrank p={out['genes'][g]['logrank_p']:.3g}"
                      for g in ['TACSTD2', 'CLDN4']))


def analyze_gse72094():
    # map probes -> gene using GPL15048 family annotation
    ann = os.path.join(DATA, "GSE72094", "GPL15048_family.txt")
    probe2gene = {}
    with open(ann, errors="replace") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 4 and p[3] in ("TACSTD2", "CLDN4", "CD274"):
                probe2gene[p[0]] = p[3]
    # read matrix expression rows for those probes
    mat = os.path.join(DATA, "GSE72094", "GSE72094_series_matrix.txt.gz")
    samples = None
    gene_vals = {"TACSTD2": [], "CLDN4": [], "CD274": []}
    with gzip.open(mat, "rt", errors="replace") as f:
        in_tab = False
        for ln in f:
            if ln.startswith("!series_matrix_table_begin"):
                in_tab = True
                continue
            if ln.startswith("!series_matrix_table_end"):
                break
            if in_tab:
                parts = ln.rstrip("\n").split("\t")
                if samples is None:
                    samples = [x.strip().strip('"') for x in parts[1:]]
                    continue
                pid = parts[0].strip().strip('"')
                if pid in probe2gene:
                    gene_vals[probe2gene[pid]].append(
                        np.array([float(x) if x not in ("", "null") else np.nan
                                  for x in parts[1:]], float))
    expr = {g: np.nanmean(np.vstack(v), axis=0) for g, v in gene_vals.items()}
    clin = {r["gsm"]: r for r in read_clinical("GSE72094")}
    times, events, keep = [], [], []
    for i, gsm in enumerate(samples):
        r = clin.get(gsm, {})
        st = r.get("survival_time_in_days", "")
        vs = r.get("vital_status", "").strip().lower()
        try:
            t = float(st)
        except ValueError:
            continue
        if vs not in ("alive", "dead"):
            continue
        times.append(t); events.append(1 if vs == "dead" else 0); keep.append(i)
    keep = np.array(keep)
    expr_k = {g: expr[g][keep] for g in expr}
    survival_analysis("GSE72094",
                      {g: np.asarray(v) for g, v in expr_k.items()},
                      [samples[i] for i in keep],
                      np.array(times, float), np.array(events, int),
                      "LUAD resected cohort — OVERALL SURVIVAL (not ICI). log2 array intensity")


def analyze_gse81089():
    fpkm = os.path.join(DATA, "GSE81089", "GSE81089_FPKM_cufflinks.tsv.gz")
    want = {info["ensembl"]: g for g, info in GENES.items()}
    with gzip.open(fpkm, "rt") as f:
        cols = f.readline().rstrip("\n").split("\t")[1:]
        expr = {}
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if p[0] in want:
                expr[want[p[0]]] = np.array([float(x) for x in p[1:]], float)
    expr = {g: np.log2(v + 1) for g, v in expr.items()}
    # clinical keyed by title (L###T); match FPKM columns
    clin = {r["title"]: r for r in read_clinical("GSE81089")}
    times, events, keep = [], [], []
    for i, c in enumerate(cols):
        r = clin.get(c, {})
        dead = r.get("dead", "").strip()
        sdate = r.get("surgery_date", "").strip()
        vdate = r.get("vital_date", "").strip()
        if dead not in ("0", "1"):
            continue
        try:
            d0 = datetime.strptime(sdate, "%Y-%m-%d")
            d1 = datetime.strptime(vdate, "%Y-%m-%d")
        except ValueError:
            continue
        t = (d1 - d0).days
        if t < 0:
            continue
        times.append(float(t)); events.append(int(dead)); keep.append(i)
    keep = np.array(keep)
    expr_k = {g: expr[g][keep] for g in expr}
    survival_analysis("GSE81089", expr_k, [cols[i] for i in keep],
                      np.array(times, float), np.array(events, int),
                      "NSCLC resected cohort — OVERALL SURVIVAL (not ICI). log2(FPKM+1)")


def main():
    audit_gene_availability()
    analyze_gse91061()
    analyze_gse72094()
    analyze_gse81089()
    with open(os.path.join(ANA, "analysis_results.json"), "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote analysis outputs -> {ANA}")


if __name__ == "__main__":
    main()
