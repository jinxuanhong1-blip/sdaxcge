#!/usr/bin/env python3
"""C-CPE pharmacologic CLDN4-axis analog: GSE22421 ± GSE22493.

Primary numbers are author-deposited series-matrix VALUE columns:
  GSE22421 VALUE = log2(C-CPE treated / untreated)
  GSE22493 VALUE = log2(CLDN4 siRNA / CLDN4 overexpression)

Both series are two-color Operon 60-mer arrays on GPL10555, n=3.
Gene symbols are mapped from ORF or DESCRIPTION prefix (same rule as
the C4 GSE22493 slice). TACSTD2 is requested but is not annotated on
this platform; that absence is reported, not filled.

ScanArray Cy5/Cy3 is a sensitivity check only.
One-sample t-tests and gene-set tests with n=3 are descriptive.

Usage:
  python3 methods/ccpe_cldn4_analog/scripts/download_data.py
  python3 methods/ccpe_cldn4_analog/scripts/analyze.py
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from gene_sets import (
    ALIASES,
    APM,
    AXIS,
    HALLMARK_FOCUS,
    IFN_IMMUNE,
    PRIORITY,
    TJ,
)
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
OUT = os.path.join(REPO, "methods", "ccpe_cldn4_analog")
TABLES = os.path.join(OUT, "tables")
FIGS = os.path.join(OUT, "figures")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)

DATA_DIR = os.environ.get("CCPE_CLDN4_DATA", "/tmp/ccpe_cldn4_analog")
GMT_PATH = os.path.join(HERE, "h.all.v2023.2.Hs.symbols.gmt")

STUDIES = {
    "GSE22421": {
        "gsms": ["GSM557400", "GSM557401", "GSM557402"],
        "matrix": os.path.join(DATA_DIR, "gse22421/GSE22421_series_matrix.txt.gz"),
        "soft": os.path.join(DATA_DIR, "gse22421/GSE22421_family.soft.gz"),
        "scan_dir": os.path.join(DATA_DIR, "gse22421/scanarray"),
        "contrast": "C-CPE 5ug/ml 72h (Cy5) / untreated (Cy3)",
        "value_def": "log2(C-CPE / untreated)",
        "model": "SKOV-3 ovarian cancer cell line",
        "perturbation": "C-CPE (CLDN3/4 binder)",
        "control": "untreated",
        "n_arrays": 3,
        "pmid": "21123456",
    },
    "GSE22493": {
        "gsms": ["GSM558700", "GSM558701", "GSM558702"],
        "matrix": os.path.join(DATA_DIR, "gse22493/GSE22493_series_matrix.txt.gz"),
        "soft": os.path.join(DATA_DIR, "gse22493/GSE22493_family.soft.gz"),
        "scan_dir": os.path.join(DATA_DIR, "gse22493/scanarray"),
        "contrast": "CLDN4 lentiviral siRNA (Cy5) / CLDN4 overexpression (Cy3)",
        "value_def": "log2(CLDN4 siRNA / CLDN4 OE)",
        "model": "SKOV-3-IP-Luc ovarian cancer cell line",
        "perturbation": "CLDN4 lentiviral siRNA",
        "control": "CLDN4 overexpression (not scramble/WT)",
        "n_arrays": 3,
        "pmid": "",
    },
}


def _split_fields(line: str) -> list[str]:
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]


def load_soft_annotations(path: str) -> dict[str, dict]:
    ann: dict[str, dict] = {}
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                next(fh)
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table or line.startswith("!"):
                continue
            fields = _split_fields(line)
            if not fields or not fields[0]:
                continue
            probe = fields[0]
            desc = fields[1] if len(fields) > 1 else ""
            orf = fields[2].strip() if len(fields) > 2 else ""
            raw = orf.upper()
            mapped_via = "ORF"
            if not raw and "--" in desc:
                raw = desc.split("--", 1)[0].strip().upper()
                mapped_via = "DESCRIPTION_PREFIX"
            if not raw:
                continue
            symbol = ALIASES.get(raw, raw)
            alias_note = f"{raw}->{symbol}" if symbol != raw else ""
            ann[probe] = {
                "probe": probe,
                "symbol": symbol,
                "raw_symbol": raw,
                "mapped_via": mapped_via,
                "alias_note": alias_note,
                "description": desc,
            }
    return ann


def load_series_matrix(path: str) -> tuple[list[str], dict[str, list[float | None]]]:
    values: dict[str, list[float | None]] = {}
    samples: list[str] = []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            fields = _split_fields(line)
            if fields[0] == "ID_REF":
                samples = fields[1:]
                continue
            probe = fields[0]
            row: list[float | None] = []
            for cell in fields[1:]:
                if cell == "":
                    row.append(None)
                else:
                    row.append(float(cell))
            while len(row) < len(samples):
                row.append(None)
            values[probe] = row[: len(samples)]
    return samples, values


def parse_scanarray(path: str) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        in_data = False
        header: list[str] = []
        for line in fh:
            if line.startswith("BEGIN DATA"):
                in_data = True
                continue
            if not in_data:
                continue
            if not header:
                header = line.rstrip("\n").split("\t")
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0].isdigit():
                continue
            cols = {name: i for i, name in enumerate(header)}
            probe = fields[0]
            try:
                ch1 = float(fields[cols["Ch1 Median - B"]])
                ch2 = float(fields[cols["Ch2 Median - B"]])
            except (ValueError, KeyError):
                continue
            rec = {
                "ch1_cy3": ch1,
                "ch2_cy5": ch2,
                "log2_cy5_over_cy3": None,
            }
            if ch1 > 0 and ch2 > 0:
                rec["log2_cy5_over_cy3"] = math.log2(ch2 / ch1)
            out[probe] = rec
    return out


def _finite(xs: list[float | None]) -> list[float]:
    return [x for x in xs if x is not None and np.isfinite(x)]


def collapse_gene(
    probe_ids: list[str],
    matrix: dict[str, list[float | None]],
    n_samples: int,
) -> tuple[list[float | None], int]:
    per_array: list[float | None] = []
    for j in range(n_samples):
        vals = []
        for pid in probe_ids:
            row = matrix.get(pid)
            if row is None or j >= len(row) or row[j] is None:
                continue
            vals.append(row[j])
        per_array.append(float(np.median(vals)) if vals else None)
    return per_array, len(probe_ids)


def onesamp(vals: list[float]) -> tuple[float | None, float | None]:
    if len(vals) < 2:
        return None, None
    t, p = stats.ttest_1samp(vals, 0.0)
    if not np.isfinite(t) or not np.isfinite(p):
        return None, None
    return float(t), float(p)


def direction_of(median: float | None) -> str:
    if median is None:
        return "NA"
    if median > 0:
        return "UP"
    if median < 0:
        return "DOWN"
    return "ZERO"


def fmt(x: float | None, nd: int = 4) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return ""
    return f"{x:.{nd}f}"


def write_tsv(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore", delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})


def load_gmt(path: str) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sets[parts[0]] = [g.upper() for g in parts[2:] if g]
    return sets


def gsea_preranked(
    ranked: list[tuple[str, float]],
    gene_set: set[str],
    n_perm: int = 1000,
    seed: int = 1,
) -> dict:
    """Classic weighted KS enrichment on a pre-ranked list. Descriptive."""
    genes = [g for g, _ in ranked]
    scores = np.array([s for _, s in ranked], dtype=float)
    n = len(genes)
    in_set = np.array([g in gene_set for g in genes], dtype=bool)
    n_hit = int(in_set.sum())
    if n_hit < 5 or n_hit > n - 5:
        return {
            "n_in_set_measured": n_hit,
            "es": None,
            "nes": None,
            "p_perm": None,
        }
    abs_s = np.abs(scores)
    hit_w = np.where(in_set, abs_s, 0.0)
    hit_den = hit_w.sum()
    if hit_den == 0:
        return {
            "n_in_set_measured": n_hit,
            "es": None,
            "nes": None,
            "p_perm": None,
        }
    miss_step = 1.0 / (n - n_hit)
    walk = np.cumsum(np.where(in_set, hit_w / hit_den, -miss_step))
    es = float(walk[np.argmax(np.abs(walk))])

    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_perm):
        perm = rng.permutation(in_set)
        hw = np.where(perm, abs_s, 0.0)
        hd = hw.sum()
        if hd == 0:
            continue
        w = np.cumsum(np.where(perm, hw / hd, -miss_step))
        null.append(float(w[np.argmax(np.abs(w))]))
    null = np.array(null)
    if es >= 0:
        pos = null[null >= 0]
        nes = es / pos.mean() if len(pos) else None
        p = float((null >= es).mean()) if len(null) else None
    else:
        neg = null[null < 0]
        nes = es / abs(neg.mean()) if len(neg) else None
        p = float((null <= es).mean()) if len(null) else None
    return {
        "n_in_set_measured": n_hit,
        "es": es,
        "nes": float(nes) if nes is not None else None,
        "p_perm": p,
    }


def analyze_study(name: str, spec: dict, ann: dict[str, dict]) -> dict:
    samples, matrix = load_series_matrix(spec["matrix"])
    # keep declared GSM order
    idx = [samples.index(g) for g in spec["gsms"]]
    reordered: dict[str, list[float | None]] = {}
    for pid, row in matrix.items():
        reordered[pid] = [row[j] if j < len(row) else None for j in idx]
    matrix = reordered
    samples = spec["gsms"]
    n = len(samples)

    by_symbol: dict[str, list[str]] = defaultdict(list)
    for pid, a in ann.items():
        by_symbol[a["symbol"]].append(pid)

    # background: genes with >=2 arrays
    gene_means: dict[str, float] = {}
    gene_medians: dict[str, float] = {}
    gene_n: dict[str, int] = {}
    for sym, pids in by_symbol.items():
        per, _ = collapse_gene(pids, matrix, n)
        vals = _finite(per)
        if len(vals) >= 2:
            gene_means[sym] = float(np.mean(vals))
            gene_medians[sym] = float(np.median(vals))
            gene_n[sym] = len(vals)

    scan: dict[str, dict[str, dict]] = {}
    for gsm in samples:
        sp = os.path.join(spec["scan_dir"], f"{gsm}.txt.gz")
        if os.path.exists(sp):
            scan[gsm] = parse_scanarray(sp)

    def scan_collapse(pids: list[str]) -> list[float | None]:
        out: list[float | None] = []
        for gsm in samples:
            vals = []
            recs = scan.get(gsm, {})
            for pid in pids:
                rec = recs.get(pid)
                if rec and rec["log2_cy5_over_cy3"] is not None:
                    vals.append(rec["log2_cy5_over_cy3"])
            out.append(float(np.median(vals)) if vals else None)
        return out

    panels = {
        "AXIS": AXIS,
        "PRIORITY": PRIORITY,
        "APM": APM,
        "IFN_IMMUNE": IFN_IMMUNE,
        "TJ": TJ,
    }

    gene_rows: list[dict] = []
    probe_rows: list[dict] = []
    for panel, genes in panels.items():
        for sym in genes:
            pids = by_symbol.get(sym, [])
            per, n_probes = collapse_gene(pids, matrix, n) if pids else ([None] * n, 0)
            vals = _finite(per)
            mean_v = float(np.mean(vals)) if vals else None
            med_v = float(np.median(vals)) if vals else None
            t, p = onesamp(vals)
            scan_per = scan_collapse(pids) if pids else [None] * n
            scan_vals = _finite(scan_per)
            scan_med = float(np.median(scan_vals)) if scan_vals else None
            notes = []
            if pids:
                vias = sorted({ann[p]["mapped_via"] for p in pids if p in ann})
                aliases = sorted({ann[p]["alias_note"] for p in pids if p in ann and ann[p]["alias_note"]})
                if aliases:
                    notes.extend(aliases)
                notes.extend(vias)
            if not pids:
                notes.append("NOT_ON_GPL10555")
            row = {
                "study": name,
                "panel": panel,
                "symbol": sym,
                "measured": "yes" if vals else "no",
                "n_probes_annotated": n_probes,
                "n_arrays_with_value": len(vals),
                "mean_log2": fmt(mean_v),
                "median_log2": fmt(med_v),
                "t_1samp": fmt(t),
                "p_1samp": fmt(p, 4) if p is None else (f"{p:.4g}" if p < 0.001 else fmt(p, 4)),
                "direction": direction_of(med_v),
                "scanarray_median_log2": fmt(scan_med),
                "scanarray_direction": direction_of(scan_med),
                "scanarray_n": len(scan_vals),
                "mapping_note": ";".join(notes),
            }
            for i, gsm in enumerate(samples):
                row[gsm] = fmt(per[i])
                row[f"{gsm}_scan"] = fmt(scan_per[i])
            gene_rows.append(row)
            for pid in pids:
                prow = {
                    "study": name,
                    "panel": panel,
                    "symbol": sym,
                    "probe": pid,
                    "mapped_via": ann[pid]["mapped_via"],
                    "raw_symbol": ann[pid]["raw_symbol"],
                    "alias_note": ann[pid]["alias_note"],
                    "description": ann[pid]["description"],
                }
                dep = matrix.get(pid, [None] * n)
                for i, gsm in enumerate(samples):
                    prow[gsm] = fmt(dep[i] if i < len(dep) else None)
                    rec = scan.get(gsm, {}).get(pid, {})
                    prow[f"{gsm}_scan"] = fmt(rec.get("log2_cy5_over_cy3"))
                    prow[f"{gsm}_ch1"] = fmt(rec.get("ch1_cy3"))
                    prow[f"{gsm}_ch2"] = fmt(rec.get("ch2_cy5"))
                probe_rows.append(prow)

    # BH within each panel among measured genes
    for panel in panels:
        idxs = [i for i, r in enumerate(gene_rows) if r["panel"] == panel and r["p_1samp"] != ""]
        if not idxs:
            continue
        ps = [float(gene_rows[i]["p_1samp"]) for i in idxs]
        _, q, _, _ = multipletests(ps, method="fdr_bh")
        for i, qi in zip(idxs, q):
            gene_rows[i]["q_BH_within_panel"] = fmt(float(qi), 4)

    # gene-set stats vs background
    bg = np.array(list(gene_means.values()))
    set_rows = []
    for set_name, genes in panels.items():
        measured = [g for g in genes if g in gene_means]
        missing = [g for g in genes if g not in gene_means]
        set_vals = np.array([gene_means[g] for g in measured]) if measured else np.array([])
        if len(set_vals) >= 3 and len(bg) > len(set_vals):
            u, p = stats.mannwhitneyu(set_vals, bg, alternative="two-sided")
            u = float(u)
            p = float(p)
        else:
            u, p = None, None
        set_rows.append(
            {
                "study": name,
                "gene_set": set_name,
                "n_in_list": len(genes),
                "n_measured_ge2_arrays": len(measured),
                "n_missing_or_lt2": len(missing),
                "missing_or_lt2": ",".join(missing),
                "n_mean_up": int((set_vals > 0).sum()) if len(set_vals) else 0,
                "n_mean_down": int((set_vals < 0).sum()) if len(set_vals) else 0,
                "median_set_mean_log2": fmt(float(np.median(set_vals)) if len(set_vals) else None),
                "mean_set_mean_log2": fmt(float(np.mean(set_vals)) if len(set_vals) else None),
                "median_background_mean_log2": fmt(float(np.median(bg))),
                "n_background_genes": len(bg),
                "mannwhitney_U": fmt(u, 1) if u is not None else "",
                "mannwhitney_p": f"{p:.4g}" if p is not None and p < 0.001 else fmt(p, 4),
                "set_shift": (
                    "UP"
                    if len(set_vals) and float(np.median(set_vals)) > float(np.median(bg))
                    else ("DOWN" if len(set_vals) else "NA")
                ),
            }
        )

    # Hallmark preranked GSEA (descriptive)
    gmt = load_gmt(GMT_PATH) if os.path.exists(GMT_PATH) else {}
    ranked = sorted(gene_means.items(), key=lambda kv: kv[1], reverse=True)
    gsea_rows = []
    for hs in HALLMARK_FOCUS:
        genes = set(gmt.get(hs, []))
        res = gsea_preranked(ranked, genes)
        gsea_rows.append(
            {
                "study": name,
                "gene_set": hs,
                "n_in_gmt": len(genes),
                "n_in_set_measured": res["n_in_set_measured"],
                "es": fmt(res["es"]),
                "nes": fmt(res["nes"]),
                "p_perm_1000": fmt(res["p_perm"], 4),
                "note": "preranked on gene-mean deposited log2; n=3 arrays; descriptive",
            }
        )

    # CLDN4 diagnostic
    cldn4_pids = by_symbol.get("CLDN4", [])
    cldn4_dep, _ = collapse_gene(cldn4_pids, matrix, n) if cldn4_pids else ([None] * n, 0)
    cldn4_scan = scan_collapse(cldn4_pids) if cldn4_pids else [None] * n
    diag_rows = []
    for i, gsm in enumerate(samples):
        rec = {}
        if cldn4_pids:
            rec = scan.get(gsm, {}).get(cldn4_pids[0], {})
        diag_rows.append(
            {
                "study": name,
                "gsm": gsm,
                "probe": ",".join(cldn4_pids),
                "deposited_log2": fmt(cldn4_dep[i]),
                "scanarray_log2": fmt(cldn4_scan[i]),
                "ch1_cy3_median_minus_B": fmt(rec.get("ch1_cy3")),
                "ch2_cy5_median_minus_B": fmt(rec.get("ch2_cy5")),
            }
        )

    return {
        "name": name,
        "spec": spec,
        "samples": samples,
        "gene_rows": gene_rows,
        "probe_rows": probe_rows,
        "set_rows": set_rows,
        "gsea_rows": gsea_rows,
        "diag_rows": diag_rows,
        "gene_means": gene_means,
        "n_background": len(gene_means),
        "n_mapped_symbols": len(by_symbol),
        "n_mapped_probes": len(ann),
        "tacstd2_on_platform": "TACSTD2" in by_symbol,
        "cldn4_median": (
            float(np.median(_finite(cldn4_dep))) if _finite(cldn4_dep) else None
        ),
        "cldn4_per_array": cldn4_dep,
        "cldn4_n": len(_finite(cldn4_dep)),
    }


def panel_summary(rows: list[dict], panel: str) -> dict:
    sub = [r for r in rows if r["panel"] == panel]
    meas = [r for r in sub if r["measured"] == "yes"]
    up = [r for r in meas if r["direction"] == "UP"]
    down = [r for r in meas if r["direction"] == "DOWN"]
    return {
        "n_list": len(sub),
        "n_measured": len(meas),
        "n_up": len(up),
        "n_down": len(down),
        "up": [r["symbol"] for r in up],
        "down": [r["symbol"] for r in down],
    }


def make_figures(results: dict[str, dict]) -> None:
    # Fig 1: key axis + priority genes, GSE22421
    r21 = results["GSE22421"]
    gsms = r21["samples"]
    key_syms = ["CLDN4", "CLDN3"] + PRIORITY
    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    xs = []
    labels = []
    for i, sym in enumerate(key_syms):
        row = next(
            (x for x in r21["gene_rows"] if x["panel"] in ("AXIS", "PRIORITY") and x["symbol"] == sym),
            None,
        )
        if row is None or row["measured"] != "yes":
            continue
        vals = [float(row[g]) for g in gsms if row[g] != ""]
        jitter = np.linspace(-0.12, 0.12, len(vals)) if len(vals) > 1 else [0]
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=36, c="#1f4e79", zorder=3)
        if vals:
            ax.hlines(float(np.median(vals)), i - 0.28, i + 0.28, colors="#c0392b", lw=2)
        xs.append(i)
        labels.append(sym)
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("deposited log2 (C-CPE / untreated)")
    ax.set_title("GSE22421 SKOV-3 C-CPE vs untreated  (n=3 two-color arrays)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "fig1_gse22421_axis_priority.png"), dpi=160)
    fig.savefig(os.path.join(FIGS, "fig1_gse22421_axis_priority.pdf"))
    plt.close(fig)

    # Fig 2: APM + TJ medians for GSE22421
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.0), sharey=True)
    for ax, panel, title in (
        (axes[0], "APM", "MHC-I / APM"),
        (axes[1], "TJ", "Tight junction"),
    ):
        rows = [x for x in r21["gene_rows"] if x["panel"] == panel and x["measured"] == "yes"]
        rows = sorted(rows, key=lambda r: float(r["median_log2"]))
        y = np.arange(len(rows))
        meds = [float(r["median_log2"]) for r in rows]
        colors = ["#c0392b" if m < 0 else "#1f4e79" for m in meds]
        ax.barh(y, meds, color=colors, height=0.7)
        ax.axvline(0, color="0.3", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([r["symbol"] for r in rows], fontsize=8)
        ax.set_title(f"GSE22421 {title}\nmedian log2, n=3")
        ax.set_xlabel("median deposited log2")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "fig2_gse22421_apm_tj.png"), dpi=160)
    fig.savefig(os.path.join(FIGS, "fig2_gse22421_apm_tj.pdf"))
    plt.close(fig)

    # Fig 3: paired GSE22421 vs GSE22493 directions for shared measured genes
    if "GSE22493" in results:
        r93 = results["GSE22493"]
        pairs = []
        seen = set()
        for panel in ("AXIS", "PRIORITY", "APM", "TJ"):
            for row in r21["gene_rows"]:
                if row["panel"] != panel or row["measured"] != "yes":
                    continue
                if (panel, row["symbol"]) in seen:
                    continue
                o = next(
                    (
                        x
                        for x in r93["gene_rows"]
                        if x["panel"] == panel and x["symbol"] == row["symbol"] and x["measured"] == "yes"
                    ),
                    None,
                )
                if o is None:
                    continue
                seen.add((panel, row["symbol"]))
                pairs.append(
                    {
                        "panel": panel,
                        "symbol": row["symbol"],
                        "x": float(row["median_log2"]),
                        "y": float(o["median_log2"]),
                    }
                )
        fig, ax = plt.subplots(figsize=(6.4, 6.2))
        colors = {
            "AXIS": "#8e44ad",
            "PRIORITY": "#1f4e79",
            "APM": "#16a085",
            "TJ": "#d35400",
        }
        for panel, c in colors.items():
            sub = [p for p in pairs if p["panel"] == panel]
            if not sub:
                continue
            ax.scatter(
                [p["x"] for p in sub],
                [p["y"] for p in sub],
                s=42,
                c=c,
                label=panel,
                zorder=3,
            )
            for p in sub:
                if p["symbol"] in ("CLDN4", "CLDN3", "HLA-A", "IFI27", "IFIT1", "ISG15", "B2M", "TAP1"):
                    ax.annotate(p["symbol"], (p["x"], p["y"]), fontsize=7, xytext=(4, 3), textcoords="offset points")
        ax.axhline(0, color="0.5", lw=0.7)
        ax.axvline(0, color="0.5", lw=0.7)
        ax.set_xlabel("GSE22421 median log2 (C-CPE / untreated)")
        ax.set_ylabel("GSE22493 median log2 (CLDN4 siRNA / OE)")
        ax.set_title("Paired deposited medians  (same GPL10555 mapping)")
        ax.legend(frameon=False, fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.tight_layout()
        fig.savefig(os.path.join(FIGS, "fig3_paired_gse22421_gse22493.png"), dpi=160)
        fig.savefig(os.path.join(FIGS, "fig3_paired_gse22421_gse22493.pdf"))
        plt.close(fig)

    # Fig 4: gene-set median vs background
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    sets_order = ["PRIORITY", "APM", "IFN_IMMUNE", "TJ"]
    studies = [k for k in ("GSE22421", "GSE22493") if k in results]
    width = 0.35
    x = np.arange(len(sets_order))
    for i, st in enumerate(studies):
        vals = []
        for sname in sets_order:
            row = next(r for r in results[st]["set_rows"] if r["gene_set"] == sname)
            vals.append(float(row["median_set_mean_log2"]) if row["median_set_mean_log2"] else 0)
        ax.bar(x + (i - 0.5) * width, vals, width=width, label=st)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(sets_order)
    ax.set_ylabel("median of gene-mean deposited log2")
    ax.set_title("Gene-set location vs 0  (n=3 arrays each; descriptive)")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "fig4_geneset_medians.png"), dpi=160)
    fig.savefig(os.path.join(FIGS, "fig4_geneset_medians.pdf"))
    plt.close(fig)

    # Fig 5: IFN_IMMUNE per-gene for GSE22421
    rows = [x for x in r21["gene_rows"] if x["panel"] == "IFN_IMMUNE" and x["measured"] == "yes"]
    rows = sorted(rows, key=lambda r: float(r["median_log2"]))
    fig, ax = plt.subplots(figsize=(7.2, 9.5))
    y = np.arange(len(rows))
    meds = [float(r["median_log2"]) for r in rows]
    colors = ["#c0392b" if m < 0 else "#1f4e79" for m in meds]
    ax.barh(y, meds, color=colors, height=0.75)
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([r["symbol"] for r in rows], fontsize=7)
    ax.set_xlabel("median deposited log2 (C-CPE / untreated)")
    ax.set_title("GSE22421 IFN_IMMUNE panel  (measured genes, n=3)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "fig5_gse22421_ifn_panel.png"), dpi=160)
    fig.savefig(os.path.join(FIGS, "fig5_gse22421_ifn_panel.pdf"))
    plt.close(fig)


def main() -> None:
    # Same platform; load annotations once from GSE22421 family SOFT.
    ann21 = load_soft_annotations(STUDIES["GSE22421"]["soft"])
    ann93 = load_soft_annotations(STUDIES["GSE22493"]["soft"])
    # Mapping audit: both should be GPL10555 with identical probe IDs.
    keys21 = set(ann21)
    keys93 = set(ann93)
    mapping_rows = [
        {
            "item": "platform",
            "GSE22421": "GPL10555",
            "GSE22493": "GPL10555",
            "note": "same custom BWH Human Release 3.0 Operon 60-mer",
        },
        {
            "item": "n_probes_with_symbol",
            "GSE22421": str(len(ann21)),
            "GSE22493": str(len(ann93)),
            "note": "ORF or DESCRIPTION prefix",
        },
        {
            "item": "n_shared_mapped_probes",
            "GSE22421": str(len(keys21 & keys93)),
            "GSE22493": str(len(keys21 & keys93)),
            "note": "probe ID intersection",
        },
        {
            "item": "TACSTD2_annotated",
            "GSE22421": "no",
            "GSE22493": "no",
            "note": "only TACSTD1 (EPCAM-family) is on GPL10555; no TROP2/GA733-1/M1S1",
        },
        {
            "item": "ISG15_mapping",
            "GSE22421": "G1P2->ISG15",
            "GSE22493": "G1P2->ISG15",
            "note": "historical IFI-15K symbol",
        },
    ]
    write_tsv(
        os.path.join(TABLES, "mapping_audit.tsv"),
        mapping_rows,
        ["item", "GSE22421", "GSE22493", "note"],
    )

    results = {}
    results["GSE22421"] = analyze_study("GSE22421", STUDIES["GSE22421"], ann21)
    results["GSE22493"] = analyze_study("GSE22493", STUDIES["GSE22493"], ann93)

    # sample table
    sample_rows = [
        {
            "study": "GSE22421",
            "gsm": "GSM557400",
            "title": "CCPE Replicate 1",
            "cell_line": "SKOV-3",
            "ch1_Cy3": "untreated 72h",
            "ch2_Cy5": "C-CPE 5ug/ml 72h",
            "deposited_value": "log2 ratio treated/control",
        },
        {
            "study": "GSE22421",
            "gsm": "GSM557401",
            "title": "CCPE Replicate 2",
            "cell_line": "SKOV-3",
            "ch1_Cy3": "untreated 72h",
            "ch2_Cy5": "C-CPE 5ug/ml 72h",
            "deposited_value": "log2 ratio treated/control",
        },
        {
            "study": "GSE22421",
            "gsm": "GSM557402",
            "title": "CCPE Replicate 3",
            "cell_line": "SKOV-3",
            "ch1_Cy3": "untreated 72h",
            "ch2_Cy5": "C-CPE 5ug/ml 72h",
            "deposited_value": "log2 ratio treated/control",
        },
        {
            "study": "GSE22493",
            "gsm": "GSM558700",
            "title": "CLDN4 Replicate 1",
            "cell_line": "SKOV-3-IP-Luc",
            "ch1_Cy3": "CLDN4 overexpression (control)",
            "ch2_Cy5": "CLDN4 knockdown (lentiviral siRNA)",
            "deposited_value": "log2 ratio knockdown/control",
        },
        {
            "study": "GSE22493",
            "gsm": "GSM558701",
            "title": "CLDN4 Replicate 2",
            "cell_line": "SKOV-3-IP-Luc",
            "ch1_Cy3": "CLDN4 overexpression (control)",
            "ch2_Cy5": "CLDN4 knockdown (lentiviral siRNA)",
            "deposited_value": "log2 ratio knockdown/control",
        },
        {
            "study": "GSE22493",
            "gsm": "GSM558702",
            "title": "CLDN4 Replicate 3",
            "cell_line": "SKOV-3-IP-Luc",
            "ch1_Cy3": "CLDN4 overexpression (control)",
            "ch2_Cy5": "CLDN4 knockdown (lentiviral siRNA)",
            "deposited_value": "log2 ratio knockdown/control",
        },
    ]
    write_tsv(
        os.path.join(TABLES, "sample_table.tsv"),
        sample_rows,
        ["study", "gsm", "title", "cell_line", "ch1_Cy3", "ch2_Cy5", "deposited_value"],
    )

    all_genes = results["GSE22421"]["gene_rows"] + results["GSE22493"]["gene_rows"]
    gsm_cols_21 = STUDIES["GSE22421"]["gsms"]
    gsm_cols_93 = STUDIES["GSE22493"]["gsms"]
    gene_fields = [
        "study",
        "panel",
        "symbol",
        "measured",
        "n_probes_annotated",
        "n_arrays_with_value",
        *gsm_cols_21,
        *gsm_cols_93,
        "mean_log2",
        "median_log2",
        "t_1samp",
        "p_1samp",
        "q_BH_within_panel",
        "direction",
        "scanarray_median_log2",
        "scanarray_direction",
        "scanarray_n",
        "mapping_note",
    ]
    write_tsv(os.path.join(TABLES, "all_panel_genes.tsv"), all_genes, gene_fields)
    for panel, fname in (
        ("AXIS", "axis_genes.tsv"),
        ("PRIORITY", "priority_genes.tsv"),
        ("APM", "apm_genes.tsv"),
        ("IFN_IMMUNE", "ifn_genes.tsv"),
        ("TJ", "tj_genes.tsv"),
    ):
        write_tsv(
            os.path.join(TABLES, fname),
            [r for r in all_genes if r["panel"] == panel],
            gene_fields,
        )

    # compact key table for the PR
    key_syms = ["CLDN4", "CLDN3", "TACSTD2"] + PRIORITY + ["B2M", "TAP1", "TAP2", "HLA-B", "HLA-C"]
    key_rows = []
    for st in ("GSE22421", "GSE22493"):
        for sym in key_syms:
            hits = [r for r in results[st]["gene_rows"] if r["symbol"] == sym]
            if not hits:
                continue
            # prefer AXIS/PRIORITY over later duplicate panels
            row = hits[0]
            key_rows.append(row)
    write_tsv(os.path.join(TABLES, "key_genes.tsv"), key_rows, gene_fields)

    set_rows = results["GSE22421"]["set_rows"] + results["GSE22493"]["set_rows"]
    write_tsv(
        os.path.join(TABLES, "geneset_stats.tsv"),
        set_rows,
        [
            "study",
            "gene_set",
            "n_in_list",
            "n_measured_ge2_arrays",
            "n_missing_or_lt2",
            "missing_or_lt2",
            "n_mean_up",
            "n_mean_down",
            "median_set_mean_log2",
            "mean_set_mean_log2",
            "median_background_mean_log2",
            "n_background_genes",
            "mannwhitney_U",
            "mannwhitney_p",
            "set_shift",
        ],
    )
    gsea_rows = results["GSE22421"]["gsea_rows"] + results["GSE22493"]["gsea_rows"]
    write_tsv(
        os.path.join(TABLES, "hallmark_gsea_preranked.tsv"),
        gsea_rows,
        ["study", "gene_set", "n_in_gmt", "n_in_set_measured", "es", "nes", "p_perm_1000", "note"],
    )
    write_tsv(
        os.path.join(TABLES, "cldn4_diagnostic.tsv"),
        results["GSE22421"]["diag_rows"] + results["GSE22493"]["diag_rows"],
        [
            "study",
            "gsm",
            "probe",
            "deposited_log2",
            "scanarray_log2",
            "ch1_cy3_median_minus_B",
            "ch2_cy5_median_minus_B",
        ],
    )

    probe_fields = [
        "study",
        "panel",
        "symbol",
        "probe",
        "mapped_via",
        "raw_symbol",
        "alias_note",
        "description",
        *gsm_cols_21,
        *gsm_cols_93,
    ]
    write_tsv(
        os.path.join(TABLES, "probe_level.tsv"),
        results["GSE22421"]["probe_rows"] + results["GSE22493"]["probe_rows"],
        [
            "study",
            "panel",
            "symbol",
            "probe",
            "mapped_via",
            "raw_symbol",
            "alias_note",
            "description",
            *gsm_cols_21,
            *[f"{g}_scan" for g in gsm_cols_21],
            *gsm_cols_93,
            *[f"{g}_scan" for g in gsm_cols_93],
        ],
    )

    # paired table
    pair_rows = []
    seen = set()
    for row in results["GSE22421"]["gene_rows"]:
        key = (row["panel"], row["symbol"])
        if key in seen:
            continue
        seen.add(key)
        o = next(
            (
                x
                for x in results["GSE22493"]["gene_rows"]
                if x["panel"] == row["panel"] and x["symbol"] == row["symbol"]
            ),
            None,
        )
        pair_rows.append(
            {
                "panel": row["panel"],
                "symbol": row["symbol"],
                "GSE22421_measured": row["measured"],
                "GSE22421_n": row["n_arrays_with_value"],
                "GSE22421_median_log2": row["median_log2"],
                "GSE22421_direction": row["direction"],
                "GSE22421_p_1samp": row["p_1samp"],
                "GSE22493_measured": o["measured"] if o else "no",
                "GSE22493_n": o["n_arrays_with_value"] if o else 0,
                "GSE22493_median_log2": o["median_log2"] if o else "",
                "GSE22493_direction": o["direction"] if o else "NA",
                "GSE22493_p_1samp": o["p_1samp"] if o else "",
                "same_sign": (
                    "yes"
                    if o
                    and row["direction"] in ("UP", "DOWN")
                    and o["direction"] == row["direction"]
                    else (
                        "no"
                        if o
                        and row["direction"] in ("UP", "DOWN")
                        and o["direction"] in ("UP", "DOWN")
                        else "NA"
                    )
                ),
                "mapping_note": row["mapping_note"],
            }
        )
    write_tsv(
        os.path.join(TABLES, "paired_gse22421_vs_gse22493.tsv"),
        pair_rows,
        [
            "panel",
            "symbol",
            "GSE22421_measured",
            "GSE22421_n",
            "GSE22421_median_log2",
            "GSE22421_direction",
            "GSE22421_p_1samp",
            "GSE22493_measured",
            "GSE22493_n",
            "GSE22493_median_log2",
            "GSE22493_direction",
            "GSE22493_p_1samp",
            "same_sign",
            "mapping_note",
        ],
    )

    make_figures(results)

    def set_get(st, name, field):
        row = next(r for r in results[st]["set_rows"] if r["gene_set"] == name)
        return row[field]

    def gene_get(st, panel, sym, field):
        row = next(r for r in results[st]["gene_rows"] if r["panel"] == panel and r["symbol"] == sym)
        return row[field]

    s21 = panel_summary(results["GSE22421"]["gene_rows"], "PRIORITY")
    a21 = panel_summary(results["GSE22421"]["gene_rows"], "APM")
    i21 = panel_summary(results["GSE22421"]["gene_rows"], "IFN_IMMUNE")
    t21 = panel_summary(results["GSE22421"]["gene_rows"], "TJ")
    s93 = panel_summary(results["GSE22493"]["gene_rows"], "PRIORITY")
    a93 = panel_summary(results["GSE22493"]["gene_rows"], "APM")

    key_stats = {
        "GSE22421": {
            "accession": "GSE22421",
            "pmid": "21123456",
            "platform": "GPL10555",
            "model": STUDIES["GSE22421"]["model"],
            "contrast": STUDIES["GSE22421"]["contrast"],
            "n_arrays": 3,
            "primary_metric": "deposited series-matrix VALUE = log2(C-CPE/untreated); gene = median of probes per array",
            "cldn4_per_array": results["GSE22421"]["cldn4_per_array"],
            "cldn4_median_log2": results["GSE22421"]["cldn4_median"],
            "cldn4_n_arrays": results["GSE22421"]["cldn4_n"],
            "tacstd2_on_platform": False,
            "priority": s21,
            "apm": a21,
            "ifn": i21,
            "tj": t21,
            "ifn_set_mannwhitney_p": set_get("GSE22421", "IFN_IMMUNE", "mannwhitney_p"),
            "ifn_set_median_log2": set_get("GSE22421", "IFN_IMMUNE", "median_set_mean_log2"),
            "apm_set_mannwhitney_p": set_get("GSE22421", "APM", "mannwhitney_p"),
            "apm_set_median_log2": set_get("GSE22421", "APM", "median_set_mean_log2"),
            "n_background_genes": results["GSE22421"]["n_background"],
        },
        "GSE22493": {
            "accession": "GSE22493",
            "platform": "GPL10555",
            "model": STUDIES["GSE22493"]["model"],
            "contrast": STUDIES["GSE22493"]["contrast"],
            "n_arrays": 3,
            "symbols_mapped": True,
            "mapping_rule": "ORF or DESCRIPTION prefix; G1P2->ISG15",
            "cldn4_per_array": results["GSE22493"]["cldn4_per_array"],
            "cldn4_median_log2": results["GSE22493"]["cldn4_median"],
            "priority": s93,
            "apm": a93,
            "ifn_set_mannwhitney_p": set_get("GSE22493", "IFN_IMMUNE", "mannwhitney_p"),
            "ifn_set_median_log2": set_get("GSE22493", "IFN_IMMUNE", "median_set_mean_log2"),
        },
        "pairing": {
            "paired": True,
            "reason": "Both series use GPL10555; gene symbols map by the same ORF/DESCRIPTION rule.",
        },
        "honesty": [
            "Not lung; SKOV-3 / SKOV-3-IP-Luc ovarian lines.",
            "Not SKB264; C-CPE is a CLDN3/4-binding C-terminal CPE fragment.",
            "GSE22421 n=3 two-color arrays; GSE22493 n=3 two-color arrays.",
            "TACSTD2 is not annotated on GPL10555.",
            "GSE22493 control is CLDN4 overexpression, not scramble/WT.",
            "One-sample t and MW/GSEA with n=3 are descriptive.",
            "Dye-swap is claimed in text; deposited GSMs all have Cy3=control, Cy5=treated.",
        ],
    }
    with open(os.path.join(OUT, "key_stats.json"), "w") as fh:
        json.dump(key_stats, fh, indent=2, default=str)
        fh.write("\n")

    print("wrote", TABLES)
    print("wrote", FIGS)
    print("GSE22421 CLDN4 median", results["GSE22421"]["cldn4_median"], "n", results["GSE22421"]["cldn4_n"])
    print("GSE22421 PRIORITY", s21)
    print("GSE22421 APM", a21)
    print("GSE22421 IFN MW p", set_get("GSE22421", "IFN_IMMUNE", "mannwhitney_p"))
    print("GSE22493 PRIORITY", s93)
    print("TACSTD2 on platform", results["GSE22421"]["tacstd2_on_platform"])


if __name__ == "__main__":
    main()
