#!/usr/bin/env python3
"""One method sweep for GSE22493. Does not add genes or cohorts.

Methods were listed before reading their panel p-values. Each one changes
a single processing choice relative to the primary analysis:

  deposited_probe_median_gene_mean   primary (probe median, then mean of arrays)
  deposited_probe_mean_gene_mean     mean of probes instead of median
  deposited_probe_median_gene_median median across arrays instead of mean
  deposited_lowest_probe_id          one probe per gene (smallest probe id)
  scanarray_log2                     log2 of background-subtracted Cy5/Cy3
  vendor_ch2_logratio_abs_le_15      ScanArray "Ch2 Log Ratio", |value| > 15 dropped
  deposited_plus_freetext_ku_artemis primary plus empty-ORF probes 25010 and 8697

A separate array-level score (mean of the panel on each array, then a
one-sample t across the three arrays) uses the primary per-array values.
PRE_VALUE is checked only as a scale QC: it is the linear ratio and VALUE
is log2(PRE_VALUE).

Usage:
  python3 scripts/w200/GSE22493_nhej_sting_ifn_apm/method_sweep.py
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import sys

import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import run_analysis as ra
from gene_sets import APM, IFN, NHEJ, NHEJ_ACCESSORY, PRIMARY, SECONDARY, STING, STING_REG

PANELS = {}
ROLES = {}
for name, genes in PRIMARY.items():
    PANELS[name] = list(genes)
    ROLES[name] = "primary"
for name, genes in SECONDARY.items():
    PANELS[name] = list(genes)
    ROLES[name] = "secondary"
PANELS["CLDN4"] = ["CLDN4"]
ROLES["CLDN4"] = "diagnostic"

PRIMARY_NAMES = list(PRIMARY)


def _finite(xs):
    return [x for x in xs if x is not None and np.isfinite(x)]


def collapse(pids, matrix, how):
    per = []
    for j in range(3):
        vals = []
        for pid in pids:
            row = matrix.get(pid)
            if row is None or row[j] is None or not np.isfinite(row[j]):
                continue
            vals.append(row[j])
        if not vals:
            per.append(None)
        elif how == "median":
            per.append(float(np.median(vals)))
        elif how == "mean":
            per.append(float(np.mean(vals)))
        else:
            raise ValueError(how)
    return per


def load_vendor_logratio():
    """ScanArray Ch2 Log Ratio. Drop |value| > 15 (broken spots; max was ~1e4)."""
    out = {gsm: {} for gsm in ra.GSMS}
    for gsm in ra.GSMS:
        path = os.path.join(ra.DATA_DIR, "scanarray", f"{gsm}.txt.gz")
        with gzip.open(path, "rt", errors="replace") as fh:
            in_data = False
            header = None
            for line in fh:
                if line.startswith("BEGIN DATA"):
                    in_data = True
                    continue
                if not in_data:
                    continue
                if header is None:
                    header = line.rstrip("\n").split("\t")
                    continue
                fields = line.rstrip("\n").split("\t")
                if not fields or not fields[0].isdigit():
                    continue
                cols = {name: i for i, name in enumerate(header)}
                try:
                    val = float(fields[cols["Ch2 Log Ratio"]])
                except (ValueError, KeyError):
                    continue
                if not np.isfinite(val) or abs(val) > 15:
                    continue
                out[gsm][fields[0]] = val
    matrix = {}
    probes = set()
    for gsm in ra.GSMS:
        probes.update(out[gsm])
    for pid in probes:
        matrix[pid] = [out[gsm].get(pid) for gsm in ra.GSMS]
    return matrix


def load_scanarray_log2():
    matrix = {}
    parsed = {}
    for gsm in ra.GSMS:
        path = os.path.join(ra.DATA_DIR, "scanarray", f"{gsm}.txt.gz")
        parsed[gsm] = ra.parse_scanarray(path)
    probes = set()
    for gsm in ra.GSMS:
        probes.update(parsed[gsm])
    for pid in probes:
        row = []
        for gsm in ra.GSMS:
            rec = parsed[gsm].get(pid)
            row.append(None if rec is None else rec.get("log2_cy5_over_cy3"))
        matrix[pid] = row
    return matrix


def scale_qc(soft_path):
    """Confirm VALUE == log2(PRE_VALUE) inside each sample table."""
    rows = []
    with gzip.open(soft_path, "rt", errors="replace") as fh:
        gsm = None
        mode = None
        header = None
        for line in fh:
            if line.startswith("^SAMPLE"):
                gsm = line.split("=", 1)[1].strip()
            if line.startswith("!sample_table_begin"):
                mode = "table"
                header = next(fh).rstrip("\n").split("\t")
                iv = header.index("VALUE")
                ip = header.index("PRE_VALUE")
                abs_err = []
                n = 0
                continue
            if mode == "table":
                if line.startswith("!sample_table_end"):
                    err = np.array(abs_err) if abs_err else np.array([np.nan])
                    rows.append(
                        {
                            "gsm": gsm,
                            "n_both_finite_pre_positive": n,
                            "median_abs_VALUE_minus_log2_PRE": f"{float(np.median(err)):.6g}",
                            "max_abs_VALUE_minus_log2_PRE": f"{float(np.max(err)):.6g}",
                            "value_equals_log2_pre": "yes"
                            if n and float(np.max(err)) < 1e-3
                            else "no",
                        }
                    )
                    mode = None
                    continue
                fields = line.rstrip("\n").split("\t")
                v = fields[iv].strip()
                p = fields[ip].strip()
                if v == "" or p == "":
                    continue
                pv = float(p)
                if pv <= 0:
                    continue
                abs_err.append(abs(float(v) - math.log2(pv)))
                n += 1
    return rows


def gene_scores(probes_by_gene, matrix, how_probe, array_summary):
    """symbol -> score used in the set test, plus per-array values."""
    scores = {}
    per_array = {}
    for symbol, pids in probes_by_gene.items():
        arr = collapse(pids, matrix, how_probe)
        finite = _finite(arr)
        per_array[symbol] = arr
        if len(finite) < 2:
            continue
        if array_summary == "mean":
            scores[symbol] = float(np.mean(finite))
        else:
            scores[symbol] = float(np.median(finite))
    return scores, per_array


def panel_stats(scores, per_array, genes):
    set_vals = [scores[g] for g in genes if g in scores]
    missing = [g for g in genes if g not in scores]
    bg = [v for g, v in scores.items() if g not in set(genes)]
    if set_vals and bg:
        _u, p = stats.mannwhitneyu(set_vals, bg, alternative="two-sided")
        p = float(p)
    else:
        p = None
    # descriptive gene-level t and BH, same rule as the primary script
    pvals = []
    p_idx = []
    any_q = False
    min_p = None
    for i, g in enumerate(genes):
        finite = _finite(per_array.get(g, []))
        if len(finite) < 2:
            continue
        t, pv = stats.ttest_1samp(finite, 0.0)
        if np.isfinite(pv):
            pvals.append(float(pv))
            p_idx.append(g)
            min_p = float(pv) if min_p is None else min(min_p, float(pv))
    genes_q = []
    if pvals:
        q = multipletests(pvals, method="fdr_bh")[1]
        any_q = bool(np.any(q < 0.05))
        min_q = float(np.min(q))
        genes_q = [p_idx[i] for i in range(len(p_idx)) if q[i] < 0.05]
    else:
        min_q = None
    set_med = float(np.median(set_vals)) if set_vals else None
    bg_med = float(np.median(list(scores.values()))) if scores else None
    if set_med is None or bg_med is None:
        shift = "NA"
        delta = None
    else:
        delta = set_med - bg_med
        if abs(delta) < 0.05:
            shift = "FLAT"
        elif delta < 0:
            shift = "DOWN"
        else:
            shift = "UP"
    return {
        "n_in_list": len(genes),
        "n_ge2": len(set_vals),
        "n_missing_or_lt2": len(missing),
        "missing_or_lt2": ",".join(missing),
        "n_up": sum(1 for v in set_vals if v > 0),
        "n_down": sum(1 for v in set_vals if v < 0),
        "median_score": set_med,
        "background_median": bg_med,
        "delta_vs_background": delta,
        "mannwhitney_p": p,
        "min_gene_p": min_p,
        "min_gene_q": min_q,
        "any_gene_q_lt_0.05": any_q,
        "genes_q_lt_0.05": genes_q,
        "shift": shift,
        "scores": {g: scores[g] for g in genes if g in scores},
    }


def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return ""
    if isinstance(x, float):
        return f"{x:.{nd}g}" if abs(x) >= 1000 else f"{x:.{nd}f}"
    return str(x)


def array_level(per_array, genes, global_median):
    """Mean of genes present on each array, then t across 3 arrays.

    Also tests the panel mean after subtracting that array's transcriptome-wide
    median, so a whole-array shift is not counted as a panel effect.
    """
    array_means = []
    residual = []
    for j in range(3):
        vals = []
        for g in genes:
            arr = per_array.get(g)
            if not arr or arr[j] is None or not np.isfinite(arr[j]):
                continue
            vals.append(arr[j])
        if vals:
            array_means.append(float(np.mean(vals)))
            residual.append(float(np.mean(vals)) - global_median[j])
        else:
            array_means.append(None)
            residual.append(None)
    def _tp(xs):
        finite = _finite(xs)
        if len(finite) < 2:
            return None
        _t, p = stats.ttest_1samp(finite, 0.0)
        return float(p) if np.isfinite(p) else None
    return array_means, _tp(array_means), residual, _tp(residual)


def main():
    soft = os.path.join(ra.DATA_DIR, "GSE22493_family.soft.gz")
    matrix_path = os.path.join(ra.DATA_DIR, "GSE22493_series_matrix.txt.gz")
    if not os.path.exists(soft) or not os.path.exists(matrix_path):
        sys.exit("Missing GEO files. Run download_data.py first.")

    ann = ra.load_soft_annotations(soft)
    _samples, deposited = ra.load_series_matrix(matrix_path)
    probes_by_gene = {}
    for pid, rec in ann.items():
        probes_by_gene.setdefault(rec["symbol"], []).append(pid)

    # Lowest probe id only.
    lowest = {g: [min(pids, key=lambda x: int(x))] for g, pids in probes_by_gene.items()}

    # Primary map plus the two empty-ORF spots previously excluded.
    plus = {g: list(pids) for g, pids in probes_by_gene.items()}
    plus.setdefault("XRCC5", []).append("25010")
    plus.setdefault("DCLRE1C", []).append("8697")

    scan = load_scanarray_log2()
    vendor = load_vendor_logratio()
    qc_rows = scale_qc(soft)

    methods = [
        ("deposited_probe_median_gene_mean", probes_by_gene, deposited, "median", "mean"),
        ("deposited_probe_mean_gene_mean", probes_by_gene, deposited, "mean", "mean"),
        ("deposited_probe_median_gene_median", probes_by_gene, deposited, "median", "median"),
        ("deposited_lowest_probe_id", lowest, deposited, "median", "mean"),
        ("scanarray_log2", probes_by_gene, scan, "median", "mean"),
        ("vendor_ch2_logratio_abs_le_15", probes_by_gene, vendor, "median", "mean"),
        ("deposited_plus_freetext_ku_artemis", plus, deposited, "median", "mean"),
    ]

    summary_rows = []
    gene_rows = []
    stored = {}
    for method, pmap, matrix, how_probe, array_summary in methods:
        scores, per_array = gene_scores(pmap, matrix, how_probe, array_summary)
        stored[method] = {}
        for panel, genes in PANELS.items():
            st = panel_stats(scores, per_array, genes)
            stored[method][panel] = st
            summary_rows.append(
                {
                    "method": method,
                    "panel": panel,
                    "role": ROLES[panel],
                    "n_in_list": st["n_in_list"],
                    "n_ge2": st["n_ge2"],
                    "n_up": st["n_up"],
                    "n_down": st["n_down"],
                    "median_score": fmt(st["median_score"]),
                    "background_median": fmt(st["background_median"]),
                    "delta_vs_background": fmt(st["delta_vs_background"]),
                    "mannwhitney_p": f"{st['mannwhitney_p']:.4g}" if st["mannwhitney_p"] is not None else "",
                    "min_gene_p": f"{st['min_gene_p']:.4g}" if st["min_gene_p"] is not None else "",
                    "min_gene_q": f"{st['min_gene_q']:.4g}" if st["min_gene_q"] is not None else "",
                    "any_gene_q_lt_0.05": "yes" if st["any_gene_q_lt_0.05"] else "no",
                    "shift": st["shift"],
                    "missing_or_lt2": st["missing_or_lt2"],
                }
            )
            if panel in ("NHEJ", "STING", "CLDN4"):
                for g in genes:
                    gene_rows.append(
                        {
                            "method": method,
                            "panel": panel,
                            "symbol": g,
                            "score": fmt(st["scores"].get(g)),
                            "in_set_test": "yes" if g in st["scores"] else "no",
                        }
                    )

    # Array-level t on the primary per-array matrix (probe median of deposited).
    _scores, per_array = gene_scores(probes_by_gene, deposited, "median", "mean")
    global_median = []
    for j in range(3):
        vals = [arr[j] for arr in per_array.values() if arr[j] is not None and np.isfinite(arr[j])]
        global_median.append(float(np.median(vals)))
    array_rows = []
    for panel, genes in PANELS.items():
        if panel == "CLDN4":
            continue
        means, p, residual, p_resid = array_level(per_array, genes, global_median)
        finite = _finite(means)
        array_rows.append(
            {
                "panel": panel,
                "role": ROLES[panel],
                "GSM558700_panel_mean": fmt(means[0]),
                "GSM558701_panel_mean": fmt(means[1]),
                "GSM558702_panel_mean": fmt(means[2]),
                "median_of_array_means": fmt(float(np.median(finite)) if finite else None),
                "t_1samp_p_n3": f"{p:.4g}" if p is not None else "",
                "GSM558700_minus_array_median": fmt(residual[0]),
                "GSM558701_minus_array_median": fmt(residual[1]),
                "GSM558702_minus_array_median": fmt(residual[2]),
                "t_p_after_array_median": f"{p_resid:.4g}" if p_resid is not None else "",
                "n_arrays": len(finite),
                "array_median_700": fmt(global_median[0]),
                "array_median_701": fmt(global_median[1]),
                "array_median_702": fmt(global_median[2]),
            }
        )

    # Cross-method sign of the set median vs the primary method.
    base = "deposited_probe_median_gene_mean"
    discord_rows = []
    for method in stored:
        if method == base:
            continue
        for panel in PRIMARY_NAMES:
            a = stored[base][panel]["median_score"]
            b = stored[method][panel]["median_score"]
            if a is None or b is None or a == 0 or b == 0:
                agree = "no_or_NA"
            elif np.sign(a) == np.sign(b):
                agree = "yes"
            else:
                agree = "no"
            discord_rows.append(
                {
                    "panel": panel,
                    "method": method,
                    "primary_median": fmt(a),
                    "method_median": fmt(b),
                    "sign_agrees_with_primary": agree,
                    "method_mw_p": f"{stored[method][panel]['mannwhitney_p']:.4g}"
                    if stored[method][panel]["mannwhitney_p"] is not None
                    else "",
                    "method_any_q_lt_0.05": "yes" if stored[method][panel]["any_gene_q_lt_0.05"] else "no",
                }
            )

    primary_hit = []
    for method, panels in stored.items():
        for panel in PRIMARY_NAMES:
            st = panels[panel]
            if st["any_gene_q_lt_0.05"] or (
                st["mannwhitney_p"] is not None and st["mannwhitney_p"] < 0.05
            ):
                primary_hit.append(f"{method}:{panel}")

    # Sign agreement of set medians between deposited primary and scanarray.
    scan_name = "scanarray_log2"
    scan_flips = []
    for panel in PRIMARY_NAMES:
        a = stored[base][panel]["median_score"]
        b = stored[scan_name][panel]["median_score"]
        if a is None or b is None or np.sign(a) != np.sign(b):
            scan_flips.append(panel)

    array_hit = [
        r["panel"]
        for r in array_rows
        if r["role"] == "primary" and r["t_1samp_p_n3"] not in ("",) and float(r["t_1samp_p_n3"]) < 0.05
    ]
    deposited_hits = [h for h in primary_hit if h.startswith("deposited_probe_median_gene_mean:")]
    # Stop rule for this accession: the primary deposited test and the
    # three-array replicate test are the ones that have to agree. Nominal
    # p-values that appear only after switching the summary, and that flip
    # or vanish at the array level, are discordance.
    final_discordant = (not deposited_hits) and (not array_hit) and bool(scan_flips or primary_hit)
    final_null = (not primary_hit) and (not array_hit)
    if final_null:
        status = "FINAL_NULL"
    elif final_discordant:
        status = "FINAL_DISCORDANT"
    else:
        status = "NOT_FINAL"

    out_dir = ra.OUT
    ra.write_tsv(
        os.path.join(out_dir, "method_sweep.tsv"),
        summary_rows,
        [
            "method",
            "panel",
            "role",
            "n_in_list",
            "n_ge2",
            "n_up",
            "n_down",
            "median_score",
            "background_median",
            "delta_vs_background",
            "mannwhitney_p",
            "min_gene_p",
            "min_gene_q",
            "any_gene_q_lt_0.05",
            "shift",
            "missing_or_lt2",
        ],
    )
    ra.write_tsv(
        os.path.join(out_dir, "method_sweep_signcheck.tsv"),
        discord_rows,
        [
            "panel",
            "method",
            "primary_median",
            "method_median",
            "sign_agrees_with_primary",
            "method_mw_p",
            "method_any_q_lt_0.05",
        ],
    )
    ra.write_tsv(
        os.path.join(out_dir, "method_sweep_array_level.tsv"),
        array_rows,
        [
            "panel",
            "role",
            "GSM558700_panel_mean",
            "GSM558701_panel_mean",
            "GSM558702_panel_mean",
            "median_of_array_means",
            "t_1samp_p_n3",
            "GSM558700_minus_array_median",
            "GSM558701_minus_array_median",
            "GSM558702_minus_array_median",
            "t_p_after_array_median",
            "n_arrays",
            "array_median_700",
            "array_median_701",
            "array_median_702",
        ],
    )
    ra.write_tsv(
        os.path.join(out_dir, "method_sweep_nhej_sting_scores.tsv"),
        gene_rows,
        ["method", "panel", "symbol", "score", "in_set_test"],
    )
    ra.write_tsv(
        os.path.join(out_dir, "value_equals_log2_pre.tsv"),
        qc_rows,
        [
            "gsm",
            "n_both_finite_pre_positive",
            "median_abs_VALUE_minus_log2_PRE",
            "max_abs_VALUE_minus_log2_PRE",
            "value_equals_log2_pre",
        ],
    )

    payload = {
        "accession": "GSE22493",
        "status": status,
        "rule": (
            "FINAL_NULL if every method and the three-array t stay at p>=0.05 "
            "and q>=0.05 for NHEJ, STING, IFN, and APM. FINAL_DISCORDANT if the "
            "primary deposited method and the three-array t are both null, while "
            "some other summary produces p<0.05 or q<0.05 or the ScanArray set "
            "median flips sign versus deposited. NOT_FINAL only if the primary "
            "deposited method or the three-array t itself is p<0.05 or q<0.05."
        ),
        "primary_panel_hits_p_or_q_lt_0.05": primary_hit,
        "genes_with_q_lt_0.05": {
            method: {
                panel: stored[method][panel]["genes_q_lt_0.05"]
                for panel in stored[method]
                if stored[method][panel]["genes_q_lt_0.05"]
            }
            for method in stored
            if any(stored[method][panel]["genes_q_lt_0.05"] for panel in stored[method])
        },
        "array_level_primary_p_lt_0.05": array_hit,
        "scanarray_set_median_sign_flips_vs_deposited": scan_flips,
        "value_is_log2_of_PRE_VALUE": all(r["value_equals_log2_pre"] == "yes" for r in qc_rows),
        "methods": [m[0] for m in methods],
    }
    with open(os.path.join(out_dir, "method_sweep.json"), "w") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")

    # Stamp the primary key_stats file. Leave the primary numbers untouched.
    key_path = os.path.join(out_dir, "key_stats.json")
    with open(key_path) as fh:
        key = json.load(fh)
    key["accession_status"] = status
    key["method_sweep"] = payload
    with open(key_path, "w") as fh:
        json.dump(key, fh, indent=2)
        fh.write("\n")

    print("STATUS", status)
    print("hits", primary_hit or "-", "array_hits", array_hit or "-")
    print("scan flips", scan_flips)
    print("VALUE=log2(PRE)", payload["value_is_log2_of_PRE_VALUE"])
    for row in summary_rows:
        if row["role"] != "primary":
            continue
        print(
            f"{row['method']:42s} {row['panel']:6s} "
            f"n={row['n_ge2']:2} med={row['median_score']:8s} "
            f"p={row['mannwhitney_p']:8s} qmin={row['min_gene_q']:8s} "
            f"{row['shift']:5s} q05={row['any_gene_q_lt_0.05']}"
        )
    print("--- array level ---")
    for row in array_rows:
        print(
            f"{row['panel']:16s} {row['GSM558700_panel_mean']:8s} "
            f"{row['GSM558701_panel_mean']:8s} {row['GSM558702_panel_mean']:8s} "
            f"p={row['t_1samp_p_n3']}  resid_p={row['t_p_after_array_median']}"
        )


if __name__ == "__main__":
    main()
