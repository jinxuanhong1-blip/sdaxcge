#!/usr/bin/env python3
"""GSE22493 SKOV-3: CLDN4 siRNA vs CLDN4-overexpression — NHEJ, STING, IFN, APM.

Primary numbers are the author-deposited series-matrix VALUE. The numeric
distribution is symmetric about 0 (median about -0.04), and this slice reads
VALUE as log2(knockdown/control), the same reading as the earlier C4 slice
of this accession. GEO text says "normalized sample to control ratios" and
does not print the word log2.

ScanArray Cy5/Cy3 background-subtracted medians are a sensitivity check only.
This slice does not re-analyze other CLDN4-loss datasets and does not invent
values for genes absent from GPL10555.

Usage:
  python3 scripts/w200/GSE22493_nhej_sting_ifn_apm/download_data.py
  python3 scripts/w200/GSE22493_nhej_sting_ifn_apm/run_analysis.py
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

from gene_sets import APM, IFN, NHEJ, NHEJ_ACCESSORY, PERTURBATION, PRIMARY, SECONDARY, STING, STING_REG

DATA_DIR = os.environ.get("W200_GSE22493_PANEL_DATA", "/tmp/gse22493")
REPO = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
OUT = os.path.join(REPO, "results", "w200", "GSE22493_nhej_sting_ifn_apm")
os.makedirs(OUT, exist_ok=True)

GSMS = ["GSM558700", "GSM558701", "GSM558702"]

# Explicit historical symbols only. Near-misses are NOT mapped:
# TTBK1 is tau-tubulin kinase, not TBK1; KUB3 is not Ku70; DNTTIP1 is not
# DNTT; WRNIP1 is not WRN.
ALIASES = {
    "G22P1": "XRCC6",  # Ku70 thyroid autoantigen, ORF on this Operon array
    "C6ORF150": "CGAS",  # historical cGAS / MB21D1
    "G1P2": "ISG15",
    "IFI15": "ISG15",
    "G1P3": "IFI6",
}

# Free-text platform rows with an empty ORF. Recorded, not used.
EXCLUDED_FREETEXT = [
    {
        "probe": "8697",
        "why_not_used": "empty ORF; description names Artemis / DNA cross-link repair 1C",
        "would_have_been": "DCLRE1C",
    },
    {
        "probe": "25010",
        "why_not_used": "empty ORF; description names Ku80 / Ku86 / ATP-dependent DNA helicase II 80 kDa",
        "would_have_been": "XRCC5",
    },
]


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
                "ch1_ctrl_cy3": ch1,
                "ch2_kd_cy5": ch2,
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
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})


def bh_within(gene_rows: list[dict], panel: str) -> None:
    idx = [
        i
        for i, r in enumerate(gene_rows)
        if r["panel"] == panel and r["p_num"] is not None
    ]
    qmap: dict[int, float] = {}
    if idx:
        pvals = [gene_rows[i]["p_num"] for i in idx]
        qvals = multipletests(pvals, method="fdr_bh")[1]
        qmap = {idx[j]: float(qvals[j]) for j in range(len(idx))}
    for i, r in enumerate(gene_rows):
        if r["panel"] != panel:
            continue
        r["q_BH_within_panel"] = fmt(qmap[i], 4) if i in qmap else ""


def main() -> None:
    soft = os.path.join(DATA_DIR, "GSE22493_family.soft.gz")
    matrix_path = os.path.join(DATA_DIR, "GSE22493_series_matrix.txt.gz")
    if not os.path.exists(soft) or not os.path.exists(matrix_path):
        sys.exit("Missing GEO files. Run download_data.py first.")

    ann = load_soft_annotations(soft)
    samples, matrix = load_series_matrix(matrix_path)
    if samples != GSMS:
        raise SystemExit(f"Unexpected sample order: {samples}")

    probes_by_gene: dict[str, list[str]] = {}
    for pid, rec in ann.items():
        probes_by_gene.setdefault(rec["symbol"], []).append(pid)

    scan: dict[str, dict[str, dict]] = {}
    for gsm in GSMS:
        sp = os.path.join(DATA_DIR, "scanarray", f"{gsm}.txt.gz")
        if os.path.exists(sp):
            scan[gsm] = parse_scanarray(sp)
    if len(scan) != 3:
        raise SystemExit("ScanArray files missing. Run download_data.py first.")

    bg_logfc: dict[str, float] = {}
    for symbol, pids in probes_by_gene.items():
        per_array, _ = collapse_gene(pids, matrix, len(samples))
        finite = _finite(per_array)
        if len(finite) >= 2:
            bg_logfc[symbol] = float(np.mean(finite))

    panels: dict[str, list[str]] = {}
    roles: dict[str, str] = {}
    for name, genes in PRIMARY.items():
        panels[name] = genes
        roles[name] = "primary"
    for name, genes in SECONDARY.items():
        panels[name] = genes
        roles[name] = "secondary"
    panels["PERTURBATION"] = PERTURBATION
    roles["PERTURBATION"] = "diagnostic"

    gene_rows: list[dict] = []
    for panel, genes in panels.items():
        for symbol in genes:
            pids = sorted(probes_by_gene.get(symbol, []), key=lambda x: int(x))
            notes = []
            if pids:
                notes = sorted({ann[p]["alias_note"] for p in pids if ann[p]["alias_note"]})
                vias = sorted({ann[p]["mapped_via"] for p in pids})
                notes = notes + vias
            per_array, n_probes = collapse_gene(pids, matrix, len(samples))
            finite = _finite(per_array)
            t_stat, p_val = onesamp(finite)
            med = float(np.median(finite)) if finite else None
            mean = float(np.mean(finite)) if finite else None
            gene_rows.append(
                {
                    "panel": panel,
                    "role": roles[panel],
                    "symbol": symbol,
                    "measured": "yes" if finite else "no",
                    "n_probes_annotated": n_probes,
                    "n_arrays_with_value": len(finite),
                    "probe_ids": ",".join(pids),
                    "GSM558700": fmt(per_array[0]) if pids else "",
                    "GSM558701": fmt(per_array[1]) if pids else "",
                    "GSM558702": fmt(per_array[2]) if pids else "",
                    "mean_log2_KD_over_ctrl": fmt(mean),
                    "median_log2_KD_over_ctrl": fmt(med),
                    "t_1samp": fmt(t_stat),
                    "p_1samp": f"{p_val:.4g}" if p_val is not None else "",
                    "direction_median": direction_of(med),
                    "direction_mean": direction_of(mean),
                    "mapping_note": ";".join(notes),
                    "mean_log2_num": mean,
                    "median_log2_num": med,
                    "p_num": p_val,
                    "per_array": per_array,
                }
            )

    for panel in panels:
        bh_within(gene_rows, panel)

    # Probe-level table for the new panels plus CLDN4.
    probe_focus = set(NHEJ + NHEJ_ACCESSORY + STING + STING_REG + PERTURBATION)
    probe_rows: list[dict] = []
    for symbol in sorted(probe_focus):
        for pid in sorted(probes_by_gene.get(symbol, []), key=lambda x: int(x)):
            row = matrix.get(pid, [None, None, None])
            rec = ann[pid]
            scan_bits = {}
            for gsm in GSMS:
                s = scan.get(gsm, {}).get(pid, {})
                scan_bits[f"{gsm}_scan_log2"] = fmt(s.get("log2_cy5_over_cy3")) if s else ""
            probe_rows.append(
                {
                    "probe": pid,
                    "symbol": symbol,
                    "raw_symbol": rec["raw_symbol"],
                    "mapped_via": rec["mapped_via"],
                    "alias_note": rec["alias_note"],
                    "description": rec["description"],
                    "GSM558700_deposited": fmt(row[0]),
                    "GSM558701_deposited": fmt(row[1]),
                    "GSM558702_deposited": fmt(row[2]),
                    **scan_bits,
                }
            )

    cldn4_pid = "17169"
    cldn4_dep = matrix.get(cldn4_pid, [None, None, None])
    cldn4_rows = []
    for j, gsm in enumerate(GSMS):
        s = scan.get(gsm, {}).get(cldn4_pid, {})
        cldn4_rows.append(
            {
                "gsm": gsm,
                "array_title": f"CLDN4 Replicate {j + 1}",
                "channel1": "Cy3 = CLDN4 overexpression (control)",
                "channel2": "Cy5 = CLDN4 lentiviral siRNA",
                "deposited_log2_KD_over_ctrl": fmt(cldn4_dep[j]),
                "scan_ch1_median_minus_B_Cy3": fmt(s.get("ch1_ctrl_cy3"), 2),
                "scan_ch2_median_minus_B_Cy5": fmt(s.get("ch2_kd_cy5"), 2),
                "scan_log2_Cy5_over_Cy3": fmt(s.get("log2_cy5_over_cy3")),
                "scan_usable": "yes"
                if s.get("log2_cy5_over_cy3") is not None
                else "no_nonpositive_BG_subtracted",
            }
        )

    bg_vals = list(bg_logfc.values())
    bg_med = float(np.median(bg_vals))
    set_rows = []
    for panel, genes in panels.items():
        if panel == "PERTURBATION":
            continue
        set_vals = [bg_logfc[g] for g in genes if g in bg_logfc]
        missing = [g for g in genes if g not in bg_logfc]
        if set_vals:
            others = [bg_logfc[g] for g in bg_logfc if g not in set(genes)]
            u, p = stats.mannwhitneyu(set_vals, others, alternative="two-sided")
            # Descriptive signed-rank of gene means vs 0. Needs variation.
            if len(set_vals) >= 6 and np.unique(np.sign(set_vals)).size > 1:
                wstat, wp = stats.wilcoxon(set_vals, alternative="two-sided", zero_method="wilcox")
            else:
                wstat, wp = None, None
        else:
            u, p, wstat, wp = None, None, None, None
        n_up = sum(1 for v in set_vals if v > 0)
        n_down = sum(1 for v in set_vals if v < 0)
        set_med = float(np.median(set_vals)) if set_vals else None
        # A gap under 0.05 log2 is the same neighborhood as the background
        # median. Do not label that gap UP or DOWN.
        if set_med is None:
            shift = "NA"
            delta = None
        else:
            delta = float(set_med - bg_med)
            if abs(delta) < 0.05:
                shift = "FLAT"
            elif delta < 0:
                shift = "DOWN"
            else:
                shift = "UP"
        set_rows.append(
            {
                "gene_set": panel,
                "role": roles[panel],
                "n_in_list": len(genes),
                "n_measured_ge2_arrays": len(set_vals),
                "n_missing_or_lt2": len(missing),
                "missing_or_lt2": ",".join(missing),
                "n_mean_up": n_up,
                "n_mean_down": n_down,
                "median_set_mean_log2": fmt(set_med),
                "mean_set_mean_log2": fmt(float(np.mean(set_vals)) if set_vals else None),
                "median_background_mean_log2": fmt(bg_med),
                "delta_vs_background": fmt(delta),
                "mannwhitney_U": fmt(float(u), 1) if u is not None else "",
                "mannwhitney_p": f"{p:.4g}" if p is not None else "",
                "wilcoxon_vs0_p": f"{wp:.4g}" if wp is not None else "",
                "set_shift_vs_background": shift,
                "set_med_num": set_med,
                "mw_p_num": float(p) if p is not None else None,
            }
        )

    # Sensitivity: deposited vs ScanArray for primary-panel genes that exist,
    # plus CLDN4. Gene-level ScanArray is the median of usable probes.
    sens_symbols = []
    for name in ("NHEJ", "NHEJ_ACCESSORY", "STING", "STING_REG", "APM"):
        sens_symbols.extend(panels[name])
    # IFN is long; sensitivity is restricted to genes that are also in the
    # short axis lists or are canonical ISGs already used as the C4 priority.
    ifn_sens = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "STAT1", "IRF7", "IFNB1"]
    sens_symbols.extend(ifn_sens)
    sens_symbols.append("CLDN4")
    # unique, keep order
    seen_s = set()
    sens_symbols = [g for g in sens_symbols if not (g in seen_s or seen_s.add(g))]

    sens_rows = []
    for symbol in sens_symbols:
        pids = sorted(probes_by_gene.get(symbol, []), key=lambda x: int(x))
        dep_per, _ = collapse_gene(pids, matrix, len(samples))
        scan_per: list[float | None] = []
        for gsm in GSMS:
            vals = []
            for pid in pids:
                rec = scan.get(gsm, {}).get(pid, {})
                lr = rec.get("log2_cy5_over_cy3")
                if lr is not None:
                    vals.append(lr)
            scan_per.append(float(np.median(vals)) if vals else None)
        dep_f = _finite(dep_per)
        sc_f = _finite(scan_per)
        dep_med = float(np.median(dep_f)) if dep_f else None
        sc_med = float(np.median(sc_f)) if sc_f else None
        if dep_med is None or sc_med is None or dep_med == 0 or sc_med == 0:
            agree = "no_or_NA"
        elif np.sign(dep_med) == np.sign(sc_med):
            agree = "yes"
        else:
            agree = "no"
        memberships = [r["panel"] for r in gene_rows if r["symbol"] == symbol]
        sens_rows.append(
            {
                "symbol": symbol,
                "panels": ",".join(memberships),
                "deposited_median": fmt(dep_med),
                "deposited_direction": direction_of(dep_med),
                "scanarray_median": fmt(sc_med),
                "scanarray_direction": direction_of(sc_med),
                "GSM558700_deposited": fmt(dep_per[0]) if pids else "",
                "GSM558701_deposited": fmt(dep_per[1]) if pids else "",
                "GSM558702_deposited": fmt(dep_per[2]) if pids else "",
                "GSM558700_scan": fmt(scan_per[0]),
                "GSM558701_scan": fmt(scan_per[1]),
                "GSM558702_scan": fmt(scan_per[2]),
                "sign_agreement": agree,
            }
        )

    def panel_rows(name: str) -> list[dict]:
        return [r for r in gene_rows if r["panel"] == name]

    def count_dir(rows: list[dict], direction: str) -> int:
        return sum(1 for r in rows if r["measured"] == "yes" and r["direction_median"] == direction)

    cldn4 = next(r for r in gene_rows if r["symbol"] == "CLDN4" and r["panel"] == "PERTURBATION")

    def set_lookup(name: str) -> dict:
        return next(r for r in set_rows if r["gene_set"] == name)

    key = {
        "accession": "GSE22493",
        "platform": "GPL10555",
        "model": "SKOV-3-IP-Luc ovarian cancer cell line",
        "contrast": "CLDN4 lentiviral siRNA (Cy5) / CLDN4 overexpression control (Cy3)",
        "n_arrays": 3,
        "primary_metric": (
            "deposited series-matrix VALUE read as log2(KD/control); "
            "gene = median of ORF-mapped probes per array; "
            "set test uses the mean of arrays with a value, genes with >=2 arrays"
        ),
        "n_background_genes": len(bg_logfc),
        "background_median_mean_log2": fmt(bg_med),
        "cldn4_deposited_per_array": [fmt(x) for x in cldn4["per_array"]],
        "cldn4_median_log2": cldn4["median_log2_KD_over_ctrl"],
        "cldn4_n_arrays": cldn4["n_arrays_with_value"],
        "cldn4_knockdown_confirmed_on_array": False,
        "aliases_applied": ALIASES,
        "excluded_freetext_probes": EXCLUDED_FREETEXT,
        "panels": {},
        "honesty": [
            "Ovarian SKOV-3-IP-Luc, not lung.",
            "Control is CLDN4 overexpression, not scramble or WT.",
            "CLDN4 knockdown is not cleanly confirmed on the array.",
            "n=3 two-color Operon arrays from 2010. GEO text mentions dye-swap, but all three GSMs label Cy3 as overexpression and Cy5 as knockdown.",
            "STING1/TMEM173 and TBK1 are not annotated on GPL10555, so the STING axis is incomplete.",
            "NHEJ1/XLF and PAXX are not annotated on GPL10555.",
            "XRCC6 is present only as historical ORF G22P1. CGAS is present only as C6orf150.",
            "TTBK1, KUB3, DNTTIP1, and WRNIP1 were not remapped onto TBK1, XRCC6, DNTT, or WRN.",
            "Empty-ORF free-text Ku80 (probe 25010) and Artemis (probe 8697) were not used.",
            "One-sample t tests at n=2-3 and set Wilcoxon tests are descriptive.",
            "IFN and APM lists match the earlier C4 slice of this accession. This file does not replace that slice and does not claim a lung or ICI result.",
        ],
    }
    for name in list(PRIMARY) + list(SECONDARY):
        rows = panel_rows(name)
        measured = [r for r in rows if r["measured"] == "yes"]
        st = set_lookup(name)
        key["panels"][name] = {
            "role": roles[name],
            "n_in_list": len(rows),
            "n_measured": len(measured),
            "n_median_up": count_dir(rows, "UP"),
            "n_median_down": count_dir(rows, "DOWN"),
            "genes_median_up": [r["symbol"] for r in measured if r["direction_median"] == "UP"],
            "genes_median_down": [r["symbol"] for r in measured if r["direction_median"] == "DOWN"],
            "genes_absent": [r["symbol"] for r in rows if r["measured"] == "no"],
            "n_ge2_arrays": int(st["n_measured_ge2_arrays"]),
            "n_mean_up": int(st["n_mean_up"]),
            "n_mean_down": int(st["n_mean_down"]),
            "median_of_gene_means": st["median_set_mean_log2"],
            "delta_vs_background": st["delta_vs_background"],
            "mannwhitney_p_vs_background": st["mannwhitney_p"],
            "wilcoxon_vs0_p": st["wilcoxon_vs0_p"],
            "shift_vs_background": st["set_shift_vs_background"],
            "any_within_panel_q_lt_0.05": any(
                r["q_BH_within_panel"] not in ("", None) and float(r["q_BH_within_panel"]) < 0.05
                for r in measured
            ),
        }

    gene_fields = [
        "panel",
        "role",
        "symbol",
        "measured",
        "n_probes_annotated",
        "n_arrays_with_value",
        "probe_ids",
        "GSM558700",
        "GSM558701",
        "GSM558702",
        "mean_log2_KD_over_ctrl",
        "median_log2_KD_over_ctrl",
        "t_1samp",
        "p_1samp",
        "q_BH_within_panel",
        "direction_median",
        "direction_mean",
        "mapping_note",
    ]
    for name, filename in (
        ("NHEJ", "nhej_genes.tsv"),
        ("NHEJ_ACCESSORY", "nhej_accessory_genes.tsv"),
        ("STING", "sting_genes.tsv"),
        ("STING_REG", "sting_reg_genes.tsv"),
        ("IFN", "ifn_genes.tsv"),
        ("APM", "apm_genes.tsv"),
    ):
        write_tsv(os.path.join(OUT, filename), panel_rows(name), gene_fields)

    write_tsv(
        os.path.join(OUT, "probe_level_nhej_sting.tsv"),
        probe_rows,
        [
            "probe",
            "symbol",
            "raw_symbol",
            "mapped_via",
            "alias_note",
            "description",
            "GSM558700_deposited",
            "GSM558701_deposited",
            "GSM558702_deposited",
            "GSM558700_scan_log2",
            "GSM558701_scan_log2",
            "GSM558702_scan_log2",
        ],
    )
    write_tsv(
        os.path.join(OUT, "cldn4_diagnostic.tsv"),
        cldn4_rows,
        [
            "gsm",
            "array_title",
            "channel1",
            "channel2",
            "deposited_log2_KD_over_ctrl",
            "scan_ch1_median_minus_B_Cy3",
            "scan_ch2_median_minus_B_Cy5",
            "scan_log2_Cy5_over_Cy3",
            "scan_usable",
        ],
    )
    write_tsv(
        os.path.join(OUT, "geneset_stats.tsv"),
        set_rows,
        [
            "gene_set",
            "role",
            "n_in_list",
            "n_measured_ge2_arrays",
            "n_missing_or_lt2",
            "missing_or_lt2",
            "n_mean_up",
            "n_mean_down",
            "median_set_mean_log2",
            "mean_set_mean_log2",
            "median_background_mean_log2",
            "delta_vs_background",
            "mannwhitney_U",
            "mannwhitney_p",
            "wilcoxon_vs0_p",
            "set_shift_vs_background",
        ],
    )
    write_tsv(
        os.path.join(OUT, "sensitivity_deposited_vs_scanarray.tsv"),
        sens_rows,
        [
            "symbol",
            "panels",
            "deposited_median",
            "deposited_direction",
            "scanarray_median",
            "scanarray_direction",
            "sign_agreement",
            "GSM558700_deposited",
            "GSM558701_deposited",
            "GSM558702_deposited",
            "GSM558700_scan",
            "GSM558701_scan",
            "GSM558702_scan",
        ],
    )
    write_tsv(
        os.path.join(OUT, "excluded_freetext_probes.tsv"),
        EXCLUDED_FREETEXT,
        ["probe", "would_have_been", "why_not_used"],
    )
    write_tsv(
        os.path.join(OUT, "sample_table.tsv"),
        [
            {
                "gsm": gsm,
                "title": f"CLDN4 Replicate {i + 1}",
                "cell_line": "SKOV-3-IP-Luc",
                "tissue": "ovarian cancer cell line",
                "ch1_Cy3": "CLDN4 overexpression (control)",
                "ch2_Cy5": "CLDN4 knockdown (lentiviral siRNA)",
                "deposited_value": "series-matrix VALUE, read as log2(knockdown/control)",
            }
            for i, gsm in enumerate(GSMS)
        ],
        [
            "gsm",
            "title",
            "cell_line",
            "tissue",
            "ch1_Cy3",
            "ch2_Cy5",
            "deposited_value",
        ],
    )

    key["verdict"] = (
        "null_to_discordant_on_this_array: no within-panel q<0.05 in NHEJ, "
        "STING, IFN, or APM; STING1 and TBK1 are absent; CLDN4 knockdown "
        "is not confirmed; control is overexpression"
    )
    key["scanarray_sign_agreement"] = {}
    for name in ("NHEJ", "STING", "NHEJ_ACCESSORY", "STING_REG", "APM"):
        rows = [r for r in sens_rows if name in r["panels"].split(",")]
        key["scanarray_sign_agreement"][name] = {
            "n_with_both": sum(1 for r in rows if r["sign_agreement"] in ("yes", "no")),
            "n_agree": sum(1 for r in rows if r["sign_agreement"] == "yes"),
            "n_disagree": sum(1 for r in rows if r["sign_agreement"] == "no"),
            "agree_down": [
                r["symbol"]
                for r in rows
                if r["sign_agreement"] == "yes" and r["deposited_direction"] == "DOWN"
            ],
            "agree_up": [
                r["symbol"]
                for r in rows
                if r["sign_agreement"] == "yes" and r["deposited_direction"] == "UP"
            ],
            "disagree": [r["symbol"] for r in rows if r["sign_agreement"] == "no"],
        }

    with open(os.path.join(OUT, "key_stats.json"), "w") as fh:
        json.dump(key, fh, indent=2)
        fh.write("\n")

    _figure(gene_rows, set_rows, bg_med)
    print("OUT", OUT)
    print("background", len(bg_logfc), "median", fmt(bg_med))
    print("CLDN4", key["cldn4_deposited_per_array"], "median", key["cldn4_median_log2"])
    for name in list(PRIMARY) + list(SECONDARY):
        block = key["panels"][name]
        print(
            f"{name:16s} ge2={block['n_ge2_arrays']}/{block['n_in_list']}  "
            f"meanUp={block['n_mean_up']} meanDown={block['n_mean_down']}  "
            f"setMed={block['median_of_gene_means']}  "
            f"MW={block['mannwhitney_p_vs_background']}  {block['shift_vs_background']}  "
            f"absent={','.join(block['genes_absent']) or '-'}"
        )


def _figure(gene_rows: list[dict], set_rows: list[dict], bg_med: float) -> None:
    order_panels = ["NHEJ", "NHEJ_ACCESSORY", "STING", "STING_REG", "APM"]
    ordered = []
    for panel in order_panels:
        for r in gene_rows:
            if r["panel"] == panel:
                ordered.append(r)

    fig, axes = plt.subplots(
        1, 2, figsize=(11.2, 8.4), gridspec_kw={"width_ratios": [1.35, 1.0]}
    )
    ax = axes[0]
    y = np.arange(len(ordered))[::-1]
    ax.axvline(0, color="0.35", lw=0.8)
    ax.axvline(bg_med, color="0.65", lw=0.7, ls="--")
    for yi, r in zip(y, ordered):
        xs = _finite(r["per_array"])
        if xs:
            ax.scatter(xs, [yi] * len(xs), c="0.55", s=12, zorder=2)
        if r["median_log2_num"] is not None:
            color = "#1f4e79" if r["direction_median"] == "UP" else "#9b1b30"
            ax.scatter(
                [r["median_log2_num"]],
                [yi],
                c=color,
                s=36,
                zorder=3,
                edgecolors="white",
                linewidths=0.4,
            )
    labels = []
    for r in ordered:
        tag = {"NHEJ": "N", "NHEJ_ACCESSORY": "n", "STING": "S", "STING_REG": "s", "APM": "A"}[
            r["panel"]
        ]
        miss = "" if r["measured"] == "yes" else "  absent"
        labels.append(f"{r['symbol']} [{tag}]{miss}")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("deposited log2 (CLDN4 siRNA / CLDN4-high control)")
    ax.set_title("Gene medians\nN NHEJ  n accessory  S STING axis  s regulators  A APM")
    finite_x = [r["median_log2_num"] for r in ordered if r["median_log2_num"] is not None]
    finite_x += [x for r in ordered for x in _finite(r["per_array"])]
    pad = 0.4
    ax.set_xlim(min(finite_x) - pad, max(finite_x) + pad)

    ax2 = axes[1]
    show = [r for r in set_rows if r["gene_set"] in ("NHEJ", "STING", "IFN", "APM", "NHEJ_ACCESSORY", "STING_REG")]
    # stable display order
    want = ["NHEJ", "STING", "IFN", "APM", "NHEJ_ACCESSORY", "STING_REG"]
    show = sorted(show, key=lambda r: want.index(r["gene_set"]))
    y2 = np.arange(len(show))[::-1]
    ax2.axvline(0, color="0.35", lw=0.8)
    ax2.axvline(bg_med, color="0.65", lw=0.7, ls="--")
    for yi, r in zip(y2, show):
        val = r["set_med_num"]
        if val is None:
            continue
        color = "#1f4e79" if val > 0 else "#9b1b30"
        ax2.barh(yi, val, color=color, height=0.55, alpha=0.9)
        ax2.text(
            val + (0.04 if val >= 0 else -0.04),
            yi,
            f"p={r['mannwhitney_p']}",
            va="center",
            ha="left" if val >= 0 else "right",
            fontsize=7,
            color="0.2",
        )
    ax2.set_yticks(y2)
    ax2.set_yticklabels(
        [
            f"{r['gene_set']}  {r['n_measured_ge2_arrays']}/{r['n_in_list']}"
            for r in show
        ],
        fontsize=8,
    )
    ax2.set_xlabel("median of gene-mean log2")
    ax2.set_title("Set vs background\nMann–Whitney p; dashed = background median")
    fig.suptitle(
        "GSE22493 SKOV-3 ovarian CLDN4 siRNA vs CLDN4 overexpression (n=3)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_nhej_sting_ifn_apm.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
