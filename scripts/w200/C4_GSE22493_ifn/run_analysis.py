#!/usr/bin/env python3
"""C4 / GSE22493: CLDN4 siRNA vs CLDN4-high control — IFN / MHC-I / APM.

Primary numbers are the author-deposited series-matrix VALUE column
(GEO: "log2 ratio knockdown/control"). ScanArray Cy5/Cy3 intensities
are a sensitivity check only.

This slice does not re-analyze other CLDN4-loss datasets. It does not
invent p-values for missing genes. Direction is KD / overexpression
control, not KD / scramble WT.

Usage:
  python3 scripts/w200/C4_GSE22493_ifn/download_data.py
  python3 scripts/w200/C4_GSE22493_ifn/run_analysis.py
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

from gene_sets import APM, IFN_IMMUNE, PERTURBATION, PRIORITY

DATA_DIR = os.environ.get("W200_C4_GSE22493_DATA", "/tmp/w200_c4_gse22493_ifn")
REPO = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
OUT = os.path.join(REPO, "results", "w200", "C4_GSE22493_ifn")
os.makedirs(OUT, exist_ok=True)

GSMS = ["GSM558700", "GSM558701", "GSM558702"]
ALIASES = {
    "G1P2": "ISG15",  # historical symbol; IFI-15K / ISG15
    "IFI15": "ISG15",
    "G1P3": "IFI6",
}

# Conservative: only ORF or "SYMBOL--description". No free-text HLA salvage.
# G1P2 is remapped explicitly and labelled in the tables.


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
            # pad if a trailing empty was dropped
            while len(row) < len(samples):
                row.append(None)
            values[probe] = row[: len(samples)]
    return samples, values


def parse_scanarray(path: str) -> dict[str, dict]:
    """Return probe-index -> raw channel medians and log2(Cy5/Cy3) if both >0."""
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
            name = fields[cols["Name"]].strip().strip('"')
            try:
                ch1 = float(fields[cols["Ch1 Median - B"]])
                ch2 = float(fields[cols["Ch2 Median - B"]])
            except (ValueError, KeyError):
                continue
            rec = {
                "name": name,
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
    """Per-array median across probes. Missing if no probe has a value."""
    per_array: list[float | None] = []
    n_probes = len(probe_ids)
    for j in range(n_samples):
        vals = []
        for pid in probe_ids:
            row = matrix.get(pid)
            if row is None or j >= len(row) or row[j] is None:
                continue
            vals.append(row[j])
        per_array.append(float(np.median(vals)) if vals else None)
    return per_array, n_probes


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

    # Background: every mapped symbol with ≥1 deposited value on ≥2 arrays.
    bg_logfc: dict[str, float] = {}
    bg_n: dict[str, int] = {}
    for symbol, pids in probes_by_gene.items():
        per_array, _ = collapse_gene(pids, matrix, len(samples))
        finite = _finite(per_array)
        if len(finite) >= 2:
            bg_logfc[symbol] = float(np.mean(finite))
            bg_n[symbol] = len(finite)

    panels = {
        "PRIORITY": PRIORITY,
        "APM": APM,
        "IFN_IMMUNE": IFN_IMMUNE,
        "PERTURBATION": PERTURBATION,
    }

    gene_rows: list[dict] = []
    for panel, genes in panels.items():
        for symbol in genes:
            pids = sorted(probes_by_gene.get(symbol, []), key=lambda x: int(x))
            mapping_note = ""
            if pids:
                notes = sorted(
                    {
                        ann[p]["alias_note"]
                        for p in pids
                        if ann[p]["alias_note"]
                    }
                )
                vias = sorted({ann[p]["mapped_via"] for p in pids})
                mapping_note = ";".join(notes + vias)
            per_array, n_probes = collapse_gene(pids, matrix, len(samples))
            finite = _finite(per_array)
            t_stat, p_val = onesamp(finite)
            med = float(np.median(finite)) if finite else None
            mean = float(np.mean(finite)) if finite else None
            gene_rows.append(
                {
                    "panel": panel,
                    "symbol": symbol,
                    "measured": "yes" if finite else "no",
                    "n_probes_annotated": n_probes,
                    "n_arrays_with_value": len(finite),
                    "GSM558700": fmt(per_array[0]) if pids else "",
                    "GSM558701": fmt(per_array[1]) if pids else "",
                    "GSM558702": fmt(per_array[2]) if pids else "",
                    "mean_log2_KD_over_ctrl": fmt(mean),
                    "median_log2_KD_over_ctrl": fmt(med),
                    "t_1samp": fmt(t_stat),
                    "p_1samp": f"{p_val:.4g}" if p_val is not None else "",
                    "direction": direction_of(med),
                    "mapping_note": mapping_note,
                    "mean_log2_num": mean,
                    "median_log2_num": med,
                    "p_num": p_val,
                    "per_array": per_array,
                }
            )

    # BH-FDR within each panel among measured genes only.
    for panel in panels:
        idx = [
            i
            for i, r in enumerate(gene_rows)
            if r["panel"] == panel and r["p_num"] is not None
        ]
        if not idx:
            for r in gene_rows:
                if r["panel"] == panel:
                    r["q_BH_within_panel"] = ""
            continue
        pvals = [gene_rows[i]["p_num"] for i in idx]
        qvals = multipletests(pvals, method="fdr_bh")[1]
        qmap = {idx[j]: float(qvals[j]) for j in range(len(idx))}
        for i, r in enumerate(gene_rows):
            if r["panel"] != panel:
                continue
            r["q_BH_within_panel"] = fmt(qmap[i], 4) if i in qmap else ""

    # Probe-level table for priority + APM + CLDN4 (honesty about HLA-A).
    probe_focus = set(PRIORITY + APM + PERTURBATION)
    probe_rows: list[dict] = []
    for symbol in sorted(probe_focus):
        for pid in sorted(probes_by_gene.get(symbol, []), key=lambda x: int(x)):
            row = matrix.get(pid, [None, None, None])
            rec = ann[pid]
            scan_bits = {}
            for gsm in GSMS:
                s = scan.get(gsm, {}).get(pid, {})
                scan_bits[f"{gsm}_ch1_cy3"] = fmt(s.get("ch1_ctrl_cy3"), 2) if s else ""
                scan_bits[f"{gsm}_ch2_cy5"] = fmt(s.get("ch2_kd_cy5"), 2) if s else ""
                scan_bits[f"{gsm}_scan_log2"] = (
                    fmt(s.get("log2_cy5_over_cy3")) if s else ""
                )
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

    # CLDN4 diagnostic (deposited + ScanArray).
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

    # Gene-set tests vs background (primary deposited, gene-collapsed means).
    set_rows = []
    bg_vals = list(bg_logfc.values())
    for panel, genes in panels.items():
        if panel == "PERTURBATION":
            continue
        set_vals = [bg_logfc[g] for g in genes if g in bg_logfc]
        missing = [g for g in genes if g not in bg_logfc]
        if set_vals and bg_vals:
            # Mann-Whitney: set vs genes not in the set.
            others = [bg_logfc[g] for g in bg_logfc if g not in set(genes)]
            u, p = stats.mannwhitneyu(set_vals, others, alternative="two-sided")
        else:
            u, p = None, None
        n_up = sum(1 for v in set_vals if v > 0)
        n_down = sum(1 for v in set_vals if v < 0)
        set_rows.append(
            {
                "gene_set": panel,
                "n_in_list": len(genes),
                "n_measured_ge2_arrays": len(set_vals),
                "n_missing_or_lt2": len(missing),
                "missing_or_lt2": ",".join(missing),
                "n_mean_up": n_up,
                "n_mean_down": n_down,
                "median_set_mean_log2": fmt(float(np.median(set_vals)) if set_vals else None),
                "mean_set_mean_log2": fmt(float(np.mean(set_vals)) if set_vals else None),
                "median_background_mean_log2": fmt(float(np.median(bg_vals)) if bg_vals else None),
                "mannwhitney_U": fmt(float(u), 1) if u is not None else "",
                "mannwhitney_p": f"{p:.4g}" if p is not None else "",
                "set_shift": (
                    "DOWN"
                    if set_vals and np.median(set_vals) < np.median(bg_vals)
                    else "UP"
                    if set_vals and np.median(set_vals) > np.median(bg_vals)
                    else "NA"
                ),
            }
        )

    # Sample table.
    sample_rows = [
        {
            "gsm": "GSM558700",
            "title": "CLDN4 Replicate 1",
            "cell_line": "SKOV-3-IP-Luc",
            "ch1_Cy3": "CLDN4 overexpression (control)",
            "ch2_Cy5": "CLDN4 knockdown (lentiviral siRNA)",
            "deposited_value": "log2 ratio knockdown/control",
        },
        {
            "gsm": "GSM558701",
            "title": "CLDN4 Replicate 2",
            "cell_line": "SKOV-3-IP-Luc",
            "ch1_Cy3": "CLDN4 overexpression (control)",
            "ch2_Cy5": "CLDN4 knockdown (lentiviral siRNA)",
            "deposited_value": "log2 ratio knockdown/control",
        },
        {
            "gsm": "GSM558702",
            "title": "CLDN4 Replicate 3",
            "cell_line": "SKOV-3-IP-Luc",
            "ch1_Cy3": "CLDN4 overexpression (control)",
            "ch2_Cy5": "CLDN4 knockdown (lentiviral siRNA)",
            "deposited_value": "log2 ratio knockdown/control",
        },
    ]

    # Verdict helpers.
    pri = [r for r in gene_rows if r["panel"] == "PRIORITY"]
    apm = [r for r in gene_rows if r["panel"] == "APM"]
    pri_meas = [r for r in pri if r["measured"] == "yes"]
    apm_meas = [r for r in apm if r["measured"] == "yes"]
    pri_up = [r for r in pri_meas if r["direction"] == "UP"]
    apm_up = [r for r in apm_meas if r["direction"] == "UP"]
    cldn4 = next(r for r in gene_rows if r["symbol"] == "CLDN4")
    cldn4_finite = _finite(cldn4["per_array"])

    key = {
        "accession": "GSE22493",
        "platform": "GPL10555",
        "model": "SKOV-3-IP-Luc ovarian cancer cell line",
        "contrast": "CLDN4 lentiviral siRNA (Cy5) / CLDN4 overexpression control (Cy3)",
        "n_arrays": 3,
        "primary_metric": "deposited series-matrix VALUE = log2(KD/control); gene = median of probes per array, then mean/median across arrays",
        "cldn4_deposited_per_array": [fmt(x) for x in cldn4["per_array"]],
        "cldn4_median_log2": cldn4["median_log2_KD_over_ctrl"],
        "cldn4_n_arrays": cldn4["n_arrays_with_value"],
        "cldn4_knockdown_confirmed_on_array": False,
        "cldn4_note": (
            "Probe 17169 is missing on GSM558700; the two remaining deposited "
            "log-ratios are negative (KD < overexpression control) but the "
            "ScanArray Cy5 background-subtracted median is non-positive on "
            "GSM558700 and opposite-signed on GSM558701 vs GSM558702. "
            "Control is CLDN4 overexpression, not scramble WT."
        ),
        "priority_measured": len(pri_meas),
        "priority_up": len(pri_up),
        "priority_down": sum(1 for r in pri_meas if r["direction"] == "DOWN"),
        "priority_genes_up": [r["symbol"] for r in pri_up],
        "priority_genes_down": [
            r["symbol"] for r in pri_meas if r["direction"] == "DOWN"
        ],
        "apm_measured": len(apm_meas),
        "apm_up": len(apm_up),
        "apm_down": sum(1 for r in apm_meas if r["direction"] == "DOWN"),
        "ifn_set_mannwhitney_p": next(
            r["mannwhitney_p"] for r in set_rows if r["gene_set"] == "IFN_IMMUNE"
        ),
        "ifn_set_median_log2": next(
            r["median_set_mean_log2"] for r in set_rows if r["gene_set"] == "IFN_IMMUNE"
        ),
        "n_background_genes": len(bg_logfc),
        "verdict": (
            "does_not_support_CLDN4_loss_opening_IFN_MHCI_APM"
        ),
        "honesty": [
            "Not lung; SKOV-3 ovarian line.",
            "Control is CLDN4 overexpression, not WT/scramble.",
            "CLDN4 knockdown is not cleanly confirmed on the array.",
            "n=3 two-color arrays; dye-swap is claimed in text but all three GSMs have the same Cy3/Cy5 assignment.",
            "Old custom Operon 60-mer array; HLA-A has 11 discordant probes.",
            "ISG15 is present only as historical symbol G1P2.",
            "NLRC5, ERAP1, ERAP2, PDIA3 are not annotated on GPL10555.",
            "One-sample t-tests with n=2–3 are descriptive, not confirmatory.",
            "This accession does not reproduce an IFN/MHC-I/APM up-signature after CLDN4 siRNA.",
        ],
    }

    # Write tables.
    gene_fields = [
        "panel",
        "symbol",
        "measured",
        "n_probes_annotated",
        "n_arrays_with_value",
        "GSM558700",
        "GSM558701",
        "GSM558702",
        "mean_log2_KD_over_ctrl",
        "median_log2_KD_over_ctrl",
        "t_1samp",
        "p_1samp",
        "q_BH_within_panel",
        "direction",
        "mapping_note",
    ]
    write_tsv(
        os.path.join(OUT, "priority_genes.tsv"),
        [r for r in gene_rows if r["panel"] == "PRIORITY"],
        gene_fields,
    )
    write_tsv(
        os.path.join(OUT, "apm_genes.tsv"),
        [r for r in gene_rows if r["panel"] == "APM"],
        gene_fields,
    )
    write_tsv(
        os.path.join(OUT, "ifn_genes.tsv"),
        [r for r in gene_rows if r["panel"] == "IFN_IMMUNE"],
        gene_fields,
    )
    write_tsv(
        os.path.join(OUT, "probe_level.tsv"),
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
            "GSM558700_ch1_cy3",
            "GSM558700_ch2_cy5",
            "GSM558700_scan_log2",
            "GSM558701_ch1_cy3",
            "GSM558701_ch2_cy5",
            "GSM558701_scan_log2",
            "GSM558702_ch1_cy3",
            "GSM558702_ch2_cy5",
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
            "n_in_list",
            "n_measured_ge2_arrays",
            "n_missing_or_lt2",
            "missing_or_lt2",
            "n_mean_up",
            "n_mean_down",
            "median_set_mean_log2",
            "mean_set_mean_log2",
            "median_background_mean_log2",
            "mannwhitney_U",
            "mannwhitney_p",
            "set_shift",
        ],
    )
    write_tsv(
        os.path.join(OUT, "sample_table.tsv"),
        sample_rows,
        [
            "gsm",
            "title",
            "cell_line",
            "ch1_Cy3",
            "ch2_Cy5",
            "deposited_value",
        ],
    )
    # Sensitivity: ScanArray raw log2(Cy5/Cy3) vs deposited VALUE.
    # Only spots with both background-subtracted medians > 0.
    sens_rows = []
    for symbol in PRIORITY + ["CLDN4"]:
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
        sens_rows.append(
            {
                "symbol": symbol,
                "deposited_median": fmt(float(np.median(dep_f)) if dep_f else None),
                "deposited_direction": direction_of(
                    float(np.median(dep_f)) if dep_f else None
                ),
                "scanarray_median": fmt(float(np.median(sc_f)) if sc_f else None),
                "scanarray_direction": direction_of(
                    float(np.median(sc_f)) if sc_f else None
                ),
                "GSM558700_deposited": fmt(dep_per[0]),
                "GSM558701_deposited": fmt(dep_per[1]),
                "GSM558702_deposited": fmt(dep_per[2]),
                "GSM558700_scan": fmt(scan_per[0]),
                "GSM558701_scan": fmt(scan_per[1]),
                "GSM558702_scan": fmt(scan_per[2]),
                "sign_agreement": (
                    "yes"
                    if dep_f
                    and sc_f
                    and np.sign(np.median(dep_f)) == np.sign(np.median(sc_f))
                    and np.median(dep_f) != 0
                    else "no_or_NA"
                ),
            }
        )
    write_tsv(
        os.path.join(OUT, "sensitivity_deposited_vs_scanarray.tsv"),
        sens_rows,
        [
            "symbol",
            "deposited_median",
            "deposited_direction",
            "scanarray_median",
            "scanarray_direction",
            "GSM558700_deposited",
            "GSM558701_deposited",
            "GSM558702_deposited",
            "GSM558700_scan",
            "GSM558701_scan",
            "GSM558702_scan",
            "sign_agreement",
        ],
    )
    key["priority_scanarray_vs_deposited"] = {
        r["symbol"]: {
            "deposited": r["deposited_direction"],
            "scanarray": r["scanarray_direction"],
            "agree": r["sign_agreement"],
        }
        for r in sens_rows
    }

    with open(os.path.join(OUT, "key_stats.json"), "w") as fh:
        json.dump(key, fh, indent=2)
        fh.write("\n")

    # Figure: priority + APM medians with per-array points.
    plot_rows = [
        r
        for r in gene_rows
        if r["panel"] in ("PRIORITY", "APM") and r["symbol"] != "HLA-A" or r["panel"] == "PRIORITY"
    ]
    # unique: priority first, then remaining APM
    seen = set()
    ordered = []
    for r in gene_rows:
        if r["panel"] == "PRIORITY":
            ordered.append(r)
            seen.add(r["symbol"])
    for r in gene_rows:
        if r["panel"] == "APM" and r["symbol"] not in seen:
            ordered.append(r)
            seen.add(r["symbol"])

    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    y = np.arange(len(ordered))[::-1]
    meds = [
        r["median_log2_num"] if r["median_log2_num"] is not None else np.nan
        for r in ordered
    ]
    colors = []
    for r in ordered:
        if r["panel"] == "PRIORITY":
            colors.append("#1f4e79" if r["direction"] == "UP" else "#9b1b30" if r["direction"] == "DOWN" else "#888888")
        else:
            colors.append("#4c78a8" if r["direction"] == "UP" else "#c44e52" if r["direction"] == "DOWN" else "#888888")
    ax.axvline(0, color="0.4", lw=0.8)
    ax.scatter(meds, y, c=colors, s=42, zorder=3, edgecolors="white", linewidths=0.4)
    for yi, r in zip(y, ordered):
        xs = _finite(r["per_array"])
        ax.scatter(xs, [yi] * len(xs), c="0.45", s=14, alpha=0.7, zorder=2)
    labels = []
    for r in ordered:
        tag = "P" if r["panel"] == "PRIORITY" else "A"
        miss = "" if r["measured"] == "yes" else " (absent)"
        labels.append(f"{r['symbol']} [{tag}]{miss}")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("deposited log2 (CLDN4 siRNA / CLDN4-high control)")
    ax.set_title(
        "GSE22493 SKOV-3: IFN priority [P] and MHC-I/APM [A]\n"
        "median of probes/array (large); arrays as small points"
    )
    ax.set_xlim(-4.2, 2.6)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_priority_apm_log2.png"), dpi=140)
    plt.close(fig)

    # Print a short console summary for the writeup author.
    print("OUT", OUT)
    print("CLDN4 deposited", [fmt(x) for x in cldn4["per_array"]], "median", cldn4["median_log2_KD_over_ctrl"])
    print("PRIORITY")
    for r in pri:
        print(
            f"  {r['symbol']:8s} n={r['n_arrays_with_value']}  "
            f"median={r['median_log2_KD_over_ctrl'] or 'NA':8s}  "
            f"{r['direction']:5s}  p={r['p_1samp'] or 'NA'}"
        )
    print("APM up", [r["symbol"] for r in apm_up], "down", [r["symbol"] for r in apm if r["direction"] == "DOWN"])
    for r in set_rows:
        print(
            f"{r['gene_set']:12s} meas={r['n_measured_ge2_arrays']}  "
            f"med={r['median_set_mean_log2']}  p={r['mannwhitney_p']}  {r['set_shift']}"
        )
    print("verdict", key["verdict"])


if __name__ == "__main__":
    main()
