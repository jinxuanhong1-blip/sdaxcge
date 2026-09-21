#!/usr/bin/env python3
"""GSE22493: remap GPL10555 symbols, then score NHEJ, STING, and IFN.

Primary numbers are the deposited series-matrix VALUE
(log2 CLDN4-siRNA / CLDN4-overexpression control). ScanArray Cy5/Cy3 is
a sensitivity check for the NHEJ, STING, and break-marker genes.

Symbol rule (NCBI Homo_sapiens.gene_info, protein-coding, official):
  - token is a current symbol: keep it, even if some other gene lists that
    token as a retired synonym (the GPL description is the tie-break; see
    symbol_collisions_kept.tsv for the panel cases)
  - token is not a current symbol and is an unambiguous synonym of exactly
    one current gene: remap
  - token matches several genes and none of them is the token itself:
    unresolved, not in a panel and not in the background
  - empty ORF with no SYMBOL-- prefix: unresolved. Free-text descriptions
    are not mined for gene names.

Usage:
  python3 scripts/nhej_sting_gse22493/download_data.py
  python3 scripts/nhej_sting_gse22493/run_analysis.py
"""

from __future__ import annotations

import csv
import gzip
import hashlib
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

from gene_sets import (
    BREAK_MARKERS,
    IFN,
    NHEJ,
    PERTURBATION,
    REQUIRED_REMAPS,
    STING,
    STING_OUTPUTS,
)

DATA_DIR = os.environ.get("NHEJ_STING_GSE22493_DATA", "/tmp/nhej_sting_gse22493")
REPO = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
OUT = os.path.join(REPO, "results", "nhej_sting_gse22493")
os.makedirs(OUT, exist_ok=True)

GSMS = ["GSM558700", "GSM558701", "GSM558702"]

# Free-text probes with no ORF. Reported, not folded into gene medians.
FREETEXT_PROBES = {
    "8697": "DCLRE1C",
    "25010": "XRCC5",
}


def _split_fields(line: str) -> list[str]:
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]


def md5_of(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_symbol_table(path: str) -> tuple[dict[str, str], dict[str, set[str]]]:
    """Return official-symbol lookup and token -> set of official symbols."""
    current: dict[str, str] = {}
    syn_to: dict[str, set[str]] = defaultdict(set)
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if fields[idx["type_of_gene"]] != "protein-coding":
                continue
            if fields[idx["Nomenclature_status"]] != "O":
                continue
            auth = fields[idx["Symbol_from_nomenclature_authority"]]
            if not auth or auth == "-":
                continue
            syns: list[str] = []
            if fields[idx["Synonyms"]] != "-":
                syns = [
                    s
                    for s in fields[idx["Synonyms"]].split("|")
                    if s and s != "-"
                ]
            current[auth.upper()] = auth
            for token in [auth] + syns:
                token = token.strip()
                if token:
                    syn_to[token.upper()].add(auth)
    return current, syn_to


def resolve_token(
    raw: str,
    current: dict[str, str],
    syn_to: dict[str, set[str]],
) -> tuple[str | None, str]:
    key = raw.upper()
    cands = syn_to.get(key, set())
    if key in current:
        return current[key], "current_symbol"
    if len(cands) == 1:
        return next(iter(cands)), "synonym"
    if len(cands) > 1:
        return None, "ambiguous_synonym"
    return None, "unmapped"


def load_platform(path: str, current: dict[str, str], syn_to: dict[str, set[str]]):
    ann: dict[str, dict] = {}
    counts = defaultdict(int)
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                next(fh)
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            fields = _split_fields(line)
            if not fields or not fields[0]:
                continue
            counts["probes"] += 1
            desc = fields[1] if len(fields) > 1 else ""
            orf = fields[2].strip() if len(fields) > 2 else ""
            via = "ORF"
            raw = orf
            if not raw and "--" in desc:
                raw = desc.split("--", 1)[0].strip()
                via = "DESCRIPTION_PREFIX"
            if not raw:
                counts["empty"] += 1
                ann[fields[0]] = {
                    "probe": fields[0],
                    "raw": "",
                    "via": "",
                    "symbol": None,
                    "how": "empty",
                    "description": desc,
                }
                continue
            symbol, how = resolve_token(raw, current, syn_to)
            counts[how] += 1
            ann[fields[0]] = {
                "probe": fields[0],
                "raw": raw,
                "via": via,
                "symbol": symbol,
                "how": how,
                "description": desc,
            }
    return ann, counts


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
            row: list[float | None] = []
            for cell in fields[1:]:
                row.append(None if cell == "" else float(cell))
            while len(row) < len(samples):
                row.append(None)
            values[fields[0]] = row[: len(samples)]
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
            try:
                ch1 = float(fields[cols["Ch1 Median - B"]])
                ch2 = float(fields[cols["Ch2 Median - B"]])
            except (ValueError, KeyError, IndexError):
                continue
            rec = {
                "ch1_ctrl_cy3": ch1,
                "ch2_kd_cy5": ch2,
                "log2_cy5_over_cy3": None,
            }
            if ch1 > 0 and ch2 > 0:
                rec["log2_cy5_over_cy3"] = math.log2(ch2 / ch1)
            out[fields[0]] = rec
    return out


def collapse_gene(
    probe_ids: list[str],
    matrix: dict[str, list[float | None]],
    n_samples: int,
) -> list[float | None]:
    per_array: list[float | None] = []
    for j in range(n_samples):
        vals = []
        for pid in probe_ids:
            row = matrix.get(pid)
            if row is None or j >= len(row) or row[j] is None:
                continue
            if not np.isfinite(row[j]):
                continue
            vals.append(row[j])
        per_array.append(float(np.median(vals)) if vals else None)
    return per_array


def finite(xs: list[float | None]) -> list[float]:
    return [x for x in xs if x is not None and np.isfinite(x)]


def onesamp(vals: list[float]) -> tuple[float | None, float | None]:
    if len(vals) < 2:
        return None, None
    if float(np.std(vals)) == 0.0:
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
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def panel_defs() -> list[tuple[str, list[tuple[str, str]], str]]:
    return [
        ("NHEJ", NHEJ, "down"),
        ("STING", STING, "up"),
        ("IFN", [(g, "IFN program") for g in IFN], "up"),
        ("BREAK_MARKERS", BREAK_MARKERS, "up"),
        ("PERTURBATION", [(g, "CLDN4 siRNA target") for g in PERTURBATION], "down"),
    ]


def main() -> None:
    soft = os.path.join(DATA_DIR, "GSE22493_family.soft.gz")
    matrix_path = os.path.join(DATA_DIR, "GSE22493_series_matrix.txt.gz")
    gene_info = os.path.join(DATA_DIR, "Homo_sapiens.gene_info.gz")
    for path in (soft, matrix_path, gene_info):
        if not os.path.exists(path):
            sys.exit(f"Missing {path}. Run download_data.py first.")

    current, syn_to = load_symbol_table(gene_info)
    for raw, expected in REQUIRED_REMAPS.items():
        got, how = resolve_token(raw, current, syn_to)
        if got != expected or how != "synonym":
            sys.exit(f"Remap check failed: {raw} -> {got} ({how}), expected {expected}")

    ann, counts = load_platform(soft, current, syn_to)
    samples, matrix = load_series_matrix(matrix_path)
    if samples != GSMS:
        sys.exit(f"Unexpected sample order: {samples}")

    probes_by_gene: dict[str, list[str]] = defaultdict(list)
    for pid, rec in ann.items():
        if rec["symbol"]:
            probes_by_gene[rec["symbol"]].append(pid)

    scan: dict[str, dict[str, dict]] = {}
    for gsm in GSMS:
        sp = os.path.join(DATA_DIR, "scanarray", f"{gsm}.txt.gz")
        if os.path.exists(sp):
            scan[gsm] = parse_scanarray(sp)

    bg_mean: dict[str, float] = {}
    for symbol, pids in probes_by_gene.items():
        per = collapse_gene(pids, matrix, len(samples))
        vals = finite(per)
        if len(vals) >= 2:
            bg_mean[symbol] = float(np.mean(vals))

    gene_rows: list[dict] = []
    for panel, genes, _expected in panel_defs():
        for symbol, role in genes:
            if panel == "IFN" and symbol in STING_OUTPUTS:
                role = "STING transcriptional output, scored in IFN"
            pids = sorted(probes_by_gene.get(symbol, []), key=lambda x: int(x))
            raws = sorted({ann[p]["raw"] for p in pids})
            hows = sorted({ann[p]["how"] for p in pids})
            vias = sorted({ann[p]["via"] for p in pids})
            per = collapse_gene(pids, matrix, len(samples)) if pids else [None] * 3
            vals = finite(per)
            t_stat, p_val = onesamp(vals)
            med = float(np.median(vals)) if vals else None
            mean = float(np.mean(vals)) if vals else None
            mapping = "absent"
            if hows == ["synonym"]:
                mapping = "synonym:" + ",".join(raws)
            elif hows == ["current_symbol"]:
                mapping = "current_symbol"
            elif hows:
                mapping = "+".join(hows) + ":" + ",".join(raws)
            gene_rows.append(
                {
                    "panel": panel,
                    "symbol": symbol,
                    "role": role,
                    "on_platform": "yes" if pids else "no",
                    "scored_ge2_arrays": "yes" if len(vals) >= 2 else "no",
                    "n_probes": len(pids),
                    "n_arrays": len(vals),
                    "raw_symbols": ",".join(raws),
                    "mapping": mapping,
                    "mapped_via": ",".join(vias),
                    "GSM558700": fmt(per[0]) if pids else "",
                    "GSM558701": fmt(per[1]) if pids else "",
                    "GSM558702": fmt(per[2]) if pids else "",
                    "mean_log2_KD_over_ctrl": fmt(mean),
                    "median_log2_KD_over_ctrl": fmt(med),
                    "t_1samp": fmt(t_stat),
                    "p_1samp": f"{p_val:.4g}" if p_val is not None else "",
                    "direction_median": direction_of(med) if len(vals) >= 2 else "NA",
                    "mean_num": mean,
                    "median_num": med if len(vals) >= 2 else None,
                    "p_num": p_val,
                    "per_array": per,
                }
            )

    for panel, _genes, _expected in panel_defs():
        idx = [
            i
            for i, row in enumerate(gene_rows)
            if row["panel"] == panel and row["p_num"] is not None
        ]
        qmap: dict[int, float] = {}
        if idx:
            qvals = multipletests([gene_rows[i]["p_num"] for i in idx], method="fdr_bh")[1]
            qmap = {idx[j]: float(qvals[j]) for j in range(len(idx))}
        for i, row in enumerate(gene_rows):
            if row["panel"] != panel:
                continue
            row["q_BH_within_panel"] = fmt(qmap[i], 4) if i in qmap else ""

    set_rows = []
    bg_vals = list(bg_mean.values())
    bg_median = float(np.median(bg_vals)) if bg_vals else None
    for panel, genes, expected in panel_defs():
        if panel == "PERTURBATION":
            continue
        symbols = [g for g, _role in genes]
        set_vals = [bg_mean[g] for g in symbols if g in bg_mean]
        scored_rows = [
            r
            for r in gene_rows
            if r["panel"] == panel and r["scored_ge2_arrays"] == "yes"
        ]
        n_up = sum(1 for r in scored_rows if r["direction_median"] == "UP")
        n_down = sum(1 for r in scored_rows if r["direction_median"] == "DOWN")
        n_zero = sum(1 for r in scored_rows if r["direction_median"] == "ZERO")
        missing = [g for g in symbols if g not in probes_by_gene]
        lt2 = [
            g
            for g in symbols
            if g in probes_by_gene and g not in bg_mean
        ]
        if set_vals and bg_vals:
            others = [bg_mean[g] for g in bg_mean if g not in set(symbols)]
            _u, p = stats.mannwhitneyu(set_vals, others, alternative="two-sided")
        else:
            p = None
        med_set = float(np.median(set_vals)) if set_vals else None
        on_expected = ""
        if med_set is not None:
            on_expected = (
                "yes"
                if (expected == "up" and med_set > 0)
                or (expected == "down" and med_set < 0)
                else "no"
            )
        set_rows.append(
            {
                "gene_set": panel,
                "expected_if_nhej_sting_ifn_logic": expected,
                "n_in_list": len(symbols),
                "n_on_platform": len(symbols) - len(missing),
                "n_scored_ge2_arrays": len(set_vals),
                "n_median_up": n_up,
                "n_median_down": n_down,
                "n_median_zero": n_zero,
                "median_of_mean_log2": fmt(med_set),
                "background_median_mean_log2": fmt(bg_median),
                "mannwhitney_p_two_sided": f"{p:.4g}" if p is not None else "",
                "median_on_expected_side_of_zero": on_expected,
                "absent_from_platform": ",".join(missing),
                "on_platform_but_lt2_arrays": ",".join(lt2),
            }
        )

    # Probe-level audit for NHEJ, STING, markers, CLDN4, and synonym-mapped IFN genes.
    focus_symbols = set()
    for row in gene_rows:
        if row["panel"] in {"NHEJ", "STING", "BREAK_MARKERS", "PERTURBATION"}:
            focus_symbols.add(row["symbol"])
        elif row["panel"] == "IFN" and row["mapping"].startswith("synonym:"):
            focus_symbols.add(row["symbol"])
    probe_rows = []
    for symbol in sorted(focus_symbols):
        for pid in sorted(probes_by_gene.get(symbol, []), key=lambda x: int(x)):
            rec = ann[pid]
            row = matrix.get(pid, [None, None, None])
            probe_rows.append(
                {
                    "probe": pid,
                    "symbol": symbol,
                    "raw_symbol": rec["raw"],
                    "mapping": rec["how"],
                    "mapped_via": rec["via"],
                    "description": rec["description"],
                    "GSM558700_deposited": fmt(row[0]),
                    "GSM558701_deposited": fmt(row[1]),
                    "GSM558702_deposited": fmt(row[2]),
                }
            )

    panel_symbols = {r["symbol"] for r in gene_rows}
    collision_rows = []
    for pid, rec in sorted(ann.items(), key=lambda kv: int(kv[0])):
        if rec["symbol"] not in panel_symbols or rec["how"] != "current_symbol":
            continue
        others = sorted(
            c for c in syn_to.get(rec["raw"].upper(), set()) if c != rec["symbol"]
        )
        if not others:
            continue
        collision_rows.append(
            {
                "probe": pid,
                "kept_symbol": rec["symbol"],
                "platform_token": rec["raw"],
                "also_retired_synonym_of": ",".join(others),
                "description": rec["description"],
                "decision": "kept as the current symbol; GPL description names that gene",
            }
        )

    freetext_rows = []
    for pid, symbol in FREETEXT_PROBES.items():
        rec = ann.get(pid)
        row = matrix.get(pid, [None, None, None])
        if rec is None:
            continue
        freetext_rows.append(
            {
                "probe": pid,
                "likely_gene": symbol,
                "orf": rec["raw"],
                "how": rec["how"],
                "description": rec["description"][:180],
                "used_in_gene_median": "no",
                "GSM558700": fmt(row[0]),
                "GSM558701": fmt(row[1]),
                "GSM558702": fmt(row[2]),
                "reason": "no ORF and description is not SYMBOL--prefix; excluded",
            }
        )

    sens_rows = []
    sens_symbols = [
        r["symbol"]
        for r in gene_rows
        if r["panel"] in {"NHEJ", "STING", "BREAK_MARKERS", "PERTURBATION"}
        and r["on_platform"] == "yes"
    ]
    for symbol in sens_symbols:
        pids = probes_by_gene.get(symbol, [])
        dep = next(
            r
            for r in gene_rows
            if r["symbol"] == symbol and r["panel"] != "IFN"
        )
        scan_per: list[float | None] = []
        for gsm in GSMS:
            vals = []
            for pid in pids:
                rec = scan.get(gsm, {}).get(pid)
                if rec and rec["log2_cy5_over_cy3"] is not None:
                    vals.append(rec["log2_cy5_over_cy3"])
            scan_per.append(float(np.median(vals)) if vals else None)
        scan_vals = finite(scan_per)
        scan_med = float(np.median(scan_vals)) if scan_vals else None
        dep_dir = dep["direction_median"]
        scan_dir = direction_of(scan_med) if len(scan_vals) >= 2 else "NA"
        if dep_dir == "NA" or scan_dir == "NA":
            agree = "NA"
        elif dep_dir == scan_dir:
            agree = "yes"
        else:
            agree = "no"
        sens_rows.append(
            {
                "symbol": symbol,
                "panel": dep["panel"],
                "mapping": dep["mapping"],
                "deposited_median_log2": dep["median_log2_KD_over_ctrl"],
                "deposited_direction": dep_dir,
                "scan_GSM558700": fmt(scan_per[0]),
                "scan_GSM558701": fmt(scan_per[1]),
                "scan_GSM558702": fmt(scan_per[2]),
                "scan_median_log2": fmt(scan_med) if len(scan_vals) >= 2 else "",
                "scan_direction": scan_dir,
                "sign_agrees": agree,
            }
        )

    cldn = next(r for r in gene_rows if r["symbol"] == "CLDN4")
    cldn_rows = []
    for j, gsm in enumerate(GSMS):
        s = scan.get(gsm, {}).get("17169", {})
        cldn_rows.append(
            {
                "gsm": gsm,
                "probe": "17169",
                "channel1": "Cy3 = CLDN4 overexpression (control)",
                "channel2": "Cy5 = CLDN4 lentiviral siRNA",
                "deposited_log2_KD_over_ctrl": cldn[gsm],
                "scan_ch1_median_minus_B_Cy3": fmt(s.get("ch1_ctrl_cy3"), 2),
                "scan_ch2_median_minus_B_Cy5": fmt(s.get("ch2_kd_cy5"), 2),
                "scan_log2_Cy5_over_Cy3": fmt(s.get("log2_cy5_over_cy3")),
            }
        )

    gene_fields = [
        "panel",
        "symbol",
        "role",
        "on_platform",
        "scored_ge2_arrays",
        "n_probes",
        "n_arrays",
        "raw_symbols",
        "mapping",
        "mapped_via",
        "GSM558700",
        "GSM558701",
        "GSM558702",
        "mean_log2_KD_over_ctrl",
        "median_log2_KD_over_ctrl",
        "t_1samp",
        "p_1samp",
        "q_BH_within_panel",
        "direction_median",
    ]
    write_tsv(os.path.join(OUT, "panel_genes.tsv"), gene_rows, gene_fields)
    for panel in ("NHEJ", "STING", "IFN", "BREAK_MARKERS"):
        write_tsv(
            os.path.join(OUT, f"{panel.lower()}_genes.tsv"),
            [r for r in gene_rows if r["panel"] == panel],
            gene_fields,
        )
    write_tsv(
        os.path.join(OUT, "geneset_stats.tsv"),
        set_rows,
        list(set_rows[0].keys()),
    )
    write_tsv(
        os.path.join(OUT, "probe_level_focus.tsv"),
        probe_rows,
        [
            "probe",
            "symbol",
            "raw_symbol",
            "mapping",
            "mapped_via",
            "description",
            "GSM558700_deposited",
            "GSM558701_deposited",
            "GSM558702_deposited",
        ],
    )
    write_tsv(
        os.path.join(OUT, "symbol_collisions_kept.tsv"),
        collision_rows,
        [
            "probe",
            "kept_symbol",
            "platform_token",
            "also_retired_synonym_of",
            "description",
            "decision",
        ],
    )
    write_tsv(
        os.path.join(OUT, "freetext_probes_excluded.tsv"),
        freetext_rows,
        list(freetext_rows[0].keys()) if freetext_rows else ["probe"],
    )
    write_tsv(
        os.path.join(OUT, "sensitivity_deposited_vs_scanarray.tsv"),
        sens_rows,
        list(sens_rows[0].keys()) if sens_rows else ["symbol"],
    )
    write_tsv(
        os.path.join(OUT, "cldn4_diagnostic.tsv"),
        cldn_rows,
        list(cldn_rows[0].keys()),
    )
    write_tsv(
        os.path.join(OUT, "symbol_remap_counts.tsv"),
        [
            {
                "probes": counts["probes"],
                "kept_current_symbol": counts["current_symbol"],
                "remapped_unambiguous_synonym": counts["synonym"],
                "ambiguous_synonym": counts["ambiguous_synonym"],
                "unmapped_token": counts["unmapped"],
                "empty_no_symbol_token": counts["empty"],
                "resolved_genes": len(probes_by_gene),
                "background_genes_ge2_arrays": len(bg_mean),
                "gene_info_md5": md5_of(gene_info),
                "series_matrix_md5": md5_of(matrix_path),
                "family_soft_md5": md5_of(soft),
            }
        ],
        [
            "probes",
            "kept_current_symbol",
            "remapped_unambiguous_synonym",
            "ambiguous_synonym",
            "unmapped_token",
            "empty_no_symbol_token",
            "resolved_genes",
            "background_genes_ge2_arrays",
            "gene_info_md5",
            "series_matrix_md5",
            "family_soft_md5",
        ],
    )

    def set_lookup(name: str) -> dict:
        return next(r for r in set_rows if r["gene_set"] == name)

    nhej_s = set_lookup("NHEJ")
    sting_s = set_lookup("STING")
    ifn_s = set_lookup("IFN")
    key = {
        "accession": "GSE22493",
        "platform": "GPL10555",
        "model": "SKOV-3-IP-Luc ovarian cancer cell line",
        "contrast": "CLDN4 lentiviral siRNA (Cy5) / CLDN4 overexpression control (Cy3)",
        "n_arrays": 3,
        "primary_metric": (
            "deposited series-matrix VALUE = log2(KD/control); "
            "gene = median of probes per array; set score uses the mean across arrays with n>=2"
        ),
        "symbol_remap": {
            "source": "NCBI Homo_sapiens.gene_info protein-coding official symbols and synonyms",
            "gene_info_md5": md5_of(gene_info),
            "probes_remapped_synonym": counts["synonym"],
            "probes_kept_current": counts["current_symbol"],
            "probes_ambiguous_excluded": counts["ambiguous_synonym"],
            "background_genes_ge2": len(bg_mean),
            "required_remaps_ok": True,
        },
        "cldn4_median_log2": cldn["median_log2_KD_over_ctrl"],
        "cldn4_n_arrays": cldn["n_arrays"],
        "cldn4_knockdown_confirmed_on_array": False,
        "nhej": nhej_s,
        "sting": sting_s,
        "ifn": ifn_s,
        "break_markers": set_lookup("BREAK_MARKERS"),
        "logic": (
            "Pre-specified pattern was NHEJ down, STING up, IFN up "
            "versus the overexpression control. Score what was measured. "
            "STING1 and TBK1 are not on GPL10555."
        ),
        "verdict": "panels_do_not_show_coordinated_NHEJ_down_STING_up_IFN_up",
        "honesty": [
            "Not lung; SKOV-3 ovarian line.",
            "Control is CLDN4 overexpression, not scramble or WT.",
            "CLDN4 knockdown is not cleanly confirmed on the array.",
            "n=3 two-color arrays from 2010.",
            "STING1 and TBK1 have no probe, so the STING hub is unscored.",
            "CGAS is present only as C6orf150; XRCC6 only as G22P1; APLF only as C2orf13.",
            "NHEJ1 and PAXX are absent.",
            "Ambiguous synonyms were dropped rather than assigned.",
            "Free-text Artemis and Ku80 probes were not folded into gene medians.",
            "One-sample t-tests at n=2-3 are descriptive.",
        ],
    }
    with open(os.path.join(OUT, "key_stats.json"), "w") as fh:
        json.dump(key, fh, indent=2)
        fh.write("\n")

    plot_panels(gene_rows, bg_median if bg_median is not None else 0.0)
    print(json.dumps({"nhej": nhej_s, "sting": sting_s, "ifn": ifn_s}, indent=2))
    print(f"wrote {OUT}")


def plot_panels(gene_rows: list[dict], bg_median: float) -> None:
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(11.4, 6.4),
        gridspec_kw={"width_ratios": [1.2, 1.15, 0.85]},
    )
    fig.suptitle(
        "GSE22493  log2 (CLDN4 siRNA / CLDN4-overexpression control)",
        fontsize=12,
    )

    def bars(ax, panel: str, title: str) -> None:
        rows = [
            r
            for r in gene_rows
            if r["panel"] == panel and r["scored_ge2_arrays"] == "yes"
        ]
        rows = sorted(rows, key=lambda r: r["median_num"] if r["median_num"] is not None else 0)
        labels = []
        for r in rows:
            label = r["symbol"]
            if str(r["mapping"]).startswith("synonym:"):
                label = f"{r['symbol']}  [{r['raw_symbols']}]"
            labels.append(label)
        vals = [r["median_num"] for r in rows]
        colors = ["#b85c38" if v > 0 else "#3b6d9a" if v < 0 else "#888888" for v in vals]
        y = np.arange(len(rows))
        ax.barh(y, vals, color=colors, height=0.72)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("median log2 across arrays")
        absent = [
            r["symbol"] for r in gene_rows if r["panel"] == panel and r["on_platform"] == "no"
        ]
        thin = [
            r["symbol"]
            for r in gene_rows
            if r["panel"] == panel and r["on_platform"] == "yes" and r["scored_ge2_arrays"] == "no"
        ]
        note = title
        if absent:
            note += "\nnot on array: " + ", ".join(absent)
        if thin:
            note += "\nn<2 arrays: " + ", ".join(thin)
        ax.set_title(note, fontsize=9, loc="left")

    bars(axes[0], "NHEJ", "NHEJ")
    bars(axes[1], "STING", "STING axis")

    ifn_meds = [
        r["median_num"]
        for r in gene_rows
        if r["panel"] == "IFN" and r["median_num"] is not None
    ]
    axes[2].boxplot(
        ifn_meds,
        orientation="vertical",
        widths=0.45,
        showfliers=False,
        medianprops={"color": "black"},
    )
    jitter = np.random.default_rng(0).uniform(-0.08, 0.08, size=len(ifn_meds))
    axes[2].scatter(
        1 + jitter,
        ifn_meds,
        s=12,
        c="#4d4d4d",
        alpha=0.75,
        zorder=3,
    )
    axes[2].axhline(0, color="black", lw=0.8)
    axes[2].axhline(bg_median, color="#888888", lw=0.8, ls="--")
    axes[2].set_xticks([])
    axes[2].set_ylabel("per-gene median log2")
    axes[2].set_title(
        f"IFN  n={len(ifn_meds)} scored\ndashed = background median",
        fontsize=9,
        loc="left",
    )

    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_nhej_sting_ifn.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
