#!/usr/bin/env python3
"""GSE50927 baseline (no VILI): NHEJ, Sting1/Cgas, IFN/chemokine.

Public non-cancer analog. Author EdgeR table only
(GSE50927_Cldn4lungWTvsKOgenes.csv). logFC is KO minus WT, locked by
Cldn4 logFC ≈ −6.06. Honest n is 1 vs 1 deposited GSM. VILI contrasts
are not scored in this wave.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import random
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, wilcoxon

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
DATA = Path(os.environ.get("GSE50927_NHEJ_STING_DATA", "/tmp/gse50927_nhej_sting"))
URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/"
    "GSE50927_Cldn4lungWTvsKOgenes.csv.gz"
)
SOURCE_FILE = "GSE50927_Cldn4lungWTvsKOgenes.csv.gz"
SEED = 42
N_PERM = 10000
# CPM = 1. Below this, the gene is reported and plotted, and kept out of means.
MIN_LOGCPM = 0.0

# Canonical c-NHEJ machinery. Current alias in parentheses when the mm9
# symbol in this 2014 table is the old one. Paxx was 1110057K04Rik.
NHEJ_CORE = [
    ("Xrcc6", "Ku70", "c-NHEJ"),
    ("Xrcc5", "Ku80", "c-NHEJ"),
    ("Prkdc", "DNA-PKcs", "c-NHEJ"),
    ("Dclre1c", "Artemis", "c-NHEJ"),
    ("Nhej1", "XLF", "c-NHEJ"),
    ("Xrcc4", "XRCC4", "c-NHEJ"),
    ("Lig4", "LIG4", "c-NHEJ"),
    ("Poll", "Pol-lambda", "c-NHEJ"),
    ("Polm", "Pol-mu", "c-NHEJ"),
    ("1110057K04Rik", "Paxx", "c-NHEJ"),
]

# End-recognition / 53BP1-shieldin neighborhood. Not folded into the core mean.
NHEJ_EXTENDED = [
    ("Trp53bp1", "53BP1", "NHEJ-extended"),
    ("Rnf8", "RNF8", "NHEJ-extended"),
    ("Rnf168", "RNF168", "NHEJ-extended"),
    ("Rif1", "RIF1", "NHEJ-extended"),
    ("Mad2l2", "REV7", "NHEJ-extended"),
    ("Atm", "ATM", "NHEJ-extended"),
    ("H2afx", "H2AX", "NHEJ-extended"),
    ("Mre11a", "MRE11", "NHEJ-extended"),
    ("Rad50", "RAD50", "NHEJ-extended"),
    ("Nbn", "NBS1", "NHEJ-extended"),
    ("Parp1", "PARP1", "NHEJ-extended"),
    ("Pnkp", "PNKP", "NHEJ-extended"),
    ("Aplf", "APLF", "NHEJ-extended"),
    ("Aptx", "APTX", "NHEJ-extended"),
]

# Asked pair. mm9 symbols: Cgas = Mb21d1, Sting1 = Tmem173.
STING_PAIR = [
    ("Mb21d1", "Cgas", "STING-pair"),
    ("Tmem173", "Sting1", "STING-pair"),
]

# Immediate pathway neighbors. Scored separately so they do not dilute the pair.
STING_NEIGHBOR = [
    ("Tbk1", "TBK1", "STING-neighbor"),
    ("Irf3", "IRF3", "STING-neighbor"),
    ("Ifnb1", "IFNB1", "STING-neighbor"),
    ("Zbp1", "ZBP1", "STING-neighbor"),
    ("Aim2", "AIM2", "STING-neighbor"),
    ("Trex1", "TREX1", "STING-neighbor"),
]

# Pre-specified IFN/chemokine output. Includes genes that move both ways
# (Ifit2, Ifit3) so the mean is not the FDR-only subset.
IFN_CHEMOKINE = [
    ("Ccl5", "CCL5", "IFN-chemokine"),
    ("Cxcl9", "CXCL9", "IFN-chemokine"),
    ("Cxcl10", "CXCL10", "IFN-chemokine"),
    ("Cxcl11", "CXCL11", "IFN-chemokine"),
    ("Ifnb1", "IFNB1", "IFN-chemokine"),
    ("Ifng", "IFNG", "IFN-chemokine"),
    ("Isg15", "ISG15", "IFN-chemokine"),
    ("Ifit1", "IFIT1", "IFN-chemokine"),
    ("Ifit2", "IFIT2", "IFN-chemokine"),
    ("Ifit3", "IFIT3", "IFN-chemokine"),
    ("Rsad2", "RSAD2", "IFN-chemokine"),
    ("Stat1", "STAT1", "IFN-chemokine"),
    ("Irf7", "IRF7", "IFN-chemokine"),
    ("Oasl2", "OASL2", "IFN-chemokine"),
    ("Mx1", "MX1", "IFN-chemokine"),
    ("Gbp4", "GBP4", "IFN-chemokine"),
]

# Other junction genes. Control that the table is not globally up in the KO.
TJ_CONTROL = [
    ("Cldn3", "CLDN3", "TJ-control"),
    ("Cldn5", "CLDN5", "TJ-control"),
    ("Cldn7", "CLDN7", "TJ-control"),
    ("Cldn18", "CLDN18", "TJ-control"),
    ("Tjp1", "TJP1", "TJ-control"),
    ("Ocln", "OCLN", "TJ-control"),
]

ANCHOR = ("Cldn4", "CLDN4", "anchor")


def download() -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    dest = DATA / SOURCE_FILE
    if not dest.exists() or dest.stat().st_size < 1000:
        urllib.request.urlretrieve(URL, dest)
    return dest


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_table(path: Path) -> list[dict]:
    rows = []
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = ["EntrezID", "logFC", "logCPM", "PValue", "FDR", "Marker.Symbol", "Marker.Name"]
        if reader.fieldnames != expected:
            raise SystemExit(f"unexpected header: {reader.fieldnames}")
        for row in reader:
            rows.append(
                {
                    "entrez": row["EntrezID"],
                    "symbol": row["Marker.Symbol"],
                    "name": row["Marker.Name"],
                    "logFC": float(row["logFC"]),
                    "logCPM": float(row["logCPM"]),
                    "PValue": float(row["PValue"]),
                    "FDR": float(row["FDR"]),
                }
            )
    if len(rows) < 20000:
        raise SystemExit(f"table too short: {len(rows)}")
    return rows


def index_symbols(rows: list[dict]) -> dict[str, dict]:
    """First row wins for the two Excel-mangled March symbols (1-Mar, 2-Mar).

    Those symbols are not in any scored panel. A duplicate inside a panel
    still aborts later, because lookup would be ambiguous.
    """
    out = {}
    dup = []
    for row in rows:
        sym = row["symbol"]
        if sym in out:
            dup.append(sym)
            continue
        out[sym] = row
    allowed = {"1-Mar", "2-Mar"}
    extra = sorted(set(dup) - allowed)
    if extra:
        raise SystemExit(f"unexpected duplicate symbols: {extra}")
    return out


def expressed(row: dict) -> bool:
    return row["logCPM"] >= MIN_LOGCPM


def lookup(panel: list[tuple[str, str, str]], by_sym: dict[str, dict]) -> list[dict]:
    found = []
    for symbol, alias, module in panel:
        row = by_sym.get(symbol)
        if row is None:
            found.append(
                {
                    "module": module,
                    "symbol": symbol,
                    "alias": alias,
                    "entrez": "",
                    "name": "",
                    "logFC": math.nan,
                    "logCPM": math.nan,
                    "PValue": math.nan,
                    "FDR": math.nan,
                    "in_table": "no",
                    "expressed": "no",
                    "used_in_mean": "no",
                    "direction": "missing",
                    "fdr_lt_0.05": "no",
                }
            )
            continue
        use = expressed(row)
        direction = "up" if row["logFC"] > 0 else ("down" if row["logFC"] < 0 else "zero")
        found.append(
            {
                "module": module,
                "symbol": symbol,
                "alias": alias,
                "entrez": row["entrez"],
                "name": row["name"],
                "logFC": row["logFC"],
                "logCPM": row["logCPM"],
                "PValue": row["PValue"],
                "FDR": row["FDR"],
                "in_table": "yes",
                "expressed": "yes" if use else "no",
                "used_in_mean": "yes" if use else "no",
                "direction": direction,
                "fdr_lt_0.05": "yes" if row["FDR"] < 0.05 and use else "no",
            }
        )
    return found


def _fmt(value: float) -> str:
    if isinstance(value, float) and math.isnan(value):
        return ""
    return f"{value:.10g}"


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row[key]) if isinstance(row[key], float) else row[key] for key in fields})


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else math.nan


def median(values: list[float]) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def permutation_mean(values: list[float], universe: list[float], n_perm: int, seed: int) -> dict:
    """Gene-label permutation of the mean logFC. Not a mouse-level test."""
    if not values:
        return {"perm_p_twosided": math.nan, "n_perm": 0, "null_mean": math.nan}
    obs = mean(values)
    k = len(values)
    rng = random.Random(seed)
    exceed = 0
    null_sum = 0.0
    for _ in range(n_perm):
        stat = mean(rng.sample(universe, k))
        null_sum += stat
        if abs(stat) >= abs(obs) - 1e-15:
            exceed += 1
    return {
        "perm_p_twosided": (exceed + 1) / (n_perm + 1),
        "n_perm": n_perm,
        "null_mean": null_sum / n_perm,
        "obs_mean": obs,
    }


def module_summary(name: str, genes: list[dict], universe: list[tuple[str, float]], seed: int) -> dict:
    used = [g for g in genes if g["used_in_mean"] == "yes"]
    values = [g["logFC"] for g in used]
    used_symbols = {g["symbol"] for g in used}
    present = [g for g in genes if g["in_table"] == "yes"]
    fdr_hits = [g for g in used if g["fdr_lt_0.05"] == "yes"]
    up = sum(1 for v in values if v > 0)
    down = sum(1 for v in values if v < 0)
    background = [logfc for symbol, logfc in universe if symbol not in used_symbols]
    mw_p = math.nan
    wx_p = math.nan
    if values and background:
        mw = mannwhitneyu(values, background, alternative="two-sided", method="asymptotic")
        mw_p = float(mw.pvalue)
    if len(values) >= 6 and any(v != 0 for v in values):
        # Signed-rank against 0. Gene-level, one contrast.
        try:
            wx = wilcoxon(values, alternative="two-sided", zero_method="wilcox")
            wx_p = float(wx.pvalue)
        except ValueError:
            wx_p = math.nan
    perm = permutation_mean(values, [logfc for _, logfc in universe], N_PERM, seed)
    return {
        "module": name,
        "n_panel": len(genes),
        "n_in_table": len(present),
        "n_expressed_logCPM_ge_0": len(used),
        "n_low_count_excluded": len(present) - len(used),
        "n_missing": len(genes) - len(present),
        "n_up": up,
        "n_down": down,
        "n_fdr_lt_0.05": len(fdr_hits),
        "fdr_hits": ",".join(g["alias"] for g in fdr_hits),
        "mean_logFC_KO_minus_WT": mean(values),
        "median_logFC_KO_minus_WT": median(values),
        "mannwhitney_p_vs_expressed_background": mw_p,
        "wilcoxon_p_vs_0": wx_p,
        "gene_label_perm_p": perm["perm_p_twosided"],
        "gene_label_perm_n": perm["n_perm"],
        "note": "gene-label tests on one deposited contrast; n_mice = 1 vs 1",
    }


def universe_logfc(rows: list[dict]) -> list[tuple[str, float]]:
    return [(row["symbol"], row["logFC"]) for row in rows if expressed(row)]


def plot(genes: list[dict], path_png: Path, path_pdf: Path) -> None:
    panels = [
        ("c-NHEJ", [g for g in genes if g["module"] == "c-NHEJ"], "#3C6E8F"),
        ("Sting1 / Cgas", [g for g in genes if g["module"] == "STING-pair"], "#C47B2B"),
        ("IFN / chemokine", [g for g in genes if g["module"] == "IFN-chemokine"], "#A33B3B"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 6.4), sharex=True)
    for ax, (title, subset, color) in zip(axes, panels):
        subset = list(reversed(subset))
        y = list(range(len(subset)))
        vals = [0.0 if math.isnan(g["logFC"]) else g["logFC"] for g in subset]
        colors = []
        for g, v in zip(subset, vals):
            if g["expressed"] != "yes":
                colors.append("#B0B0B0")
            elif g["fdr_lt_0.05"] == "yes":
                colors.append(color)
            else:
                colors.append(color)
        bars = ax.barh(y, vals, color=colors, edgecolor="white", height=0.72)
        for bar, g in zip(bars, subset):
            if g["expressed"] != "yes":
                bar.set_hatch("///")
                bar.set_facecolor("#D0D0D0")
            elif g["fdr_lt_0.05"] != "yes":
                bar.set_alpha(0.45)
        labels = []
        for g in subset:
            if g["symbol"] == g["alias"] or g["alias"].upper() == g["symbol"].upper():
                labels.append(g["alias"])
            else:
                labels.append(f"{g['alias']} ({g['symbol']})")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.axvline(0, color="#222222", linewidth=0.8)
        ax.set_title(title, fontsize=11, loc="left")
        ax.set_xlabel("logFC (KO − WT)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for i, g in enumerate(subset):
            if g["fdr_lt_0.05"] == "yes":
                x = g["logFC"]
                ax.text(x + (0.08 if x >= 0 else -0.08), i, "FDR<0.05", va="center", ha="left" if x >= 0 else "right", fontsize=6.5, color="#333333")
    fig.suptitle(
        "GSE50927 naive lung, no VILI  ·  Cldn4 KO vs WT  ·  n = 1 vs 1 GSM",
        fontsize=12,
        x=0.01,
        ha="left",
    )
    fig.text(
        0.01,
        0.01,
        "Solid = author FDR < 0.05. Pale = in the table, FDR ≥ 0.05. Hatched = logCPM < 0, excluded from module means. Author EdgeR, mm9.",
        fontsize=7.5,
        color="#333333",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    fig.savefig(path_png, dpi=160)
    fig.savefig(path_pdf)
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = download()
    digest = sha256(path)
    rows = load_table(path)
    by_sym = index_symbols(rows)
    anchor = by_sym["Cldn4"]
    if not (-6.2 < anchor["logFC"] < -5.9):
        raise SystemExit(f"Cldn4 logFC {anchor['logFC']} is not the locked KO-minus-WT sign")

    panels = {
        "c-NHEJ": lookup(NHEJ_CORE, by_sym),
        "NHEJ-extended": lookup(NHEJ_EXTENDED, by_sym),
        "STING-pair": lookup(STING_PAIR, by_sym),
        "STING-neighbor": lookup(STING_NEIGHBOR, by_sym),
        "IFN-chemokine": lookup(IFN_CHEMOKINE, by_sym),
        "TJ-control": lookup(TJ_CONTROL, by_sym),
    }
    # Ifnb1 is in two panels; keep both rows (module membership differs).
    gene_rows = []
    for name in ["c-NHEJ", "NHEJ-extended", "STING-pair", "STING-neighbor", "IFN-chemokine", "TJ-control"]:
        gene_rows.extend(panels[name])
    anchor_row = lookup([ANCHOR], by_sym)
    gene_rows.extend(anchor_row)

    universe = universe_logfc(rows)
    summaries = []
    seeds = {
        "c-NHEJ": SEED + 1,
        "NHEJ-extended": SEED + 2,
        "STING-pair": SEED + 3,
        "STING-neighbor": SEED + 4,
        "IFN-chemokine": SEED + 5,
        "TJ-control": SEED + 6,
    }
    for name, seed in seeds.items():
        summaries.append(module_summary(name, panels[name], universe, seed))

    gene_fields = [
        "module",
        "symbol",
        "alias",
        "entrez",
        "name",
        "logFC",
        "logCPM",
        "PValue",
        "FDR",
        "in_table",
        "expressed",
        "used_in_mean",
        "direction",
        "fdr_lt_0.05",
    ]
    write_tsv(TABLES / "gene_level.tsv", gene_rows, gene_fields)
    sum_fields = list(summaries[0].keys())
    write_tsv(TABLES / "module_summary.tsv", summaries, sum_fields)

    sting = {g["alias"]: g for g in panels["STING-pair"]}
    ifn = next(s for s in summaries if s["module"] == "IFN-chemokine")
    nhej = next(s for s in summaries if s["module"] == "c-NHEJ")
    one = {
        "dataset": "GSE50927",
        "contrast": "naive Cldn4 KO vs WT, no VILI",
        "gsm_high_Cldn4": "GSM1232580 WT no VILI",
        "gsm_low_Cldn4": "GSM1232581 Cldn4 KO no VILI",
        "n_high_vs_low": "1 vs 1",
        "logFC_definition": "KO minus WT",
        "Cldn4_logFC": anchor["logFC"],
        "Cldn4_FDR": anchor["FDR"],
        "NHEJ_core_mean_logFC": nhej["mean_logFC_KO_minus_WT"],
        "NHEJ_core_n_fdr": nhej["n_fdr_lt_0.05"],
        "NHEJ_core_perm_p": nhej["gene_label_perm_p"],
        "Cgas_Mb21d1_logFC": sting["Cgas"]["logFC"],
        "Cgas_FDR": sting["Cgas"]["FDR"],
        "Sting1_Tmem173_logFC": sting["Sting1"]["logFC"],
        "Sting1_FDR": sting["Sting1"]["FDR"],
        "IFN_chemokine_mean_logFC": ifn["mean_logFC_KO_minus_WT"],
        "IFN_chemokine_n_expressed": ifn["n_expressed_logCPM_ge_0"],
        "IFN_chemokine_n_fdr": ifn["n_fdr_lt_0.05"],
        "IFN_chemokine_perm_p": ifn["gene_label_perm_p"],
        "call": "IFN/chemokine output up after Cldn4 loss; NHEJ and Sting1/Cgas mRNA flat",
        "cancer": "no",
        "VILI_scored": "no",
    }
    write_tsv(TABLES / "one_row.tsv", [one], list(one.keys()))

    inventory = [
        {"item": "GEO series", "n": "1", "used": "yes", "note": "GSE50927 Kage/Borok 2014 PMID 25106430"},
        {"item": "baseline GSM WT no VILI", "n": "1", "used": "yes", "note": "GSM1232580; Cldn4-intact"},
        {"item": "baseline GSM Cldn4 KO no VILI", "n": "1", "used": "yes", "note": "GSM1232581; Cldn4-null"},
        {"item": "deposited pairwise n", "n": "1 vs 1", "used": "yes", "note": "author EdgeR collapses the design-text duplicates"},
        {"item": "SRA runs for these two GSM", "n": "4", "used": "no", "note": "SRX352050 and SRX352051; not a processed count matrix"},
        {"item": "VILI GSM", "n": "3", "used": "no", "note": "WT VILI, KO VILIlow, KO VILIhigh; out of this wave"},
        {"item": "series-matrix expression rows", "n": "0", "used": "no", "note": "no per-sample matrix"},
        {"item": "tumour / LUAD / ICI", "n": "0", "used": "no", "note": "whole lung, mixed 129S6 / C57BL/6 / BALB/c"},
    ]
    write_tsv(TABLES / "label_inventory.tsv", inventory, ["item", "n", "used", "note"])

    def json_safe(obj):
        if isinstance(obj, float):
            if math.isnan(obj) or math.isinf(obj):
                return None
            return obj
        if isinstance(obj, dict):
            return {key: json_safe(value) for key, value in obj.items()}
        if isinstance(obj, list):
            return [json_safe(value) for value in obj]
        return obj

    summary = {
        "source_url": URL,
        "source_sha256": digest,
        "n_genes_in_table": len(rows),
        "n_expressed_logCPM_ge_0": len(universe),
        "min_logCPM_for_mean": MIN_LOGCPM,
        "seed": SEED,
        "n_perm": N_PERM,
        "cldn4_logFC": anchor["logFC"],
        "cldn4_FDR": anchor["FDR"],
        "modules": summaries,
        "one_row": one,
    }
    (TABLES / "summary.json").write_text(json.dumps(json_safe(summary), indent=2) + "\n")
    samples = [
        {
            "gsm": "GSM1232580",
            "title": "WT no VILI",
            "genotype": "WT",
            "vili": "none",
            "role": "Cldn4-intact baseline",
            "in_this_contrast": "yes",
        },
        {
            "gsm": "GSM1232581",
            "title": "Cldn4 KO no VILI",
            "genotype": "Cldn4 KO",
            "vili": "none",
            "role": "Cldn4-null baseline",
            "in_this_contrast": "yes",
        },
        {
            "gsm": "GSM1232582",
            "title": "WT VILI",
            "genotype": "WT",
            "vili": "40 cmH2O, 2 h",
            "role": "injury companion",
            "in_this_contrast": "no",
        },
        {
            "gsm": "GSM1232583",
            "title": "Cldn4 KO VILIlow",
            "genotype": "Cldn4 KO",
            "vili": "40 cmH2O, 2 h",
            "role": "injury stratum, not a Cldn4-expression class",
            "in_this_contrast": "no",
        },
        {
            "gsm": "GSM1232584",
            "title": "Cldn4 KO VILIhigh",
            "genotype": "Cldn4 KO",
            "vili": "40 cmH2O, 2 h",
            "role": "injury stratum, not a Cldn4-expression class",
            "in_this_contrast": "no",
        },
    ]
    write_tsv(
        TABLES / "sample_annotation.tsv",
        samples,
        ["gsm", "title", "genotype", "vili", "role", "in_this_contrast"],
    )
    plot(
        panels["c-NHEJ"] + panels["STING-pair"] + panels["IFN-chemokine"],
        FIGURES / "fig1_baseline_nhej_sting_ifn.png",
        FIGURES / "fig1_baseline_nhej_sting_ifn.pdf",
    )
    print(json.dumps({"one_row": one, "modules": summaries}, indent=2))


if __name__ == "__main__":
    main()
