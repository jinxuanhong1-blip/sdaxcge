#!/usr/bin/env python3
"""GSE22493 CLDN4 siRNA vs CLDN4-overexpression, with HGNC symbol repair.

GPL10555 (BWH / Operon-style 60-mer) stores 2006-era ORF tokens. The prior
C4 and prerank runs treated those tokens as current HGNC symbols and only
rewrote G1P2->ISG15 and G1P3->IFI6. Genes that exist on the array only under
an obsolete token (NOD27, GRP58, ARTS-1, LRAP, cig5, PRKR, C1orf29, ...)
were scored as absent.

This script maps each probe token through the HGNC complete set
(approved symbol, else a unique previous symbol, else a unique alias).
It then tests the user KD thesis on the deposited log2(siRNA / CLDN4-OE)
values:

  IFN and APM go UP after CLDN4 loss
  other tight-junction genes go DOWN (CLDN4-high is the barrier)

CLDN4 itself is the perturbation check and is not inside the TJ test.

Primary gene statistic: median across arrays of the per-array probe median.
Primary set test: one-sided Mann-Whitney of that statistic vs other mapped
genes. n=3 two-colour arrays; ovarian SKOV-3; control is overexpression.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import re
import urllib.request
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "gse22493_cldn4_remap"
TAB = OUT / "tables"
FIG = OUT / "figures"
CACHE = Path(os.environ.get("GSE22493_CACHE", "/tmp/gse22493"))

SERIES_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493/"
    "matrix/GSE22493_series_matrix.txt.gz"
)
PLATFORM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10555/"
    "soft/GPL10555_family.soft.gz"
)
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/"
    "hgnc_complete_set.txt"
)

GSMS = ["GSM558700", "GSM558701", "GSM558702"]
TOKEN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9.\-]{0,19}$")

# User-named priority panel (private KD: these go up).
PRIORITY = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

# MHC-I / antigen-presentation machinery. Same 16 genes as the C4 slice.
APM = [
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "NLRC5",
    "TAP1",
    "TAP2",
    "TAPBP",
    "PSMB8",
    "PSMB9",
    "PSMB10",
    "ERAP1",
    "ERAP2",
    "CALR",
    "CANX",
    "PDIA3",
]

# Broader IFN / ISG set used in the prior CLDN4-loss slices.
IFN_IMMUNE = [
    "ISG15",
    "MX1",
    "MX2",
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL",
    "RSAD2",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IFIT5",
    "IFITM1",
    "IFITM2",
    "IFITM3",
    "IFI6",
    "IFI27",
    "IFI35",
    "IFI44",
    "IFI44L",
    "IFI16",
    "DDX58",
    "IFIH1",
    "DDX60",
    "XAF1",
    "HERC5",
    "USP18",
    "CMPK2",
    "BST2",
    "GBP1",
    "GBP2",
    "GBP4",
    "GBP5",
    "EIF2AK2",
    "ZBP1",
    "STAT1",
    "STAT2",
    "IRF1",
    "IRF7",
    "IRF9",
    "JAK2",
    "SOCS1",
    "SOCS3",
    "B2M",
    "TAP1",
    "TAP2",
    "PSMB8",
    "PSMB9",
    "NLRC5",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "HLA-F",
    "TAPBP",
    "CXCL10",
    "CXCL11",
    "CXCL9",
    "CCL5",
    "CCL2",
    "IL6",
    "TNF",
    "NFKB1",
    "RELA",
    "IL15",
    "IL32",
    "IFNB1",
    "IFNG",
    "IFNAR1",
    "IFNAR2",
    "IFNGR1",
    "IFNGR2",
    "IFNL1",
]

# Tight-junction core. CLDN4 is the perturbed gene and is tested apart.
TJ_CORE = [
    "CLDN1",
    "CLDN3",
    "CLDN7",
    "CLDN18",
    "OCLN",
    "TJP1",
    "TJP2",
    "TJP3",
    "MARVELD2",
    "MARVELD3",
    "F11R",
    "JAM2",
    "JAM3",
    "CGN",
    "CGNL1",
]

# Thesis direction for log2(siRNA / CLDN4-OE).
THESIS_UP = {"PRIORITY", "APM", "IFN_IMMUNE"}
THESIS_DOWN = {"TJ_CORE"}

# Tokens the prior run hard-coded. Everything else was left as the ORF string.
LEGACY_ALIAS = {"G1P2": "ISG15", "IFI15": "ISG15", "G1P3": "IFI6"}

# Remaps that must succeed or the run stops. These are the symbols the
# prior ORF-as-HGNC map missed or only caught by a two-gene hard-code.
REQUIRED_REMAPS = {
    "G1P2": "ISG15",
    "G1P3": "IFI6",
    "NOD27": "NLRC5",
    "GRP58": "PDIA3",
    "ARTS-1": "ERAP1",
    "LRAP": "ERAP2",
    "CIG5": "RSAD2",
    "PRKR": "EIF2AK2",
    "C1ORF29": "IFI44L",
    "TACSTD1": "EPCAM",
    "CLDN4": "CLDN4",
    "HLA-A": "HLA-A",
}


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    print(f"download {url}", flush=True)
    urllib.request.urlretrieve(url, dest)


def is_symbol_token(text: str) -> bool:
    s = text.strip()
    if not s or " " in s or len(s) > 20:
        return False
    return TOKEN_RE.match(s) is not None


def extract_token(orf: str, desc: str) -> tuple[str, str]:
    orf = (orf or "").strip()
    desc = (desc or "").strip()
    if is_symbol_token(orf):
        return orf, "ORF"
    if "--" in desc:
        prefix = desc.split("--", 1)[0].strip()
        if is_symbol_token(prefix):
            return prefix, "DESCRIPTION_PREFIX"
    return "", ""


def load_hgnc(path: Path) -> dict:
    approved: dict[str, str] = {}
    prev_map: dict[str, list[str]] = defaultdict(list)
    alias_map: dict[str, list[str]] = defaultdict(list)
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        n = 0
        for row in reader:
            n += 1
            sym = (row.get("symbol") or "").strip()
            if not sym:
                continue
            approved[sym.upper()] = sym
            for raw in (row.get("prev_symbol") or "").split("|"):
                tok = raw.strip()
                if tok:
                    prev_map[tok.upper()].append(sym)
            for raw in (row.get("alias_symbol") or "").split("|"):
                tok = raw.strip()
                if tok:
                    alias_map[tok.upper()].append(sym)
    return {
        "n_rows": n,
        "approved": approved,
        "prev": prev_map,
        "alias": alias_map,
    }


def resolve_token(token: str, hgnc: dict) -> tuple[str, str]:
    key = token.upper()
    approved = hgnc["approved"]
    if key in approved:
        return approved[key], "current"
    prevs = list(dict.fromkeys(hgnc["prev"].get(key, [])))
    aliases = list(dict.fromkeys(hgnc["alias"].get(key, [])))
    targets = list(dict.fromkeys(prevs + aliases))
    if len(targets) == 1:
        how = "prev_symbol" if targets[0] in prevs else "alias_symbol"
        return targets[0], how
    if len(targets) > 1:
        return "", "ambiguous"
    return "", "unmapped"


def legacy_symbol(token: str) -> str:
    return LEGACY_ALIAS.get(token.upper(), token.upper())


def load_platform(path: Path) -> list[dict]:
    rows = []
    with gzip.open(path, "rt", errors="replace") as fh:
        started = False
        header = None
        for line in fh:
            if line.startswith("!platform_table_begin"):
                started = True
                header = None
                continue
            if line.startswith("!platform_table_end"):
                break
            if not started:
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                fields = fields + [""] * (len(header) - len(fields))
            rec = dict(zip(header, fields))
            rows.append(rec)
    return rows


def load_matrix(path: Path) -> tuple[list[str], dict[str, list[float | None]]]:
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
            fields = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if fields[0] == "ID_REF":
                samples = fields[1:]
                continue
            row: list[float | None] = []
            for cell in fields[1:]:
                if cell == "":
                    row.append(None)
                else:
                    row.append(float(cell))
            while len(row) < len(samples):
                row.append(None)
            values[fields[0]] = row[: len(samples)]
    return samples, values


def annotate_probes(platform: list[dict], hgnc: dict) -> list[dict]:
    out = []
    for rec in platform:
        probe = (rec.get("ID") or "").strip()
        if not probe:
            continue
        token, token_from = extract_token(rec.get("ORF", ""), rec.get("DESCRIPTION", ""))
        if token:
            symbol, how = resolve_token(token, hgnc)
            legacy = legacy_symbol(token)
        else:
            symbol, how, legacy = "", "no_token", ""
        out.append(
            {
                "probe": probe,
                "token": token,
                "token_from": token_from,
                "symbol": symbol,
                "map_class": how,
                "legacy_symbol": legacy,
                "description": (rec.get("DESCRIPTION") or "").strip(),
            }
        )
    return out


def collapse_genes(
    probes: list[dict],
    matrix: dict[str, list[float | None]],
    symbol_key: str,
) -> pd.DataFrame:
    by_gene: dict[str, list[dict]] = defaultdict(list)
    for rec in probes:
        sym = rec[symbol_key]
        if sym:
            by_gene[sym].append(rec)
    rows = []
    for sym, plist in by_gene.items():
        per_array: list[float | None] = []
        for j in range(len(GSMS)):
            vals = []
            for rec in plist:
                row = matrix.get(rec["probe"])
                if row is None or row[j] is None or not math.isfinite(row[j]):
                    continue
                vals.append(row[j])
            per_array.append(float(np.median(vals)) if vals else None)
        finite = [v for v in per_array if v is not None]
        tokens = sorted({p["token"].upper() for p in plist if p["token"]})
        rows.append(
            {
                "symbol": sym,
                "n_probes": len(plist),
                "n_arrays": len(finite),
                "tokens": "|".join(tokens),
                "GSM558700": per_array[0],
                "GSM558701": per_array[1],
                "GSM558702": per_array[2],
                "median_log2": float(np.median(finite)) if finite else np.nan,
                "mean_log2": float(np.mean(finite)) if finite else np.nan,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values("symbol").reset_index(drop=True)


def measured(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    return df[df["n_arrays"] >= 2].copy()


def gene_stat(df: pd.DataFrame, symbol: str, col: str = "median_log2") -> float:
    hit = df[df["symbol"] == symbol]
    if hit.empty:
        return np.nan
    return float(hit.iloc[0][col])


def panel_values_from_rows(rows: list[dict], panel: str, col: str) -> list[float]:
    vals = []
    seen = set()
    for r in rows:
        if r["panel"] != panel or r["symbol"] in seen:
            continue
        seen.add(r["symbol"])
        v = r.get(col, np.nan)
        if v is None or not isinstance(v, (int, float)) or not np.isfinite(v):
            continue
        # Primary columns are already blank unless the gene is measured.
        if col in ("median_log2", "mean_log2") and r["measured"] != "yes":
            continue
        vals.append(float(v))
    return vals


def safe_wilcoxon(vals: list[float], alternative: str) -> float:
    x = np.asarray(vals, dtype=float)
    x = x[np.isfinite(x)]
    x = x[x != 0]
    if len(x) < 5:
        return np.nan
    try:
        res = stats.wilcoxon(x, alternative=alternative, zero_method="wilcox")
    except ValueError:
        return np.nan
    return float(res.pvalue)


def safe_mw(panel: list[float], background: list[float], alternative: str) -> float:
    a = np.asarray(panel, dtype=float)
    b = np.asarray(background, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 5 or len(b) < 20:
        return np.nan
    res = stats.mannwhitneyu(a, b, alternative=alternative)
    return float(res.pvalue)


def bh(pvals: list[float]) -> list[float]:
    arr = np.asarray(pvals, dtype=float)
    out = np.full(len(arr), np.nan)
    ok = np.isfinite(arr)
    if ok.sum() == 0:
        return out.tolist()
    order = np.argsort(arr[ok])
    ranked = arr[ok][order]
    m = len(ranked)
    q = ranked * m / (np.arange(1, m + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out[np.flatnonzero(ok)[order]] = q
    return out.tolist()


def direction(x: float) -> str:
    if not np.isfinite(x):
        return "NA"
    if x > 0:
        return "UP"
    if x < 0:
        return "DOWN"
    return "ZERO"


def fmt(x: float, nd: int = 3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{nd}f}"


def fmt_p(x: float) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def call_panel(median_effect: float, p_mw: float, thesis: str) -> str:
    if not np.isfinite(median_effect):
        return "not_testable"
    sign_ok = (thesis == "UP" and median_effect > 0) or (
        thesis == "DOWN" and median_effect < 0
    )
    if not np.isfinite(p_mw):
        return "same_direction_untested" if sign_ok else "opposite_thesis"
    if sign_ok and p_mw < 0.05:
        return "matches_thesis"
    if sign_ok:
        return "same_direction_not_significant"
    return "opposite_thesis"


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            out = {}
            for k in fields:
                v = row.get(k, "")
                if isinstance(v, float):
                    out[k] = "" if not math.isfinite(v) else f"{v:.6g}"
                else:
                    out[k] = v
            w.writerow(out)


def array_diagnostics(genes: pd.DataFrame) -> dict:
    m = measured(genes)
    cols = GSMS
    mat = m[cols].to_numpy(dtype=float)
    # Spearman on genes present on all three arrays.
    complete = mat[np.isfinite(mat).all(axis=1)]
    corr = {}
    for i, a in enumerate(cols):
        for j, b in enumerate(cols):
            if j <= i:
                continue
            rho, p = stats.spearmanr(complete[:, i], complete[:, j])
            corr[f"{a}_vs_{b}"] = {"rho": float(rho), "p": float(p), "n": int(complete.shape[0])}
    medians = {}
    for j, g in enumerate(cols):
        col = mat[:, j]
        col = col[np.isfinite(col)]
        medians[g] = float(np.median(col)) if len(col) else np.nan
    return {"spearman_complete_genes": corr, "array_median_log2": medians, "n_genes_ge2": int(len(m))}


def make_figure(panel_rows: list[dict], tests: list[dict]) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    show = []
    seen = set()
    preferred = ["PRIORITY", "APM", "TJ_CORE", "PERTURB"]
    for panel_name in preferred:
        for r in panel_rows:
            if r["panel"] != panel_name or r["measured"] != "yes" or r["symbol"] in seen:
                continue
            seen.add(r["symbol"])
            show.append(r)
    order = {s: i for i, s in enumerate(PRIORITY + [g for g in APM if g not in PRIORITY] + TJ_CORE + ["CLDN4"])}
    show.sort(key=lambda r: (0 if r["panel"] == "PRIORITY" else 1 if r["panel"] == "APM" else 2 if r["panel"] == "TJ_CORE" else 3, order.get(r["symbol"], 99)))

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 8.2), gridspec_kw={"width_ratios": [1.45, 1]})

    ax = axes[0]
    y = np.arange(len(show))
    colors = []
    for r in show:
        if r["panel"] == "PERTURB":
            colors.append("#444444")
        elif r["panel"] == "TJ_CORE":
            colors.append("#2c7fb8")
        elif r["panel"] == "APM" and r["symbol"] not in PRIORITY:
            colors.append("#e6550d")
        else:
            colors.append("#d95f0e")
    vals = [r["median_log2"] for r in show]
    ax.barh(y, vals, color=colors, height=0.72, zorder=2)
    for i, r in enumerate(show):
        xs = [r["GSM558700"], r["GSM558701"], r["GSM558702"]]
        ax.scatter(xs, [i] * 3, s=12, c="black", zorder=3)
    ax.axvline(0, color="black", lw=0.8)
    labels = []
    for r in show:
        lab = r["symbol"]
        if r.get("hgnc_symbol") and r["hgnc_symbol"] != r["symbol"]:
            lab = f"{r['symbol']} ({r['hgnc_symbol']})"
        elif r["remapped"] == "yes":
            lab = f"{r['symbol']}  [{r['tokens']}]"
        labels.append(lab)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("median log2 (CLDN4 siRNA / CLDN4-OE)")
    ax.set_title("Gene-level ratios after HGNC remap")
    ax.tick_params(axis="x", labelsize=8)

    ax2 = axes[1]
    # One row per primary test.
    tshow = [t for t in tests if t["stat"] == "median_log2" and t["panel"] in ("PRIORITY", "APM", "IFN_IMMUNE", "TJ_CORE")]
    yy = np.arange(len(tshow))
    for i, t in enumerate(tshow):
        col = "#1b9e77" if t["call"].startswith("matches") else ("#7570b3" if "same_direction" in t["call"] else "#d95f02")
        ax2.scatter([t["set_median"]], [i], s=60, c=col, zorder=3)
        ax2.plot([t["set_median"]], [i], marker="o", color=col)
    ax2.axvline(0, color="black", lw=0.8)
    ax2.set_yticks(yy)
    ylabels = []
    for t in tshow:
        ylabels.append(
            f"{t['panel']}\n{t['thesis']}  {t['n_up']}/{t['n_measured']} up\nMW {fmt_p(t['mw_p'])}"
        )
    ax2.set_yticklabels(ylabels, fontsize=8)
    ax2.invert_yaxis()
    ax2.set_xlabel("set median of gene medians")
    ax2.set_title("Thesis test vs background")
    fig.tight_layout()
    fig.savefig(FIG / "fig_ifn_apm_tj.png", dpi=160)
    fig.savefig(FIG / "fig_ifn_apm_tj.pdf")
    plt.close(fig)


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    series = CACHE / "GSE22493_series_matrix.txt.gz"
    platform = CACHE / "GPL10555_family.soft.gz"
    hgnc_path = CACHE / "hgnc_complete_set.txt"
    download(SERIES_URL, series)
    download(PLATFORM_URL, platform)
    download(HGNC_URL, hgnc_path)

    hgnc = load_hgnc(hgnc_path)
    for token, expect in REQUIRED_REMAPS.items():
        got, how = resolve_token(token, hgnc)
        if got != expect:
            raise SystemExit(f"HGNC remap failed: {token} -> {got} ({how}), expected {expect}")

    plat = load_platform(platform)
    samples, matrix = load_matrix(series)
    if samples != GSMS:
        raise SystemExit(f"Unexpected samples: {samples}")

    probes = annotate_probes(plat, hgnc)
    n_token = sum(1 for p in probes if p["token"])
    n_current = sum(1 for p in probes if p["map_class"] == "current")
    n_prev = sum(1 for p in probes if p["map_class"] == "prev_symbol")
    n_alias = sum(1 for p in probes if p["map_class"] == "alias_symbol")
    n_amb = sum(1 for p in probes if p["map_class"] == "ambiguous")
    n_unmap = sum(1 for p in probes if p["map_class"] == "unmapped")
    n_notok = sum(1 for p in probes if p["map_class"] == "no_token")
    approved_upper = set(hgnc["approved"])
    n_legacy_not_current = sum(
        1 for p in probes if p["legacy_symbol"] and p["legacy_symbol"].upper() not in approved_upper
    )

    genes = collapse_genes(probes, matrix, "symbol")
    legacy_genes = collapse_genes(probes, matrix, "legacy_symbol")
    genes_m = measured(genes)
    legacy_m = measured(legacy_genes)
    gene_ix = genes.set_index("symbol", drop=False)
    legacy_ix = set(legacy_m["symbol"])

    def hgnc_symbol(query: str) -> tuple[str, str]:
        got, how = resolve_token(query, hgnc)
        if got:
            return got, how
        return query, "as_written"

    diag = array_diagnostics(genes)

    # Probe rows for audited genes.
    audit_symbols = set(PRIORITY + APM + IFN_IMMUNE + TJ_CORE + ["CLDN4", "EPCAM", "TACSTD2", "CD274"])
    token_to_symbol = {}
    for p in probes:
        if p["symbol"] in audit_symbols or p["legacy_symbol"] in audit_symbols:
            token_to_symbol.setdefault(p["symbol"], set()).add(p["token"].upper())

    probe_rows = []
    for p in probes:
        if p["symbol"] not in audit_symbols and p["legacy_symbol"] not in audit_symbols:
            continue
        rowv = matrix.get(p["probe"], [None, None, None])
        probe_rows.append(
            {
                "probe": p["probe"],
                "token": p["token"],
                "token_from": p["token_from"],
                "symbol": p["symbol"],
                "map_class": p["map_class"],
                "legacy_symbol": p["legacy_symbol"],
                "remapped": "yes" if p["symbol"] and p["symbol"].upper() != p["token"].upper() else "no",
                "GSM558700": rowv[0],
                "GSM558701": rowv[1],
                "GSM558702": rowv[2],
                "description": p["description"][:180],
            }
        )

    panels = {
        "PRIORITY": PRIORITY,
        "APM": APM,
        "IFN_IMMUNE": IFN_IMMUNE,
        "TJ_CORE": TJ_CORE,
        "PERTURB": ["CLDN4"],
    }
    def subset_median(r: pd.Series, cols: list[str]) -> float:
        vals = [r[c] for c in cols if pd.notna(r[c]) and np.isfinite(r[c])]
        if len(vals) < len(cols):
            return np.nan
        return float(np.median(vals))

    panel_gene_rows = []
    for name, glist in panels.items():
        thesis = "UP" if name in THESIS_UP else ("DOWN" if name in THESIS_DOWN else "DOWN_expected_for_CLDN4")
        for g in glist:
            canon, canon_how = hgnc_symbol(g)
            if canon not in gene_ix.index:
                panel_gene_rows.append(
                    {
                        "panel": name,
                        "symbol": g,
                        "hgnc_symbol": canon,
                        "hgnc_class": canon_how,
                        "measured": "no",
                        "n_probes": 0,
                        "n_arrays": 0,
                        "tokens": "",
                        "remapped": "no",
                        "legacy_measured": "yes" if g in legacy_ix else "no",
                        "GSM558700": np.nan,
                        "GSM558701": np.nan,
                        "GSM558702": np.nan,
                        "median_log2": np.nan,
                        "mean_log2": np.nan,
                        "median_log2_700_702": np.nan,
                        "direction_median": "NA",
                        "thesis": thesis,
                    }
                )
                continue
            r = gene_ix.loc[canon]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            toks = [t for t in str(r["tokens"]).split("|") if t]
            remapped = "yes" if any(t.upper() != g.upper() for t in toks) or canon.upper() != g.upper() else "no"
            n_arr = int(r["n_arrays"])
            panel_gene_rows.append(
                {
                    "panel": name,
                    "symbol": g,
                    "hgnc_symbol": canon,
                    "hgnc_class": canon_how,
                    "measured": "yes" if n_arr >= 2 else "no",
                    "n_probes": int(r["n_probes"]),
                    "n_arrays": n_arr,
                    "tokens": r["tokens"],
                    "remapped": remapped,
                    "legacy_measured": "yes" if g in legacy_ix else "no",
                    "GSM558700": r["GSM558700"],
                    "GSM558701": r["GSM558701"],
                    "GSM558702": r["GSM558702"],
                    "median_log2": r["median_log2"] if n_arr >= 2 else np.nan,
                    "mean_log2": r["mean_log2"] if n_arr >= 2 else np.nan,
                    "median_log2_700_702": subset_median(r, ["GSM558700", "GSM558702"]),
                    "direction_median": direction(float(r["median_log2"])) if n_arr >= 2 else "NA",
                    "thesis": thesis,
                }
            )

    # Claudins present after remap (descriptive, not a second primary test).
    claudin_rows = []
    for _, r in genes_m.sort_values("symbol").iterrows():
        if str(r["symbol"]).startswith("CLDN"):
            claudin_rows.append(
                {
                    "symbol": r["symbol"],
                    "n_probes": int(r["n_probes"]),
                    "n_arrays": int(r["n_arrays"]),
                    "tokens": r["tokens"],
                    "median_log2": r["median_log2"],
                    "mean_log2": r["mean_log2"],
                    "direction_median": direction(float(r["median_log2"])),
                    "GSM558700": r["GSM558700"],
                    "GSM558701": r["GSM558701"],
                    "GSM558702": r["GSM558702"],
                }
            )

    # Background for the two-array sensitivity: gene median of GSM558700 and
    # GSM558702 only. GSM558701 anti-correlates with both (see diagnostics).
    sens_rows = []
    for _, r in genes.iterrows():
        v = subset_median(r, ["GSM558700", "GSM558702"])
        if np.isfinite(v):
            sens_rows.append((r["symbol"], v))
    sens_map = dict(sens_rows)

    def run_test(name: str, glist: list[str], stat_col: str, vals: list[float], bg: list[float], legacy_vals: list[float]) -> dict:
        thesis = "UP" if name in THESIS_UP else "DOWN"
        alt = "greater" if thesis == "UP" else "less"
        n_up = sum(v > 0 for v in vals)
        n_down = sum(v < 0 for v in vals)
        med = float(np.median(vals)) if vals else np.nan
        mean = float(np.mean(vals)) if vals else np.nan
        p_mw = safe_mw(vals, bg, alt)
        p_mw_two = safe_mw(vals, bg, "two-sided")
        p_w = safe_wilcoxon(vals, alt)
        n_nz = n_up + n_down
        if thesis == "UP":
            p_sign = float(stats.binomtest(n_up, n_nz, 0.5, alternative="greater").pvalue) if n_nz else np.nan
        else:
            p_sign = float(stats.binomtest(n_down, n_nz, 0.5, alternative="greater").pvalue) if n_nz else np.nan
        return {
            "panel": name,
            "stat": stat_col,
            "thesis": thesis,
            "n_in_list": len(glist),
            "n_measured": len(vals),
            "n_up": n_up,
            "n_down": n_down,
            "set_median": med,
            "set_mean": mean,
            "bg_median": float(np.median(bg)) if bg else np.nan,
            "mw_p": p_mw,
            "mw_p_two_sided": p_mw_two,
            "wilcoxon_p": p_w,
            "sign_p": p_sign,
            "legacy_n_measured": len(legacy_vals),
            "legacy_set_median": float(np.median(legacy_vals)) if legacy_vals else np.nan,
            "call": call_panel(med, p_mw, thesis),
        }

    tests = []
    for stat_col in ("median_log2", "mean_log2"):
        for name, glist in panels.items():
            if name == "PERTURB":
                continue
            vals = panel_values_from_rows(panel_gene_rows, name, stat_col)
            canon_set = {hgnc_symbol(g)[0] for g in glist}
            bg = [
                float(r[stat_col])
                for _, r in genes_m.iterrows()
                if r["symbol"] not in canon_set and np.isfinite(r[stat_col])
            ]
            # Legacy coverage still uses the query symbol, which is what the
            # old ORF-as-symbol map would have called the gene.
            leg_vals = []
            leg_ix = legacy_m.set_index("symbol")
            for g in glist:
                if g in leg_ix.index and np.isfinite(leg_ix.loc[g, stat_col]):
                    cell = leg_ix.loc[g, stat_col]
                    leg_vals.append(float(cell if not isinstance(cell, pd.Series) else cell.iloc[0]))
            tests.append(run_test(name, glist, stat_col, vals, bg, leg_vals))

    for name, glist in panels.items():
        if name == "PERTURB":
            continue
        vals = panel_values_from_rows(panel_gene_rows, name, "median_log2_700_702")
        canon_set = {hgnc_symbol(g)[0] for g in glist}
        bg = [v for s, v in sens_map.items() if s not in canon_set]
        tests.append(run_test(name, glist, "median_log2_700_702", vals, bg, []))

    # BH across the four primary tests on the pre-specified median statistic.
    for i, t in enumerate(tests):
        t["mw_q_bh4"] = np.nan
    for stat_name in ("median_log2", "median_log2_700_702"):
        primary_idx = [
            i
            for i, t in enumerate(tests)
            if t["stat"] == stat_name and t["panel"] in ("PRIORITY", "APM", "IFN_IMMUNE", "TJ_CORE")
        ]
        qvals = bh([tests[i]["mw_p"] for i in primary_idx])
        for i, q in zip(primary_idx, qvals):
            tests[i]["mw_q_bh4"] = q

    # Recovered panel genes: measured now, absent when the ORF string was the symbol.
    recovered = []
    for row in panel_gene_rows:
        if row["panel"] == "PERTURB":
            continue
        if row["measured"] == "yes" and row["legacy_measured"] == "no":
            recovered.append(row)

    # All-gene table is the analysis object; keep it.
    gene_out = genes_m.copy()
    gene_out["direction_median"] = gene_out["median_log2"].map(lambda x: direction(float(x)))

    cldn = genes[genes["symbol"] == "CLDN4"]
    cldn_row = cldn.iloc[0].to_dict() if not cldn.empty else {}

    summary = {
        "accession": "GSE22493",
        "platform": "GPL10555",
        "contrast": "CLDN4 lentiviral siRNA (Cy5) / CLDN4 overexpression control (Cy3)",
        "value": "deposited series-matrix VALUE, treated as log2(siRNA/OE) because GEO says normalized sample-to-control ratios and CLDN4 is the orientation check",
        "n_platform_rows": len(plat),
        "n_probes_in_matrix": len(matrix),
        "n_probes_with_token": n_token,
        "n_mapped_current": n_current,
        "n_mapped_prev_symbol": n_prev,
        "n_mapped_alias": n_alias,
        "n_ambiguous": n_amb,
        "n_unmapped_token": n_unmap,
        "n_no_token": n_notok,
        "n_legacy_token_not_current_hgnc": n_legacy_not_current,
        "n_genes_fixed_ge2": int(len(genes_m)),
        "n_genes_legacy_ge2": int(len(legacy_m)),
        "hgnc_rows": hgnc["n_rows"],
        "cldn4": {
            "n_probes": int(cldn_row.get("n_probes", 0) or 0),
            "n_arrays": int(cldn_row.get("n_arrays", 0) or 0),
            "GSM558700": cldn_row.get("GSM558700", np.nan),
            "GSM558701": cldn_row.get("GSM558701", np.nan),
            "GSM558702": cldn_row.get("GSM558702", np.nan),
            "median_log2": cldn_row.get("median_log2", np.nan),
            "mean_log2": cldn_row.get("mean_log2", np.nan),
        },
        "array_diagnostics": diag,
        "primary_stat": "median across arrays of per-array probe median; genes need >=2 arrays",
        "primary_test": "one-sided Mann-Whitney vs other mapped genes; BH across 4 panels",
        "tacstd2_on_platform": bool((genes["symbol"] == "TACSTD2").any()),
        "epcam_on_platform": bool((genes_m["symbol"] == "EPCAM").any()),
    }

    TAB.mkdir(parents=True, exist_ok=True)
    write_tsv(
        TAB / "panel_genes.tsv",
        panel_gene_rows,
        [
            "panel",
            "symbol",
            "hgnc_symbol",
            "hgnc_class",
            "measured",
            "n_probes",
            "n_arrays",
            "tokens",
            "remapped",
            "legacy_measured",
            "GSM558700",
            "GSM558701",
            "GSM558702",
            "median_log2",
            "mean_log2",
            "median_log2_700_702",
            "direction_median",
            "thesis",
        ],
    )
    write_tsv(
        TAB / "panel_tests.tsv",
        tests,
        [
            "panel",
            "stat",
            "thesis",
            "n_in_list",
            "n_measured",
            "n_up",
            "n_down",
            "set_median",
            "set_mean",
            "bg_median",
            "mw_p",
            "mw_p_two_sided",
            "mw_q_bh4",
            "wilcoxon_p",
            "sign_p",
            "legacy_n_measured",
            "legacy_set_median",
            "call",
        ],
    )
    write_tsv(
        TAB / "recovered_panel_genes.tsv",
        recovered,
        [
            "panel",
            "symbol",
            "tokens",
            "n_probes",
            "n_arrays",
            "median_log2",
            "mean_log2",
            "direction_median",
            "GSM558700",
            "GSM558701",
            "GSM558702",
        ],
    )
    write_tsv(
        TAB / "panel_probes.tsv",
        probe_rows,
        [
            "probe",
            "token",
            "token_from",
            "symbol",
            "map_class",
            "legacy_symbol",
            "remapped",
            "GSM558700",
            "GSM558701",
            "GSM558702",
            "description",
        ],
    )
    write_tsv(
        TAB / "claudins.tsv",
        claudin_rows,
        [
            "symbol",
            "n_probes",
            "n_arrays",
            "tokens",
            "median_log2",
            "mean_log2",
            "direction_median",
            "GSM558700",
            "GSM558701",
            "GSM558702",
        ],
    )
    gene_out.to_csv(TAB / "gene_median_log2.tsv", sep="\t", index=False, float_format="%.6g")

    def _jsonable(obj):
        if isinstance(obj, dict):
            return {k: _jsonable(v) for k, v in obj.items()}
        if isinstance(obj, float):
            return None if not math.isfinite(obj) else obj
        if isinstance(obj, (np.floating,)):
            v = float(obj)
            return None if not math.isfinite(v) else v
        if isinstance(obj, (np.integer,)):
            return int(obj)
        return obj

    (TAB / "summary.json").write_text(json.dumps(_jsonable(summary), indent=2) + "\n")
    make_figure(panel_gene_rows, tests)

    print(json.dumps(_jsonable(summary), indent=2))
    print("--- tests (median) ---")
    for t in tests:
        if t["stat"] == "median_log2":
            print(
                f"{t['panel']:12} n={t['n_measured']:3} up={t['n_up']:3} down={t['n_down']:3} "
                f"med={fmt(t['set_median'])} MW={fmt_p(t['mw_p'])} q={fmt_p(t['mw_q_bh4'])} "
                f"W={fmt_p(t['wilcoxon_p'])} sign={fmt_p(t['sign_p'])} "
                f"legacy_n={t['legacy_n_measured']} {t['call']}"
            )
    print("--- recovered ---")
    for r in recovered:
        print(f"{r['panel']:12} {r['symbol']:10} tokens={r['tokens']:20} med={fmt(r['median_log2'])} {r['direction_median']}")


if __name__ == "__main__":
    main()
