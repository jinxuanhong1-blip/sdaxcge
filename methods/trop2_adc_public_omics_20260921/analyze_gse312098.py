#!/usr/bin/env python3
"""Weak-analogy reanalysis of public GSE312098 (colon CX-1, not lung).

Contrast is IMMU132 (sacituzumab govitecan) vs control at 2 days, n=3 vs 3.
A second contrast, IMMU132+GSK2606414 vs control, is reported beside it and
is not an ADC-only comparison.

Expression is author FPKM. Statistic is Welch's t on log2(FPKM+1).
Duplicate gene symbols are collapsed to the protein-coding row with the
highest mean FPKM. Hallmark sets are MSigDB 2024.1.Hs.
"""

from __future__ import annotations

import csv
import gzip
import math
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
FPKM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312098/suppl/"
    "GSE312098_gene_fpkm.txt.gz"
)
GMT_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
    "h.all.v2024.1.Hs.symbols.gmt"
)

# Library names in the FPKM file match GEO sample descriptions X_1..X_12.
ARMS = {
    "control": ["X_1", "X_2", "X_3"],
    "IMMU132": ["X_4", "X_5", "X_6"],
    "GSK2606414": ["X_7", "X_8", "X_9"],
    "combination": ["X_10", "X_11", "X_12"],
}
CONTRASTS = [
    ("IMMU132_vs_control", "IMMU132", "control"),
    ("combination_vs_control", "combination", "control"),
]
KEY_GENES = [
    "TACSTD2",
    "CLDN4",
    "CLDN1",
    "CLDN3",
    "CLDN7",
    "OCLN",
    "TJP1",
    "F11R",
    "IFNG",
    "IFNGR1",
    "STAT1",
    "IRF1",
    "IRF7",
    "ISG15",
    "IFI6",
    "IFI27",
    "IFIT1",
    "IFITM1",
    "MX1",
    "OAS1",
    "OAS2",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "TAP1",
    "TAP2",
    "PSMB8",
    "PSMB9",
    "CD274",
]
HALLMARKS = [
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
]


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research-hunt"})
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as out:
        out.write(resp.read())


def load_fpkm(path: Path):
    raw = gzip.open(path, "rb").read()
    text = raw.decode("utf-16")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    return rows


def collapse_protein_coding(rows):
    """One symbol: protein-coding row with the highest mean FPKM."""
    cols = [c for arm in ARMS.values() for c in arm]
    best = {}
    for row in rows:
        if row.get("gene_biotype") != "protein_coding":
            continue
        symbol = (row.get("gene_name") or "").strip()
        if not symbol or symbol == "-":
            continue
        vals = []
        ok = True
        for c in cols:
            try:
                vals.append(float(row[c]))
            except ValueError:
                ok = False
                break
        if not ok:
            continue
        mean_fpkm = float(np.mean(vals))
        prev = best.get(symbol)
        if prev is None or mean_fpkm > prev[0]:
            best[symbol] = (mean_fpkm, row["gene_id"], np.asarray(vals, dtype=float))
    symbols = sorted(best)
    mat = np.vstack([best[s][2] for s in symbols])
    ensembl = {s: best[s][1] for s in symbols}
    col_index = {c: i for i, c in enumerate(cols)}
    return symbols, mat, ensembl, col_index


def arm_matrix(mat, col_index, arm: str) -> np.ndarray:
    idx = [col_index[c] for c in ARMS[arm]]
    return mat[:, idx]


def welch_log2(treat: np.ndarray, ctrl: np.ndarray):
    lt = np.log2(treat + 1.0)
    lc = np.log2(ctrl + 1.0)
    log2fc = lt.mean(axis=1) - lc.mean(axis=1)
    # SciPy returns nan when a gene has zero variance on both sides.
    res = stats.ttest_ind(lt, lc, axis=1, equal_var=False, alternative="two-sided")
    p = np.asarray(res.pvalue, dtype=float)
    return log2fc, p, lt.mean(axis=1), lc.mean(axis=1), treat.mean(axis=1), ctrl.mean(axis=1)


def bh(pvals: np.ndarray) -> np.ndarray:
    q = np.full(pvals.shape, np.nan)
    valid = np.isfinite(pvals)
    idx = np.flatnonzero(valid)
    if idx.size == 0:
        return q
    order = idx[np.argsort(pvals[idx])]
    m = order.size
    ranked = pvals[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.clip(ranked, 0, 1)
    q[order] = ranked
    return q


def load_gmt(path: Path) -> dict[str, set[str]]:
    sets = {}
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        sets[parts[0]] = set(parts[2:])
    return sets


def fmt(x, digits=4):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "NA"
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    fpkm_path = DATA / "GSE312098_gene_fpkm.txt.gz"
    gmt_path = DATA / "h.all.v2024.1.Hs.symbols.gmt"
    download(FPKM_URL, fpkm_path)
    download(GMT_URL, gmt_path)

    rows = load_fpkm(fpkm_path)
    symbols, mat, ensembl, col_index = collapse_protein_coding(rows)
    symbol_index = {s: i for i, s in enumerate(symbols)}
    hallmark = load_gmt(gmt_path)

    sample_rows = []
    for gsm, title, arm, col in [
        ("GSM9337711", "CX-1 cells, Control, 2Day, rep1", "control", "X_1"),
        ("GSM9337712", "CX-1 cells, Control, 2Day, rep2", "control", "X_2"),
        ("GSM9337713", "CX-1 cells, Control, 2Day, rep3", "control", "X_3"),
        ("GSM9337714", "CX-1 cells, IMMU132, 2Day, rep1", "IMMU132", "X_4"),
        ("GSM9337715", "CX-1 cells, IMMU132, 2Day, rep2", "IMMU132", "X_5"),
        ("GSM9337716", "CX-1 cells, IMMU132, 2Day, rep3", "IMMU132", "X_6"),
        ("GSM9337717", "CX-1 cells, GSK2606414, 2Day, rep1", "GSK2606414", "X_7"),
        ("GSM9337718", "CX-1 cells, GSK2606414, 2Day, rep2", "GSK2606414", "X_8"),
        ("GSM9337719", "CX-1 cells, GSK2606414, 2Day, rep3", "GSK2606414", "X_9"),
        ("GSM9337720", "CX-1 cells, Combination, 2Day, rep1", "combination", "X_10"),
        ("GSM9337721", "CX-1 cells, Combination, 2Day, rep2", "combination", "X_11"),
        ("GSM9337722", "CX-1 cells, Combination, 2Day, rep3", "combination", "X_12"),
    ]:
        sample_rows.append(
            {
                "gsm": gsm,
                "title": title,
                "arm": arm,
                "fpkm_column": col,
                "cell_line": "CX-1",
                "tissue": "colorectal cancer cell line",
                "time": "48 h (series title: 2Day)",
                "dose_protocol": "3 ug/ml IMMU132 and/or 3 uM GSK2606414 for 48 h",
            }
        )
    write_tsv(TABLES / "gse312098_samples.tsv", sample_rows)

    key_rows = []
    set_rows = []
    contrast_cache = {}
    for contrast, treat_arm, ctrl_arm in CONTRASTS:
        treat = arm_matrix(mat, col_index, treat_arm)
        ctrl = arm_matrix(mat, col_index, ctrl_arm)
        log2fc, p, mean_lt, mean_lc, mean_t, mean_c = welch_log2(treat, ctrl)
        q_all = bh(p)
        contrast_cache[contrast] = (log2fc, p, q_all)
        for gene in KEY_GENES:
            i = symbol_index.get(gene)
            if i is None:
                key_rows.append(
                    {
                        "contrast": contrast,
                        "gene": gene,
                        "present": "no",
                        "ensembl": "",
                        "n_treat": treat.shape[1],
                        "n_control": ctrl.shape[1],
                        "mean_fpkm_treat": "",
                        "mean_fpkm_control": "",
                        "log2FC_log2fpkm1": "",
                        "welch_p": "",
                        "panel_bh_q": "",
                    }
                )
                continue
            key_rows.append(
                {
                    "contrast": contrast,
                    "gene": gene,
                    "present": "yes",
                    "ensembl": ensembl[gene],
                    "n_treat": treat.shape[1],
                    "n_control": ctrl.shape[1],
                    "mean_fpkm_treat": mean_t[i],
                    "mean_fpkm_control": mean_c[i],
                    "log2FC_log2fpkm1": log2fc[i],
                    "welch_p": p[i],
                    "panel_bh_q": "",
                }
            )
        # BH only inside the reported key-gene panel for this contrast.
        panel_p = []
        panel_ix = []
        for j, row in enumerate(key_rows):
            if row["contrast"] != contrast or row["present"] != "yes":
                continue
            panel_p.append(float(row["welch_p"]))
            panel_ix.append(j)
        panel_q = bh(np.asarray(panel_p, dtype=float))
        for j, qv in zip(panel_ix, panel_q):
            key_rows[j]["panel_bh_q"] = qv

        for set_name in HALLMARKS:
            members = hallmark[set_name]
            in_set = np.array([s in members for s in symbols])
            n_in = int(in_set.sum())
            n_out = int((~in_set).sum())
            a = log2fc[in_set]
            b = log2fc[~in_set]
            a = a[np.isfinite(a)]
            b = b[np.isfinite(b)]
            mw = stats.mannwhitneyu(a, b, alternative="two-sided")
            set_rows.append(
                {
                    "contrast": contrast,
                    "geneset": set_name,
                    "source": "MSigDB 2024.1.Hs hallmark",
                    "n_genes_in_matrix": n_in,
                    "n_background": n_out,
                    "median_log2FC": float(np.median(a)),
                    "mean_log2FC": float(np.mean(a)),
                    "fraction_log2FC_gt_0": float(np.mean(a > 0)),
                    "mw_p_vs_other_protein_coding": float(mw.pvalue),
                    "universe": "protein_coding symbols, one row each",
                }
            )

    # stringify numbers
    for row in key_rows:
        for k in ("mean_fpkm_treat", "mean_fpkm_control", "log2FC_log2fpkm1", "welch_p", "panel_bh_q"):
            if row[k] != "":
                row[k] = fmt(float(row[k]))
    for row in set_rows:
        for k in ("median_log2FC", "mean_log2FC", "fraction_log2FC_gt_0", "mw_p_vs_other_protein_coding"):
            row[k] = fmt(float(row[k]))

    write_tsv(TABLES / "gse312098_key_genes.tsv", key_rows)
    write_tsv(TABLES / "gse312098_genesets.tsv", set_rows)

    plot_key_genes(key_rows)
    print(f"protein_coding symbols: {len(symbols)}")
    print("wrote", TABLES)
    for row in key_rows:
        if row["gene"] in {"TACSTD2", "CLDN4", "IFNG", "STAT1", "CD274"}:
            print(row["contrast"], row["gene"], "log2FC", row["log2FC_log2fpkm1"], "p", row["welch_p"], "fpkm", row["mean_fpkm_control"], "->", row["mean_fpkm_treat"])
    for row in set_rows:
        print(row["contrast"], row["geneset"], "median", row["median_log2FC"], "mw", row["mw_p_vs_other_protein_coding"], "frac", row["fraction_log2FC_gt_0"])


def plot_key_genes(key_rows):
    contrast = "IMMU132_vs_control"
    genes = []
    vals = []
    ps = []
    for gene in KEY_GENES:
        row = next(r for r in key_rows if r["contrast"] == contrast and r["gene"] == gene)
        if row["present"] != "yes":
            continue
        genes.append(gene)
        vals.append(float(row["log2FC_log2fpkm1"]))
        try:
            ps.append(float(row["welch_p"]))
        except ValueError:
            ps.append(math.nan)
    order = np.argsort(vals)
    genes = [genes[i] for i in order]
    vals = [vals[i] for i in order]
    ps = [ps[i] for i in order]
    colors = ["#b45309" if g in {"TACSTD2", "CLDN4"} else "#1f4e79" for g in genes]
    fig, ax = plt.subplots(figsize=(7.2, 8.2))
    y = np.arange(len(genes))
    ax.barh(y, vals, color=colors, height=0.72)
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(genes, fontsize=8)
    ax.set_xlabel("log2FC  log2(FPKM+1), IMMU132 − control")
    ax.set_title("GSE312098 CX-1 colorectal line, 48 h\nIMMU132 vs control (n=3 vs 3), weak analogy")
    for yi, (v, p) in enumerate(zip(vals, ps)):
        label = f"p={p:.2g}" if math.isfinite(p) else ""
        x = v + (0.02 if v >= 0 else -0.02)
        ha = "left" if v >= 0 else "right"
        ax.text(x, yi, label, va="center", ha=ha, fontsize=6, color="#333333")
    fig.tight_layout()
    fig.savefig(FIGS / "gse312098_immu132_keygenes.png", dpi=160)
    fig.savefig(FIGS / "gse312098_immu132_keygenes.pdf")
    plt.close(fig)


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
