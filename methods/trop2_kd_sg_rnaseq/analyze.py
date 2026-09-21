#!/usr/bin/env python3
"""Score public TACSTD2-knockdown and sacituzumab RNA-seq for CLDN4, NHEJ, STING, IFN.

Downloads depositor-processed matrices from GEO, keeps the pre-specified
panel, and writes tables under results/trop2_kd_sg_rnaseq/. Raw matrices are
cached under results/.../data and are gitignored.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import sys
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import ALIASES, MOUSE_ENSEMBL, MOUSE_SYMBOL, PANEL, SETS  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "trop2_kd_sg_rnaseq"
DATA = OUT / "data"
TABLES = OUT / "tables"
FIGS = OUT / "figures"

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"cached {dest.name}", flush=True)
        return dest
    print(f"GET {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "trop2-public-rnaseq/1.0"})
    with urllib.request.urlopen(req, timeout=300) as resp, dest.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    print(f"  wrote {dest.stat().st_size} bytes", flush=True)
    return dest


def open_text(path: Path):
    raw = path.open("rb")
    if path.suffix == ".gz":
        raw = gzip.GzipFile(fileobj=raw)
    start = raw.read(4)
    # UTF-16 LE BOM, or a NUL in the first bytes of a UTF-16 header.
    if start.startswith(b"\xff\xfe") or (len(start) >= 2 and start[1] == 0):
        raw.close()
        if path.suffix == ".gz":
            fh = gzip.open(path, "rt", encoding="utf-16")
        else:
            fh = path.open("r", encoding="utf-16")
        return fh
    raw.close()
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def canon_symbol(name: str) -> str | None:
    if name is None:
        return None
    sym = name.strip().strip('"').split(".")[0]
    if not sym or sym.upper() in {"NA", "NAN", "NONE"}:
        return None
    key = sym.upper()
    key = ALIASES.get(key, key)
    if key in {g.upper() for g in PANEL}:
        # return the panel's canonical casing
        for g in PANEL:
            if g.upper() == key:
                return g
    return None


def norm_col(name: str) -> str:
    s = name.strip().strip('"')
    if s.endswith("_1"):
        s = s[:-2]
    return s


def extract_by_symbol(path: Path, symbol_col: str, sample_cols: list[str]) -> dict[str, dict[str, float]]:
    """Return panel gene -> {sample_col: value}. Highest-mean row wins ties."""
    want = set(sample_cols)
    best: dict[str, tuple[float, dict[str, float]]] = {}
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        header = [h.strip().strip('"') for h in header]
        idx = {h: i for i, h in enumerate(header)}
        missing = [c for c in sample_cols if c not in idx]
        if missing:
            # try normalized names (drop trailing _1)
            idx_norm = {norm_col(h): i for h, i in idx.items()}
            missing = [c for c in sample_cols if c not in idx and norm_col(c) not in idx_norm]
            if missing:
                raise SystemExit(f"{path.name}: missing columns {missing}; header={header[:20]}")
            col_index = {}
            for c in sample_cols:
                col_index[c] = idx[c] if c in idx else idx_norm[norm_col(c)]
        else:
            col_index = {c: idx[c] for c in sample_cols}
        if symbol_col not in idx:
            raise SystemExit(f"{path.name}: no symbol column {symbol_col}")
        si = idx[symbol_col]
        biotype_i = idx.get("Biotype", idx.get("gene_biotype"))
        for row in reader:
            if len(row) <= si:
                continue
            if biotype_i is not None and biotype_i < len(row):
                bt = row[biotype_i].strip().strip('"')
                if bt and bt not in {"protein_coding", "protein_coding_gene"}:
                    # keep protein-coding when the column is present; TACSTD2 etc. are coding
                    if "protein_coding" in bt or bt == "protein_coding":
                        pass
                    else:
                        continue
            sym = canon_symbol(row[si])
            if sym is None:
                continue
            vals = {}
            ok = True
            for c, i in col_index.items():
                if i >= len(row) or row[i].strip() in {"", "NA", "NaN"}:
                    ok = False
                    break
                vals[c] = float(row[i])
            if not ok:
                continue
            mean = float(np.mean(list(vals.values())))
            prev = best.get(sym)
            if prev is None or mean > prev[0]:
                best[sym] = (mean, vals)
    return {g: v[1] for g, v in best.items()}


def extract_mouse(path: Path, sample_cols: list[str]) -> dict[str, dict[str, float]]:
    ens_to_human = {}
    for human, mouse in MOUSE_SYMBOL.items():
        ens = MOUSE_ENSEMBL[mouse]
        ens_to_human[ens] = human
    found: dict[str, dict[str, float]] = {}
    with open_text(path) as fh:
        reader = csv.reader(fh)
        header = [h.strip().strip('"') for h in next(reader)]
        # first field may be an empty name
        idx = {h: i for i, h in enumerate(header)}
        missing = [c for c in sample_cols if c not in idx]
        if missing:
            raise SystemExit(f"{path.name}: missing {missing}; header={header}")
        for row in reader:
            if not row:
                continue
            ens = row[0].strip().strip('"').split(".")[0]
            human = ens_to_human.get(ens)
            if human is None:
                continue
            found[human] = {c: float(row[idx[c]]) for c in sample_cols}
    return found


def log2p(x: np.ndarray) -> np.ndarray:
    return np.log2(x + 1.0)


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    finite = np.where(np.isfinite(p))[0]
    if len(finite) == 0:
        return out.tolist()
    order = finite[np.argsort(p[finite])]
    m = len(order)
    adj = np.empty(m)
    running = 1.0
    for rank in range(m - 1, -1, -1):
        val = p[order[rank]] * m / (rank + 1)
        running = min(running, val)
        adj[rank] = min(running, 1.0)
    for i, idx in enumerate(order):
        out[idx] = adj[i]
    return out.tolist()


def welch(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    if np.allclose(a, a[0]) and np.allclose(b, b[0]):
        return 1.0 if np.isclose(a[0], b[0]) else float("nan")
    res = stats.ttest_ind(a, b, equal_var=False, alternative="two-sided")
    return float(res.pvalue)


def mwu(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 1 or len(b) < 1:
        return float("nan")
    try:
        res = stats.mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        return float("nan")
    return float(res.pvalue)


def paired_t(d: np.ndarray) -> float:
    if len(d) < 2 or np.allclose(d, d[0]) and np.isclose(d[0], 0):
        return 1.0 if len(d) and np.allclose(d, 0) else float("nan")
    if np.allclose(d, d[0]) and not np.isclose(d[0], 0):
        # zero variance, nonzero mean: t is infinite
        return 0.0
    res = stats.ttest_1samp(d, popmean=0.0, alternative="two-sided")
    return float(res.pvalue)


def wilcoxon(d: np.ndarray) -> float:
    if len(d) < 1 or np.allclose(d, 0):
        return 1.0 if len(d) else float("nan")
    try:
        res = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox")
    except ValueError:
        return float("nan")
    return float(res.pvalue)


def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return ""
    return f"{x:.{nd}g}"


def analyze_unpaired(expr: dict[str, dict[str, float]], pert: list[str], ctrl: list[str], genes: list[str]):
    """expr values are linear (FPKM or normalized counts)."""
    gene_rows = []
    log_samples = {}  # gene -> array aligned to pert+ctrl? we'll build per group
    for g in genes:
        if g not in expr:
            gene_rows.append({"gene": g, "present": 0})
            continue
        a = log2p(np.array([expr[g][s] for s in pert], dtype=float))
        b = log2p(np.array([expr[g][s] for s in ctrl], dtype=float))
        gene_rows.append(
            {
                "gene": g,
                "present": 1,
                "mean_log2_pert": float(np.mean(a)),
                "mean_log2_ctrl": float(np.mean(b)),
                "log2fc": float(np.mean(a) - np.mean(b)),
                "welch_p": welch(a, b),
                "mwu_p": mwu(a, b),
                "pert_values_log2": ";".join(f"{v:.4f}" for v in a),
                "ctrl_values_log2": ";".join(f"{v:.4f}" for v in b),
            }
        )
        log_samples[g] = (a, b)
    pvals = [r.get("welch_p", float("nan")) for r in gene_rows]
    qvals = bh(pvals)
    for r, q in zip(gene_rows, qvals):
        r["welch_q_bh"] = q
    return gene_rows, log_samples


def analyze_paired(expr, pairs: list[tuple[str, str]], genes: list[str]):
    """pairs are (pert_sample, ctrl_sample)."""
    gene_rows = []
    log_samples = {}
    for g in genes:
        if g not in expr:
            gene_rows.append({"gene": g, "present": 0})
            continue
        a = log2p(np.array([expr[g][p] for p, _ in pairs], dtype=float))
        b = log2p(np.array([expr[g][c] for _, c in pairs], dtype=float))
        d = a - b
        gene_rows.append(
            {
                "gene": g,
                "present": 1,
                "mean_log2_pert": float(np.mean(a)),
                "mean_log2_ctrl": float(np.mean(b)),
                "log2fc": float(np.mean(d)),
                "welch_p": paired_t(d),
                "mwu_p": wilcoxon(d),
                "pert_values_log2": ";".join(f"{v:.4f}" for v in a),
                "ctrl_values_log2": ";".join(f"{v:.4f}" for v in b),
                "pair_delta_log2": ";".join(f"{v:.4f}" for v in d),
            }
        )
        log_samples[g] = (a, b)
    pvals = [r.get("welch_p", float("nan")) for r in gene_rows]
    qvals = bh(pvals)
    for r, q in zip(gene_rows, qvals):
        r["welch_q_bh"] = q
    return gene_rows, log_samples


def set_summary(gene_rows, log_samples, set_name, members, paired: bool):
    present = [g for g in members if g in log_samples]
    if not present:
        return {
            "set": set_name,
            "n_genes": len(members),
            "n_present": 0,
        }
    # sample x gene log matrix
    A = np.vstack([log_samples[g][0] for g in present])  # genes x pert samples
    B = np.vstack([log_samples[g][1] for g in present])
    # z-score each gene across the contrast samples, then mean across genes
    both = np.concatenate([A, B], axis=1)
    mu = both.mean(axis=1, keepdims=True)
    sd = both.std(axis=1, keepdims=True, ddof=1)
    sd[sd == 0] = np.nan
    Z = (both - mu) / sd
    # if a gene has no variance, drop it from the score
    keep = np.isfinite(Z).all(axis=1)
    if not keep.any():
        score_p = score_c = None
        score_fc = float("nan")
        p_s = float("nan")
        p_u = float("nan")
    else:
        Zk = Z[keep]
        n_pert = A.shape[1]
        score_p = np.nanmean(Zk[:, :n_pert], axis=0)
        score_c = np.nanmean(Zk[:, n_pert:], axis=0)
        if paired:
            d = score_p - score_c
            score_fc = float(np.mean(d))
            p_s = paired_t(d)
            p_u = wilcoxon(d)
        else:
            score_fc = float(np.mean(score_p) - np.mean(score_c))
            p_s = welch(score_p, score_c)
            p_u = mwu(score_p, score_c)
    fcs = [r["log2fc"] for r in gene_rows if r["gene"] in present and r.get("present")]
    n_up = sum(1 for v in fcs if v > 0)
    n_down = sum(1 for v in fcs if v < 0)
    return {
        "set": set_name,
        "n_genes": len(members),
        "n_present": len(present),
        "n_scored": int(keep.sum()) if not (isinstance(keep, bool)) else 0,
        "mean_log2fc": float(np.mean(fcs)),
        "median_log2fc": float(np.median(fcs)),
        "n_up": n_up,
        "n_down": n_down,
        "score_delta": score_fc,
        "score_p": p_s,
        "score_rank_p": p_u,
        "genes": ",".join(present),
    }


def write_tsv(path: Path, rows: list[dict], cols: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore", delimiter="\t")
        w.writeheader()
        for r in rows:
            out = {}
            for c in cols:
                v = r.get(c, "")
                if isinstance(v, float):
                    out[c] = "" if math.isnan(v) else f"{v:.6g}"
                else:
                    out[c] = v
            w.writerow(out)


def heatmap(contrasts_gene_rows, path: Path):
    # rows = panel genes that are not only QC, columns = contrasts
    genes = [g for g in PANEL if g not in {"TACSTD2", "CDKN1A"}]
    # still show TACSTD2 and CDKN1A at the bottom as QC
    genes = genes + ["TACSTD2", "CDKN1A"]
    names = [c["id"] for c in contrasts_gene_rows]
    mat = np.full((len(genes), len(names)), np.nan)
    for j, c in enumerate(contrasts_gene_rows):
        by = {r["gene"]: r for r in c["genes"]}
        for i, g in enumerate(genes):
            rec = by.get(g)
            if rec and rec.get("present"):
                mat[i, j] = rec["log2fc"]
    fig_w = max(8, 1.3 * len(names) + 3)
    fig_h = max(8, 0.28 * len(genes) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    finite = mat[np.isfinite(mat)]
    lim = 2.0 if finite.size == 0 else float(np.nanpercentile(np.abs(finite), 95))
    lim = max(lim, 0.5)
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                ax.text(j, i, f"{mat[i, j]:+.2f}", ha="center", va="center", fontsize=6, color="black")
    ax.set_title("log2FC (perturbation − control), log2(x+1)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="log2FC")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def set_bars(set_rows, path: Path):
    # grouped bars of mean log2FC for the four asked sets
    contrasts = []
    for r in set_rows:
        if r["contrast"] not in contrasts:
            contrasts.append(r["contrast"])
    sets = ["CLDN4", "NHEJ", "STING", "IFN"]
    x = np.arange(len(contrasts))
    width = 0.18
    fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(contrasts) + 2), 4.8))
    for k, s in enumerate(sets):
        vals = []
        for c in contrasts:
            hit = next(r for r in set_rows if r["contrast"] == c and r["set"] == s)
            vals.append(hit.get("mean_log2fc", np.nan))
        ax.bar(x + (k - 1.5) * width, vals, width, label=s)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(contrasts, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("mean log2FC of set genes")
    ax.set_title("CLDN4 / NHEJ / STING / IFN")
    ax.legend(frameon=False, ncol=4)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    files = {
        "GSE245459": download(
            f"{GEO}/GSE245nnn/GSE245459/suppl/GSE245459_fpkm.anno.txt.gz",
            DATA / "GSE245459_fpkm.anno.txt.gz",
        ),
        "GSE334497": download(
            f"{GEO}/GSE334nnn/GSE334497/suppl/GSE334497_normalized_counts.csv.gz",
            DATA / "GSE334497_normalized_counts.csv.gz",
        ),
        "GSE312098": download(
            f"{GEO}/GSE312nnn/GSE312098/suppl/GSE312098_gene_fpkm.txt.gz",
            DATA / "GSE312098_gene_fpkm.txt.gz",
        ),
        "GSE311016": download(
            f"{GEO}/GSE311nnn/GSE311016/suppl/GSE311016_gene_fpkm.txt.gz",
            DATA / "GSE311016_gene_fpkm.txt.gz",
        ),
        "GSE304294": download(
            f"{GEO}/GSE304nnn/GSE304294/suppl/GSE304294_gene_fpkm.txt.gz",
            DATA / "GSE304294_gene_fpkm.txt.gz",
        ),
    }

    print("extract GSE245459", flush=True)
    skov = extract_by_symbol(
        files["GSE245459"],
        "GeneName",
        ["shNC1", "shNC2", "shNC3", "sh1", "sh2", "sh3", "shNCDDP1", "shNCDDP2", "shNCDDP3", "shDDP1", "shDDP2", "shDDP3"],
    )
    print("  genes", sorted(skov))

    print("extract GSE334497", flush=True)
    ko_libs = ["KO162", "KO164", "KO165", "KO172", "RESUB-KO163R"]
    wt_libs = ["RESUB-171R", "RESUB-170R", "RESUB-169R", "RESUB-168R", "control170"]
    mouse = extract_mouse(files["GSE334497"], ko_libs + wt_libs)
    print("  genes", sorted(mouse))

    print("extract GSE312098", flush=True)
    cx_ctrl = ["X_1", "X_2", "X_3"]
    cx_sg = ["X_4", "X_5", "X_6"]
    cx_combo = ["X_10", "X_11", "X_12"]
    cx = extract_by_symbol(files["GSE312098"], "gene_name", cx_ctrl + cx_sg + cx_combo)
    print("  genes", sorted(cx))

    print("extract GSE304294", flush=True)
    ky_ctrl = ["OX1_1", "OX1_2", "OX1_3"]
    ky_sg = ["OX2_1", "OX2_2"]
    ky_combo = ["OX4_1", "OX4_2", "OX4_3"]
    ky = extract_by_symbol(files["GSE304294"], "gene_name", ky_ctrl + ky_sg + ky_combo)
    print("  genes", sorted(ky))

    print("extract GSE311016", flush=True)
    pdx_pairs = [("T_36", "C_36"), ("T_82", "C_82"), ("T_83", "C_83"), ("T_114", "C_114"), ("T_196", "C_196")]
    pdx_cols = [c for pair in pdx_pairs for c in pair]
    pdx = extract_by_symbol(files["GSE311016"], "gene_name", pdx_cols)
    print("  genes", sorted(pdx))

    contrasts = [
        {
            "id": "SKOV3_shTACSTD2",
            "accession": "GSE245459",
            "class": "TACSTD2_KD",
            "system": "SKOV3 ovarian cell line",
            "contrast": "shTACSTD2 vs shNC (no cisplatin)",
            "n_pert": 3,
            "n_ctrl": 3,
            "unit": "library",
            "test": "unpaired",
            "matrix": "FPKM",
            "bulk_tumor": 0,
            "primary": 1,
            "expr": skov,
            "pert": ["sh1", "sh2", "sh3"],
            "ctrl": ["shNC1", "shNC2", "shNC3"],
        },
        {
            "id": "SKOV3_shTACSTD2_DDP",
            "accession": "GSE245459",
            "class": "TACSTD2_KD",
            "system": "SKOV3 ovarian cell line + cisplatin",
            "contrast": "shTACSTD2+cisplatin vs shNC+cisplatin",
            "n_pert": 3,
            "n_ctrl": 3,
            "unit": "library",
            "test": "unpaired",
            "matrix": "FPKM",
            "bulk_tumor": 0,
            "primary": 0,
            "expr": skov,
            "pert": ["shDDP1", "shDDP2", "shDDP3"],
            "ctrl": ["shNCDDP1", "shNCDDP2", "shNCDDP3"],
        },
        {
            "id": "4T1_Trop2KO",
            "accession": "GSE334497",
            "class": "TACSTD2_KO",
            "system": "4T1 mouse mammary tumor, bulk",
            "contrast": "Trop2 KO vs WT tumor",
            "n_pert": 5,
            "n_ctrl": 5,
            "unit": "tumor",
            "test": "unpaired",
            "matrix": "normalized counts",
            "bulk_tumor": 1,
            "primary": 1,
            "expr": mouse,
            "pert": ko_libs,
            "ctrl": wt_libs,
        },
        {
            "id": "CX1_IMMU132",
            "accession": "GSE312098",
            "class": "sacituzumab",
            "system": "CX-1 CRC cell line, 2 day",
            "contrast": "IMMU132 vs control",
            "n_pert": 3,
            "n_ctrl": 3,
            "unit": "library",
            "test": "unpaired",
            "matrix": "FPKM",
            "bulk_tumor": 0,
            "primary": 1,
            "expr": cx,
            "pert": cx_sg,
            "ctrl": cx_ctrl,
        },
        {
            "id": "CX1_IMMU132_GSK",
            "accession": "GSE312098",
            "class": "sacituzumab_combo",
            "system": "CX-1 CRC cell line, 2 day",
            "contrast": "IMMU132+GSK2606414 vs control",
            "n_pert": 3,
            "n_ctrl": 3,
            "unit": "library",
            "test": "unpaired",
            "matrix": "FPKM",
            "bulk_tumor": 0,
            "primary": 0,
            "expr": cx,
            "pert": cx_combo,
            "ctrl": cx_ctrl,
        },
        {
            "id": "KYSE30_IMMU132",
            "accession": "GSE304294",
            "class": "sacituzumab",
            "system": "KYSE30 ESCC cell line, 1 day",
            "contrast": "IMMU132 vs control",
            "n_pert": 2,
            "n_ctrl": 3,
            "unit": "library",
            "test": "unpaired",
            "matrix": "FPKM",
            "bulk_tumor": 0,
            "primary": 1,
            "expr": ky,
            "pert": ky_sg,
            "ctrl": ky_ctrl,
        },
        {
            "id": "KYSE30_IMMU132_IACS",
            "accession": "GSE304294",
            "class": "sacituzumab_combo",
            "system": "KYSE30 ESCC cell line, 1 day",
            "contrast": "IMMU132+IACS-010759 vs control",
            "n_pert": 3,
            "n_ctrl": 3,
            "unit": "library",
            "test": "unpaired",
            "matrix": "FPKM",
            "bulk_tumor": 0,
            "primary": 0,
            "expr": ky,
            "pert": ky_combo,
            "ctrl": ky_ctrl,
        },
        {
            "id": "CRC_PDX_IMMU132",
            "accession": "GSE311016",
            "class": "sacituzumab",
            "system": "CRC PDX, 29 day, paired",
            "contrast": "IMMU132 vs matched control",
            "n_pert": 5,
            "n_ctrl": 5,
            "unit": "PDX model",
            "test": "paired",
            "matrix": "FPKM",
            "bulk_tumor": 1,
            "primary": 1,
            "expr": pdx,
            "pairs": pdx_pairs,
        },
    ]

    gene_out = []
    set_out = []
    for c in contrasts:
        print("contrast", c["id"], flush=True)
        if c["test"] == "paired":
            grows, logs = analyze_paired(c["expr"], c["pairs"], PANEL)
        else:
            grows, logs = analyze_unpaired(c["expr"], c["pert"], c["ctrl"], PANEL)
        c["genes"] = grows
        for r in grows:
            gene_out.append(
                {
                    "contrast": c["id"],
                    "accession": c["accession"],
                    "class": c["class"],
                    "primary": c["primary"],
                    "bulk_tumor": c["bulk_tumor"],
                    **r,
                }
            )
        for set_name, members in SETS.items():
            sm = set_summary(grows, logs, set_name, members, paired=c["test"] == "paired")
            sm.update(
                {
                    "contrast": c["id"],
                    "accession": c["accession"],
                    "class": c["class"],
                    "primary": c["primary"],
                    "bulk_tumor": c["bulk_tumor"],
                    "n_pert": c["n_pert"],
                    "n_ctrl": c["n_ctrl"],
                    "test": c["test"],
                }
            )
            set_out.append(sm)
            print(
                f"  {set_name}: mean log2FC={sm.get('mean_log2fc')} scoreΔ={sm.get('score_delta')} p={sm.get('score_p')}",
                flush=True,
            )

    gene_cols = [
        "contrast",
        "accession",
        "class",
        "primary",
        "bulk_tumor",
        "gene",
        "present",
        "mean_log2_pert",
        "mean_log2_ctrl",
        "log2fc",
        "welch_p",
        "welch_q_bh",
        "mwu_p",
        "pert_values_log2",
        "ctrl_values_log2",
        "pair_delta_log2",
    ]
    set_cols = [
        "contrast",
        "accession",
        "class",
        "primary",
        "bulk_tumor",
        "test",
        "n_pert",
        "n_ctrl",
        "set",
        "n_genes",
        "n_present",
        "n_scored",
        "mean_log2fc",
        "median_log2fc",
        "n_up",
        "n_down",
        "score_delta",
        "score_p",
        "score_rank_p",
        "genes",
    ]
    write_tsv(TABLES / "gene_log2fc.tsv", gene_out, gene_cols)
    write_tsv(TABLES / "set_scores.tsv", set_out, set_cols)

    # compact primary-only pivot for the four sets
    primary_ids = [c["id"] for c in contrasts if c["primary"]]
    pivot = []
    for g in ["CLDN4", "TACSTD2", "CDKN1A"] + [x for s in ("NHEJ", "STING", "IFN") for x in SETS[s]]:
        row = {"gene": g}
        for cid in primary_ids:
            rec = next(r for r in gene_out if r["contrast"] == cid and r["gene"] == g)
            row[f"{cid}_log2fc"] = rec.get("log2fc", float("nan")) if rec.get("present") else float("nan")
            row[f"{cid}_p"] = rec.get("welch_p", float("nan")) if rec.get("present") else float("nan")
        pivot.append(row)
    write_tsv(
        TABLES / "primary_gene_pivot.tsv",
        pivot,
        ["gene"] + [f"{cid}_{k}" for cid in primary_ids for k in ("log2fc", "p")],
    )

    heatmap(contrasts, FIGS / "fig1_log2fc_heatmap.png")
    set_bars(set_out, FIGS / "fig2_set_mean_log2fc.png")

    summary = {
        "contrasts": [
            {k: c[k] for k in ("id", "accession", "class", "system", "contrast", "n_pert", "n_ctrl", "test", "matrix", "bulk_tumor", "primary")}
            for c in contrasts
        ],
        "sets": {k: v for k, v in SETS.items()},
        "note_oas1_mouse": "Mouse OAS1 ortholog in the 4T1 score is Oas1a (ENSMUSG00000052776).",
        "note_stats": "Unpaired Welch t and Mann-Whitney on log2(x+1). Paired CRC PDX uses one-sample t and Wilcoxon on within-model deltas. BH q is within the pre-specified panel of each contrast, not genome-wide. n=3 vs 3 Mann-Whitney cannot reach p<0.05.",
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2))
    print("done", flush=True)


if __name__ == "__main__":
    main()
